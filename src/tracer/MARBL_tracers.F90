!> A tracer package for tracers computed in the MARBL library
module MARBL_tracers

! This file is part of MOM6. See LICENSE.md for the license.

! Currently configured for use with marbl0.36.0
! https://github.com/marbl-ecosys/MARBL/releases/tag/marbl0.36.0
! (clone entire repo into pkg/MARBL)
#ifdef _USE_MARBL_TRACERS
use MOM_coms,            only : root_PE, broadcast
use MOM_diag_mediator,   only : diag_ctrl
use MOM_error_handler,   only : is_root_PE, MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_grid,            only : ocean_grid_type
use MOM_hor_index,       only : hor_index_type
use MOM_io,              only : file_exists, read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : query_initialized, MOM_restart_CS, register_restart_field
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init,   only : tracer_Z_init
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type
use MOM_diag_mediator,   only : register_diag_field, post_data!, safe_alloc_ptr

use MARBL_interface,              only : MARBL_interface_class
use MARBL_interface_public_types, only : marbl_diagnostics_type, marbl_saved_state_type

use coupler_types_mod,      only : coupler_type_set_data, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_MARBL_tracers, initialize_MARBL_tracers
public MARBL_tracers_column_physics, MARBL_tracers_surface_state
public MARBL_tracer_stock, MARBL_tracers_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Temporary type for diagnostic variables coming from MARBL
!! Allocate exactly one of field_[23]d
type :: temp_MARBL_diag
  integer :: id !< index into MOM diagnostic structure
  real, allocatable :: field_2d(:,:) !< memory for 2D field
  real, allocatable :: field_3d(:,:,:) !< memory for 3D field
end type temp_MARBL_diag

type :: saved_state_for_MARBL_type
  character(len=200) :: short_name !< name of variable being saved
  character(len=200) :: file_varname !< name of variable in restart file
  character(len=200) :: units !< variable units
  real, pointer :: field_2d(:,:) !< memory for 2D field
  real, pointer :: field_3d(:,:,:) !< memory for 3D field
end type saved_state_for_MARBL_type

!> All calls to MARBL are done via the interface class
type(MARBL_interface_class) :: MARBL_instances

!> The control structure for the MARBL tracer package
type, public :: MARBL_tracers_CS ; private
  integer :: ntr    !< The number of tracers that are actually used.
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the coupler.
  character(len=200) :: IC_file !< The file in which the age-tracer initial values cam be found.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL() !< The array of tracers used in this subroutine, in g m-3?

  integer, allocatable, dimension(:) :: ind_tr !< Indices returned by aof_set_coupler_flux if it is used and the
                                               !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(vardesc), allocatable :: tr_desc(:) !< Descriptions and metadata for the tracers
  logical :: tracers_may_reinit = .false. !< If true the tracers may be initialized if not found in a restart file

  !> Driver-specific parameters
  character(len=35) :: marbl_settings_file
  real :: atm_co2_const, atm_alt_co2_const
  real :: ndep_scale_factor

  !> Indices to the registered diagnostics and saved state match the indices used in MARBL
  type(temp_MARBL_diag), allocatable :: surface_flux_diags(:), interior_tendency_diags(:)
  type(saved_state_for_MARBL_type), allocatable :: surface_flux_saved_state(:), interior_tendency_saved_state(:)

  !> Surface fluxes returned from MARBL and passed to tracer_vertdiff
  !! tracer_vertdiff expects units concentrations times meters per second,
  !! but MARBL will return in cgs so we need to remember to convert
  real, allocatable :: STF(:,:,:) ! i, j, tracer

  integer :: u10_sqr_ind, sss_ind, sst_ind, ifrac_ind, dust_dep_ind, fe_dep_ind
  integer :: nox_flux_ind, nhy_flux_ind, atmpress_ind, xco2_ind, xco2_alt_ind
end type MARBL_tracers_CS

!> If we can post data column by column, all we need are integer
!! arrays for ids
! integer, allocatable :: id_surface_flux_diags(:) !, id_interior_tendency_diags(:)

!> Module parameters
  real, parameter :: cm_per_m = 100.  !< convert from m -> cm (MARBL is cgs)
  real, parameter :: g_per_kg = 1000. !< convert from kg -> g (MARBL is cgs)
  real, parameter :: m_per_cm = 0.01  !< convert from cm -> m
  real, parameter :: atm_per_Pa = 1./101325.  !< convert from Pa -> atm

contains

!> This subroutine is used to read marbl_in, configure MARBL accordingly, and then
!! call MARBL's initialization routine
subroutine configure_MARBL_tracers(GV, param_file, CS)
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(MARBL_tracers_CS),     pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module

#include "version_variable.h"
  character(len=40)  :: mdl = "MARBL_tracers" ! This module's name.
  character(len=256) :: log_message
  character(len=256) :: marbl_in_line(1)
  integer :: m, nz, marbl_settings_in, read_error
  nz = GV%ke
  marbl_settings_in = 615

  ! (1) Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MARBL_SETTINGS_FILE", CS%marbl_settings_file, &
                 "The name of a file from which to read the run-time "//&
                 "settings for MARBL.", default="marbl_in")
  call get_param(param_file, mdl, "ATM_CO2_CONST", CS%atm_co2_const, &
                 "Value to send to MARBL as xco2", default=284.317)
  call get_param(param_file, mdl, "ATM_ALT_CO2_CONST", CS%atm_alt_co2_const, &
                 "Value to send to MARBL as xco2_alt_co2", default=284.317)
  call get_param(param_file, mdl, "NDEP_SCALE_FACTOR", CS%ndep_scale_factor, &
                 "Scale factor applied to nitrogen deposition terms", default=1e5)

  ! (2) Read marbl settings file and call put_setting()

  ! (2a) only master task opens file
  if (is_root_PE()) then
     ! read the marbl_in into buffer
     open(unit=marbl_settings_in, file=CS%marbl_settings_file, iostat=read_error)
     if (read_error .ne. 0) then
        write(log_message, '(A, I0, 2A)') "IO ERROR ", read_error, &
              "opening namelist file : ", trim(CS%marbl_settings_file)
        call MOM_error(FATAL, log_message)
     end if
  end if

  ! (2b) master task reads file and broadcasts line-by-line
  marbl_in_line = ''
  do
    ! i. Read next line on master, iostat value out
    !    (Exit loop if read is not successful; either read error or end of file)
    if (is_root_PE()) read(marbl_settings_in, "(A)", iostat=read_error) marbl_in_line(1)
    call broadcast(read_error, root_PE())
    if (read_error .ne. 0) exit

    ! ii. Broadcast line just read in on root PE to all tasks
    call broadcast(marbl_in_line, 256, root_PE())

    ! iii. All tasks call put_setting (TODO: openMP blocks?)
    call marbl_instances%put_setting(marbl_in_line(1))
  end do

  ! (2c) we should always reach the EOF to capture the entire file...
  if (.not. is_iostat_end(read_error)) then
     write(log_message, '(3A, I0)') "IO ERROR reading ", trim(CS%marbl_settings_file), ": ", read_error
     call MOM_error(FATAL, log_message)
  else
     if (is_root_PE()) then
       write(log_message, '(3A)') "Read '", trim(CS%marbl_settings_file), "' until EOF."
       call MOM_error(NOTE, log_message)
     end if
  end if
  if (is_root_PE()) close(marbl_settings_in)

  ! (3) call marbl%init()
  call MARBL_instances%init(&
                            gcm_num_levels = nz, &
                            gcm_num_PAR_subcols = 1, &
                            gcm_num_elements_surface_flux = 1, & ! FIXME: change to number of grid cells on MPI task
                            gcm_delta_z = GV%sInterface(2:nz+1) - GV%sInterface(1:nz), &
                            gcm_zw = GV%sInterface(2:nz+1), &
                            gcm_zt = GV%sLayer, &
                            lgcm_has_global_ops = .true. &
                           )
  if (marbl_instances%StatusLog%labort_marbl) &
    call marbl_instances%StatusLog%log_error_trace("marbl_instances%init", "configure_MARBL_tracers")
  call print_marbl_log(MARBL_instances%StatusLog)
  call MARBL_instances%StatusLog%erase()

  ! (4) Initialize forcing fields
  !     i. store all surface forcing indices
  CS%u10_sqr_ind = -1
  CS%sss_ind = -1
  CS%sst_ind = -1
  CS%ifrac_ind = -1
  CS%dust_dep_ind = -1
  CS%fe_dep_ind = -1
  CS%nox_flux_ind = -1
  CS%nhy_flux_ind = -1
  CS%atmpress_ind = -1
  CS%xco2_ind = -1
  CS%xco2_alt_ind = -1
  do m=1,size(marbl_instances%surface_flux_forcings)
    select case (trim(marbl_instances%surface_flux_forcings(m)%metadata%varname))
      case('u10_sqr')
        CS%u10_sqr_ind = m
      case('sss')
        CS%sss_ind = m
      case('sst')
        CS%sst_ind = m
      case('Ice Fraction')
        CS%ifrac_ind = m
      case('Dust Flux')
        CS%dust_dep_ind = m
      case('Iron Flux')
        CS%fe_dep_ind = m
      case('NOx Flux')
        CS%nox_flux_ind = m
      case('NHy Flux')
        CS%nhy_flux_ind = m
      case('Atmospheric Pressure')
        CS%atmpress_ind = m
      case('xco2')
        CS%xco2_ind = m
      case('xco2_alt_co2')
        CS%xco2_alt_ind = m
      case DEFAULT
        write(log_message, "(A,1X,A)") trim(marbl_instances%surface_flux_forcings(m)%metadata%varname), &
                                   'is not a valid surface flux forcing field name.'
        call MOM_error(FATAL, log_message)
    end select
  end do

end subroutine configure_MARBL_tracers

!> This subroutine is used to register tracer fields and subroutines
!! to be used with MOM.
function register_MARBL_tracers(HI, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(MARBL_tracers_CS),     pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                 !! structure for the tracer advection and diffusion module.
  type(MOM_restart_CS),       pointer    :: restart_CS !< A pointer to the restart control structure.

! Local variables
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MARBL_tracers" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  character(len=128) :: desc_name ! The variable's descriptor.
  character(len=48)  :: units ! The variable's units.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_MARBL_tracers
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_MARBL_tracers called with an "// &
                             "associated control structure.")
    return
  endif
  allocate(CS)

  call configure_MARBL_tracers(GV, param_file, CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MARBL_TRACERS_IC_FILE", CS%IC_file, &
                 "The file in which the MARBL tracers initial values can be found.", &
                 default="ecosys_jan_IC_omip_MOM_tx0.66v1_c190925.nc")
  if (scan(CS%IC_file,'/') == 0) then
    ! Add the directory if CS%IC_file is not already a complete path.
    call get_param(param_file, mdl, "MARBL_TRACERS_INPUTDIR", inputdir, default="/glade/work/mlevy/cesm_inputdata")
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mdl, "INPUTDIR/MARBL_TRACERS_IC_FILE", CS%IC_file)
  endif

  CS%ntr = size(marbl_instances%tracer_metadata)
  allocate(CS%ind_tr(CS%ntr))
  allocate(CS%tr_desc(CS%ntr))
  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr)) ; CS%tr(:,:,:,:) = 0.0

  do m = 1, CS%ntr
    write(var_name(:),'(A)') trim(marbl_instances%tracer_metadata(m)%short_name)
    write(desc_name(:),'(A)') trim(marbl_instances%tracer_metadata(m)%long_name)
    write(units(:),'(A)') trim(marbl_instances%tracer_metadata(m)%units)
    CS%tr_desc(m) = var_desc(trim(var_name), trim(units), trim(desc_name), caller=mdl)

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=var_name, &
                       caller="register_MARBL_tracers")
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, units = units, &
                         tr_desc=CS%tr_desc(m), registry_diags=.true., &
                         restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_MARBL_tracers")
  enddo

  ! Set up memory for saved state
  call setup_saved_state(marbl_instances%surface_flux_saved_state, HI, GV, restart_CS, CS%tracers_may_reinit, CS%surface_flux_saved_state)
  call setup_saved_state(marbl_instances%interior_tendency_saved_state, HI, GV, restart_CS, CS%tracers_may_reinit, CS%interior_tendency_saved_state)

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_MARBL_tracers = .true.
end function register_MARBL_tracers

!> This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
!! and it sets up the tracer output.
subroutine initialize_MARBL_tracers(restart, day, G, GV, US, h, diag, OBC, CS, sponge_CSp)
  logical,                            intent(in) :: restart !< .true. if the fields have already been
                                                            !! read from a restart file.
  type(time_type), target,            intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: h !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl), target,            intent(in) :: diag !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(MARBL_tracers_CS),                pointer    :: CS   !< The control structure returned by a previous
                                                         !! call to register_MARBL_tracers.
  type(sponge_CS),                    pointer    :: sponge_CSp    !< A pointer to the control structure
                                                                  !! for the sponges, if they are in use.

! Local variables
  character(len=48) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                ! years m3 s-1 or years kg s-1.
  logical :: OK
  integer :: i, j, k, m, diag_size

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  CS%diag => diag

  ! Allocate memory for surface tracer fluxes
  allocate(CS%STF(SZI_(G), SZJ_(G), CS%ntr))
  CS%STF(:,:,:) = 0.

  ! Register diagnostics (surface flux first, then interior tendency)
  call register_MARBL_diags(marbl_instances%surface_flux_diags, diag, day, G, CS%surface_flux_diags)
  call register_MARBL_diags(marbl_instances%interior_tendency_diags, diag, day, G, CS%interior_tendency_diags)

  ! Initialize tracers from file (unless they were initialized by restart file)
  if (.not.restart) then
    do m= 1, CS%ntr
      write(name(:),'(A)') trim(marbl_instances%tracer_metadata(m)%short_name)
      OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, name, G, US, -1e34)
      if (.not.OK) call MOM_error(FATAL,"initialize_MARBL_tracers: "//&
                                  "Unable to read "//trim(name)//" from "//&
                                  trim(CS%IC_file)//".")
    end do
  end if

end subroutine initialize_MARBL_tracers

!> This subroutine is used to register tracer fields and subroutines
!! to be used with MOM.
subroutine register_MARBL_diags(MARBL_diags, diag, day, G, id_diags)

  type(marbl_diagnostics_type), intent(in)    :: MARBL_diags !< MARBL diagnostics from marbl_instances
  type(time_type), target,      intent(in)    :: day  !< Time of the start of the run.
  type(diag_ctrl), target,      intent(in)    :: diag !< Structure used to regulate diagnostic output.
  !integer, allocatable,         intent(inout) :: id_diags(:) !< allocatable array storing diagnostic index number
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(temp_marbl_diag), allocatable, intent(inout) :: id_diags(:) !< allocatable array storing diagnostic index
                                                                   !! number and buffer space for collecting diags
                                                                   !! from all columns

  integer :: m, diag_size

  diag_size = size(MARBL_diags%diags)
  allocate(id_diags(diag_size))
  do m = 1, diag_size
    id_diags(m)%id = -1
    if (trim(MARBL_diags%diags(m)%vertical_grid) .eq. "none") then ! 2D field
      id_diags(m)%id = register_diag_field("ocean_model", &
        trim(MARBL_diags%diags(m)%short_name), &
        diag%axesT1, & ! T => tracer grid? 1 => no vertical grid
        day, &
        trim(MARBL_diags%diags(m)%long_name), &
        trim(MARBL_diags%diags(m)%units))
      allocate(id_diags(m)%field_2d(SZI_(G),SZJ_(G)))
      id_diags(m)%field_2d(:,:) = 0.
    else ! 3D field
      id_diags(m)%id = register_diag_field("ocean_model", &
        trim(MARBL_diags%diags(m)%short_name), &
        diag%axesTL, & ! T=> tracer grid? L => layer center
        day, &
        trim(MARBL_diags%diags(m)%long_name), &
        trim(MARBL_diags%diags(m)%units))
      allocate(id_diags(m)%field_3d(SZI_(G),SZJ_(G), SZK_(G)))
      id_diags(m)%field_3d(:,:,:) = 0.
    end if
  end do

end subroutine register_MARBL_diags

!> This subroutine allocates memory for saved state fields and registers them in the restart files
subroutine setup_saved_state(MARBL_saved_state, HI, GV, restart_CS, tracers_may_reinit, local_saved_state)

  type(marbl_saved_state_type),                  intent(in)    :: MARBL_saved_state !< MARBL saved state from marbl_instances
  type(hor_index_type),                          intent(in)    :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),                       intent(in)    :: GV   !< The ocean's vertical grid structure
  type(MOM_restart_CS), pointer,                 intent(in)    :: restart_CS !< control structure to add saved state to restarts
  logical,                                       intent(in)    :: tracers_may_reinit  !< used to determine mandatory flag in restart
  type(saved_state_for_MARBL_type), allocatable, intent(inout) :: local_saved_state(:) !< allocatable array storing saved state locally

  integer :: num_fields, m
  character(len=200) :: log_message, varname

  num_fields = MARBL_saved_state%saved_state_cnt
  allocate(local_saved_state(num_fields))

  do m=1,num_fields
    write(varname, "(2A)") "MARBL_", trim(MARBL_saved_state%state(m)%short_name)
    select case (MARBL_saved_state%state(m)%rank)
      case (2)
        allocate(local_saved_state(m)%field_2d(SZI_(HI),SZJ_(HI)))
        local_saved_state(m)%field_2d(:,:) = 0.
        call register_restart_field(local_saved_state(m)%field_2d, varname, .not.tracers_may_reinit, restart_CS)
      case (3)
        if (trim(MARBL_saved_state%state(m)%vertical_grid).eq."layer_avg") then
          allocate(local_saved_state(m)%field_3d(SZI_(HI),SZJ_(HI), SZK_(GV)))
          local_saved_state(m)%field_3d(:,:,:) = 0.
          call register_restart_field(local_saved_state(m)%field_3d, varname, .not.tracers_may_reinit, restart_CS)
        else
          write(log_message, "(3A, I0, A)") "'", trim(MARBL_saved_state%state(m)%vertical_grid), &
                "' is an invalid vertical grid for saved state (ind = ", m, ")"
          call MOM_error(FATAL, log_message)
        end if
      case DEFAULT
        write(log_message, "(I0, A, I0, A)") MARBL_saved_state%state(m)%rank, &
              " is an invalid rank for saved state (ind = ", m, ")"
        call MOM_error(FATAL, log_message)
    end select
    local_saved_state(m)%short_name = trim(MARBL_saved_state%state(m)%short_name)
    write(local_saved_state(m)%file_varname, "(2A)") "MARBL_", trim(local_saved_state(m)%short_name)
    local_saved_state(m)%units = trim(MARBL_saved_state%state(m)%units)
  end do

end subroutine setup_saved_state

!> This subroutine applies diapycnal diffusion and any other column
!! tracer physics or chemistry to the tracers from this file.
subroutine MARBL_tracers_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, tv, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(MARBL_tracers_CS),     pointer :: CS   !< The control structure returned by a previous
                                              !! call to register_MARBL_tracers.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [m]

! Local variables
  character(len=256) :: log_message
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified
  real :: sfc_val  ! The surface value for the tracers.
  real :: Isecs_per_year  ! The number of seconds in a year.
  real :: year            ! The time in years.
  integer :: secs, days   ! Integer components of the time type.
  real, dimension(0:GV%ke) :: zi  ! z-coordinate interface depth
  real, dimension(GV%ke) :: zc, dz  ! z-coordinate layer center depth and cell thickness
  integer :: i, j, k, is, ie, js, je, nz, m
  real :: kg_m2_s_conversion, ndep_conversion

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  kg_m2_s_conversion = US%RZ_T_to_kg_m2s
  ndep_conversion = (US%m_to_L)**2 * US%s_to_T

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  ! (1) Compute surface fluxes
  ! FIXME: MARBL can handle computing surface fluxes for all columns simultaneously
  !        I was just thinking going column-by-column at first might be easier
  do i=is,ie
    do j=js,je
      ! i. only want ocean points in this loop
      if (G%mask2dT(i,j) == 0) cycle

      ! ii. Load proper column data
      !     * surface flux forcings
      !       These fields are getting the correct data
      if (CS%sss_ind > 0) marbl_instances%surface_flux_forcings(CS%sss_ind)%field_0d(1) = tv%S(i,j,1)
      if (CS%sst_ind > 0) marbl_instances%surface_flux_forcings(CS%sst_ind)%field_0d(1) = tv%T(i,j,1)
      if (CS%ifrac_ind > 0) marbl_instances%surface_flux_forcings(CS%ifrac_ind)%field_0d(1) = fluxes%ice_fraction(i,j)
      ! MARBL wants u10_sqr in (cm/s)^2
      if (CS%u10_sqr_ind > 0) marbl_instances%surface_flux_forcings(CS%u10_sqr_ind)%field_0d(1) = fluxes%u10_sqr(i,j) * (US%L_t_to_m_s * cm_per_m)**2
      ! mct_driver/ocn_cap_methods:93 -- ice_ocean_boundary%p(i,j) comes from coupler
      ! We may need a new ice_ocean_boundary%p_atm because %p includes ice in GFDL driver
      if (CS%atmpress_ind > 0) marbl_instances%surface_flux_forcings(CS%atmpress_ind)%field_0d(1) = fluxes%p_surf_full(i,j) * &
                                                                                                    ((US%R_to_kg_m3 * US%L_T_to_m_s**2) * atm_per_Pa)

      !       These are okay, but need option to come in from coupler
      if (CS%xco2_ind > 0) marbl_instances%surface_flux_forcings(CS%xco2_ind)%field_0d(1) = CS%atm_co2_const
      if (CS%xco2_alt_ind > 0) marbl_instances%surface_flux_forcings(CS%xco2_alt_ind)%field_0d(1) = CS%atm_alt_co2_const

      !       These are okay, but need option to read in from file
      if (CS%dust_dep_ind > 0) marbl_instances%surface_flux_forcings(CS%dust_dep_ind)%field_0d(1) = fluxes%dust_flux(i,j) * (kg_m2_s_conversion * g_per_kg * m_per_cm**2)
      if (CS%fe_dep_ind > 0) marbl_instances%surface_flux_forcings(CS%fe_dep_ind)%field_0d(1) = fluxes%iron_flux(i,j) * (kg_m2_s_conversion * g_per_kg * m_per_cm**2)

      !       These are read from /glade/work/mlevy/cesm_inputdata/ndep_ocn_1850_w_nhx_emis_MOM_tx0.66v1_c200827.nc
      if (CS%nox_flux_ind > 0) marbl_instances%surface_flux_forcings(CS%nox_flux_ind)%field_0d(1) = fluxes%noy_dep(i,j) * (ndep_conversion * CS%ndep_scale_factor)
      if (CS%nhy_flux_ind > 0) marbl_instances%surface_flux_forcings(CS%nhy_flux_ind)%field_0d(1) = fluxes%nhx_dep(i,j) * (ndep_conversion * CS%ndep_scale_factor)

      !     * tracers at surface
      do m=1,CS%ntr
        marbl_instances%tracers_at_surface(1,m) = CS%tr(i,j,1,m)
      end do

      !     * surface flux saved state
      do m=1,size(marbl_instances%surface_flux_saved_state%state)
      !       (currently only 2D fields are saved from surface_flux_compute())
        marbl_instances%surface_flux_saved_state%state(m)%field_2d(1) = CS%surface_flux_saved_state(m)%field_2d(i,j)
      end do

      ! iii. Compute surface fluxes in MARBL
      call marbl_instances%surface_flux_compute()
      if (marbl_instances%StatusLog%labort_marbl) then
        call marbl_instances%StatusLog%log_error_trace("marbl_instances%surface_flux_compute()", &
                                                       "MARBL_tracers_column_physics")
        call print_marbl_log(MARBL_instances%StatusLog)
      end if

      ! iv. Copy output that MOM6 needs to hold on to
      !     * saved state
      do m=1,size(marbl_instances%surface_flux_saved_state%state)
        CS%surface_flux_saved_state(m)%field_2d(i,j) = marbl_instances%surface_flux_saved_state%state(m)%field_2d(1)
      end do

      !     * diagnostics
      do m=1,size(marbl_instances%surface_flux_diags%diags)
        ! All diags are 2D coming from surface
        CS%surface_flux_diags(m)%field_2d(i,j) = real(marbl_instances%surface_flux_diags%diags(m)%field_2d(1))
      end do

      !     * Surface tracer flux
      CS%STF(i,j,:) = marbl_instances%surface_fluxes(1,:) * m_per_cm
    end do
  end do

  ! (2) Post surface flux diagnostics (currently all 2D)
  do m=1,size(CS%surface_flux_diags)
    if (CS%surface_flux_diags(m)%id > 0) &
      call post_data(CS%surface_flux_diags(m)%id, CS%surface_flux_diags(m)%field_2d(:,:), CS%diag)
  end do

  ! (3) Apply surface fluxes via vertical diffusion
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,CS%ntr
      do k=1,nz ;do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m) , dt, fluxes, h_work, &
          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV, sfc_flux=CS%STF(:,:,m))
    enddo
  else
    do m=1,CS%ntr
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV, sfc_flux=CS%STF(:,:,m))
    enddo
  endif

  ! (4) Compute interior tendencies
  do i=is,ie
    do j=js,je
      ! i. only want ocean points in this loop
      if (G%mask2dT(i,j) == 0) cycle

      ! ii. Set up vertical domain
      marbl_instances%domain%kmt = GV%ke
      ! Calculate depth of interface by building up thicknesses from the bottom (top interface is always 0)
      ! MARBL wants this to be positive-down
      zi(GV%ke) = G%bathyT(i,j)
      do k = GV%ke, 1, -1
        ! TODO: h_new or h_old?
        dz(k) = h_new(i,j,k)*GV%H_to_Z ! cell thickness
        zc(k) = zi(k) - 0.5 * dz(k)
        zi(k-1) = zi(k) - dz(k)
      enddo
      ! mks -> cgs
      ! zw(1:nz) is bottom cell depth so no element of zw = 0, it is assumed to be top layer depth
      marbl_instances%domain%zw(:) = zi(1:GV%ke) * m_per_cm
      marbl_instances%domain%zt(:) = zc(:) * m_per_cm
      marbl_instances%domain%delta_z(:) = dz(:) * m_per_cm

      ! iii. Load proper column data
      !      * Forcing Fields
      !      * Column Tracers
      !        NOTE: POP averages previous two timesteps, should we do that too?
      do m=1,CS%ntr
        marbl_instances%tracers(m, :) = CS%tr(i,j,:,m)
      end do

      !     * interior tendency saved state
      !       (currently only 3D fields are saved from interior_tendency_compute())
      do m=1,size(marbl_instances%interior_tendency_saved_state%state)
        marbl_instances%interior_tendency_saved_state%state(m)%field_3d(:,1) = CS%interior_tendency_saved_state(m)%field_3d(i,j,:)
      end do

      ! iv. Compute interior tendencies in MARBL
      ! call marbl_instances%interior_tendency_compute()
      ! if (marbl_instances%StatusLog%labort_marbl) then
      !   call marbl_instances%StatusLog%log_error_trace("marbl_instances%interior_tendency_compute()", &
      !                                                  "MARBL_tracers_column_physics")
      !   call print_marbl_log(MARBL_instances%StatusLog)
      ! end if

      ! v. Apply tendencies immediately
      do k=1,GV%ke
        CS%tr(i,j,k,:) = CS%tr(i,j,k,:) + G%mask2dT(i,j)*dt*marbl_instances%interior_tendencies(:, k)
      end do

      ! vi. Copy output that MOM6 needs to hold on to
      !     * saved state
      do m=1,size(marbl_instances%interior_tendency_saved_state%state)
        CS%interior_tendency_saved_state(m)%field_3d(i,j,:) = marbl_instances%interior_tendency_saved_state%state(m)%field_3d(:,1)
      end do

      !     * diagnostics
      do m=1,size(marbl_instances%interior_tendency_diags%diags)
        if (allocated(CS%interior_tendency_diags(m)%field_2d)) then
          CS%interior_tendency_diags(m)%field_2d(i,j) = real(marbl_instances%interior_tendency_diags%diags(m)%field_2d(1))
        else
          CS%interior_tendency_diags(m)%field_3d(i,j,:) = real(marbl_instances%interior_tendency_diags%diags(m)%field_3d(:,1))
        end if
      end do

    end do
  end do

  ! (5) Post diagnostics from our buffer
  !     i. Surface flux diagnostics (currently all 2D)
  !     ii. Interior tendency diagnostics (mix of 2D and 3D)
  do m=1,size(CS%interior_tendency_diags)
    if (CS%interior_tendency_diags(m)%id > 0) then
      if (allocated(CS%interior_tendency_diags(m)%field_2d)) then
        call post_data(CS%interior_tendency_diags(m)%id, CS%interior_tendency_diags(m)%field_2d(:,:), CS%diag)
      else
        call post_data(CS%interior_tendency_diags(m)%id, CS%interior_tendency_diags(m)%field_3d(:,:,:), CS%diag)
      end if
    end if
  end do

end subroutine MARBL_tracers_column_physics

!> This function calculates the mass-weighted integral of all tracer stocks,
!! returning the number of stocks it has calculated.  If the stock_index
!! is present, only the stock corresponding to that coded index is returned.
function MARBL_tracer_stock(h, stocks, G, GV, CS, names, units, stock_index)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(:),                 intent(out)   :: stocks !< the mass-weighted integrated amount of
                                                            !! each tracer, in kg times concentration units [kg conc].
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  type(MARBL_tracers_CS),                pointer       :: CS   !< The control structure returned by a
                                                            !! previous call to register_MARBL_tracers.
  character(len=*), dimension(:),     intent(out)   :: names !< the names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units !< the units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< the coded index of a specific stock
                                                                   !! being sought.
  integer                                           :: MARBL_tracer_stock   !< Return value: the number of stocks
                                                                   !! calculated here.

! Local variables
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  MARBL_tracer_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="MARBL_tracer_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      stocks(m) = stocks(m) + CS%tr(i,j,k,m) * &
                             (G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k))
    enddo ; enddo ; enddo
    stocks(m) = GV%H_to_kg_m2 * stocks(m)
  enddo
  MARBL_tracer_stock = CS%ntr

end function MARBL_tracer_stock

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
!! This particular tracer package does not report anything back to the coupler.
subroutine MARBL_tracers_surface_state(state, h, G, CS)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure.
  type(surface),          intent(inout) :: state !< A structure containing fields that
                                              !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                          intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(MARBL_tracers_CS),    pointer       :: CS !< The control structure returned by a previous
                                              !! call to register_MARBL_tracers.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

  integer :: m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,CS%ntr
      !   This call loads the surface values into the appropriate array in the
      ! coupler-type structure.
      call coupler_type_set_data(CS%tr(:,:,1,m), CS%ind_tr(m), ind_csurf, &
                   state%tr_fields, idim=(/isd, is, ie, ied/), &
                   jdim=(/jsd, js, je, jed/) )
    enddo
  endif

end subroutine MARBL_tracers_surface_state

!> Clean up any allocated memory after the run.
subroutine MARBL_tracers_end(CS)
  type(MARBL_tracers_CS), pointer :: CS !< The control structure returned by a previous
                                        !! call to register_MARBL_tracers.
  integer :: m

  call print_marbl_log(MARBL_instances%StatusLog)
  call MARBL_instances%StatusLog%erase()
  call marbl_instances%shutdown()
  ! TODO: print MARBL timers to stdout as well

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine MARBL_tracers_end

!> This subroutine writes the contents of the MARBL log using MOM_error(NOTE, ...).
!! TODO: some log messages come from a specific grid point, and this routine
!!       needs to include the location in the preamble
subroutine print_marbl_log(log_to_print)

  use marbl_logging, only : marbl_status_log_entry_type
  use marbl_logging, only : marbl_log_type
  use MOM_coms,      only : PE_here

  class(marbl_log_type), intent(in) :: log_to_print

  character(len=*), parameter :: subname = 'MARBL_tracers:print_marbl_log'
  character(len=256)          :: message_prefix, message_location, log_message
  type(marbl_status_log_entry_type), pointer :: tmp
  integer :: msg_lev

  write(message_prefix, "(A,I0,A)") '(Task ', PE_here(), ')'

  tmp => log_to_print%FullLog
  do while (associated(tmp))
    ! 1) Do I need to write this message? Yes, if all tasks should write this
    !    or if I am master_task
    if ((.not. tmp%lonly_master_writes) .or. is_root_PE()) then
      ! master task does not need prefix
      if (.not. is_root_PE()) then
        write(log_message, "(A,1X,A)") trim(message_prefix), trim(tmp%LogMessage)
        msg_lev = WARNING
      else
        write(log_message, "(A)") trim(tmp%LogMessage)
        msg_lev = NOTE
      end if     ! print message prefix?
      call MOM_error(msg_lev, log_message, all_print=.true.)
    end if       ! write the message?
    tmp => tmp%next
  end do

  if (log_to_print%labort_marbl) then
    call MOM_error(WARNING, 'ERROR reported from MARBL library', all_print=.true.)
    call MOM_error(FATAL, 'Stopping in ' // subname)
  end if

end subroutine print_marbl_log

#endif /* _USE_MARBL_TRACERS */
!> \namespace MARBL_tracers
!!
!!    This file contains an example of the code that is needed to set
!!  up and use a set (in this case two) of dynamically passive tracers
!!  for diagnostic purposes.  The tracers here are dye tracers which
!!  are set to 1 within the geographical region specified. The depth
!!  which a tracer is set is determined by calculating the depth from
!!  the seafloor upwards through the column.

end module MARBL_tracers
