!> A tracer package for tracers computed in the MARBL library
!!
!! Currently configured for use with marbl0.36.0
!! https://github.com/marbl-ecosys/MARBL/releases/tag/marbl0.36.0
!! (clone entire repo into pkg/MARBL)
module MARBL_tracers

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,            only : root_PE, broadcast
use MOM_diag_mediator,   only : diag_ctrl
use MOM_diag_vkernels,   only : reintegrate_column
use MOM_error_handler,   only : is_root_PE, MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_grid,            only : ocean_grid_type
use MOM_CVMix_KPP,       only : KPP_NonLocalTransport, KPP_CS
use MOM_hor_index,       only : hor_index_type
use MOM_io,              only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : query_initialized, MOM_restart_CS, register_restart_field
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_registry, only : register_tracer
use MOM_tracer_types,    only : tracer_type, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_initialization_from_Z, only : MOM_initialize_tracer_from_Z
use MOM_tracer_Z_init,   only : read_Z_edges
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
public MARBL_tracer_stock, MARBL_tracers_get_output_for_GCM, MARBL_tracers_end

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

!> MOM6 needs to know the index of some MARBL tracers to properly apply river fluxes
type :: tracer_ind_type
  integer :: no3_ind  !< NO3 index
  integer :: po4_ind  !< PO4 index
  integer :: don_ind  !< DON index
  integer :: donr_ind  !< DONr index
  integer :: dop_ind  !< DOP index
  integer :: dopr_ind  !< DOPr index
  integer :: sio3_ind  !< SiO3 index
  integer :: fe_ind  !< Fe index
  integer :: doc_ind  !< DOC index
  integer :: docr_ind  !< DOCr index
  integer :: alk_ind  !< ALK index
  integer :: alk_alt_co2_ind  !< ALK_ALT_CO2 index
  integer :: dic_ind  !< DIC index
  integer :: dic_alt_co2_ind  !< DIC_ALT_CO2 index
end type tracer_ind_type

!> MOM needs to store some information about saved_state; besides providing these
!! fields to MARBL, they are also written to restart files
type :: saved_state_for_MARBL_type
  character(len=200) :: short_name !< name of variable being saved
  character(len=200) :: file_varname !< name of variable in restart file
  character(len=200) :: units !< variable units
  real, pointer :: field_2d(:,:) !< memory for 2D field
  real, pointer :: field_3d(:,:,:) !< memory for 3D field
end type saved_state_for_MARBL_type

!> All calls to MARBL are done via the interface class
type(MARBL_interface_class) :: MARBL_instances

!> Pointer to tracer concentration and to tracer_type in tracer registry
type, private :: MARBL_tracer_data
  real, pointer              :: tr(:,:,:) => NULL() !< The array of tracers used in this subroutine, in g m-3?
  type(tracer_type), pointer :: tr_ptr    => NULL() !< pointer to tracer inside Tr_reg
end type MARBL_tracer_data

!> The control structure for the MARBL tracer package
type, public :: MARBL_tracers_CS ; private
  integer :: ntr    !< The number of tracers that are actually used.
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the coupler.
  logical :: use_ice_category_fields    !< Forcing will include multiple ice categories for ice_frac and shortwave
  logical :: request_Chl_from_MARBL     !< MARBL can provide Chl to use in set_pen_shortwave()
  integer :: ice_ncat                   !< Number of ice categories when use_ice_category_fields = True
  character(len=200) :: IC_file !< The file in which the age-tracer initial values cam be found.
  logical :: ongrid                     !< True if IC_file is already interpolated to MOM grid
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  type(MARBL_tracer_data), dimension(:), allocatable :: tracer_data  !< type containing tracer data and pointer
                                                                     !! into tracer registry

  integer, allocatable, dimension(:) :: ind_tr !< Indices returned by aof_set_coupler_flux if it is used and the
                                               !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(vardesc), allocatable :: tr_desc(:) !< Descriptions and metadata for the tracers
  logical :: tracers_may_reinit = .false. !< If true the tracers may be initialized if not found in a restart file

  character(len=200) :: fesedflux_file  !< name of [netCDF] file containing iron sediment flux
  character(len=200) :: feventflux_file  !< name of [netCDF] file containing iron vent flux
  character(len=35) :: marbl_settings_file  !< name of [text] file containing MARBL settings

  real :: bot_flux_mix_thickness !< for bottom flux -> tendency conversion, assume uniform mixing over
                                 !! bottom layer of prescribed thickness
  real :: bfmt_r                 !< Reciprocal of above

  type(temp_MARBL_diag), allocatable :: surface_flux_diags(:)  !< collect surface flux diagnostics from all columns
                                                               !! before posting
  type(temp_MARBL_diag), allocatable :: interior_tendency_diags(:)  !< collect tendency diagnostics from all columns
                                                                    !! before posting
  type(saved_state_for_MARBL_type), allocatable :: surface_flux_saved_state(:)  !< surface_flux saved state
  type(saved_state_for_MARBL_type), allocatable :: interior_tendency_saved_state(:)  !< interior_tendency saved state

  ! TODO: If we can post data column by column, all we need are integer arrays for ids
  ! integer, allocatable :: id_surface_flux_diags(:)  !< array of indices for surface_flux diagnostics
  ! integer, allocatable :: id_interior_tendency_diags(:)  !< array of indices for interior_tendency diagnostics

  type(tracer_ind_type) :: tracer_inds  !< Indices to tracers that will have river fluxes added to STF

  !> Need to store global output from both marbl_instance%surface_flux_compute() and
  !! marbl_instance%interior_tendency_compute(). For the former, just need id to register
  !! because we already copy data into CS%STF; latter requires copying data and indices
  !! so currently using temp_MARBL_diag for that.
  integer, allocatable :: id_surface_flux_out(:)  !< register_diag indices for surface_flux output
  type(temp_MARBL_diag), allocatable :: interior_tendency_out(:)  !< collect interior tendencies for diagnostic output
  type(temp_MARBL_diag), allocatable :: interior_tendency_out_zint(:)  !< vertical integral of interior tendencies
                                                                       !! (full column)
  type(temp_MARBL_diag), allocatable :: interior_tendency_out_zint_100m(:)  !< vertical integral of interior tendencies
                                                                            !! (top 100m)
  integer :: bot_flux_to_tend_id  !< register_diag index for BOT_FLUX_TO_TEND
  integer, allocatable :: fracr_cat_id(:) !< register_diag index for per-category ice fraction
  integer, allocatable :: qsw_cat_id(:)   !< register_diag index for per-category shortwave

  real, allocatable :: STF(:,:,:) !< surface fluxes returned from MARBL to use in tracer_vertdiff (dims: i, j, tracer)
                                  !! [conc m/s]
  real, allocatable :: RIV_FLUXES(:,:,:) !< time-integrate river flux forcing for applyTracerBoundaryFluxesInOut
                                         !! (dims: i, j, tracer) [conc m]
  real, allocatable :: SFO(:,:,:)  !< surface flux output returned from MARBL for use in GCM
                                   !! e.g. CO2 flux to pass to atmosphere (dims: i, j, num_sfo)

  integer :: u10_sqr_ind  !< index of MARBL forcing field array to copy 10-m wind (squared) into
  integer :: sss_ind  !< index of MARBL forcing field array to copy sea surface salinity into
  integer :: sst_ind  !< index of MARBL forcing field array to copy sea surface temperature into
  integer :: ifrac_ind  !< index of MARBL forcing field array to copy ice fraction into
  integer :: dust_dep_ind  !< index of MARBL forcing field array to copy dust flux into
  integer :: fe_dep_ind  !< index of MARBL forcing field array to copy iron flux into
  integer :: nox_flux_ind  !< index of MARBL forcing field array to copy NOx flux into
  integer :: nhy_flux_ind  !< index of MARBL forcing field array to copy NHy flux into
  integer :: atmpress_ind  !< index of MARBL forcing field array to copy atmospheric pressure into
  integer :: xco2_ind  !< index of MARBL forcing field array to copy CO2 flux into
  integer :: xco2_alt_ind  !< index of MARBL forcing field array to copy CO2 flux (alternate CO2) into

  !> Indices for forcing fields required to compute interior tendencies
  integer :: dustflux_ind  !< index of MARBL forcing field array to copy dust flux into
  integer :: PAR_col_frac_ind  !< index of MARBL forcing field array to copy PAR column fraction into
  integer :: surf_shortwave_ind  !< index of MARBL forcing field array to copy surface shortwave into
  integer :: potemp_ind  !< index of MARBL forcing field array to copy potential temperature into
  integer :: salinity_ind  !< index of MARBL forcing field array to copy salinity into
  integer :: pressure_ind  !< index of MARBL forcing field array to copy pressure into
  integer :: fesedflux_ind  !< index of MARBL forcing field array to copy iron sediment flux into
  integer :: o2_scalef_ind  !< index of MARBL forcing field array to copy O2 scale length into
  integer :: remin_scalef_ind  !< index of MARBL forcing field array to copy remin scale length into

  !> Number of surface flux outputs as well as specific indices for each one
  integer :: sfo_cnt
  integer :: flux_co2_ind

  ! TODO: create generic 3D forcing input type to read z coordinate + values
  real    :: fesedflux_scale_factor !< scale factor for iron sediment flux
  integer :: fesedflux_nz  !< number of levels in iron sediment flux file
  real, allocatable, dimension(:,:,:) :: fesedflux_in  !< Field to read iron sediment flux into
  real, allocatable, dimension(:,:,:) :: feventflux_in  !< Field to read iron vent flux into
  real, allocatable, dimension(:) :: &
    fesedflux_z_edges  !< The depths of the cell interfaces in the input data [Z ~> m]
  ! TODO: this thickness does not need to be 3D, but that's a problem for future Mike
  real, allocatable, dimension(:,:,:) :: &
    fesedflux_dz  !< The thickness of the cell layers in the input data [Z ~> m]
end type MARBL_tracers_CS

! Module parameters
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
  logical :: chl_from_file
  nz = GV%ke
  marbl_settings_in = 615

  ! (1) Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MARBL_SETTINGS_FILE", CS%marbl_settings_file, &
                 "The name of a file from which to read the run-time "//&
                 "settings for MARBL.", default="marbl_in")
  call get_param(param_file, mdl, "BOT_FLUX_MIX_THICKNESS", CS%bot_flux_mix_thickness, &
                 "Bottom fluxes are uniformly mixed over layer of this thickness", &
                 default=1., units="m")
  call get_param(param_file, mdl, "CHL_FROM_FILE", chl_from_file, &
                 "If true, chl_a is read from a file.", default=.true.)
  CS%request_Chl_from_MARBL = (.not. chl_from_file)
  call get_param(param_file, mdl, "USE_ICE_CATEGORIES", CS%use_ice_category_fields, &
       "If true, allocate memory for shortwave and ice fraction split by ice thickness category.", &
       default=.false.)
  call get_param(param_file, mdl, "ICE_NCAT", CS%ice_ncat, &
       "Number of ice thickness categories in shortwave and ice fraction forcings.", &
       default=0)
  CS%bfmt_r = 1. / CS%bot_flux_mix_thickness

  if (CS%use_ice_category_fields .and. (CS%ice_ncat == 0)) &
    call MOM_error(FATAL, "Can not configure MARBL to use multiple ice categories without ice_ncat present")

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
    call MARBL_instances%put_setting(marbl_in_line(1))
  end do
  ! iv. (TEMPORARY) don't set tracer restoring
  !     TODO: fix this
  call MARBL_instances%put_setting("tracer_restore_vars(1) = ''")
  call MARBL_instances%put_setting("tracer_restore_vars(2) = ''")
  call MARBL_instances%put_setting("tracer_restore_vars(3) = ''")
  call MARBL_instances%put_setting("tracer_restore_vars(4) = ''")
  call MARBL_instances%put_setting("tracer_restore_vars(5) = ''")

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
  ! TODO: the units in gcm_delta_z, gcm_zw, and gcm_zt are wrong, but we want to strip these values
  !       out of init anyway because MOM updates them every time step / every column
  call MARBL_instances%init(&
                            gcm_num_levels = nz, &
                            gcm_num_PAR_subcols = CS%ice_ncat + 1, &
                            gcm_num_elements_surface_flux = 1, & ! FIXME: change to number of grid cells on MPI task
                            gcm_delta_z = GV%sInterface(2:nz+1) - GV%sInterface(1:nz), &
                            gcm_zw = GV%sInterface(2:nz+1), &
                            gcm_zt = GV%sLayer, &
                            lgcm_has_global_ops = .true. &
                           )
  if (MARBL_instances%StatusLog%labort_marbl) &
    call MARBL_instances%StatusLog%log_error_trace("MARBL_instances%init", "configure_MARBL_tracers")
  call print_marbl_log(MARBL_instances%StatusLog)
  call MARBL_instances%StatusLog%erase()

  ! (4) Request fields needed by MOM6
  CS%sfo_cnt = 0

  ! CO2 Flux to the atmosphere
  CS%sfo_cnt = CS%sfo_cnt +1
  call MARBL_instances%surface_flux_output%add_output(num_elements=1, &
                                                      field_name="flux_co2", &
                                                      output_id=CS%flux_co2_ind, &
                                                      marbl_status_log=MARBL_instances%StatusLog)

  ! (5) Initialize forcing fields
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
  do m=1,size(MARBL_instances%surface_flux_forcings)
    select case (trim(MARBL_instances%surface_flux_forcings(m)%metadata%varname))
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
        write(log_message, "(A,1X,A)") trim(MARBL_instances%surface_flux_forcings(m)%metadata%varname), &
                                   'is not a valid surface flux forcing field name.'
        call MOM_error(FATAL, log_message)
    end select
  end do

  !     ii. store all interior forcing indices
  CS%dustflux_ind = -1
  CS%PAR_col_frac_ind = -1
  CS%surf_shortwave_ind = -1
  CS%potemp_ind = -1
  CS%salinity_ind = -1
  CS%pressure_ind = -1
  CS%fesedflux_ind = -1
  CS%o2_scalef_ind = -1
  CS%remin_scalef_ind = -1
  do m=1,size(MARBL_instances%interior_tendency_forcings)
    ! i. Check to see if this is a tracer restoring field or timescale
    ! ii. If not, should be one of the following:
    select case (trim(MARBL_instances%interior_tendency_forcings(m)%metadata%varname))
      case('Dust Flux')
        CS%dustflux_ind = m
      case('PAR Column Fraction')
        CS%PAR_col_frac_ind = m
      case('Surface Shortwave')
        CS%surf_shortwave_ind = m
      case('Potential Temperature')
        CS%potemp_ind = m
      case('Salinity')
        CS%salinity_ind = m
      case('Pressure')
        CS%pressure_ind = m
      case('Iron Sediment Flux')
        CS%fesedflux_ind = m
      case('O2 Consumption Scale Factor')
        CS%o2_scalef_ind = m
      case('Particulate Remin Scale Factor')
        CS%remin_scalef_ind = m
      case DEFAULT
        write(log_message, "(A,1X,A)") trim(MARBL_instances%interior_tendency_forcings(m)%metadata%varname), &
                                   'is not a valid interior tendency forcing field name.'
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
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct

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
  ! ** Input directory
  ! TODO: just use DIN_LOC_ROOT
  call get_param(param_file, mdl, "CESM_INPUTDIR", inputdir, default="/glade/work/mlevy/cesm_inputdata")
  ! ** Tracer initial conditions
  call get_param(param_file, mdl, "MARBL_TRACERS_IC_FILE", CS%IC_file, &
                 "The file in which the MARBL tracers initial values can be found.", &
                 default="ecosys_jan_IC_omip_latlon_1x1_180W_c230331.nc")
  if (scan(CS%IC_file,'/') == 0) then
    ! Add the directory if CS%IC_file is not already a complete path.
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mdl, "INPUTDIR/MARBL_TRACERS_IC_FILE", CS%IC_file)
  endif
  call get_param(param_file, mdl, "MARBL_TRACERS_INIT_VERTICAL_REMAP_ONLY", CS%ongrid, &
                 "If true, initial conditions are on the model horizontal grid. " //&
                 "Extrapolation over missing ocean values is done using an ICE-9 "//&
                 "procedure with vertical ALE remapping .", &
                 default=.false.)
  ! ** FESEDFLUX
  call get_param(param_file, mdl, "MARBL_FESEDFLUX_FILE", CS%fesedflux_file, &
                 "The file in which the iron sediment flux forcing field can be found.", &
                 default="fesedflux_total_reduce_oxic_tx0.66v1.c211109.nc")
  if (scan(CS%fesedflux_file,'/') == 0) then
    ! Add the directory if CS%fesedflux_file is not already a complete path.
    CS%fesedflux_file = trim(slasher(inputdir))//trim(CS%fesedflux_file)
    call log_param(param_file, mdl, "INPUTDIR/MARBL_TRACERS_FESEDFLUX_FILE", CS%fesedflux_file)
  endif
  ! ** FEVENTFLUX
  call get_param(param_file, mdl, "MARBL_FEVENTFLUX_FILE", CS%feventflux_file, &
                 "The file in which the iron vent flux forcing field can be found.", &
                 default="feventflux_5gmol_tx0.66v1.c211109.nc")
  if (scan(CS%feventflux_file,'/') == 0) then
    ! Add the directory if CS%feventflux_file is not already a complete path.
    CS%feventflux_file = trim(slasher(inputdir))//trim(CS%feventflux_file)
    call log_param(param_file, mdl, "INPUTDIR/MARBL_TRACERS_FEVENTFLUX_FILE", CS%feventflux_file)
  endif
  ! ** Scale factor for FESEDFLUX
  call get_param(param_file, mdl, "MARBL_FESEDFLUX_SCALE_FACTOR", CS%fesedflux_scale_factor, &
                 "Conversion factor between FESEDFLUX file and MARBL units (umol / m^2 / d -> nmol / cm^2 / s)", &
                 default=1000. / 10000. / 86400.)

  CS%ntr = size(MARBL_instances%tracer_metadata)
  allocate(CS%ind_tr(CS%ntr))
  allocate(CS%tr_desc(CS%ntr))
  allocate(CS%tracer_data(CS%ntr))

  do m = 1, CS%ntr
    allocate(CS%tracer_data(m)%tr(isd:ied,jsd:jed,nz)) ; CS%tracer_data(m)%tr(:,:,:) = 0.0
    write(var_name(:),'(A)') trim(MARBL_instances%tracer_metadata(m)%short_name)
    write(desc_name(:),'(A)') trim(MARBL_instances%tracer_metadata(m)%long_name)
    write(units(:),'(A)') trim(MARBL_instances%tracer_metadata(m)%units)
    CS%tr_desc(m) = var_desc(trim(var_name), trim(units), trim(desc_name), caller=mdl)

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tracer_data(m)%tr(:,:,:)
    call query_vardesc(CS%tr_desc(m), name=var_name, &
                       caller="register_MARBL_tracers")
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, units = units, &
                         tr_desc=CS%tr_desc(m), registry_diags=.true., &
                         restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit, &
                         Tr_out=CS%tracer_data(m)%tr_ptr)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_MARBL_tracers")
  enddo

  ! Set up memory for saved state
  call setup_saved_state(MARBL_instances%surface_flux_saved_state, HI, GV, restart_CS, CS%tracers_may_reinit, &
                         CS%surface_flux_saved_state)
  call setup_saved_state(MARBL_instances%interior_tendency_saved_state, HI, GV, restart_CS, CS%tracers_may_reinit, &
                         CS%interior_tendency_saved_state)

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS

  call set_riv_flux_tracer_inds(CS)
  register_MARBL_tracers = .true.

end function register_MARBL_tracers

!> This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
!! and it sets up the tracer output.
subroutine initialize_MARBL_tracers(restart, day, G, GV, US, h, param_file, diag, OBC, CS, sponge_CSp)
  logical,                            intent(in) :: restart !< .true. if the fields have already been
                                                            !! read from a restart file.
  type(time_type), target,            intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: h !< Layer thicknesses [H ~> m or kg m-2]
  type(param_file_type),              intent(in) :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,            intent(in) :: diag !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(MARBL_tracers_CS),                pointer    :: CS   !< The control structure returned by a previous
                                                         !! call to register_MARBL_tracers.
  type(sponge_CS),                    pointer    :: sponge_CSp    !< A pointer to the control structure
                                                                  !! for the sponges, if they are in use.

! Local variables
  character(len=200) :: log_message
  character(len=48) :: name       ! A variable's name in a NetCDF file.
  character(len=100) :: longname   ! The long name of that variable.
  character(len=48) :: units      ! The units of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                  ! years m3 s-1 or years kg s-1.
  character(len=48) :: tracer_name
  logical :: fesedflux_has_edges, fesedflux_use_missing
  real    :: fesedflux_missing
  integer :: i, j, k, kbot, m, diag_size

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  CS%diag => diag

  ! Allocate memory for surface tracer fluxes
  allocate(CS%STF(SZI_(G), SZJ_(G), CS%ntr))
  allocate(CS%RIV_FLUXES(SZI_(G), SZJ_(G), CS%ntr))
  allocate(CS%SFO(SZI_(G), SZJ_(G), CS%sfo_cnt))
  CS%STF(:,:,:) = 0.
  CS%RIV_FLUXES(:,:,:) = 0.
  CS%SFO(:,:,:) = 0.

  ! Register diagnostics returned from MARBL (surface flux first, then interior tendency)
  call register_MARBL_diags(MARBL_instances%surface_flux_diags, diag, day, G, CS%surface_flux_diags)
  call register_MARBL_diags(MARBL_instances%interior_tendency_diags, diag, day, G, CS%interior_tendency_diags)

  ! Register per-tracer diagnostics computed from MARBL surface flux / interior tendency values
  allocate(CS%id_surface_flux_out(CS%ntr))
  allocate(CS%interior_tendency_out(CS%ntr))
  allocate(CS%interior_tendency_out_zint(CS%ntr))
  allocate(CS%interior_tendency_out_zint_100m(CS%ntr))
  do m=1,CS%ntr
    write(name, "(2A)") "STF_", trim(MARBL_instances%tracer_metadata(m)%short_name)
    write(longname, "(2A)") trim(MARBL_instances%tracer_metadata(m)%long_name), " Surface Flux"
    write(units, "(2A)") trim(MARBL_instances%tracer_metadata(m)%units), " m/s"
    CS%id_surface_flux_out(m) = register_diag_field("ocean_model", &
                                                    trim(name), &
                                                    diag%axesT1, & ! T => tracer grid? 1 => no vertical grid
                                                    day, &
                                                    trim(longname), &
                                                    trim(units))

    write(name, "(2A)") "J_", trim(MARBL_instances%tracer_metadata(m)%short_name)
    write(longname, "(2A)") trim(MARBL_instances%tracer_metadata(m)%long_name), " Source Sink Term"
    write(units, "(2A)") trim(MARBL_instances%tracer_metadata(m)%units), "/s"
    CS%interior_tendency_out(m)%id = register_diag_field("ocean_model", &
                                                         trim(name), &
                                                         diag%axesTL, & ! T=> tracer grid? L => layer center
                                                         day, &
                                                         trim(longname), &
                                                         trim(units))
    if (CS%interior_tendency_out(m)%id > 0) then
      allocate(CS%interior_tendency_out(m)%field_3d(SZI_(G),SZJ_(G), SZK_(G)))
      CS%interior_tendency_out(m)%field_3d(:,:,:) = 0.
    end if

    write(name, "(2A)") "Jint_", trim(MARBL_instances%tracer_metadata(m)%short_name)
    write(longname, "(2A)") trim(MARBL_instances%tracer_metadata(m)%long_name), " Source Sink Term Vertical Integral"
    write(units, "(2A)") trim(MARBL_instances%tracer_metadata(m)%units), " m/s"
    CS%interior_tendency_out_zint(m)%id = register_diag_field("ocean_model", &
                                                              trim(name), &
                                                              diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                              day, &
                                                              trim(longname), &
                                                              trim(units))
    if (CS%interior_tendency_out_zint(m)%id > 0) then
      allocate(CS%interior_tendency_out_zint(m)%field_2d(SZI_(G),SZJ_(G)))
      CS%interior_tendency_out_zint(m)%field_2d(:,:) = 0.
    end if

    write(name, "(2A)") "Jint_100m_", trim(MARBL_instances%tracer_metadata(m)%short_name)
    write(longname, "(2A)") trim(MARBL_instances%tracer_metadata(m)%long_name), &
                            " Source Sink Term Vertical Integral, 0-100m"
    write(units, "(2A)") trim(MARBL_instances%tracer_metadata(m)%units), " m/s"
    CS%interior_tendency_out_zint_100m(m)%id = register_diag_field("ocean_model", &
                                                                   trim(name), &
                                                                   diag%axesT1, &
                                                                   day, &
                                                                   trim(longname), &
                                                                   trim(units))
    if (CS%interior_tendency_out_zint_100m(m)%id > 0) then
      allocate(CS%interior_tendency_out_zint_100m(m)%field_2d(SZI_(G),SZJ_(G)))
      CS%interior_tendency_out_zint_100m(m)%field_2d(:,:) = 0.
    end if

  end do

  ! Register diagnostics for MOM to report that are not tracer specific
  CS%bot_flux_to_tend_id = register_diag_field("ocean_model", &
                                               "BOT_FLUX_TO_TEND", &
                                               diag%axesTL, & ! T=> tracer grid? L => layer center
                                               day, &
                                               "Conversion Factor for Bottom Flux -> Tend", &
                                               "1/m")

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=name, caller="initialize_MARBL_tracers")
    if ((.not. restart) .or. &
        (CS%tracers_may_reinit .and. &
         .not. query_initialized(CS%tracer_data(m)%tr(:,:,:), name, CS%restart_CSp))) then
      ! TODO: added the ongrid optional argument, but is there a good way to detect if the file is on grid?
      call MOM_initialize_tracer_from_Z(h, CS%tracer_data(m)%tr, G, GV, US, param_file, &
                                        CS%IC_file, name, ongrid=CS%ongrid)
      do k=1,GV%ke
        do j=G%jsc, G%jec
          do i=G%isc, G%iec
            ! Set negative tracer concentrations to 0
            if (CS%tracer_data(m)%tr(i,j,k) < 0) CS%tracer_data(m)%tr(i,j,k) = 0.
          end do
        end do
      end do
    end if
  end do

  ! Register diagnostics for per-category forcing fields
  if (CS%ice_ncat > 0) then
    allocate(CS%fracr_cat_id(CS%ice_ncat+1))
    allocate(CS%qsw_cat_id(CS%ice_ncat+1))
    do m=1,CS%ice_ncat+1
      write(name, "(A,I0)") "FRACR_CAT_", m
      write(longname, "(A,I0)") "Fraction of area in ice category ", m
      units = "fraction"
      CS%fracr_cat_id(m) = register_diag_field("ocean_model", &
                                               trim(name), &
                                               diag%axesT1, & ! T => tracer grid? 1 => no vertical grid
                                               day, &
                                               trim(longname), &
                                               trim(units))
      write(name, "(A,I0)") "QSW_CAT_", m
      write(longname, "(A,I0)") "Shortwave penetrating through ice category ", m
      units = "TODO: set units"
      CS%qsw_cat_id(m) = register_diag_field("ocean_model", &
                                             trim(name), &
                                             diag%axesT1, & ! T => tracer grid? 1 => no vertical grid
                                             day, &
                                             trim(longname), &
                                             trim(units))
    enddo
  endif

  ! Read initial fesedflux and feventflux fields
  ! (1) get vertical dimension
  !     -- comes from fesedflux_file, assume same dimension in feventflux
  !        (maybe these fields should be combined?)
  !     -- note: read_Z_edges treats depth as positive UP => 0 at surface, negative at depth
  fesedflux_use_missing = .false.
  call read_Z_edges(CS%fesedflux_file, "FESEDFLUXIN", CS%fesedflux_z_edges, CS%fesedflux_nz, &
                    fesedflux_has_edges, fesedflux_use_missing, fesedflux_missing, &
                    scale=US%m_to_Z)

  ! (2) Allocate memory for fesedflux and feventflux
  allocate(CS%fesedflux_in(SZI_(G), SZJ_(G), CS%fesedflux_nz))
  allocate(CS%feventflux_in(SZI_(G), SZJ_(G), CS%fesedflux_nz))
  allocate(CS%fesedflux_dz(SZI_(G), SZJ_(G), CS%fesedflux_nz))

  ! (3) Read data
  !     TODO: Add US term to scale
  call MOM_read_data(CS%fesedflux_file, "FESEDFLUXIN", CS%fesedflux_in(:,:,:), G%Domain, &
                     scale=CS%fesedflux_scale_factor)
  call MOM_read_data(CS%feventflux_file, "FESEDFLUXIN", CS%feventflux_in(:,:,:), G%Domain, &
                     scale=CS%fesedflux_scale_factor)

  ! (4) Relocate values that are below ocean bottom to layer that intersects bathymetry
  !     Remember, fesedflux_z_edges = 0 at surface and is < 0 below surface

  do k=CS%fesedflux_nz, 1, -1
    kbot = k + 1 ! level k is between z(k) and z(k+1)
    do j=G%jsc, G%jec
      do i=G%isc, G%iec
        if (G%mask2dT(i,j) == 0) cycle
        if (G%bathyT(i,j) + CS%fesedflux_z_edges(1) < 1e-8) then
          write(log_message, *) "Current implementation of fesedflux assumes G%bathyT >= first edge;", &
                                "first edge =", -CS%fesedflux_z_edges(1), &
                                "bathyT =", G%bathyT(i,j)
          call MOM_error(FATAL, log_message)
        end if
        ! Also figure out layer thickness while we're here
        CS%fesedflux_dz(i,j,k) = CS%fesedflux_z_edges(k) - CS%fesedflux_z_edges(kbot)
        ! If top interface is at or below ocean bottom, move flux in current layer up one
        ! and set thickness of current level to 0
        if (G%bathyT(i,j) + CS%fesedflux_z_edges(k) < 1e-8) then
          CS%fesedflux_in(i,j,k-1) = CS%fesedflux_in(i,j,k-1) + CS%fesedflux_in(i,j,k)
          CS%fesedflux_in(i,j,k) = 0.
          CS%feventflux_in(i,j,k-1) = CS%feventflux_in(i,j,k-1) + CS%feventflux_in(i,j,k)
          CS%feventflux_in(i,j,k) = 0.
          CS%fesedflux_dz(i,j,k) = 0.
        elseif (G%bathyT(i,j) + CS%fesedflux_z_edges(kbot) < 1e-8) then
          ! Otherwise, if lower interface is below bathymetry move interface to ocean bottom
          CS%fesedflux_dz(i,j,k) = G%bathyT(i,j) + CS%fesedflux_z_edges(k)
        end if
      end do
    end do
  end do

end subroutine initialize_MARBL_tracers

!> This subroutine is used to register tracer fields and subroutines
!! to be used with MOM.
subroutine register_MARBL_diags(MARBL_diags, diag, day, G, id_diags)

  type(marbl_diagnostics_type), intent(in)    :: MARBL_diags !< MARBL diagnostics from MARBL_instances
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
      if (id_diags(m)%id > 0) then
        allocate(id_diags(m)%field_2d(SZI_(G),SZJ_(G)))
        id_diags(m)%field_2d(:,:) = 0.
      end if
    else ! 3D field
      ! TODO: MARBL should provide v_extensive through MARBL_diags
      !       (for now, FESEDFLUX is the only one that should be true)
      id_diags(m)%id = register_diag_field("ocean_model", &
        trim(MARBL_diags%diags(m)%short_name), &
        diag%axesTL, & ! T=> tracer grid? L => layer center
        day, &
        trim(MARBL_diags%diags(m)%long_name), &
        trim(MARBL_diags%diags(m)%units), &
        v_extensive=(trim(MARBL_diags%diags(m)%short_name).eq."FESEDFLUX"))
      if (id_diags(m)%id > 0) then
        allocate(id_diags(m)%field_3d(SZI_(G),SZJ_(G), SZK_(G)))
        id_diags(m)%field_3d(:,:,:) = 0.
      end if
    end if
  end do

end subroutine register_MARBL_diags

!> This subroutine allocates memory for saved state fields and registers them in the restart files
subroutine setup_saved_state(MARBL_saved_state, HI, GV, restart_CS, tracers_may_reinit, local_saved_state)

  type(marbl_saved_state_type),                  intent(in)    :: MARBL_saved_state !< MARBL saved state from
                                                                                    !! MARBL_instances
  type(hor_index_type),                          intent(in)    :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),                       intent(in)    :: GV   !< The ocean's vertical grid structure
  type(MOM_restart_CS), pointer,                 intent(in)    :: restart_CS !< control structure to add saved state
                                                                             !! to restarts
  logical,                                       intent(in)    :: tracers_may_reinit  !< used to determine mandatory
                                                                                      !! flag in restart
  type(saved_state_for_MARBL_type), allocatable, intent(inout) :: local_saved_state(:) !< allocatable array for local
                                                                                       !! saved state

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
              KPP_CSp, nonLocalTrans, evap_CFL_limit, minimum_forcing_depth)
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
  type(KPP_CS),  optional, pointer    :: KPP_CSp  !< KPP control structure
  real,          optional, intent(in) :: nonLocalTrans(:,:,:) !< Non-local transport [nondim]
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [m]

! Local variables
  character(len=256) :: log_message
  real, dimension(SZI_(G),SZJ_(G)) :: ref_mask ! Mask for 2D MARBL diags using ref_depth
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: bot_flux_to_tend
  real :: cum_bftt_dz     ! sum of bot_flux_to_tend * dz from the bottom layer to current layer
  real :: sfc_val  ! The surface value for the tracers.
  real :: Isecs_per_year  ! The number of seconds in a year.
  real :: year            ! The time in years.
  integer :: secs, days   ! Integer components of the time type.
  real, dimension(0:GV%ke) :: zi  ! z-coordinate interface depth
  real, dimension(GV%ke) :: zc, dz  ! z-coordinate layer center depth and cell thickness
  integer :: i, j, k, is, ie, js, je, nz, m
  real :: ndep_conversion ! mol L-2 T-1 -> nmol cm-2 s-1

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  ! First two terms get us to mol m-2 s-1; 1e5 gets us nmol cm-2 s-1, which is what MARBL wants
  ndep_conversion = (US%m_to_L)**2 * US%s_to_T * 1.e5

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  ! (1) Compute surface fluxes
  ! FIXME: MARBL can handle computing surface fluxes for all columns simultaneously
  !        I was just thinking going column-by-column at first might be easier
  do j=js,je
    do i=is,ie
      ! i. only want ocean points in this loop
      if (G%mask2dT(i,j) == 0) cycle

      ! ii. Load proper column data
      !     * surface flux forcings
      !       These fields are getting the correct data
      !       TODO: if top layer is vanishly thin, do we actually want (e.g.) top 5m average temp / salinity?
      !             How does MOM pass SST and SSS to GFDL coupler? (look in core.F90?)
      if (CS%sss_ind > 0) MARBL_instances%surface_flux_forcings(CS%sss_ind)%field_0d(1) = tv%S(i,j,1)
      if (CS%sst_ind > 0) MARBL_instances%surface_flux_forcings(CS%sst_ind)%field_0d(1) = tv%T(i,j,1)
      if (CS%ifrac_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%ifrac_ind)%field_0d(1) = fluxes%ice_fraction(i,j)
      ! MARBL wants u10_sqr in (cm/s)^2
      if (CS%u10_sqr_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%u10_sqr_ind)%field_0d(1) = fluxes%u10_sqr(i,j) * &
                                                                            (US%L_t_to_m_s * cm_per_m)**2
      ! mct_driver/ocn_cap_methods:93 -- ice_ocean_boundary%p(i,j) comes from coupler
      ! We may need a new ice_ocean_boundary%p_atm because %p includes ice in GFDL driver
      if (CS%atmpress_ind > 0) then
        if (associated(fluxes%p_surf_full)) then
          MARBL_instances%surface_flux_forcings(CS%atmpress_ind)%field_0d(1) = fluxes%p_surf_full(i,j) * &
                                                                               ((US%R_to_kg_m3 * US%L_T_to_m_s**2) * &
                                                                                atm_per_Pa)
        else
          ! hardcode value of 1 atm (can't figure out how to get this from solo_driver)
          MARBL_instances%surface_flux_forcings(CS%atmpress_ind)%field_0d(1) = 1 * (US%R_to_kg_m3 * US%L_T_to_m_s**2)
        endif
      endif
      !       These are okay, but need option to come in from coupler
      if (CS%xco2_ind > 0) MARBL_instances%surface_flux_forcings(CS%xco2_ind)%field_0d(1) = fluxes%atm_co2(i,j)
      if (CS%xco2_alt_ind > 0) MARBL_instances%surface_flux_forcings(CS%xco2_alt_ind)%field_0d(1) = fluxes%atm_alt_co2(i,j)

      !       These are okay, but need option to read in from file
      if (CS%dust_dep_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%dust_dep_ind)%field_0d(1) = fluxes%dust_flux(i,j)

      if (CS%fe_dep_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%fe_dep_ind)%field_0d(1) = fluxes%iron_flux(i,j)

      !       MARBL wants ndep in (nmol/cm^2/s)
      if (CS%nox_flux_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%nox_flux_ind)%field_0d(1) = fluxes%noy_dep(i,j) * &
                                                                             ndep_conversion
      if (CS%nhy_flux_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%nhy_flux_ind)%field_0d(1) = fluxes%nhx_dep(i,j) * &
                                                                             ndep_conversion

      !     * tracers at surface
      !       TODO: average over some shallow depth (e.g. 5m)
      do m=1,CS%ntr
        MARBL_instances%tracers_at_surface(1,m) = CS%tracer_data(m)%tr(i,j,1)
      end do

      !     * surface flux saved state
      do m=1,size(MARBL_instances%surface_flux_saved_state%state)
      !       (currently only 2D fields are saved from surface_flux_compute())
        MARBL_instances%surface_flux_saved_state%state(m)%field_2d(1) = CS%surface_flux_saved_state(m)%field_2d(i,j)
      end do

      ! iii. Compute surface fluxes in MARBL
      call MARBL_instances%surface_flux_compute()
      if (MARBL_instances%StatusLog%labort_marbl) then
        call MARBL_instances%StatusLog%log_error_trace("MARBL_instances%surface_flux_compute()", &
                                                       "MARBL_tracers_column_physics")
      end if
      call print_marbl_log(MARBL_instances%StatusLog)
      call MARBL_instances%StatusLog%erase()

      ! iv. Copy output that MOM6 needs to hold on to
      !     * saved state
      do m=1,size(MARBL_instances%surface_flux_saved_state%state)
        CS%surface_flux_saved_state(m)%field_2d(i,j) = MARBL_instances%surface_flux_saved_state%state(m)%field_2d(1)
      end do

      !     * diagnostics
      do m=1,size(MARBL_instances%surface_flux_diags%diags)
        ! All diags are 2D coming from surface
        if (CS%surface_flux_diags(m)%id > 0) &
          CS%surface_flux_diags(m)%field_2d(i,j) = real(MARBL_instances%surface_flux_diags%diags(m)%field_2d(1))
      end do

      !     * Surface tracer flux
      CS%STF(i,j,:) = MARBL_instances%surface_fluxes(1,:) * m_per_cm

      !     * Surface flux output
      do m=1,CS%sfo_cnt
        CS%SFO(i,j,m) = MARBL_instances%surface_flux_output%outputs_for_GCM(m)%forcing_field_0d(1)
      end do

    end do
  end do
  ! Add River Fluxes to STF
  if (CS%tracer_inds%no3_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%no3_ind) = fluxes%no3_riv_flux(:,:)
  if (CS%tracer_inds%po4_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%po4_ind) = fluxes%po4_riv_flux(:,:)
  if (CS%tracer_inds%don_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%don_ind) = fluxes%don_riv_flux(:,:)
  if (CS%tracer_inds%donr_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%donr_ind) = fluxes%donr_riv_flux(:,:)
  if (CS%tracer_inds%dop_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%dop_ind) = fluxes%dop_riv_flux(:,:)
  if (CS%tracer_inds%dopr_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%dopr_ind) = fluxes%dopr_riv_flux(:,:)
  if (CS%tracer_inds%sio3_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%sio3_ind) = fluxes%sio3_riv_flux(:,:)
  if (CS%tracer_inds%fe_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%fe_ind) = fluxes%fe_riv_flux(:,:)
  if (CS%tracer_inds%doc_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%doc_ind) = fluxes%doc_riv_flux(:,:)
  if (CS%tracer_inds%docr_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%docr_ind) = fluxes%docr_riv_flux(:,:)
  if (CS%tracer_inds%alk_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%alk_ind) = fluxes%alk_riv_flux(:,:)
  if (CS%tracer_inds%alk_alt_co2_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%alk_alt_co2_ind) = fluxes%alk_alt_co2_riv_flux(:,:)
  if (CS%tracer_inds%dic_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%dic_ind) = fluxes%dic_riv_flux(:,:)
  if (CS%tracer_inds%dic_alt_co2_ind > 0) &
    CS%RIV_FLUXES(:,:,CS%tracer_inds%dic_alt_co2_ind) = fluxes%dic_alt_co2_riv_flux(:,:)
  ! fluxes%*_riv_flux are conc m/s, in_flux_optional expects time-integrated flux (conc m)
  CS%RIV_FLUXES(:,:,:) = CS%RIV_FLUXES(:,:,:) * dt

  ! (2) Post surface fluxes and their diagnostics (currently all 2D)
  do m=1,CS%ntr
    if (CS%id_surface_flux_out(m) > 0) &
      call post_data(CS%id_surface_flux_out(m), CS%STF(:,:,m), CS%diag)
  enddo
  do m=1,size(CS%surface_flux_diags)
    if (CS%surface_flux_diags(m)%id > 0) &
      call post_data(CS%surface_flux_diags(m)%id, CS%surface_flux_diags(m)%field_2d(:,:), CS%diag)
  enddo

  ! (3) Apply surface fluxes via vertical diffusion
  ! Compute KPP nonlocal term if necessary
  if (present(KPP_CSp)) then
    if (associated(KPP_CSp) .and. present(nonLocalTrans)) then
      do m=1,CS%ntr
        call KPP_NonLocalTransport(KPP_CSp, G, GV, h_old, nonLocalTrans, CS%STF(:,:,m), &
                                   dt, CS%diag, CS%tracer_data(m)%tr_ptr, &
                                   CS%tracer_data(m)%tr(:,:,:), flux_scale=GV%Z_to_H * US%T_to_s)
      enddo
    endif
  endif

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,CS%ntr
      do k=1,nz ;do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tracer_data(m)%tr(:,:,:) , dt, fluxes, h_work, &
          evap_CFL_limit, minimum_forcing_depth, in_flux_optional=CS%RIV_FLUXES(:,:,m))
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tracer_data(m)%tr(:,:,:), G, GV, &
                           sfc_flux=GV%Rho0 * CS%STF(:,:,m) * US%T_to_s)
    enddo
  else
    do m=1,CS%ntr
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tracer_data(m)%tr(:,:,:), G, GV, &
                           sfc_flux=GV%Rho0 * CS%STF(:,:,m) * US%T_to_s)
    enddo
  endif

  ! (4) Compute interior tendencies
  bot_flux_to_tend(:, :, :) = 0.
  do j=js,je
    do i=is,ie
      ! i. only want ocean points in this loop
      if (G%mask2dT(i,j) == 0) cycle

      ! ii. Set up vertical domain and bot_flux_to_tend
      MARBL_instances%domain%kmt = GV%ke
      ! Calculate depth of interface by building up thicknesses from the bottom (top interface is always 0)
      ! MARBL wants this to be positive-down
      zi(GV%ke) = G%bathyT(i,j)
      MARBL_instances%bot_flux_to_tend(:) = 0.
      cum_bftt_dz = 0.
      do k = GV%ke, 1, -1
        ! TODO: if we move this above vertical mixing, use h_old
        dz(k) = h_new(i,j,k)*GV%H_to_Z ! cell thickness
        zc(k) = zi(k) - 0.5 * dz(k)
        zi(k-1) = zi(k) - dz(k)
        if (G%bathyT(i,j) - zi(k-1) <= CS%bot_flux_mix_thickness) then
          MARBL_instances%bot_flux_to_tend(k) = CS%bfmt_r
          cum_bftt_dz = cum_bftt_dz + MARBL_instances%bot_flux_to_tend(k) * dz(k)
        elseif (G%bathyT(i,j) - zi(k) < CS%bot_flux_mix_thickness) then
          ! MARBL_instances%bot_flux_to_tend(k) = (1. - (G%bathyT(i,j) - zi(k)) * CS%bfmt_r) / dz(k)
          MARBL_instances%bot_flux_to_tend(k) = (1. - cum_bftt_dz) / dz(k)
        end if
      enddo
      if (G%bathyT(i,j) - zi(0) < CS%bot_flux_mix_thickness) &
        MARBL_instances%bot_flux_to_tend(:) = MARBL_instances%bot_flux_to_tend(:) * &
                                              CS%bot_flux_mix_thickness / (G%bathyT(i,j) - zi(0))
      if (CS%bot_flux_to_tend_id > 0) &
        bot_flux_to_tend(i, j, :) = MARBL_instances%bot_flux_to_tend(:)
      ! When MARBL is mks, we can drop the m_per_cm conversion
      MARBL_instances%bot_flux_to_tend(:) = MARBL_instances%bot_flux_to_tend(:) * m_per_cm

      ! mks -> cgs
      ! zw(1:nz) is bottom cell depth so no element of zw = 0, it is assumed to be top layer depth
      MARBL_instances%domain%zw(:) = zi(1:GV%ke) * cm_per_m
      MARBL_instances%domain%zt(:) = zc(:) * cm_per_m
      MARBL_instances%domain%delta_z(:) = dz(:) * cm_per_m

      ! iii. Load proper column data
      !      * Forcing Fields
      !       These fields are getting the correct data
      if (CS%potemp_ind > 0) MARBL_instances%interior_tendency_forcings(CS%potemp_ind)%field_1d(1,:) = tv%T(i,j,:)
      if (CS%salinity_ind > 0) MARBL_instances%interior_tendency_forcings(CS%salinity_ind)%field_1d(1,:) = tv%S(i,j,:)

      !       This are okay, but need option to read in from file
      !       (Same as dust_dep_ind for surface_flux_forcings)
      if (CS%dustflux_ind > 0) &
        MARBL_instances%interior_tendency_forcings(CS%dustflux_ind)%field_0d(1) = &
            fluxes%dust_flux(i,j)

      !        TODO: Support PAR (currently just using single subcolumn)
      !              (Look for Pen_sw_bnd?)
      if (CS%PAR_col_frac_ind > 0) then
        ! second index is num_subcols, not depth
        !MARBL_instances%interior_tendency_forcings(CS%PAR_col_frac_ind)%field_1d(1,:) = fluxes%fracr_cat(i,j,:)
        if (CS%use_ice_category_fields) then
          MARBL_instances%interior_tendency_forcings(CS%PAR_col_frac_ind)%field_1d(1,:) = fluxes%fracr_cat(i,j,:)
        else
          MARBL_instances%interior_tendency_forcings(CS%PAR_col_frac_ind)%field_1d(1,1) = 1
        end if
      end if
      if (CS%surf_shortwave_ind > 0) then
        ! second index is num_subcols, not depth
        if (CS%use_ice_category_fields) then
          MARBL_instances%interior_tendency_forcings(CS%surf_shortwave_ind)%field_1d(1,:) = fluxes%qsw_cat(i,j,:)
        else
          MARBL_instances%interior_tendency_forcings(CS%surf_shortwave_ind)%field_1d(1,1) = fluxes%sw(i,j)
        end if
      end if

      !        TODO: In POP, pressure comes from a function in state_mod.F90; I don't see a similar function here
      !              This formulation is from Levitus 1994, and I think it belongs in MOM_EOS.F90?
      !              Converts depth [m] -> pressure [bars]
      !        NOTE: Andrew recommends using GV%H_to_Pa
      if (CS%pressure_ind > 0) MARBL_instances%interior_tendency_forcings(CS%pressure_ind)%field_1d(1,:) = &
          0.0598088*(exp(-0.025*zc(:)) - 1) + 0.100766*zc(:) + 2.28405e-7*(zc(:)**2)

      if (CS%fesedflux_ind > 0) then
        MARBL_instances%interior_tendency_forcings(CS%fesedflux_ind)%field_1d(1,:) = 0.
        call reintegrate_column(CS%fesedflux_nz, &
                                CS%fesedflux_dz(i,j,:) * (sum(dz(:)) / G%bathyT(i,j)), &
                                CS%fesedflux_in(i,j,:) + CS%feventflux_in(i,j,:), &
                                GV%ke, &
                                dz(:), &
                                0., &
                                MARBL_instances%interior_tendency_forcings(CS%fesedflux_ind)%field_1d(1,:))
      end if

      !        TODO: add ability to read these fields from file
      !              also, add constant values to CS
      if (CS%o2_scalef_ind > 0) MARBL_instances%interior_tendency_forcings(CS%o2_scalef_ind)%field_1d(1,:) = 1
      if (CS%remin_scalef_ind > 0) MARBL_instances%interior_tendency_forcings(CS%remin_scalef_ind)%field_1d(1,:) = 1

      !      * Column Tracers
      !        NOTE: POP averages previous two timesteps, should we do that too?
      do m=1,CS%ntr
        MARBL_instances%tracers(m, :) = CS%tracer_data(m)%tr(i,j,:)
      end do

      !     * interior tendency saved state
      !       (currently only 3D fields are saved from interior_tendency_compute())
      do m=1,size(MARBL_instances%interior_tendency_saved_state%state)
        MARBL_instances%interior_tendency_saved_state%state(m)%field_3d(:,1) = &
            CS%interior_tendency_saved_state(m)%field_3d(i,j,:)
      end do

      ! iv. Compute interior tendencies in MARBL
      call MARBL_instances%interior_tendency_compute()
      if (MARBL_instances%StatusLog%labort_marbl) then
        call MARBL_instances%StatusLog%log_error_trace("MARBL_instances%interior_tendency_compute()", &
                                                       "MARBL_tracers_column_physics")
      end if
      call print_marbl_log(MARBL_instances%StatusLog, G, i, j)
      call MARBL_instances%StatusLog%erase()

      ! v. Apply tendencies immediately
      !    First pass - Euler step; if stability issues, we can do something different (subcycle?)
      do m=1,CS%ntr
        CS%tracer_data(m)%tr(i,j,:) = CS%tracer_data(m)%tr(i,j,:) + &
                                      G%mask2dT(i,j)*dt*MARBL_instances%interior_tendencies(m,:)
      end do

      ! vi. Copy output that MOM6 needs to hold on to
      !     * saved state
      do m=1,size(MARBL_instances%interior_tendency_saved_state%state)
        CS%interior_tendency_saved_state(m)%field_3d(i,j,:) = &
            MARBL_instances%interior_tendency_saved_state%state(m)%field_3d(:,1)
      end do

      !     * diagnostics
      do m=1,size(MARBL_instances%interior_tendency_diags%diags)
        if (CS%interior_tendency_diags(m)%id > 0) then
          if (allocated(CS%interior_tendency_diags(m)%field_2d)) then
            ! Only copy values if ref_depth < bathyT
            if (G%bathyT(i,j) > real(MARBL_instances%interior_tendency_diags%diags(m)%ref_depth)) then
              CS%interior_tendency_diags(m)%field_2d(i,j) = &
                  real(MARBL_instances%interior_tendency_diags%diags(m)%field_2d(1))
            end if
          else ! not a 2D diagnostic
            CS%interior_tendency_diags(m)%field_3d(i,j,:) = &
                real(MARBL_instances%interior_tendency_diags%diags(m)%field_3d(:,1))
          end if
        end if
      end do
      !     * tendency values themselves (and vertical integrals of them)
      do m=1,CS%ntr
        if (allocated(CS%interior_tendency_out(m)%field_3d)) then
          CS%interior_tendency_out(m)%field_3d(i,j,:) = G%mask2dT(i,j)*MARBL_instances%interior_tendencies(m,:)
        end if

        if (allocated(CS%interior_tendency_out_zint(m)%field_2d)) then
          CS%interior_tendency_out_zint(m)%field_2d(i,j) = G%mask2dT(i,j) * &
                                                           sum(dz(:) * MARBL_instances%interior_tendencies(m,:))
        end if

        if (allocated(CS%interior_tendency_out_zint_100m(m)%field_2d)) then
          CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) = 0.
          do k=1,GV%ke
            if (zi(k) < 100.) then
              CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) = &
                  CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) + dz(k) * MARBL_instances%interior_tendencies(m,k)
            elseif (zi(k-1) < 100.) then
              CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) = &
                  CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) + dz(k) * &
                                                                        ((100. - zi(k-1)) / (zi(k) - zi(k-1))) * &
                                                                        MARBL_instances%interior_tendencies(m,k)
            else
              exit
            end if
          end do
          CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) = &
              G%mask2dT(i,j)*CS%interior_tendency_out_zint_100m(m)%field_2d(i,j)
        end if
      end do

    end do
  end do

  ! (5) Post diagnostics from our buffer
  !     i. Interior tendency diagnostics (mix of 2D and 3D)
  !     ii. Interior tendencies themselves
  !     iii. Forcing fields
  if (CS%bot_flux_to_tend_id > 0) &
    call post_data(CS%bot_flux_to_tend_id, bot_flux_to_tend(:, :, :), CS%diag)

  do m=1,size(CS%interior_tendency_diags)
    if (CS%interior_tendency_diags(m)%id > 0) then
      if (allocated(CS%interior_tendency_diags(m)%field_2d)) then
        if (real(MARBL_instances%interior_tendency_diags%diags(m)%ref_depth) == 0.) then
          call post_data(CS%interior_tendency_diags(m)%id, CS%interior_tendency_diags(m)%field_2d(:,:), CS%diag)
        else ! non-zero ref-depth
          ref_mask(:, :) = 0.
          do j=js,je ; do i=is,ie
            if (G%bathyT(i,j) > real(MARBL_instances%interior_tendency_diags%diags(m)%ref_depth)) &
              ref_mask(i,j) = 1.
          end do ; end do
          call post_data(CS%interior_tendency_diags(m)%id, CS%interior_tendency_diags(m)%field_2d(:,:), &
                         CS%diag, mask=ref_mask(:,:))
        end if
      elseif (allocated(CS%interior_tendency_diags(m)%field_3d)) then
        call post_data(CS%interior_tendency_diags(m)%id, CS%interior_tendency_diags(m)%field_3d(:,:,:), CS%diag)
      else
        write(log_message, "(A, I0, A, I0, A)") "Diagnostic number ", m, " post id ", &
                                                CS%interior_tendency_diags(m)%id," did not allocate 2D or 3D array"
        call MOM_error(FATAL, log_message)
      end if
    end if
  end do
  do m=1,CS%ntr
    if (allocated(CS%interior_tendency_out(m)%field_3d)) &
      call post_data(CS%interior_tendency_out(m)%id, CS%interior_tendency_out(m)%field_3d(:,:,:), CS%diag)
    if (allocated(CS%interior_tendency_out_zint(m)%field_2d)) &
      call post_data(CS%interior_tendency_out_zint(m)%id, CS%interior_tendency_out_zint(m)%field_2d(:,:), CS%diag)
    if (allocated(CS%interior_tendency_out_zint_100m(m)%field_2d)) &
      call post_data(CS%interior_tendency_out_zint_100m(m)%id, CS%interior_tendency_out_zint_100m(m)%field_2d(:,:), &
                     CS%diag)
  end do
  if (CS%ice_ncat > 0) then
    do m=1,CS%ice_ncat+1
      if (CS%fracr_cat_id(m) > 0) call post_data(CS%fracr_cat_id(m), fluxes%fracr_cat(:,:,m), CS%diag)
      if (CS%qsw_cat_id(m) > 0) call post_data(CS%qsw_cat_id(m), fluxes%qsw_cat(:,:,m), CS%diag)
    end do
  endif


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
      stocks(m) = stocks(m) + CS%tracer_data(m)%tr(i,j,k) * &
                             (G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k))
    enddo ; enddo ; enddo
    stocks(m) = GV%H_to_kg_m2 * stocks(m)
  enddo
  MARBL_tracer_stock = CS%ntr

end function MARBL_tracer_stock

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
subroutine MARBL_tracers_surface_state(sfc_state, G, US, CS)
  type(ocean_grid_type),   intent(in)    :: G   !< The ocean's grid structure.
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                                      !! describe the surface state of the ocean.
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(MARBL_tracers_CS),  pointer       :: CS  !< The control structure returned by a previous
                                                !! call to register_MARBL_tracers.

  ! Local variables
  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(CS)) return

  if (allocated(sfc_state%sfc_co2)) then
    do j=js,je ; do i=is,ie
      !  nmol/cm^2/s (positive down) to kg CO2/m^2/s (positive down)
      sfc_state%sfc_co2(i,j) = US%kg_m2s_to_RZ_T * (44.0e-8*CS%SFO(i,j,CS%flux_co2_ind))
    enddo ; enddo
  endif

end subroutine MARBL_tracers_surface_state

!> Copy the requested interior tendency output field into an array.
subroutine MARBL_tracers_get_output_for_GCM(name, G, GV, array, CS)

  use marbl_settings_mod, only : output_for_GCM_iopt_total_Chl_3d

  character(len=*),         intent(in)    :: name   !< Name of requested tracer.
  type(ocean_grid_type),    intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(inout) :: array  !< Array filled by this routine.
  type(MARBL_tracers_CS),   pointer       :: CS     !< Pointer to the control structure for this module.

  character(len=128), parameter :: sub_name = 'MARBL_tracers_get_output_for_GCM'
  character(len=128) :: log_message
  integer :: i, j, is, ie, js, je, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  array(:,:,:) = 0.0
  select case(trim(name))
    case ('Chl')
      do j=js,je ; do i=is,ie
        do m=1,size(CS%tracer_data)
          MARBL_instances%tracers(m,:) = CS%tracer_data(m)%tr(i,j,:)
        end do
        if (G%mask2dT(i,j) /= 0) then
          call MARBL_instances%get_output_for_GCM(output_for_GCM_iopt_total_Chl_3d, array(i,j,:))
          if (MARBL_instances%StatusLog%labort_marbl) &
            call print_marbl_log(MARBL_instances%StatusLog, G, i, j)
        end if
      end do ; end do
    case DEFAULT
      write(log_message, "(3A)") "'", trim(name), "' is not a valid interior tendency output field name"
      call MOM_error(FATAL, log_message)
  end select

end subroutine MARBL_tracers_get_output_for_GCM

!> Clean up any allocated memory after the run.
subroutine MARBL_tracers_end(CS)
  type(MARBL_tracers_CS), pointer, intent(inout) :: CS !< The control structure returned by a previous
                                                       !! call to register_MARBL_tracers.

  integer :: m

  call print_marbl_log(MARBL_instances%StatusLog)
  call MARBL_instances%StatusLog%erase()
  call MARBL_instances%shutdown()
  ! TODO: print MARBL timers to stdout as well

  if (associated(CS)) then
    if (allocated(CS%tracer_data)) then
      do m=1,CS%ntr
        if (associated(CS%tracer_data(m)%tr)) deallocate(CS%tracer_data(m)%tr)
      end do
      deallocate(CS%tracer_data)
    end if
    if (allocated(CS%ind_tr)) &
      deallocate(CS%ind_tr)
    if (allocated(CS%id_surface_flux_out)) &
      deallocate(CS%id_surface_flux_out)
    if (allocated(CS%fracr_cat_id)) &
      deallocate(CS%fracr_cat_id)
    if (allocated(CS%qsw_cat_id)) &
      deallocate(CS%qsw_cat_id)
    if (allocated(CS%STF)) &
      deallocate(CS%STF)
    if (allocated(CS%RIV_FLUXES)) &
      deallocate(CS%RIV_FLUXES)
    if (allocated(CS%SFO)) &
      deallocate(CS%SFO)
    deallocate(CS)
  endif
end subroutine MARBL_tracers_end

subroutine set_riv_flux_tracer_inds(CS)

  type(MARBL_tracers_CS), pointer, intent(inout) :: CS   !< The MARBL tracers control structure

  character(len=256) :: log_message
  character(len=48) :: name       ! A variable's name in a NetCDF file.
  integer :: m

  ! Initialize tracers from file (unless they were initialized by restart file)
  ! Also save indices of tracers that have river fluxes
  CS%tracer_inds%no3_ind = 0
  CS%tracer_inds%po4_ind = 0
  CS%tracer_inds%don_ind = 0
  CS%tracer_inds%donr_ind = 0
  CS%tracer_inds%dop_ind = 0
  CS%tracer_inds%dopr_ind = 0
  CS%tracer_inds%sio3_ind = 0
  CS%tracer_inds%fe_ind = 0
  CS%tracer_inds%doc_ind = 0
  CS%tracer_inds%docr_ind = 0
  CS%tracer_inds%alk_ind = 0
  CS%tracer_inds%alk_alt_co2_ind = 0
  CS%tracer_inds%dic_ind = 0
  CS%tracer_inds%dic_alt_co2_ind = 0
  do m=1,CS%ntr
    name = MARBL_instances%tracer_metadata(m)%short_name
    if (trim(name) == "NO3") then
      CS%tracer_inds%no3_ind = m
    elseif (trim(name) == "PO4") then
       CS%tracer_inds%po4_ind = m
    elseif (trim(name) == "DON") then
       CS%tracer_inds%don_ind = m
    elseif (trim(name) == "DONr") then
       CS%tracer_inds%donr_ind = m
    elseif (trim(name) == "DOP") then
       CS%tracer_inds%dop_ind = m
    elseif (trim(name) == "DOPr") then
       CS%tracer_inds%dopr_ind = m
    elseif (trim(name) == "SiO3") then
       CS%tracer_inds%sio3_ind = m
    elseif (trim(name) == "Fe") then
       CS%tracer_inds%fe_ind = m
    elseif (trim(name) == "DOC") then
       CS%tracer_inds%doc_ind = m
    elseif (trim(name) == "DOCr") then
       CS%tracer_inds%docr_ind = m
    elseif (trim(name) == "ALK") then
       CS%tracer_inds%alk_ind = m
    elseif (trim(name) == "ALK_ALT_CO2") then
       CS%tracer_inds%alk_alt_co2_ind = m
    elseif (trim(name) == "DIC") then
       CS%tracer_inds%dic_ind = m
    elseif (trim(name) == "DIC_ALT_CO2") then
       CS%tracer_inds%dic_alt_co2_ind = m
    end if
  end do

  ! Log indices for each tracer to ensure we set them all correctly
  write(log_message, "(A,I0)") "NO3 index: ", CS%tracer_inds%no3_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "PO4 index: ", CS%tracer_inds%po4_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DON index: ", CS%tracer_inds%don_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DONr index: ", CS%tracer_inds%donr_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DOP index: ", CS%tracer_inds%dop_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DOPr index: ", CS%tracer_inds%dopr_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "SiO3 index: ", CS%tracer_inds%sio3_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "Fe index: ", CS%tracer_inds%fe_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DOC index: ", CS%tracer_inds%doc_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DOCr index: ", CS%tracer_inds%docr_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "ALK index: ", CS%tracer_inds%alk_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "ALK_ALT_CO2 index: ", CS%tracer_inds%alk_alt_co2_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DIC index: ", CS%tracer_inds%dic_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DIC_ALT_CO2 index: ", CS%tracer_inds%dic_alt_co2_ind
  call MOM_error(NOTE, log_message)

end subroutine set_riv_flux_tracer_inds

! TODO: some log messages come from a specific grid point, and this routine
!       needs to include the location in the preamble
!> This subroutine writes the contents of the MARBL log using MOM_error(NOTE, ...).
subroutine print_marbl_log(log_to_print, G, i, j)

  use marbl_logging, only : marbl_status_log_entry_type
  use marbl_logging, only : marbl_log_type
  use MOM_coms,      only : PE_here

  class(marbl_log_type),           intent(in) :: log_to_print  !< MARBL log to include in MOM6 logfile
  type(ocean_grid_type), optional, intent(in) :: G             !< The ocean's grid structure
  integer,               optional, intent(in) :: i             !< i of (i,j) index of column providing the log
  integer,               optional, intent(in) :: j             !< j of (i,j) index of column providing the log

  character(len=*), parameter :: subname = 'MARBL_tracers:print_marbl_log'
  character(len=256)          :: message_prefix, message_location, log_message
  type(marbl_status_log_entry_type), pointer :: tmp
  integer :: msg_lev, elem_old

  ! elem_old is used to keep track of whether all messages are coming from the same point
  elem_old = -1
  write(message_prefix, "(A,I0,A)") '(Task ', PE_here(), ')'

  tmp => log_to_print%FullLog
  do while (associated(tmp))
    ! 1) Do I need to write this message? Yes, if all tasks should write this
    !    or if I am master_task
    if ((.not. tmp%lonly_master_writes) .or. is_root_PE()) then
      ! 2) Print message location? (only if ElementInd changed and is positive; requires G)
      if ((present(G)) .and. (tmp%ElementInd .ne. elem_old)) then
        if (tmp%ElementInd .gt. 0) then
          if (present(i) .and. present(j)) then
            write(message_location, "(A,F8.3,A,F7.3,A,I0,A,I0,A,I0)") &
                 'Message from (lon, lat) (', G%geoLonT(i,j), ', ', &
                 G%geoLatT(i,j), '), which is global (i,j) (', &
                 i + G%HI%idg_offset, ', ', j + G%HI%jdg_offset, &
                 '). Level: ', tmp%ElementInd
          else
            write(message_location, "(A)") "Grid cell responsible for message is unknown"
          !   i_loc = marbl_col_to_pop_i(tmp%ElementInd, iblock)
          !   j_loc = marbl_col_to_pop_j(tmp%ElementInd, iblock)
          !   write(message_location, "(A,F8.3,A,F7.3,A,I0,A,I0,A)") &
          !        'Message from (lon, lat) (', TLOND(i_loc, j_loc, iblock), &
          !        ", ", TLATD(i_loc, j_loc, iblock), '), which is global (i,j) (', &
          !        this_block%i_glob(i_loc), ', ', this_block%j_glob(j_loc), ')'
          end if ! i,j present
          ! master task does not need prefix
          if (is_root_PE()) then
            write(log_message, "(A)") trim(message_location)
            msg_lev = NOTE
          else
            write(log_message, "(A,1X,A)") trim(message_prefix), trim(message_location)
            msg_lev = WARNING
          end if ! print message prefix?
          call MOM_error(msg_lev, log_message, all_print=.true.)
        end if   ! ElementInd > 0
        elem_old = tmp%ElementInd
      end if     ! ElementInd /= elem_old

      ! 3) Write message from the log
      ! master task does not need prefix
      if (is_root_PE()) then
        write(log_message, "(A)") trim(tmp%LogMessage)
        msg_lev = NOTE
      else
        write(log_message, "(A,1X,A)") trim(message_prefix), trim(tmp%LogMessage)
        msg_lev = WARNING
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

! TODO: fix the comment below
!> \namespace MARBL_tracers
!!
!!    This file contains an example of the code that is needed to set
!!  up and use a set (in this case two) of dynamically passive tracers
!!  for diagnostic purposes.  The tracers here are dye tracers which
!!  are set to 1 within the geographical region specified. The depth
!!  which a tracer is set is determined by calculating the depth from
!!  the seafloor upwards through the column.

end module MARBL_tracers
