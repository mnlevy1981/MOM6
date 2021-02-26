module marbl_forcing_type_main

use MOM_error_handler,        only : MOM_error, WARNING
use MOM_file_parser,          only : get_param, param_file_type
use MOM_grid,                 only : ocean_grid_type
use time_interp_external_mod, only : init_external_field

implicit none ; private

!> Contains pointers to the forcing fields needed to drive MARBL
type, public :: marbl_forcing_type
  real, pointer, dimension(:,:) :: noy_dep => NULL() !> NOy Deposition [R Z T-1 ~> kgN m-2 s-1]
  real, pointer, dimension(:,:) :: nhx_dep => NULL() !> NHx Deposition [R Z T-1 ~> kgN m-2 s-1]
  real, pointer, dimension(:,:) :: dust_flux => NULL() !< Flux of dust into the ocean [m2 m-2]
  real, pointer, dimension(:,:) :: iron_flux => NULL() !< Flux of dust into the ocean [m2 m-2]
  real, pointer, dimension(:,:) :: ice_fraction => NULL() !< Fraction of ocean cell under seaice [m2 m-2]
  real, pointer, dimension(:,:) :: u10_sqr => NULL() !< 10m wind speed squared [L2 T-2 ~> m2 s-2]
end type marbl_forcing_type

type, public :: marbl_forcing_CS
  logical :: read_ndep                      !< If true, use nitrogen deposition supplied from an input file.
                                            !! This is temporary, we will always read NDEP
  character(len=200) :: ndep_file           !< If read_ndep, then this is the file from which to read
  real    :: dust_ratio_to_fe_bioavail_frac !< TODO: Add description
  real    :: fe_bioavail_frac_offset        !< TODO: Add description
  real    :: atm_fe_to_bc_ratio             !< TODO: Add description
  real    :: atm_bc_fe_bioavail_frac        !< TODO: Add description
  real    :: seaice_fe_to_bc_ratio          !< TODO: Add description
  real    :: seaice_bc_fe_bioavail_frac     !< TODO: Add description
  real    :: iron_frac_in_atm_fine_dust     !< Fraction of fine dust from the atmosphere that is iron
  real    :: iron_frac_in_atm_coarse_dust   !< Fraction of coarse dust from the atmosphere that is iron
  real    :: iron_frac_in_seaice_dust       !< Fraction of dust from the sea ice that is iron

  integer :: id_noydep   = -1     !< id number for time_interp_external.
  integer :: id_nhxdep   = -1     !< id number for time_interp_external.

end type marbl_forcing_CS

public :: marbl_forcing_init

contains

  subroutine marbl_forcing_init(G, param_file, inputdir, CS)
    type(ocean_grid_type),           intent(in)    :: G           !< The ocean's grid structure
    type(param_file_type),           intent(in)    :: param_file  !< A structure to parse for run-time parameters
    character(len=*),                intent(in)    :: inputdir    !< Directory containing input files
    type(marbl_forcing_CS), pointer, intent(inout) :: CS          !< A pointer that is set to point to control
                                                                  !! structure for MARBL forcing

    character(len=40)  :: mdl = "MOM_forcing_type"  ! This module's name.

    if (associated(CS)) then
      call MOM_error(WARNING, "marbl_forcing_init called with an associated "// &
                              "control structure.")
      return
    endif

    allocate(CS)
    call get_param(param_file, mdl, "DUST_RATIO_TO_FE_BIOAVAIL_FRAC", CS%dust_ratio_to_fe_bioavail_frac, &
    "TODO: Add description", default=1./170.)
    call get_param(param_file, mdl, "FE_BIOAVAIL_FRAC_OFFSET", CS%fe_bioavail_frac_offset, &
        "TODO: Add description", default=0.01)
    call get_param(param_file, mdl, "ATM_FE_TO_BC_RATIO", CS%atm_fe_to_bc_ratio, &
        "TODO: Add description", default=1.)
    call get_param(param_file, mdl, "ATM_BC_FE_BIOAVAIL_FRAC", CS%atm_bc_fe_bioavail_frac, &
        "TODO: Add description", default=0.06)
    call get_param(param_file, mdl, "SEAICE_FE_TO_BC_RATIO", CS%seaice_fe_to_bc_ratio, &
        "TODO: Add description", default=1.)
    call get_param(param_file, mdl, "SEAICE_BC_FE_BIOAVAIL_FRAC", CS%seaice_bc_fe_bioavail_frac, &
        "TODO: Add description", default=0.06)
    call get_param(param_file, mdl, "IRON_FRAC_IN_ATM_FINE_DUST", CS%iron_frac_in_atm_fine_dust, &
        "Fraction of fine dust from the atmosphere that is iron", default=0.035)
    call get_param(param_file, mdl, "IRON_FRAC_IN_ATM_COARSE_DUST", CS%iron_frac_in_atm_coarse_dust, &
        "Fraction of coarse dust from the atmosphere that is iron", default=0.035)
    call get_param(param_file, mdl, "IRON_FRAC_IN_SEAICE_DUST", CS%iron_frac_in_seaice_dust, &
        "Fraction of dust from sea ice that is iron", default=0.035)

    call get_param(param_file, mdl, "READ_NDEP", CS%read_ndep, &
        "If true, use nitrogen deposition supplied from "//&
        "an input file", default=.true.)
    if (CS%read_ndep) then
      ! TODO: we only want to read this variable in when running with MARBL
      call get_param(param_file, mdl, "NDEP_FILE", CS%ndep_file, &
            "The file in which the nitrogen deposition is found in "//&
            "variables NOy_deposition and NHx_deposition.", &
            default='ndep_ocn_1850_w_nhx_emis_MOM_tx0.66v1_c210222.nc')
      ! CS%ndep_file = trim(inputdir) // trim(CS%ndep_file)
      CS%ndep_file = trim('/glade/work/mlevy/cesm_inputdata/') // trim(CS%ndep_file)
      CS%id_noydep = init_external_field(CS%ndep_file, 'NDEP_NOy_month', domain=G%Domain%mpp_domain)
      CS%id_nhxdep = init_external_field(CS%ndep_file, 'NDEP_NHx_month', domain=G%Domain%mpp_domain)
    end if
  end subroutine marbl_forcing_init

end module marbl_forcing_type_main