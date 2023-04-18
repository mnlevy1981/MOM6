!> This module provides a common datatype to provide forcing for MARBL tracers
!! regardless of driver
module marbl_forcing_mod

!! This module exists to house code used by multiple drivers in config_src/
!! for passing forcing fields to MARBL
!! (This comment can go in the wiki on the NCAR fork?)

use MOM_diag_mediator,        only : safe_alloc_ptr, diag_ctrl, register_diag_field, post_data
use MOM_time_manager,         only : time_type
use MOM_error_handler,        only : MOM_error, WARNING
use MOM_file_parser,          only : get_param, log_param, param_file_type
use MOM_grid,                 only : ocean_grid_type
use MOM_unit_scaling,         only : unit_scale_type
use MOM_interpolate,          only : init_external_field, time_interp_external
use MOM_io,                   only : slasher
use tracer_forcing_utils_mod, only : forcing_timeseries_dataset
use tracer_forcing_utils_mod, only : forcing_timeseries_set_time_type_vars
use tracer_forcing_utils_mod, only : map_model_time_to_forcing_time
use marbl_constants_mod,      only : molw_Fe
use MOM_forcing_type,         only : forcing

implicit none ; private

#include <MOM_memory.h>

!> Data type used to store diagnostic index returned from register_diag_field()
!! For the forcing fields that can be written via post_data()
type, private :: marbl_forcing_diag_ids
  integer :: atm_fine_dust   !< Atmospheric fine dust component of dust_flux
  integer :: atm_coarse_dust !< Atmospheric coarse dust component of dust_flux
  integer :: atm_bc          !< Atmospheric black carbon component of iron_flux
  integer :: ice_dust        !< Sea-ice dust component of dust_flux
  integer :: ice_bc          !< Sea-ice black carbon component of iron_flux
  ! River fluxes
  integer :: no3_riv_flux          !< NO3 riverine flux
  integer :: po4_riv_flux          !< PO4 riverine flux
  integer :: don_riv_flux          !< DON riverine flux
  integer :: donr_riv_flux         !< DONr riverine flux
  integer :: dop_riv_flux          !< DOP riverine flux
  integer :: dopr_riv_flux         !< DOPr riverine flux
  integer :: sio3_riv_flux         !< SiO3 riverine flux
  integer :: fe_riv_flux           !< Fe riverine flux
  integer :: doc_riv_flux          !< DOC riverine flux
  integer :: docr_riv_flux         !< DOCr riverine flux
  integer :: alk_riv_flux          !< ALK riverine flux
  integer :: alk_alt_co2_riv_flux  !< ALK (alternate CO2) riverine flux
  integer :: dic_riv_flux          !< DIC riverine flux
  integer :: dic_alt_co2_riv_flux  !< DIC (alternate CO2) riverine flux
end type marbl_forcing_diag_ids

!> Control structure for this module
type, public :: marbl_forcing_CS
  logical :: read_riv_fluxes                !< If true, use river fluxes supplied from an input file.
                                            !! This is temporary, we will always read river fluxes
  type(forcing_timeseries_dataset) :: riv_flux_dataset !< File and time axis information for river fluxes
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.

  real    :: dust_ratio_thres               !< TODO: Add description
  real    :: dust_ratio_to_fe_bioavail_frac !< TODO: Add description
  real    :: fe_bioavail_frac_offset        !< TODO: Add description
  real    :: atm_fe_to_bc_ratio             !< TODO: Add description
  real    :: atm_bc_fe_bioavail_frac        !< TODO: Add description
  real    :: seaice_fe_to_bc_ratio          !< TODO: Add description
  real    :: seaice_bc_fe_bioavail_frac     !< TODO: Add description
  real    :: iron_frac_in_atm_fine_dust     !< Fraction of fine dust from the atmosphere that is iron
  real    :: iron_frac_in_atm_coarse_dust   !< Fraction of coarse dust from the atmosphere that is iron
  real    :: iron_frac_in_seaice_dust       !< Fraction of dust from the sea ice that is iron
  real    :: atm_co2_const                  !< atmospheric CO2 (if specifying a constant value) [ppm]
  real    :: atm_alt_co2_const              !< alternate atmospheric CO2 for _ALT_CO2 tracers
                                            !! (if specifying a constant value) [ppm]

  type(marbl_forcing_diag_ids) :: diag_ids  !< used for registering and posting some MARBL forcing fields as diagnostics


  integer :: id_din_riv  = -1     !< id number for time_interp_external.
  integer :: id_don_riv  = -1     !< id number for time_interp_external.
  integer :: id_dip_riv  = -1     !< id number for time_interp_external.
  integer :: id_dop_riv  = -1     !< id number for time_interp_external.
  integer :: id_dsi_riv  = -1     !< id number for time_interp_external.
  integer :: id_dfe_riv  = -1     !< id number for time_interp_external.
  integer :: id_dic_riv  = -1     !< id number for time_interp_external.
  integer :: id_alk_riv  = -1     !< id number for time_interp_external.
  integer :: id_doc_riv  = -1     !< id number for time_interp_external.

  logical :: use_marbl_tracers    !< most functions can return immediately
                                  !! MARBL tracers are turned off

end type marbl_forcing_CS

public :: marbl_forcing_init
public :: convert_marbl_IOB_to_forcings

contains

  subroutine marbl_forcing_init(G, param_file, diag, day, inputdir, use_marbl, CS)
    type(ocean_grid_type),           intent(in)    :: G           !< The ocean's grid structure
    type(param_file_type),           intent(in)    :: param_file  !< A structure to parse for run-time parameters
    type(diag_ctrl), target,         intent(in)    :: diag        !< Structure used to regulate diagnostic output.
    type(time_type), target,         intent(in)    :: day         !< Time of the start of the run.
    character(len=*),                intent(in)    :: inputdir    !< Directory containing input files
    logical,                         intent(in)    :: use_marbl   !< Is MARBL tracer package active?
    type(marbl_forcing_CS), pointer, intent(inout) :: CS          !< A pointer that is set to point to control
                                                                  !! structure for MARBL forcing

    character(len=40)  :: mdl = "MOM_forcing_type"  ! This module's name.
    character(len=200) :: inputdir2 ! The directory where the input files are.
    integer :: riv_flux_file_start_year
    integer :: riv_flux_file_end_year
    integer :: riv_flux_file_data_ref_year
    integer :: riv_flux_file_model_ref_year
    integer :: riv_flux_forcing_year

    if (associated(CS)) then
      call MOM_error(WARNING, "marbl_forcing_init called with an associated "// &
                              "control structure.")
      return
    endif

    allocate(CS)
    CS%diag => diag

    CS%use_marbl_tracers = .true.
    if (.not. use_marbl) then
      CS%use_marbl_tracers = .false.
      return
    end if

    ! TODO: just use DIN_LOC_ROOT
    call get_param(param_file, mdl, "CESM_INPUTDIR", inputdir2, default="/glade/work/mlevy/cesm_inputdata")

    call get_param(param_file, mdl, "DUST_RATIO_THRES", CS%dust_ratio_thres, &
    "TODO: Add description", default=69.00594)
    call get_param(param_file, mdl, "DUST_RATIO_TO_FE_BIOAVAIL_FRAC", CS%dust_ratio_to_fe_bioavail_frac, &
    "TODO: Add description", default=1./366.314)
    call get_param(param_file, mdl, "FE_BIOAVAIL_FRAC_OFFSET", CS%fe_bioavail_frac_offset, &
        "TODO: Add description", default=0.0146756)
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
    call get_param(param_file, mdl, "ATM_CO2_CONST", CS%atm_co2_const, &
        "Value to send to MARBL as xco2", &
        default=284.317, units="ppm")
    call get_param(param_file, mdl, "ATM_ALT_CO2_CONST", CS%atm_alt_co2_const, &
        "Value to send to MARBL as xco2_alt_co2", &
        default=284.317, units="ppm")

    ! ** River fluxes
    call get_param(param_file, mdl, "READ_RIV_FLUXES", CS%read_riv_fluxes, &
        "If true, use nitrogen deposition supplied from "//&
        "an input file", default=.true.)
    if (CS%read_riv_fluxes) then
      call get_param(param_file, mdl, "RIV_FLUX_FILE", CS%riv_flux_dataset%file_name, &
                    "The file in which the river fluxes can be found", &
                    default="riv_nut.gnews_gnm.JRA025m_to_tx0.66v1_nnsm_e333r100_190910.20210405.nc")
      ! call get_param(param_file, mdl, "RIV_FLUX_OFFSET_YEAR", CS%riv)
      if (scan(CS%riv_flux_dataset%file_name,'/') == 0) then
        ! CS%riv_flux_dataset%file_name = trim(inputdir) // trim(CS%riv_flux_dataset%file_name)
        CS%riv_flux_dataset%file_name = trim(slasher(inputdir2)) // trim(CS%riv_flux_dataset%file_name)
        call log_param(param_file, mdl, "INPUTDIR/RIV_FLUX_FILE", CS%riv_flux_dataset%file_name)
      end if
      call get_param(param_file, mdl, "RIV_FLUX_L_TIME_VARYING", CS%riv_flux_dataset%l_time_varying, &
                    ".true. for time-varying forcing, .false. for static forcing", default=.false.)
      if (CS%riv_flux_dataset%l_time_varying) then
        call get_param(param_file, mdl, "RIV_FLUX_FILE_START_YEAR", riv_flux_file_start_year, &
                      "Time coordinate of earliest date in RIV_FLUX_FILE", default=1900)
        call get_param(param_file, mdl, "RIV_FLUX_FILE_END_YEAR", riv_flux_file_end_year, &
                      "Time coordinate of earliest date in RIV_FLUX_FILE", default=1999)
        call get_param(param_file, mdl, "RIV_FLUX_FILE_DATA_REF_YEAR", riv_flux_file_data_ref_year, &
                      "Time coordinate of latest date in RIV_FLUX_FILE", default=1900)
        call get_param(param_file, mdl, "RIV_FLUX_FILE_MODEL_REF_YEAR", riv_flux_file_model_ref_year, &
                      "Time coordinate of latest date in RIV_FLUX_FILE",  default=1900)
      else
        call get_param(param_file, mdl, "RIV_FLUX_FORCING_YEAR", riv_flux_forcing_year, &
                      "Year from RIV_FLUX_FILE to use for forcing",  default=1900)
      end if
      call forcing_timeseries_set_time_type_vars(riv_flux_file_start_year, &
                                                riv_flux_file_end_year, &
                                                riv_flux_file_data_ref_year, &
                                                riv_flux_file_model_ref_year, &
                                                riv_flux_forcing_year, &
                                                CS%riv_flux_dataset)

      CS%id_din_riv = init_external_field(CS%riv_flux_dataset%file_name, 'din_riv_flux', domain=G%Domain%mpp_domain)
      CS%id_don_riv = init_external_field(CS%riv_flux_dataset%file_name, 'don_riv_flux', domain=G%Domain%mpp_domain)
      CS%id_dip_riv = init_external_field(CS%riv_flux_dataset%file_name, 'dip_riv_flux', domain=G%Domain%mpp_domain)
      CS%id_dop_riv = init_external_field(CS%riv_flux_dataset%file_name, 'dop_riv_flux', domain=G%Domain%mpp_domain)
      CS%id_dsi_riv = init_external_field(CS%riv_flux_dataset%file_name, 'dsi_riv_flux', domain=G%Domain%mpp_domain)
      CS%id_dfe_riv = init_external_field(CS%riv_flux_dataset%file_name, 'dfe_riv_flux', domain=G%Domain%mpp_domain)
      CS%id_dic_riv = init_external_field(CS%riv_flux_dataset%file_name, 'dic_riv_flux', domain=G%Domain%mpp_domain)
      CS%id_alk_riv = init_external_field(CS%riv_flux_dataset%file_name, 'alk_riv_flux', domain=G%Domain%mpp_domain)
      CS%id_doc_riv = init_external_field(CS%riv_flux_dataset%file_name, 'doc_riv_flux', domain=G%Domain%mpp_domain)
    endif

    ! Register diagnostic fields for outputing forcing values
    CS%diag_ids%atm_fine_dust = register_diag_field("ocean_model", &
                                                    "ATM_FINE_DUST_FLUX_CPL", &
                                                    CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                    day, &
                                                    "ATM_FINE_DUST_FLUX from cpl", &
                                                    "kg/m^2/s")
    CS%diag_ids%atm_coarse_dust = register_diag_field("ocean_model", &
                                                      "ATM_COARSE_DUST_FLUX_CPL", &
                                                      CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                      day, &
                                                      "ATM_COARSE_DUST_FLUX from cpl", &
                                                      "kg/m^2/s")
    CS%diag_ids%atm_bc = register_diag_field("ocean_model", &
                                             "ATM_BLACK_CARBON_FLUX_CPL", &
                                             CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                             day, &
                                             "ATM_BLACK_CARBON_FLUX from cpl", &
                                             "kg/m^2/s")

    CS%diag_ids%ice_dust = register_diag_field("ocean_model", &
                                               "SEAICE_DUST_FLUX_CPL", &
                                               CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                               day, &
                                               "SEAICE_DUST_FLUX from cpl", &
                                               "kg/m^2/s")
    CS%diag_ids%ice_bc = register_diag_field("ocean_model", &
                                             "SEAICE_BLACK_CARBON_FLUX_CPL", &
                                             CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                             day, &
                                             "SEAICE_BLACK_CARBON_FLUX from cpl", &
                                             "kg/m^2/s")

    CS%diag_ids%no3_riv_flux = register_diag_field("ocean_model", &
                                                   "NO3_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Inorganic Nitrate Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%po4_riv_flux = register_diag_field("ocean_model", &
                                                   "PO4_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Inorganic Phosphate Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%don_riv_flux = register_diag_field("ocean_model", &
                                                   "DON_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Organic Nitrogen Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%donr_riv_flux = register_diag_field("ocean_model", &
                                                    "DONR_RIV_FLUX", &
                                                    CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                    day, &
                                                    "Refractory DON Riverine Flux", &
                                                    "mmol/m^3 m/s")
    CS%diag_ids%dop_riv_flux = register_diag_field("ocean_model", &
                                                   "DOP_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Organic Phosphorus Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%dopr_riv_flux = register_diag_field("ocean_model", &
                                                    "DOPR_RIV_FLUX", &
                                                    CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                    day, &
                                                    "Refractory DOP Riverine Flux", &
                                                    "mmol/m^3 m/s")
    CS%diag_ids%sio3_riv_flux = register_diag_field("ocean_model", &
                                                    "SiO3_RIV_FLUX", &
                                                    CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                    day, &
                                                    "Dissolved Inorganic Silicate Riverine Flux", &
                                                    "mmol/m^3 m/s")
    CS%diag_ids%fe_riv_flux = register_diag_field("ocean_model", &
                                                  "Fe_RIV_FLUX", &
                                                  CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                  day, &
                                                  "Dissolved Inorganic Iron Riverine Flux", &
                                                  "mmol/m^3 m/s")
    CS%diag_ids%doc_riv_flux = register_diag_field("ocean_model", &
                                                   "DOC_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Organic Carbon Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%docr_riv_flux = register_diag_field("ocean_model", &
                                                    "DOCR_RIV_FLUX", &
                                                    CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                    day, &
                                                    "Refractory DOC Riverine Flux", &
                                                    "mmol/m^3 m/s")
    CS%diag_ids%alk_riv_flux = register_diag_field("ocean_model", &
                                                   "ALK_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Alkalinity Riverine Flux", &
                                                   "meq/m^3 m/s")
    CS%diag_ids%alk_alt_co2_riv_flux = register_diag_field("ocean_model", &
                                                          "ALK_ALT_CO2_RIV_FLUX", &
                                                          CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                          day, &
                                                          "Alkalinity Riverine Flux, Alternative CO2", &
                                                          "meq/m^3 m/s")
    CS%diag_ids%dic_riv_flux = register_diag_field("ocean_model", &
                                                   "DIC_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Inorganic Carbon Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%dic_alt_co2_riv_flux = register_diag_field("ocean_model", &
                                                          "DIC_ALT_CO2_RIV_FLUX", &
                                                          CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                          day, &
                                                          "Dissolved Inorganic Carbon Riverine Flux, Alternative CO2", &
                                                          "mmol/m^3 m/s")

  end subroutine marbl_forcing_init

  ! Note: ice fraction and u10_sqr are handled in mom_surface_forcing because of CFCs
  subroutine convert_marbl_IOB_to_forcings(atm_fine_dust_flux, atm_coarse_dust_flux, &
                                           seaice_dust_flux, atm_bc_flux, seaice_bc_flux, &
                                           nhx_dep, noy_dep, afracr, swnet_afracr, ifrac_n, &
                                           swpen_ifrac_n, Time, G, US, i0, j0, fluxes, CS)

    real, dimension(:,:),   pointer, intent(in)    :: atm_fine_dust_flux   !< atmosphere fine dust flux from IOB
    real, dimension(:,:),   pointer, intent(in)    :: atm_coarse_dust_flux !< atmosphere coarse dust flux from IOB
    real, dimension(:,:),   pointer, intent(in)    :: seaice_dust_flux     !< sea ice dust flux from IOB
    real, dimension(:,:),   pointer, intent(in)    :: atm_bc_flux          !< atmosphere black carbon flux from IOB
    real, dimension(:,:),   pointer, intent(in)    :: seaice_bc_flux       !< sea ice black carbon flux from IOB
    real, dimension(:,:),   pointer, intent(in)    :: afracr               !< open ocean fraction
    real, dimension(:,:),   pointer, intent(in)    :: nhx_dep              !< NHx flux from atmosphere
    real, dimension(:,:),   pointer, intent(in)    :: noy_dep              !< NOy flux from atmosphere
    real, dimension(:,:),   pointer, intent(in)    :: swnet_afracr         !< shortwave flux * open ocean fraction
    real, dimension(:,:,:), pointer, intent(in)    :: ifrac_n              !< per-category ice fraction
    real, dimension(:,:,:), pointer, intent(in)    :: swpen_ifrac_n        !< per-category shortwave flux * ice fraction
    type(time_type),                 intent(in)    :: Time                 !< The time of the fluxes, used for
                                                                           !! interpolating the salinity to the
                                                                           !! right time, when it is being
                                                                           !! restored.
    type(ocean_grid_type),           intent(in)    :: G                    !< The ocean's grid structure
    type(unit_scale_type),           intent(in)    :: US                   !< A dimensional unit scaling type
    integer,                         intent(in)    :: i0                   !< i index offset
    integer,                         intent(in)    :: j0                   !< j index offset
    type(forcing),                   intent(inout) :: fluxes               !< MARBL-specific forcing fields
    type(marbl_forcing_CS), pointer, intent(inout) :: CS                   !< A pointer that is set to point to
                                                                           !! control structure for MARBL forcing

    ! These are Fortran parameters in POP
    real, parameter :: DONriv_refract = 0.1
    real, parameter :: DOCriv_refract = 0.2
    real, parameter :: DOPriv_refract = 0.025

    real, dimension(SZI_(G),SZJ_(G)) :: time_varying_data  !< The field read in from forcing file with time dimension
    type(time_type) :: Time_riv_flux  !< For reading river flux fields, we use a modified version of Time
    integer :: i, j, is, ie, js, je, m
    real :: atm_fe_bioavail_frac     !< TODO: define this (local) term
    real :: seaice_fe_bioavail_frac  !< TODO: define this (local) term
    real :: dust_flux_conversion     !< TODO: define this (local) term
    real :: iron_flux_conversion     !< TODO: define this (local) term
    real :: ndep_conversion          !< Combination of unit conversion factors for rescaling
                                     !! nitrogen deposition [kg(N) m-2 s-1 ~> mol L-2 T-1]

    if (.not. CS%use_marbl_tracers) return

    is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
    ndep_conversion = (1000./14.) * ((US%L_to_m)**2 * US%T_to_s)
    dust_flux_conversion = US%kg_m2s_to_RZ_T * 0.1            ! kg / m^2 / s -> g / cm^2 / s
    iron_flux_conversion = US%kg_m2s_to_RZ_T * 1.e8 / molw_Fe ! kg / m^2 / s -> nmol / cm^2 / s

    ! Post fields from coupler to diagnostics
    ! TODO: units from diag register are incorrect; we should be converting these in the cap, I think
    time_varying_data(:,:) = 1.
    if (CS%diag_ids%atm_fine_dust > 0) &
      call post_data(CS%diag_ids%atm_fine_dust, &
                     US%kg_m2s_to_RZ_T * atm_fine_dust_flux(is-i0:ie-i0,js-j0:je-j0), &
                     CS%diag, mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%atm_coarse_dust > 0) &
      call post_data(CS%diag_ids%atm_coarse_dust, &
                     US%kg_m2s_to_RZ_T * atm_coarse_dust_flux(is-i0:ie-i0,js-j0:je-j0), &
                     CS%diag, mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%atm_bc > 0) &
      call post_data(CS%diag_ids%atm_bc, US%kg_m2s_to_RZ_T * atm_bc_flux(is-i0:ie-i0,js-j0:je-j0), &
                     CS%diag, mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%ice_dust > 0) &
      call post_data(CS%diag_ids%ice_dust, US%kg_m2s_to_RZ_T * seaice_dust_flux(is-i0:ie-i0,js-j0:je-j0), &
                     CS%diag, mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%ice_bc > 0) &
      call post_data(CS%diag_ids%ice_bc, US%kg_m2s_to_RZ_T * seaice_bc_flux(is-i0:ie-i0,js-j0:je-j0), &
                     CS%diag, mask=G%mask2dT(is:ie,js:je))

    do j=js,je ; do i=is,ie
      ! Nitrogen Deposition
      fluxes%nhx_dep(i,j) = (G%mask2dT(i,j) * ndep_conversion) * nhx_dep(i-i0,j-j0)
      fluxes%noy_dep(i,j) = (G%mask2dT(i,j) * ndep_conversion) * noy_dep(i-i0,j-j0)

      ! Atmospheric CO2
      fluxes%atm_co2(i,j) = G%mask2dT(i,j) * CS%atm_co2_const
      fluxes%atm_alt_co2(i,j) = G%mask2dT(i,j) * CS%atm_alt_co2_const

      if (associated(atm_fine_dust_flux)) then
        ! TODO: MARBL wants g/cm^2/s; we should convert to RZ_T in ocn_cap_methods then back to MARBL units here
        fluxes%dust_flux(i,j) = (G%mask2dT(i,j) * dust_flux_conversion) * &
                                       (atm_fine_dust_flux(i-i0,j-j0) + &
                                        atm_coarse_dust_flux(i-i0,j-j0) + &
                                        seaice_dust_flux(i-i0,j-j0))
      end if

      if (associated(atm_bc_flux)) then
        ! TODO: abort if atm_fine_dust_flux and atm_coarse_dust_flux are not associated?
        ! Contribution of atmospheric dust to iron flux
        if (atm_coarse_dust_flux(i-i0,j-j0) < &
            CS%dust_ratio_thres * atm_fine_dust_flux(i-i0,j-j0)) then
          atm_fe_bioavail_frac = CS%fe_bioavail_frac_offset + CS%dust_ratio_to_fe_bioavail_frac * &
            (CS%dust_ratio_thres - atm_coarse_dust_flux(i-i0,j-j0) / atm_fine_dust_flux(i-i0,j-j0))
        else
          atm_fe_bioavail_frac = CS%fe_bioavail_frac_offset
        end if
        ! Contribution of atmospheric dust to iron flux
        fluxes%iron_flux(i,j) = (atm_fe_bioavail_frac * &
                                 (CS%iron_frac_in_atm_fine_dust * atm_fine_dust_flux(i-i0,j-j0) + &
                                  CS%iron_frac_in_atm_coarse_dust * atm_coarse_dust_flux(i-i0,j-j0)))
        ! Contribution of atmospheric black carbon to iron flux
        fluxes%iron_flux(i,j) = fluxes%iron_flux(i,j) + (CS%atm_bc_fe_bioavail_frac * &
                                 (CS%atm_fe_to_bc_ratio * atm_bc_flux(i-i0,j-j0)))

        seaice_fe_bioavail_frac = atm_fe_bioavail_frac
        ! Contribution of seaice dust to iron flux
        fluxes%iron_flux(i,j) = fluxes%iron_flux(i,j) + (seaice_fe_bioavail_frac * &
                                 (CS%iron_frac_in_seaice_dust * seaice_dust_flux(i-i0,j-j0)))
        ! Contribution of seaice black carbon to iron flux
        fluxes%iron_flux(i,j) = fluxes%iron_flux(i,j) + (CS%seaice_bc_fe_bioavail_frac * &
                                 (CS%seaice_fe_to_bc_ratio * seaice_bc_flux(i-i0,j-j0)))

        ! Unit conversion (kg / m^2 / s -> nmol / cm^2 / s)
        fluxes%iron_flux(i,j) = (G%mask2dT(i,j) * iron_flux_conversion) * fluxes%iron_flux(i,j)

      end if

      ! Per ice-category forcings
      ! If the cap receives per-category fields, memory should be allocated in fluxes
      if (associated(ifrac_n)) then
        fluxes%fracr_cat(i,j,1) = min(1., afracr(i-i0,j-j0))
        fluxes%qsw_cat(i,j,1) = swnet_afracr(i-i0,j-j0)
        do m=1,size(ifrac_n, 3)
          fluxes%fracr_cat(i,j,m+1) = min(1., ifrac_n(i-i0,j-j0,m))
          fluxes%qsw_cat(i,j,m+1)   = swpen_ifrac_n(i-i0,j-j0,m)
        end do
        where (fluxes%fracr_cat(i,j,:) > 0.)
          fluxes%qsw_cat(i,j,:) = fluxes%qsw_cat(i,j,:) / fluxes%fracr_cat(i,j,:)
        elsewhere
          fluxes%fracr_cat(i,j,:) = 0.
          fluxes%qsw_cat(i,j,:) = 0.
        endwhere
        fluxes%fracr_cat(i,j,:) = G%mask2dT(i,j) * fluxes%fracr_cat(i,j,:)
        fluxes%qsw_cat(i,j,:)   = G%mask2dT(i,j) * fluxes%qsw_cat(i,j,:)
      endif

    enddo; enddo

    ! River fluxes
    if (CS%read_riv_fluxes) then
      fluxes%alk_riv_flux(:,:) = 0.
      fluxes%alk_alt_co2_riv_flux(:,:) = 0.

      ! DIN river flux affects NO3, ALK, and ALK_ALT_CO2
      Time_riv_flux = map_model_time_to_forcing_time(Time, CS%riv_flux_dataset)

      call time_interp_external(CS%id_din_riv,Time_riv_flux,time_varying_data)
      fluxes%no3_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)
      fluxes%alk_riv_flux(:,:) = fluxes%alk_riv_flux(:,:) - time_varying_data(:,:)
      fluxes%alk_alt_co2_riv_flux(:,:) = fluxes%alk_alt_co2_riv_flux(:,:) - time_varying_data(:,:)

      call time_interp_external(CS%id_dip_riv,Time_riv_flux,time_varying_data)
      fluxes%po4_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)

      call time_interp_external(CS%id_don_riv,Time_riv_flux,time_varying_data)
      fluxes%don_riv_flux(:,:) = G%mask2dT(:,:) * (1. - DONriv_refract) * &
                                        time_varying_data(:,:)
      fluxes%donr_riv_flux(:,:) = G%mask2dT(:,:) * DONriv_refract * &
                                        time_varying_data(:,:)

      call time_interp_external(CS%id_dop_riv,Time_riv_flux,time_varying_data)
      fluxes%dop_riv_flux(:,:) = G%mask2dT(:,:) * (1. - DOPriv_refract) * &
                                        time_varying_data(:,:)
      fluxes%dopr_riv_flux(:,:) = G%mask2dT(:,:) * DOPriv_refract * &
                                        time_varying_data(:,:)

      call time_interp_external(CS%id_dsi_riv,Time_riv_flux,time_varying_data)
      fluxes%sio3_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)

      call time_interp_external(CS%id_dfe_riv,Time_riv_flux,time_varying_data)
      fluxes%fe_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)

      call time_interp_external(CS%id_dic_riv,Time_riv_flux,time_varying_data)
      fluxes%dic_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)
      fluxes%dic_alt_co2_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)

      call time_interp_external(CS%id_alk_riv,Time_riv_flux,time_varying_data)
      fluxes%alk_riv_flux(:,:) = fluxes%alk_riv_flux(:,:) + time_varying_data(:,:)
      fluxes%alk_alt_co2_riv_flux(:,:) = fluxes%alk_alt_co2_riv_flux(:,:) + time_varying_data(:,:)

      call time_interp_external(CS%id_doc_riv,Time_riv_flux,time_varying_data)
      fluxes%doc_riv_flux(:,:) = G%mask2dT(:,:) * (1. - DOCriv_refract) * time_varying_data(:,:)
      fluxes%docr_riv_flux(:,:) = G%mask2dT(:,:) * DOCriv_refract * time_varying_data(:,:)
    else
      fluxes%no3_riv_flux(:,:) = 0.
      fluxes%po4_riv_flux(:,:) = 0.
      fluxes%don_riv_flux(:,:) = 0.
      fluxes%donr_riv_flux(:,:) = 0.
      fluxes%dop_riv_flux(:,:) = 0.
      fluxes%dopr_riv_flux(:,:) = 0.
      fluxes%sio3_riv_flux(:,:) = 0.
      fluxes%fe_riv_flux(:,:) = 0.
      fluxes%doc_riv_flux(:,:) = 0.
      fluxes%docr_riv_flux(:,:) = 0.
      fluxes%alk_riv_flux(:,:) = 0.
      fluxes%alk_alt_co2_riv_flux(:,:) = 0.
      fluxes%dic_riv_flux(:,:) = 0.
      fluxes%dic_alt_co2_riv_flux(:,:) = 0.
    end if

    ! Post to diags
    if (CS%diag_ids%no3_riv_flux > 0) &
      call post_data(CS%diag_ids%no3_riv_flux, fluxes%no3_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%po4_riv_flux > 0) &
      call post_data(CS%diag_ids%po4_riv_flux, fluxes%po4_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%don_riv_flux > 0) &
      call post_data(CS%diag_ids%don_riv_flux, fluxes%don_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%donr_riv_flux > 0) &
      call post_data(CS%diag_ids%donr_riv_flux, fluxes%donr_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%dop_riv_flux > 0) &
      call post_data(CS%diag_ids%dop_riv_flux, fluxes%dop_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%dopr_riv_flux > 0) &
      call post_data(CS%diag_ids%dopr_riv_flux, fluxes%dopr_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%sio3_riv_flux > 0) &
      call post_data(CS%diag_ids%sio3_riv_flux, fluxes%sio3_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%fe_riv_flux > 0) &
      call post_data(CS%diag_ids%fe_riv_flux, fluxes%fe_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%doc_riv_flux > 0) &
      call post_data(CS%diag_ids%doc_riv_flux, fluxes%doc_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%docr_riv_flux > 0) &
      call post_data(CS%diag_ids%docr_riv_flux, fluxes%docr_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%alk_riv_flux > 0) &
      call post_data(CS%diag_ids%alk_riv_flux, fluxes%alk_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%alk_alt_co2_riv_flux > 0) &
      call post_data(CS%diag_ids%alk_alt_co2_riv_flux, fluxes%alk_alt_co2_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%dic_riv_flux > 0) &
      call post_data(CS%diag_ids%dic_riv_flux, fluxes%dic_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%dic_alt_co2_riv_flux > 0) &
      call post_data(CS%diag_ids%dic_alt_co2_riv_flux, fluxes%dic_alt_co2_riv_flux(:,:), CS%diag)

  end subroutine convert_marbl_IOB_to_forcings

end module marbl_forcing_mod
