! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!!SKILLS: 0.3.2

#include "do_concurrent_compat.h"

!> Accelerations due to the Coriolis force and momentum advection
module MOM_CoriolisAdv

!> \author Robert Hallberg, April 1994 - June 2002

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : post_product_u, post_product_sum_u
use MOM_diag_mediator, only : post_product_v, post_product_sum_v
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_grid_containers, only : grid_core_type, grid_OBC_type, grid_core, grid_OBC
use MOM_open_boundary, only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary, only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_open_boundary, only : OBC_VORTICITY_ZERO, OBC_VORTICITY_FREESLIP
use MOM_open_boundary, only : OBC_VORTICITY_COMPUTED, OBC_VORTICITY_SPECIFIED
use MOM_string_functions, only : uppercase
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : accel_diag_ptrs, porous_barrier_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS
use box_mod,           only : Box_t, Box_C
use array_mod,         only : RealArray_t, RealArray_C
use iso_c_binding,     only : c_int, c_bool, c_double, c_null_char, c_ptr, c_loc, c_null_ptr
use posix,             only : mkdir_posix
use turbotmp_helperF,  only : getenv_mode, io_recorder, already_recorded, mark_recorded
use turbotmp_helperF,  only : TIMH_runAMREX, TIMH_capture, TIMH_runFORTRAN

implicit none ; private

public CorAdCalc, CoriolisAdv_init, CoriolisAdv_end, CoriolisAdv_stencil

#include <MOM_memory.h>

!> Options selecting the discretization of the Coriolis terms.
type, public :: Coriolis_scheme_CS
  integer :: Coriolis_Scheme !< Selects the discretization for the Coriolis terms.
                             !! Valid values are:
                             !! - SADOURNY75_ENERGY - Sadourny, 1975
                             !! - ARAKAWA_HSU90     - Arakawa & Hsu, 1990, Energy & non-div. Enstrophy
                             !! - ROBUST_ENSTRO     - Pseudo-enstrophy scheme
                             !! - SADOURNY75_ENSTRO - Sadourny, JAS 1975, Enstrophy
                             !! - ARAKAWA_LAMB81    - Arakawa & Lamb, MWR 1981, Energy & Enstrophy
                             !! - ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with Arakawa & Hsu and Sadourny energy.
                             !! - WENOVI3RD_PV_ENSTRO    - 3rd-order WENO scheme for PV reconstruction
                             !! - WENOVI5TH_PV_ENSTRO    - 5th-order WENO scheme for PV reconstruction
                             !! - WENOVI7TH_PV_ENSTRO    - 7th-order WENO scheme for PV reconstruction
                             !! The default, SADOURNY75_ENERGY, is the safest choice then the
                             !! deformation radius is poorly resolved.
  logical :: Coriolis_En_Dis !< If CORIOLIS_EN_DIS is defined, two estimates of
                             !! the thickness fluxes are used to estimate the
                             !! Coriolis term, and the one that dissipates energy
                             !! relative to the other one is used.  This is only
                             !! available at present if Coriolis scheme is
                             !! SADOURNY75_ENERGY.
end type Coriolis_scheme_CS

!> bind(C) mirror of Coriolis_scheme_CS, field-for-field, same order.
type, bind(C) :: Coriolis_scheme_CS_C
  integer(c_int)  :: Coriolis_Scheme !< Selects the discretization for the Coriolis terms.
  logical(c_bool) :: Coriolis_En_Dis !< If true, use the energy-dissipating biased Coriolis scheme.
end type Coriolis_scheme_CS_C

!> Options used in the Coriolis and momentum advection calculations.
type, public :: CorAdv_opts_CS
  type(Coriolis_scheme_CS) :: Coriolis_scheme_CS !< Options selecting the discretization of the
                             !! Coriolis terms.
  integer :: KE_Scheme       !< KE_SCHEME selects the discretization for
                             !! the kinetic energy. Valid values are:
                             !!  KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV
  logical :: KE_use_limiter  !< If true, use the Koren limiter for KE_UP3 scheme
  integer :: PV_Adv_Scheme   !< PV_ADV_SCHEME selects the discretization for PV advection
                             !! Valid values are:
                             !! - PV_ADV_CENTERED - centered (aka Sadourny, 75)
                             !! - PV_ADV_UPWIND1  - upwind, first order
  integer :: nkblock         !< The k block size used in Coriolis/advection calculations [nondim].
  real    :: F_eff_max_blend !< The factor by which the maximum effective Coriolis
                             !! acceleration from any point can be increased when
                             !! blending different discretizations with the
                             !! ARAKAWA_LAMB_BLEND Coriolis scheme [nondim].
                             !! This must be greater than 2.0, and is 4.0 by default.
  real    :: wt_lin_blend    !< A weighting value beyond which the blending between
                             !! Sadourny and Arakawa & Hsu goes linearly to 0 [nondim].
                             !! This must be between 1 and 1e-15, often 1/8.
  logical :: no_slip         !< If true, no slip boundary conditions are used.
                             !! Otherwise free slip boundary conditions are assumed.
                             !! The implementation of the free slip boundary
                             !! conditions on a C-grid is much cleaner than the
                             !! no slip boundary conditions. The use of free slip
                             !! b.c.s is strongly encouraged. The no slip b.c.s
                             !! are not implemented with the biharmonic viscosity.
  logical :: bound_Coriolis  !< If true, the Coriolis terms at u points are
                             !! bounded by the four estimates of (f+rv)v from the
                             !! four neighboring v points, and similarly at v
                             !! points.  This option would have no effect on the
                             !! SADOURNY75_ENERGY scheme if it were possible to
                             !! use centered difference thickness fluxes.
  logical :: weno_velocity_smooth !< If true, use velocity to compute the smoothness indicator for WENO
  logical :: do_CAS_diag = .false. !< If true, the Stokes contribution to CAu or CAv is needed
                                   !! for diagnostics.
  logical :: do_RV_diag = .false. !< If true, the relative vorticity is needed for diagnostics.
  logical :: do_PV_diag = .false. !< If true, the potential vorticity is needed for diagnostics.
end type CorAdv_opts_CS

!> bind(C) mirror of CorAdv_opts_CS, field-for-field, same order.
type, bind(C) :: CorAdv_opts_CS_C
  type(Coriolis_scheme_CS_C) :: Coriolis_scheme_CS !< Options selecting the discretization of the
                               !! Coriolis terms.
  integer(c_int) :: KE_Scheme !< KE_SCHEME selects the discretization for
                              !! the kinetic energy. Valid values are:
                              !!  KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV
  logical(c_bool) :: KE_use_limiter !< If true, use the Koren limiter for KE_UP3 scheme
  integer(c_int) :: PV_Adv_Scheme !< PV_ADV_SCHEME selects the discretization for PV advection
                                  !! Valid values are:
                                  !! - PV_ADV_CENTERED - centered (aka Sadourny, 75)
                                  !! - PV_ADV_UPWIND1  - upwind, first order
  integer(c_int) :: nkblock  !< The k block size used in Coriolis/advection calculations [nondim].
  real(c_double) :: F_eff_max_blend !< The factor by which the maximum effective Coriolis
                                    !! acceleration from any point can be increased when
                                    !! blending different discretizations with the
                                    !! ARAKAWA_LAMB_BLEND Coriolis scheme [nondim].
                                    !! This must be greater than 2.0, and is 4.0 by default.
  real(c_double) :: wt_lin_blend !< A weighting value beyond which the blending between
                                 !! Sadourny and Arakawa & Hsu goes linearly to 0 [nondim].
                                 !! This must be between 1 and 1e-15, often 1/8.
  logical(c_bool) :: no_slip !< If true, no slip boundary conditions are used.
                             !! Otherwise free slip boundary conditions are assumed.
                             !! The implementation of the free slip boundary
                             !! conditions on a C-grid is much cleaner than the
                             !! no slip boundary conditions. The use of free slip
                             !! b.c.s is strongly encouraged. The no slip b.c.s
                             !! are not implemented with the biharmonic viscosity.
  logical(c_bool) :: bound_Coriolis !< If true, the Coriolis terms at u points are
                                    !! bounded by the four estimates of (f+rv)v from the
                                    !! four neighboring v points, and similarly at v
                                    !! points.  This option would have no effect on the
                                    !! SADOURNY75_ENERGY scheme if it were possible to
                                    !! use centered difference thickness fluxes.
  logical(c_bool) :: weno_velocity_smooth !< If true, use velocity to compute the smoothness
                                          !! indicator for WENO
  logical(c_bool) :: do_CAS_diag !< If true, the Stokes contribution to CAu or CAv is needed
                                 !! for diagnostics.
  logical(c_bool) :: do_RV_diag !< If true, the relative vorticity is needed for diagnostics.
  logical(c_bool) :: do_PV_diag !< If true, the potential vorticity is needed for diagnostics.
end type CorAdv_opts_CS_C

!> Copies of the acceleration diagnostics that the Coriolis and momentum advection
!! calculations fill, each allocated only if its accel_diag_ptrs counterpart is associated.
type, public :: CorAdv_AD_type
  type(RealArray_t) :: gradKEu_a !< gradKEu = - d/dx(u2) [L T-2 ~> m s-2]
  type(RealArray_t) :: gradKEv_a !< gradKEv = - d/dy(u2) [L T-2 ~> m s-2]
  type(RealArray_t) :: rv_x_v_a  !< rv_x_v = rv * v at u [L T-2 ~> m s-2]
  type(RealArray_t) :: rv_x_u_a  !< rv_x_u = rv * u at v [L T-2 ~> m s-2]
contains
  procedure :: build_from => CorAdv_AD_build_from !< Copy in the associated diagnostics
  procedure :: copy_back => CorAdv_AD_copy_back   !< Copy out and free the diagnostics
end type CorAdv_AD_type

!> bind(C) mirror of CorAdv_AD_type, field-for-field, same order.
type, bind(C) :: CorAdv_AD_C
  type(RealArray_C) :: gradKEu !< gradKEu = - d/dx(u2) [L T-2 ~> m s-2]
  type(RealArray_C) :: gradKEv !< gradKEv = - d/dy(u2) [L T-2 ~> m s-2]
  type(RealArray_C) :: rv_x_v  !< rv_x_v = rv * v at u [L T-2 ~> m s-2]
  type(RealArray_C) :: rv_x_u  !< rv_x_u = rv * u at v [L T-2 ~> m s-2]
end type CorAdv_AD_C

!> Control structure for mom_coriolisadv
type, public :: CoriolisAdv_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  type(CorAdv_opts_CS) :: CorAdv_opts_CS !< Options used in the Coriolis and momentum advection
                                         !! calculations.
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the timing of diagnostic output.
  !>@{ Diagnostic IDs
  integer :: id_rv = -1, id_PV = -1, id_gKEu = -1, id_gKEv = -1
  integer :: id_rvxu = -1, id_rvxv = -1
  ! integer :: id_hf_gKEu    = -1, id_hf_gKEv    = -1
  integer :: id_hf_gKEu_2d = -1, id_hf_gKEv_2d = -1
  integer :: id_intz_gKEu_2d = -1, id_intz_gKEv_2d = -1
  ! integer :: id_hf_rvxu    = -1, id_hf_rvxv    = -1
  integer :: id_hf_rvxu_2d = -1, id_hf_rvxv_2d = -1
  integer :: id_h_gKEu = -1, id_h_gKEv = -1
  integer :: id_h_rvxu = -1, id_h_rvxv = -1
  integer :: id_intz_rvxu_2d = -1, id_intz_rvxv_2d = -1
  integer :: id_CAuS = -1, id_CAvS = -1
  !>@}
end type CoriolisAdv_CS

!>@{ Enumeration values for Coriolis_Scheme
integer, parameter :: SADOURNY75_ENERGY = 1
integer, parameter :: ARAKAWA_HSU90     = 2
integer, parameter :: ROBUST_ENSTRO     = 3
integer, parameter :: SADOURNY75_ENSTRO = 4
integer, parameter :: ARAKAWA_LAMB81    = 5
integer, parameter :: AL_BLEND          = 6
integer, parameter :: wenovi7th_PV_ENSTRO = 7
integer, parameter :: wenovi5th_PV_ENSTRO = 8
integer, parameter :: wenovi3rd_PV_ENSTRO = 9
character*(20), parameter :: SADOURNY75_ENERGY_STRING = "SADOURNY75_ENERGY"
character*(20), parameter :: ARAKAWA_HSU_STRING = "ARAKAWA_HSU90"
character*(20), parameter :: ROBUST_ENSTRO_STRING = "ROBUST_ENSTRO"
character*(20), parameter :: SADOURNY75_ENSTRO_STRING = "SADOURNY75_ENSTRO"
character*(20), parameter :: ARAKAWA_LAMB_STRING = "ARAKAWA_LAMB81"
character*(20), parameter :: AL_BLEND_STRING = "ARAKAWA_LAMB_BLEND"
character*(20), parameter :: WENOVI7TH_PV_ENSTRO_STRING = "WENOVI7TH_PV_ENSTRO"
character*(20), parameter :: WENOVI5TH_PV_ENSTRO_STRING = "WENOVI5TH_PV_ENSTRO"
character*(20), parameter :: WENOVI3RD_PV_ENSTRO_STRING = "WENOVI3RD_PV_ENSTRO"
!>@}
!>@{ Enumeration values for KE_Scheme
integer, parameter :: KE_ARAKAWA        = 10
integer, parameter :: KE_SIMPLE_GUDONOV = 11
integer, parameter :: KE_GUDONOV        = 12
integer, parameter :: KE_UP3            = 13
character*(20), parameter :: KE_ARAKAWA_STRING = "KE_ARAKAWA"
character*(20), parameter :: KE_SIMPLE_GUDONOV_STRING = "KE_SIMPLE_GUDONOV"
character*(20), parameter :: KE_GUDONOV_STRING = "KE_GUDONOV"
character*(20), parameter :: KE_UP3_STRING = "KE_UP3"
!>@}
!>@{ Enumeration values for PV_Adv_Scheme
integer, parameter :: PV_ADV_CENTERED   = 21
integer, parameter :: PV_ADV_UPWIND1    = 22
character*(20), parameter :: PV_ADV_CENTERED_STRING = "PV_ADV_CENTERED"
character*(20), parameter :: PV_ADV_UPWIND1_STRING = "PV_ADV_UPWIND1"
!>@}

  !----------------------------------------
  ! C interface (bridge to C++)
  !----------------------------------------
  interface
    !> Bridge for the gradKE subroutine
    subroutine turbotmp_gradke_bridge(bxQ, bxU, bxV, u, v, KE, KEx, KEy, Gcore, Gobc, KE_Scheme, &
                                      KE_use_limiter) bind(C)
      use iso_c_binding, only : c_int, c_bool
      use array_mod, only : RealArray_C
      use box_mod, only : Box_C
      use MOM_grid_containers, only : grid_core_C, grid_OBC_C
      implicit none
      type(Box_C),       intent(in)    :: bxQ   !< The B-grid-index iteration box of this tile
      type(Box_C),       intent(in)    :: bxU   !< The u-point iteration box of this tile
      type(Box_C),       intent(in)    :: bxV   !< The v-point iteration box of this tile
      type(RealArray_C), intent(in)    :: u     !< Zonal velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)    :: v     !< Meridional velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(inout) :: KE    !< Kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
      type(RealArray_C), intent(inout) :: KEx   !< Zonal acceleration due to kinetic
                                                !! energy gradient [L T-2 ~> m s-2]
      type(RealArray_C), intent(inout) :: KEy   !< Meridional acceleration due to kinetic
                                                !! energy gradient [L T-2 ~> m s-2]
      type(grid_core_C), intent(in)    :: Gcore !< Persistent copies of the core
                                                !! ocean grid fields
      type(grid_OBC_C),  intent(in)    :: Gobc  !< Persistent copies of the ocean
                                                !! grid OBC mask fields
      integer(c_int),  intent(in), value :: KE_Scheme !< Selects the discretization
                                                !! for the kinetic energy: KE_ARAKAWA,
                                                !! KE_SIMPLE_GUDONOV, KE_GUDONOV or KE_UP3
      logical(c_bool), intent(in), value :: KE_use_limiter !< If true, use the
                                                !! Koren limiter for the KE_UP3 scheme
    end subroutine turbotmp_gradke_bridge

    !> Bridge for the CorAdv_setup subroutine
    subroutine turbotmp_coradv_setup_bridge(bxH, bxQ, bxU, bxV, bxQs, u, v, h, OBC, AD, &
                                            por_face_areaU, por_face_areaV, Waves, Gcore, &
                                            Area_h, Area_q, Stokes_VF, use_weno, &
                                            vol_neglect, area_neglect, q, Ih_q, abs_vort, &
                                            h_q, qS, q2, uh_center, vh_center, RV, PV, CS) bind(C)
      use iso_c_binding, only : c_bool, c_double, c_ptr
      use array_mod, only : RealArray_C
      use box_mod, only : Box_C
      use MOM_grid_containers, only : grid_core_C
      import :: CorAdv_AD_C, CorAdv_opts_CS_C
      implicit none
      type(Box_C), intent(in)            :: bxH            !< The h-point iteration box of this
                                                           !! tile, [isc:iec, jsc:jec, ksc:kec]
      type(Box_C), intent(in)            :: bxQ            !< The B-grid-index iteration box of this
                                                           !! tile, [IscB:IecB, JscB:JecB, ksc:kec]
      type(Box_C), intent(in)            :: bxU            !< The u-point iteration box of this
                                                           !! tile, [IscB:IecB, jsc:jec, ksc:kec]
      type(Box_C), intent(in)            :: bxV            !< The v-point iteration box of this
                                                           !! tile, [isc:iec, JscB:JecB, ksc:kec]
      type(Box_C), intent(in)            :: bxQs           !< The scheme-dependent vorticity-point
                                                           !! box, [Is_q:Ie_q, Js_q:Je_q, ksc:kec]
      type(RealArray_C), intent(in)      :: u              !< Zonal velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)      :: v              !< Meridional velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)      :: h              !< Layer thickness [H ~> m or kg m-2]
      type(c_ptr), intent(in), value     :: OBC            !< Open boundary control structure
      type(CorAdv_AD_C), intent(in)      :: AD             !< Copies of the acceleration diagnostics
      type(RealArray_C), intent(in)      :: por_face_areaU !< fractional open area of U-faces
                                                           !! [nondim]
      type(RealArray_C), intent(in)      :: por_face_areaV !< fractional open area of V-faces
                                                           !! [nondim]
      type(c_ptr), intent(in), value     :: Waves          !< An optional pointer to Stokes drift CS
      type(grid_core_C), intent(in)      :: Gcore          !< Persistent copies of the core ocean
                                                           !! grid fields
      type(RealArray_C), intent(in)      :: Area_h         !< The ocean area at h points [L2 ~> m2]
      type(RealArray_C), intent(in)      :: Area_q         !< The sum of the ocean areas at the 4
                                                           !! adjacent thickness points [L2 ~> m2]
      logical(c_bool), intent(in), value :: Stokes_VF      !< If true, include the Stokes drift
      logical(c_bool), intent(in), value :: use_weno       !< True if using one of the WENO schemes
      real(c_double), intent(in), value  :: vol_neglect    !< A negligible volume [H L2 ~> m3 or kg]
      real(c_double), intent(in), value  :: area_neglect   !< A negligible area [L2 ~> m2]
      type(RealArray_C), intent(inout)   :: q              !< Layer potential vorticity
                                                           !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
      type(RealArray_C), intent(inout)   :: Ih_q           !< Inverse of thickness interpolated to q
                                                           !! points [H-1 ~> m-1 or m2 kg-1]
      type(RealArray_C), intent(inout)   :: abs_vort       !< Absolute vorticity at q-points
                                                           !! [T-1 ~> s-1]
      type(RealArray_C), intent(inout)   :: h_q            !< Thickness interpolated to q points
                                                           !! [H ~> m or kg m-2]
      type(RealArray_C), intent(inout)   :: qS             !< Layer Stokes vorticity
                                                           !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
      type(RealArray_C), intent(inout)   :: q2             !< Relative vorticity over thickness
                                                           !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
      type(RealArray_C), intent(inout)   :: uh_center      !< Transport based on arithmetic mean h
                                                           !! at u-points
                                                           !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(inout)   :: vh_center      !< Transport based on arithmetic mean h
                                                           !! at v-points
                                                           !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(inout)   :: RV             !< Diagnostic relative vorticity
                                                           !! [T-1 ~> s-1]
      type(RealArray_C), intent(inout)   :: PV             !< Diagnostic potential vorticity
                                                           !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
      type(CorAdv_opts_CS_C), intent(in) :: CS             !< Options used in the Coriolis and
                                                           !! momentum advection calculations
    end subroutine turbotmp_coradv_setup_bridge

    !> Bridge for the CorAdv_sadourny subroutine
    subroutine turbotmp_coradv_sadourny_bridge(bxH, bxU, bxV, u, v, uh, vh, q, uh_center, &
                                               vh_center, Gcore, CAu, CAv, CS) bind(C)
      use array_mod, only : RealArray_C
      use box_mod, only : Box_C
      use MOM_grid_containers, only : grid_core_C
      import :: Coriolis_scheme_CS_C
      implicit none
      type(Box_C), intent(in)                :: bxH       !< The h-point iteration box of this tile,
                                                          !! [isc:iec, jsc:jec, ksc:kec]
      type(Box_C), intent(in)                :: bxU       !< The u-point iteration box of this tile,
                                                          !! [IscB:IecB, jsc:jec, ksc:kec]
      type(Box_C), intent(in)                :: bxV       !< The v-point iteration box of this tile,
                                                          !! [isc:iec, JscB:JecB, ksc:kec]
      type(RealArray_C), intent(in)          :: u         !< Zonal velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)          :: v         !< Meridional velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)          :: uh        !< Zonal transport u*h*dy
                                                          !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)          :: vh        !< Meridional transport v*h*dx
                                                          !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)          :: q         !< Layer potential vorticity
                                                          !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
      type(RealArray_C), intent(in)          :: uh_center !< Transport based on arithmetic mean h at
                                                          !! u-points [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)          :: vh_center !< Transport based on arithmetic mean h at
                                                          !! v-points [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(grid_core_C), intent(in)          :: Gcore     !< Persistent copies of the core ocean
                                                          !! grid fields
      type(RealArray_C), intent(inout)       :: CAu       !< Zonal acceleration due to Coriolis and
                                                          !! momentum advection [L T-2 ~> m s-2].
      type(RealArray_C), intent(inout)       :: CAv       !< Meridional acceleration due to Coriolis
                                                          !! and momentum advection
                                                          !! [L T-2 ~> m s-2].
      type(Coriolis_scheme_CS_C), intent(in) :: CS        !< Options selecting the Coriolis
                                                          !! discretization
    end subroutine turbotmp_coradv_sadourny_bridge

    !> Bridge for the CorAdv_arakawa subroutine
    subroutine turbotmp_coradv_arakawa_bridge(bxH, bxQ, bxU, bxV, uh, vh, q, Ih_q, Fe_m2, &
                                              rat_lin, Gcore, CAu, CAv, wt_lin_blend, CS) bind(C)
      use iso_c_binding, only : c_double
      use array_mod, only : RealArray_C
      use box_mod, only : Box_C
      use MOM_grid_containers, only : grid_core_C
      import :: Coriolis_scheme_CS_C
      implicit none
      type(Box_C), intent(in)                :: bxH          !< The h-point iteration box of this
                                                             !! tile, [isc:iec, jsc:jec, ksc:kec]
      type(Box_C), intent(in)                :: bxQ          !< The B-grid-index iteration box of
                                                             !! this tile,
                                                             !! [IscB:IecB, JscB:JecB, ksc:kec]
      type(Box_C), intent(in)                :: bxU          !< The u-point iteration box of this
                                                             !! tile, [IscB:IecB, jsc:jec, ksc:kec]
      type(Box_C), intent(in)                :: bxV          !< The v-point iteration box of this
                                                             !! tile, [isc:iec, JscB:JecB, ksc:kec]
      type(RealArray_C), intent(in)          :: uh           !< Zonal transport u*h*dy
                                                             !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)          :: vh           !< Meridional transport v*h*dx
                                                             !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)          :: q            !< Layer potential vorticity
                                                             !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
      type(RealArray_C), intent(in)          :: Ih_q         !< Inverse of thickness interpolated to
                                                             !! q points [H-1 ~> m-1 or m2 kg-1]
      real(c_double), intent(in), value      :: Fe_m2        !< ARAKAWA_LAMB_BLEND blending
                                                             !! parameter [nondim]
      real(c_double), intent(in), value      :: rat_lin      !< ARAKAWA_LAMB_BLEND blending
                                                             !! parameter [nondim]
      type(grid_core_C), intent(in)          :: Gcore        !< Persistent copies of the core ocean
                                                             !! grid fields
      type(RealArray_C), intent(inout)       :: CAu          !< Zonal acceleration due to Coriolis
                                                             !! and momentum advection
                                                             !! [L T-2 ~> m s-2].
      type(RealArray_C), intent(inout)       :: CAv          !< Meridional acceleration due to
                                                             !! Coriolis and momentum advection
                                                             !! [L T-2 ~> m s-2].
      real(c_double), intent(in), value      :: wt_lin_blend !< A weighting value beyond which the
                                                             !! blending between Sadourny and
                                                             !! Arakawa & Hsu goes linearly to 0
                                                             !! [nondim].
      type(Coriolis_scheme_CS_C), intent(in) :: CS           !< Options selecting the Coriolis
                                                             !! discretization
    end subroutine turbotmp_coradv_arakawa_bridge

    !> Bridge for the CorAdv_robust_enstro subroutine
    subroutine turbotmp_coradv_robust_enstro_bridge(bxU, bxV, u, v, h, uh, vh, abs_vort, &
                                                    eps_vel, h_tiny, Gcore, CAu, CAv, &
                                                    PV_Adv_Scheme) bind(C)
      use iso_c_binding, only : c_double, c_int
      use array_mod, only : RealArray_C
      use box_mod, only : Box_C
      use MOM_grid_containers, only : grid_core_C
      implicit none
      type(Box_C), intent(in)           :: bxU           !< The u-point iteration box of this tile,
                                                         !! [IscB:IecB, jsc:jec, ksc:kec]
      type(Box_C), intent(in)           :: bxV           !< The v-point iteration box of this tile,
                                                         !! [isc:iec, JscB:JecB, ksc:kec]
      type(RealArray_C), intent(in)     :: u             !< Zonal velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)     :: v             !< Meridional velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)     :: h             !< Layer thickness [H ~> m or kg m-2]
      type(RealArray_C), intent(in)     :: uh            !< Zonal transport u*h*dy
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)     :: vh            !< Meridional transport v*h*dx
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)     :: abs_vort      !< Absolute vorticity at q-points
                                                         !! [T-1 ~> s-1]
      real(c_double), intent(in), value :: eps_vel       !< A tiny, positive velocity
                                                         !! [L T-1 ~> m s-1]
      real(c_double), intent(in), value :: h_tiny        !< A very small thickness
                                                         !! [H ~> m or kg m-2]
      type(grid_core_C), intent(in)     :: Gcore         !< Persistent copies of the core ocean grid
                                                         !! fields
      type(RealArray_C), intent(inout)  :: CAu           !< Zonal acceleration due to Coriolis and
                                                         !! momentum advection [L T-2 ~> m s-2].
      type(RealArray_C), intent(inout)  :: CAv           !< Meridional acceleration due to Coriolis
                                                         !! and momentum advection [L T-2 ~> m s-2].
      integer(c_int), intent(in), value :: PV_Adv_Scheme !< Selects the discretization for PV
                                                         !! advection: PV_ADV_CENTERED or
                                                         !! PV_ADV_UPWIND1
    end subroutine turbotmp_coradv_robust_enstro_bridge

    !> Bridge for the CorAdv_weno subroutine
    subroutine turbotmp_coradv_weno_bridge(bxU, bxV, u, v, uh, vh, q, abs_vort, h_q, &
                                           Gcore, H_subroundoff, CAu, CAv, &
                                           weno_velocity_smooth, CS) bind(C)
      use iso_c_binding, only : c_bool, c_double
      use array_mod, only : RealArray_C
      use box_mod, only : Box_C
      use MOM_grid_containers, only : grid_core_C
      import :: Coriolis_scheme_CS_C
      implicit none
      type(Box_C), intent(in)                :: bxU           !< The u-point iteration box of this
                                                              !! tile, [IscB:IecB, jsc:jec, ksc:kec]
      type(Box_C), intent(in)                :: bxV           !< The v-point iteration box of this
                                                              !! tile, [isc:iec, JscB:JecB, ksc:kec]
      type(RealArray_C), intent(in)          :: u             !< Zonal velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)          :: v             !< Meridional velocity
                                                              !! [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)          :: uh            !< Zonal transport u*h*dy
                                                              !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)          :: vh            !< Meridional transport v*h*dx
                                                              !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)          :: q             !< Layer potential vorticity
                                                              !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
      type(RealArray_C), intent(in)          :: abs_vort      !< Absolute vorticity at q-points
                                                              !! [T-1 ~> s-1]
      type(RealArray_C), intent(in)          :: h_q           !< Thickness interpolated to q points
                                                              !! [H ~> m or kg m-2]
      type(grid_core_C), intent(in)          :: Gcore         !< Persistent copies of the core ocean
                                                              !! grid fields
      real(c_double), intent(in), value      :: H_subroundoff !< A thickness that is so small that
                                                              !! it can be added to a thickness of
                                                              !! Angstrom or larger without changing
                                                              !! it at the bit level
                                                              !! [H ~> m or kg m-2].
      type(RealArray_C), intent(inout)       :: CAu           !< Zonal acceleration due to Coriolis
                                                              !! and momentum advection
                                                              !! [L T-2 ~> m s-2].
      type(RealArray_C), intent(inout)       :: CAv           !< Meridional acceleration due to
                                                              !! Coriolis and momentum advection
                                                              !! [L T-2 ~> m s-2].
      logical(c_bool), intent(in), value     :: weno_velocity_smooth !< If true, use velocity to
                                                                     !! compute the smoothness
                                                                     !! indicator for WENO
      type(Coriolis_scheme_CS_C), intent(in) :: CS            !< Options selecting the Coriolis
                                                              !! discretization
    end subroutine turbotmp_coradv_weno_bridge

    !> Bridge for the CorAdv_common_terms subroutine
    subroutine turbotmp_coradv_common_terms_bridge(bxU, bxV, u, v, uh, vh, abs_vort, qS, &
                                                   q2, KEx, KEy, Stokes_VF, Gcore, CAu, &
                                                   CAv, CAuS, CAvS, AD, CS) bind(C)
      use iso_c_binding, only : c_bool
      use array_mod, only : RealArray_C
      use box_mod, only : Box_C
      use MOM_grid_containers, only : grid_core_C
      import :: CorAdv_AD_C, CorAdv_opts_CS_C
      implicit none
      type(Box_C), intent(in)            :: bxU       !< The u-point iteration box of this tile,
                                                      !! [IscB:IecB, jsc:jec, ksc:kec]
      type(Box_C), intent(in)            :: bxV       !< The v-point iteration box of this tile,
                                                      !! [isc:iec, JscB:JecB, ksc:kec]
      type(RealArray_C), intent(in)      :: u         !< Zonal velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)      :: v         !< Meridional velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)      :: uh        !< Zonal transport u*h*dy
                                                      !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)      :: vh        !< Meridional transport v*h*dx
                                                      !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)      :: abs_vort  !< Absolute vorticity at q-points [T-1 ~> s-1]
      type(RealArray_C), intent(in)      :: qS        !< Layer Stokes vorticity
                                                      !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
      type(RealArray_C), intent(in)      :: q2        !< Relative vorticity over thickness
                                                      !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
      type(RealArray_C), intent(in)      :: KEx       !< Zonal gradient of kinetic energy per unit
                                                      !! mass [L T-2 ~> m s-2]
      type(RealArray_C), intent(in)      :: KEy       !< Meridional gradient of kinetic energy per
                                                      !! unit mass [L T-2 ~> m s-2]
      logical(c_bool), intent(in), value :: Stokes_VF !< If true, include the Stokes drift
      type(grid_core_C), intent(in)      :: Gcore     !< Persistent copies of the core ocean grid
                                                      !! fields
      type(RealArray_C), intent(inout)   :: CAu       !< Zonal acceleration due to Coriolis and
                                                      !! momentum advection [L T-2 ~> m s-2].
      type(RealArray_C), intent(inout)   :: CAv       !< Meridional acceleration due to Coriolis and
                                                      !! momentum advection [L T-2 ~> m s-2].
      type(RealArray_C), intent(inout)   :: CAuS      !< Stokes contribution to CAu [L T-2 ~> m s-2]
      type(RealArray_C), intent(inout)   :: CAvS      !< Stokes contribution to CAv [L T-2 ~> m s-2]
      type(CorAdv_AD_C), intent(inout)   :: AD        !< Copies of the acceleration diagnostics
      type(CorAdv_opts_CS_C), intent(in) :: CS        !< Options used in the Coriolis and momentum
                                                      !! advection calculations
    end subroutine turbotmp_coradv_common_terms_bridge

    !> Bridge for the CorAdCalc_TR subroutine
    subroutine turbotmp_coradcalc_tr_bridge(bxH0, bxQ0, u, v, h, uh, vh, CAu, CAv, RV, PV, &
                                            CAuS, CAvS, OBC, AD, Gcore, Gobc, nz, &
                                            H_subroundoff, Angstrom_H, m_to_L, m_s_to_L_T, &
                                            CS, por_face_areaU, por_face_areaV, Waves, &
                                            Stokes_VF) bind(C)
      use iso_c_binding, only : c_bool, c_double, c_int, c_ptr
      use array_mod, only : RealArray_C
      use box_mod, only : Box_C
      use MOM_grid_containers, only : grid_core_C, grid_OBC_C
      import :: CorAdv_AD_C, CorAdv_opts_CS_C
      implicit none
      type(Box_C), intent(in)            :: bxH0           !< The h-point box of the whole column,
                                                           !! [isc:iec, jsc:jec, 1:nz]
      type(Box_C), intent(in)            :: bxQ0           !< The B-grid-index box of the whole
                                                           !! column, [IscB:IecB, JscB:JecB, 1:nz]
      type(RealArray_C), intent(in)      :: u              !< Zonal velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)      :: v              !< Meridional velocity [L T-1 ~> m s-1]
      type(RealArray_C), intent(in)      :: h              !< Layer thickness [H ~> m or kg m-2]
      type(RealArray_C), intent(in)      :: uh             !< Zonal transport u*h*dy
                                                           !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(in)      :: vh             !< Meridional transport v*h*dx
                                                           !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      type(RealArray_C), intent(inout)   :: CAu            !< Zonal acceleration due to Coriolis and
                                                           !! momentum advection [L T-2 ~> m s-2].
      type(RealArray_C), intent(inout)   :: CAv            !< Meridional acceleration due to
                                                           !! Coriolis and momentum advection
                                                           !! [L T-2 ~> m s-2].
      type(RealArray_C), intent(inout)   :: RV             !< Diagnostic relative vorticity
                                                           !! [T-1 ~> s-1]
      type(RealArray_C), intent(inout)   :: PV             !< Diagnostic potential vorticity
                                                           !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
      type(RealArray_C), intent(inout)   :: CAuS           !< Stokes contribution to CAu
                                                           !! [L T-2 ~> m s-2]
      type(RealArray_C), intent(inout)   :: CAvS           !< Stokes contribution to CAv
                                                           !! [L T-2 ~> m s-2]
      type(c_ptr), intent(in), value     :: OBC            !< Open boundary control structure
      type(CorAdv_AD_C), intent(inout)   :: AD             !< Copies of the acceleration diagnostics
      type(grid_core_C), intent(in)      :: Gcore          !< Persistent copies of the core ocean
                                                           !! grid fields
      type(grid_OBC_C), intent(in)       :: Gobc           !< Persistent copies of the ocean grid
                                                           !! OBC mask fields
      integer(c_int), intent(in), value  :: nz             !< The number of layers/levels in the
                                                           !! vertical
      real(c_double), intent(in), value  :: H_subroundoff  !< A thickness that is so small that it
                                                           !! can be added to a thickness of
                                                           !! Angstrom or larger without changing it
                                                           !! at the bit level [H ~> m or kg m-2].
      real(c_double), intent(in), value  :: Angstrom_H     !< A one-Angstrom thickness in the model
                                                           !! thickness units [H ~> m or kg m-2].
      real(c_double), intent(in), value  :: m_to_L         !< A constant that translates lengths in
                                                           !! meters to the units of horizontal
                                                           !! lengths [L m-1 ~> 1]
      real(c_double), intent(in), value  :: m_s_to_L_T     !< Convert lateral velocities from m s-1
                                                           !! to L T-1 [L s T-1 m-1 ~> 1]
      type(CorAdv_opts_CS_C), intent(in) :: CS             !< Options used in the Coriolis and
                                                           !! momentum advection calculations
      type(RealArray_C), intent(in)      :: por_face_areaU !< fractional open area of U-faces
                                                           !! [nondim]
      type(RealArray_C), intent(in)      :: por_face_areaV !< fractional open area of V-faces
                                                           !! [nondim]
      type(c_ptr), intent(in), value     :: Waves          !< An optional pointer to Stokes drift CS
      logical(c_bool), intent(in), value :: Stokes_VF      !< If true, include the Stokes drift
    end subroutine turbotmp_coradcalc_tr_bridge
  end interface

contains

!> Converts a Coriolis_scheme_CS to its bind(C) mirror.
function Coriolis_scheme_CS_to_c(opts) result(cdesc)
  type(Coriolis_scheme_CS), intent(in) :: opts !< Options to convert
  type(Coriolis_scheme_CS_C) :: cdesc          !< bind(C) mirror of opts
  cdesc%Coriolis_Scheme = opts%Coriolis_Scheme
  cdesc%Coriolis_En_Dis = opts%Coriolis_En_Dis
end function Coriolis_scheme_CS_to_c

!> Converts a CorAdv_opts_CS to its bind(C) mirror.
function CorAdv_opts_CS_to_c(opts) result(cdesc)
  type(CorAdv_opts_CS), intent(in) :: opts !< Options to convert
  type(CorAdv_opts_CS_C) :: cdesc          !< bind(C) mirror of opts
  cdesc%Coriolis_scheme_CS = Coriolis_scheme_CS_to_c(opts%Coriolis_scheme_CS)
  cdesc%KE_Scheme = opts%KE_Scheme
  cdesc%KE_use_limiter = opts%KE_use_limiter
  cdesc%PV_Adv_Scheme = opts%PV_Adv_Scheme
  cdesc%nkblock = opts%nkblock
  cdesc%F_eff_max_blend = opts%F_eff_max_blend
  cdesc%wt_lin_blend = opts%wt_lin_blend
  cdesc%no_slip = opts%no_slip
  cdesc%bound_Coriolis = opts%bound_Coriolis
  cdesc%weno_velocity_smooth = opts%weno_velocity_smooth
  cdesc%do_CAS_diag = opts%do_CAS_diag
  cdesc%do_RV_diag = opts%do_RV_diag
  cdesc%do_PV_diag = opts%do_PV_diag
end function CorAdv_opts_CS_to_c

!> Copies each associated acceleration diagnostic in AD into this container set.
subroutine CorAdv_AD_build_from(this, AD)
  class(CorAdv_AD_type),  intent(inout) :: this !< Copies of the acceleration diagnostics
  type(accel_diag_ptrs),  intent(in)    :: AD   !< Storage for acceleration diagnostics

  if (associated(AD%gradKEu)) &
    call this%gradKEu_a%alloc(lb=LBOUND(AD%gradKEu), ub=UBOUND(AD%gradKEu), source=AD%gradKEu)
  if (associated(AD%gradKEv)) &
    call this%gradKEv_a%alloc(lb=LBOUND(AD%gradKEv), ub=UBOUND(AD%gradKEv), source=AD%gradKEv)
  if (associated(AD%rv_x_v)) &
    call this%rv_x_v_a%alloc(lb=LBOUND(AD%rv_x_v), ub=UBOUND(AD%rv_x_v), source=AD%rv_x_v)
  if (associated(AD%rv_x_u)) &
    call this%rv_x_u_a%alloc(lb=LBOUND(AD%rv_x_u), ub=UBOUND(AD%rv_x_u), source=AD%rv_x_u)
end subroutine CorAdv_AD_build_from

!> Copies each allocated acceleration diagnostic back into AD and frees it.
subroutine CorAdv_AD_copy_back(this, AD)
  class(CorAdv_AD_type),  intent(inout) :: this !< Copies of the acceleration diagnostics
  type(accel_diag_ptrs),  intent(inout) :: AD   !< Storage for acceleration diagnostics

  if (this%gradKEu_a%associated()) then
    call this%gradKEu_a%copy2F(AD%gradKEu) ; call this%gradKEu_a%free()
  endif
  if (this%gradKEv_a%associated()) then
    call this%gradKEv_a%copy2F(AD%gradKEv) ; call this%gradKEv_a%free()
  endif
  if (this%rv_x_v_a%associated()) then
    call this%rv_x_v_a%copy2F(AD%rv_x_v) ; call this%rv_x_v_a%free()
  endif
  if (this%rv_x_u_a%associated()) then
    call this%rv_x_u_a%copy2F(AD%rv_x_u) ; call this%rv_x_u_a%free()
  endif
end subroutine CorAdv_AD_copy_back

!> Converts a CorAdv_AD_type to its bind(C) mirror.
#ifdef _TIM
function CorAdv_AD_to_c(AD_a) result(cdesc)
  type(CorAdv_AD_type), intent(in) :: AD_a !< Copies of the acceleration diagnostics
  type(CorAdv_AD_C) :: cdesc               !< bind(C) mirror of AD_a
  cdesc%gradKEu = AD_a%gradKEu_a%to_c()
  cdesc%gradKEv = AD_a%gradKEv_a%to_c()
  cdesc%rv_x_v = AD_a%rv_x_v_a%to_c()
  cdesc%rv_x_u = AD_a%rv_x_u_a%to_c()
end function CorAdv_AD_to_c
#endif

!> Records the members of a Coriolis_scheme_CS in a capture file, with names starting with
!! prefix.
subroutine record_Coriolis_scheme_CS(rec, prefix, opts)
  type(io_recorder),        intent(inout) :: rec    !< The capture file being written
  character(len=*),         intent(in)    :: prefix !< The start of each record name
  type(Coriolis_scheme_CS), intent(in)    :: opts   !< The options to record

  call rec%add(prefix//"_Coriolis_Scheme", opts%Coriolis_Scheme)
  call rec%add(prefix//"_Coriolis_En_Dis", opts%Coriolis_En_Dis)
end subroutine record_Coriolis_scheme_CS

!> Records the members of a CorAdv_opts_CS in a capture file, with names starting with prefix.
subroutine record_CorAdv_opts_CS(rec, prefix, opts)
  type(io_recorder),     intent(inout) :: rec    !< The capture file being written
  character(len=*),      intent(in)    :: prefix !< The start of each record name
  type(CorAdv_opts_CS),  intent(in)    :: opts   !< The options to record

  call record_Coriolis_scheme_CS(rec, prefix//"_Coriolis_scheme_CS", opts%Coriolis_scheme_CS)
  call rec%add(prefix//"_KE_Scheme", opts%KE_Scheme)
  call rec%add(prefix//"_KE_use_limiter", opts%KE_use_limiter)
  call rec%add(prefix//"_PV_Adv_Scheme", opts%PV_Adv_Scheme)
  call rec%add(prefix//"_nkblock", opts%nkblock)
  call rec%add(prefix//"_F_eff_max_blend", opts%F_eff_max_blend)
  call rec%add(prefix//"_wt_lin_blend", opts%wt_lin_blend)
  call rec%add(prefix//"_no_slip", opts%no_slip)
  call rec%add(prefix//"_bound_Coriolis", opts%bound_Coriolis)
  call rec%add(prefix//"_weno_velocity_smooth", opts%weno_velocity_smooth)
  call rec%add(prefix//"_do_CAS_diag", opts%do_CAS_diag)
  call rec%add(prefix//"_do_RV_diag", opts%do_RV_diag)
  call rec%add(prefix//"_do_PV_diag", opts%do_PV_diag)
end subroutine record_CorAdv_opts_CS

!> Records which acceleration diagnostics in a CorAdv_AD_type are allocated (with an empty
!! suffix), or the values of the allocated ones (with a non-empty suffix), in a capture file.
subroutine record_CorAdv_AD(rec, prefix, AD_a, suffix)
  type(io_recorder),     intent(inout) :: rec    !< The capture file being written
  character(len=*),      intent(in)    :: prefix !< The start of each record name
  type(CorAdv_AD_type),  intent(in)    :: AD_a   !< Copies of the acceleration diagnostics
  character(len=*),      intent(in)    :: suffix !< The end of each value record name, or empty
                                                 !! to record which diagnostics are allocated

  if (len(suffix) == 0) then
    call rec%add(prefix//"_gradKEu_allocated", AD_a%gradKEu_a%associated())
  elseif (AD_a%gradKEu_a%associated()) then
    call rec%add(prefix//"_gradKEu"//suffix, AD_a%gradKEu_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_gradKEv_allocated", AD_a%gradKEv_a%associated())
  elseif (AD_a%gradKEv_a%associated()) then
    call rec%add(prefix//"_gradKEv"//suffix, AD_a%gradKEv_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_rv_x_v_allocated", AD_a%rv_x_v_a%associated())
  elseif (AD_a%rv_x_v_a%associated()) then
    call rec%add(prefix//"_rv_x_v"//suffix, AD_a%rv_x_v_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_rv_x_u_allocated", AD_a%rv_x_u_a%associated())
  elseif (AD_a%rv_x_u_a%associated()) then
    call rec%add(prefix//"_rv_x_u"//suffix, AD_a%rv_x_u_a)
  endif
end subroutine record_CorAdv_AD

!> Calculates the Coriolis and momentum advection contributions to the acceleration.
subroutine CorAdCalc(u, v, h, uh, vh, CAu, CAv, OBC, AD, G, GV, US, CS, pbv, Waves)
  type(ocean_grid_type),                      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v  !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: uh !< Zonal transport u*h*dy
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: vh !< Meridional transport v*h*dx
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out)   :: CAu !< Zonal acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out)   :: CAv !< Meridional acceleration due to Coriolis
                                                                  !! and momentum advection [L T-2 ~> m s-2].
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary control structure
  type(accel_diag_ptrs),                      intent(inout) :: AD  !< Storage for acceleration diagnostics
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  type(CoriolisAdv_CS),                       intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv
  type(porous_barrier_type),                  intent(in)    :: pbv !< porous barrier fractional cell metrics
  type(Wave_parameters_CS),         optional, pointer       :: Waves !< An optional pointer to Stokes drift CS

  type(RealArray_t) :: u_a, v_a, h_a, uh_a, vh_a, CAu_a, CAv_a
  type(RealArray_t) :: RV_a, PV_a, CAuS_a, CAvS_a
  type(RealArray_t) :: por_face_areaU_a, por_face_areaV_a
  type(CorAdv_AD_type) :: AD_a ! Copies of the acceleration diagnostics in AD
  logical :: Stokes_VF
  type(grid_core_type), pointer :: Gcore ! Persistent copies of the core ocean grid fields
  type(grid_OBC_type),  pointer :: Gobc  ! Persistent copies of the ocean grid OBC mask fields
  type(Box_t) :: bxH0  ! The h-point box of the whole column, [isc:iec, jsc:jec, 1:nz]
  type(Box_t) :: bxQ0  ! The B-grid-index box of the whole column, [IscB:IecB, JscB:JecB, 1:nz]

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_CoriolisAdv: Module must be initialized before it is used.")

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uh_a%alloc(lb=LBOUND(uh), ub=UBOUND(uh), source=uh)
  call vh_a%alloc(lb=LBOUND(vh), ub=UBOUND(vh), source=vh)
  call CAu_a%alloc(lb=LBOUND(CAu), ub=UBOUND(CAu), source=CAu)
  call CAv_a%alloc(lb=LBOUND(CAv), ub=UBOUND(CAv), source=CAv)
  call RV_a%alloc(lb=[u_a%lb(1),v_a%lb(2),u_a%lb(3)], ub=[u_a%ub(1),v_a%ub(2),u_a%ub(3)])
  call PV_a%alloc(lb=[u_a%lb(1),v_a%lb(2),u_a%lb(3)], ub=[u_a%ub(1),v_a%ub(2),u_a%ub(3)])
  call CAuS_a%alloc(lb=[u_a%lb(1),u_a%lb(2),1], ub=[u_a%ub(1),u_a%ub(2),G%ke])
  call CAvS_a%alloc(lb=[v_a%lb(1),v_a%lb(2),1], ub=[v_a%ub(1),v_a%ub(2),G%ke])
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)
  call AD_a%build_from(AD)

  Gcore => grid_core(G)
  Gobc => grid_OBC(G)

  Stokes_VF = .false.
  if (present(Waves)) then ; if (associated(Waves)) then
    Stokes_VF = Waves%Stokes_VF
  endif ; endif

  call bxH0%safe_alloc(ndims=3) ; call bxQ0%safe_alloc(ndims=3)
  call bxH0%set(idxS=[G%isc,G%jsc,1], idxE=[G%iec,G%jec,GV%ke])
  call bxQ0%set(idxS=[G%IscB,G%JscB,1], idxE=[G%IecB,G%JecB,GV%ke])

  call CorAdCalc_TR(bxH0, bxQ0, u_a, v_a, h_a, uh_a, vh_a, CAu_a, CAv_a, RV_a, PV_a, CAuS_a, &
                    CAvS_a, OBC, AD_a, Gcore, Gobc, GV%ke, GV%H_subroundoff, GV%Angstrom_H, &
                    US%m_to_L, US%m_s_to_L_T, CS%CorAdv_opts_CS, por_face_areaU_a, &
                    por_face_areaV_a, Waves, Stokes_VF)
  call bxH0%free() ; call bxQ0%free()
  call AD_a%copy_back(AD)

  call CorAdv_finalize_diagnostics(RV_a, PV_a, CAuS_a, CAvS_a, Stokes_VF, AD, G, GV%ke, CS)

  call CAu_a%copy2F(CAu)
  call CAv_a%copy2F(CAv)

  call u_a%free() ; call v_a%free() ; call h_a%free()
  call uh_a%free() ; call vh_a%free()
  call CAu_a%free() ; call CAv_a%free()
  call RV_a%free() ; call PV_a%free() ; call CAuS_a%free() ; call CAvS_a%free()
  call por_face_areaU_a%free() ; call por_face_areaV_a%free()

end subroutine CorAdCalc

!> The implementation of CorAdCalc (the root of its call tree): sets up the k-invariant fields,
!! then loops over iteration tiles, calling the setup, the selected scheme and the common terms
!! for each one.
subroutine CorAdCalc_TR_fortran(bxH0, bxQ0, u_a, v_a, h_a, uh_a, vh_a, CAu_a, CAv_a, RV_a, PV_a, &
                                CAuS_a, CAvS_a, OBC, AD_a, Gcore, Gobc, nz, H_subroundoff, &
                                Angstrom_H, m_to_L, m_s_to_L_T, CS, por_face_areaU_a, &
                                por_face_areaV_a, Waves, Stokes_VF)
  type(Box_t),                                intent(in)    :: bxH0 !< The h-point box of the
                                                                    !! whole column,
                                                                    !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),                                intent(in)    :: bxQ0 !< The B-grid-index box of the
                                                                    !! whole column,
                                                                    !! [IscB:IecB, JscB:JecB, 1:nz]
  type(grid_core_type),                       intent(in)    :: Gcore !< Persistent copies of the
                                                                     !! core ocean grid fields
  type(grid_OBC_type),                        intent(in)    :: Gobc !< Persistent copies of the
                                                                    !! ocean grid OBC mask fields
  integer, intent(in) :: nz !< The number of layers/levels in the vertical
  real, intent(in) :: H_subroundoff !< A thickness that is so small that it can be added to a
                                    !! thickness of Angstrom or larger without changing it at the
                                    !! bit level [H ~> m or kg m-2].
  real, intent(in) :: Angstrom_H !< A one-Angstrom thickness in the model thickness units
                                 !! [H ~> m or kg m-2].
  type(RealArray_t),         intent(in)    :: u_a   !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t),         intent(in)    :: v_a   !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t),         intent(in)    :: h_a   !< Layer thickness [H ~> m or kg m-2]
  type(RealArray_t),         intent(in)    :: uh_a  !< Zonal transport u*h*dy
                                                    !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t),         intent(in)    :: vh_a  !< Meridional transport v*h*dx
                                                    !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t),         intent(inout) :: CAu_a !< Zonal acceleration due to Coriolis
                                                    !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t),         intent(inout) :: CAv_a !< Meridional acceleration due to Coriolis
                                                    !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: RV_a     !< Diagnostic relative vorticity [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: PV_a     !< Diagnostic potential vorticity
                                               !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(inout) :: CAuS_a   !< Stokes contribution to CAu [L T-2 ~> m s-2]
  type(RealArray_t), intent(inout) :: CAvS_a   !< Stokes contribution to CAv [L T-2 ~> m s-2]
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary control structure
  type(CorAdv_AD_type),                       intent(inout) :: AD_a !< Copies of the acceleration
                                                                    !! diagnostics
  real, intent(in) :: m_to_L !< A constant that translates lengths in meters to the units of
                             !! horizontal lengths [L m-1 ~> 1]
  real, intent(in) :: m_s_to_L_T !< Convert lateral velocities from m s-1 to L T-1
                                 !! [L s T-1 m-1 ~> 1]
  type(CorAdv_opts_CS),                       intent(in)    :: CS  !< Options used in the Coriolis
                                                                   !! and momentum advection calculations
  type(RealArray_t), intent(in) :: por_face_areaU_a !< fractional open area of U-faces [nondim]
  type(RealArray_t), intent(in) :: por_face_areaV_a !< fractional open area of V-faces [nondim]
  type(Wave_parameters_CS),         optional, pointer       :: Waves !< An optional pointer to Stokes drift CS
  logical,                   intent(in)    :: Stokes_VF !< If true, include the Stokes drift

  ! Local variables
  type(RealArray_t) :: Area_h_a, Area_q_a
  real, dimension(:,:,:), contiguous, pointer :: u, v, h, uh, vh, CAu, CAv
  real, dimension(:,:), contiguous, pointer :: &
    Area_q      ! The sum of the ocean areas at the 4 adjacent thickness points [L2 ~> m2].
  real, dimension(:,:), contiguous, pointer :: &
    Area_h      ! The ocean area at h points [L2 ~> m2].  Area_h is used to find the
                ! average thickness in the denominator of q.  0 for land points.
  real, dimension(:,:,:), contiguous, pointer :: &
    PV, &       ! A diagnostic array of the potential vorticities [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    RV          ! A diagnostic array of the relative vorticities [T-1 ~> s-1].
  real, dimension(:,:,:), contiguous, pointer :: CAuS ! Stokes contribution to CAu [L T-2 ~> m s-2]
  real, dimension(:,:,:), contiguous, pointer :: CAvS ! Stokes contribution to CAv [L T-2 ~> m s-2]
  type(RealArray_t) :: KE_a, KEx_a, KEy_a
  type(RealArray_t) :: q_a, Ih_q_a, abs_vort_a, h_q_a, qS_a, q2_a, uh_center_a, vh_center_a
  real, dimension(:,:,:), contiguous, pointer :: &
    q, &        ! Layer potential vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    qS, &       ! Layer Stokes vorticity [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
    Ih_q, &     ! The inverse of thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1].
    h_q, &      ! The thickness interpolated to q points [H-1 ~> m-1 or m2 kg-1].
    abs_vort, & ! Absolute vorticity at q-points [T-1 ~> s-1].
    q2          ! Relative vorticity over thickness [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
  real, dimension(:,:,:), contiguous, pointer :: &
    KEx, &      ! The zonal gradient of Kinetic energy per unit mass [L T-2 ~> m s-2],
                ! KEx = d/dx KE.
    uh_center   ! Transport based on arithmetic mean h at u-points [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(:,:,:), contiguous, pointer :: &
    KE          ! Kinetic energy per unit mass [L2 T-2 ~> m2 s-2], KE = (u^2 + v^2)/2.
  real, dimension(:,:,:), contiguous, pointer :: &
    KEy, &      ! The meridional gradient of Kinetic energy per unit mass [L T-2 ~> m s-2],
                ! KEy = d/dy KE.
    vh_center   ! Transport based on arithmetic mean h at v-points [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(:,:), contiguous, pointer :: mask2dT, areaT
  real :: vol_neglect            ! A volume so small that is expected to be
                                 ! lost in roundoff [H L2 ~> m3 or kg].
  real :: area_neglect           ! An area so small that is expected to be
                                 ! lost in roundoff [L2 ~> m2].
  real :: eps_vel                ! A tiny, positive velocity [L T-1 ~> m s-1].
  real :: Fe_m2         ! Temporary variable associated with the ARAKAWA_LAMB_BLEND scheme [nondim]
  real :: rat_lin       ! Temporary variable associated with the ARAKAWA_LAMB_BLEND scheme [nondim]
  real :: h_tiny        ! A very small thickness [H ~> m or kg m-2].
  type(Box_t) :: bxH    ! The h-point iteration box of one tile, [isc:iec, jsc:jec, k_start:k_end]
  type(Box_t) :: bxQ    ! The B-grid-index iteration box of one tile, [IscB:IecB, JscB:JecB, k_start:k_end]
  type(Box_t) :: bxU    ! The u-point iteration box of one tile, [IscB:IecB, jsc:jec, k_start:k_end]
  type(Box_t) :: bxV    ! The v-point iteration box of one tile, [isc:iec, JscB:JecB, k_start:k_end]
  type(Box_t) :: bxQs   ! The scheme-dependent vorticity-point box of one tile,
                        ! [Is_q:Ie_q, Js_q:Je_q, k_start:k_end]
  type(Box_t) :: bxQs_pij ! bxQs grown by one point at the end in i and j, before the tile loop,
                         ! [Is_q:Ie_q+1, Js_q:Je_q+1, 1:nz]
  integer :: i, j, n, is, ie, js, je, nkblock, k_start, k_end
  integer :: Is_q, Ie_q, Js_q, Je_q  ! The scheme-dependent range of values at which vorticity is set.
  logical :: use_weno   ! True if using one of the WENO schemes
  integer :: stencil    ! Stencil size of WENO scheme

! To work, the following fields must be set outside of the usual
! is to ie range before this subroutine is called:
!   v(is-1:ie+2,js-1:je+1), u(is-1:ie+1,js-1:je+2), h(is-1:ie+2,js-1:je+2),
!   uh(is-1,ie,js:je+1) and vh(is:ie+1,js-1:je).

  is = bxH0%idxS(1) ; ie = bxH0%idxE(1) ; js = bxH0%idxS(2) ; je = bxH0%idxE(2)
  nkblock = merge(nz, CS%nkblock, CS%nkblock==0)
  vol_neglect = H_subroundoff * (1e-4 * m_to_L)**2
  area_neglect = (1e-4 * m_to_L)**2
  eps_vel = 1.0e-10*m_s_to_L_T
  h_tiny = Angstrom_H  ! Perhaps this should be set to h_neglect instead.

  stencil = Coriolis_scheme_stencil(CS%Coriolis_scheme_CS)

  use_weno = CS%Coriolis_scheme_CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO &
      .or. CS%Coriolis_scheme_CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO &
      .or. CS%Coriolis_scheme_CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO

  if (use_weno) then
    Is_q = is - stencil ; Ie_q = ie + stencil - 1 ; Js_q = js - stencil ; Je_q = je + stencil - 1
  else
    Is_q = bxQ0%idxS(1) - 1 ; Ie_q = bxQ0%idxE(1) + 1
    Js_q = bxQ0%idxS(2) - 1 ; Je_q = bxQ0%idxE(2) + 1
  endif

  call Area_q_a%alloc(lb=[u_a%lb(1),v_a%lb(2)], ub=[u_a%ub(1),v_a%ub(2)])
  call Area_h_a%alloc(lb=[h_a%lb(1),h_a%lb(2)], ub=[h_a%ub(1),h_a%ub(2)])
  call Area_q_a%view(Area_q) ; call Area_h_a%view(Area_h) ; call RV_a%view(RV)
  call PV_a%view(PV) ; call CAuS_a%view(CAuS) ; call CAvS_a%view(CAvS)
  call Gcore%mask2dT%view(mask2dT) ; call Gcore%areaT%view(areaT)

  call bxQs%safe_alloc(ndims=3)
  call bxQs%set(idxS=[Is_q,Js_q,1], idxE=[Ie_q,Je_q,nz])
  bxQs_pij = bxQs%growBy(lo=[0,0,0], hi=[1,1,0])

  !$omp target enter data map(alloc: Area_h, Area_q)

  do concurrent (j=bxQs_pij%idxS(2):bxQs_pij%idxE(2), &
                 I=bxQs_pij%idxS(1):bxQs_pij%idxE(1))
    Area_h(i,j) = mask2dT(i,j) * areaT(i,j)
  enddo

  if (associated(OBC)) then
    !$omp target update from(Area_h)

    do n=1,OBC%number_of_segments
      if (.not. OBC%segment(n)%on_pe) cycle
      I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
      if (OBC%segment(n)%is_N_or_S .and. (J >= Js_q) .and. (J <= Je_q)) then
        do i = max(Is_q,OBC%segment(n)%HI%isd), min(Ie_q+1,OBC%segment(n)%HI%ied)
          if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
            Area_h(i,j+1) = Area_h(i,j)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
            Area_h(i,j) = Area_h(i,j+1)
          endif
        enddo
      elseif (OBC%segment(n)%is_E_or_W .and. (I >= Is_q) .and. (I <= Ie_q)) then
        do j = max(Js_q,OBC%segment(n)%HI%jsd), min(Je_q+1,OBC%segment(n)%HI%jed)
          if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
            Area_h(i+1,j) = Area_h(i,j)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
            Area_h(i,j) = Area_h(i+1,j)
          endif
        enddo
      endif
    enddo

    !$omp target update to(Area_h)
  endif

  do concurrent (J=bxQs%idxS(2):bxQs%idxE(2), &
                 I=bxQs%idxS(1):bxQs%idxE(1))
    Area_q(i,j) = (Area_h(i,j) + Area_h(i+1,j+1)) + &
                  (Area_h(i+1,j) + Area_h(i,j+1))
  enddo

  !$omp target enter data map(alloc: CAuS, CAvS) if (Stokes_VF)

  ! Diagnostics
  !$omp target enter data map(alloc: RV) if (CS%do_RV_diag)
  !$omp target enter data map(alloc: PV) if (CS%do_PV_diag)

  ! TODO: Do this outside of the function

  ! Hoist AL_BLEND k-independent scalars out of the tile loop.
  Fe_m2 = 0.0 ; rat_lin = 0.0
  if (CS%Coriolis_scheme_CS%Coriolis_Scheme == AL_BLEND) then
    Fe_m2 = CS%F_eff_max_blend - 2.0
    rat_lin = 1.5 * Fe_m2 / max(CS%wt_lin_blend, 1.0e-16)
    if (CS%F_eff_max_blend <= 2.0) then ; Fe_m2 = -1. ; rat_lin = -1.0 ; endif
  endif

  ! Loop over iteration tiles.  Each tile covers the whole horizontal compute domain and
  ! nkblock layers; the tile boxes carry all of the index ranges used within the loop.
  call u_a%view(u) ; call v_a%view(v) ; call h_a%view(h)
  call uh_a%view(uh) ; call vh_a%view(vh)
  call CAu_a%view(CAu) ; call CAv_a%view(CAv)

  call bxH%safe_alloc(ndims=3) ; call bxQ%safe_alloc(ndims=3)
  call bxU%safe_alloc(ndims=3) ; call bxV%safe_alloc(ndims=3)
  call q_a%alloc(lb=[Area_q_a%lb(1),Area_q_a%lb(2),1], &
                 ub=[Area_q_a%ub(1),Area_q_a%ub(2),min(nkblock,nz)])
  call qS_a%alloc(lb=[Area_q_a%lb(1),Area_q_a%lb(2),1], &
                  ub=[Area_q_a%ub(1),Area_q_a%ub(2),min(nkblock,nz)])
  call Ih_q_a%alloc(lb=[Area_q_a%lb(1),Area_q_a%lb(2),1], &
                    ub=[Area_q_a%ub(1),Area_q_a%ub(2),min(nkblock,nz)])
  call h_q_a%alloc(lb=[Area_q_a%lb(1),Area_q_a%lb(2),1], &
                   ub=[Area_q_a%ub(1),Area_q_a%ub(2),min(nkblock,nz)])
  call abs_vort_a%alloc(lb=[Area_q_a%lb(1),Area_q_a%lb(2),1], &
                        ub=[Area_q_a%ub(1),Area_q_a%ub(2),min(nkblock,nz)])
  call q2_a%alloc(lb=[Area_q_a%lb(1),Area_q_a%lb(2),1], &
                  ub=[Area_q_a%ub(1),Area_q_a%ub(2),min(nkblock,nz)])
  call KEx_a%alloc(lb=[u_a%lb(1),u_a%lb(2),1], ub=[u_a%ub(1),u_a%ub(2),min(nkblock,nz)])
  call uh_center_a%alloc(lb=[u_a%lb(1),u_a%lb(2),1], &
                         ub=[u_a%ub(1),u_a%ub(2),min(nkblock,nz)])
  call KE_a%alloc(lb=[Area_h_a%lb(1),Area_h_a%lb(2),1], &
                  ub=[Area_h_a%ub(1),Area_h_a%ub(2),min(nkblock,nz)])
  call KEy_a%alloc(lb=[v_a%lb(1),v_a%lb(2),1], ub=[v_a%ub(1),v_a%ub(2),min(nkblock,nz)])
  call vh_center_a%alloc(lb=[v_a%lb(1),v_a%lb(2),1], &
                         ub=[v_a%ub(1),v_a%ub(2),min(nkblock,nz)])
  call q_a%view(q) ; call qS_a%view(qS) ; call Ih_q_a%view(Ih_q)
  call h_q_a%view(h_q) ; call abs_vort_a%view(abs_vort) ; call q2_a%view(q2)
  call KEx_a%view(KEx) ; call uh_center_a%view(uh_center) ; call KE_a%view(KE)
  call KEy_a%view(KEy) ; call vh_center_a%view(vh_center)

  !$omp target enter data map(alloc: abs_vort, q, Ih_q)
  !$omp target enter data map(alloc: h_q) if (use_weno)
  !$omp target enter data map(alloc: KE, KEx, KEy)
  ! TODO: These Stokes_VF fields seem associated with diagnostics
  !$omp target enter data map(alloc: qS) if (Stokes_VF)
  !$omp target enter data map(alloc: uh_center, vh_center) &
  !$omp   if (CS%Coriolis_scheme_CS%Coriolis_En_Dis)
  !$omp target enter data map(alloc: q2) &
  !$omp   if(AD_a%rv_x_u_a%associated() .or. AD_a%rv_x_v_a%associated())

  do k_start=1,nz,nkblock
    k_end = min(k_start+nkblock-1, nz)
    call bxH%set(idxS=[bxH0%idxS(1),bxH0%idxS(2),k_start], &
                 idxE=[bxH0%idxE(1),bxH0%idxE(2),k_end])
    call bxQ%set(idxS=[bxQ0%idxS(1),bxQ0%idxS(2),k_start], &
                 idxE=[bxQ0%idxE(1),bxQ0%idxE(2),k_end])
    call bxQs%set(idxS=[Is_q,Js_q,k_start], idxE=[Ie_q,Je_q,k_end])
    call bxU%set(idxS=[bxQ0%idxS(1),bxH0%idxS(2),k_start], &
                 idxE=[bxQ0%idxE(1),bxH0%idxE(2),k_end])
    call bxV%set(idxS=[bxH0%idxS(1),bxQ0%idxS(2),k_start], &
                 idxE=[bxH0%idxE(1),bxQ0%idxE(2),k_end])

    call q_a%rebound(lb=[Area_q_a%lb(1),Area_q_a%lb(2),k_start], &
                     ub=[Area_q_a%ub(1),Area_q_a%ub(2),k_end])
    call qS_a%rebound(lb=[Area_q_a%lb(1),Area_q_a%lb(2),k_start], &
                      ub=[Area_q_a%ub(1),Area_q_a%ub(2),k_end])
    call Ih_q_a%rebound(lb=[Area_q_a%lb(1),Area_q_a%lb(2),k_start], &
                        ub=[Area_q_a%ub(1),Area_q_a%ub(2),k_end])
    call h_q_a%rebound(lb=[Area_q_a%lb(1),Area_q_a%lb(2),k_start], &
                       ub=[Area_q_a%ub(1),Area_q_a%ub(2),k_end])
    call abs_vort_a%rebound(lb=[Area_q_a%lb(1),Area_q_a%lb(2),k_start], &
                            ub=[Area_q_a%ub(1),Area_q_a%ub(2),k_end])
    call q2_a%rebound(lb=[Area_q_a%lb(1),Area_q_a%lb(2),k_start], &
                      ub=[Area_q_a%ub(1),Area_q_a%ub(2),k_end])
    call KEx_a%rebound(lb=[u_a%lb(1),u_a%lb(2),k_start], ub=[u_a%ub(1),u_a%ub(2),k_end])
    call uh_center_a%rebound(lb=[u_a%lb(1),u_a%lb(2),k_start], ub=[u_a%ub(1),u_a%ub(2),k_end])
    call KE_a%rebound(lb=[Area_h_a%lb(1),Area_h_a%lb(2),k_start], &
                      ub=[Area_h_a%ub(1),Area_h_a%ub(2),k_end])
    call KEy_a%rebound(lb=[v_a%lb(1),v_a%lb(2),k_start], ub=[v_a%ub(1),v_a%ub(2),k_end])
    call vh_center_a%rebound(lb=[v_a%lb(1),v_a%lb(2),k_start], ub=[v_a%ub(1),v_a%ub(2),k_end])
    call q_a%view(q) ; call qS_a%view(qS) ; call Ih_q_a%view(Ih_q)
    call h_q_a%view(h_q) ; call abs_vort_a%view(abs_vort) ; call q2_a%view(q2)
    call KEx_a%view(KEx) ; call uh_center_a%view(uh_center) ; call KE_a%view(KE)
    call KEy_a%view(KEy) ; call vh_center_a%view(vh_center)

    ! Potential vorticity and the related quantities at q points used by every scheme.
    call CorAdv_setup(bxH, bxQ, bxU, bxV, bxQs, u_a, v_a, h_a, OBC, AD_a, por_face_areaU_a, &
                      por_face_areaV_a, Waves, Gcore, Area_h_a, Area_q_a, Stokes_VF, use_weno, &
                      vol_neglect, area_neglect, q_a, Ih_q_a, abs_vort_a, h_q_a, qS_a, q2_a, &
                      uh_center_a, vh_center_a, RV_a, PV_a, CS)

    ! Calculate KE and the gradient of KE
    call gradKE(bxQ, bxU, bxV, u_a, v_a, KE_a, KEx_a, KEy_a, Gcore, Gobc, CS%KE_Scheme, &
                CS%KE_use_limiter)
    ! TODO: Can KE be removed from this function?

    ! Calculate the Coriolis and momentum advection accelerations, CAu and CAv, with the selected
    ! scheme.  On a Cartesian grid, CAu =  q * vh - d(KE)/dx and CAv = - q * uh - d(KE)/dy.
    select case (CS%Coriolis_scheme_CS%Coriolis_Scheme)
    case (SADOURNY75_ENERGY, SADOURNY75_ENSTRO)
      call CorAdv_sadourny(bxH, bxU, bxV, u_a, v_a, uh_a, vh_a, q_a, uh_center_a, vh_center_a, &
                           Gcore, CAu_a, CAv_a, CS%Coriolis_scheme_CS)
    case (ARAKAWA_HSU90, ARAKAWA_LAMB81, AL_BLEND)
      call CorAdv_arakawa(bxH, bxQ, bxU, bxV, uh_a, vh_a, q_a, Ih_q_a, Fe_m2, rat_lin, Gcore, &
                          CAu_a, CAv_a, CS%wt_lin_blend, CS%Coriolis_scheme_CS)
    case (ROBUST_ENSTRO)
      call CorAdv_robust_enstro(bxU, bxV, u_a, v_a, h_a, uh_a, vh_a, abs_vort_a, eps_vel, h_tiny, &
                                Gcore, CAu_a, CAv_a, CS%PV_Adv_Scheme)
    case (wenovi7th_PV_ENSTRO, wenovi5th_PV_ENSTRO, wenovi3rd_PV_ENSTRO)
      call CorAdv_weno(bxU, bxV, u_a, v_a, uh_a, vh_a, q_a, abs_vort_a, h_q_a, Gcore, &
                       H_subroundoff, CAu_a, CAv_a, CS%weno_velocity_smooth, &
                       CS%Coriolis_scheme_CS)
    end select

    ! The Stokes-drift diagnostic, bounding, the kinetic energy gradient, and diagnostics.
    call CorAdv_common_terms(bxU, bxV, u_a, v_a, uh_a, vh_a, abs_vort_a, qS_a, q2_a, KEx_a, KEy_a, &
                             Stokes_VF, Gcore, CAu_a, CAv_a, CAuS_a, CAvS_a, AD_a, CS)
  enddo ! end of tile loop.
  !$omp target exit data map(delete: abs_vort, q, Ih_q)
  !$omp target exit data map(delete: h_q) if (use_weno)
  !$omp target exit data map(delete: KE, KEx, KEy)
  !$omp target exit data map(delete: qS) if (Stokes_VF)
  !$omp target exit data map(delete: uh_center, vh_center) &
  !$omp   if (CS%Coriolis_scheme_CS%Coriolis_En_Dis)
  !$omp target exit data map(delete: q2) &
  !$omp     if(AD_a%rv_x_u_a%associated() .or. AD_a%rv_x_v_a%associated())

  call q_a%free() ; call qS_a%free() ; call Ih_q_a%free()
  call h_q_a%free() ; call abs_vort_a%free() ; call q2_a%free()
  call KEx_a%free() ; call uh_center_a%free() ; call KE_a%free()
  call KEy_a%free() ; call vh_center_a%free()
  call bxH%free() ; call bxQ%free() ; call bxQs%free() ; call bxQs_pij%free()
  call bxU%free() ; call bxV%free()

  !$omp target exit data map(delete: Area_h, Area_q)

  ! TODO: Move outside function

  ! Diagnostics
  !$omp target exit data map(from: RV) if (CS%do_RV_diag)
  !$omp target exit data map(from: PV) if (CS%do_PV_diag)
  !$omp target exit data map(from: CAuS, CAvS) if (Stokes_VF)

  call Area_q_a%free() ; call Area_h_a%free()
end subroutine CorAdCalc_TR_fortran

!> Runs CorAdCalc_TR_fortran, a capture of its inputs and outputs, or the C++/AMReX bridge,
!! as selected by the CORADCALC_TR_MODE environment variable.
subroutine CorAdCalc_TR(bxH0, bxQ0, u_a, v_a, h_a, uh_a, vh_a, CAu_a, CAv_a, RV_a, PV_a, CAuS_a, &
                        CAvS_a, OBC, AD_a, Gcore, Gobc, nz, H_subroundoff, Angstrom_H, m_to_L, &
                        m_s_to_L_T, CS, por_face_areaU_a, por_face_areaV_a, Waves, Stokes_VF)
  type(Box_t),                                intent(in)    :: bxH0 !< The h-point box of the
                                                                    !! whole column,
                                                                    !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),                                intent(in)    :: bxQ0 !< The B-grid-index box of the
                                                                    !! whole column,
                                                                    !! [IscB:IecB, JscB:JecB, 1:nz]
  type(grid_core_type),                       intent(in)    :: Gcore !< Persistent copies of the
                                                                     !! core ocean grid fields
  type(grid_OBC_type),                        intent(in)    :: Gobc !< Persistent copies of the
                                                                    !! ocean grid OBC mask fields
  integer, intent(in) :: nz !< The number of layers/levels in the vertical
  real, intent(in) :: H_subroundoff !< A thickness that is so small that it can be added to a
                                    !! thickness of Angstrom or larger without changing it at the
                                    !! bit level [H ~> m or kg m-2].
  real, intent(in) :: Angstrom_H !< A one-Angstrom thickness in the model thickness units
                                 !! [H ~> m or kg m-2].
  type(RealArray_t),         intent(in)    :: u_a   !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t),         intent(in)    :: v_a   !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t),         intent(in)    :: h_a   !< Layer thickness [H ~> m or kg m-2]
  type(RealArray_t),         intent(in)    :: uh_a  !< Zonal transport u*h*dy
                                                    !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t),         intent(in)    :: vh_a  !< Meridional transport v*h*dx
                                                    !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t),         intent(inout) :: CAu_a !< Zonal acceleration due to Coriolis
                                                    !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t),         intent(inout) :: CAv_a !< Meridional acceleration due to Coriolis
                                                    !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: RV_a     !< Diagnostic relative vorticity [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: PV_a     !< Diagnostic potential vorticity
                                               !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(inout) :: CAuS_a   !< Stokes contribution to CAu [L T-2 ~> m s-2]
  type(RealArray_t), intent(inout) :: CAvS_a   !< Stokes contribution to CAv [L T-2 ~> m s-2]
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary control structure
  type(CorAdv_AD_type),                       intent(inout) :: AD_a !< Copies of the acceleration
                                                                    !! diagnostics
  real, intent(in) :: m_to_L !< A constant that translates lengths in meters to the units of
                             !! horizontal lengths [L m-1 ~> 1]
  real, intent(in) :: m_s_to_L_T !< Convert lateral velocities from m s-1 to L T-1
                                 !! [L s T-1 m-1 ~> 1]
  type(CorAdv_opts_CS),                       intent(in)    :: CS  !< Options used in the Coriolis
                                                                   !! and momentum advection calculations
  type(RealArray_t), intent(in) :: por_face_areaU_a !< fractional open area of U-faces [nondim]
  type(RealArray_t), intent(in) :: por_face_areaV_a !< fractional open area of V-faces [nondim]
  type(Wave_parameters_CS),         optional, pointer       :: Waves !< An optional pointer to Stokes drift CS
  logical,                   intent(in)    :: Stokes_VF !< If true, include the Stokes drift
  ! Local variables
  integer            :: mode, rc
  type(Box_C)        :: bxH0_c, bxQ0_c
  type(RealArray_C)  :: u_c, v_c, h_c, uh_c, vh_c, CAu_c, CAv_c, RV_c, PV_c, CAuS_c, CAvS_c, &
    por_face_areaU_c, por_face_areaV_c
  type(c_ptr) :: OBC_c
  type(CorAdv_AD_C) :: AD_a_c
  integer(c_int) :: nz_c
  type(CorAdv_opts_CS_C) :: CS_c
  type(c_ptr) :: Waves_c
  logical(c_bool) :: Stokes_VF_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "coradcalc_tr"
  mode = getenv_mode("CORADCALC_TR_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxH0", bxH0)
        call rec%add("_bxQ0", bxQ0)
        call rec%add("_u", u_a)
        call rec%add("_v", v_a)
        call rec%add("_h", h_a)
        call rec%add("_uh", uh_a)
        call rec%add("_vh", vh_a)
        call rec%add("_CAu_before", CAu_a)
        call rec%add("_CAv_before", CAv_a)
        call rec%add("_RV_before", RV_a)
        call rec%add("_PV_before", PV_a)
        call rec%add("_CAuS_before", CAuS_a)
        call rec%add("_CAvS_before", CAvS_a)
        call record_CorAdv_AD(rec, "_AD", AD_a, "")
        call record_CorAdv_AD(rec, "_AD", AD_a, "_before")
        call rec%add("_mask2dT", Gcore%mask2dT)
        call rec%add("_areaT", Gcore%areaT)
        call rec%add("_areaCu", Gcore%areaCu)
        call rec%add("_areaCv", Gcore%areaCv)
        call rec%add("_IareaT", Gcore%IareaT)
        call rec%add("_mask2dCu", Gcore%mask2dCu)
        call rec%add("_mask2dCv", Gcore%mask2dCv)
        call rec%add("_dyCv", Gcore%dyCv)
        call rec%add("_dxCu", Gcore%dxCu)
        call rec%add("_dy_Cu", Gcore%dy_Cu)
        call rec%add("_dx_Cv", Gcore%dx_Cv)
        call rec%add("_dyBu", Gcore%dyBu)
        call rec%add("_dxBu", Gcore%dxBu)
        call rec%add("_mask2dBu", Gcore%mask2dBu)
        call rec%add("_IareaBu", Gcore%IareaBu)
        call rec%add("_CoriolisBu", Gcore%CoriolisBu)
        call rec%add("_IdxCu", Gcore%IdxCu)
        call rec%add("_IdyCv", Gcore%IdyCv)
        call rec%add("_IdxCv", Gcore%IdxCv)
        call rec%add("_IdyCu", Gcore%IdyCu)
        call rec%add("_IdxCu_OBCmask", Gobc%IdxCu_OBCmask)
        call rec%add("_IdyCv_OBCmask", Gobc%IdyCv_OBCmask)
        call rec%add("_nz", nz)
        call rec%add("_H_subroundoff", H_subroundoff)
        call rec%add("_Angstrom_H", Angstrom_H)
        call rec%add("_m_to_L", m_to_L)
        call rec%add("_m_s_to_L_T", m_s_to_L_T)
        call record_CorAdv_opts_CS(rec, "_CS", CS)
        call rec%add("_por_face_areaU", por_face_areaU_a)
        call rec%add("_por_face_areaV", por_face_areaV_a)
        call rec%add("_Stokes_VF", Stokes_VF)
      endif

      call CorAdCalc_TR_fortran(bxH0, bxQ0, u_a, v_a, h_a, uh_a, vh_a, CAu_a, CAv_a, RV_a, PV_a, &
                                CAuS_a, CAvS_a, OBC, AD_a, Gcore, Gobc, nz, H_subroundoff, &
                                Angstrom_H, m_to_L, m_s_to_L_T, CS, por_face_areaU_a, &
                                por_face_areaV_a, Waves, Stokes_VF)

      if (capture) then
        call rec%add("_CAu_after", CAu_a)
        call rec%add("_CAv_after", CAv_a)
        call rec%add("_RV_after", RV_a)
        call rec%add("_PV_after", PV_a)
        call rec%add("_CAuS_after", CAuS_a)
        call rec%add("_CAvS_after", CAvS_a)
        call record_CorAdv_AD(rec, "_AD", AD_a, "_after")
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxH0_c = bxH0%to_c()
      bxQ0_c = bxQ0%to_c()
      u_c = u_a%to_c()
      v_c = v_a%to_c()
      h_c = h_a%to_c()
      uh_c = uh_a%to_c()
      vh_c = vh_a%to_c()
      CAu_c = CAu_a%to_c()
      CAv_c = CAv_a%to_c()
      RV_c = RV_a%to_c()
      PV_c = PV_a%to_c()
      CAuS_c = CAuS_a%to_c()
      CAvS_c = CAvS_a%to_c()
      if (associated(OBC)) then
        OBC_c = c_loc(OBC)
      else
        OBC_c = c_null_ptr
      endif
      AD_a_c = CorAdv_AD_to_c(AD_a)
      nz_c = nz
      CS_c = CorAdv_opts_CS_to_c(CS)
      por_face_areaU_c = por_face_areaU_a%to_c()
      por_face_areaV_c = por_face_areaV_a%to_c()
      if (Stokes_VF) then
        Waves_c = c_loc(Waves)
      else
        Waves_c = c_null_ptr
      endif
      Stokes_VF_c = Stokes_VF
      call turbotmp_coradcalc_tr_bridge(bxH0_c, bxQ0_c, u_c, v_c, h_c, uh_c, vh_c, CAu_c, CAv_c, &
                                        RV_c, PV_c, CAuS_c, CAvS_c, OBC_c, AD_a_c, Gcore%C, &
                                        Gobc%C, nz_c, H_subroundoff, Angstrom_H, m_to_L, &
                                        m_s_to_L_T, CS_c, por_face_areaU_c, por_face_areaV_c, &
                                        Waves_c, Stokes_VF_c)
#endif

    case default
      call CorAdCalc_TR_fortran(bxH0, bxQ0, u_a, v_a, h_a, uh_a, vh_a, CAu_a, CAv_a, RV_a, PV_a, &
                                CAuS_a, CAvS_a, OBC, AD_a, Gcore, Gobc, nz, H_subroundoff, &
                                Angstrom_H, m_to_L, m_s_to_L_T, CS, por_face_areaU_a, &
                                por_face_areaV_a, Waves, Stokes_VF)

  end select

end subroutine CorAdCalc_TR


!> Calculates the potential vorticity and the related quantities at q points that every
!! Coriolis scheme uses, over one iteration tile.
subroutine CorAdv_setup_fortran(bxH, bxQ, bxU, bxV, bxQs, u_a, v_a, h_a, OBC, AD_a, &
                                por_face_areaU_a, por_face_areaV_a, Waves, Gcore, Area_h_a, &
                                Area_q_a, Stokes_VF, use_weno, vol_neglect, area_neglect, q_a, &
                                Ih_q_a, abs_vort_a, h_q_a, qS_a, q2_a, uh_center_a, vh_center_a, &
                                RV_a, PV_a, CS)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxQ !< The B-grid-index iteration box of this tile,
                                                !! [IscB:IecB, JscB:JecB, ksc:kec]
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(Box_t),             intent(in)    :: bxQs !< The scheme-dependent vorticity-point box,
                                                 !! [Is_q:Ie_q, Js_q:Je_q, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a         !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a         !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: h_a         !< Layer thickness [H ~> m or kg m-2]
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary control structure
  type(CorAdv_AD_type),    intent(in)    :: AD_a !< Copies of the acceleration diagnostics
  type(RealArray_t), intent(in)    :: por_face_areaU_a !< fractional open area of U-faces [nondim]
  type(RealArray_t), intent(in)    :: por_face_areaV_a !< fractional open area of V-faces [nondim]
  type(Wave_parameters_CS), optional, pointer :: Waves !< An optional pointer to Stokes drift CS
  type(RealArray_t), intent(in)    :: Area_h_a    !< The ocean area at h points [L2 ~> m2]
  type(RealArray_t), intent(in)    :: Area_q_a    !< The sum of the ocean areas at the 4
                                                  !! adjacent thickness points [L2 ~> m2]
  logical,                 intent(in)    :: Stokes_VF !< If true, include the Stokes drift
  logical,                 intent(in)    :: use_weno  !< True if using one of the WENO schemes
  real,                    intent(in)    :: vol_neglect !< A negligible volume [H L2 ~> m3 or kg]
  real,                    intent(in)    :: area_neglect !< A negligible area [L2 ~> m2]
  type(RealArray_t), intent(inout) :: q_a         !< Layer potential vorticity
                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(inout) :: Ih_q_a      !< Inverse of thickness interpolated to q points
                                                  !! [H-1 ~> m-1 or m2 kg-1]
  type(RealArray_t), intent(inout) :: abs_vort_a  !< Absolute vorticity at q-points [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: h_q_a       !< Thickness interpolated to q points
                                                  !! [H ~> m or kg m-2]
  type(RealArray_t), intent(inout) :: qS_a        !< Layer Stokes vorticity
                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(inout) :: q2_a        !< Relative vorticity over thickness
                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(inout) :: uh_center_a !< Transport based on arithmetic mean h at
                                                  !! u-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(inout) :: vh_center_a !< Transport based on arithmetic mean h at
                                                  !! v-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(inout) :: RV_a        !< Diagnostic relative vorticity [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: PV_a        !< Diagnostic potential vorticity
                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(CorAdv_opts_CS),    intent(in)    :: CS  !< Options used in the Coriolis and momentum
                                                !! advection calculations

  ! Local variables
  type(Box_t) :: bxH_mij  ! bxH grown by one point at the start in i and j,
                          ! [isc-1:iec, jsc-1:jec, ksc:kec]
  type(Box_t) :: bxQ_gij  ! bxQ grown by one point on both sides in i and j,
                          ! [IscB-1:IecB+1, JscB-1:JecB+1, ksc:kec]
  type(Box_t) :: bxQs_pi  ! bxQs grown by one point at the end in i,
                          ! [Is_q:Ie_q+1, Js_q:Je_q, ksc:kec]
  type(Box_t) :: bxQs_pj  ! bxQs grown by one point at the end in j,
                          ! [Is_q:Ie_q, Js_q:Je_q+1, ksc:kec]
  type(Box_t) :: bxUx  ! bxU grown by one point at the end in i and the start in j,
                       ! [IscB:IecB+1, jsc-1:jec, ksc:kec]
  type(Box_t) :: bxVx  ! bxV grown by one point at the start in i and the end in j,
                       ! [isc-1:iec, JscB:JecB+1, ksc:kec]
  real, dimension(:,:,:), contiguous, pointer :: u, v, h, q, Ih_q, abs_vort, h_q, qS, q2, &
    uh_center, vh_center, RV, PV
  real, dimension(:,:), contiguous, pointer :: Area_h, Area_q
  real, dimension(:,:), contiguous, pointer :: dyCv, dxCu, dy_Cu, dx_Cv, dyBu, dxBu, mask2dBu, &
    IareaBu, CoriolisBu
  real, dimension(:,:,:), contiguous, pointer :: por_face_areaU, por_face_areaV
  real, dimension(q_a%lb(1):q_a%ub(1),q_a%lb(2):q_a%ub(2),bxH%idxS(3):bxH%idxE(3)) :: &
    dvdx, dudy, &   ! Contributions to the circulation around q-points [L2 T-1 ~> m2 s-1]
    dvSdx, duSdy, & ! idem. for Stokes drift [L2 T-1 ~> m2 s-1]
    rel_vort, &     ! Relative vorticity at q-points [T-1 ~> s-1].
    stk_vort        ! Stokes vorticity at q-points [T-1 ~> s-1].
  real, dimension(uh_center_a%lb(1):uh_center_a%ub(1), &
                  uh_center_a%lb(2):uh_center_a%ub(2),bxH%idxS(3):bxH%idxE(3)) :: &
    hArea_u         ! The cell area weighted thickness interpolated to u points
                    ! times the effective areas [H L2 ~> m3 or kg].
  real, dimension(vh_center_a%lb(1):vh_center_a%ub(1), &
                  vh_center_a%lb(2):vh_center_a%ub(2),bxH%idxS(3):bxH%idxE(3)) :: &
    hArea_v         ! The cell area weighted thickness interpolated to v points
                    ! times the effective areas [H L2 ~> m3 or kg].
  real :: hArea_q   ! The sum of area times thickness of the cells
                    ! surrounding a q point [H L2 ~> m3 or kg].
  integer :: i, j, k, Isq, Ieq, Jsq, Jeq, ksc, kec, n
  integer :: Is_q, Ie_q, Js_q, Je_q  ! The scheme-dependent range of values at which vorticity is set.

  call u_a%view(u) ; call v_a%view(v) ; call h_a%view(h)
  call Area_h_a%view(Area_h) ; call Area_q_a%view(Area_q) ; call q_a%view(q)
  call Ih_q_a%view(Ih_q) ; call abs_vort_a%view(abs_vort) ; call h_q_a%view(h_q)
  call qS_a%view(qS) ; call q2_a%view(q2) ; call uh_center_a%view(uh_center)
  call vh_center_a%view(vh_center) ; call RV_a%view(RV) ; call PV_a%view(PV)
  call Gcore%dyCv%view(dyCv) ; call Gcore%dxCu%view(dxCu) ; call Gcore%dy_Cu%view(dy_Cu)
  call Gcore%dx_Cv%view(dx_Cv) ; call Gcore%dyBu%view(dyBu) ; call Gcore%dxBu%view(dxBu)
  call Gcore%mask2dBu%view(mask2dBu) ; call Gcore%IareaBu%view(IareaBu)
  call Gcore%CoriolisBu%view(CoriolisBu)
  call por_face_areaU_a%view(por_face_areaU) ; call por_face_areaV_a%view(por_face_areaV)

  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  ksc = bxH%idxS(3) ; kec = bxH%idxE(3)
  Is_q = bxQs%idxS(1) ; Ie_q = bxQs%idxE(1) ; Js_q = bxQs%idxS(2) ; Je_q = bxQs%idxE(2)
  bxUx = bxU%growBy(lo=[0,1,0], hi=[1,0,0])
  bxVx = bxV%growBy(lo=[1,0,0], hi=[0,1,0])
  bxH_mij = bxH%growBy(lo=[1,1,0], hi=[0,0,0])
  bxQ_gij = bxQ%growBy(lo=[1,1,0], hi=[1,1,0])
  bxQs_pi = bxQs%growHi(dim=1, n=1)
  bxQs_pj = bxQs%growHi(dim=2, n=1)

  !$omp target enter data map(alloc: dvdx, dudy)
  !$omp target enter data map(alloc: hArea_u, hArea_v)
  !$omp target enter data map(alloc: rel_vort)
  !$omp target enter data map(alloc: dvSdx, duSdy, stk_vort) if (Stokes_VF)

  ! Here the second order accurate layer potential vorticities, q,
  ! are calculated.  hq is  second order accurate in space.  Relative
  ! vorticity is second order accurate everywhere with free slip b.c.s,
  ! but only first order accurate at boundaries with no slip b.c.s.
  ! First calculate the contributions to the circulation around the q-point.
  if (Stokes_VF) then
    if (CS%do_CAS_diag) then
      do concurrent (k=bxQs%idxS(3):bxQs%idxE(3), J=bxQs%idxS(2):bxQs%idxE(2), &
                     I=bxQs%idxS(1):bxQs%idxE(1))
        dvSdx(I,J,k) = (-Waves%us_y(i+1,J,k)*dyCv(i+1,J)) - &
                         (-Waves%us_y(i,J,k)*dyCv(i,J))
        duSdy(I,J,k) = (-Waves%us_x(I,j+1,k)*dxCu(I,j+1)) - &
                         (-Waves%us_x(I,j,k)*dxCu(I,j))
      enddo
    endif
    if (.not. Waves%Passive_Stokes_VF) then
      do concurrent (k=bxQs%idxS(3):bxQs%idxE(3), J=bxQs%idxS(2):bxQs%idxE(2), &
                     I=bxQs%idxS(1):bxQs%idxE(1))
        dvdx(I,J,k) = ((v(i+1,J,k)-Waves%us_y(i+1,J,k))*dyCv(i+1,J)) - &
                        ((v(i,J,k)-Waves%us_y(i,J,k))*dyCv(i,J))
        dudy(I,J,k) = ((u(I,j+1,k)-Waves%us_x(I,j+1,k))*dxCu(I,j+1)) - &
                        ((u(I,j,k)-Waves%us_x(I,j,k))*dxCu(I,j))
      enddo
    else
      do concurrent (k=bxQs%idxS(3):bxQs%idxE(3), J=bxQs%idxS(2):bxQs%idxE(2), &
                     I=bxQs%idxS(1):bxQs%idxE(1))
        dvdx(I,J,k) = (v(i+1,J,k)*dyCv(i+1,J)) - (v(i,J,k)*dyCv(i,J))
        dudy(I,J,k) = (u(I,j+1,k)*dxCu(I,j+1)) - (u(I,j,k)*dxCu(I,j))
      enddo
    endif
  else
    do concurrent (k=bxQs%idxS(3):bxQs%idxE(3), J=bxQs%idxS(2):bxQs%idxE(2), &
                   I=bxQs%idxS(1):bxQs%idxE(1))
      dvdx(I,J,k) = (v(i+1,J,k)*dyCv(i+1,J)) - (v(i,J,k)*dyCv(i,J))
      dudy(I,J,k) = (u(I,j+1,k)*dxCu(I,j+1)) - (u(I,j,k)*dxCu(I,j))
    enddo
  endif

  do concurrent (k=bxQs_pi%idxS(3):bxQs_pi%idxE(3), J=bxQs_pi%idxS(2):bxQs_pi%idxE(2), &
                 i=bxQs_pi%idxS(1):bxQs_pi%idxE(1))
    hArea_v(i,J,k) = 0.5*((Area_h(i,j) * h(i,j,k)) + (Area_h(i,j+1) * h(i,j+1,k)))
  enddo

  do concurrent (k=bxQs_pj%idxS(3):bxQs_pj%idxE(3), j=bxQs_pj%idxS(2):bxQs_pj%idxE(2), &
                 I=bxQs_pj%idxS(1):bxQs_pj%idxE(1))
    hArea_u(I,j,k) = 0.5*((Area_h(i,j) * h(i,j,k)) + (Area_h(i+1,j) * h(i+1,j,k)))
  enddo

  if (CS%Coriolis_scheme_CS%Coriolis_En_Dis) then
    do concurrent (k=bxVx%idxS(3):bxVx%idxE(3), J=bxVx%idxS(2):bxVx%idxE(2), &
                   I=bxVx%idxS(1):bxVx%idxE(1))
      uh_center(I,j,k) = 0.5 * ((dy_Cu(I,j)*por_face_areaU(I,j,k)) * u(I,j,k)) * (h(i,j,k) + h(i+1,j,k))
    enddo

    do concurrent (k=bxUx%idxS(3):bxUx%idxE(3), J=bxUx%idxS(2):bxUx%idxE(2), &
                   i=bxUx%idxS(1):bxUx%idxE(1))
      vh_center(i,J,k) = 0.5 * ((dx_Cv(i,J)*por_face_areaV(i,J,k)) * v(i,J,k)) * (h(i,j,k) + h(i,j+1,k))
    enddo
  endif

  ! Adjust circulation components to relative vorticity and thickness projected onto
  ! velocity points on open boundaries.
  if (associated(OBC)) then
    !$omp target update from(Area_h)
    do k=ksc,kec ! TODO: port (OBC GPU path not yet implemented)
      !$omp target update from(dvdx(:,:,k), dudy(:,:,k))
      !$omp target update from(hArea_u(:,:,k), hArea_v(:,:,k))
      !$omp target update from(uh_center(:,:,k), vh_center(:,:,k))

      do n=1,OBC%number_of_segments
        if (.not. OBC%segment(n)%on_pe) cycle

        I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB

        if (OBC%segment(n)%is_N_or_S .and. (J >= Js_q) .and. (J <= Je_q)) then
          select case (OBC%vorticity_config)
          case (OBC_VORTICITY_ZERO)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              dvdx(I,J,k) = 0. ; dudy(I,J,k) = 0.
            enddo
          case (OBC_VORTICITY_FREESLIP)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              dudy(I,J,k) = 0.
            enddo
          case (OBC_VORTICITY_COMPUTED)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                dudy(I,J,k) = 2.0*(OBC%segment(n)%tangential_vel(I,J,k) - u(I,j,k))*dxCu(I,j)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
                dudy(I,J,k) = 2.0*(u(I,j+1,k) - OBC%segment(n)%tangential_vel(I,J,k))*dxCu(I,j+1)
              endif
            enddo
          case (OBC_VORTICITY_SPECIFIED)
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                dudy(I,J,k) = OBC%segment(n)%tangential_grad(I,J,k)*dxCu(I,j)*dyBu(I,J)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
                dudy(I,J,k) = OBC%segment(n)%tangential_grad(I,J,k)*dxCu(I,j+1)*dyBu(I,J)
              endif
            enddo
          end select

          ! Project thicknesses across OBC points with a no-gradient condition.
          do i=max(Is_q,OBC%segment(n)%HI%isd), min(Ie_q+1,OBC%segment(n)%HI%ied)
            if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
              hArea_v(i,J,k) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j,k)
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
              hArea_v(i,J,k) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j+1,k)
            endif
          enddo

          if (CS%Coriolis_scheme_CS%Coriolis_En_Dis) then
            do i=max(Isq,OBC%segment(n)%HI%isd), min(Ieq+1,OBC%segment(n)%HI%ied)
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                vh_center(i,J,k) = (dx_Cv(i,J)*por_face_areaV(i,J,k)) * v(i,J,k) * h(i,j,k)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
                vh_center(i,J,k) = (dx_Cv(i,J)*por_face_areaV(i,J,k)) * v(i,J,k) * h(i,j+1,k)
              endif
            enddo
          endif
        elseif (OBC%segment(n)%is_E_or_W .and. (I >= Is_q) .and. (I <= Ie_q)) then
          select case (OBC%vorticity_config)
          case (OBC_VORTICITY_ZERO)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              dvdx(I,J,k) = 0. ; dudy(I,J,k) = 0.
            enddo
          case (OBC_VORTICITY_FREESLIP)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              dvdx(I,J,k) = 0.
            enddo
          case (OBC_VORTICITY_COMPUTED)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                dvdx(I,J,k) = 2.0*(OBC%segment(n)%tangential_vel(I,J,k) - v(i,J,k))*dyCv(i,J)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
                dvdx(I,J,k) = 2.0*(v(i+1,J,k) - OBC%segment(n)%tangential_vel(I,J,k))*dyCv(i+1,J)
              endif
            enddo
          case (OBC_VORTICITY_SPECIFIED)
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                dvdx(I,J,k) = OBC%segment(n)%tangential_grad(I,J,k)*dyCv(i,J)*dxBu(I,J)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
                dvdx(I,J,k) = OBC%segment(n)%tangential_grad(I,J,k)*dyCv(i+1,J)*dxBu(I,J)
              endif
            enddo
          end select

          ! Project thicknesses across OBC points with a no-gradient condition.
          do j=max(Js_q,OBC%segment(n)%HI%jsd), min(Je_q+1,OBC%segment(n)%HI%jed)
            if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
              hArea_u(I,j,k) = 0.5*(Area_h(i,j) + Area_h(i+1,j)) * h(i,j,k)
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
              hArea_u(I,j,k) = 0.5*(Area_h(i,j) + Area_h(i+1,j)) * h(i+1,j,k)
            endif
          enddo

          if (CS%Coriolis_scheme_CS%Coriolis_En_Dis) then
            do j=max(Jsq,OBC%segment(n)%HI%jsd), min(Jeq+1,OBC%segment(n)%HI%jed)
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                uh_center(I,j,k) = (dy_Cu(I,j)*por_face_areaU(I,j,k)) * u(I,j,k) * h(i,j,k)
              else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
                uh_center(I,j,k) = (dy_Cu(I,j)*por_face_areaU(I,j,k)) * u(I,j,k) * h(i+1,j,k)
              endif
            enddo
          endif
        endif
      enddo

      ! Now project thicknesses across cell-corner points in the OBCs.  The two
      ! projections have to occur in sequence and can not be combined easily.
      do n=1,OBC%number_of_segments
        if (.not. OBC%segment(n)%on_pe) cycle
        I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%is_N_or_S .and. (J >= Js_q) .and. (J <= Je_q)) then
          do I = max(Is_q,OBC%segment(n)%HI%IsdB), min(Ie_q,OBC%segment(n)%HI%IedB)
            if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
              if (Area_h(i,j) + Area_h(i+1,j) > 0.0) then
                hArea_u(I,j+1,k) = hArea_u(I,j,k) * ((Area_h(i,j+1) + Area_h(i+1,j+1)) / &
                                                         (Area_h(i,j) + Area_h(i+1,j)))
              else ; hArea_u(I,j+1,k) = 0.0 ; endif
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
              if (Area_h(i,j+1) + Area_h(i+1,j+1) > 0.0) then
                hArea_u(I,j,k) = hArea_u(I,j+1,k) * ((Area_h(i,j) + Area_h(i+1,j)) / &
                                                         (Area_h(i,j+1) + Area_h(i+1,j+1)))
              else ; hArea_u(I,j,k) = 0.0 ; endif
            endif
          enddo
        elseif (OBC%segment(n)%is_E_or_W .and. (I >= Is_q) .and. (I <= Ie_q)) then
          do J = max(Js_q,OBC%segment(n)%HI%JsdB), min(Je_q,OBC%segment(n)%HI%JedB)
            if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
              if (Area_h(i,j) + Area_h(i,j+1) > 0.0) then
                hArea_v(i+1,J,k) = hArea_v(i,J,k) * ((Area_h(i+1,j) + Area_h(i+1,j+1)) / &
                                                         (Area_h(i,j) + Area_h(i,j+1)))
              else ; hArea_v(i+1,J,k) = 0.0 ; endif
            else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
              hArea_v(i,J,k) = 0.5 * (Area_h(i,j) + Area_h(i,j+1)) * h(i,j+1,k)
              if (Area_h(i+1,j) + Area_h(i+1,j+1) > 0.0) then
                hArea_v(i,J,k) = hArea_v(i+1,J,k) * ((Area_h(i,j) + Area_h(i,j+1)) / &
                                                         (Area_h(i+1,j) + Area_h(i+1,j+1)))
              else ; hArea_v(i,J,k) = 0.0 ; endif
            endif
          enddo
        endif
      enddo

      !$omp target update to(dvdx(:,:,k), dudy(:,:,k))
      !$omp target update to(hArea_u(:,:,k), hArea_v(:,:,k))
      !$omp target update to(uh_center(:,:,k), vh_center(:,:,k))
    enddo ! k=ksc,kec OBC
  endif

  if (CS%no_slip) then
    do concurrent (k=bxQs%idxS(3):bxQs%idxE(3), J=bxQs%idxS(2):bxQs%idxE(2), &
                   I=bxQs%idxS(1):bxQs%idxE(1))
      rel_vort(I,J,k) = (2.0 - mask2dBu(I,J)) * (dvdx(I,J,k) - dudy(I,J,k)) * IareaBu(I,J)
    enddo

    if (Stokes_VF) then
      if (CS%do_CAS_diag) then
        do concurrent (k=bxQ_gij%idxS(3):bxQ_gij%idxE(3), J=bxQ_gij%idxS(2):bxQ_gij%idxE(2), &
                       I=bxQ_gij%idxS(1):bxQ_gij%idxE(1))
          stk_vort(I,J,k) = (2.0 - mask2dBu(I,J)) * (dvSdx(I,J,k) - duSdy(I,J,k)) * IareaBu(I,J)
        enddo
      endif
    endif
  else
    do concurrent (k=bxQs%idxS(3):bxQs%idxE(3), J=bxQs%idxS(2):bxQs%idxE(2), &
                   I=bxQs%idxS(1):bxQs%idxE(1))
      rel_vort(I,J,k) = mask2dBu(I,J) * (dvdx(I,J,k) - dudy(I,J,k)) * IareaBu(I,J)
    enddo

    if (Stokes_VF) then
      if (CS%do_CAS_diag) then
        do concurrent (k=bxQ_gij%idxS(3):bxQ_gij%idxE(3), J=bxQ_gij%idxS(2):bxQ_gij%idxE(2), &
                       I=bxQ_gij%idxS(1):bxQ_gij%idxE(1))
          stk_vort(I,J,k) = (2.0 - mask2dBu(I,J)) * (dvSdx(I,J,k) - duSdy(I,J,k)) * IareaBu(I,J)
        enddo
      endif
    endif
  endif

  do concurrent (k=bxQs%idxS(3):bxQs%idxE(3), J=bxQs%idxS(2):bxQs%idxE(2), &
                 I=bxQs%idxS(1):bxQs%idxE(1))
    abs_vort(I,J,k) = CoriolisBu(I,J) + rel_vort(I,J,k)
  enddo

  do concurrent (k=bxQs%idxS(3):bxQs%idxE(3), J=bxQs%idxS(2):bxQs%idxE(2), &
                 I=bxQs%idxS(1):bxQs%idxE(1)) DO_LOCALITY(local(hArea_q))
    hArea_q = (hArea_u(I,j,k) + hArea_u(I,j+1,k)) + (hArea_v(i,J,k) + hArea_v(i+1,J,k))
    Ih_q(I,J,k) = Area_q(I,J) / (hArea_q + vol_neglect)
    q(I,J,k) = abs_vort(I,J,k) * Ih_q(I,J,k)
  enddo

  ! NOTE: `h_q` is only used by WENO and was pulled out of the above loop to
  !   improve GPU performance, but it may need to be moved back.
  if (use_weno) then
    do concurrent (k=bxQs%idxS(3):bxQs%idxE(3), J=bxQs%idxS(2):bxQs%idxE(2), &
                   I=bxQs%idxS(1):bxQs%idxE(1)) DO_LOCALITY(local(hArea_q))
      hArea_q = (hArea_u(I,j,k) + hArea_u(I,j+1,k)) + (hArea_v(i,J,k) + hArea_v(i+1,J,k))
      h_q(I,J,k) = hArea_q / max(Area_q(I,J), area_neglect)
    enddo
  endif

  if (Stokes_VF) then
    if (CS%do_CAS_diag) then
      do concurrent (k=bxH_mij%idxS(3):bxH_mij%idxE(3), J=bxH_mij%idxS(2):bxH_mij%idxE(2), &
                     I=bxH_mij%idxS(1):bxH_mij%idxE(1))
        qS(I,J,k) = stk_vort(I,J,k) * Ih_q(I,J,k)
      enddo
    endif
  endif

  if (CS%do_RV_diag) then
    do concurrent (k=bxQ_gij%idxS(3):bxQ_gij%idxE(3), J=bxQ_gij%idxS(2):bxQ_gij%idxE(2), &
                   I=bxQ_gij%idxS(1):bxQ_gij%idxE(1))
      RV(I,J,k) = rel_vort(I,J,k)
    enddo
  endif

  if (CS%do_PV_diag) then
    do concurrent (k=bxQ_gij%idxS(3):bxQ_gij%idxE(3), J=bxQ_gij%idxS(2):bxQ_gij%idxE(2), &
                   I=bxQ_gij%idxS(1):bxQ_gij%idxE(1))
      PV(I,J,k) = q(I,J,k)
    enddo
  endif

  if (AD_a%rv_x_v_a%associated() .or. AD_a%rv_x_u_a%associated()) then
    do concurrent (k=bxQ_gij%idxS(3):bxQ_gij%idxE(3), J=bxQ_gij%idxS(2):bxQ_gij%idxE(2), &
                   I=bxQ_gij%idxS(1):bxQ_gij%idxE(1))
      q2(I,J,k) = rel_vort(I,J,k) * Ih_q(I,J,k)
    enddo
  endif

  !$omp target exit data map(delete: dvdx, dudy)
  !$omp target exit data map(delete: hArea_u, hArea_v)
  !$omp target exit data map(delete: rel_vort)
  !$omp target exit data map(delete: dvSdx, duSdy, stk_vort) if (Stokes_VF)
  call bxUx%free() ; call bxVx%free()
  call bxH_mij%free() ; call bxQ_gij%free() ; call bxQs_pi%free() ; call bxQs_pj%free()
end subroutine CorAdv_setup_fortran

!> Runs CorAdv_setup_fortran, a capture of its inputs and outputs, or the C++/AMReX bridge,
!! as selected by the CORADV_SETUP_MODE environment variable.
subroutine CorAdv_setup(bxH, bxQ, bxU, bxV, bxQs, u_a, v_a, h_a, OBC, AD_a, por_face_areaU_a, &
                        por_face_areaV_a, Waves, Gcore, Area_h_a, Area_q_a, Stokes_VF, use_weno, &
                        vol_neglect, area_neglect, q_a, Ih_q_a, abs_vort_a, h_q_a, qS_a, q2_a, &
                        uh_center_a, vh_center_a, RV_a, PV_a, CS)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxQ !< The B-grid-index iteration box of this tile,
                                                !! [IscB:IecB, JscB:JecB, ksc:kec]
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(Box_t),             intent(in)    :: bxQs !< The scheme-dependent vorticity-point box,
                                                 !! [Is_q:Ie_q, Js_q:Je_q, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a         !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a         !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: h_a         !< Layer thickness [H ~> m or kg m-2]
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary control structure
  type(CorAdv_AD_type),    intent(in)    :: AD_a !< Copies of the acceleration diagnostics
  type(RealArray_t), intent(in)    :: por_face_areaU_a !< fractional open area of U-faces [nondim]
  type(RealArray_t), intent(in)    :: por_face_areaV_a !< fractional open area of V-faces [nondim]
  type(Wave_parameters_CS), optional, pointer :: Waves !< An optional pointer to Stokes drift CS
  type(RealArray_t), intent(in)    :: Area_h_a    !< The ocean area at h points [L2 ~> m2]
  type(RealArray_t), intent(in)    :: Area_q_a    !< The sum of the ocean areas at the 4
                                                  !! adjacent thickness points [L2 ~> m2]
  logical,                 intent(in)    :: Stokes_VF !< If true, include the Stokes drift
  logical,                 intent(in)    :: use_weno  !< True if using one of the WENO schemes
  real,                    intent(in)    :: vol_neglect !< A negligible volume [H L2 ~> m3 or kg]
  real,                    intent(in)    :: area_neglect !< A negligible area [L2 ~> m2]
  type(RealArray_t), intent(inout) :: q_a         !< Layer potential vorticity
                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(inout) :: Ih_q_a      !< Inverse of thickness interpolated to q points
                                                  !! [H-1 ~> m-1 or m2 kg-1]
  type(RealArray_t), intent(inout) :: abs_vort_a  !< Absolute vorticity at q-points [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: h_q_a       !< Thickness interpolated to q points
                                                  !! [H ~> m or kg m-2]
  type(RealArray_t), intent(inout) :: qS_a        !< Layer Stokes vorticity
                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(inout) :: q2_a        !< Relative vorticity over thickness
                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(inout) :: uh_center_a !< Transport based on arithmetic mean h at
                                                  !! u-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(inout) :: vh_center_a !< Transport based on arithmetic mean h at
                                                  !! v-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(inout) :: RV_a        !< Diagnostic relative vorticity [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: PV_a        !< Diagnostic potential vorticity
                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(CorAdv_opts_CS),    intent(in)    :: CS  !< Options used in the Coriolis and momentum
                                                !! advection calculations
  ! Local variables
  integer            :: mode, rc
  type(Box_C)        :: bxH_c, bxQ_c, bxU_c, bxV_c, bxQs_c
  type(RealArray_C)  :: u_c, v_c, h_c, por_face_areaU_c, por_face_areaV_c, Area_h_c, Area_q_c, &
    q_c, Ih_q_c, abs_vort_c, h_q_c, qS_c, q2_c, uh_center_c, vh_center_c, RV_c, PV_c
  type(c_ptr) :: OBC_c
  type(CorAdv_AD_C) :: AD_a_c
  type(c_ptr) :: Waves_c
  logical(c_bool) :: Stokes_VF_c
  logical(c_bool) :: use_weno_c
  type(CorAdv_opts_CS_C) :: CS_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "coradv_setup"
  mode = getenv_mode("CORADV_SETUP_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxH", bxH)
        call rec%add("_bxQ", bxQ)
        call rec%add("_bxU", bxU)
        call rec%add("_bxV", bxV)
        call rec%add("_bxQs", bxQs)
        call rec%add("_u", u_a)
        call rec%add("_v", v_a)
        call rec%add("_h", h_a)
        call record_CorAdv_AD(rec, "_AD", AD_a, "")
        call rec%add("_por_face_areaU", por_face_areaU_a)
        call rec%add("_por_face_areaV", por_face_areaV_a)
        call rec%add("_dyCv", Gcore%dyCv)
        call rec%add("_dxCu", Gcore%dxCu)
        call rec%add("_dy_Cu", Gcore%dy_Cu)
        call rec%add("_dx_Cv", Gcore%dx_Cv)
        call rec%add("_dyBu", Gcore%dyBu)
        call rec%add("_dxBu", Gcore%dxBu)
        call rec%add("_mask2dBu", Gcore%mask2dBu)
        call rec%add("_IareaBu", Gcore%IareaBu)
        call rec%add("_CoriolisBu", Gcore%CoriolisBu)
        call rec%add("_Area_h", Area_h_a)
        call rec%add("_Area_q", Area_q_a)
        call rec%add("_Stokes_VF", Stokes_VF)
        call rec%add("_use_weno", use_weno)
        call rec%add("_vol_neglect", vol_neglect)
        call rec%add("_area_neglect", area_neglect)
        call rec%add("_q_before", q_a)
        call rec%add("_Ih_q_before", Ih_q_a)
        call rec%add("_abs_vort_before", abs_vort_a)
        call rec%add("_h_q_before", h_q_a)
        call rec%add("_qS_before", qS_a)
        call rec%add("_q2_before", q2_a)
        call rec%add("_uh_center_before", uh_center_a)
        call rec%add("_vh_center_before", vh_center_a)
        call rec%add("_RV_before", RV_a)
        call rec%add("_PV_before", PV_a)
        call record_CorAdv_opts_CS(rec, "_CS", CS)
      endif

      call CorAdv_setup_fortran(bxH, bxQ, bxU, bxV, bxQs, u_a, v_a, h_a, OBC, AD_a, &
                                por_face_areaU_a, por_face_areaV_a, Waves, Gcore, Area_h_a, &
                                Area_q_a, Stokes_VF, use_weno, vol_neglect, area_neglect, q_a, &
                                Ih_q_a, abs_vort_a, h_q_a, qS_a, q2_a, uh_center_a, vh_center_a, &
                                RV_a, PV_a, CS)

      if (capture) then
        call rec%add("_q_after", q_a)
        call rec%add("_Ih_q_after", Ih_q_a)
        call rec%add("_abs_vort_after", abs_vort_a)
        call rec%add("_h_q_after", h_q_a)
        call rec%add("_qS_after", qS_a)
        call rec%add("_q2_after", q2_a)
        call rec%add("_uh_center_after", uh_center_a)
        call rec%add("_vh_center_after", vh_center_a)
        call rec%add("_RV_after", RV_a)
        call rec%add("_PV_after", PV_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxH_c = bxH%to_c()
      bxQ_c = bxQ%to_c()
      bxU_c = bxU%to_c()
      bxV_c = bxV%to_c()
      bxQs_c = bxQs%to_c()
      u_c = u_a%to_c()
      v_c = v_a%to_c()
      h_c = h_a%to_c()
      if (associated(OBC)) then
        OBC_c = c_loc(OBC)
      else
        OBC_c = c_null_ptr
      endif
      AD_a_c = CorAdv_AD_to_c(AD_a)
      por_face_areaU_c = por_face_areaU_a%to_c()
      por_face_areaV_c = por_face_areaV_a%to_c()
      if (Stokes_VF) then
        Waves_c = c_loc(Waves)
      else
        Waves_c = c_null_ptr
      endif
      Area_h_c = Area_h_a%to_c()
      Area_q_c = Area_q_a%to_c()
      Stokes_VF_c = Stokes_VF
      use_weno_c = use_weno
      q_c = q_a%to_c()
      Ih_q_c = Ih_q_a%to_c()
      abs_vort_c = abs_vort_a%to_c()
      h_q_c = h_q_a%to_c()
      qS_c = qS_a%to_c()
      q2_c = q2_a%to_c()
      uh_center_c = uh_center_a%to_c()
      vh_center_c = vh_center_a%to_c()
      RV_c = RV_a%to_c()
      PV_c = PV_a%to_c()
      CS_c = CorAdv_opts_CS_to_c(CS)
      call turbotmp_coradv_setup_bridge(bxH_c, bxQ_c, bxU_c, bxV_c, bxQs_c, u_c, v_c, h_c, OBC_c, &
                                        AD_a_c, por_face_areaU_c, por_face_areaV_c, Waves_c, &
                                        Gcore%C, Area_h_c, Area_q_c, Stokes_VF_c, use_weno_c, &
                                        vol_neglect, area_neglect, q_c, Ih_q_c, abs_vort_c, h_q_c, &
                                        qS_c, q2_c, uh_center_c, vh_center_c, RV_c, PV_c, CS_c)
#endif

    case default
      call CorAdv_setup_fortran(bxH, bxQ, bxU, bxV, bxQs, u_a, v_a, h_a, OBC, AD_a, &
                                por_face_areaU_a, por_face_areaV_a, Waves, Gcore, Area_h_a, &
                                Area_q_a, Stokes_VF, use_weno, vol_neglect, area_neglect, q_a, &
                                Ih_q_a, abs_vort_a, h_q_a, qS_a, q2_a, uh_center_a, vh_center_a, &
                                RV_a, PV_a, CS)

  end select

end subroutine CorAdv_setup


!> Calculates the Coriolis and momentum advection accelerations with the Sadourny (1975)
!! energy- or enstrophy-conserving schemes (SADOURNY75_ENERGY, SADOURNY75_ENSTRO), over one tile.
subroutine CorAdv_sadourny_fortran(bxH, bxU, bxV, u_a, v_a, uh_a, vh_a, q_a, uh_center_a, &
                                   vh_center_a, Gcore, CAu_a, CAv_a, CS)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a         !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a         !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: uh_a        !< Zonal transport u*h*dy
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_a        !< Meridional transport v*h*dx
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: q_a         !< Layer potential vorticity
                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(in)    :: uh_center_a !< Transport based on arithmetic mean h at
                                                  !! u-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_center_a !< Transport based on arithmetic mean h at
                                                  !! v-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(inout) :: CAu_a       !< Zonal acceleration due to Coriolis
                                                  !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAv_a       !< Meridional acceleration due to Coriolis
                                                  !! and momentum advection [L T-2 ~> m s-2].
  type(Coriolis_scheme_CS), intent(in)  :: CS  !< Options selecting the Coriolis discretization

  ! Local variables
  type(Box_t) :: bxUx  ! bxU grown by one point at the end in i and the start in j,
                       ! [IscB:IecB+1, jsc-1:jec, ksc:kec]
  type(Box_t) :: bxVx  ! bxV grown by one point at the start in i and the end in j,
                       ! [isc-1:iec, JscB:JecB+1, ksc:kec]
  real, dimension(:,:,:), contiguous, pointer :: u, v, uh, vh, q, uh_center, vh_center, CAu, CAv
  real, dimension(:,:), contiguous, pointer :: dy_Cu, dx_Cv, IdxCu, IdyCv
  real, dimension(v_a%lb(1):v_a%ub(1),u_a%lb(2):u_a%ub(2),bxH%idxS(3):bxH%idxE(3)) :: &
    uh_min, uh_max, &   ! The smallest and largest estimates of the zonal volume fluxes through
                        ! the faces (i.e. u*h*dy) [H L2 T-1 ~> m3 s-1 or kg s-1]
    vh_min, vh_max      ! The smallest and largest estimates of the meridional volume fluxes through
                        ! the faces (i.e. v*h*dx) [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: temp1, temp2           ! Temporary variables [L2 T-2 ~> m2 s-2].
  real :: uhc, vhc               ! Centered estimates of uh and vh [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: uhm, vhm               ! The input estimates of uh and vh [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: c1, c2, c3, slope      ! Nondimensional parameters for the Coriolis limiter scheme [nondim]
  integer :: i, j, k

  call u_a%view(u) ; call v_a%view(v) ; call uh_a%view(uh)
  call vh_a%view(vh) ; call q_a%view(q) ; call uh_center_a%view(uh_center)
  call vh_center_a%view(vh_center) ; call CAu_a%view(CAu) ; call CAv_a%view(CAv)
  call Gcore%dy_Cu%view(dy_Cu) ; call Gcore%dx_Cv%view(dx_Cv) ; call Gcore%IdxCu%view(IdxCu)
  call Gcore%IdyCv%view(IdyCv)

  bxUx = bxU%growBy(lo=[0,1,0], hi=[1,0,0])
  bxVx = bxV%growBy(lo=[1,0,0], hi=[0,1,0])

  ! TODO: May also need SADOURNEY75_ENERGY
  !$omp target enter data map(alloc: uh_min, vh_min) if (CS%Coriolis_En_Dis)
  !$omp target enter data map(alloc: uh_max, vh_max) if (CS%Coriolis_En_Dis)

  ! The energy-dissipating biased limiter (Coriolis_En_Dis) is only used by these schemes.
  ! .and. SADOURNEY75_ENERGY ??
  if (CS%Coriolis_En_Dis) then
  !  c1 = 1.0-1.5*RANGE ; c2 = 1.0-RANGE ; c3 = 2.0 ; slope = 0.5
    c1 = 1.0-1.5*0.5 ; c2 = 1.0-0.5 ; c3 = 2.0 ; slope = 0.5

    do concurrent (k=bxVx%idxS(3):bxVx%idxE(3), j=bxVx%idxS(2):bxVx%idxE(2), &
                   I=bxVx%idxS(1):bxVx%idxE(1)) DO_LOCALITY(local(uhc, uhm))
      uhc = uh_center(I,j,k)
      uhm = uh(I,j,k)
      ! This sometimes matters with some types of open boundary conditions.
      if (dy_Cu(I,j) == 0.0) uhc = uhm

      if (abs(uhc) < 0.1*abs(uhm)) then
        uhm = 10.0*uhc
      elseif (abs(uhc) > c1*abs(uhm)) then
        if (abs(uhc) < c2*abs(uhm)) then ; uhc = (3.0*uhc+(1.0-c2*3.0)*uhm)
        elseif (abs(uhc) <= c3*abs(uhm)) then ; uhc = uhm
        else ; uhc = slope*uhc+(1.0-c3*slope)*uhm
        endif
      endif

      if (uhc > uhm) then
        uh_min(I,j,k) = uhm ; uh_max(I,j,k) = uhc
      else
        uh_max(I,j,k) = uhm ; uh_min(I,j,k) = uhc
      endif
    enddo

    do concurrent (k=bxUx%idxS(3):bxUx%idxE(3), J=bxUx%idxS(2):bxUx%idxE(2), &
                   i=bxUx%idxS(1):bxUx%idxE(1)) DO_LOCALITY(local(vhc, vhm))
      vhc = vh_center(i,J,k)
      vhm = vh(i,J,k)
      ! This sometimes matters with some types of open boundary conditions.
      if (dx_Cv(i,J) == 0.0) vhc = vhm

      if (abs(vhc) < 0.1*abs(vhm)) then
        vhm = 10.0*vhc
      elseif (abs(vhc) > c1*abs(vhm)) then
        if (abs(vhc) < c2*abs(vhm)) then ; vhc = (3.0*vhc+(1.0-c2*3.0)*vhm)
        elseif (abs(vhc) <= c3*abs(vhm)) then ; vhc = vhm
        else ; vhc = slope*vhc+(1.0-c3*slope)*vhm
        endif
      endif

      if (vhc > vhm) then
        vh_min(i,J,k) = vhm ; vh_max(i,J,k) = vhc
      else
        vh_max(i,J,k) = vhm ; vh_min(i,J,k) = vhc
      endif
    enddo
  endif

  ! Calculate the tendencies of zonal velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAu =  q * vh - d(KE)/dx.
  if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
    if (CS%Coriolis_En_Dis) then
      ! Energy dissipating biased scheme, Hallberg 200x
      do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                     I=bxU%idxS(1):bxU%idxE(1)) DO_LOCALITY(local(temp1, temp2))
        if (q(I,J,k)*u(I,j,k) == 0.0) then
          temp1 = q(I,J,k) * ( (vh_max(i,j,k)+vh_max(i+1,j,k)) &
                           + (vh_min(i,j,k)+vh_min(i+1,j,k)) )*0.5
        elseif (q(I,J,k)*u(I,j,k) < 0.0) then
          temp1 = q(I,J,k) * (vh_max(i,j,k)+vh_max(i+1,j,k))
        else
          temp1 = q(I,J,k) * (vh_min(i,j,k)+vh_min(i+1,j,k))
        endif
        if (q(I,J-1,k)*u(I,j,k) == 0.0) then
          temp2 = q(I,J-1,k) * ( (vh_max(i,j-1,k)+vh_max(i+1,j-1,k)) &
                             + (vh_min(i,j-1,k)+vh_min(i+1,j-1,k)) )*0.5
        elseif (q(I,J-1,k)*u(I,j,k) < 0.0) then
          temp2 = q(I,J-1,k) * (vh_max(i,j-1,k)+vh_max(i+1,j-1,k))
        else
          temp2 = q(I,J-1,k) * (vh_min(i,j-1,k)+vh_min(i+1,j-1,k))
        endif
        CAu(I,j,k) = 0.25 * IdxCu(I,j) * (temp1 + temp2)
      enddo
    else
      ! Energy conserving scheme, Sadourny 1975
      do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                     I=bxU%idxS(1):bxU%idxE(1))
        CAu(I,j,k) = 0.25 * &
          ((q(I,J,k) * (vh(i+1,J,k) + vh(i,J,k))) + &
           (q(I,J-1,k) * (vh(i,J-1,k) + vh(i+1,J-1,k)))) * IdxCu(I,j)
      enddo
    endif
  elseif (CS%Coriolis_Scheme == SADOURNY75_ENSTRO) then
    do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                   I=bxU%idxS(1):bxU%idxE(1))
      CAu(I,j,k) = 0.125 * (IdxCu(I,j) * (q(I,J,k) + q(I,J-1,k))) * &
                   ((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
    enddo
  endif

  ! Calculate the tendencies of meridional velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAv = - q * uh - d(KE)/dy.
  if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
    if (CS%Coriolis_En_Dis) then
      ! Energy dissipating biased scheme, Hallberg 200x
      do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                     i=bxV%idxS(1):bxV%idxE(1)) DO_LOCALITY(local(temp1, temp2))
        if (q(I-1,J,k)*v(i,J,k) == 0.0) then
          temp1 = q(I-1,J,k) * ( (uh_max(i-1,j,k)+uh_max(i-1,j+1,k)) &
                             + (uh_min(i-1,j,k)+uh_min(i-1,j+1,k)) )*0.5
        elseif (q(I-1,J,k)*v(i,J,k) > 0.0) then
          temp1 = q(I-1,J,k) * (uh_max(i-1,j,k)+uh_max(i-1,j+1,k))
        else
          temp1 = q(I-1,J,k) * (uh_min(i-1,j,k)+uh_min(i-1,j+1,k))
        endif
        if (q(I,J,k)*v(i,J,k) == 0.0) then
          temp2 = q(I,J,k) * ( (uh_max(i,j,k)+uh_max(i,j+1,k)) &
                           + (uh_min(i,j,k)+uh_min(i,j+1,k)) )*0.5
        elseif (q(I,J,k)*v(i,J,k) > 0.0) then
          temp2 = q(I,J,k) * (uh_max(i,j,k)+uh_max(i,j+1,k))
        else
          temp2 = q(I,J,k) * (uh_min(i,j,k)+uh_min(i,j+1,k))
        endif
        CAv(i,J,k) = -0.25 * IdyCv(i,J) * (temp1 + temp2)
      enddo
    else
      ! Energy conserving scheme, Sadourny 1975
      do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                     i=bxV%idxS(1):bxV%idxE(1))
        CAv(i,J,k) = - 0.25* &
            ((q(I-1,J,k)*(uh(I-1,j,k) + uh(I-1,j+1,k))) + &
             (q(I,J,k)*(uh(I,j,k) + uh(I,j+1,k)))) * IdyCv(i,J)
      enddo
    endif
  elseif (CS%Coriolis_Scheme == SADOURNY75_ENSTRO) then
    do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                   i=bxV%idxS(1):bxV%idxE(1))
      CAv(i,J,k) = -0.125 * (IdyCv(i,J) * (q(I-1,J,k) + q(I,J,k))) * &
                   ((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))
    enddo
  endif

  !$omp target exit data map(delete: uh_min, vh_min) if (CS%Coriolis_En_Dis)
  !$omp target exit data map(delete: uh_max, vh_max) if (CS%Coriolis_En_Dis)
  call bxUx%free() ; call bxVx%free()
end subroutine CorAdv_sadourny_fortran

!> Runs CorAdv_sadourny_fortran, a capture of its inputs and outputs, or the C++/AMReX bridge,
!! as selected by the CORADV_SADOURNY_MODE environment variable.
subroutine CorAdv_sadourny(bxH, bxU, bxV, u_a, v_a, uh_a, vh_a, q_a, uh_center_a, vh_center_a, &
                           Gcore, CAu_a, CAv_a, CS)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a         !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a         !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: uh_a        !< Zonal transport u*h*dy
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_a        !< Meridional transport v*h*dx
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: q_a         !< Layer potential vorticity
                                                  !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(in)    :: uh_center_a !< Transport based on arithmetic mean h at
                                                  !! u-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_center_a !< Transport based on arithmetic mean h at
                                                  !! v-points
                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(inout) :: CAu_a       !< Zonal acceleration due to Coriolis
                                                  !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAv_a       !< Meridional acceleration due to Coriolis
                                                  !! and momentum advection [L T-2 ~> m s-2].
  type(Coriolis_scheme_CS), intent(in)  :: CS  !< Options selecting the Coriolis discretization
  ! Local variables
  integer            :: mode, rc
  type(Box_C)        :: bxH_c, bxU_c, bxV_c
  type(RealArray_C)  :: u_c, v_c, uh_c, vh_c, q_c, uh_center_c, vh_center_c, CAu_c, CAv_c
  type(Coriolis_scheme_CS_C) :: CS_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "coradv_sadourny"
  mode = getenv_mode("CORADV_SADOURNY_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxH", bxH)
        call rec%add("_bxU", bxU)
        call rec%add("_bxV", bxV)
        call rec%add("_u", u_a)
        call rec%add("_v", v_a)
        call rec%add("_uh", uh_a)
        call rec%add("_vh", vh_a)
        call rec%add("_q", q_a)
        call rec%add("_uh_center", uh_center_a)
        call rec%add("_vh_center", vh_center_a)
        call rec%add("_dy_Cu", Gcore%dy_Cu)
        call rec%add("_dx_Cv", Gcore%dx_Cv)
        call rec%add("_IdxCu", Gcore%IdxCu)
        call rec%add("_IdyCv", Gcore%IdyCv)
        call rec%add("_CAu_before", CAu_a)
        call rec%add("_CAv_before", CAv_a)
        call record_Coriolis_scheme_CS(rec, "_CS", CS)
      endif

      call CorAdv_sadourny_fortran(bxH, bxU, bxV, u_a, v_a, uh_a, vh_a, q_a, uh_center_a, &
                                   vh_center_a, Gcore, CAu_a, CAv_a, CS)

      if (capture) then
        call rec%add("_CAu_after", CAu_a)
        call rec%add("_CAv_after", CAv_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxH_c = bxH%to_c()
      bxU_c = bxU%to_c()
      bxV_c = bxV%to_c()
      u_c = u_a%to_c()
      v_c = v_a%to_c()
      uh_c = uh_a%to_c()
      vh_c = vh_a%to_c()
      q_c = q_a%to_c()
      uh_center_c = uh_center_a%to_c()
      vh_center_c = vh_center_a%to_c()
      CAu_c = CAu_a%to_c()
      CAv_c = CAv_a%to_c()
      CS_c = Coriolis_scheme_CS_to_c(CS)
      call turbotmp_coradv_sadourny_bridge(bxH_c, bxU_c, bxV_c, u_c, v_c, uh_c, vh_c, q_c, &
                                           uh_center_c, vh_center_c, Gcore%C, CAu_c, CAv_c, CS_c)
#endif

    case default
      call CorAdv_sadourny_fortran(bxH, bxU, bxV, u_a, v_a, uh_a, vh_a, q_a, uh_center_a, &
                                   vh_center_a, Gcore, CAu_a, CAv_a, CS)

  end select

end subroutine CorAdv_sadourny


!> Calculates the Coriolis and momentum advection accelerations with the Arakawa & Hsu (1990),
!! Arakawa & Lamb (1981) or blended (ARAKAWA_LAMB_BLEND) schemes, over one tile.
subroutine CorAdv_arakawa_fortran(bxH, bxQ, bxU, bxV, uh_a, vh_a, q_a, Ih_q_a, Fe_m2, rat_lin, &
                                  Gcore, CAu_a, CAv_a, wt_lin_blend, CS)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxQ !< The B-grid-index iteration box of this tile,
                                                !! [IscB:IecB, JscB:JecB, ksc:kec]
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: uh_a   !< Zonal transport u*h*dy
                                             !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_a   !< Meridional transport v*h*dx
                                             !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: q_a    !< Layer potential vorticity
                                             !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(in)    :: Ih_q_a !< Inverse of thickness interpolated to q points
                                             !! [H-1 ~> m-1 or m2 kg-1]
  real,                    intent(in)    :: Fe_m2   !< ARAKAWA_LAMB_BLEND blending parameter [nondim]
  real,                    intent(in)    :: rat_lin !< ARAKAWA_LAMB_BLEND blending parameter [nondim]
  type(RealArray_t), intent(inout) :: CAu_a  !< Zonal acceleration due to Coriolis
                                             !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAv_a  !< Meridional acceleration due to Coriolis
                                             !! and momentum advection [L T-2 ~> m s-2].
  real,                    intent(in)    :: wt_lin_blend !< A weighting value beyond which the
                                          !! blending between Sadourny and Arakawa & Hsu goes
                                          !! linearly to 0 [nondim].
  type(Coriolis_scheme_CS), intent(in)  :: CS  !< Options selecting the Coriolis discretization

  ! Local variables
  type(Box_t) :: bxQ_pij  ! bxQ grown by one point at the end in i and j,
                          ! [IscB:IecB+1, JscB:JecB+1, ksc:kec]
  type(Box_t) :: bxQ_pj   ! bxQ grown by one point at the end in j,
                          ! [IscB:IecB, JscB:JecB+1, ksc:kec]
  type(Box_t) :: bxVx  ! bxV grown by one point at the start in i and the end in j,
                       ! [isc-1:iec, JscB:JecB+1, ksc:kec]
  real, dimension(:,:,:), contiguous, pointer :: uh, vh, q, Ih_q, CAu, CAv
  real, dimension(:,:), contiguous, pointer :: IdxCu, IdyCv
  real, dimension(uh_a%lb(1):uh_a%ub(1),uh_a%lb(2):uh_a%ub(2),bxH%idxS(3):bxH%idxE(3)) :: &
    a, b, c, d    ! a, b, c, & d are combinations of the potential vorticities
                  ! surrounding an h grid point.  At small scales, a = q/4,
                  ! b = q/4, etc.  All are in [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1],
                  ! and use the indexing of the corresponding u point.
  real, dimension(vh_a%lb(1):vh_a%ub(1),uh_a%lb(2):uh_a%ub(2),bxH%idxS(3):bxH%idxE(3)) :: &
    ep_u, ep_v    ! Additional pseudo-Coriolis terms in the Arakawa and Lamb
                  ! discretization [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1].
  real, parameter :: C1_12 = 1.0 / 12.0 ! C1_12 = 1/12 [nondim]
  real, parameter :: C1_24 = 1.0 / 24.0 ! C1_24 = 1/24 [nondim]
  real :: max_Ihq, min_Ihq       ! The maximum and minimum of the nearby Ihq [H-1 ~> m-1 or m2 kg-1].
  real :: rat_m1        ! The ratio of the maximum neighboring inverse thickness
                        ! to the minimum inverse thickness minus 1 [nondim]. rat_m1 >= 0.
  real :: AL_wt         ! The relative weight of the Arakawa & Lamb scheme to the
                        ! Arakawa & Hsu scheme [nondim], between 0 and 1.
  real :: Sad_wt        ! The relative weight of the Sadourny energy scheme to
                        ! the other two with the ARAKAWA_LAMB_BLEND scheme [nondim],
                        ! between 0 and 1.
  integer :: i, j, k

  call uh_a%view(uh) ; call vh_a%view(vh) ; call q_a%view(q)
  call Ih_q_a%view(Ih_q) ; call CAu_a%view(CAu) ; call CAv_a%view(CAv)
  call Gcore%IdxCu%view(IdxCu) ; call Gcore%IdyCv%view(IdyCv)

  bxVx = bxV%growBy(lo=[1,0,0], hi=[0,1,0])
  bxQ_pij = bxQ%growBy(lo=[0,0,0], hi=[1,1,0])
  bxQ_pj = bxQ%growHi(dim=2, n=1)

  !$omp target enter data map(alloc: a, b, c, d, ep_u, ep_v)

  !   a, b, c, and d are combinations of neighboring potential
  ! vorticities which form the Arakawa and Hsu vorticity advection
  ! scheme.  All are defined at u grid points.

  if (CS%Coriolis_Scheme == ARAKAWA_HSU90) then
    do concurrent (k=bxVx%idxS(3):bxVx%idxE(3), j=bxVx%idxS(2):bxVx%idxE(2), &
                   I=bxVx%idxS(1):bxVx%idxE(1))
      a(I,j,k) = (q(I,J,k) + (q(I+1,J,k) + q(I,J-1,k))) * C1_12
      d(I,j,k) = ((q(I,J,k) + q(I+1,J-1,k)) + q(I,J-1,k)) * C1_12
    enddo

    do concurrent (k=bxQ_pj%idxS(3):bxQ_pj%idxE(3), j=bxQ_pj%idxS(2):bxQ_pj%idxE(2), &
                   I=bxQ_pj%idxS(1):bxQ_pj%idxE(1))
      b(I,j,k) = (q(I,J,k) + (q(I-1,J,k) + q(I,J-1,k))) * C1_12
      c(I,j,k) = ((q(I,J,k) + q(I-1,J-1,k)) + q(I,J-1,k)) * C1_12
    enddo
  elseif (CS%Coriolis_Scheme == ARAKAWA_LAMB81) then
    do concurrent (k=bxQ_pij%idxS(3):bxQ_pij%idxE(3), j=bxQ_pij%idxS(2):bxQ_pij%idxE(2), &
                   I=bxQ_pij%idxS(1):bxQ_pij%idxE(1))
      a(I-1,j,k) = (2.0*(q(I,J,k) + q(I-1,J-1,k)) + (q(I-1,J,k) + q(I,J-1,k))) * C1_24
      d(I-1,j,k) = ((q(I,j,k) + q(I-1,J-1,k)) + 2.0*(q(I-1,J,k) + q(I,J-1,k))) * C1_24
      b(I,j,k) =   ((q(I,J,k) + q(I-1,J-1,k)) + 2.0*(q(I-1,J,k) + q(I,J-1,k))) * C1_24
      c(I,j,k) =   (2.0*(q(I,J,k) + q(I-1,J-1,k)) + (q(I-1,J,k) + q(I,J-1,k))) * C1_24
      ep_u(i,j,k) = ((q(I,J,k) - q(I-1,J-1,k)) + (q(I-1,J,k) - q(I,J-1,k))) * C1_24
      ep_v(i,j,k) = (-(q(I,J,k) - q(I-1,J-1,k)) + (q(I-1,J,k) - q(I,J-1,k))) * C1_24
    enddo
  elseif (CS%Coriolis_Scheme == AL_BLEND) then
    ! Fe_m2 and rat_lin are k-independent; computed in CorAdCalc_TR before the tile loop.
    do concurrent (k=bxQ_pij%idxS(3):bxQ_pij%idxE(3), j=bxQ_pij%idxS(2):bxQ_pij%idxE(2), &
                   I=bxQ_pij%idxS(1):bxQ_pij%idxE(1)) &
        DO_LOCALITY(local(min_Ihq, max_Ihq, rat_m1, AL_wt, Sad_wt))
      min_Ihq = MIN(Ih_q(I-1,J-1,k), Ih_q(I,J-1,k), Ih_q(I-1,J,k), Ih_q(I,J,k))
      max_Ihq = MAX(Ih_q(I-1,J-1,k), Ih_q(I,J-1,k), Ih_q(I-1,J,k), Ih_q(I,J,k))
      rat_m1 = 1.0e15
      if (max_Ihq < 1.0e15*min_Ihq) rat_m1 = max_Ihq / min_Ihq - 1.0
      ! The weights used here are designed to keep the effective Coriolis
      ! acceleration from any one point on its neighbors within a factor
      ! of F_eff_max.  The minimum permitted value is 2 (the factor for
      ! Sadourny's energy conserving scheme).

      ! Determine the relative weights of Arakawa & Lamb vs. Arakawa and Hsu.
      if (rat_m1 <= Fe_m2) then ; AL_wt = 1.0
      elseif (rat_m1 < 1.5*Fe_m2) then ; AL_wt = 3.0*Fe_m2 / rat_m1 - 2.0
      else ; AL_wt = 0.0 ; endif

      ! Determine the relative weights of Sadourny Energy vs. the other two.
      if (rat_m1 <= 1.5*Fe_m2) then ; Sad_wt = 0.0
      elseif (rat_m1 <= rat_lin) then
        Sad_wt = 1.0 - (1.5*Fe_m2) / rat_m1
      elseif (rat_m1 < 2.0*rat_lin) then
        Sad_wt = 1.0 - (wt_lin_blend / rat_lin) * (rat_m1 - 2.0*rat_lin)
      else ; Sad_wt = 1.0 ; endif

      a(I-1,j,k) = Sad_wt * 0.25 * q(I-1,J,k) + (1.0 - Sad_wt) * &
                    ( ((2.0-AL_wt)* q(I-1,J,k) + AL_wt*q(I,J-1,k)) + &
                       2.0 * (q(I,J,k) + q(I-1,J-1,k)) ) * C1_24
      d(I-1,j,k) = Sad_wt * 0.25 * q(I-1,J-1,k) + (1.0 - Sad_wt) * &
                    ( ((2.0-AL_wt)* q(I-1,J-1,k) + AL_wt*q(I,J,k)) + &
                       2.0 * (q(I-1,J,k) + q(I,J-1,k)) ) * C1_24
      b(I,j,k) =   Sad_wt * 0.25 * q(I,J,k) + (1.0 - Sad_wt) * &
                    ( ((2.0-AL_wt)* q(I,J,k) + AL_wt*q(I-1,J-1,k)) + &
                       2.0 * (q(I-1,J,k) + q(I,J-1,k)) ) * C1_24
      c(I,j,k) =   Sad_wt * 0.25 * q(I,J-1,k) + (1.0 - Sad_wt) * &
                    ( ((2.0-AL_wt)* q(I,J-1,k) + AL_wt*q(I-1,J,k)) + &
                       2.0 * (q(I,J,k) + q(I-1,J-1,k)) ) * C1_24
      ep_u(i,j,k) = AL_wt  * ((q(I,J,k) - q(I-1,J-1,k)) + (q(I-1,J,k) - q(I,J-1,k))) * C1_24
      ep_v(i,j,k) = AL_wt * (-(q(I,J,k) - q(I-1,J-1,k)) + (q(I-1,J,k) - q(I,J-1,k))) * C1_24
    enddo
  endif

  ! Calculate the tendencies of zonal velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAu =  q * vh - d(KE)/dx.
    ! (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
    do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                   I=bxU%idxS(1):bxU%idxE(1))
      CAu(I,j,k) = (((a(I,j,k) * vh(i+1,J,k)) +  (c(I,j,k) * vh(i,J-1,k)))  + &
                    ((b(I,j,k) * vh(i,J,k)) +  (d(I,j,k) * vh(i+1,J-1,k)))) * IdxCu(I,j)
    enddo

  ! Add in the additional terms with Arakawa & Lamb.
  if ((CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
      (CS%Coriolis_Scheme == AL_BLEND)) then
    do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                   I=bxU%idxS(1):bxU%idxE(1))
      CAu(I,j,k) = CAu(I,j,k) + &
            ((ep_u(i,j,k)*uh(I-1,j,k)) - (ep_u(i+1,j,k)*uh(I+1,j,k))) * IdxCu(I,j)
    enddo
  endif

  ! Calculate the tendencies of meridional velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAv = - q * uh - d(KE)/dy.
    ! (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
    do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                   i=bxV%idxS(1):bxV%idxE(1))
      CAv(i,J,k) = - (((a(I-1,j,k)   * uh(I-1,j,k)) + &
                       (c(I,j+1,k)   * uh(I,j+1,k)))  &
                    + ((b(I,j,k)     * uh(I,j,k)) +   &
                       (d(I-1,j+1,k) * uh(I-1,j+1,k)))) * IdyCv(i,J)
    enddo
  ! Add in the additonal terms with Arakawa & Lamb.
  if ((CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
      (CS%Coriolis_Scheme == AL_BLEND)) then
    do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                   i=bxV%idxS(1):bxV%idxE(1))
      CAv(i,J,k) = CAv(i,J,k) + &
            ((ep_v(i,j,k)*vh(i,J-1,k)) - (ep_v(i,j+1,k)*vh(i,J+1,k))) * IdyCv(i,J)
    enddo
  endif

  !$omp target exit data map(delete: a, b, c, d, ep_u, ep_v)
  call bxVx%free()
  call bxQ_pij%free() ; call bxQ_pj%free()
end subroutine CorAdv_arakawa_fortran

!> Runs CorAdv_arakawa_fortran, a capture of its inputs and outputs, or the C++/AMReX bridge,
!! as selected by the CORADV_ARAKAWA_MODE environment variable.
subroutine CorAdv_arakawa(bxH, bxQ, bxU, bxV, uh_a, vh_a, q_a, Ih_q_a, Fe_m2, rat_lin, Gcore, &
                          CAu_a, CAv_a, wt_lin_blend, CS)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),             intent(in)    :: bxH !< The h-point iteration box of this tile,
                                                !! [isc:iec, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxQ !< The B-grid-index iteration box of this tile,
                                                !! [IscB:IecB, JscB:JecB, ksc:kec]
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: uh_a   !< Zonal transport u*h*dy
                                             !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_a   !< Meridional transport v*h*dx
                                             !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: q_a    !< Layer potential vorticity
                                             !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(in)    :: Ih_q_a !< Inverse of thickness interpolated to q points
                                             !! [H-1 ~> m-1 or m2 kg-1]
  real,                    intent(in)    :: Fe_m2   !< ARAKAWA_LAMB_BLEND blending parameter [nondim]
  real,                    intent(in)    :: rat_lin !< ARAKAWA_LAMB_BLEND blending parameter [nondim]
  type(RealArray_t), intent(inout) :: CAu_a  !< Zonal acceleration due to Coriolis
                                             !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAv_a  !< Meridional acceleration due to Coriolis
                                             !! and momentum advection [L T-2 ~> m s-2].
  real,                    intent(in)    :: wt_lin_blend !< A weighting value beyond which the
                                          !! blending between Sadourny and Arakawa & Hsu goes
                                          !! linearly to 0 [nondim].
  type(Coriolis_scheme_CS), intent(in)  :: CS  !< Options selecting the Coriolis discretization
  ! Local variables
  integer            :: mode, rc
  type(Box_C)        :: bxH_c, bxQ_c, bxU_c, bxV_c
  type(RealArray_C)  :: uh_c, vh_c, q_c, Ih_q_c, CAu_c, CAv_c
  type(Coriolis_scheme_CS_C) :: CS_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "coradv_arakawa"
  mode = getenv_mode("CORADV_ARAKAWA_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxH", bxH)
        call rec%add("_bxQ", bxQ)
        call rec%add("_bxU", bxU)
        call rec%add("_bxV", bxV)
        call rec%add("_uh", uh_a)
        call rec%add("_vh", vh_a)
        call rec%add("_q", q_a)
        call rec%add("_Ih_q", Ih_q_a)
        call rec%add("_Fe_m2", Fe_m2)
        call rec%add("_rat_lin", rat_lin)
        call rec%add("_IdxCu", Gcore%IdxCu)
        call rec%add("_IdyCv", Gcore%IdyCv)
        call rec%add("_CAu_before", CAu_a)
        call rec%add("_CAv_before", CAv_a)
        call rec%add("_wt_lin_blend", wt_lin_blend)
        call record_Coriolis_scheme_CS(rec, "_CS", CS)
      endif

      call CorAdv_arakawa_fortran(bxH, bxQ, bxU, bxV, uh_a, vh_a, q_a, Ih_q_a, Fe_m2, rat_lin, &
                                  Gcore, CAu_a, CAv_a, wt_lin_blend, CS)

      if (capture) then
        call rec%add("_CAu_after", CAu_a)
        call rec%add("_CAv_after", CAv_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxH_c = bxH%to_c()
      bxQ_c = bxQ%to_c()
      bxU_c = bxU%to_c()
      bxV_c = bxV%to_c()
      uh_c = uh_a%to_c()
      vh_c = vh_a%to_c()
      q_c = q_a%to_c()
      Ih_q_c = Ih_q_a%to_c()
      CAu_c = CAu_a%to_c()
      CAv_c = CAv_a%to_c()
      CS_c = Coriolis_scheme_CS_to_c(CS)
      call turbotmp_coradv_arakawa_bridge(bxH_c, bxQ_c, bxU_c, bxV_c, uh_c, vh_c, q_c, Ih_q_c, &
                                          Fe_m2, rat_lin, Gcore%C, CAu_c, CAv_c, wt_lin_blend, &
                                          CS_c)
#endif

    case default
      call CorAdv_arakawa_fortran(bxH, bxQ, bxU, bxV, uh_a, vh_a, q_a, Ih_q_a, Fe_m2, rat_lin, &
                                  Gcore, CAu_a, CAv_a, wt_lin_blend, CS)

  end select

end subroutine CorAdv_arakawa


!> Calculates the Coriolis and momentum advection accelerations with the pseudo-enstrophy
!! scheme that is robust to vanishing layers (ROBUST_ENSTRO), over one tile.
subroutine CorAdv_robust_enstro_fortran(bxU, bxV, u_a, v_a, h_a, uh_a, vh_a, abs_vort_a, eps_vel, &
                                        h_tiny, Gcore, CAu_a, CAv_a, PV_Adv_Scheme)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a        !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a        !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: h_a        !< Layer thickness [H ~> m or kg m-2]
  type(RealArray_t), intent(in)    :: uh_a       !< Zonal transport u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_a       !< Meridional transport v*h*dx
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: abs_vort_a !< Absolute vorticity at q-points [T-1 ~> s-1]
  real,                    intent(in)    :: eps_vel !< A tiny, positive velocity [L T-1 ~> m s-1]
  real,                    intent(in)    :: h_tiny  !< A very small thickness [H ~> m or kg m-2]
  type(RealArray_t), intent(inout) :: CAu_a      !< Zonal acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAv_a      !< Meridional acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  integer,                 intent(in)    :: PV_Adv_Scheme !< Selects the discretization for PV
                                          !! advection: PV_ADV_CENTERED or PV_ADV_UPWIND1

  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: u, v, h, uh, vh, abs_vort, CAu, CAv
  real, dimension(:,:), contiguous, pointer :: IdxCv, IdxCu, IdyCu, IdyCv
  real :: Heff1, Heff2  ! Temporary effective H at U or V points [H ~> m or kg m-2].
  real :: Heff3, Heff4  ! Temporary effective H at U or V points [H ~> m or kg m-2].
  real :: UHeff, VHeff  ! More temporary variables [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: QUHeff,QVHeff ! More temporary variables [H L2 T-2 ~> m3 s-2 or kg s-2].
  integer :: i, j, k

  call u_a%view(u) ; call v_a%view(v) ; call h_a%view(h)
  call uh_a%view(uh) ; call vh_a%view(vh) ; call abs_vort_a%view(abs_vort)
  call CAu_a%view(CAu) ; call CAv_a%view(CAv)
  call Gcore%IdxCv%view(IdxCv) ; call Gcore%IdxCu%view(IdxCu) ; call Gcore%IdyCu%view(IdyCu)
  call Gcore%IdyCv%view(IdyCv)

  ! Calculate the tendencies of zonal velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAu =  q * vh - d(KE)/dx.
    ! An enstrophy conserving scheme robust to vanishing layers
    ! Note: Heffs are in lieu of h_at_v that should be returned by the
    !       continuity solver. AJA
    do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                   I=bxU%idxS(1):bxU%idxE(1)) &
        DO_LOCALITY(local(Heff1, Heff2, Heff3, Heff4, VHeff, QVHeff))
      Heff1 = abs(vh(i,J,k) * IdxCv(i,J)) / (eps_vel+abs(v(i,J,k)))
      Heff1 = max(Heff1, min(h(i,j,k),h(i,j+1,k)))
      Heff1 = min(Heff1, max(h(i,j,k),h(i,j+1,k)))
      Heff2 = abs(vh(i,J-1,k) * IdxCv(i,J-1)) / (eps_vel+abs(v(i,J-1,k)))
      Heff2 = max(Heff2, min(h(i,j-1,k),h(i,j,k)))
      Heff2 = min(Heff2, max(h(i,j-1,k),h(i,j,k)))
      Heff3 = abs(vh(i+1,J,k) * IdxCv(i+1,J)) / (eps_vel+abs(v(i+1,J,k)))
      Heff3 = max(Heff3, min(h(i+1,j,k),h(i+1,j+1,k)))
      Heff3 = min(Heff3, max(h(i+1,j,k),h(i+1,j+1,k)))
      Heff4 = abs(vh(i+1,J-1,k) * IdxCv(i+1,J-1)) / (eps_vel+abs(v(i+1,J-1,k)))
      Heff4 = max(Heff4, min(h(i+1,j-1,k),h(i+1,j,k)))
      Heff4 = min(Heff4, max(h(i+1,j-1,k),h(i+1,j,k)))
      if (PV_Adv_Scheme == PV_ADV_CENTERED) then
        CAu(I,j,k) = 0.5*(abs_vort(I,J,k)+abs_vort(I,J-1,k)) * &
                     ((vh(i,J,k) + vh(i+1,J-1,k)) + (vh(i,J-1,k) + vh(i+1,J,k)) ) /  &
                     (h_tiny + ((Heff1+Heff4) + (Heff2+Heff3)) ) * IdxCu(I,j)
      elseif (PV_Adv_Scheme == PV_ADV_UPWIND1) then
        VHeff = ((vh(i,J,k) + vh(i+1,J-1,k)) + (vh(i,J-1,k) + vh(i+1,J,k)) )
        QVHeff = 0.5*( ((abs_vort(I,J,k)+abs_vort(I,J-1,k))*VHeff) &
                     - ((abs_vort(I,J,k)-abs_vort(I,J-1,k))*abs(VHeff)) )
        CAu(I,j,k) = (QVHeff / ( h_tiny + ((Heff1+Heff4) + (Heff2+Heff3)) ) ) * IdxCu(I,j)
      endif
    enddo

  ! Calculate the tendencies of meridional velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAv = - q * uh - d(KE)/dy.
    ! An enstrophy conserving scheme robust to vanishing layers
    ! Note: Heffs are in lieu of h_at_u that should be returned by the
    !       continuity solver. AJA
    do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                   i=bxV%idxS(1):bxV%idxE(1)) &
        DO_LOCALITY(local(Heff1, Heff2, Heff3, Heff4, UHeff, QUHeff))
      Heff1 = abs(uh(I,j,k) * IdyCu(I,j)) / (eps_vel+abs(u(I,j,k)))
      Heff1 = max(Heff1, min(h(i,j,k),h(i+1,j,k)))
      Heff1 = min(Heff1, max(h(i,j,k),h(i+1,j,k)))
      Heff2 = abs(uh(I-1,j,k) * IdyCu(I-1,j)) / (eps_vel+abs(u(I-1,j,k)))
      Heff2 = max(Heff2, min(h(i-1,j,k),h(i,j,k)))
      Heff2 = min(Heff2, max(h(i-1,j,k),h(i,j,k)))
      Heff3 = abs(uh(I,j+1,k) * IdyCu(I,j+1)) / (eps_vel+abs(u(I,j+1,k)))
      Heff3 = max(Heff3, min(h(i,j+1,k),h(i+1,j+1,k)))
      Heff3 = min(Heff3, max(h(i,j+1,k),h(i+1,j+1,k)))
      Heff4 = abs(uh(I-1,j+1,k) * IdyCu(I-1,j+1)) / (eps_vel+abs(u(I-1,j+1,k)))
      Heff4 = max(Heff4, min(h(i-1,j+1,k),h(i,j+1,k)))
      Heff4 = min(Heff4, max(h(i-1,j+1,k),h(i,j+1,k)))
      if (PV_Adv_Scheme == PV_ADV_CENTERED) then
        CAv(i,J,k) = - 0.5*(abs_vort(I,J,k)+abs_vort(I-1,J,k)) * &
                       ((uh(I  ,j  ,k)+uh(I-1,j+1,k)) +      &
                        (uh(I-1,j  ,k)+uh(I  ,j+1,k)) ) /    &
                    (h_tiny + ((Heff1+Heff4) +(Heff2+Heff3)) ) * IdyCv(i,J)
      elseif (PV_Adv_Scheme == PV_ADV_UPWIND1) then
        UHeff = ((uh(I  ,j  ,k)+uh(I-1,j+1,k)) +      &
                 (uh(I-1,j  ,k)+uh(I  ,j+1,k)) )
        QUHeff = 0.5*( ((abs_vort(I,J,k)+abs_vort(I-1,J,k))*UHeff) &
                     - ((abs_vort(I,J,k)-abs_vort(I-1,J,k))*abs(UHeff)) )
        CAv(i,J,k) = - QUHeff / &
                     (h_tiny + ((Heff1+Heff4) +(Heff2+Heff3)) ) * IdyCv(i,J)
      endif
    enddo
end subroutine CorAdv_robust_enstro_fortran

!> Runs CorAdv_robust_enstro_fortran, a capture of its inputs and outputs, or the C++/AMReX bridge,
!! as selected by the CORADV_ROBUST_ENSTRO_MODE environment variable.
subroutine CorAdv_robust_enstro(bxU, bxV, u_a, v_a, h_a, uh_a, vh_a, abs_vort_a, eps_vel, h_tiny, &
                                Gcore, CAu_a, CAv_a, PV_Adv_Scheme)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a        !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a        !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: h_a        !< Layer thickness [H ~> m or kg m-2]
  type(RealArray_t), intent(in)    :: uh_a       !< Zonal transport u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_a       !< Meridional transport v*h*dx
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: abs_vort_a !< Absolute vorticity at q-points [T-1 ~> s-1]
  real,                    intent(in)    :: eps_vel !< A tiny, positive velocity [L T-1 ~> m s-1]
  real,                    intent(in)    :: h_tiny  !< A very small thickness [H ~> m or kg m-2]
  type(RealArray_t), intent(inout) :: CAu_a      !< Zonal acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAv_a      !< Meridional acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  integer,                 intent(in)    :: PV_Adv_Scheme !< Selects the discretization for PV
                                          !! advection: PV_ADV_CENTERED or PV_ADV_UPWIND1
  ! Local variables
  integer            :: mode, rc
  type(Box_C)        :: bxU_c, bxV_c
  type(RealArray_C)  :: u_c, v_c, h_c, uh_c, vh_c, abs_vort_c, CAu_c, CAv_c
  integer(c_int) :: PV_Adv_Scheme_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "coradv_robust_enstro"
  mode = getenv_mode("CORADV_ROBUST_ENSTRO_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxU", bxU)
        call rec%add("_bxV", bxV)
        call rec%add("_u", u_a)
        call rec%add("_v", v_a)
        call rec%add("_h", h_a)
        call rec%add("_uh", uh_a)
        call rec%add("_vh", vh_a)
        call rec%add("_abs_vort", abs_vort_a)
        call rec%add("_eps_vel", eps_vel)
        call rec%add("_h_tiny", h_tiny)
        call rec%add("_IdxCv", Gcore%IdxCv)
        call rec%add("_IdxCu", Gcore%IdxCu)
        call rec%add("_IdyCu", Gcore%IdyCu)
        call rec%add("_IdyCv", Gcore%IdyCv)
        call rec%add("_CAu_before", CAu_a)
        call rec%add("_CAv_before", CAv_a)
        call rec%add("_PV_Adv_Scheme", PV_Adv_Scheme)
      endif

      call CorAdv_robust_enstro_fortran(bxU, bxV, u_a, v_a, h_a, uh_a, vh_a, abs_vort_a, eps_vel, &
                                        h_tiny, Gcore, CAu_a, CAv_a, PV_Adv_Scheme)

      if (capture) then
        call rec%add("_CAu_after", CAu_a)
        call rec%add("_CAv_after", CAv_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxU_c = bxU%to_c()
      bxV_c = bxV%to_c()
      u_c = u_a%to_c()
      v_c = v_a%to_c()
      h_c = h_a%to_c()
      uh_c = uh_a%to_c()
      vh_c = vh_a%to_c()
      abs_vort_c = abs_vort_a%to_c()
      CAu_c = CAu_a%to_c()
      CAv_c = CAv_a%to_c()
      PV_Adv_Scheme_c = PV_Adv_Scheme
      call turbotmp_coradv_robust_enstro_bridge(bxU_c, bxV_c, u_c, v_c, h_c, uh_c, vh_c, &
                                                abs_vort_c, eps_vel, h_tiny, Gcore%C, CAu_c, &
                                                CAv_c, PV_Adv_Scheme_c)
#endif

    case default
      call CorAdv_robust_enstro_fortran(bxU, bxV, u_a, v_a, h_a, uh_a, vh_a, abs_vort_a, eps_vel, &
                                        h_tiny, Gcore, CAu_a, CAv_a, PV_Adv_Scheme)

  end select

end subroutine CorAdv_robust_enstro


!> Calculates the Coriolis and momentum advection accelerations with the WENO PV-reconstruction
!! schemes (WENOVI7TH/5TH/3RD_PV_ENSTRO), over one tile.  Near land the reconstruction falls back to
!! narrower stencils, down to first-order upwind.
subroutine CorAdv_weno_fortran(bxU, bxV, u_a, v_a, uh_a, vh_a, q_a, abs_vort_a, h_q_a, Gcore, &
                               H_subroundoff, CAu_a, CAv_a, weno_velocity_smooth, CS)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  real, intent(in) :: H_subroundoff !< A thickness that is so small that it can be added to a
                                    !! thickness of Angstrom or larger without changing it at the
                                    !! bit level [H ~> m or kg m-2].
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a        !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a        !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: uh_a       !< Zonal transport u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_a       !< Meridional transport v*h*dx
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: q_a        !< Layer potential vorticity
                                                 !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(in)    :: abs_vort_a !< Absolute vorticity at q-points [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: h_q_a      !< Thickness interpolated to q points
                                                 !! [H ~> m or kg m-2]
  type(RealArray_t), intent(inout) :: CAu_a      !< Zonal acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAv_a      !< Meridional acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  logical,                 intent(in)    :: weno_velocity_smooth !< If true, use velocity to
                                          !! compute the smoothness indicator for WENO
  type(Coriolis_scheme_CS), intent(in)  :: CS  !< Options selecting the Coriolis discretization

  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: u, v, uh, vh, q, abs_vort, h_q, CAu, CAv
  real, dimension(:,:), contiguous, pointer :: IdxCu, mask2dCu, IdyCv, mask2dCv
  real :: u_v, v_u      ! u_v is the u velocity at v point, v_u is the v velocity at u point [L T-1 ~> m s-1]
  real :: q_v, q_u      ! PV at the u and v points [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  integer :: seventh_order, fifth_order, third_order ! Order of accuracy for the WENO calculations
  real :: u_q8(8) ! Eight-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
  real :: u_q6(6) ! Six-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
  real :: u_q4(4) ! Four-point zonal velocity at WENO stencils [L T-1 ~> m s-1]
  real :: v_q8(8) ! Eight-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
  real :: v_q6(6) ! Six-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
  real :: v_q4(4) ! Four-point meridional velocity at WENO stencils [L T-1 ~> m s-1]
  integer :: i, j, k

  call u_a%view(u) ; call v_a%view(v) ; call uh_a%view(uh)
  call vh_a%view(vh) ; call q_a%view(q) ; call abs_vort_a%view(abs_vort)
  call h_q_a%view(h_q) ; call CAu_a%view(CAu) ; call CAv_a%view(CAv)
  call Gcore%IdxCu%view(IdxCu) ; call Gcore%mask2dCu%view(mask2dCu) ; call Gcore%IdyCv%view(IdyCv)
  call Gcore%mask2dCv%view(mask2dCv)

  ! Calculate the tendencies of zonal velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAu =  q * vh - d(KE)/dx.
  if (CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO) then
    do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                   I=bxU%idxS(1):bxU%idxE(1)) &
        DO_LOCALITY(local(v_u, q_u, third_order, fifth_order, seventh_order, u_q8, u_q6, u_q4))
      v_u = 0.25*IdxCu(I,j)*((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
      ! check whether there is masked land points in the stencil
      third_order = (mask2dCu(I,j-2) * mask2dCu(I,j-1) * mask2dCu(I,j) * &
                     mask2dCu(I,j+1) * mask2dCu(I,j+2))

      fifth_order   = third_order * mask2dCu(I,j-3) * mask2dCu(I,j+3)
      seventh_order = fifth_order * mask2dCu(I,j-4) * mask2dCu(I,j+4)

      ! compute the masking to make sure that inland values are not used
      if (seventh_order == 1) then
        ! all values are valid, we use seventh order reconstruction
        u_q8(:) = (u(I,j-4:j+3,k) + u(I,j-3:j+4,k)) * 0.5
        call weno_seven_h_weight_reconstruction(abs_vort(I,J-4:J+3,k), &
                                       h_q(I,J-4:J+3,k), &
                                       u_q8, &
                                       H_subroundoff, v_u, q_u, weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)

      elseif (fifth_order == 1) then
        ! all values are valid, we use fifth order reconstruction
        u_q6(:) = (u(I,j-3:j+2,k) + u(I,j-2:j+3,k)) * 0.5
        call weno_five_h_weight_reconstruction(abs_vort(I,J-3:J+2,k), &
                                      h_q(I,J-3:J+2,k), &
                                      u_q6, &
                                      H_subroundoff, v_u, q_u, weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)

      elseif (third_order == 1) then
        ! only the middle values are valid, we use third order reconstruction
        u_q4(:) = (u(I,j-2:j+1,k) + u(I,j-1:j+2,k)) * 0.5
        call weno_three_h_weight_reconstruction(abs_vort(I,J-2:J+1,k), &
                                       h_q(I,J-2:J+1,k), &
                                       u_q4, &
                                       H_subroundoff, v_u, q_u, weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)
      else ! Upwind first order
        if (v_u>0.) then
            q_u = q(I,J-1,k)
        else
            q_u = q(I,J,k)
        endif
        CAu(I,j,k) = (q_u * v_u)

      endif
    enddo
  elseif (CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO) then
    do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                   I=bxU%idxS(1):bxU%idxE(1)) &
        DO_LOCALITY(local(v_u, q_u, third_order, fifth_order, u_q6, u_q4))
      v_u = 0.25*IdxCu(I,j)*((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
      third_order = (mask2dCu(I,j-2) * mask2dCu(I,j-1) * mask2dCu(I,j) * &
                     mask2dCu(I,j+1) * mask2dCu(I,j+2))

      fifth_order   = third_order * mask2dCu(I,j-3) * mask2dCu(I,j+3)

      if (fifth_order == 1) then
        ! all values are valid, we use fifth order reconstruction
        u_q6(:) = (u(I,j-3:j+2,k) + u(I,j-2:j+3,k)) * 0.5
        call weno_five_h_weight_reconstruction(abs_vort(I,J-3:J+2,k), &
                                      h_q(I,J-3:J+2,k), &
                                      u_q6, &
                                      H_subroundoff, v_u, q_u, weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)

      elseif (third_order == 1) then
        ! only the middle values are valid, we use third order reconstruction
        u_q4(:) = (u(I,j-2:j+1,k) + u(I,j-1:j+2,k)) * 0.5
        call weno_three_h_weight_reconstruction(abs_vort(I,J-2:J+1,k), &
                                       h_q(I,J-2:J+1,k), &
                                       u_q4, &
                                       H_subroundoff, v_u, q_u, weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)

      else ! Upwind first order
        if (v_u>0.) then
            q_u = q(I,J-1,k)
        else
            q_u = q(I,J,k)
        endif
        CAu(I,j,k) = (q_u * v_u)
      endif
    enddo
  elseif (CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO) then
    do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                   I=bxU%idxS(1):bxU%idxE(1)) &
        DO_LOCALITY(local(v_u, q_u, third_order, u_q4))
      v_u = 0.25*IdxCu(I,j)*((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
      third_order = (mask2dCu(I,j-2) * mask2dCu(I,j-1) * mask2dCu(I,j) * &
                     mask2dCu(I,j+1) * mask2dCu(I,j+2))

      if (third_order == 1) then
        ! only the middle values are valid, we use third order reconstruction
        u_q4(:) = (u(I,j-2:j+1,k) + u(I,j-1:j+2,k)) * 0.5
        call weno_three_h_weight_reconstruction(abs_vort(I,J-2:J+1,k), &
                                       h_q(I,J-2:J+1,k), &
                                       u_q4, &
                                       H_subroundoff, v_u, q_u, weno_velocity_smooth)
        CAu(I,j,k) = (q_u * v_u)

      else ! Upwind first order
        if (v_u>0.) then
            q_u = q(I,J-1,k)
        else
            q_u = q(I,J,k)
        endif
        CAu(I,j,k) = (q_u * v_u)
      endif
    enddo
  endif

  ! Calculate the tendencies of meridional velocity due to the Coriolis
  ! force and momentum advection.  On a Cartesian grid, this is
  !     CAv = - q * uh - d(KE)/dy.
  if (CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO) then
    do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                   i=bxV%idxS(1):bxV%idxE(1)) &
        DO_LOCALITY(local(u_v, q_v, third_order, fifth_order, seventh_order, v_q8, v_q6, v_q4))
      u_v = 0.25*IdyCv(i,J)*((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))

      ! check whether there is any masked land values within the stencils
      third_order = (mask2dCv(i-2,J) * mask2dCv(i-1,J) * mask2dCv(i,J) * mask2dCv(i+1,J) * &
                     mask2dCv(i+2,J))
      fifth_order   = third_order * mask2dCv(i-3,J) * mask2dCv(i+3,J)
      seventh_order = fifth_order * mask2dCv(i-4,J) * mask2dCv(i+4,J)

      ! compute the masking to make sure that inland values are not used
      if (seventh_order == 1) then
        v_q8(:) = (v(i-4:i+3,J,k) + v(i-3:i+4,J,k)) * 0.5
        ! all values are valid, we use seventh order reconstruction
        call weno_seven_h_weight_reconstruction(abs_vort(I-4:I+3,J,k), &
                                       h_q(I-4:I+3,J,k), &
                                       v_q8, &
                                       H_subroundoff, u_v, q_v, weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)

      elseif (fifth_order == 1) then
        v_q6(:) = (v(i-3:i+2,J,k) + v(i-2:i+3,J,k)) * 0.5
        ! all values are valid, we use fifth order reconstruction
        call weno_five_h_weight_reconstruction(abs_vort(I-3:I+2,J,k), &
                                      h_q(I-3:I+2,J,k), &
                                      v_q6, &
                                      H_subroundoff, u_v, q_v, weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)

      elseif (third_order == 1) then
        v_q4(:) = (v(i-2:i+1,J,k) + v(i-1:i+2,J,k)) * 0.5
!          ! only the middle values are valid, we use third order reconstruction
        call weno_three_h_weight_reconstruction(abs_vort(I-2:I+1,J,k), &
                                               h_q(I-2:I+1,J,k), &
                                               v_q4, &
                                               H_subroundoff, u_v, q_v, weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)
      else ! Upwind first order!
        if (u_v>0.) then
            q_v = q(I-1,J,k)
        else
            q_v = q(I,J,k)
        endif
        CAv(i,J,k) = - (q_v * u_v)
      endif

    enddo
  elseif (CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO) then
    do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                   i=bxV%idxS(1):bxV%idxE(1)) &
        DO_LOCALITY(local(u_v, q_v, third_order, fifth_order, v_q6, v_q4))
      u_v = 0.25*IdyCv(i,J)*((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))

      third_order = (mask2dCv(i-2,J) * mask2dCv(i-1,J) * mask2dCv(i,J) * mask2dCv(i+1,J) * &
                     mask2dCv(i+2,J))
      fifth_order   = third_order * mask2dCv(i-3,J) * mask2dCv(i+3,J)

      ! compute the masking to make sure that inland values are not used
      if (fifth_order == 1) then
        v_q6(:) = (v(i-3:i+2,J,k) + v(i-2:i+3,J,k)) * 0.5
        ! all values are valid, we use fifth order reconstruction
        call weno_five_h_weight_reconstruction(abs_vort(I-3:I+2,J,k), &
                                      h_q(I-3:I+2,J,k), &
                                      v_q6, &
                                      H_subroundoff, u_v, q_v, weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)

      elseif (third_order == 1) then
        v_q4(:) = (v(i-2:i+1,J,k) + v(i-1:i+2,J,k)) * 0.5
!          ! only the middle values are valid, we use third order reconstruction
        call weno_three_h_weight_reconstruction(abs_vort(I-2:I+1,J,k), &
                                               h_q(I-2:I+1,J,k), &
                                               v_q4, &
                                               H_subroundoff, u_v, q_v, weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)

      else
        if (u_v>0.) then
            q_v = q(I-1,J,k)
        else
            q_v = q(I,J,k)
        endif
        CAv(i,J,k) = - (q_v * u_v)
      endif

    enddo
  elseif (CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO) then
    do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                   i=bxV%idxS(1):bxV%idxE(1)) &
        DO_LOCALITY(local(u_v, q_v, third_order, v_q4))
      u_v = 0.25*IdyCv(i,J)*((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))

      third_order = (mask2dCv(i-2,J) * mask2dCv(i-1,J) * mask2dCv(i,J) * mask2dCv(i+1,J) * &
                     mask2dCv(i+2,J))

      ! compute the masking to make sure that inland values are not used
      if (third_order == 1) then
        v_q4(:) = (v(i-2:i+1,J,k) + v(i-1:i+2,J,k)) * 0.5
!          ! only the middle values are valid, we use third order reconstruction
        call weno_three_h_weight_reconstruction(abs_vort(I-2:I+1,J,k), &
                                               h_q(I-2:I+1,J,k), &
                                               v_q4, &
                                               H_subroundoff, u_v, q_v, weno_velocity_smooth)
        CAv(i,J,k) = - (q_v * u_v)

      else
        if (u_v>0.) then
            q_v = q(I-1,J,k)
        else
            q_v = q(I,J,k)
        endif
        CAv(i,J,k) = - (q_v * u_v)
      endif

    enddo
  endif
end subroutine CorAdv_weno_fortran

!> Runs CorAdv_weno_fortran, a capture of its inputs and outputs, or the C++/AMReX bridge,
!! as selected by the CORADV_WENO_MODE environment variable.
subroutine CorAdv_weno(bxU, bxV, u_a, v_a, uh_a, vh_a, q_a, abs_vort_a, h_q_a, Gcore, &
                       H_subroundoff, CAu_a, CAv_a, weno_velocity_smooth, CS)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  real, intent(in) :: H_subroundoff !< A thickness that is so small that it can be added to a
                                    !! thickness of Angstrom or larger without changing it at the
                                    !! bit level [H ~> m or kg m-2].
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a        !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a        !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: uh_a       !< Zonal transport u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_a       !< Meridional transport v*h*dx
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: q_a        !< Layer potential vorticity
                                                 !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(in)    :: abs_vort_a !< Absolute vorticity at q-points [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: h_q_a      !< Thickness interpolated to q points
                                                 !! [H ~> m or kg m-2]
  type(RealArray_t), intent(inout) :: CAu_a      !< Zonal acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAv_a      !< Meridional acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  logical,                 intent(in)    :: weno_velocity_smooth !< If true, use velocity to
                                          !! compute the smoothness indicator for WENO
  type(Coriolis_scheme_CS), intent(in)  :: CS  !< Options selecting the Coriolis discretization
  ! Local variables
  integer            :: mode, rc
  type(Box_C)        :: bxU_c, bxV_c
  type(RealArray_C)  :: u_c, v_c, uh_c, vh_c, q_c, abs_vort_c, h_q_c, CAu_c, CAv_c
  logical(c_bool) :: weno_velocity_smooth_c
  type(Coriolis_scheme_CS_C) :: CS_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "coradv_weno"
  mode = getenv_mode("CORADV_WENO_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxU", bxU)
        call rec%add("_bxV", bxV)
        call rec%add("_u", u_a)
        call rec%add("_v", v_a)
        call rec%add("_uh", uh_a)
        call rec%add("_vh", vh_a)
        call rec%add("_q", q_a)
        call rec%add("_abs_vort", abs_vort_a)
        call rec%add("_h_q", h_q_a)
        call rec%add("_IdxCu", Gcore%IdxCu)
        call rec%add("_mask2dCu", Gcore%mask2dCu)
        call rec%add("_IdyCv", Gcore%IdyCv)
        call rec%add("_mask2dCv", Gcore%mask2dCv)
        call rec%add("_H_subroundoff", H_subroundoff)
        call rec%add("_CAu_before", CAu_a)
        call rec%add("_CAv_before", CAv_a)
        call rec%add("_weno_velocity_smooth", weno_velocity_smooth)
        call record_Coriolis_scheme_CS(rec, "_CS", CS)
      endif

      call CorAdv_weno_fortran(bxU, bxV, u_a, v_a, uh_a, vh_a, q_a, abs_vort_a, h_q_a, Gcore, &
                               H_subroundoff, CAu_a, CAv_a, weno_velocity_smooth, CS)

      if (capture) then
        call rec%add("_CAu_after", CAu_a)
        call rec%add("_CAv_after", CAv_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxU_c = bxU%to_c()
      bxV_c = bxV%to_c()
      u_c = u_a%to_c()
      v_c = v_a%to_c()
      uh_c = uh_a%to_c()
      vh_c = vh_a%to_c()
      q_c = q_a%to_c()
      abs_vort_c = abs_vort_a%to_c()
      h_q_c = h_q_a%to_c()
      CAu_c = CAu_a%to_c()
      CAv_c = CAv_a%to_c()
      weno_velocity_smooth_c = weno_velocity_smooth
      CS_c = Coriolis_scheme_CS_to_c(CS)
      call turbotmp_coradv_weno_bridge(bxU_c, bxV_c, u_c, v_c, uh_c, vh_c, q_c, abs_vort_c, h_q_c, &
                                       Gcore%C, H_subroundoff, CAu_c, CAv_c, &
                                       weno_velocity_smooth_c, CS_c)
#endif

    case default
      call CorAdv_weno_fortran(bxU, bxV, u_a, v_a, uh_a, vh_a, q_a, abs_vort_a, h_q_a, Gcore, &
                               H_subroundoff, CAu_a, CAv_a, weno_velocity_smooth, CS)

  end select

end subroutine CorAdv_weno


!> Adds the terms that are common to every Coriolis scheme to the accelerations over one tile:
!! the Stokes-drift diagnostic, the optional bounding of the Coriolis terms and the kinetic energy
!! gradient, and the diagnostics of the kinetic energy gradient and relative-vorticity terms.
subroutine CorAdv_common_terms_fortran(bxU, bxV, u_a, v_a, uh_a, vh_a, abs_vort_a, qS_a, q2_a, &
                                       KEx_a, KEy_a, Stokes_VF, Gcore, CAu_a, CAv_a, CAuS_a, &
                                       CAvS_a, AD_a, CS)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a        !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a        !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: uh_a       !< Zonal transport u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_a       !< Meridional transport v*h*dx
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: abs_vort_a !< Absolute vorticity at q-points [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: qS_a       !< Layer Stokes vorticity
                                                 !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(in)    :: q2_a       !< Relative vorticity over thickness
                                                 !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(in)    :: KEx_a      !< Zonal gradient of kinetic energy per unit mass
                                                 !! [L T-2 ~> m s-2]
  type(RealArray_t), intent(in)    :: KEy_a      !< Meridional gradient of kinetic energy per unit
                                                 !! mass
                                                 !! [L T-2 ~> m s-2]
  logical,                 intent(in)    :: Stokes_VF !< If true, include the Stokes drift
  type(RealArray_t), intent(inout) :: CAu_a      !< Zonal acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAv_a      !< Meridional acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAuS_a     !< Stokes contribution to CAu [L T-2 ~> m s-2]
  type(RealArray_t), intent(inout) :: CAvS_a     !< Stokes contribution to CAv [L T-2 ~> m s-2]
  type(CorAdv_AD_type),    intent(inout) :: AD_a !< Copies of the acceleration diagnostics
  type(CorAdv_opts_CS),    intent(in)    :: CS  !< Options used in the Coriolis and momentum
                                                !! advection calculations

  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: u, v, uh, vh, abs_vort, qS, q2, KEx, KEy, CAu, &
    CAv, CAuS, CAvS
  real, dimension(:,:), contiguous, pointer :: IdxCu, IdyCv
  real, dimension(:,:,:), contiguous, pointer :: gradKEu, gradKEv, rv_x_u, rv_x_v
  real :: fv1, fv2, fv3, fv4   ! (f+rv)*v at the 4 points surrounding a u points[L T-2 ~> m s-2]
  real :: fu1, fu2, fu3, fu4   ! -(f+rv)*u at the 4 points surrounding a v point [L T-2 ~> m s-2]
  real :: max_fv, max_fu       ! The maximum of the neighboring Coriolis accelerations [L T-2 ~> m s-2]
  real :: min_fv, min_fu       ! The minimum of the neighboring Coriolis accelerations [L T-2 ~> m s-2]
  real, parameter :: C1_12 = 1.0 / 12.0 ! C1_12 = 1/12 [nondim]
  integer :: i, j, k

  call u_a%view(u) ; call v_a%view(v) ; call uh_a%view(uh)
  call vh_a%view(vh) ; call abs_vort_a%view(abs_vort) ; call qS_a%view(qS)
  call q2_a%view(q2) ; call KEx_a%view(KEx) ; call KEy_a%view(KEy)
  call CAu_a%view(CAu) ; call CAv_a%view(CAv) ; call CAuS_a%view(CAuS)
  call CAvS_a%view(CAvS)
  call Gcore%IdxCu%view(IdxCu) ; call Gcore%IdyCv%view(IdyCv)
  if (AD_a%gradKEu_a%associated()) call AD_a%gradKEu_a%view(gradKEu)
  if (AD_a%gradKEv_a%associated()) call AD_a%gradKEv_a%view(gradKEv)
  if (AD_a%rv_x_u_a%associated()) call AD_a%rv_x_u_a%view(rv_x_u)
  if (AD_a%rv_x_v_a%associated()) call AD_a%rv_x_v_a%view(rv_x_v)

  if (Stokes_VF) then
    if (CS%do_CAS_diag) then
      ! Computing the diagnostic Stokes contribution to CAu
      do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                     I=bxU%idxS(1):bxU%idxE(1))
        CAuS(I,j,k) = 0.25 * &
              ((qS(I,J,k) * (vh(i+1,J,k) + vh(i,J,k))) + &
               (qS(I,J-1,k) * (vh(i,J-1,k) + vh(i+1,J-1,k)))) * IdxCu(I,j)
      enddo
    endif
  endif

  if (CS%bound_Coriolis) then
    do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                   I=bxU%idxS(1):bxU%idxE(1)) DO_LOCALITY(local(fv1, fv2, fv3, fv4, max_fv, min_fv))
      fv1 = abs_vort(I,J,k) * v(i+1,J,k)
      fv2 = abs_vort(I,J,k) * v(i,J,k)
      fv3 = abs_vort(I,J-1,k) * v(i+1,J-1,k)
      fv4 = abs_vort(I,J-1,k) * v(i,J-1,k)

      max_fv = max(fv1, fv2, fv3, fv4)
      min_fv = min(fv1, fv2, fv3, fv4)

      CAu(I,j,k) = min(CAu(I,j,k), max_fv)
      CAu(I,j,k) = max(CAu(I,j,k), min_fv)
    enddo
  endif

  ! Term - d(KE)/dx.
  do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                 I=bxU%idxS(1):bxU%idxE(1))
    CAu(I,j,k) = CAu(I,j,k) - KEx(I,j,k)
  enddo

  if (AD_a%gradKEu_a%associated()) then
    do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                   I=bxU%idxS(1):bxU%idxE(1))
      gradKEu(I,j,k) = -KEx(I,j,k)
    enddo
  endif

  if (Stokes_VF) then
    if (CS%do_CAS_diag) then
      ! Computing the diagnostic Stokes contribution to CAv
      do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                     i=bxV%idxS(1):bxV%idxE(1))
        CAvS(i,J,k) = 0.25 * &
              ((qS(I,J,k) * (uh(I,j+1,k) + uh(I,j,k))) + &
               (qS(I-1,J,k) * (uh(I-1,j,k) + uh(I-1,j+1,k)))) * IdyCv(i,J)
      enddo
    endif
  endif

  if (CS%bound_Coriolis) then
    do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                   i=bxV%idxS(1):bxV%idxE(1)) DO_LOCALITY(local(fu1, fu2, fu3, fu4, max_fu, min_fu))
      fu1 = -abs_vort(I,J,k) * u(I,j+1,k)
      fu2 = -abs_vort(I,J,k) * u(I,j,k)
      fu3 = -abs_vort(I-1,J,k) * u(I-1,j+1,k)
      fu4 = -abs_vort(I-1,J,k) * u(I-1,j,k)

      max_fu = max(fu1, fu2, fu3, fu4)
      min_fu = min(fu1, fu2, fu3, fu4)

      CAv(I,j,k) = min(CAv(I,j,k), max_fu)
      CAv(I,j,k) = max(CAv(I,j,k), min_fu)
    enddo
  endif

  ! Term - d(KE)/dy.
  do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                 i=bxV%idxS(1):bxV%idxE(1))
    CAv(i,J,k) = CAv(i,J,k) - KEy(i,J,k)
  enddo
  if (AD_a%gradKEv_a%associated()) then
    do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                   i=bxV%idxS(1):bxV%idxE(1))
      gradKEv(i,J,k) = -KEy(i,J,k)
    enddo
  endif

  if (AD_a%rv_x_u_a%associated() .or. AD_a%rv_x_v_a%associated()) then
    ! Calculate the Coriolis-like acceleration due to relative vorticity.
    if (CS%Coriolis_scheme_CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
      if (AD_a%rv_x_u_a%associated()) then
        do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                       i=bxV%idxS(1):bxV%idxE(1))
          rv_x_u(i,J,k) = - 0.25* &
            ((q2(I-1,j,k)*(uh(I-1,j,k) + uh(I-1,j+1,k))) + &
             (q2(I,j,k)*(uh(I,j,k) + uh(I,j+1,k)))) * IdyCv(i,J)
        enddo
      endif

      if (AD_a%rv_x_v_a%associated()) then
        do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                       I=bxU%idxS(1):bxU%idxE(1))
          rv_x_v(I,j,k) = 0.25 * &
            ((q2(I,j,k) * (vh(i+1,J,k) + vh(i,J,k))) + &
             (q2(I,j-1,k) * (vh(i,J-1,k) + vh(i+1,J-1,k)))) * IdxCu(I,j)
        enddo
      endif
    else
      if (AD_a%rv_x_u_a%associated()) then
        do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                       i=bxV%idxS(1):bxV%idxE(1))
          rv_x_u(i,J,k) = -IdyCv(i,J) * C1_12 * &
            (((((q2(I,J,k) + q2(I-1,J-1,k)) + q2(I-1,J,k)) * uh(I-1,j,k)) + &
              (((q2(I-1,J,k) + q2(I,J+1,k)) + q2(I,J,k)) * uh(I,j+1,k))) + &
             ((((q2(I-1,J,k) + q2(I,J-1,k)) + q2(I,J,k)) * uh(I,j,k))+ &
              (((q2(I,J,k) + q2(I-1,J+1,k)) + q2(I-1,J,k)) * uh(I-1,j+1,k))))
        enddo
      endif

      if (AD_a%rv_x_v_a%associated()) then
        do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                       I=bxU%idxS(1):bxU%idxE(1))
          rv_x_v(I,j,k) = IdxCu(I,j) * C1_12 * &
            (((((q2(I+1,J,k) + q2(I,J-1,k)) + q2(I,J,k)) * vh(i+1,J,k)) + &
              (((q2(I-1,J-1,k) + q2(I,J,k)) + q2(I,J-1,k)) * vh(i,J-1,k))) + &
             ((((q2(I-1,J,k) + q2(I,J-1,k)) + q2(I,J,k)) * vh(i,J,k)) + &
              (((q2(I+1,J-1,k) + q2(I,J,k)) + q2(I,J-1,k)) * vh(i+1,J-1,k))))
        enddo
      endif
    endif
  endif
end subroutine CorAdv_common_terms_fortran

!> Runs CorAdv_common_terms_fortran, a capture of its inputs and outputs, or the C++/AMReX bridge,
!! as selected by the CORADV_COMMON_TERMS_MODE environment variable.
subroutine CorAdv_common_terms(bxU, bxV, u_a, v_a, uh_a, vh_a, abs_vort_a, qS_a, q2_a, KEx_a, &
                               KEy_a, Stokes_VF, Gcore, CAu_a, CAv_a, CAuS_a, CAvS_a, AD_a, CS)
  type(grid_core_type),    intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),             intent(in)    :: bxU !< The u-point iteration box of this tile,
                                                !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),             intent(in)    :: bxV !< The v-point iteration box of this tile,
                                                !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a        !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a        !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: uh_a       !< Zonal transport u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: vh_a       !< Meridional transport v*h*dx
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t), intent(in)    :: abs_vort_a !< Absolute vorticity at q-points [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: qS_a       !< Layer Stokes vorticity
                                                 !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(in)    :: q2_a       !< Relative vorticity over thickness
                                                 !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(in)    :: KEx_a      !< Zonal gradient of kinetic energy per unit mass
                                                 !! [L T-2 ~> m s-2]
  type(RealArray_t), intent(in)    :: KEy_a      !< Meridional gradient of kinetic energy per unit
                                                 !! mass
                                                 !! [L T-2 ~> m s-2]
  logical,                 intent(in)    :: Stokes_VF !< If true, include the Stokes drift
  type(RealArray_t), intent(inout) :: CAu_a      !< Zonal acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAv_a      !< Meridional acceleration due to Coriolis
                                                 !! and momentum advection [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: CAuS_a     !< Stokes contribution to CAu [L T-2 ~> m s-2]
  type(RealArray_t), intent(inout) :: CAvS_a     !< Stokes contribution to CAv [L T-2 ~> m s-2]
  type(CorAdv_AD_type),    intent(inout) :: AD_a !< Copies of the acceleration diagnostics
  type(CorAdv_opts_CS),    intent(in)    :: CS  !< Options used in the Coriolis and momentum
                                                !! advection calculations
  ! Local variables
  integer            :: mode, rc
  type(Box_C)        :: bxU_c, bxV_c
  type(RealArray_C)  :: u_c, v_c, uh_c, vh_c, abs_vort_c, qS_c, q2_c, KEx_c, KEy_c, CAu_c, CAv_c, &
    CAuS_c, CAvS_c
  logical(c_bool) :: Stokes_VF_c
  type(CorAdv_AD_C) :: AD_a_c
  type(CorAdv_opts_CS_C) :: CS_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "coradv_common_terms"
  mode = getenv_mode("CORADV_COMMON_TERMS_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxU", bxU)
        call rec%add("_bxV", bxV)
        call rec%add("_u", u_a)
        call rec%add("_v", v_a)
        call rec%add("_uh", uh_a)
        call rec%add("_vh", vh_a)
        call rec%add("_abs_vort", abs_vort_a)
        call rec%add("_qS", qS_a)
        call rec%add("_q2", q2_a)
        call rec%add("_KEx", KEx_a)
        call rec%add("_KEy", KEy_a)
        call rec%add("_Stokes_VF", Stokes_VF)
        call rec%add("_IdxCu", Gcore%IdxCu)
        call rec%add("_IdyCv", Gcore%IdyCv)
        call rec%add("_CAu_before", CAu_a)
        call rec%add("_CAv_before", CAv_a)
        call rec%add("_CAuS_before", CAuS_a)
        call rec%add("_CAvS_before", CAvS_a)
        call record_CorAdv_AD(rec, "_AD", AD_a, "")
        call record_CorAdv_AD(rec, "_AD", AD_a, "_before")
        call record_CorAdv_opts_CS(rec, "_CS", CS)
      endif

      call CorAdv_common_terms_fortran(bxU, bxV, u_a, v_a, uh_a, vh_a, abs_vort_a, qS_a, q2_a, &
                                       KEx_a, KEy_a, Stokes_VF, Gcore, CAu_a, CAv_a, CAuS_a, &
                                       CAvS_a, AD_a, CS)

      if (capture) then
        call rec%add("_CAu_after", CAu_a)
        call rec%add("_CAv_after", CAv_a)
        call rec%add("_CAuS_after", CAuS_a)
        call rec%add("_CAvS_after", CAvS_a)
        call record_CorAdv_AD(rec, "_AD", AD_a, "_after")
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxU_c = bxU%to_c()
      bxV_c = bxV%to_c()
      u_c = u_a%to_c()
      v_c = v_a%to_c()
      uh_c = uh_a%to_c()
      vh_c = vh_a%to_c()
      abs_vort_c = abs_vort_a%to_c()
      qS_c = qS_a%to_c()
      q2_c = q2_a%to_c()
      KEx_c = KEx_a%to_c()
      KEy_c = KEy_a%to_c()
      Stokes_VF_c = Stokes_VF
      CAu_c = CAu_a%to_c()
      CAv_c = CAv_a%to_c()
      CAuS_c = CAuS_a%to_c()
      CAvS_c = CAvS_a%to_c()
      AD_a_c = CorAdv_AD_to_c(AD_a)
      CS_c = CorAdv_opts_CS_to_c(CS)
      call turbotmp_coradv_common_terms_bridge(bxU_c, bxV_c, u_c, v_c, uh_c, vh_c, abs_vort_c, &
                                               qS_c, q2_c, KEx_c, KEy_c, Stokes_VF_c, Gcore%C, &
                                               CAu_c, CAv_c, CAuS_c, CAvS_c, AD_a_c, CS_c)
#endif

    case default
      call CorAdv_common_terms_fortran(bxU, bxV, u_a, v_a, uh_a, vh_a, abs_vort_a, qS_a, q2_a, &
                                       KEx_a, KEy_a, Stokes_VF, Gcore, CAu_a, CAv_a, CAuS_a, &
                                       CAvS_a, AD_a, CS)

  end select

end subroutine CorAdv_common_terms


!> Offers the Coriolis-related derived quantities for averaging, once the accelerations over
!! every tile have been calculated.
subroutine CorAdv_finalize_diagnostics(RV_a, PV_a, CAuS_a, CAvS_a, Stokes_VF, AD, G, ke, CS)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  integer, intent(in) :: ke !< The number of layers/levels in the vertical
  type(RealArray_t), intent(in)    :: RV_a   !< Diagnostic relative vorticity [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: PV_a   !< Diagnostic potential vorticity
                                             !! [H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]
  type(RealArray_t), intent(in)    :: CAuS_a !< Stokes contribution to CAu [L T-2 ~> m s-2]
  type(RealArray_t), intent(in)    :: CAvS_a !< Stokes contribution to CAv [L T-2 ~> m s-2]
  logical,                 intent(in)    :: Stokes_VF !< If true, include the Stokes drift
  type(accel_diag_ptrs),   intent(in)    :: AD  !< Storage for acceleration diagnostics
  type(CoriolisAdv_CS),    intent(in)    :: CS  !< Control structure for MOM_CoriolisAdv

  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: RV, PV, CAuS, CAvS
  integer :: nz

  call RV_a%view(RV) ; call PV_a%view(PV) ; call CAuS_a%view(CAuS)
  call CAvS_a%view(CAvS)

  nz = ke

  ! Here the various Coriolis-related derived quantities are offered for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_rv > 0) call post_data(CS%id_rv, RV, CS%diag)
    if (CS%id_PV > 0) call post_data(CS%id_PV, PV, CS%diag)
    if (CS%id_gKEu>0) call post_data(CS%id_gKEu, AD%gradKEu, CS%diag)
    if (CS%id_gKEv>0) call post_data(CS%id_gKEv, AD%gradKEv, CS%diag)
    if (CS%id_rvxu > 0) call post_data(CS%id_rvxu, AD%rv_x_u, CS%diag)
    if (CS%id_rvxv > 0) call post_data(CS%id_rvxv, AD%rv_x_v, CS%diag)
    if (Stokes_VF) then
      if (CS%id_CAuS > 0) call post_data(CS%id_CAuS, CAuS, CS%diag)
      if (CS%id_CAvS > 0) call post_data(CS%id_CAvS, CAvS, CS%diag)
    endif

    ! Diagnostics for terms multiplied by fractional thicknesses

    ! 3D diagnostics hf_gKEu etc. are commented because there is no clarity on proper remapping grid option.
    ! The code is retained for debugging purposes in the future.
    ! if (CS%id_hf_gKEu > 0) call post_product_u(CS%id_hf_gKEu, AD%gradKEu, AD%diag_hfrac_u, G, nz, CS%diag)
    ! if (CS%id_hf_gKEv > 0) call post_product_v(CS%id_hf_gKEv, AD%gradKEv, AD%diag_hfrac_v, G, nz, CS%diag)
    ! if (CS%id_hf_rvxv > 0) call post_product_u(CS%id_hf_rvxv, AD%rv_x_v, AD%diag_hfrac_u, G, nz, CS%diag)
    ! if (CS%id_hf_rvxu > 0) call post_product_v(CS%id_hf_rvxu, AD%rv_x_u, AD%diag_hfrac_v, G, nz, CS%diag)

    if (CS%id_hf_gKEu_2d > 0) call post_product_sum_u(CS%id_hf_gKEu_2d, AD%gradKEu, AD%diag_hfrac_u, G, nz, CS%diag)
    if (CS%id_hf_gKEv_2d > 0) call post_product_sum_v(CS%id_hf_gKEv_2d, AD%gradKEv, AD%diag_hfrac_v, G, nz, CS%diag)
    if (CS%id_intz_gKEu_2d > 0) call post_product_sum_u(CS%id_intz_gKEu_2d, AD%gradKEu, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_intz_gKEv_2d > 0) call post_product_sum_v(CS%id_intz_gKEv_2d, AD%gradKEv, AD%diag_hv, G, nz, CS%diag)

    if (CS%id_hf_rvxv_2d > 0) call post_product_sum_u(CS%id_hf_rvxv_2d, AD%rv_x_v, AD%diag_hfrac_u, G, nz, CS%diag)
    if (CS%id_hf_rvxu_2d > 0) call post_product_sum_v(CS%id_hf_rvxu_2d, AD%rv_x_u, AD%diag_hfrac_v, G, nz, CS%diag)

    if (CS%id_h_gKEu > 0) call post_product_u(CS%id_h_gKEu, AD%gradKEu, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_h_gKEv > 0) call post_product_v(CS%id_h_gKEv, AD%gradKEv, AD%diag_hv, G, nz, CS%diag)
    if (CS%id_h_rvxv > 0) call post_product_u(CS%id_h_rvxv, AD%rv_x_v, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_h_rvxu > 0) call post_product_v(CS%id_h_rvxu, AD%rv_x_u, AD%diag_hv, G, nz, CS%diag)

    if (CS%id_intz_rvxv_2d > 0) call post_product_sum_u(CS%id_intz_rvxv_2d, AD%rv_x_v, AD%diag_hu, G, nz, CS%diag)
    if (CS%id_intz_rvxu_2d > 0) call post_product_sum_v(CS%id_intz_rvxu_2d, AD%rv_x_u, AD%diag_hv, G, nz, CS%diag)
  endif
end subroutine CorAdv_finalize_diagnostics


!> Calculates the acceleration due to the gradient of kinetic energy over one iteration tile.
subroutine gradKE_fortran(bxQ, bxU, bxV, u_a, v_a, KE_a, KEx_a, KEy_a, Gcore, Gobc, KE_Scheme, &
                          KE_use_limiter)
  type(grid_core_type),                       intent(in)  :: Gcore !< Persistent copies of the core
                                                                   !! ocean grid fields
  type(grid_OBC_type),                        intent(in)  :: Gobc !< Persistent copies of the ocean
                                                                  !! grid OBC mask fields
  type(Box_t),                                intent(in)  :: bxQ !< The B-grid-index iteration box of
                                                                 !! this tile, [IscB:IecB, JscB:JecB, ksc:kec]
  type(Box_t),                   intent(in)  :: bxU !< The u-point iteration box of this tile,
                                                    !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),                   intent(in)  :: bxV !< The v-point iteration box of this tile,
                                                    !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a   !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a   !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(inout) :: KE_a  !< Kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
  type(RealArray_t), intent(inout) :: KEx_a !< Zonal acceleration due to kinetic
                                            !! energy gradient [L T-2 ~> m s-2]
  type(RealArray_t), intent(inout) :: KEy_a !< Meridional acceleration due to kinetic
                                            !! energy gradient [L T-2 ~> m s-2]
  integer,                                    intent(in)  :: KE_Scheme !< Selects the discretization
                                                   !! for the kinetic energy: KE_ARAKAWA,
                                                   !! KE_SIMPLE_GUDONOV, KE_GUDONOV or KE_UP3
  logical,                                    intent(in)  :: KE_use_limiter !< If true, use the
                                                   !! Koren limiter for the KE_UP3 scheme
  ! Local variables
  type(Box_t) :: bxQ_pij  ! bxQ grown by one point at the end in i and j,
                          ! [IscB:IecB+1, JscB:JecB+1, ksc:kec]
  real, dimension(:,:,:), contiguous, pointer :: u, v, KE, KEx, KEy
  real, dimension(:,:), contiguous, pointer :: areaCu, areaCv, IareaT, mask2dCu, mask2dCv, &
    IdxCu_OBCmask, IdyCv_OBCmask
  real :: um, up, vm, vp         ! Temporary variables [L T-1 ~> m s-1].
  real :: um2, up2, vm2, vp2     ! Temporary variables [L2 T-2 ~> m2 s-2].
  real :: um2a, up2a, vm2a, vp2a ! Temporary variables [L4 T-2 ~> m4 s-2].
  real :: third_order_u, third_order_v  ! Product of mask values to determine the boundary
  integer :: i, j, k
  real, parameter     :: C1_12 = 1.0/12.0   ! The ratio of 1/12 [nondim]

  call u_a%view(u) ; call v_a%view(v)
  call KE_a%view(KE) ; call KEx_a%view(KEx) ; call KEy_a%view(KEy)
  call Gcore%areaCu%view(areaCu) ; call Gcore%areaCv%view(areaCv) ; call Gcore%IareaT%view(IareaT)
  call Gcore%mask2dCu%view(mask2dCu) ; call Gcore%mask2dCv%view(mask2dCv)
  call Gobc%IdxCu_OBCmask%view(IdxCu_OBCmask) ; call Gobc%IdyCv_OBCmask%view(IdyCv_OBCmask)

  bxQ_pij = bxQ%growBy(lo=[0,0,0], hi=[1,1,0])

  ! Calculate KE (Kinetic energy for use in the -grad(KE) acceleration term).
  if (KE_Scheme == KE_ARAKAWA) then
    ! The following calculation of Kinetic energy includes the metric terms
    ! identified in Arakawa & Lamb 1982 as important for KE conservation.  It
    ! also includes the possibility of partially-blocked tracer cell faces.
    do concurrent (k=bxQ_pij%idxS(3):bxQ_pij%idxE(3), j=bxQ_pij%idxS(2):bxQ_pij%idxE(2), &
                   i=bxQ_pij%idxS(1):bxQ_pij%idxE(1))
      KE(i,j,k) = ( ( (areaCu( I ,j)*(u( I ,j,k)*u( I ,j,k))) + &
                       (areaCu(I-1,j)*(u(I-1,j,k)*u(I-1,j,k))) ) + &
                     ( (areaCv(i, J )*(v(i, J ,k)*v(i, J ,k))) + &
                       (areaCv(i,J-1)*(v(i,J-1,k)*v(i,J-1,k))) ) )*0.25*IareaT(i,j)
    enddo
  elseif (KE_Scheme == KE_SIMPLE_GUDONOV) then
    ! The following discretization of KE is based on the one-dimensional Gudonov
    ! scheme which does not take into account any geometric factors
    do concurrent (k=bxQ_pij%idxS(3):bxQ_pij%idxE(3), j=bxQ_pij%idxS(2):bxQ_pij%idxE(2), &
                   i=bxQ_pij%idxS(1):bxQ_pij%idxE(1)) &
        DO_LOCALITY(local(up, um, vp, vm, up2, um2, vp2, vm2))
      up = 0.5*( u(I-1,j,k) + ABS( u(I-1,j,k) ) ) ; up2 = up*up
      um = 0.5*( u( I ,j,k) - ABS( u( I ,j,k) ) ) ; um2 = um*um
      vp = 0.5*( v(i,J-1,k) + ABS( v(i,J-1,k) ) ) ; vp2 = vp*vp
      vm = 0.5*( v(i, J ,k) - ABS( v(i, J ,k) ) ) ; vm2 = vm*vm
      KE(i,j,k) = ( max(up2,um2) + max(vp2,vm2) ) *0.5
    enddo
  elseif (KE_Scheme == KE_GUDONOV) then
    ! The following discretization of KE is based on the one-dimensional Gudonov
    ! scheme but has been adapted to take horizontal grid factors into account
    do concurrent (k=bxQ_pij%idxS(3):bxQ_pij%idxE(3), j=bxQ_pij%idxS(2):bxQ_pij%idxE(2), &
                   i=bxQ_pij%idxS(1):bxQ_pij%idxE(1)) &
        DO_LOCALITY(local(up, um, vp, vm, up2a, um2a, vp2a, vm2a))
      up = 0.5*( u(I-1,j,k) + ABS( u(I-1,j,k) ) ) ; up2a = up*up*areaCu(I-1,j)
      um = 0.5*( u( I ,j,k) - ABS( u( I ,j,k) ) ) ; um2a = um*um*areaCu( I ,j)
      vp = 0.5*( v(i,J-1,k) + ABS( v(i,J-1,k) ) ) ; vp2a = vp*vp*areaCv(i,J-1)
      vm = 0.5*( v(i, J ,k) - ABS( v(i, J ,k) ) ) ; vm2a = vm*vm*areaCv(i, J )
      KE(i,j,k) = ( max(um2a,up2a) + max(vm2a,vp2a) )*0.5*IareaT(i,j)
    enddo
  elseif (KE_Scheme == KE_UP3) then
    ! The following discretization of KE is based on the one-dimensional third-order
    ! upwind scheme which does not take horizontal grid factors into account

    if (KE_use_limiter) then
      do concurrent (k=bxQ_pij%idxS(3):bxQ_pij%idxE(3), j=bxQ_pij%idxS(2):bxQ_pij%idxE(2), &
                     i=bxQ_pij%idxS(1):bxQ_pij%idxE(1)) &
          DO_LOCALITY(local(up, um, vp, vm, third_order_u, third_order_v))
        ! compute the masking to make sure that inland values are not used
        third_order_u = (mask2dCu(I-2,j) * mask2dCu(I-1,j)* &
                       mask2dCu(I,j) * mask2dCu(I+1,j))

        if (third_order_u == 1) then
          up = (7.0 * (u(I-1,j,k) + u(I,j,k)) - (u(I-2,j,k) + u(I+1,j,k))) * C1_12
          call UP3_Koren_limiter_reconstruction(u(I-2:I+1,j,k), up, um)
        else
          up = (u(I-1,j,k) + u(I,j,k))*0.5
          if (up>0.) then
            um = u(I-1,j,k)
          elseif (up<0.) then
            um = u(I,j,k)
          else
            um = up
          endif
        endif

        third_order_v = (mask2dCv(i,J-2) * mask2dCv(i,J-1)* &
                       mask2dCv(i,J) * mask2dCv(i,J+1))
        if (third_order_v ==1) then
          vp = (7.0 * (v(i,J-1,k) + v(i,J,k)) - (v(i,J-2,k) + v(i,J+1,k))) * C1_12
          call UP3_Koren_limiter_reconstruction(v(i,J-2:J+1,k), vp, vm)
        else
          vp = (v(i,J-1,k) + v(i,J,k))*0.5
          if (vp>0.) then
            vm = v(i,J-1,k)
          elseif (vp<0.) then
            vm = v(i,J,k)
          else
            vm = vp
          endif
        endif

        KE(i,j,k) = ( (um*um) + (vm*vm) )*0.5
      enddo
    else
      do concurrent (k=bxQ_pij%idxS(3):bxQ_pij%idxE(3), j=bxQ_pij%idxS(2):bxQ_pij%idxE(2), &
                     i=bxQ_pij%idxS(1):bxQ_pij%idxE(1)) &
          DO_LOCALITY(local(up, um, vp, vm, third_order_u, third_order_v))
        ! compute the masking to make sure that inland values are not used
        third_order_u = (mask2dCu(I-2,j) * mask2dCu(I-1,j)* &
                       mask2dCu(I,j) * mask2dCu(I+1,j))

        if (third_order_u == 1) then
          up = (7.0 * (u(I-1,j,k) + u(I,j,k)) - (u(I-2,j,k) + u(I+1,j,k))) * C1_12
          call UP3_reconstruction(u(I-2:I+1,j,k), up, um)
        else
          up = (u(I-1,j,k) + u(I,j,k))*0.5
          if (up>0.) then
            um = u(I-1,j,k)
          elseif (up<0.) then
            um = u(I,j,k)
          else
            um = up
          endif
        endif

        third_order_v = (mask2dCv(i,J-2) * mask2dCv(i,J-1)* &
                       mask2dCv(i,J) * mask2dCv(i,J+1))
        if (third_order_v ==1) then
          vp = (7.0 * (v(i,J-1,k) + v(i,J,k)) - (v(i,J-2,k) + v(i,J+1,k))) * C1_12
          call UP3_reconstruction(v(i,J-2:J+1,k), vp, vm)
        else
          vp = (v(i,J-1,k) + v(i,J,k))*0.5
          if (vp>0.) then
            vm = v(i,J-1,k)
          elseif (vp<0.) then
            vm = v(i,J,k)
          else
            vm = vp
          endif
        endif

        KE(i,j,k) = ( (um*um) + (vm*vm) )*0.5
      enddo
    endif
  endif

  ! Term - d(KE)/dx.
  do concurrent (k=bxU%idxS(3):bxU%idxE(3), j=bxU%idxS(2):bxU%idxE(2), &
                 I=bxU%idxS(1):bxU%idxE(1))
    KEx(I,j,k) = (KE(i+1,j,k) - KE(i,j,k)) * IdxCu_OBCmask(I,j)
  enddo

  ! Term - d(KE)/dy.
  do concurrent (k=bxV%idxS(3):bxV%idxE(3), J=bxV%idxS(2):bxV%idxE(2), &
                 i=bxV%idxS(1):bxV%idxE(1))
    KEy(i,J,k) = (KE(i,j+1,k) - KE(i,j,k)) * IdyCv_OBCmask(i,J)
  enddo
  call bxQ_pij%free()
end subroutine gradKE_fortran

!> Calculates the acceleration due to the gradient of kinetic energy over one iteration tile,
!! with the Fortran implementation, a capture of its inputs and outputs, or the C++/AMReX
!! bridge, as selected by the GRADKE_MODE environment variable.
subroutine gradKE(bxQ, bxU, bxV, u_a, v_a, KE_a, KEx_a, KEy_a, Gcore, Gobc, KE_Scheme, &
                  KE_use_limiter)
  type(grid_core_type),                       intent(in)  :: Gcore !< Persistent copies of the core
                                                                   !! ocean grid fields
  type(grid_OBC_type),                        intent(in)  :: Gobc !< Persistent copies of the ocean
                                                                  !! grid OBC mask fields
  type(Box_t),                                intent(in)  :: bxQ !< The B-grid-index iteration box of
                                                                 !! this tile, [IscB:IecB, JscB:JecB, ksc:kec]
  type(Box_t),                   intent(in)  :: bxU !< The u-point iteration box of this tile,
                                                    !! [IscB:IecB, jsc:jec, ksc:kec]
  type(Box_t),                   intent(in)  :: bxV !< The v-point iteration box of this tile,
                                                    !! [isc:iec, JscB:JecB, ksc:kec]
  type(RealArray_t), intent(in)    :: u_a   !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_a   !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t), intent(inout) :: KE_a  !< Kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
  type(RealArray_t), intent(inout) :: KEx_a !< Zonal acceleration due to kinetic
                                            !! energy gradient [L T-2 ~> m s-2]
  type(RealArray_t), intent(inout) :: KEy_a !< Meridional acceleration due to kinetic
                                            !! energy gradient [L T-2 ~> m s-2]
  integer,                                    intent(in)  :: KE_Scheme !< Selects the discretization
                                                   !! for the kinetic energy: KE_ARAKAWA,
                                                   !! KE_SIMPLE_GUDONOV, KE_GUDONOV or KE_UP3
  logical,                                    intent(in)  :: KE_use_limiter !< If true, use the
                                                   !! Koren limiter for the KE_UP3 scheme
  ! Local variables
  integer            :: mode, rc
  type(Box_C)        :: bxQ_c, bxU_c, bxV_c
  type(RealArray_C)  :: u_c, v_c, KE_c, KEx_c, KEy_c
  integer(c_int)     :: KE_Scheme_c
  logical(c_bool)    :: KE_use_limiter_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "gradke"
  mode = getenv_mode("GRADKE_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxQ",            bxQ)
        call rec%add("_bxU",            bxU)
        call rec%add("_bxV",            bxV)
        call rec%add("_u",              u_a)
        call rec%add("_v",              v_a)
        call rec%add("_KE_before",      KE_a)
        call rec%add("_KEx_before",     KEx_a)
        call rec%add("_KEy_before",     KEy_a)
        call rec%add("_areaCu",         Gcore%areaCu)
        call rec%add("_areaCv",         Gcore%areaCv)
        call rec%add("_IareaT",         Gcore%IareaT)
        call rec%add("_mask2dCu",       Gcore%mask2dCu)
        call rec%add("_mask2dCv",       Gcore%mask2dCv)
        call rec%add("_IdxCu_OBCmask",  Gobc%IdxCu_OBCmask)
        call rec%add("_IdyCv_OBCmask",  Gobc%IdyCv_OBCmask)
        call rec%add("_KE_Scheme",      KE_Scheme)
        call rec%add("_KE_use_limiter", KE_use_limiter)
      endif

      call gradKE_fortran(bxQ, bxU, bxV, u_a, v_a, KE_a, KEx_a, KEy_a, Gcore, Gobc, KE_Scheme, &
                          KE_use_limiter)

      if (capture) then
        call rec%add("_KE_after",  KE_a)
        call rec%add("_KEx_after", KEx_a)
        call rec%add("_KEy_after", KEy_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxQ_c = bxQ%to_c() ; bxU_c = bxU%to_c() ; bxV_c = bxV%to_c()
      u_c   = u_a%to_c() ; v_c = v_a%to_c()
      KE_c  = KE_a%to_c() ; KEx_c = KEx_a%to_c() ; KEy_c = KEy_a%to_c()
      KE_Scheme_c      = KE_Scheme
      KE_use_limiter_c = KE_use_limiter
      call turbotmp_gradke_bridge(bxQ_c, bxU_c, bxV_c, u_c, v_c, KE_c, KEx_c, KEy_c, Gcore%C, &
                                  Gobc%C, KE_Scheme_c, KE_use_limiter_c)
#endif

    case default
      call gradKE_fortran(bxQ, bxU, bxV, u_a, v_a, KE_a, KEx_a, KEy_a, Gcore, Gobc, KE_Scheme, &
                          KE_use_limiter)

  end select

end subroutine gradKE

!> Reconstruct the scalar (e.g., pv, vorticity) onto point i-1/2 using a third-order upwind scheme
pure subroutine UP3_reconstruction(q4,u,qr)
  real, intent(in)    :: q4(4)            !< Tracer values on points i-2, i-1, i, i+1 [A ~> a]
  real, intent(in)    :: u                !< Velocity or thickness flux on point i-1/2
                                          !! [l t-1 ~> m s-1] or [l2 t-1 ~> m2 s-1]
  real, intent(inout) :: qr               !< Reconstruction of tracer q at point i-1/2 [A ~> a]
  real, parameter :: C1_6 = 1.0/6.0       ! The ratio of 1/6 [nondim]

  if (u>0.) then
    qr = ((2.*q4(3) + 5.*q4(2)) - q4(1)) * C1_6
  else
    qr = ((2.*q4(2) + 5.*q4(3)) - q4(4)) * C1_6
  endif

end subroutine UP3_reconstruction


!> Reconstruct the scalar (e.g., PV, vorticity) onto point i-1/2
!! using a third-order upwind scheme with the Koren flux limiter
pure subroutine UP3_Koren_limiter_reconstruction(q4,u,qr)
  real, intent(in)    :: q4(4)            !< Tracer values on points i-2, i-1, i, i+1 [A ~> a]
  real, intent(in)    :: u                !< Velocity or thickness flux on point i-1/2
                                          !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
  real, intent(inout) :: qr               !< Reconstruction of tracer q on point i-1/2 [A ~> a]
  real                :: theta            ! Ratio of gradient [nondim]
  real                :: psi              ! Limiter function [nondim]
  real, parameter     :: C1_3 = 1.0/3.0   ! The ratio of 1/3 [nondim]
  real, parameter     :: C1_6 = 1.0/6.0   ! The ratio of 1/6 [nondim]

  if (u>0.) then
    if (q4(3) == q4(2)) then
      qr = q4(2)
    else
      theta = (q4(2) - q4(1))/(q4(3) - q4(2))
      psi = max(0., min(1., C1_3 + C1_6*theta, theta)) ! limiter introduced by Koren (1993)
      qr = q4(2) + psi*(q4(3) - q4(2))
    endif
  else
    if (q4(3) == q4(2)) then
      qr = q4(3)
    else
      theta = (q4(4) - q4(3))/(q4(3) - q4(2))
      psi = max(0., min(1., C1_3 + C1_6*theta, theta))
      qr = q4(3) + psi*(q4(2) - q4(3))
    endif
  endif

end subroutine UP3_Koren_limiter_reconstruction

!> Compute the factor for the WENO weights
pure function fac_fn(tau, b) result(fac)
  real, intent(in)  :: tau  !< Difference of the smoothness indicator [A ~> a]
  real, intent(in)  :: b    !< The smoothness indicator [A ~> a]
  real :: fac               !< The factor for the weight [nondim]

  fac = 1.0e40 ; if (abs(b) > 1.0e-20*tau) fac = (1 + tau / b)**2

end function fac_fn


!> Reconstruct the tracer (e.g., PV, vorticity) onto the point i-1/2 using a third-order WENO scheme
!! This reconstruction is thickness-weighted
pure subroutine weno_three_h_weight_reconstruction(q4, h4, u4, &
                                              h_tiny, u, qr, velocity_smoothing)
    real, intent(in)    :: q4(4)   !< Tracer value times thickness on points i-2, i-1, i, i+1 [A ~> a]
    real, intent(in)    :: h4(4)   !< Thickness values on points i-2, i-1, i, i+1 [L ~> m]
    real, optional, intent(in)    :: u4(4) !< Velocity values on points i-2, i-1, i, i+1
                                                    !![L T-1 ~> m s-1]
    real, intent(in)    :: h_tiny  !< A tiny thickness to prevent division by zero [L ~> m]
    real, intent(in)    :: u              !< Velocity or thickness flux on point i-1/2
                                          !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
    real, intent(inout) :: qr             !< Reconstruction of tracer q on point i-1/2 [A ~> a]
    logical, intent(in) :: velocity_smoothing !< If true, use velocity to compute smoothness indicator
    real :: vr                            ! Reconstruction of hq [A ~> a]
    real :: hr                            ! Reconstruction of h [L ~> m]
    real :: c0, c1                        ! Intermediate reconstruction of q [A ~> a]
    real :: d0, d1                        ! Intermediate reconstruction of h [L ~> m]
    real :: b0, b1                        ! Smoothness indicator [A ~> a]
    real :: tau                           ! Difference of smoothness indicator [A ~> a]
    real :: w0, w1                        ! Weights [nondim]
    real :: s                             ! Temporary variables [nondim]
    real, parameter :: C2_3 = 2.0/3.0     ! The ratio of 2/3 [nondim]
    real, parameter :: C1_3 = 1.0/3.0     ! The ratio of 1/3 [nondim]

    if (u>0.) then
      call weno_three_reconstruction_0(q4(2:3), c0) ! Reconstruction in the second upwind stencil
      call weno_three_reconstruction_1(q4(1:2), c1) ! Reconstruction in the first upwind stencil

      call weno_three_reconstruction_0(h4(2:3), d0)
      call weno_three_reconstruction_1(h4(1:2), d1)
      if (velocity_smoothing) then
        call weno_three_weight(u4(2:3), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(u4(1:2), b1)  ! Smoothness indicator the first upwind stencil
      else
        call weno_three_weight(q4(2:3), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(q4(1:2), b1)  ! Smoothness indicator the first upwind stencil
      endif
    else
      call weno_three_reconstruction_0(q4(3:2:-1), c0) ! Reconstruction in the second upwind stencil
      call weno_three_reconstruction_1(q4(4:3:-1), c1) ! Reconstruction in the first upwind stencil

      call weno_three_reconstruction_0(h4(3:2:-1), d0)
      call weno_three_reconstruction_1(h4(4:3:-1), d1)
      if (velocity_smoothing) then
        call weno_three_weight(u4(3:2:-1), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(u4(4:3:-1), b1)  ! Smoothness indicator the first upwind stencil
      else
        call weno_three_weight(q4(3:2:-1), b0)  ! Smoothness indicator the second upwind stencil
        call weno_three_weight(q4(4:3:-1), b1)  ! Smoothness indicator the first upwind stencil
      endif
    endif

    tau = abs(b0-b1)
    w0  = C2_3 * fac_fn(tau, b0)
    w1  = C1_3 * fac_fn(tau, b1)

    s = 1. / (w0 + w1)
    w0 = w0 * s   ! Weights of stencils
    w1 = w1 * s

    vr = (w0 * c0) + (w1 * c1)
    hr = (w0 * d0) + (w1 * d1)
!    vr = min(max(q4(3), q4(2)), vr) ; vr = max(min(q4(3), q4(2)), vr) !Impose a monotonicity limiter
    hr = min(max(h4(3), h4(2)), hr) ; hr = max(min(h4(3), h4(2)), hr) ! A monotonicity limiter

    qr = vr / max(hr, h_tiny)

end subroutine weno_three_h_weight_reconstruction

!> Compute the smoothness indicator for the two-point stencil of the third-order WENO scheme
pure subroutine weno_three_weight(q2, w0)
    real, intent(in) :: q2(2)    !< Tracer values on the two-point stencil [A ~> a]
    real, intent(inout) :: w0    !< Smoothness indicator for this stencil [A2 ~> a2]

    w0 = (q2(1) - q2(2))**2

end subroutine weno_three_weight

!> Reconstruction in the second upwind stencil of the third-order WENO scheme
pure subroutine weno_three_reconstruction_0(q2, w0)
    real, intent(in) :: q2(2)    !< Tracer values on the two-point stencil [A ~> a]
    real, intent(inout) :: w0    !< Reconstruction of the quantity [A2 ~> a2]

    w0 = (q2(1) + q2(2)) * 0.5

end subroutine weno_three_reconstruction_0

!> Reconstruction in the first upwind stencil for third-order WENO scheme
pure subroutine weno_three_reconstruction_1(q2, w0)
    real, intent(in) :: q2(2)    !< Tracer values on the two-point stencil [A ~> a]
    real, intent(inout) :: w0    !< Reconstruction of the quantity [A ~> a]

    w0 = (- q2(1) + 3 * q2(2)) * 0.5

end subroutine weno_three_reconstruction_1


!> Reconstruct the tracer (e.g., PV, vorticity) onto point i-1/2 using a fifth-order WENO scheme
!! The reconstruction is weighted by the thickness
pure subroutine weno_five_h_weight_reconstruction(q6, h6, u6, &
                                             h_tiny, u, qr, velocity_smoothing)
    real, intent(in)    :: q6(6)
    !< Tracer values on points i-3, i-2, i-1, i, i+1, i+2 [A ~> a]
    real, intent(in)    :: h6(6)
    !< Thickness values on points i-3, i-2, i-1, i, i+1, i+2 [L ~> m]
    real, optional, intent(in)    :: u6(6)
    !< Velocity values on points i-3, i-2, i-1, i, i+1, i+2 [L T-1 ~> m s-1]
    real, intent(in)    :: h_tiny  !< A tiny thickness to prevent division by zero [L ~> m]
    real, intent(in)    :: u                      !< Velocity or thickness flux on point i-1/2
                                                  !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
    logical, intent(in) :: velocity_smoothing     !< If ture, use velocity to compute the smoothness indicator
    real, intent(inout) :: qr                     !< Reconstruction of tracer q on point i-1/2 [A ~> a]
    real :: vr                                    ! Reconstruction of hq [A ~> a]
    real :: hr                                    ! Reconstruction of h [L ~> m]
    real :: c0, c1, c2                            ! Intermediate reconstruction of hq[A ~> a]
    real :: d0, d1, d2                            ! Intermediate reconstruction of h [L ~> m]
    real :: b0, b1, b2                            ! Smoothness indicator [A ~> a]
    real :: tau                                   ! Difference of smoothness indicators [A ~> a]
    real :: w0, w1, w2                            ! Weights [nondim]
    real :: s                                     ! Temporary variables [nondim]
    real, parameter :: C3_10 = 3.0/10.0           ! The ratio of 3/10 [nondim]
    real, parameter :: C3_5 = 3.0/5.0             ! The ratio of 3/5 [nondim]
    real, parameter :: C1_10 = 1.0/10.0           ! The ratio of 1/10 [nondim]

    if (u>0.) then
      call weno_five_reconstruction_0(q6(3:5), c0) ! Reconstruction in the third upwind stencil
      call weno_five_reconstruction_1(q6(2:4), c1) ! Reconstruction in the second upwind stencil
      call weno_five_reconstruction_2(q6(1:3), c2) ! Reconstruction in the first upwind stencil

      call weno_five_reconstruction_0(h6(3:5), d0)
      call weno_five_reconstruction_1(h6(2:4), d1)
      call weno_five_reconstruction_2(h6(1:3), d2)
      if (velocity_smoothing) then
        call weno_five_weight_0(u6(3:5), b0)   ! Smoothness indicator of the third upwind stencil
        call weno_five_weight_1(u6(2:4), b1)   ! Smoothness indicator of the second upwind stencil
        call weno_five_weight_2(u6(1:3), b2)   ! Smoothness indicator of the first upwind stencil
      else
        call weno_five_weight_0(q6(3:5), b0)
        call weno_five_weight_1(q6(2:4), b1)
        call weno_five_weight_2(q6(1:3), b2)
      endif
    else
      call weno_five_reconstruction_0(q6(4:2:-1), c0) ! Reconstruction in the third upwind stencil
      call weno_five_reconstruction_1(q6(5:3:-1), c1) ! Reconstruction in the second upwind stencil
      call weno_five_reconstruction_2(q6(6:4:-1), c2) ! Reconstruction in the first upwind stencil

      call weno_five_reconstruction_0(h6(4:2:-1), d0)
      call weno_five_reconstruction_1(h6(5:3:-1), d1)
      call weno_five_reconstruction_2(h6(6:4:-1), d2)
      if (velocity_smoothing) then
        call weno_five_weight_0(u6(4:2:-1), b0) ! Smoothness indicator of the third upwind stencil
        call weno_five_weight_1(u6(5:3:-1), b1) ! Smoothness indicator of the second upwind stencil
        call weno_five_weight_2(u6(6:4:-1), b2) ! Smoothness indicator of the first upwind stencil
      else
        call weno_five_weight_0(q6(4:2:-1), b0)
        call weno_five_weight_1(q6(5:3:-1), b1)
        call weno_five_weight_2(q6(6:4:-1), b2)
      endif
    endif

    tau = abs(b0 - b2)
    w0  = C3_10 * fac_fn(tau, b0)
    w1  = C3_5  * fac_fn(tau, b1)
    w2  = C1_10 * fac_fn(tau, b2)

    s = 1. / ((w0 + w1) + w2)
    w0 = w0 * s   ! Weights of stencils
    w1 = w1 * s
    w2 = w2 * s

    vr = ((w0 * c0) + (w1 * c1)) + (w2 * c2)
    hr = ((w0 * d0) + (w1 * d1)) + (w2 * d2)
!    vr = min(max(q6(3), q6(4)), vr) ; vr = max(min(q6(3), q6(4)), vr) !Impose a monotonicity limiter
    hr = min(max(h6(3), h6(4)), hr) ; hr = max(min(h6(3), h6(4)), hr) !Impose a monotonicity limiter

    qr = vr / max(hr, h_tiny)

end subroutine weno_five_h_weight_reconstruction

!> Compute the smoothness indicator for the third upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_weight_0(q3, w0)
  real, intent(in) :: q3(3)       !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: w0       !< Smoothness indicator for this stencil [A2 ~> a2]

  w0 = (q3(1) * ((10 * q3(1) - 31 * q3(2)) + 11 * q3(3))) + &
       ((q3(2) * (25 * q3(2) - 19 * q3(3))) + 4 * (q3(3) * q3(3)))

end subroutine weno_five_weight_0

!> Compute the smoothness indicator for the second upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_weight_1(q3, w1)
  real, intent(in) :: q3(3)        !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: w1        !< Smoothness indicator for this stencil [A2 ~> a2]

  w1 = (q3(1) * ((4 * q3(1) - 13 * q3(2)) + 5 * q3(3))) + &
       ((q3(2) * (13 * q3(2) - 13 * q3(3))) + 4 * (q3(3) * q3(3)))

end subroutine weno_five_weight_1

!> Compute the smoothness indicator for the first upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_weight_2(q3, w2)
  real, intent(in) :: q3(3)        !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: w2        !< Smoothness indicator for this stencil [A2 ~> a2]

  w2 = (q3(1) * ((4 * q3(1) - 19 * q3(2)) + 11 * q3(3))) + &
       ((q3(2) * (25 * q3(2) - 31 * q3(3))) + 10 * (q3(3) * q3(3)))

end subroutine weno_five_weight_2

!> Reconstruction in the third upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_reconstruction_0(q3, p0)
  real, intent(in) :: q3(3)        !< Tracer values on three points [A ~> a]
  real, intent(inout) :: p0        !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_6 = 1.0/6.0 ! One sixth [nondim]

  p0 = ((2*q3(1) + 5*q3(2)) - q3(3)) * C1_6

end subroutine weno_five_reconstruction_0

!> Reconstruction in the second upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_reconstruction_1(q3, p1)
  real, intent(in) :: q3(3)         !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: p1         !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_6 = 1.0/6.0 ! One sixth [nondim]

  p1 = ((-q3(1) + 5*q3(2)) + 2*q3(3)) * C1_6

end subroutine weno_five_reconstruction_1

!> Reconstruction in the first upwind stencil of the fifth-order WENO scheme
pure subroutine weno_five_reconstruction_2(q3, p2)
  real, intent(in) :: q3(3)          !< Tracer values on the three-point stencil [A ~> a]
  real, intent(inout) :: p2          !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_6 = 1.0/6.0  ! One sixth [nondim]

  p2 = ((2*q3(1) - 7*q3(2)) + 11*q3(3)) * C1_6

end subroutine weno_five_reconstruction_2


!> Reconstruct the tracer (e.g., PV, vorticity) onto point i-1/2 using a seventh-order WENO scheme
!! This reconstruction computes a thickness weighted average of PV
pure subroutine weno_seven_h_weight_reconstruction(q8, h8, u8, &
                                            h_tiny, u, qr, velocity_smoothing)
  real, intent(in)    :: q8(8)
  !< Tracer values on points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3
  real, intent(in)    :: h8(8)
  !< Thickness on the same tracer points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3 [L ~> m]
  real, optional, intent(in)    :: u8(8)
  !< Velocity values on points i-4, i-3, i-2, i-1, i, i+1, i+2, i+3 [L T-1 ~> m s-1]
  real, intent(in)    :: h_tiny  !< A tiny thickness to prevent division by zero [L ~> m]
  real, intent(in)    :: u    !< Velocity or thickness flux on point i-1/2
                              !! [L T-1 ~> m s-1] or [L2 T-1 ~> m2 s-1]
  logical, intent(in) :: velocity_smoothing !< If true, use velocity to compute the smoothness indicator
  real, intent(inout) :: qr   !< Reconstruction of tracer q on point i-1/2 [A ~> a]
  real :: vr                  ! Reconstruction of hq [A ~> a]
  real :: hr                  ! Reconstruction of h [L ~> m]
  real :: c0, c1, c2, c3      ! Intermediate reconstruction of hq [A ~> a]
  real :: d0, d1, d2, d3      ! Intermediate reconstruction of h [L ~> m]
  real :: b0, b1, b2, b3      ! Smoothness indicator [A ~> a]
  real :: tau                 ! Difference of smoothness indicators [A ~> a]
  real :: w0, w1, w2, w3      ! Weights [nondim]
  real :: s                   ! Temporary variables [nondim]
  real, parameter :: C4_35 = 4.0/35.0 ! The ratio of 4/35 [nondim]
  real, parameter :: C18_35 = 18.0/35.0 ! The ratio of 18/35 [nondim]
  real, parameter :: C12_35 = 12.0/35.0 ! The ratio of 12/35 [nondim]
  real, parameter :: C1_35 = 1.0/35.0   ! The ratio of 1/35 [nondim]

  if (u>0.) then
    call weno_seven_reconstruction_0(q8(4:7), c0) ! Reconstruction in the fourth upwind stencil
    call weno_seven_reconstruction_1(q8(3:6), c1) ! Reconstruction in the third upwind stencil
    call weno_seven_reconstruction_2(q8(2:5), c2) ! Reconstruction in the second upwind stencil
    call weno_seven_reconstruction_3(q8(1:4), c3) ! Reconstruction in the first upwind stencil

    call weno_seven_reconstruction_0(h8(4:7), d0) ! Reconstruction in the fourth upwind stencil
    call weno_seven_reconstruction_1(h8(3:6), d1) ! Reconstruction in the third upwind stencil
    call weno_seven_reconstruction_2(h8(2:5), d2) ! Reconstruction in the second upwind stencil
    call weno_seven_reconstruction_3(h8(1:4), d3) ! Reconstruction in the first upwind stencil
    if (velocity_smoothing) then
      call weno_seven_weight_0(u8(4:7), b0)       ! Smoothness indicator of the fourth upwind stencil
      call weno_seven_weight_1(u8(3:6), b1)       ! Smoothness indicator of the third upwind stencil
      call weno_seven_weight_2(u8(2:5), b2)       ! Smoothness indicator of the second upwind stencil
      call weno_seven_weight_3(u8(1:4), b3)       ! Smoothness indicator of the first upwind stencil
    else
      call weno_seven_weight_0(q8(4:7), b0)
      call weno_seven_weight_1(q8(3:6), b1)
      call weno_seven_weight_2(q8(2:5), b2)
      call weno_seven_weight_3(q8(1:4), b3)
    endif
  else
    call weno_seven_reconstruction_0(q8(5:2:-1), c0) ! Reconstruction in the fourth upwind stencil
    call weno_seven_reconstruction_1(q8(6:3:-1), c1) ! Reconstruction in the third upwind stencil
    call weno_seven_reconstruction_2(q8(7:4:-1), c2) ! Reconstruction in the second upwind stencil
    call weno_seven_reconstruction_3(q8(8:5:-1), c3) ! Reconstruction in the first upwind stencil

    call weno_seven_reconstruction_0(h8(5:2:-1), d0)
    call weno_seven_reconstruction_1(h8(6:3:-1), d1)
    call weno_seven_reconstruction_2(h8(7:4:-1), d2)
    call weno_seven_reconstruction_3(h8(8:5:-1), d3)
    if (velocity_smoothing) then
      call weno_seven_weight_0(u8(5:2:-1), b0)    ! Smoothness indicator of the fourth upwind stencil
      call weno_seven_weight_1(u8(6:3:-1), b1)    ! Smoothness indicator of the third upwind stencil
      call weno_seven_weight_2(u8(7:4:-1), b2)    ! Smoothness indicator of the second upwind stencil
      call weno_seven_weight_3(u8(8:5:-1), b3)    ! Smoothness indicator of the first upwind stencil
    else
      call weno_seven_weight_0(q8(5:2:-1), b0)
      call weno_seven_weight_1(q8(6:3:-1), b1)
      call weno_seven_weight_2(q8(7:4:-1), b2)
      call weno_seven_weight_3(q8(8:5:-1), b3)
    endif
  endif

  tau = abs((b0 - b3) + 3 * (b1 - b2))
  w0  = C4_35  * fac_fn(tau, b0)
  w1  = C18_35 * fac_fn(tau, b1)
  w2  = C12_35 * fac_fn(tau, b2)
  w3  = C1_35  * fac_fn(tau, b3)

  s = 1. / ((w0 + w1) + (w2 + w3))
  w0 = w0 * s   ! Weights of the stencils
  w1 = w1 * s
  w2 = w2 * s
  w3 = w3 * s

  vr = ((w0 * c0) + (w1 * c1)) + ((w2 * c2) + (w3 * c3))
  hr = ((w0 * d0) + (w1 * d1)) + ((w2 * d2) + (w3 * d3))

!  vr = min(max(q4, q5), vr) ; vr = max(min(q4, q5), vr)
  hr = min(max(h8(4), h8(5)), hr) ; hr = max(min(h8(4), h8(5)), hr) ! Impose a monotonicity limiter

  qr = vr / max(hr, h_tiny)

end subroutine weno_seven_h_weight_reconstruction

!> Compute the smoothness indicator for the fourth upwind stencil of the seventh-order WENO scheme
pure subroutine weno_seven_weight_0(q4, w0)
  real, intent(in) :: q4(4)          !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w0          !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
  w0 = ((q4(1) * ((2.107 * q4(1) - 9.402 * q4(2)) + (7.042 * q4(3) - 1.854 * q4(4)))) + &
     (q4(2) * ((11.003 * q4(2) - 17.246 * q4(3)) + 4.642 * q4(4)))) + &
     ((q4(3) * (7.043 * q4(3) - 3.882 * q4(4))) + 0.547 * (q4(4) * q4(4)))

end subroutine weno_seven_weight_0

!> Compute the smoothness indicator for the third upwind stencil of the seventh-order WENO scheme
pure subroutine weno_seven_weight_1(q4, w1)
  real, intent(in) :: q4(4)          !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w1          !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
  w1 = ((q4(1) * ((0.547 * q4(1) - 2.522 * q4(2)) + (1.922 * q4(3) - 0.494 * q4(4)))) + &
     (q4(2) * ((3.443 * q4(2) - 5.966 * q4(3)) + 1.602 * q4(4)))) + &
     ((q4(3) * (2.843 * q4(3) - 1.642 * q4(4))) + 0.267 * (q4(4) * q4(4)))

end subroutine weno_seven_weight_1

!> Compute the smoothness indicator for the second upwind stencil of the seventh-order WENO scheme
pure subroutine weno_seven_weight_2(q4, w2)
  real, intent(in) :: q4(4)           !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w2           !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
  w2 = ((q4(1) * ((0.267 * q4(1) - 1.642 * q4(2)) + (1.602 * q4(3) - 0.494 * q4(4)))) + &
     (q4(2) * ((2.843 * q4(2) - 5.966 * q4(3)) + 1.922 * q4(4)))) + &
     ((q4(3) * (3.443 * q4(3) - 2.522 * q4(4))) + 0.547 * (q4(4) * q4(4)))

end subroutine weno_seven_weight_2

!> Compute smoothness indicator for the first upwind stencil of the seventh-order WENO scheme
pure subroutine weno_seven_weight_3(q4, w3)
  real, intent(in) :: q4(4)           !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: w3           !< Smoothness indicator for this stencil [A2 ~> a2]

  ! Coefficients from Balsara and Shu (2000). The division by 1000 will be normalized out by fac_fn
  w3 = ((q4(1) * ((0.547  * q4(1) - 3.882 * q4(2)) + (4.642 * q4(3) - 1.854 * q4(4)))) + &
     (q4(2) * ((7.043 * q4(2) - 17.246 * q4(3)) + 7.042 * q4(4)))) + &
     ((q4(3) * (11.003 * q4(3) - 9.402 * q4(4))) + 2.107 * (q4(4) * q4(4)))

end subroutine weno_seven_weight_3

!> Reconstruction in the fourth upwind stencil for seventh-order WENO scheme
pure subroutine weno_seven_reconstruction_0(q4, p0)
  real, intent(in) :: q4(4)            !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p0            !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_24 = 1.0/24.0  ! One twenty fourth [nondim]

  p0 = (((6 * q4(1) + 26 * q4(2)) - 10 * q4(3)) + 2 * q4(4)) * C1_24

end subroutine weno_seven_reconstruction_0

!> Reconstruction in the third upwind stencil for seventh-order WENO scheme
pure subroutine weno_seven_reconstruction_1(q4, p1)
  real, intent(in) :: q4(4)            !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p1            !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_24 = 1.0/24.0  ! One twenty fourth [nondim]

  p1 = (14 * (q4(2) + q4(3)) - 2 * (q4(1) + q4(4))) * C1_24

end subroutine weno_seven_reconstruction_1

!> Reconstruction in the second upwind stencil for seventh-order WENO scheme
pure subroutine weno_seven_reconstruction_2(q4, p2)
  real, intent(in) :: q4(4)             !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p2             !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_24 = 1.0/24.0   ! One twenty fourth [nondim]

  p2 = (((2 * q4(1) - 10 * q4(2)) + 26 * q4(3)) + 6 * q4(4)) * C1_24

end subroutine weno_seven_reconstruction_2

!> Reconstruction in the first upwind stencil for seventh-order WENO scheme
pure subroutine weno_seven_reconstruction_3(q4, p3)
  real, intent(in) :: q4(4)            !< Tracer values on the four-point stencil [A ~> a]
  real, intent(inout) :: p3            !< Reconstruction of the quantity [A ~> a]
  real, parameter :: C1_24 = 1.0/24.0  ! One twenty fourth [nondim]

  p3 = (((-6 * q4(1) + 26 * q4(2)) - 46 * q4(3)) + 50 * q4(4)) * C1_24

end subroutine weno_seven_reconstruction_3

function CoriolisAdv_stencil(CS) result(stencil)
  type(CoriolisAdv_CS), intent(in)  :: CS  !< Control structure for MOM_CoriolisAdv
  integer :: stencil  !< The halo stencil size for the Coriolis advection scheme

  stencil = Coriolis_scheme_stencil(CS%CorAdv_opts_CS%Coriolis_scheme_CS)

end function CoriolisAdv_stencil

!> Returns the halo stencil size needed by the selected Coriolis scheme.
function Coriolis_scheme_stencil(CS) result(stencil)
  type(Coriolis_scheme_CS), intent(in) :: CS !< Options selecting the discretization of the
                                             !! Coriolis terms.
  integer :: stencil  !< The halo stencil size for the Coriolis advection scheme

  stencil = 2
  if (CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO) stencil = 4
  if (CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO) stencil = 3

end function Coriolis_scheme_stencil


!> Initializes the control structure for MOM_CoriolisAdv
subroutine CoriolisAdv_init(Time, G, GV, US, param_file, diag, AD, CS)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Runtime parameter handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(accel_diag_ptrs),   target, intent(inout) :: AD !< Storage for acceleration diagnostics
  type(CoriolisAdv_CS),    intent(inout) :: CS   !< Control structure for MOM_CoriolisAdv
  ! Local variables
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_CoriolisAdv" ! This module's name.
  character(len=20)  :: tmpstr
  character(len=400) :: mesg
  logical :: use_weno
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
#ifdef __NVCOMPILER_OPENMP_GPU
  integer, parameter :: default_nkblock = 0
#else
  integer, parameter :: default_nkblock = 1
#endif

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%initialized = .true.
  CS%diag => diag ; CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "CORIOLIS_ADV_NKBLOCK", CS%CorAdv_opts_CS%nkblock, &
                 "The k-direction block size used in Coriolis and momentum advection "//&
                 "calculations. The default 0 setting dynamically uses the full vertical column.", &
                 default=default_nkblock, layoutParam=.true.)
  if (CS%CorAdv_opts_CS%nkblock < 0) call MOM_error(FATAL, "CORIOLIS_ADV_NKBLOCK must be >= 0.")
  call get_param(param_file, mdl, "NOSLIP", CS%CorAdv_opts_CS%no_slip, &
                 "If true, no slip boundary conditions are used; otherwise "//&
                 "free slip boundary conditions are assumed. The "//&
                 "implementation of the free slip BCs on a C-grid is much "//&
                 "cleaner than the no slip BCs. The use of free slip BCs "//&
                 "is strongly encouraged, and no slip BCs are not used with "//&
                 "the biharmonic viscosity.", default=.false.)

  call get_param(param_file, mdl, "CORIOLIS_EN_DIS", &
                 CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_En_Dis, &
                 "If true, two estimates of the thickness fluxes are used "//&
                 "to estimate the Coriolis term, and the one that "//&
                 "dissipates energy relative to the other one is used.", &
                 default=.false.)

  ! Set %Coriolis_Scheme
  ! (Select the baseline discretization for the Coriolis term)
  call get_param(param_file, mdl, "CORIOLIS_SCHEME", tmpstr, &
                 "CORIOLIS_SCHEME selects the discretization for the "//&
                 "Coriolis terms. Valid values are: \n"//&
                 "\t SADOURNY75_ENERGY - Sadourny, 1975; energy cons. \n"//&
                 "\t ARAKAWA_HSU90     - Arakawa & Hsu, 1990 \n"//&
                 "\t SADOURNY75_ENSTRO - Sadourny, 1975; enstrophy cons. \n"//&
                 "\t ARAKAWA_LAMB81    - Arakawa & Lamb, 1981; En. + Enst.\n"//&
                 "\t ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with \n"//&
                 "\t                      Arakawa & Hsu and Sadourny energy \n"//&
                 "\t WENOVI5TH_PV_ENSTRO   - 5th-order WENO PV enstrophy \n"//&
                 "\t WENOVI3RD_PV_ENSTRO   - 3rd-order WENO PV enstrophy \n"//&
                 "\t WENOVI7TH_PV_ENSTRO  - 7th-order WENO PV enstrophy \n", &
                 default=SADOURNY75_ENERGY_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (SADOURNY75_ENERGY_STRING)
      CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme = SADOURNY75_ENERGY
    case (ARAKAWA_HSU_STRING)
      CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme = ARAKAWA_HSU90
    case (SADOURNY75_ENSTRO_STRING)
      CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme = SADOURNY75_ENSTRO
    case (ARAKAWA_LAMB_STRING)
      CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme = ARAKAWA_LAMB81
    case (AL_BLEND_STRING)
      CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme = AL_BLEND
    case (ROBUST_ENSTRO_STRING)
      CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme = ROBUST_ENSTRO
      CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_En_Dis = .false.
    case (WENOVI7TH_PV_ENSTRO_STRING)
      CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme = wenovi7th_PV_ENSTRO
    case (WENOVI5TH_PV_ENSTRO_STRING)
      CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme = wenovi5th_PV_ENSTRO
    case (WENOVI3RD_PV_ENSTRO_STRING)
      CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme = wenovi3rd_PV_ENSTRO
    case default
      call MOM_mesg('CoriolisAdv_init: Coriolis_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: Unrecognized setting "// &
            "#define CORIOLIS_SCHEME "//trim(tmpstr)//" found in input file.")
  end select

  use_weno = CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme == wenovi7th_PV_ENSTRO &
      .or. CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme == wenovi5th_PV_ENSTRO &
      .or. CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme == wenovi3rd_PV_ENSTRO

  if (use_weno) then
    call get_param(param_file, mdl, "WENO_VELOCITY_SMOOTH", &
                   CS%CorAdv_opts_CS%weno_velocity_smooth, &
            "If true, use velocity to compute weighting for WENO. ", &
                  default=.false.)
  endif

  if (CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme == AL_BLEND) then
    call get_param(param_file, mdl, "CORIOLIS_BLEND_WT_LIN", CS%CorAdv_opts_CS%wt_lin_blend, &
                 "A weighting value for the ratio of inverse thicknesses, "//&
                 "beyond which the blending between Sadourny Energy and "//&
                 "Arakawa & Hsu goes linearly to 0 when CORIOLIS_SCHEME "//&
                 "is ARAWAKA_LAMB_BLEND. This must be between 1 and 1e-16.", &
                 units="nondim", default=0.125)
    call get_param(param_file, mdl, "CORIOLIS_BLEND_F_EFF_MAX", CS%CorAdv_opts_CS%F_eff_max_blend, &
                 "The factor by which the maximum effective Coriolis "//&
                 "acceleration from any point can be increased when "//&
                 "blending different discretizations with the "//&
                 "ARAKAWA_LAMB_BLEND Coriolis scheme.  This must be "//&
                 "greater than 2.0 (the max value for Sadourny energy).", &
                 units="nondim", default=4.0)
    CS%CorAdv_opts_CS%wt_lin_blend = min(1.0, max(CS%CorAdv_opts_CS%wt_lin_blend,1e-16))
    if (CS%CorAdv_opts_CS%F_eff_max_blend < 2.0) call MOM_error(WARNING, "CoriolisAdv_init: "//&
           "CORIOLIS_BLEND_F_EFF_MAX should be at least 2.")
  endif

  mesg = "If true, the Coriolis terms at u-points are bounded by "//&
         "the four estimates of (f+rv)v from the four neighboring "//&
         "v-points, and similarly at v-points."
  if (CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_En_Dis .and. &
      (CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme == SADOURNY75_ENERGY)) then
    mesg = trim(mesg)//"  This option is "//&
                 "always effectively false with CORIOLIS_EN_DIS defined and "//&
                 "CORIOLIS_SCHEME set to "//trim(SADOURNY75_ENERGY_STRING)//"."
  else
    mesg = trim(mesg)//"  This option would "//&
                 "have no effect on the SADOURNY Coriolis scheme if it "//&
                 "were possible to use centered difference thickness fluxes."
  endif
  call get_param(param_file, mdl, "BOUND_CORIOLIS", CS%CorAdv_opts_CS%bound_Coriolis, mesg, &
                 default=.false.)
  if ((CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_En_Dis .and. &
       (CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme == SADOURNY75_ENERGY)) .or. &
      (CS%CorAdv_opts_CS%Coriolis_scheme_CS%Coriolis_Scheme == ROBUST_ENSTRO)) &
    CS%CorAdv_opts_CS%bound_Coriolis = .false.

  ! Set KE_Scheme (selects discretization of KE)
  call get_param(param_file, mdl, "KE_SCHEME", tmpstr, &
                 "KE_SCHEME selects the discretization for acceleration "//&
                 "due to the kinetic energy gradient. Valid values are: \n"//&
                 "\t KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV, KE_UP3", &
                 default=KE_ARAKAWA_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (KE_ARAKAWA_STRING); CS%CorAdv_opts_CS%KE_Scheme = KE_ARAKAWA
    case (KE_SIMPLE_GUDONOV_STRING); CS%CorAdv_opts_CS%KE_Scheme = KE_SIMPLE_GUDONOV
    case (KE_GUDONOV_STRING); CS%CorAdv_opts_CS%KE_Scheme = KE_GUDONOV
    case (KE_UP3_STRING); CS%CorAdv_opts_CS%KE_Scheme = KE_UP3
    case default
      call MOM_mesg('CoriolisAdv_init: KE_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: "// &
               "#define KE_SCHEME "//trim(tmpstr)//" in input file is invalid.")
  end select

  if (CS%CorAdv_opts_CS%KE_Scheme == KE_UP3) then
    call get_param(param_file, mdl, "KE_USE_LIMITER", CS%CorAdv_opts_CS%KE_use_limiter, &
            "If true, use Koren limiter for KE_UP3 scheme", &
                  default=.True.)
  endif

  ! Set PV_Adv_Scheme (selects discretization of PV advection)
  call get_param(param_file, mdl, "PV_ADV_SCHEME", tmpstr, &
                 "PV_ADV_SCHEME selects the discretization for PV "//&
                 "advection. Valid values are: \n"//&
                 "\t PV_ADV_CENTERED - centered (aka Sadourny, 75) \n"//&
                 "\t PV_ADV_UPWIND1  - upwind, first order", &
                 default=PV_ADV_CENTERED_STRING)
  select case (uppercase(tmpstr))
    case (PV_ADV_CENTERED_STRING)
      CS%CorAdv_opts_CS%PV_Adv_Scheme = PV_ADV_CENTERED
    case (PV_ADV_UPWIND1_STRING)
      CS%CorAdv_opts_CS%PV_Adv_Scheme = PV_ADV_UPWIND1
    case default
      call MOM_mesg('CoriolisAdv_init: PV_Adv_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: "// &
                     "#DEFINE PV_ADV_SCHEME in input file is invalid.")
  end select

  CS%id_rv = register_diag_field('ocean_model', 'RV', diag%axesBL, Time, &
     'Relative Vorticity', 's-1', conversion=US%s_to_T)

  CS%id_PV = register_diag_field('ocean_model', 'PV', diag%axesBL, Time, &
     'Potential Vorticity', 'm-1 s-1', conversion=GV%m_to_H*US%s_to_T)

  CS%id_gKEu = register_diag_field('ocean_model', 'gKEu', diag%axesCuL, Time, &
     'Zonal Acceleration from Grad. Kinetic Energy', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_gKEv = register_diag_field('ocean_model', 'gKEv', diag%axesCvL, Time, &
     'Meridional Acceleration from Grad. Kinetic Energy', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_rvxu = register_diag_field('ocean_model', 'rvxu', diag%axesCvL, Time, &
     'Meridional Acceleration from Relative Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_rvxv = register_diag_field('ocean_model', 'rvxv', diag%axesCuL, Time, &
     'Zonal Acceleration from Relative Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_CAuS = register_diag_field('ocean_model', 'CAu_Stokes', diag%axesCuL, Time, &
     'Zonal Acceleration from Stokes Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)
  ! add to AD

  CS%id_CAvS = register_diag_field('ocean_model', 'CAv_Stokes', diag%axesCvL, Time, &
     'Meridional Acceleration from Stokes Vorticity', 'm s-2', conversion=US%L_T2_to_m_s2)
  ! add to AD

  !CS%id_hf_gKEu = register_diag_field('ocean_model', 'hf_gKEu', diag%axesCuL, Time, &
  !   'Fractional Thickness-weighted Zonal Acceleration from Grad. Kinetic Energy', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_gKEu_2d = register_diag_field('ocean_model', 'hf_gKEu_2d', diag%axesCu1, Time, &
     'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Grad. Kinetic Energy', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  !CS%id_hf_gKEv = register_diag_field('ocean_model', 'hf_gKEv', diag%axesCvL, Time, &
  !   'Fractional Thickness-weighted Meridional Acceleration from Grad. Kinetic Energy', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_gKEv_2d = register_diag_field('ocean_model', 'hf_gKEv_2d', diag%axesCv1, Time, &
     'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Grad. Kinetic Energy', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_h_gKEu = register_diag_field('ocean_model', 'h_gKEu', diag%axesCuL, Time, &
     'Thickness Multiplied Zonal Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_gKEu_2d = register_diag_field('ocean_model', 'intz_gKEu_2d', diag%axesCu1, Time, &
     'Depth-integral of Zonal Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  CS%id_h_gKEv = register_diag_field('ocean_model', 'h_gKEv', diag%axesCvL, Time, &
     'Thickness Multiplied Meridional Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_gKEv_2d = register_diag_field('ocean_model', 'intz_gKEv_2d', diag%axesCv1, Time, &
     'Depth-integral of Meridional Acceleration from Grad. Kinetic Energy', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  !CS%id_hf_rvxu = register_diag_field('ocean_model', 'hf_rvxu', diag%axesCvL, Time, &
  !   'Fractional Thickness-weighted Meridional Acceleration from Relative Vorticity', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_rvxu_2d = register_diag_field('ocean_model', 'hf_rvxu_2d', diag%axesCv1, Time, &
     'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Relative Vorticity', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  !CS%id_hf_rvxv = register_diag_field('ocean_model', 'hf_rvxv', diag%axesCuL, Time, &
  !   'Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
  !   'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  CS%id_hf_rvxv_2d = register_diag_field('ocean_model', 'hf_rvxv_2d', diag%axesCu1, Time, &
     'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
     'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_h_rvxu = register_diag_field('ocean_model', 'h_rvxu', diag%axesCvL, Time, &
     'Thickness Multiplied Meridional Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_rvxu_2d = register_diag_field('ocean_model', 'intz_rvxu_2d', diag%axesCv1, Time, &
     'Depth-integral of Meridional Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  CS%id_h_rvxv = register_diag_field('ocean_model', 'h_rvxv', diag%axesCuL, Time, &
     'Thickness Multiplied Zonal Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  CS%id_intz_rvxv_2d = register_diag_field('ocean_model', 'intz_rvxv_2d', diag%axesCu1, Time, &
     'Depth-integral of Fractional Thickness-weighted Zonal Acceleration from Relative Vorticity', &
     'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  ! Allocate memory needed for the diagnostics that have been enabled.
  if ((CS%id_gKEu > 0) .or. (CS%id_hf_gKEu_2d > 0) .or. &
    ! (CS%id_hf_gKEu > 0) .or. &
      (CS%id_h_gKEu > 0) .or. (CS%id_intz_gKEu_2d > 0)) then
    call safe_alloc_ptr(AD%gradKEu, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_gKEv > 0) .or. (CS%id_hf_gKEv_2d > 0) .or. &
    ! (CS%id_hf_gKEv > 0) .or. &
      (CS%id_h_gKEv > 0) .or. (CS%id_intz_gKEv_2d > 0)) then
    call safe_alloc_ptr(AD%gradKEv, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_rvxu > 0) .or. (CS%id_hf_rvxu_2d > 0) .or. &
    ! (CS%id_hf_rvxu > 0) .or. &
      (CS%id_h_rvxu > 0) .or. (CS%id_intz_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%rv_x_u, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_rvxv > 0) .or. (CS%id_hf_rvxv_2d > 0) .or. &
    ! (CS%id_hf_rvxv > 0) .or. &
      (CS%id_h_rvxv > 0) .or. (CS%id_intz_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%rv_x_v, IsdB, IedB, jsd, jed, nz)
  endif

  if ((CS%id_hf_gKEv_2d > 0) .or. &
    ! (CS%id_hf_gKEv > 0) .or. (CS%id_hf_rvxu > 0) .or. &
      (CS%id_hf_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hfrac_v, isd, ied, JsdB, JedB, nz)
  endif
  if ((CS%id_hf_gKEu_2d > 0) .or. &
    ! (CS%id_hf_gKEu > 0) .or. (CS%id_hf_rvxv > 0) .or. &
      (CS%id_hf_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hfrac_u, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_h_gKEu > 0) .or. (CS%id_intz_gKEu_2d > 0) .or. &
      (CS%id_h_rvxv > 0) .or. (CS%id_intz_rvxv_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hu, IsdB, IedB, jsd, jed, nz)
  endif
  if ((CS%id_h_gKEv > 0) .or. (CS%id_intz_gKEv_2d > 0) .or. &
      (CS%id_h_rvxu > 0) .or. (CS%id_intz_rvxu_2d > 0)) then
    call safe_alloc_ptr(AD%diag_hv, isd, ied, JsdB, JedB, nz)
  endif

  CS%CorAdv_opts_CS%do_CAS_diag = (CS%id_CAuS > 0) .or. (CS%id_CAvS > 0)
  CS%CorAdv_opts_CS%do_RV_diag = (CS%id_rv > 0)
  CS%CorAdv_opts_CS%do_PV_diag = (CS%id_PV > 0)

end subroutine CoriolisAdv_init

!> Destructor for coriolisadv_cs
subroutine CoriolisAdv_end(CS)
  type(CoriolisAdv_CS), intent(inout) :: CS !< Control structure for MOM_CoriolisAdv
end subroutine CoriolisAdv_end

!> \namespace mom_coriolisadv
!!
!! This file contains the subroutine that calculates the time
!! derivatives of the velocities due to Coriolis acceleration and
!! momentum advection.  This subroutine uses either a vorticity
!! advection scheme from Arakawa and Hsu, Mon. Wea. Rev. 1990, or
!! Sadourny's (JAS 1975) energy conserving scheme.  Both have been
!! modified to use general orthogonal coordinates as described in
!! Arakawa and Lamb, Mon. Wea. Rev. 1981.  Both schemes are second
!! order accurate, and allow for vanishingly small layer thicknesses.
!! The Arakawa and Hsu scheme globally conserves both total energy
!! and potential enstrophy in the limit of nondivergent flow.
!! Sadourny's energy conserving scheme conserves energy if the flow
!! is nondivergent or centered difference thickness fluxes are used.
!!
!! A small fragment of the grid is shown below:
!! \verbatim
!!
!!    j+1  x ^ x ^ x   At x:  q, CoriolisBu
!!    j+1  > o > o >   At ^:  v, CAv, vh
!!    j    x ^ x ^ x   At >:  u, CAu, uh, a, b, c, d
!!    j    > o > o >   At o:  h, KE
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!! \endverbatim
!!
!! The boundaries always run through q grid points (x).

end module MOM_CoriolisAdv
