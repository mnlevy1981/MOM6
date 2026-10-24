! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!!SKILLS: 0.4.0

#include "do_concurrent_compat.h"

!> Calculates horizontal viscosity and viscous stresses
module MOM_hor_visc

use array_mod,                 only : RealArray_t, RealArray_C
use iso_c_binding,             only : c_bool, c_double, c_int, c_ptr, c_loc, c_f_pointer
use turbotmp_helperF,          only : io_recorder
use box_mod,                   only : Box_t
use MOM_checksums,             only : hchksum, Bchksum, uvchksum
use MOM_coms,                  only : min_across_PEs
use MOM_diag_mediator,         only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,         only : post_product_u, post_product_sum_u
use MOM_diag_mediator,         only : post_product_v, post_product_sum_v
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_domains,               only : pass_var, CORNER, pass_vector, AGRID, BGRID_NE
use MOM_domains,               only : To_All, Scalar_Pair
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_grid_containers,       only : grid_core_type, grid_core
use MOM_interface_heights,     only : thickness_to_dz
use MOM_lateral_mixing_coeffs, only : VarMix_CS, calc_QG_slopes, calc_QG_Leith_viscosity
use MOM_barotropic,            only : barotropic_CS, barotropic_get_tav
use MOM_thickness_diffuse,     only : thickness_diffuse_CS, thickness_diffuse_get_KH
use MOM_io,                    only : MOM_read_data, slasher
use MOM_MEKE_types,            only : MEKE_type
use MOM_open_boundary,         only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary,         only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_open_boundary,         only : OBC_STRAIN_NONE, OBC_STRAIN_ZERO, OBC_STRAIN_FREESLIP
use MOM_open_boundary,         only : OBC_STRAIN_COMPUTED, OBC_STRAIN_SPECIFIED
use MOM_stochastics,           only : stochastic_CS
use MOM_unit_scaling,          only : unit_scale_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_variables,             only : accel_diag_ptrs, thermo_var_ptrs
use MOM_Zanna_Bolton,          only : ZB2020_lateral_stress, ZB2020_init, ZB2020_end
use MOM_Zanna_Bolton,          only : ZB2020_CS, ZB2020_copy_gradient_and_thickness

implicit none ; private

#include <MOM_memory.h>

public horizontal_viscosity, hor_visc_nkblock, hor_visc_init, hor_visc_end, hor_visc_vel_stencil


!> Options and diagnostic switches used by the horizontal viscosity kernels.
type :: hv_opts_CS
  logical :: Laplacian       !< Use a Laplacian horizontal viscosity if true.
  logical :: biharmonic      !< Use a biharmonic horizontal viscosity if true.
  logical :: debug           !< If true, write verbose checksums for debugging purposes.
  logical :: no_slip         !< If true, no slip boundary conditions are used.
                             !! Otherwise free slip boundary conditions are assumed.
                             !! The implementation of the free slip boundary
                             !! conditions on a C-grid is much cleaner than the
                             !! no slip boundary conditions. The use of free slip
                             !! b.c.s is strongly encouraged. The no slip b.c.s
                             !! are not implemented with the biharmonic viscosity.
  logical :: bound_Kh        !< If true, the Laplacian coefficient is locally
                             !! limited to guarantee stability.
  logical :: EY24_EBT_BS     !! If true, use an equivalent barotropic backscatter
                             !! with a stabilizing kill switch in MEKE,
                             !< developed by Yankovsky et al. 2024
  logical :: bound_Ah        !< If true, the biharmonic coefficient is locally
                             !! limited to guarantee stability.
  real    :: Re_Ah           !! If nonzero, the biharmonic coefficient is scaled
                             !< so that the biharmonic Reynolds number is equal to this [nondim].
  real    :: KS_coef         !< A nondimensional coefficient on the biharmonic viscosity that sets the
                             !! kill switch for backscatter. Default is 1.0 [nondim].
  logical :: backscatter_underbound !< If true, the bounds on the biharmonic viscosity are allowed
                             !! to increase where the Laplacian viscosity is negative (due to
                             !! backscatter parameterizations) beyond the largest timestep-dependent
                             !! stable values of biharmonic viscosity when no Laplacian viscosity is
                             !! applied.  The default is true for historical reasons, but this option
                             !! probably should not be used as it can lead to numerical instabilities.
  logical :: Smagorinsky_Kh  !< If true, use Smagorinsky nonlinear eddy
                             !! viscosity. KH is the background value.
  logical :: Smagorinsky_Ah  !< If true, use a biharmonic form of Smagorinsky
                             !! nonlinear eddy viscosity. AH is the background.
  logical :: Leith_Kh        !< If true, use 2D Leith nonlinear eddy
                             !! viscosity. KH is the background value.
  logical :: Modified_Leith  !< If true, use extra component of Leith viscosity
                             !! to damp divergent flow. To use, still set Leith_Kh=.TRUE.
  logical :: use_beta_in_Leith !< If true, includes the beta term in the Leith viscosity
  logical :: Leith_Ah        !< If true, use a biharmonic form of 2D Leith
                             !! nonlinear eddy viscosity. AH is the background.
  logical :: use_Leithy      !< If true, use a biharmonic form of 2D Leith
                             !! nonlinear eddy viscosity with harmonic backscatter.
                             !! Ah is the background. Leithy = Leith+E
  real    :: c_K             !< Fraction of energy dissipated by the biharmonic term
                             !! that gets backscattered in the Leith+E scheme. [nondim]
  logical :: smooth_Ah       !< If true (default), then Ah and m_leithy are smoothed.
                             !! This smoothing requires a lot of blocking communication.
  logical :: use_QG_Leith_visc    !< If true, use QG Leith nonlinear eddy viscosity.
                             !! KH is the background value.
  logical :: bound_Coriolis  !< If true & SMAGORINSKY_AH is used, the biharmonic
                             !! viscosity is modified to include a term that
                             !! scales quadratically with the velocity shears.
  real    :: Kh_bg_min       !< The minimum value allowed for Laplacian horizontal
                             !! viscosity [L2 T-1 ~> m2 s-1]. The default is 0.0.
  logical :: FrictWork_bug   !< If true, retain an answer-changing bug in calculating FrictWork,
                             !! which cancels the h in thickness flux and the h at velocity point.
  logical :: OBC_strain_bug  !< If true, recover a bug that specified shear strain option at open
                             !! boundaries cannot be applied.
  logical :: use_land_mask   !< Use the land mask for the computation of thicknesses
                             !! at velocity locations. This eliminates the dependence on
                             !! arbitrary values over land or outside of the domain.
                             !! Default is False to maintain answers with legacy experiments
                             !! but should be changed to True for new experiments.
  logical :: anisotropic     !< If true, allow anisotropic component to the viscosity.
  logical :: add_LES_viscosity!< If true, adds the viscosity from Smagorinsky and Leith to
                             !! the background viscosity instead of taking the maximum.
  real    :: Kh_aniso        !< The anisotropic viscosity [L2 T-1 ~> m2 s-1].
  logical :: res_scale_MEKE  !< If true, the viscosity contribution from MEKE is scaled by
                             !! the resolution function.
  logical :: use_GME         !< If true, use GME backscatter scheme.
  integer :: answer_date     !< The vintage of the order of arithmetic and expressions in the
                             !! horizontal viscosity calculations.  Values below 20190101 recover
                             !! the answers from the end of 2018, while higher values use updated
                             !! and more robust forms of the same expressions.
  real    :: GME_h0          !< The strength of GME tapers quadratically to zero when the bathymetric
                             !! total water column thickness is less than GME_H0 [H ~> m or kg m-2]
  real    :: GME_efficiency  !< The nondimensional prefactor multiplying the GME coefficient [nondim]
  real    :: GME_limiter     !< The absolute maximum value the GME coefficient is allowed to take [L2 T-1 ~> m2 s-1].
  real    :: min_grid_Kh     !< Minimum horizontal Laplacian viscosity used to
                             !! limit the grid Reynolds number [L2 T-1 ~> m2 s-1]
  real    :: min_grid_Ah     !< Minimun horizontal biharmonic viscosity used to
                             !! limit grid Reynolds number [L4 T-1 ~> m4 s-1]
  logical :: use_cont_thick  !< If true, thickness at velocity points adopts h[uv] in BT_cont from continuity solver.
  logical :: use_ZB2020      !< If true, use Zanna-Bolton 2020 parameterization.
  logical :: use_circulation !< If true, use circulation theorem to compute vorticity (for ZB20 or Leith)
  logical :: do_grid_Re_Ah = .false. !< If true, the grid_Re_Ah diagnostic is requested.
  logical :: do_grid_Re_Kh = .false. !< If true, the grid_Re_Kh diagnostic is requested.
  logical :: do_Ah_h = .false. !< If true, the Ah_h diagnostic is requested.
  logical :: do_Ah_q = .false. !< If true, the Ah_q diagnostic is requested.
  logical :: do_Kh_h = .false. !< If true, the Kh_h diagnostic is requested.
  logical :: do_Kh_q = .false. !< If true, the Kh_q diagnostic is requested.
  logical :: do_GME_coeff_h = .false. !< If true, the GME_coeff_h diagnostic is requested.
  logical :: do_GME_coeff_q = .false. !< If true, the GME_coeff_q diagnostic is requested.
  logical :: do_vort_xy_q = .false. !< If true, the vort_xy_q diagnostic is requested.
  logical :: do_div_xx_h = .false. !< If true, the div_xx_h diagnostic is requested.
  logical :: do_sh_xy_q = .false. !< If true, the sh_xy_q diagnostic is requested.
  logical :: do_sh_xx_h = .false. !< If true, the sh_xx_h diagnostic is requested.
  logical :: do_FrictWork = .false. !< If true, the FrictWork diagnostic is requested.
  logical :: do_FrictWorkIntz = .false. !< If true, the FrictWorkIntz diagnostic is requested.
  logical :: do_FrictWork_bh = .false. !< If true, the FrictWork_bh diagnostic is requested.
  logical :: do_FrictWorkIntz_bh = .false. !< If true, the FrictWorkIntz_bh diagnostic is requested.
  logical :: do_normstress = .false. !< If true, the normstress diagnostic is requested.
  logical :: do_shearstress = .false. !< If true, the shearstress diagnostic is requested.
  logical :: do_BS_coeff_h = .false. !< If true, the BS_coeff_h diagnostic is requested.
  logical :: do_BS_coeff_q = .false. !< If true, the BS_coeff_q diagnostic is requested.
end type hv_opts_CS

!> bind(C) mirror of hv_opts_CS, field-for-field, same order.
type, bind(C) :: hv_opts_CS_C
  logical(c_bool) :: Laplacian                 !< Use a Laplacian horizontal viscosity if true.
  logical(c_bool) :: biharmonic                !< Use a biharmonic horizontal viscosity if true.
  logical(c_bool) :: debug                     !< If true, write verbose checksums for debugging purposes.
  logical(c_bool) :: no_slip                   !< If true, no slip boundary conditions are used. Otherwise free slip
                                            !! boundary conditions are assumed. The implementation of the free slip
                                            !! boundary conditions on a C-grid is much cleaner than the no slip boundary
                                            !! conditions. The use of free slip b.c.s is strongly encouraged. The no
                                            !! slip b.c.s are not implemented with the biharmonic viscosity.
  logical(c_bool) :: bound_Kh                  !< If true, the Laplacian coefficient is locally limited to guarantee
                                            !! stability.
  logical(c_bool) :: EY24_EBT_BS               !< If true, use an equivalent barotropic backscatter with a stabilizing
                                            !! kill switch in MEKE, developed by Yankovsky et al. 2024
  logical(c_bool) :: bound_Ah                  !< If true, the biharmonic coefficient is locally limited to guarantee
                                            !! stability.
  real(c_double)  :: Re_Ah                     !< If nonzero, the biharmonic coefficient is scaled so that the
                                            !! biharmonic Reynolds number is equal to this [nondim].
  real(c_double)  :: KS_coef                   !< A nondimensional coefficient on the biharmonic viscosity that sets the
                                            !! kill switch for backscatter. Default is 1.0 [nondim].
  logical(c_bool) :: backscatter_underbound    !< If true, the bounds on the biharmonic viscosity are allowed to
                                            !! increase where the Laplacian viscosity is negative (due to backscatter
                                            !! parameterizations) beyond the largest timestep-dependent stable values of
                                            !! biharmonic viscosity when no Laplacian viscosity is applied.  The default
                                            !! is true for historical reasons, but this option probably should not be
                                            !! used as it can lead to numerical instabilities.
  logical(c_bool) :: Smagorinsky_Kh            !< If true, use Smagorinsky nonlinear eddy viscosity. KH is the
                                            !! background value.
  logical(c_bool) :: Smagorinsky_Ah            !< If true, use a biharmonic form of Smagorinsky nonlinear eddy
                                            !! viscosity. AH is the background.
  logical(c_bool) :: Leith_Kh                  !< If true, use 2D Leith nonlinear eddy viscosity. KH is the background
                                            !! value.
  logical(c_bool) :: Modified_Leith            !< If true, use extra component of Leith viscosity to damp divergent
                                            !! flow. To use, still set Leith_Kh=.TRUE.
  logical(c_bool) :: use_beta_in_Leith         !< If true, includes the beta term in the Leith viscosity
  logical(c_bool) :: Leith_Ah                  !< If true, use a biharmonic form of 2D Leith nonlinear eddy viscosity.
                                            !! AH is the background.
  logical(c_bool) :: use_Leithy                !< If true, use a biharmonic form of 2D Leith nonlinear eddy viscosity
                                            !! with harmonic backscatter. Ah is the background. Leithy = Leith+E
  real(c_double)  :: c_K                       !< Fraction of energy dissipated by the biharmonic term that gets
                                            !! backscattered in the Leith+E scheme. [nondim]
  logical(c_bool) :: smooth_Ah                 !< If true (default), then Ah and m_leithy are smoothed. This smoothing
                                            !! requires a lot of blocking communication.
  logical(c_bool) :: use_QG_Leith_visc         !< If true, use QG Leith nonlinear eddy viscosity. KH is the background
                                            !! value.
  logical(c_bool) :: bound_Coriolis            !< If true & SMAGORINSKY_AH is used, the biharmonic viscosity is modified
                                            !! to include a term that scales quadratically with the velocity shears.
  real(c_double)  :: Kh_bg_min                 !< The minimum value allowed for Laplacian horizontal viscosity
                                            !! [L2 T-1 ~> m2 s-1]. The default is 0.0.
  logical(c_bool) :: FrictWork_bug             !< If true, retain an answer-changing bug in calculating FrictWork, which
                                            !! cancels the h in thickness flux and the h at velocity point.
  logical(c_bool) :: OBC_strain_bug            !< If true, recover a bug that specified shear strain option at open
                                            !! boundaries cannot be applied.
  logical(c_bool) :: use_land_mask             !< Use the land mask for the computation of thicknesses at velocity
                                            !! locations. This eliminates the dependence on arbitrary values over land
                                            !! or outside of the domain. Default is False to maintain answers with
                                            !! legacy experiments but should be changed to True for new experiments.
  logical(c_bool) :: anisotropic               !< If true, allow anisotropic component to the viscosity.
  logical(c_bool) :: add_LES_viscosity         !< If true, adds the viscosity from Smagorinsky and Leith to the
                                            !! background viscosity instead of taking the maximum.
  real(c_double)  :: Kh_aniso                  !< The anisotropic viscosity [L2 T-1 ~> m2 s-1].
  logical(c_bool) :: res_scale_MEKE            !< If true, the viscosity contribution from MEKE is scaled by the
                                            !! resolution function.
  logical(c_bool) :: use_GME                   !< If true, use GME backscatter scheme.
  integer(c_int)  :: answer_date               !< The vintage of the order of arithmetic and expressions in the
                                            !! horizontal viscosity calculations.  Values below 20190101 recover the
                                            !! answers from the end of 2018, while higher values use updated and more
                                            !! robust forms of the same expressions.
  real(c_double)  :: GME_h0                    !< The strength of GME tapers quadratically to zero when the bathymetric
                                            !! total water column thickness is less than GME_H0 [H ~> m or kg m-2]
  real(c_double)  :: GME_efficiency            !< The nondimensional prefactor multiplying the GME coefficient [nondim]
  real(c_double)  :: GME_limiter               !< The absolute maximum value the GME coefficient is allowed to take
                                            !! [L2 T-1 ~> m2 s-1].
  real(c_double)  :: min_grid_Kh               !< Minimum horizontal Laplacian viscosity used to limit the grid Reynolds
                                            !! number [L2 T-1 ~> m2 s-1]
  real(c_double)  :: min_grid_Ah               !< Minimun horizontal biharmonic viscosity used to limit grid Reynolds
                                            !! number [L4 T-1 ~> m4 s-1]
  logical(c_bool) :: use_cont_thick            !< If true, thickness at velocity points adopts h[uv] in BT_cont from
                                            !! continuity solver.
  logical(c_bool) :: use_ZB2020                !< If true, use Zanna-Bolton 2020 parameterization.
  logical(c_bool) :: use_circulation           !< If true, use circulation theorem to compute vorticity (for ZB20 or
                                            !! Leith)
  logical(c_bool) :: do_grid_Re_Ah             !< If true, the grid_Re_Ah diagnostic is requested.
  logical(c_bool) :: do_grid_Re_Kh             !< If true, the grid_Re_Kh diagnostic is requested.
  logical(c_bool) :: do_Ah_h                   !< If true, the Ah_h diagnostic is requested.
  logical(c_bool) :: do_Ah_q                   !< If true, the Ah_q diagnostic is requested.
  logical(c_bool) :: do_Kh_h                   !< If true, the Kh_h diagnostic is requested.
  logical(c_bool) :: do_Kh_q                   !< If true, the Kh_q diagnostic is requested.
  logical(c_bool) :: do_GME_coeff_h            !< If true, the GME_coeff_h diagnostic is requested.
  logical(c_bool) :: do_GME_coeff_q            !< If true, the GME_coeff_q diagnostic is requested.
  logical(c_bool) :: do_vort_xy_q              !< If true, the vort_xy_q diagnostic is requested.
  logical(c_bool) :: do_div_xx_h               !< If true, the div_xx_h diagnostic is requested.
  logical(c_bool) :: do_sh_xy_q                !< If true, the sh_xy_q diagnostic is requested.
  logical(c_bool) :: do_sh_xx_h                !< If true, the sh_xx_h diagnostic is requested.
  logical(c_bool) :: do_FrictWork              !< If true, the FrictWork diagnostic is requested.
  logical(c_bool) :: do_FrictWorkIntz          !< If true, the FrictWorkIntz diagnostic is requested.
  logical(c_bool) :: do_FrictWork_bh           !< If true, the FrictWork_bh diagnostic is requested.
  logical(c_bool) :: do_FrictWorkIntz_bh       !< If true, the FrictWorkIntz_bh diagnostic is requested.
  logical(c_bool) :: do_normstress             !< If true, the normstress diagnostic is requested.
  logical(c_bool) :: do_shearstress            !< If true, the shearstress diagnostic is requested.
  logical(c_bool) :: do_BS_coeff_h             !< If true, the BS_coeff_h diagnostic is requested.
  logical(c_bool) :: do_BS_coeff_q             !< If true, the BS_coeff_q diagnostic is requested.
end type hv_opts_CS_C

!> Persistent copies of the horizontal viscosity metric and coefficient arrays, set once
!! when the horizontal viscosity is initialized.
type :: hv_metrics_CS
  type(RealArray_t) :: Kh_bg_xx         !< The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1].
                                        !! The actual viscosity may be the larger of this
                                        !! viscosity and the Smagorinsky and Leith viscosities.
  type(RealArray_t) :: Ah_bg_xx         !< The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1].
                                        !! The actual viscosity may be the larger of this
                                        !! viscosity and the Smagorinsky and Leith viscosities.
  type(RealArray_t) :: reduction_xx     !< The amount by which stresses through h points are reduced
                                        !! due to partial barriers [nondim].
  type(RealArray_t) :: Kh_Max_xx        !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
  type(RealArray_t) :: Ah_Max_xx        !< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
  type(RealArray_t) :: Ah_Max_xx_KS     !< The maximum permitted biharmonic viscosity for
                                        !! the kill switch [L4 T-1 ~> m4 s-1].
  type(RealArray_t) :: n1n2_h           !< Factor n1*n2 in the anisotropic direction tensor at h-points [nondim]
  type(RealArray_t) :: n1n1_m_n2n2_h    !< Factor n1**2-n2**2 in the anisotropic direction tensor at h-points [nondim]
  type(RealArray_t) :: grid_sp_h2       !< Harmonic mean of the squares of the grid [L2 ~> m2]
  type(RealArray_t) :: grid_sp_h3       !< Harmonic mean of the squares of the grid^(3/2) [L3 ~> m3]
  type(RealArray_t) :: Kh_bg_xy         !< The background Laplacian viscosity at q points [L2 T-1 ~> m2 s-1].
                                        !! The actual viscosity may be the larger of this
                                        !! viscosity and the Smagorinsky and Leith viscosities.
  type(RealArray_t) :: Ah_bg_xy         !< The background biharmonic viscosity at q points [L4 T-1 ~> m4 s-1].
                                        !! The actual viscosity may be the larger of this
                                        !! viscosity and the Smagorinsky and Leith viscosities.
  type(RealArray_t) :: reduction_xy     !< The amount by which stresses through q points are reduced
                                        !! due to partial barriers [nondim].
  type(RealArray_t) :: Kh_Max_xy        !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
  type(RealArray_t) :: Ah_Max_xy        !< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
  type(RealArray_t) :: Ah_Max_xy_KS     !< The maximum permitted biharmonic viscosity for
                                        !! the  kill switch [L4 T-1 ~> m4 s-1].
  type(RealArray_t) :: n1n2_q           !< Factor n1*n2 in the anisotropic direction tensor at q-points [nondim]
  type(RealArray_t) :: n1n1_m_n2n2_q    !< Factor n1**2-n2**2 in the anisotropic direction tensor at q-points [nondim]
  type(RealArray_t) :: dx2h             !< Pre-calculated dx^2 at h points [L2 ~> m2]
  type(RealArray_t) :: dy2h             !< Pre-calculated dy^2 at h points [L2 ~> m2]
  type(RealArray_t) :: dx_dyT           !< Pre-calculated dx/dy at h points [nondim]
  type(RealArray_t) :: dy_dxT           !< Pre-calculated dy/dx at h points [nondim]
  type(RealArray_t) :: m_const_leithy   !< Pre-calculated .5*sqrt(c_K)*max{dx,dy} [L ~> m]
  type(RealArray_t) :: m_leithy_max     !< Pre-calculated 4./max(dx,dy)^2 at h points [L-2 ~> m-2]
  type(RealArray_t) :: dx2q             !< Pre-calculated dx^2 at q points [L2 ~> m2]
  type(RealArray_t) :: dy2q             !< Pre-calculated dy^2 at q points [L2 ~> m2]
  type(RealArray_t) :: dx_dyBu          !< Pre-calculated dx/dy at q points [nondim]
  type(RealArray_t) :: dy_dxBu          !< Pre-calculated dy/dx at q points [nondim]
  type(RealArray_t) :: Idx2dyCu         !< 1/(dx^2 dy) at u points [L-3 ~> m-3]
  type(RealArray_t) :: Idxdy2u          !< 1/(dx dy^2) at u points [L-3 ~> m-3]
  type(RealArray_t) :: Idx2dyCv         !< 1/(dx^2 dy) at v points [L-3 ~> m-3]
  type(RealArray_t) :: Idxdy2v          !< 1/(dx dy^2) at v points [L-3 ~> m-3]
  type(RealArray_t) :: Laplac2_const_xx !< Laplacian metric-dependent constants [L2 ~> m2]
  type(RealArray_t) :: Biharm6_const_xx !< Biharmonic metric-dependent constants [L6 ~> m6]
  type(RealArray_t) :: Laplac3_const_xx !< Laplacian metric-dependent constants [L3 ~> m3]
  type(RealArray_t) :: Biharm_const_xx  !< Biharmonic metric-dependent constants [L4 ~> m4]
  type(RealArray_t) :: Biharm_const2_xx !< Biharmonic metric-dependent constants [T L4 ~> s m4]
  type(RealArray_t) :: Re_Ah_const_xx   !< Biharmonic metric-dependent constants [L3 ~> m3]
  type(RealArray_t) :: Laplac2_const_xy !< Laplacian metric-dependent constants [L2 ~> m2]
  type(RealArray_t) :: Biharm6_const_xy !< Biharmonic metric-dependent constants [L6 ~> m6]
  type(RealArray_t) :: Laplac3_const_xy !< Laplacian metric-dependent constants [L3 ~> m3]
  type(RealArray_t) :: Biharm_const_xy  !< Biharmonic metric-dependent constants [L4 ~> m4]
  type(RealArray_t) :: Biharm_const2_xy !< Biharmonic metric-dependent constants [T L4 ~> s m4]
  type(RealArray_t) :: Re_Ah_const_xy   !< Biharmonic metric-dependent constants [L3 ~> m3]
end type hv_metrics_CS

!> bind(C) mirror of hv_metrics_CS, field-for-field, same order.
type, bind(C) :: hv_metrics_CS_C
  type(RealArray_C) :: Kh_bg_xx         !< The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1]. The actual
                                     !! viscosity may be the larger of this viscosity and the Smagorinsky and Leith
                                     !! viscosities.
  type(RealArray_C) :: Ah_bg_xx         !< The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1]. The
                                     !! actual viscosity may be the larger of this viscosity and the Smagorinsky and
                                     !! Leith viscosities.
  type(RealArray_C) :: reduction_xx     !< The amount by which stresses through h points are reduced due to partial
                                     !! barriers [nondim].
  type(RealArray_C) :: Kh_Max_xx        !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
  type(RealArray_C) :: Ah_Max_xx        !< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
  type(RealArray_C) :: Ah_Max_xx_KS     !< The maximum permitted biharmonic viscosity for the kill switch
                                     !! [L4 T-1 ~> m4 s-1].
  type(RealArray_C) :: n1n2_h           !< Factor n1*n2 in the anisotropic direction tensor at h-points [nondim]
  type(RealArray_C) :: n1n1_m_n2n2_h    !< Factor n1**2-n2**2 in the anisotropic direction tensor at h-points [nondim]
  type(RealArray_C) :: grid_sp_h2       !< Harmonic mean of the squares of the grid [L2 ~> m2]
  type(RealArray_C) :: grid_sp_h3       !< Harmonic mean of the squares of the grid^(3/2) [L3 ~> m3]
  type(RealArray_C) :: Kh_bg_xy         !< The background Laplacian viscosity at q points [L2 T-1 ~> m2 s-1]. The actual
                                     !! viscosity may be the larger of this viscosity and the Smagorinsky and Leith
                                     !! viscosities.
  type(RealArray_C) :: Ah_bg_xy         !< The background biharmonic viscosity at q points [L4 T-1 ~> m4 s-1]. The
                                     !! actual viscosity may be the larger of this viscosity and the Smagorinsky and
                                     !! Leith viscosities.
  type(RealArray_C) :: reduction_xy     !< The amount by which stresses through q points are reduced due to partial
                                     !! barriers [nondim].
  type(RealArray_C) :: Kh_Max_xy        !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
  type(RealArray_C) :: Ah_Max_xy        !< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
  type(RealArray_C) :: Ah_Max_xy_KS     !< The maximum permitted biharmonic viscosity for the  kill switch
                                     !! [L4 T-1 ~> m4 s-1].
  type(RealArray_C) :: n1n2_q           !< Factor n1*n2 in the anisotropic direction tensor at q-points [nondim]
  type(RealArray_C) :: n1n1_m_n2n2_q    !< Factor n1**2-n2**2 in the anisotropic direction tensor at q-points [nondim]
  type(RealArray_C) :: dx2h             !< Pre-calculated dx^2 at h points [L2 ~> m2]
  type(RealArray_C) :: dy2h             !< Pre-calculated dy^2 at h points [L2 ~> m2]
  type(RealArray_C) :: dx_dyT           !< Pre-calculated dx/dy at h points [nondim]
  type(RealArray_C) :: dy_dxT           !< Pre-calculated dy/dx at h points [nondim]
  type(RealArray_C) :: m_const_leithy   !< Pre-calculated .5*sqrt(c_K)*max{dx,dy} [L ~> m]
  type(RealArray_C) :: m_leithy_max     !< Pre-calculated 4./max(dx,dy)^2 at h points [L-2 ~> m-2]
  type(RealArray_C) :: dx2q             !< Pre-calculated dx^2 at q points [L2 ~> m2]
  type(RealArray_C) :: dy2q             !< Pre-calculated dy^2 at q points [L2 ~> m2]
  type(RealArray_C) :: dx_dyBu          !< Pre-calculated dx/dy at q points [nondim]
  type(RealArray_C) :: dy_dxBu          !< Pre-calculated dy/dx at q points [nondim]
  type(RealArray_C) :: Idx2dyCu         !< 1/(dx^2 dy) at u points [L-3 ~> m-3]
  type(RealArray_C) :: Idxdy2u          !< 1/(dx dy^2) at u points [L-3 ~> m-3]
  type(RealArray_C) :: Idx2dyCv         !< 1/(dx^2 dy) at v points [L-3 ~> m-3]
  type(RealArray_C) :: Idxdy2v          !< 1/(dx dy^2) at v points [L-3 ~> m-3]
  type(RealArray_C) :: Laplac2_const_xx !< Laplacian metric-dependent constants [L2 ~> m2]
  type(RealArray_C) :: Biharm6_const_xx !< Biharmonic metric-dependent constants [L6 ~> m6]
  type(RealArray_C) :: Laplac3_const_xx !< Laplacian metric-dependent constants [L3 ~> m3]
  type(RealArray_C) :: Biharm_const_xx  !< Biharmonic metric-dependent constants [L4 ~> m4]
  type(RealArray_C) :: Biharm_const2_xx !< Biharmonic metric-dependent constants [T L4 ~> s m4]
  type(RealArray_C) :: Re_Ah_const_xx   !< Biharmonic metric-dependent constants [L3 ~> m3]
  type(RealArray_C) :: Laplac2_const_xy !< Laplacian metric-dependent constants [L2 ~> m2]
  type(RealArray_C) :: Biharm6_const_xy !< Biharmonic metric-dependent constants [L6 ~> m6]
  type(RealArray_C) :: Laplac3_const_xy !< Laplacian metric-dependent constants [L3 ~> m3]
  type(RealArray_C) :: Biharm_const_xy  !< Biharmonic metric-dependent constants [L4 ~> m4]
  type(RealArray_C) :: Biharm_const2_xy !< Biharmonic metric-dependent constants [T L4 ~> s m4]
  type(RealArray_C) :: Re_Ah_const_xy   !< Biharmonic metric-dependent constants [L3 ~> m3]
end type hv_metrics_CS_C

!> Control structure for horizontal viscosity
type, public :: hor_visc_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: bound_coef      !< The nondimensional coefficient of the ratio of
                             !! the viscosity bounds to the theoretical maximum
                             !! for stability without considering other terms [nondim].
                             !! The default is 0.8.
  real    :: KS_timescale    !< A timescale for computing CFL limit for turning off backscatter [T ~> s].
  logical :: use_Kh_bg_2d    !< Read 2d background viscosity from a file.
  logical :: Kh_bg_2d_bug    !< If true, retain an answer-changing horizontal indexing bug
                             !! in setting the corner-point viscosities when USE_KH_BG_2D=True.
  logical :: dynamic_aniso   !< If true, the anisotropic viscosity is recomputed as a function
                             !! of state. This is set depending on ANISOTROPIC_MODE.
  integer :: nkblock         !< The k block size used in horizontal viscosity calculations [nondim].
  logical :: use_cont_thick_bug  !< If true, retain an answer-changing bug for thickness at velocity points.
  type(ZB2020_CS) :: ZB2020  !< Zanna-Bolton 2020 control structure.

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: Kh_bg_xx
                      !< The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real, allocatable :: Kh_bg_2d(:,:)
                      !< The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: Ah_bg_xx
                      !< The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: reduction_xx
                      !< The amount by which stresses through h points are reduced
                      !! due to partial barriers [nondim].
  real, allocatable :: Kh_Max_xx(:,:)     !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
  real, allocatable :: Ah_Max_xx(:,:)     !< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
  real, allocatable :: Ah_Max_xx_KS(:,:)  !< The maximum permitted biharmonic viscosity for
                                          !! the kill switch [L4 T-1 ~> m4 s-1].
  real, allocatable :: n1n2_h(:,:)        !< Factor n1*n2 in the anisotropic direction tensor at h-points [nondim]
  real, allocatable :: n1n1_m_n2n2_h(:,:) !< Factor n1**2-n2**2 in the anisotropic direction tensor at h-points [nondim]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    grid_sp_h2,     & !< Harmonic mean of the squares of the grid [L2 ~> m2]
    grid_sp_h3        !< Harmonic mean of the squares of the grid^(3/2) [L3 ~> m3]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: Kh_bg_xy
                      !< The background Laplacian viscosity at q points [L2 T-1 ~> m2 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: Ah_bg_xy
                      !< The background biharmonic viscosity at q points [L4 T-1 ~> m4 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: reduction_xy
                      !< The amount by which stresses through q points are reduced
                      !! due to partial barriers [nondim].
  real, allocatable :: Kh_Max_xy(:,:)  !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
  real, allocatable :: Ah_Max_xy(:,:)     !< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
  real, allocatable :: Ah_Max_xy_KS(:,:)  !< The maximum permitted biharmonic viscosity for
                                          !! the  kill switch [L4 T-1 ~> m4 s-1].
  real, allocatable :: n1n2_q(:,:)        !< Factor n1*n2 in the anisotropic direction tensor at q-points [nondim]
  real, allocatable :: n1n1_m_n2n2_q(:,:) !< Factor n1**2-n2**2 in the anisotropic direction tensor at q-points [nondim]

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    dx2h,           & !< Pre-calculated dx^2 at h points [L2 ~> m2]
    dy2h,           & !< Pre-calculated dy^2 at h points [L2 ~> m2]
    dx_dyT,         & !< Pre-calculated dx/dy at h points [nondim]
    dy_dxT            !< Pre-calculated dy/dx at h points [nondim]
  real, allocatable :: m_const_leithy(:,:) !< Pre-calculated .5*sqrt(c_K)*max{dx,dy} [L ~> m]
  real, allocatable :: m_leithy_max(:,:)   !< Pre-calculated 4./max(dx,dy)^2 at h points [L-2 ~> m-2]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    dx2q,    & !< Pre-calculated dx^2 at q points [L2 ~> m2]
    dy2q,    & !< Pre-calculated dy^2 at q points [L2 ~> m2]
    dx_dyBu, & !< Pre-calculated dx/dy at q points [nondim]
    dy_dxBu    !< Pre-calculated dy/dx at q points [nondim]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    Idx2dyCu, & !< 1/(dx^2 dy) at u points [L-3 ~> m-3]
    Idxdy2u     !< 1/(dx dy^2) at u points [L-3 ~> m-3]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    Idx2dyCv, & !< 1/(dx^2 dy) at v points [L-3 ~> m-3]
    Idxdy2v     !< 1/(dx dy^2) at v points [L-3 ~> m-3]

  ! The following variables are precalculated time-invariant combinations of
  ! parameters and metric terms.
  real, allocatable :: Laplac2_const_xx(:,:) !< Laplacian metric-dependent constants [L2 ~> m2]
  real, allocatable :: Biharm6_const_xx(:,:) !< Biharmonic metric-dependent constants [L6 ~> m6]
  real, allocatable :: Laplac3_const_xx(:,:) !< Laplacian metric-dependent constants [L3 ~> m3]
  real, allocatable :: Biharm_const_xx(:,:)  !< Biharmonic metric-dependent constants [L4 ~> m4]
  real, allocatable :: Biharm_const2_xx(:,:) !< Biharmonic metric-dependent constants [T L4 ~> s m4]
  real, allocatable :: Re_Ah_const_xx(:,:)   !< Biharmonic metric-dependent constants [L3 ~> m3]

  real, allocatable :: Laplac2_const_xy(:,:) !< Laplacian metric-dependent constants [L2 ~> m2]
  real, allocatable :: Biharm6_const_xy(:,:) !< Biharmonic metric-dependent constants [L6 ~> m6]
  real, allocatable :: Laplac3_const_xy(:,:) !< Laplacian metric-dependent constants [L3 ~> m3]
  real, allocatable :: Biharm_const_xy(:,:)  !< Biharmonic metric-dependent constants [L4 ~> m4]
  real, allocatable :: Biharm_const2_xy(:,:) !< Biharmonic metric-dependent constants [T L4 ~> s m4]
  real, allocatable :: Re_Ah_const_xy(:,:)   !< Biharmonic metric-dependent constants [L3 ~> m3]

  type(diag_ctrl), pointer :: diag => NULL() !< structure to regulate diagnostics

  ! real, allocatable :: hf_diffu(:,:,:)  ! Zonal horizontal viscous acceleleration times
  !                                       ! fractional thickness [L T-2 ~> m s-2].
  ! real, allocatable :: hf_diffv(:,:,:)  ! Meridional horizontal viscous acceleleration times
  !                                       ! fractional thickness [L T-2 ~> m s-2].
  ! 3D diagnostics hf_diffu(diffv) are commented because there is no clarity on proper remapping grid option.
  ! The code is retained for debugging purposes in the future.

  integer :: num_smooth_gme !< number of smoothing passes for the GME fluxes.
  !>@{
  !! Diagnostic id
  integer :: id_grid_Re_Ah = -1, id_grid_Re_Kh   = -1
  integer :: id_diffu     = -1, id_diffv         = -1
  ! integer :: id_hf_diffu  = -1, id_hf_diffv      = -1
  integer :: id_h_diffu  = -1, id_h_diffv      = -1
  integer :: id_hf_diffu_2d = -1, id_hf_diffv_2d = -1
  integer :: id_intz_diffu_2d = -1, id_intz_diffv_2d = -1
  integer :: id_diffu_visc_rem = -1, id_diffv_visc_rem = -1
  integer :: id_Ah_h      = -1, id_Ah_q          = -1
  integer :: id_Kh_h      = -1, id_Kh_q          = -1
  integer :: id_GME_coeff_h = -1, id_GME_coeff_q = -1
  integer :: id_dudx_bt = -1, id_dvdy_bt = -1
  integer :: id_dudy_bt = -1, id_dvdx_bt = -1
  integer :: id_vort_xy_q = -1, id_div_xx_h      = -1
  integer :: id_sh_xy_q = -1,    id_sh_xx_h      = -1
  integer :: id_FrictWork = -1, id_FrictWorkIntz = -1
  integer :: id_FrictWork_bh = -1, id_FrictWorkIntz_bh = -1
  integer :: id_FrictWork_GME = -1
  integer :: id_normstress = -1, id_shearstress = -1
  integer :: id_visc_limit_h = -1, id_visc_limit_q = -1
  integer :: id_visc_limit_h_flag = -1, id_visc_limit_q_flag = -1
  integer :: id_visc_limit_h_frac = -1, id_visc_limit_q_frac = -1
  integer :: id_BS_coeff_h = -1, id_BS_coeff_q = -1
  !>@}

  type(hv_opts_CS)    :: hv_opts_CS    !< Options and diagnostic switches used by the horizontal
                                       !! viscosity kernels
  type(hv_metrics_CS) :: hv_metrics_CS !< Persistent copies of the metric and coefficient arrays
                                       !! used by the horizontal viscosity kernels
end type hor_visc_CS

!> The Fortran objects that the horizontal_viscosity host helpers need. horizontal_viscosity
!! builds one and passes it down the tree as an opaque c_ptr.
type :: hv_host_ctx
  type(ocean_grid_type),   pointer :: G => NULL()      !< The ocean's grid structure.
  type(verticalGrid_type), pointer :: GV => NULL()     !< The ocean's vertical grid structure.
  type(unit_scale_type),   pointer :: US => NULL()     !< A dimensional unit scaling type
  type(VarMix_CS),         pointer :: VarMix => NULL() !< Variable mixing control structure
  type(hor_visc_CS),       pointer :: CS => NULL()     !< Horizontal viscosity control structure
end type hv_host_ctx

!> A container-based shadow of the MEKE_type fields used within the horizontal_viscosity call
!! tree, so that they can be bridge-ready without converting the shared MEKE_type itself.
type :: hor_visc_MEKE_type
  type(RealArray_t) :: Ku_a         !< The MEKE-derived lateral viscosity coefficient
                                    !! [L2 T-1 ~> m2 s-1]. This viscosity can be negative when representing
                                    !! backscatter from unresolved eddies (see Jansen and Held, 2014).
  type(RealArray_t) :: Au_a         !< The MEKE-derived lateral biharmonic viscosity
                                    !! coefficient [L4 T-1 ~> m4 s-1].
  type(RealArray_t) :: mom_src_a    !< MEKE source from lateral friction in the
                                    !! momentum equations [R Z L2 T-3 ~> W m-2].
  type(RealArray_t) :: mom_src_bh_a !< MEKE source from the biharmonic part of the lateral friction in the
                                    !! momentum equations [R Z L2 T-3 ~> W m-2].
  type(RealArray_t) :: GME_snk_a    !< MEKE sink from GME backscatter in the momentum equations
                                    !! [R Z L2 T-3 ~> W m-2].
  real :: backscatter_Ro_pow = 0.0 !< Power in Rossby number function for backscatter [nondim].
  real :: backscatter_Ro_c = 0.0 !< Coefficient in Rossby number function for backscatter [nondim].
contains
  procedure :: build_from => hor_visc_MEKE_build_from !< Copy in the allocated fields
  procedure :: copy_back  => hor_visc_MEKE_copy_back  !< Copy out the source and sink fields and free
end type hor_visc_MEKE_type

!> A container-based shadow of the VarMix_CS fields used within the horizontal_viscosity call
!! tree, so that they can be bridge-ready without converting the shared VarMix_CS itself.
type :: hor_visc_VarMix_type
  type(RealArray_t) :: Res_fn_h_a  !< Non-dimensional function of the ratio the first baroclinic
                                   !! deformation radius to the grid spacing at h points [nondim].
  type(RealArray_t) :: Res_fn_q_a  !< Non-dimensional function of the ratio the first baroclinic
                                   !! deformation radius to the grid spacing at q points [nondim].
  type(RealArray_t) :: BS_struct_a !< Vertical structure function used in backscatter [nondim]
  logical :: use_variable_mixing = .false. !< If true, use the variable mixing.
  logical :: Resoln_scaled_Kh = .false.    !< If true, scale away the Laplacian viscosity
                                           !! when the deformation radius is well resolved.
contains
  procedure :: build_from => hor_visc_VarMix_build_from !< Copy in the allocated fields
  procedure :: free       => hor_visc_VarMix_free       !< Free every allocated field
end type hor_visc_VarMix_type

!> A container-based shadow of the stochastic_CS fields used within the horizontal_viscosity call
!! tree, so that they can be bridge-ready without converting the shared stochastic_CS itself.
type :: hor_visc_STOCH_type
  type(RealArray_t) :: skeb_diss_a !< Dissipation rate used to set amplitude of SKEBS [L2 T-3 ~> m2 s-3]
                                   !! Index into this at h points.
  logical :: skeb_use_frict = .false. !< If true, adds viscous dissipation rate to the amplitude of SKEBS
  real    :: skeb_frict_coef = 0.0 !< If skeb_use_frict is true, then skeb_gm_coef * GM_work is added to the
                                   !! dissipation rate used to set the amplitude of SKEBS [nondim]
contains
  procedure :: build_from => hor_visc_STOCH_build_from !< Copy in the allocated fields
  procedure :: copy_back  => hor_visc_STOCH_copy_back  !< Copy out the dissipation rate and free
end type hor_visc_STOCH_type

!> bind(C) mirror of hor_visc_MEKE_type, one RealArray_C per container (same order as
!! hor_visc_MEKE_type itself), followed by its scalars.
type, bind(C) :: hor_visc_MEKE_C
  type(RealArray_C) :: Ku                 !< The MEKE-derived lateral viscosity coefficient
                                          !! [L2 T-1 ~> m2 s-1]. This viscosity can be negative when representing
                                          !! backscatter from unresolved eddies (see Jansen and Held, 2014).
  type(RealArray_C) :: Au                 !< The MEKE-derived lateral biharmonic viscosity
                                          !! coefficient [L4 T-1 ~> m4 s-1].
  type(RealArray_C) :: mom_src            !< MEKE source from lateral friction in the
                                          !! momentum equations [R Z L2 T-3 ~> W m-2].
  type(RealArray_C) :: mom_src_bh         !< MEKE source from the biharmonic part of the lateral friction in the
                                          !! momentum equations [R Z L2 T-3 ~> W m-2].
  type(RealArray_C) :: GME_snk            !< MEKE sink from GME backscatter in the momentum equations
                                          !! [R Z L2 T-3 ~> W m-2].
  real(c_double)    :: backscatter_Ro_pow !< Power in Rossby number function for backscatter [nondim].
  real(c_double)    :: backscatter_Ro_c   !< Coefficient in Rossby number function for backscatter [nondim].
end type hor_visc_MEKE_C

!> bind(C) mirror of hor_visc_VarMix_type, one RealArray_C per container (same order as
!! hor_visc_VarMix_type itself), followed by its scalars.
type, bind(C) :: hor_visc_VarMix_C
  type(RealArray_C) :: Res_fn_h            !< Non-dimensional function of the ratio the first baroclinic
                                           !! deformation radius to the grid spacing at h points [nondim].
  type(RealArray_C) :: Res_fn_q            !< Non-dimensional function of the ratio the first baroclinic
                                           !! deformation radius to the grid spacing at q points [nondim].
  type(RealArray_C) :: BS_struct           !< Vertical structure function used in backscatter [nondim]
  logical(c_bool)   :: use_variable_mixing !< If true, use the variable mixing.
  logical(c_bool)   :: Resoln_scaled_Kh    !< If true, scale away the Laplacian viscosity
                                           !! when the deformation radius is well resolved.
end type hor_visc_VarMix_C

!> bind(C) mirror of hor_visc_STOCH_type, one RealArray_C per container (same order as
!! hor_visc_STOCH_type itself), followed by its scalars.
type, bind(C) :: hor_visc_STOCH_C
  type(RealArray_C) :: skeb_diss       !< Dissipation rate used to set amplitude of SKEBS [L2 T-3 ~> m2 s-3]
                                       !! Index into this at h points.
  logical(c_bool)   :: skeb_use_frict  !< If true, adds viscous dissipation rate to the amplitude of SKEBS
  real(c_double)    :: skeb_frict_coef !< If skeb_use_frict is true, then skeb_gm_coef * GM_work is added to the
                                       !! dissipation rate used to set the amplitude of SKEBS [nondim]
end type hor_visc_STOCH_C

contains

!> Calculates the acceleration due to the horizontal viscosity.
!!
!! A combination of biharmonic and Laplacian forms can be used. The coefficient
!! may either be a constant or a shear-dependent form. The biharmonic is
!! determined by twice taking the divergence of an appropriately defined stress
!! tensor. The Laplacian is determined by doing so once.
!!
!! To work, the following fields must be set outside of the usual
!! is:ie range before this subroutine is called:
!!   u(is-2:ie+2,js-2:je+2)
!!   v(is-2:ie+2,js-2:je+2)
!!   h(is-1:ie+1,js-1:je+1) or up to h(is-2:ie+2,js-2:je+2) with some Leith options.
subroutine horizontal_viscosity(u, v, h, uh, vh, diffu, diffv, MEKE, VarMix, G, GV, US, &
                                CS, nkblock, tv, dt, OBC, BT, TD, ADp, hu_cont, hv_cont, STOCH)
  type(ocean_grid_type), target, intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type), target, intent(in)  :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: uh      !< The zonal volume transport [H L2 T-1 ~> m3 s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: vh      !< The meridional volume transport [H L2 T-1 ~> m3 s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(out) :: diffu  !< Zonal acceleration due to convergence of
                                                       !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(out) :: diffv  !< Meridional acceleration due to convergence
                                                       !! of along-coordinate stress tensor [L T-2 ~> m s-2].
  type(MEKE_type),               intent(inout) :: MEKE !< MEKE fields
                                                       !! related to Mesoscale Eddy Kinetic Energy.
  type(VarMix_CS), target,       intent(inout) :: VarMix !< Variable mixing control structure
  type(unit_scale_type), target, intent(in)    :: US   !< A dimensional unit scaling type
  type(hor_visc_CS), target,     intent(inout) :: CS   !< Horizontal viscosity control structure
  integer,                       intent(in)    :: nkblock !< The effective k-block size [nondim]
  type(thermo_var_ptrs),         intent(in)    :: tv   !< A structure pointing to various
                                                       !! thermodynamic variables
  real,                          intent(in)    :: dt   !< Time increment [T ~> s]
  type(ocean_OBC_type), optional, pointer      :: OBC  !< Pointer to an open boundary condition type
  type(barotropic_CS), optional, intent(in)    :: BT   !< Barotropic control structure
  type(thickness_diffuse_CS), optional, intent(in) :: TD !< Thickness diffusion control structure
  type(accel_diag_ptrs), optional, intent(in)  :: ADp  !< Acceleration diagnostics
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                    optional, intent(inout) :: hu_cont !< Layer thickness at u-points [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                    optional, intent(inout) :: hv_cont !< Layer thickness at v-points [H ~> m or kg m-2].
  type(stochastic_CS), intent(inout), optional :: STOCH !< Stochastic control structure

  ! Local variables
  type(RealArray_t) :: u_layer_a ! One layer of the smoothed zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t) :: v_layer_a ! One layer of the smoothed meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t) :: u_a, v_a, h_a, uh_a, vh_a, diffu_a, diffv_a, hu_cont_a, hv_cont_a
  real, dimension(:,:), contiguous, pointer :: sh_xx_bt, sh_xy_bt, GME_effic_h, GME_effic_q
  real, dimension(:,:,:), contiguous, pointer :: u_smooth, v_smooth, dz, slope_x, slope_y, &
    KH_u_GME, KH_v_GME, GME_coeff_h, GME_coeff_q, str_xx_GME, str_xy_GME, NoSt, ShSt, FrictWork, &
    FrictWork_bh, FrictWork_GME, Ah_h, Ah_q, Kh_h, Kh_q, grid_Re_Kh, grid_Re_Ah, div_xx_h, &
    sh_xx_h, vort_xy_q, sh_xy_q, visc_limit_h, visc_limit_q, visc_limit_h_flag, visc_limit_q_flag, &
    visc_limit_h_frac, visc_limit_q_frac, BS_coeff_h, BS_coeff_q
  real, dimension(SZI_(G),SZJB_(G)) :: &
    dudx_bt, dvdy_bt ! components in the barotropic horizontal tension [T-1 ~> s-1]
  type(RealArray_t) :: &
    sh_xx_bt_a    ! barotropic horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
  type(RealArray_t) :: &
    str_xx_GME_a ! smoothed diagonal term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJ_(G),nkblock) :: &
    FrictWorkIntz, & ! depth integrated energy dissipated by lateral friction [R Z L2 T-3 ~> W m-2]
    FrictWorkIntz_bh   ! depth integrated energy dissipated by biharmonic lateral friction [R Z L2 T-3 ~> W m-2]
  type(RealArray_t) :: &
    GME_effic_h_a ! The filtered efficiency of the GME terms at h points [nondim]
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    dvdx_bt, dudy_bt   ! components in the barotropic shearing strain [T-1 ~> s-1]
  type(RealArray_t) :: &
    sh_xy_bt_a    ! barotropic horizontal shearing strain (du/dy + dv/dx) inc. metric terms [T-1 ~> s-1]
  type(RealArray_t) :: &
    str_xy_GME_a  ! smoothed cross term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
  type(RealArray_t) :: &
    GME_effic_q_a ! The filtered efficiency of the GME terms at q points [nondim]
  type(RealArray_t) :: &
    Ah_q_a, &      ! biharmonic viscosity at corner points [L4 T-1 ~> m4 s-1]
    Kh_q_a, &      ! Laplacian viscosity at corner points [L2 T-1 ~> m2 s-1]
    vort_xy_q_a, & ! vertical vorticity at corner points [T-1 ~> s-1]
    sh_xy_q_a,   & ! horizontal shearing strain at corner points [T-1 ~> s-1]
    GME_coeff_q_a, &  !< GME coeff. at q-points [L2 T-1 ~> m2 s-1]
    visc_limit_q_a, &  ! used to stabilize the EY24_EBT_BS backscatter [nondim]
    visc_limit_q_flag_a, & ! determines whether backscatter is shut off [nondim]
    visc_limit_q_frac_a, & ! determines how close backscatter is to shutting off [nondim]
    BS_coeff_q_a, &  ! A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
    ShSt_a       ! A diagnostic array of shear stress [T-1 ~> s-1].
  type(RealArray_t) :: &
    KH_u_GME_a, &  !< Isopycnal height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
    slope_x_a    !< Isopycnal slope in i-direction [Z L-1 ~> nondim]
  type(RealArray_t) :: &
    KH_v_GME_a, &  !< Isopycnal height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
    slope_y_a    !< Isopycnal slope in j-direction [Z L-1 ~> nondim]
  type(RealArray_t) :: &
    Ah_h_a, &          ! biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
    Kh_h_a, &          ! Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
    dz_a, &            ! Height change across layers [Z ~> m]
    FrictWork_a, &     ! work done by MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
    FrictWork_bh_a, &  ! work done by the biharmonic MKE dissipation mechanisms [R Z L2 T-3 ~> W m-2]
    FrictWork_GME_a, & ! work done by GME [R Z L2 T-3 ~> W m-2]
    div_xx_h_a,      & ! horizontal divergence [T-1 ~> s-1]
    sh_xx_h_a,       & ! horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    NoSt_a, &          ! A diagnostic array of normal stress [T-1 ~> s-1].
    BS_coeff_h_a     ! A diagnostic array of the backscatter coefficient [L2 T-1 ~> m2 s-1]
  type(RealArray_t) :: &
    grid_Re_Kh_a, &    ! Grid Reynolds number for Laplacian horizontal viscosity at h points [nondim]
    grid_Re_Ah_a, &    ! Grid Reynolds number for Biharmonic horizontal viscosity at h points [nondim]
    GME_coeff_h_a, &   ! GME coefficient at h-points [L2 T-1 ~> m2 s-1]
    visc_limit_h_a, &  ! Used to stabilize the EY24_EBT_BS backscatter [nondim]
    visc_limit_h_flag_a, & ! determines whether backscatter is shut off [nondim]
    visc_limit_h_frac_a  ! determines how close backscatter is to shutting off [nondim]
  type(RealArray_t) :: &
    u_smooth_a       ! Zonal velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
  type(RealArray_t) :: &
    v_smooth_a       ! Meridional velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
  logical :: rescale_Kh
  logical :: use_cont_huv
  integer :: is, ie, js, je, nz
  integer :: i, j, k
  type(Box_t) :: bxH0 ! The h-point box of the whole column, [isc:iec, jsc:jec, 1:nz]
  type(Box_t) :: bxQ0 ! The B-grid-index box of the whole column, [IscB:IecB, JscB:JecB, 1:nz]
  type(hor_visc_MEKE_type)   :: MEKE_a   ! Copies of the MEKE fields used by this call tree
  type(hor_visc_VarMix_type) :: VarMix_a ! Copies of the variable mixing fields used by this call tree
  type(hor_visc_STOCH_type)  :: STOCH_a  ! Copies of the stochastic fields used by this call tree
  type(grid_core_type), pointer :: Gcore ! Persistent copies of the core ocean grid fields
  type(hv_host_ctx), target :: host_ctx ! The objects used by the host helpers

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_hor_visc: Module must be initialized before it is used.")

  if (.not.(CS%hv_opts_CS%Laplacian .or. CS%hv_opts_CS%biharmonic)) return

  use_cont_huv = CS%hv_opts_CS%use_cont_thick .and. present(hu_cont) .and. present(hv_cont)
  if (use_cont_huv .and. .not.CS%use_cont_thick_bug) then
    call pass_vector(hu_cont, hv_cont, G%domain, To_All+Scalar_Pair, halo=2)
  endif

  rescale_Kh = .false.
  if (VarMix%use_variable_mixing) then
    rescale_Kh = VarMix%Resoln_scaled_Kh
    if ((rescale_Kh .or. CS%hv_opts_CS%res_scale_MEKE) &
        .and. (.not. allocated(VarMix%Res_fn_h) .or. .not. allocated(VarMix%Res_fn_q))) &
      call MOM_error(FATAL, "MOM_hor_visc: VarMix%Res_fn_h and VarMix%Res_fn_q "//&
                     "both need to be associated with Resoln_scaled_Kh or RES_SCALE_MEKE_VISC.")
  elseif (CS%hv_opts_CS%res_scale_MEKE) then
    call MOM_error(FATAL, "MOM_hor_visc: VarMix needs to be associated if "//&
                          "RES_SCALE_MEKE_VISC is True.")
  endif

  ! Set the halo sizes used for the vorticity calculations.
  if ((CS%hv_opts_CS%Leith_Kh) .or. (CS%hv_opts_CS%Leith_Ah) .or. (CS%hv_opts_CS%use_Leithy)) then
    if ((G%isc-G%isd < 3) .or. (G%isc-G%isd < 3)) call MOM_error(FATAL, &
          "The minimum halo size is 3 when a Leith viscosity is being used.")
  endif

  Gcore => grid_core(G)
  call bxH0%safe_alloc(ndims=3) ; call bxQ0%safe_alloc(ndims=3)
  call bxH0%set(idxS=[G%isc,G%jsc,1], idxE=[G%iec,G%jec,nz])
  call bxQ0%set(idxS=[G%IscB,G%JscB,1], idxE=[G%IecB,G%JecB,nz])
  call u_smooth_a%alloc(lb=[LBOUND(u,1),LBOUND(h,2),1], ub=[UBOUND(u,1),UBOUND(h,2),SZK_(GV)])
  call v_smooth_a%alloc(lb=[LBOUND(h,1),LBOUND(v,2),1], ub=[UBOUND(h,1),UBOUND(v,2),SZK_(GV)])
  call dz_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(GV)])
  call slope_x_a%alloc(lb=[LBOUND(u,1),LBOUND(h,2),1], ub=[UBOUND(u,1),UBOUND(h,2),SZK_(GV)+1])
  call slope_y_a%alloc(lb=[LBOUND(h,1),LBOUND(v,2),1], ub=[UBOUND(h,1),UBOUND(v,2),SZK_(GV)+1])
  call sh_xx_bt_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2)], ub=[UBOUND(h,1),UBOUND(h,2)])
  call sh_xy_bt_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2)], ub=[UBOUND(u,1),UBOUND(v,2)])
  call GME_effic_h_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2)], ub=[UBOUND(h,1),UBOUND(h,2)])
  call GME_effic_q_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2)], ub=[UBOUND(u,1),UBOUND(v,2)])
  call KH_u_GME_a%alloc(lb=[LBOUND(u,1),LBOUND(h,2),1], ub=[UBOUND(u,1),UBOUND(h,2),SZK_(GV)+1])
  call KH_v_GME_a%alloc(lb=[LBOUND(h,1),LBOUND(v,2),1], ub=[UBOUND(h,1),UBOUND(v,2),SZK_(GV)+1])
  call GME_coeff_h_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(G)])
  call GME_coeff_q_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2),1], ub=[UBOUND(u,1),UBOUND(v,2),SZK_(GV)])
  call str_xx_GME_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),nkblock])
  call str_xy_GME_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2),1], ub=[UBOUND(u,1),UBOUND(v,2),nkblock])
  call NoSt_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(GV)])
  call ShSt_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2),1], ub=[UBOUND(u,1),UBOUND(v,2),SZK_(GV)])
  call FrictWork_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(GV)])
  call FrictWork_bh_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(GV)])
  call FrictWork_GME_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(GV)])
  call Ah_h_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(GV)])
  call Ah_q_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2),1], ub=[UBOUND(u,1),UBOUND(v,2),SZK_(GV)])
  call Kh_h_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(GV)])
  call Kh_q_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2),1], ub=[UBOUND(u,1),UBOUND(v,2),SZK_(GV)])
  call grid_Re_Kh_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(G)])
  call grid_Re_Ah_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(G)])
  call div_xx_h_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(GV)])
  call sh_xx_h_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(GV)])
  call vort_xy_q_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2),1], ub=[UBOUND(u,1),UBOUND(v,2),SZK_(GV)])
  call sh_xy_q_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2),1], ub=[UBOUND(u,1),UBOUND(v,2),SZK_(GV)])
  call visc_limit_h_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(G)])
  call visc_limit_q_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2),1], ub=[UBOUND(u,1),UBOUND(v,2),SZK_(GV)])
  call visc_limit_h_flag_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], &
                                 ub=[UBOUND(h,1),UBOUND(h,2),SZK_(G)])
  call visc_limit_q_flag_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2),1], &
                                 ub=[UBOUND(u,1),UBOUND(v,2),SZK_(GV)])
  call visc_limit_h_frac_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], &
                                 ub=[UBOUND(h,1),UBOUND(h,2),SZK_(G)])
  call visc_limit_q_frac_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2),1], &
                                 ub=[UBOUND(u,1),UBOUND(v,2),SZK_(GV)])
  call BS_coeff_h_a%alloc(lb=[LBOUND(h,1),LBOUND(h,2),1], ub=[UBOUND(h,1),UBOUND(h,2),SZK_(GV)])
  call BS_coeff_q_a%alloc(lb=[LBOUND(u,1),LBOUND(v,2),1], ub=[UBOUND(u,1),UBOUND(v,2),SZK_(GV)])
  call u_smooth_a%view(u_smooth) ; call v_smooth_a%view(v_smooth) ; call dz_a%view(dz)
  call slope_x_a%view(slope_x) ; call slope_y_a%view(slope_y) ; call sh_xx_bt_a%view(sh_xx_bt)
  call sh_xy_bt_a%view(sh_xy_bt) ; call GME_effic_h_a%view(GME_effic_h)
  call GME_effic_q_a%view(GME_effic_q) ; call KH_u_GME_a%view(KH_u_GME)
  call KH_v_GME_a%view(KH_v_GME) ; call GME_coeff_h_a%view(GME_coeff_h)
  call GME_coeff_q_a%view(GME_coeff_q) ; call str_xx_GME_a%view(str_xx_GME)
  call str_xy_GME_a%view(str_xy_GME) ; call NoSt_a%view(NoSt) ; call ShSt_a%view(ShSt)
  call FrictWork_a%view(FrictWork) ; call FrictWork_bh_a%view(FrictWork_bh)
  call FrictWork_GME_a%view(FrictWork_GME) ; call Ah_h_a%view(Ah_h) ; call Ah_q_a%view(Ah_q)
  call Kh_h_a%view(Kh_h) ; call Kh_q_a%view(Kh_q) ; call grid_Re_Kh_a%view(grid_Re_Kh)
  call grid_Re_Ah_a%view(grid_Re_Ah) ; call div_xx_h_a%view(div_xx_h)
  call sh_xx_h_a%view(sh_xx_h) ; call vort_xy_q_a%view(vort_xy_q) ; call sh_xy_q_a%view(sh_xy_q)
  call visc_limit_h_a%view(visc_limit_h) ; call visc_limit_q_a%view(visc_limit_q)
  call visc_limit_h_flag_a%view(visc_limit_h_flag)
  call visc_limit_q_flag_a%view(visc_limit_q_flag)
  call visc_limit_h_frac_a%view(visc_limit_h_frac)
  call visc_limit_q_frac_a%view(visc_limit_q_frac) ; call BS_coeff_h_a%view(BS_coeff_h)
  call BS_coeff_q_a%view(BS_coeff_q)

  if (CS%hv_opts_CS%use_GME) then
    call hor_visc_GME_setup(G, GV, US, CS, Gcore, bxH0, bxQ0, h, BT, TD, nkblock, dudx_bt, dvdy_bt, &
                            dvdx_bt, dudy_bt, sh_xx_bt, sh_xy_bt, GME_effic_h, GME_effic_q, &
                            KH_u_GME, KH_v_GME, GME_coeff_h, GME_coeff_q, str_xx_GME, str_xy_GME)
    ! NOTE: sh_xx_bt, sh_xy_bt, GME_effic_h, GME_effic_q, KH_u_GME and KH_v_GME
    !   are computed on the host above and are only read (never written) on the
    !   device.  dudx_bt, dvdy_bt, dvdx_bt and dudy_bt are only used on the
    !   host (post_data diagnostics), so they are not mapped at all.
    ! NOTE: GME_coeff_h, GME_coeff_q, str_xx_GME and str_xy_GME are written and
    !   read entirely on the device (aside from the smooth_GME host round-trip
    !   below), so a bare alloc is enough on entry.
  endif

  if (CS%hv_opts_CS%use_Leithy) then
    ! Smooth the velocity. Right now it happens twice. In the future
    ! one might make the number of smoothing cycles a user-specified parameter
    call u_layer_a%alloc(lb=[LBOUND(u_smooth,1),LBOUND(u_smooth,2)], ub=[UBOUND(u_smooth,1),UBOUND(u_smooth,2)])
    call v_layer_a%alloc(lb=[LBOUND(v_smooth,1),LBOUND(v_smooth,2)], ub=[UBOUND(v_smooth,1),UBOUND(v_smooth,2)])
    do k=1,nz
      ! One call applies the filter twice
      u_smooth(:,:,k) = u(:,:,k)
      v_smooth(:,:,k) = v(:,:,k)
      call u_layer_a%copy2Array(u_smooth(:,:,k)) ; call v_layer_a%copy2Array(v_smooth(:,:,k))
      call smooth_x9_uv(Gcore, bxH0, bxQ0, u_layer_a, v_layer_a, zero_land=.false.)
      call u_layer_a%copy2F(u_smooth(:,:,k)) ; call v_layer_a%copy2F(v_smooth(:,:,k))
    enddo
    call u_layer_a%free() ; call v_layer_a%free()
    call pass_vector(u_smooth, v_smooth, G%Domain)
  endif

  if (CS%hv_opts_CS%use_QG_Leith_visc .and. ((CS%hv_opts_CS%Leith_Kh) .or. (CS%hv_opts_CS%Leith_Ah))) then
    call thickness_to_dz(h, tv, dz, G, GV, US, halo_size=2)
    ! Calculate isopycnal slopes that will be used for some forms of viscosity.
    call calc_QG_slopes(h, tv, dt, G, GV, US, slope_x, slope_y, VarMix, OBC)
    ! If the following halo update is added, the calculations in calc_QG_slopes could work on just
    ! the computational domains, and some halo updates outside of this routine could be smaller.
    ! call pass_vector(slope_x, slope_y, G%Domain, halo=2)
  endif

  call MEKE_a%build_from(MEKE) ; call VarMix_a%build_from(VarMix)
  if (present(STOCH)) call STOCH_a%build_from(STOCH)
  host_ctx%G => G ; host_ctx%GV => GV ; host_ctx%US => US
  host_ctx%VarMix => VarMix ; host_ctx%CS => CS
  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), &
                 source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), &
                 source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), &
                 source=h)
  call uh_a%alloc(lb=LBOUND(uh), ub=UBOUND(uh), &
                  source=uh)
  call vh_a%alloc(lb=LBOUND(vh), ub=UBOUND(vh), &
                  source=vh)
  call diffu_a%alloc(lb=LBOUND(diffu), ub=UBOUND(diffu), &
                     source=diffu)
  call diffv_a%alloc(lb=LBOUND(diffv), ub=UBOUND(diffv), &
                     source=diffv)
  if (present(hu_cont)) &
    call hu_cont_a%alloc(lb=LBOUND(hu_cont), ub=UBOUND(hu_cont), &
                         source=hu_cont)
  if (present(hv_cont)) &
    call hv_cont_a%alloc(lb=LBOUND(hv_cont), ub=UBOUND(hv_cont), &
                         source=hv_cont)
  call horizontal_viscosity_TR(bxH0, bxQ0, u_a, v_a, h_a, uh_a, vh_a, diffu_a, diffv_a, MEKE_a, &
                               VarMix_a, Gcore, CS%hv_opts_CS, CS%hv_metrics_CS, nkblock, OBC, &
                               hu_cont_a, hv_cont_a, STOCH_a, u_smooth_a, v_smooth_a, dz_a, &
                               slope_x_a, slope_y_a, sh_xx_bt_a, sh_xy_bt_a, GME_effic_h_a, &
                               GME_effic_q_a, KH_u_GME_a, KH_v_GME_a, GME_coeff_h_a, GME_coeff_q_a, &
                               str_xx_GME_a, str_xy_GME_a, NoSt_a, ShSt_a, FrictWork_a, &
                               FrictWork_bh_a, FrictWork_GME_a, Ah_h_a, Ah_q_a, Kh_h_a, Kh_q_a, &
                               grid_Re_Kh_a, grid_Re_Ah_a, div_xx_h_a, sh_xx_h_a, vort_xy_q_a, &
                               sh_xy_q_a, visc_limit_h_a, visc_limit_q_a, visc_limit_h_flag_a, &
                               visc_limit_q_flag_a, visc_limit_h_frac_a, visc_limit_q_frac_a, &
                               BS_coeff_h_a, BS_coeff_q_a, GV%H_subroundoff, GV%H_to_RZ, US%s_to_T, &
                               c_loc(host_ctx))
  call diffu_a%copy2F(diffu) ; call diffv_a%copy2F(diffv)
  if (present(hu_cont)) call hu_cont_a%copy2F(hu_cont)
  if (present(hv_cont)) call hv_cont_a%copy2F(hv_cont)
  call u_a%free() ; call v_a%free() ; call h_a%free() ; call uh_a%free() ; call vh_a%free()
  call diffu_a%free() ; call diffv_a%free() ; call hu_cont_a%free() ; call hv_cont_a%free()
  call bxH0%free() ; call bxQ0%free()
  call MEKE_a%copy_back(MEKE) ; call VarMix_a%free()
  if (present(STOCH)) call STOCH_a%copy_back(STOCH)

  ! Offer fields for diagnostic averaging.
  if (CS%id_normstress > 0) call post_data(CS%id_normstress, NoSt, CS%diag)
  if (CS%id_shearstress > 0) call post_data(CS%id_shearstress, ShSt, CS%diag)
  if (CS%id_diffu>0)     call post_data(CS%id_diffu, diffu, CS%diag)
  if (CS%id_diffv>0)     call post_data(CS%id_diffv, diffv, CS%diag)
  if (CS%id_FrictWork>0) call post_data(CS%id_FrictWork, FrictWork, CS%diag)
  if (CS%id_FrictWork_bh>0) call post_data(CS%id_FrictWork_bh, FrictWork_bh, CS%diag)
  if (CS%id_Ah_h>0)      call post_data(CS%id_Ah_h, Ah_h, CS%diag)
  if (CS%id_grid_Re_Ah>0) call post_data(CS%id_grid_Re_Ah, grid_Re_Ah, CS%diag)
  if (CS%id_div_xx_h>0)  call post_data(CS%id_div_xx_h, div_xx_h, CS%diag)
  if (CS%id_vort_xy_q>0) call post_data(CS%id_vort_xy_q, vort_xy_q, CS%diag)
  if (CS%id_sh_xx_h>0)   call post_data(CS%id_sh_xx_h, sh_xx_h, CS%diag)
  if (CS%id_sh_xy_q>0)   call post_data(CS%id_sh_xy_q, sh_xy_q, CS%diag)
  if (CS%id_Ah_q>0)      call post_data(CS%id_Ah_q, Ah_q, CS%diag)
  if (CS%id_Kh_h>0)      call post_data(CS%id_Kh_h, Kh_h, CS%diag)
  if (CS%id_grid_Re_Kh>0) call post_data(CS%id_grid_Re_Kh, grid_Re_Kh, CS%diag)
  if (CS%id_Kh_q>0)      call post_data(CS%id_Kh_q, Kh_q, CS%diag)
  if (CS%hv_opts_CS%use_GME) then  ! post barotropic tension and strain
    if (CS%id_GME_coeff_h > 0) call post_data(CS%id_GME_coeff_h, GME_coeff_h, CS%diag)
    if (CS%id_GME_coeff_q > 0) call post_data(CS%id_GME_coeff_q, GME_coeff_q, CS%diag)
    if (CS%id_FrictWork_GME>0) call post_data(CS%id_FrictWork_GME, FrictWork_GME, CS%diag)
    if (CS%id_dudx_bt > 0) call post_data(CS%id_dudx_bt, dudx_bt, CS%diag)
    if (CS%id_dvdy_bt > 0) call post_data(CS%id_dvdy_bt, dvdy_bt, CS%diag)
    if (CS%id_dudy_bt > 0) call post_data(CS%id_dudy_bt, dudy_bt, CS%diag)
    if (CS%id_dvdx_bt > 0) call post_data(CS%id_dvdx_bt, dvdx_bt, CS%diag)
  endif
  if (CS%hv_opts_CS%EY24_EBT_BS) then
    if (CS%id_visc_limit_h>0) call post_data(CS%id_visc_limit_h, visc_limit_h, CS%diag)
    if (CS%id_visc_limit_q>0) call post_data(CS%id_visc_limit_q, visc_limit_q, CS%diag)
    if (CS%id_visc_limit_h_frac>0) call post_data(CS%id_visc_limit_h_frac, visc_limit_h_frac, CS%diag)
    if (CS%id_visc_limit_q_frac>0) call post_data(CS%id_visc_limit_q_frac, visc_limit_q_frac, CS%diag)
    if (CS%id_visc_limit_h_flag>0) call post_data(CS%id_visc_limit_h_flag, visc_limit_h_flag, CS%diag)
    if (CS%id_visc_limit_q_flag>0) call post_data(CS%id_visc_limit_q_flag, visc_limit_q_flag, CS%diag)
    if (CS%id_BS_coeff_h>0) call post_data(CS%id_BS_coeff_h, BS_coeff_h, CS%diag)
    if (CS%id_BS_coeff_q>0) call post_data(CS%id_BS_coeff_q, BS_coeff_q, CS%diag)
  endif

  if (CS%hv_opts_CS%debug) then
    if (CS%hv_opts_CS%Laplacian) then
      call hchksum(Kh_h, "Kh_h", G%HI, haloshift=1, unscale=US%L_to_m**2*US%s_to_T)
      call Bchksum(Kh_q, "Kh_q", G%HI, haloshift=0, symmetric=.true., unscale=US%L_to_m**2*US%s_to_T)
    endif
    if (CS%hv_opts_CS%biharmonic) then
      call hchksum(Ah_h, "Ah_h", G%HI, haloshift=1, unscale=US%L_to_m**4*US%s_to_T)
      call Bchksum(Ah_q, "Ah_q", G%HI, haloshift=0, symmetric=.true., unscale=US%L_to_m**4*US%s_to_T)
    endif
  endif

  if (CS%id_FrictWorkIntz > 0) then
    do j=js,je
      do i=is,ie ; FrictWorkIntz(i,j,1) = FrictWork(i,j,1) ; enddo
      do k=2,nz ; do i=is,ie
        FrictWorkIntz(i,j,1) = FrictWorkIntz(i,j,1) + FrictWork(i,j,k)
      enddo ; enddo
    enddo
    call post_data(CS%id_FrictWorkIntz, FrictWorkIntz(:,:,1), CS%diag)
  endif

  if (CS%id_FrictWorkIntz_bh > 0) then
    do j=js,je
      do i=is,ie ; FrictWorkIntz_bh(i,j,1) = FrictWork_bh(i,j,1) ; enddo
      do k=2,nz ; do i=is,ie
        FrictWorkIntz_bh(i,j,1) = FrictWorkIntz_bh(i,j,1) + FrictWork_bh(i,j,k)
      enddo ; enddo
    enddo
    call post_data(CS%id_FrictWorkIntz_bh, FrictWorkIntz_bh(:,:,1), CS%diag)
  endif

  if (present(ADp)) then
    ! Diagnostics of the fractional thicknesses times momentum budget terms
    ! 3D diagnostics of hf_diffu(diffv) are commented because there is no clarity on proper remapping grid option.
    ! The code is retained for debugging purposes in the future.
    !if (CS%id_hf_diffu > 0) call post_product_u(CS%id_hf_diffu, diffu, ADp%diag_hfrac_u, G, nz, CS%diag)
    !if (CS%id_hf_diffv > 0) call post_product_v(CS%id_hf_diffv, diffv, ADp%diag_hfrac_v, G, nz, CS%diag)

    ! Diagnostics for thickness-weighted vertically averaged momentum budget terms
    if (CS%id_hf_diffu_2d > 0) call post_product_sum_u(CS%id_hf_diffu_2d, diffu, ADp%diag_hfrac_u, G, nz, CS%diag)
    if (CS%id_hf_diffv_2d > 0) call post_product_sum_v(CS%id_hf_diffv_2d, diffv, ADp%diag_hfrac_v, G, nz, CS%diag)

    ! Diagnostics for the vertical sum of layer thickness x momentum budget terms
    if (CS%id_intz_diffu_2d > 0) call post_product_sum_u(CS%id_intz_diffu_2d, diffu, ADp%diag_hu, G, nz, CS%diag)
    if (CS%id_intz_diffv_2d > 0) call post_product_sum_v(CS%id_intz_diffv_2d, diffv, ADp%diag_hv, G, nz, CS%diag)

    ! Diagnostics for thickness x momentum budget terms
    if (CS%id_h_diffu > 0) call post_product_u(CS%id_h_diffu, diffu, ADp%diag_hu, G, nz, CS%diag)
    if (CS%id_h_diffv > 0) call post_product_v(CS%id_h_diffv, diffv, ADp%diag_hv, G, nz, CS%diag)

    ! Diagnostics for momentum budget terms multiplied by visc_rem_[uv],
    if (CS%id_diffu_visc_rem > 0) call post_product_u(CS%id_diffu_visc_rem, diffu, ADp%visc_rem_u, G, nz, CS%diag)
    if (CS%id_diffv_visc_rem > 0) call post_product_v(CS%id_diffv_visc_rem, diffv, ADp%visc_rem_v, G, nz, CS%diag)
  endif

  if (CS%hv_opts_CS%use_ZB2020) then
    call ZB2020_lateral_stress(u, v, h, diffu, diffv, G, GV, CS%ZB2020, &
                               CS%dx2h, CS%dy2h, CS%dx2q, CS%dy2q)
  endif

  call u_smooth_a%free() ; call v_smooth_a%free() ; call dz_a%free() ; call slope_x_a%free()
  call slope_y_a%free() ; call sh_xx_bt_a%free() ; call sh_xy_bt_a%free()
  call GME_effic_h_a%free() ; call GME_effic_q_a%free() ; call KH_u_GME_a%free()
  call KH_v_GME_a%free() ; call GME_coeff_h_a%free() ; call GME_coeff_q_a%free()
  call str_xx_GME_a%free() ; call str_xy_GME_a%free() ; call NoSt_a%free() ; call ShSt_a%free()
  call FrictWork_a%free() ; call FrictWork_bh_a%free() ; call FrictWork_GME_a%free()
  call Ah_h_a%free() ; call Ah_q_a%free() ; call Kh_h_a%free() ; call Kh_q_a%free()
  call grid_Re_Kh_a%free() ; call grid_Re_Ah_a%free() ; call div_xx_h_a%free()
  call sh_xx_h_a%free() ; call vort_xy_q_a%free() ; call sh_xy_q_a%free()
  call visc_limit_h_a%free() ; call visc_limit_q_a%free() ; call visc_limit_h_flag_a%free()
  call visc_limit_q_flag_a%free() ; call visc_limit_h_frac_a%free()
  call visc_limit_q_frac_a%free() ; call BS_coeff_h_a%free() ; call BS_coeff_q_a%free()

end subroutine horizontal_viscosity

!> The tree root of horizontal_viscosity: calculates the acceleration due to the horizontal
!! viscosity once the Fortran-only work (halo updates, calls into other modules, error checks
!! and diagnostics) has been done by the horizontal_viscosity wrapper.
subroutine horizontal_viscosity_TR(bxH0, bxQ0, u_a, v_a, h_a, uh_a, vh_a, diffu_a, diffv_a, MEKE_a, &
                                   VarMix_a, Gcore, CS, metrics, nkblock, OBC, hu_cont_a, &
                                   hv_cont_a, STOCH_a, u_smooth_a, v_smooth_a, dz_a, slope_x_a, &
                                   slope_y_a, sh_xx_bt_a, sh_xy_bt_a, GME_effic_h_a, GME_effic_q_a, &
                                   KH_u_GME_a, KH_v_GME_a, GME_coeff_h_a, GME_coeff_q_a, &
                                   str_xx_GME_a, str_xy_GME_a, NoSt_a, ShSt_a, FrictWork_a, &
                                   FrictWork_bh_a, FrictWork_GME_a, Ah_h_a, Ah_q_a, Kh_h_a, Kh_q_a, &
                                   grid_Re_Kh_a, grid_Re_Ah_a, div_xx_h_a, sh_xx_h_a, vort_xy_q_a, &
                                   sh_xy_q_a, visc_limit_h_a, visc_limit_q_a, visc_limit_h_flag_a, &
                                   visc_limit_q_flag_a, visc_limit_h_frac_a, visc_limit_q_frac_a, &
                                   BS_coeff_h_a, BS_coeff_q_a, H_subroundoff, H_to_RZ, s_to_T, host)
  type(Box_t),                   intent(in)    :: bxH0 !< The h-point box of the whole column,
                                                       !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),                   intent(in)    :: bxQ0 !< The B-grid-index box of the whole column,
                                                       !! [IscB:IecB, JscB:JecB, 1:nz]
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  real,                 intent(in)    :: H_subroundoff !< A thickness that is so small that it can be added
                                                       !! to a thickness of Angstrom or larger without
                                                       !! changing it at the bit level [H ~> m or kg m-2]
  real,                 intent(in)    :: H_to_RZ !< A constant that translates thickness units to the units
                                                 !! of mass per unit area [R Z H-1 ~> kg m-3 or 1]
  type(RealArray_t), intent(in)    :: u_a !< The zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: v_a !< The meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: h_a !< Layer thicknesses [H ~> m or kg m-2].
  type(RealArray_t), intent(in)    :: uh_a !< The zonal volume transport [H L2 T-1 ~> m3 s-1].
  type(RealArray_t), intent(in)    :: vh_a !< The meridional volume transport [H L2 T-1 ~> m3 s-1].
  type(RealArray_t), intent(inout) :: diffu_a !< Zonal acceleration due to convergence of
                                              !! along-coordinate stress tensor [L T-2 ~> m s-2]
  type(RealArray_t), intent(inout) :: diffv_a !< Meridional acceleration due to convergence
                                              !! of along-coordinate stress tensor [L T-2 ~> m s-2].
  type(hor_visc_MEKE_type),      intent(inout) :: MEKE_a !< MEKE fields
                                                       !! related to Mesoscale Eddy Kinetic Energy.
  type(hor_visc_VarMix_type),    intent(inout) :: VarMix_a !< Variable mixing control structure
  type(c_ptr), value,   intent(in)    :: host  !< Opaque pointer to the hv_host_ctx that holds the Fortran
                                               !! objects used by the host helpers
  real,                 intent(in)    :: s_to_T  !< A constant that translates time intervals in seconds to
                                                 !! the units of time [T s-1 ~> 1]
  type(hv_opts_CS),     intent(in)    :: CS      !< Options used by the horizontal viscosity kernels
  type(hv_metrics_CS),  intent(in)    :: metrics !< Persistent copies of the horizontal
                                                 !! viscosity metric and coefficient arrays
  integer,                       intent(in)    :: nkblock !< The effective k-block size [nondim]
  type(ocean_OBC_type), optional, pointer      :: OBC  !< Pointer to an open boundary condition type
  type(RealArray_t), intent(inout) :: hu_cont_a !< Layer thickness at u-points [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: hv_cont_a !< Layer thickness at v-points [H ~> m or kg m-2].
  type(hor_visc_STOCH_type), intent(inout)           :: STOCH_a !< Stochastic control structure
  type(RealArray_t), intent(in)    :: u_smooth_a !< Zonal velocity, smoothed with a spatial low-pass filter
                                                 !! [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_smooth_a !< Meridional velocity, smoothed with a spatial low-pass
                                                 !! filter [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: dz_a !< Height change across layers [Z ~> m]
  type(RealArray_t), intent(inout) :: slope_x_a !< Isopycnal slope in i-direction [Z L-1 ~> nondim]
  type(RealArray_t), intent(inout) :: slope_y_a !< Isopycnal slope in j-direction [Z L-1 ~> nondim]
  type(RealArray_t), intent(in)    :: sh_xx_bt_a !< barotropic horizontal tension (du/dx - dv/dy) including
                                                 !! metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: sh_xy_bt_a !< barotropic horizontal shearing strain (du/dy + dv/dx) inc.
                                                 !! metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: GME_effic_h_a !< The filtered efficiency of the GME terms at h points
                                                    !! [nondim]
  type(RealArray_t), intent(in)    :: GME_effic_q_a !< The filtered efficiency of the GME terms at q points
                                                    !! [nondim]
  type(RealArray_t), intent(in)    :: KH_u_GME_a !< Isopycnal height diffusivities in u-columns
                                                 !! [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(in)    :: KH_v_GME_a !< Isopycnal height diffusivities in v-columns
                                                 !! [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: GME_coeff_h_a !< GME coefficient at h-points [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: GME_coeff_q_a !< GME coeff. at q-points [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: str_xx_GME_a !< smoothed diagonal term in the stress tensor from GME
                                                   !! [L2 T-2 ~> m2 s-2]
  type(RealArray_t), intent(inout) :: str_xy_GME_a !< smoothed cross term in the stress tensor from GME
                                                   !! [L2 T-2 ~> m2 s-2]
  type(RealArray_t), intent(inout) :: NoSt_a !< A diagnostic array of normal stress [T-1 ~> s-1].
  type(RealArray_t), intent(inout) :: ShSt_a !< A diagnostic array of shear stress [T-1 ~> s-1].
  type(RealArray_t), intent(inout) :: FrictWork_a !< work done by MKE dissipation mechanisms
                                                  !! [R Z L2 T-3 ~> W m-2]
  type(RealArray_t), intent(inout) :: FrictWork_bh_a !< work done by the biharmonic MKE dissipation mechanisms
                                                     !! [R Z L2 T-3 ~> W m-2]
  type(RealArray_t), intent(inout) :: FrictWork_GME_a !< work done by GME [R Z L2 T-3 ~> W m-2]
  type(RealArray_t), intent(inout) :: Ah_h_a !< biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
  type(RealArray_t), intent(inout) :: Ah_q_a !< biharmonic viscosity at corner points [L4 T-1 ~> m4 s-1]
  type(RealArray_t), intent(inout) :: Kh_h_a !< Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: Kh_q_a !< Laplacian viscosity at corner points [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: grid_Re_Kh_a !< Grid Reynolds number for Laplacian horizontal viscosity
                                                   !! at h points [nondim]
  type(RealArray_t), intent(inout) :: grid_Re_Ah_a !< Grid Reynolds number for Biharmonic horizontal viscosity
                                                   !! at h points [nondim]
  type(RealArray_t), intent(inout) :: div_xx_h_a !< horizontal divergence [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xx_h_a !< horizontal tension (du/dx - dv/dy) including metric terms
                                                !! [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: vort_xy_q_a !< vertical vorticity at corner points [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xy_q_a !< horizontal shearing strain at corner points [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: visc_limit_h_a !< Used to stabilize the EY24_EBT_BS backscatter [nondim]
  type(RealArray_t), intent(inout) :: visc_limit_q_a !< used to stabilize the EY24_EBT_BS backscatter [nondim]
  type(RealArray_t), intent(inout) :: visc_limit_h_flag_a !< determines whether backscatter is shut off
                                                          !! [nondim]
  type(RealArray_t), intent(inout) :: visc_limit_q_flag_a !< determines whether backscatter is shut off
                                                          !! [nondim]
  type(RealArray_t), intent(inout) :: visc_limit_h_frac_a !< determines how close backscatter is to shutting
                                                          !! off [nondim]
  type(RealArray_t), intent(inout) :: visc_limit_q_frac_a !< determines how close backscatter is to shutting
                                                          !! off [nondim]
  type(RealArray_t), intent(inout) :: BS_coeff_h_a !< A diagnostic array of the backscatter coefficient
                                                   !! [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: BS_coeff_q_a !< A diagnostic array of the backscatter coefficient
                                                   !! [L2 T-1 ~> m2 s-1]

  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: h_u, Del2u, h_v, Del2v, sh_xx_smooth, sh_xx, dvdy, &
    dudx, grad_div_mag_h, grad_vort_mag_h_2d, grad_vort_mag_h, bhstr_xx, str_xx, m_leithy, &
    dDel2udy, dDel2vdx, dudy_smooth, dvdx_smooth, dudy, dvdx, sh_xy_smooth, sh_xy, hq, &
    grad_div_mag_q, Del2vort_q, grad_vort_mag_q_2d, grad_vort_mag_q, vort_xy_smooth, vort_xy, &
    bhstr_xy, str_xy, visc_bound_rem, hrat_min, vert_vort_mag_smooth, vert_vort_mag, Shear_mag, Kh, &
    Ah, div_xx_dy, vort_xy_dx_smooth, vort_xy_dx, div_xx_dx, vort_xy_dy_smooth, vort_xy_dy, Ah_sq
  real, dimension(:,:), contiguous, pointer :: GME_effic_q, GME_effic_h, sh_xy_bt, sh_xx_bt
  real, dimension(:,:,:), contiguous, pointer :: BS_coeff_q, BS_coeff_h, visc_limit_q_frac, &
    visc_limit_h_frac, visc_limit_q_flag, visc_limit_h_flag, visc_limit_q, visc_limit_h, sh_xy_q, &
    vort_xy_q, sh_xx_h, div_xx_h, grid_Re_Ah, grid_Re_Kh, Kh_q, Kh_h, Ah_q, Ah_h, FrictWork_GME, &
    FrictWork_bh, FrictWork, ShSt, NoSt, str_xy_GME, str_xx_GME, GME_coeff_q, GME_coeff_h, &
    KH_v_GME, KH_u_GME, slope_y, slope_x, dz, v_smooth, u_smooth, hv_cont, hu_cont, diffv, diffu, &
    vh, uh, h, v, u
  real, dimension(:,:), contiguous, pointer :: MEKE_mom_src, MEKE_mom_src_bh, MEKE_GME_snk
  type(RealArray_t) :: &
    Del2u_a, &      ! The u-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
    h_u_a           ! Thickness interpolated to u points [H ~> m or kg m-2].
  type(RealArray_t) :: &
    Del2v_a, &      ! The v-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
    h_v_a           ! Thickness interpolated to v points [H ~> m or kg m-2].
  type(RealArray_t) :: &
    sh_xx_a, &      ! horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    sh_xx_smooth_a  ! horizontal tension from smoothed velocity including metric terms [T-1 ~> s-1]
  type(RealArray_t) :: &
    str_xx_a, &     ! str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but
                  ! at some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2].
    bhstr_xx_a, &   ! A copy of str_xx that only contains the biharmonic contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
    grad_vort_mag_h_a, & ! Magnitude of vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
    grad_vort_mag_h_2d_a, & ! Magnitude of 2d vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
    grad_div_mag_h_a, &     ! Magnitude of divergence gradient at h-points [L-1 T-1 ~> m-1 s-1]
    dudx_a, dvdy_a      ! components in the horizontal tension [T-1 ~> s-1]
  type(RealArray_t) :: &
    m_leithy_a      ! Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]

  type(RealArray_t) :: &
    dvdx_a, dudy_a, & ! components in the shearing strain [T-1 ~> s-1]
    dvdx_smooth_a, dudy_smooth_a, & ! components in the shearing strain from smoothed velocity [T-1 ~> s-1]
    dDel2vdx_a, dDel2udy_a ! Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
  type(RealArray_t) :: &
    sh_xy_a,  &     ! horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
    sh_xy_smooth_a  ! horizontal shearing strain from smoothed velocity including metric terms [T-1 ~> s-1]
  type(RealArray_t) :: &
    str_xy_a, &     ! str_xy is the cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but
                  ! at some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2].
    bhstr_xy_a, &   ! A copy of str_xy that only contains the biharmonic contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
    vort_xy_a, &    ! Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
    vort_xy_smooth_a, & ! Vertical vorticity including metric terms, smoothed [T-1 ~> s-1]
    grad_vort_mag_q_a, & ! Magnitude of vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
    grad_vort_mag_q_2d_a, & ! Magnitude of 2d vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
    Del2vort_q_a, & ! Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
    grad_div_mag_q_a, &  ! Magnitude of divergence gradient at q-points [L-1 T-1 ~> m-1 s-1]
    hq_a             ! harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2]
                   ! This form guarantees that hq/hu < 4.

  real :: h_neglect  ! thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
  real :: h_neglect3 ! h_neglect^3 [H3 ~> m3 or kg3 m-6]
  real :: backscat_subround ! The ratio of f over Shear_mag that is so small that the backscatter
                    ! calculation gives the same value as if f were 0 [nondim].

  logical :: rescale_Kh
  logical :: find_FrictWork
  logical :: find_FrictWork_bh
  logical :: apply_OBC = .false.
  logical :: apply_OBC_strain
  logical :: use_MEKE_Ku
  logical :: use_MEKE_Au
  logical :: skeb_use_frict
  logical :: use_cont_huv
  logical :: use_kh_struct
  logical :: use_Leith    ! True if any Leith parameterizations are enabled
  logical :: use_vort_xy  ! True if vort_xy must be computed
  logical :: use_Smag     ! True if a Smagorinsky viscosity is enabled
  integer :: is_vort, ie_vort, js_vort, je_vort  ! Loop ranges for vorticity terms
  integer :: is_Kh, ie_Kh, js_Kh, je_Kh  ! Loop ranges for thickness point viscosities
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, kstart, kend, kmax
  real :: inv_PI3, inv_PI2, inv_PI6 ! Powers of the inverse of pi [nondim]

  ! Fields evaluated on active layers, used for constructing 3D stress fields
  ! NOTE: The position of these declarations can impact performance, due to the
  !   very large number of stack arrays in this function.  Move with caution!
  ! NOTE: Several of these are declared with the memory extent of q-points, but the
  !   same arrays are also used at h-points to reduce the memory footprint of this
  !   module, so they should never be used in halo point or checksum calls.
  type(RealArray_t) :: &
    Ah_a, &           ! biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
    Kh_a, &           ! Laplacian  viscosity (h or q) [L2 T-1 ~> m2 s-1]
    Shear_mag_a, &    ! magnitude of the shear (h or q) [T-1 ~> s-1]
    vert_vort_mag_a, &  ! magnitude of the vertical vorticity gradient (h or q) [L-1 T-1 ~> m-1 s-1]
    vert_vort_mag_smooth_a, &  ! magnitude of gradient of smoothed vertical vorticity (h or q) [L-1 T-1 ~> m-1 s-1]
    hrat_min_a, &     ! h_min divided by the thickness at the stress point (h or q) [nondim]
    visc_bound_rem_a  ! fraction of overall viscous bounds that remain to be applied (h or q) [nondim]
  type(RealArray_t) :: &
    vort_xy_dx_a, &        ! x-derivative of vertical vorticity (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
    vort_xy_dx_smooth_a, & ! x-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
    div_xx_dy_a            ! y-derivative of horizontal divergence (d/dy(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t) :: &
    vort_xy_dy_a, &        ! y-derivative of vertical vorticity (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
    vort_xy_dy_smooth_a, & ! y-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
    div_xx_dx_a            ! x-derivative of horizontal divergence (d/dx(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t) :: &
    Ah_sq_a  ! The square of the biharmonic viscosity [L8 T-2 ~> m8 s-2]
  type(Box_t) :: bxH  ! The h-point iteration box of one tile, [isc:iec, jsc:jec, kstart:kend]
  type(Box_t) :: bxQ  ! The B-grid-index iteration box of one tile, [IscB:IecB, JscB:JecB, kstart:kend]

  nullify(hv_cont, hu_cont)
  if (hv_cont_a%associated()) call hv_cont_a%view(hv_cont)
  if (hu_cont_a%associated()) call hu_cont_a%view(hu_cont)
  call BS_coeff_q_a%view(BS_coeff_q) ; call BS_coeff_h_a%view(BS_coeff_h)
  call visc_limit_q_frac_a%view(visc_limit_q_frac)
  call visc_limit_h_frac_a%view(visc_limit_h_frac)
  call visc_limit_q_flag_a%view(visc_limit_q_flag)
  call visc_limit_h_flag_a%view(visc_limit_h_flag) ; call visc_limit_q_a%view(visc_limit_q)
  call visc_limit_h_a%view(visc_limit_h) ; call sh_xy_q_a%view(sh_xy_q)
  call vort_xy_q_a%view(vort_xy_q) ; call sh_xx_h_a%view(sh_xx_h)
  call div_xx_h_a%view(div_xx_h) ; call grid_Re_Ah_a%view(grid_Re_Ah)
  call grid_Re_Kh_a%view(grid_Re_Kh) ; call Kh_q_a%view(Kh_q) ; call Kh_h_a%view(Kh_h)
  call Ah_q_a%view(Ah_q) ; call Ah_h_a%view(Ah_h) ; call FrictWork_GME_a%view(FrictWork_GME)
  call FrictWork_bh_a%view(FrictWork_bh) ; call FrictWork_a%view(FrictWork)
  call ShSt_a%view(ShSt) ; call NoSt_a%view(NoSt) ; call str_xy_GME_a%view(str_xy_GME)
  call str_xx_GME_a%view(str_xx_GME) ; call GME_coeff_q_a%view(GME_coeff_q)
  call GME_coeff_h_a%view(GME_coeff_h) ; call KH_v_GME_a%view(KH_v_GME)
  call KH_u_GME_a%view(KH_u_GME) ; call GME_effic_q_a%view(GME_effic_q)
  call GME_effic_h_a%view(GME_effic_h) ; call sh_xy_bt_a%view(sh_xy_bt)
  call sh_xx_bt_a%view(sh_xx_bt) ; call slope_y_a%view(slope_y) ; call slope_x_a%view(slope_x)
  call dz_a%view(dz) ; call v_smooth_a%view(v_smooth) ; call u_smooth_a%view(u_smooth)
  call diffv_a%view(diffv) ; call diffu_a%view(diffu) ; call vh_a%view(vh) ; call uh_a%view(uh)
  call h_a%view(h) ; call v_a%view(v) ; call u_a%view(u)

  use_Leith = CS%Leith_Kh .or. CS%Leith_Ah .or. CS%use_Leithy
  use_vort_xy = use_Leith .or. CS%do_vort_xy_q .or. CS%use_ZB2020
  use_Smag = CS%Smagorinsky_Kh .or. CS%Smagorinsky_Ah

  is  = bxH0%idxS(1) ; ie  = bxH0%idxE(1) ; js  = bxH0%idxS(2) ; je  = bxH0%idxE(2) ; nz = bxH0%idxE(3)
  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2)

  h_neglect  = H_subroundoff
  !h_neglect3 = h_neglect**3
  h_neglect3 = h_neglect*h_neglect*h_neglect
  inv_PI3 = 1.0/((4.0*atan(1.0))**3)
  inv_PI2 = 1.0/((4.0*atan(1.0))**2)
  inv_PI6 = inv_PI3 * inv_PI3

  if (CS%EY24_EBT_BS) then
    visc_limit_h(:,:,:) = 0.
    visc_limit_q(:,:,:) = 0.
    visc_limit_h_flag(:,:,:) = 0.
    visc_limit_q_flag(:,:,:) = 0.
    visc_limit_h_frac(:,:,:) = 0.
    visc_limit_q_frac(:,:,:) = 0.
  endif

  skeb_use_frict = STOCH_a%skeb_use_frict

  call h_u_a%alloc(lb=[u_a%lb(1),h_a%lb(2),1], ub=[u_a%ub(1),h_a%ub(2),nkblock])
  call Del2u_a%alloc(lb=[u_a%lb(1),h_a%lb(2),1], ub=[u_a%ub(1),h_a%ub(2),nkblock])
  call h_v_a%alloc(lb=[h_a%lb(1),v_a%lb(2),1], ub=[h_a%ub(1),v_a%ub(2),nkblock])
  call Del2v_a%alloc(lb=[h_a%lb(1),v_a%lb(2),1], ub=[h_a%ub(1),v_a%ub(2),nkblock])
  call sh_xx_smooth_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),nkblock])
  call sh_xx_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),nkblock])
  call dvdy_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),nkblock])
  call dudx_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),nkblock])
  call grad_div_mag_h_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),nkblock])
  call grad_vort_mag_h_2d_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),nkblock])
  call grad_vort_mag_h_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),nkblock])
  call bhstr_xx_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),nkblock])
  call str_xx_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),nkblock])
  call m_leithy_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),nkblock])
  call dDel2udy_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call dDel2vdx_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call dudy_smooth_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call dvdx_smooth_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call dudy_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call dvdx_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call sh_xy_smooth_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call sh_xy_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call hq_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call grad_div_mag_q_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call Del2vort_q_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call grad_vort_mag_q_2d_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call grad_vort_mag_q_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call vort_xy_smooth_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call vort_xy_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call bhstr_xy_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call str_xy_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call visc_bound_rem_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call hrat_min_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call vert_vort_mag_smooth_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call vert_vort_mag_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call Shear_mag_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call Kh_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call Ah_a%alloc(lb=[u_a%lb(1),v_a%lb(2),1], ub=[u_a%ub(1),v_a%ub(2),nkblock])
  call div_xx_dy_a%alloc(lb=[h_a%lb(1),v_a%lb(2),1], ub=[h_a%ub(1),v_a%ub(2),nkblock])
  call vort_xy_dx_smooth_a%alloc(lb=[h_a%lb(1),v_a%lb(2),1], ub=[h_a%ub(1),v_a%ub(2),nkblock])
  call vort_xy_dx_a%alloc(lb=[h_a%lb(1),v_a%lb(2),1], ub=[h_a%ub(1),v_a%ub(2),nkblock])
  call div_xx_dx_a%alloc(lb=[u_a%lb(1),h_a%lb(2),1], ub=[u_a%ub(1),h_a%ub(2),nkblock])
  call vort_xy_dy_smooth_a%alloc(lb=[u_a%lb(1),h_a%lb(2),1], ub=[u_a%ub(1),h_a%ub(2),nkblock])
  call vort_xy_dy_a%alloc(lb=[u_a%lb(1),h_a%lb(2),1], ub=[u_a%ub(1),h_a%ub(2),nkblock])
  call Ah_sq_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),nkblock])
  call h_u_a%view(h_u) ; call Del2u_a%view(Del2u) ; call h_v_a%view(h_v)
  call Del2v_a%view(Del2v) ; call sh_xx_smooth_a%view(sh_xx_smooth) ; call sh_xx_a%view(sh_xx)
  call dvdy_a%view(dvdy) ; call dudx_a%view(dudx) ; call grad_div_mag_h_a%view(grad_div_mag_h)
  call grad_vort_mag_h_2d_a%view(grad_vort_mag_h_2d)
  call grad_vort_mag_h_a%view(grad_vort_mag_h) ; call bhstr_xx_a%view(bhstr_xx)
  call str_xx_a%view(str_xx) ; call m_leithy_a%view(m_leithy) ; call dDel2udy_a%view(dDel2udy)
  call dDel2vdx_a%view(dDel2vdx) ; call dudy_smooth_a%view(dudy_smooth)
  call dvdx_smooth_a%view(dvdx_smooth) ; call dudy_a%view(dudy) ; call dvdx_a%view(dvdx)
  call sh_xy_smooth_a%view(sh_xy_smooth) ; call sh_xy_a%view(sh_xy) ; call hq_a%view(hq)
  call grad_div_mag_q_a%view(grad_div_mag_q) ; call Del2vort_q_a%view(Del2vort_q)
  call grad_vort_mag_q_2d_a%view(grad_vort_mag_q_2d)
  call grad_vort_mag_q_a%view(grad_vort_mag_q) ; call vort_xy_smooth_a%view(vort_xy_smooth)
  call vort_xy_a%view(vort_xy) ; call bhstr_xy_a%view(bhstr_xy) ; call str_xy_a%view(str_xy)
  call visc_bound_rem_a%view(visc_bound_rem) ; call hrat_min_a%view(hrat_min)
  call vert_vort_mag_smooth_a%view(vert_vort_mag_smooth)
  call vert_vort_mag_a%view(vert_vort_mag) ; call Shear_mag_a%view(Shear_mag)
  call Kh_a%view(Kh) ; call Ah_a%view(Ah) ; call div_xx_dy_a%view(div_xx_dy)
  call vort_xy_dx_smooth_a%view(vort_xy_dx_smooth) ; call vort_xy_dx_a%view(vort_xy_dx)
  call div_xx_dx_a%view(div_xx_dx) ; call vort_xy_dy_smooth_a%view(vort_xy_dy_smooth)
  call vort_xy_dy_a%view(vort_xy_dy) ; call Ah_sq_a%view(Ah_sq)

  if (CS%use_Leithy) m_leithy(:,:,:) = 0.0 ! Initialize

  if (present(OBC)) then ; if (associated(OBC)) then ; if (OBC%OBC_pe) then
    apply_OBC = OBC%Flather_u_BCs_exist_globally .or. OBC%Flather_v_BCs_exist_globally
    apply_OBC = .true.
  endif ; endif ; endif

  apply_OBC_strain = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    apply_OBC_strain = (OBC%strain_config /= OBC_STRAIN_NONE) &
      .and. ((.not. CS%OBC_strain_bug) .or. (OBC%strain_config /= OBC_STRAIN_SPECIFIED))
  endif ; endif


  find_FrictWork = CS%do_FrictWork .or. CS%do_FrictWorkIntz &
      .or. MEKE_a%mom_src_a%associated()
  find_FrictWork_bh = CS%do_FrictWork_bh .or. CS%do_FrictWorkIntz_bh &
      .or. MEKE_a%mom_src_bh_a%associated()

  use_kh_struct = VarMix_a%BS_struct_a%associated()
  backscat_subround = 0.0
  if (find_FrictWork .and. MEKE_a%mom_src_a%associated() .and. (MEKE_a%backscatter_Ro_c > 0.0) .and. &
      (MEKE_a%backscatter_Ro_Pow /= 0.0)) &
    backscat_subround = (1.0e-16/MEKE_a%backscatter_Ro_c)**(1.0/MEKE_a%backscatter_Ro_Pow)

  ! Toggle whether to use a Laplacian viscosity derived from MEKE
  use_MEKE_Ku = MEKE_a%Ku_a%associated()
  use_MEKE_Au = MEKE_a%Au_a%associated()

  use_cont_huv = CS%use_cont_thick .and. hu_cont_a%associated() .and. hv_cont_a%associated()

  rescale_Kh = .false.
  if (VarMix_a%use_variable_mixing) then
    rescale_Kh = VarMix_a%Resoln_scaled_Kh
  endif

  ! Set the halo sizes used for the thickness-point viscosities.
  if (CS%use_Leithy .or. CS%debug) then
    js_Kh = js-1 ; je_Kh = je+1 ; is_Kh = is-1 ; ie_Kh = ie+1
  else
    js_Kh = Jsq ; je_Kh = je+1 ; is_Kh = Isq ; ie_Kh = ie+1
  endif

  ! Set the halo sizes used for the vorticity calculations.
  if ((CS%Leith_Kh) .or. (CS%Leith_Ah) .or. (CS%use_Leithy)) then
    js_vort = js_Kh-2 ; je_vort = Jeq+2 ; is_vort = is_Kh-2 ; ie_vort = Ieq+2
  else
    js_vort = js-2 ; je_vort = Jeq+1 ; is_vort = is-2 ; ie_vort = Ieq+1
  endif

  nullify(MEKE_mom_src, MEKE_mom_src_bh, MEKE_GME_snk)
  if (MEKE_a%mom_src_a%associated()) call MEKE_a%mom_src_a%view(MEKE_mom_src)
  if (MEKE_a%mom_src_bh_a%associated()) call MEKE_a%mom_src_bh_a%view(MEKE_mom_src_bh)
  if (MEKE_a%GME_snk_a%associated()) call MEKE_a%GME_snk_a%view(MEKE_GME_snk)

  if (find_FrictWork .and. MEKE_a%mom_src_a%associated()) then
    do j=js,je ; do i=is,ie
      MEKE_mom_src(i,j) = 0.
    enddo ; enddo

    if (MEKE_a%mom_src_bh_a%associated()) then
      do j=js,je ; do i=is,ie
        MEKE_mom_src_bh(i,j) = 0.
      enddo ; enddo
    endif

    if (MEKE_a%GME_snk_a%associated()) then
      do j=js,je ; do i=is,ie
        MEKE_GME_snk(i,j) = 0.
      enddo ; enddo
    endif
  endif

  use_Leith = CS%Leith_Kh .or. CS%Leith_Ah .or. CS%use_Leithy
  use_vort_xy = use_Leith .or. CS%do_vort_xy_q .or. CS%use_ZB2020
  use_Smag = CS%Smagorinsky_Kh .or. CS%Smagorinsky_Ah

  ! TODO: Only needed if FrictWork_bh is true, and currently only used on CPU,
  !   but I do not yet see any benefit to breaking up the calculation.


  !!$omp target enter data map(alloc: CS%use_Leithy)

  call bxH%safe_alloc(ndims=3) ; call bxQ%safe_alloc(ndims=3)
  do kstart=1,nz,nkblock
    kend = min(kstart+nkblock-1,nz)
    kmax = kend-kstart+1
    call bxH%set(idxS=[is,js,kstart], idxE=[ie,je,kend])
    call bxQ%set(idxS=[Isq,Jsq,kstart], idxE=[Ieq,Jeq,kend])

    call hor_visc_strain(Gcore, CS, metrics, OBC, bxH, bxQ, nkblock, use_vort_xy, apply_OBC, &
                         apply_OBC_strain, use_cont_huv, is_Kh, ie_Kh, js_Kh, je_Kh, is_vort, &
                         ie_vort, js_vort, je_vort, u_a, v_a, h_a, hu_cont_a, hv_cont_a, &
                         u_smooth_a, v_smooth_a, NoSt_a, ShSt_a, Del2u_a, h_u_a, Del2v_a, h_v_a, &
                         sh_xx_a, sh_xx_smooth_a, dudx_a, dvdy_a, dvdx_a, dudy_a, dvdx_smooth_a, &
                         dudy_smooth_a, sh_xy_a, sh_xy_smooth_a, vort_xy_a, vort_xy_smooth_a)

    if (use_Leith) then
      call hor_visc_Leith_grad_a(Gcore, CS, nkblock, kmax, is, ie, js, je, is_Kh, ie_Kh, js_Kh, &
                                 je_Kh, Ieq, Jeq, dudx_a, dvdy_a, vort_xy_a, vort_xy_smooth_a, &
                                 grad_vort_mag_h_2d_a, grad_div_mag_h_a, grad_vort_mag_q_2d_a, &
                                 grad_div_mag_q_a, Del2vort_q_a, vort_xy_dx_a, vort_xy_dx_smooth_a, &
                                 div_xx_dy_a, vort_xy_dy_a, vort_xy_dy_smooth_a, div_xx_dx_a)
      if (CS%use_QG_Leith_visc) &
        call hor_visc_QG_Leith_host(host, nkblock, kstart, kmax, h_a, dz_a, slope_x_a, slope_y_a, &
                                    div_xx_dx_a, div_xx_dy_a, vort_xy_dx_a, vort_xy_dy_a)
      call hor_visc_Leith_grad_b(CS, nkblock, kmax, is, js, is_Kh, ie_Kh, js_Kh, je_Kh, Ieq, Jeq, &
                                 grad_vort_mag_h_a, grad_vort_mag_q_a, vert_vort_mag_smooth_a, &
                                 vort_xy_dx_a, vort_xy_dx_smooth_a, vort_xy_dy_a, &
                                 vort_xy_dy_smooth_a)
    endif ! CS%Leith_Kh

    call hor_visc_h_coeffs_a(MEKE_a, VarMix_a, CS, metrics, bxH, bxQ, nkblock, h_neglect, inv_PI3, &
                             inv_PI6, use_kh_struct, use_MEKE_Ku, rescale_Kh, is_Kh, ie_Kh, js_Kh, &
                             je_Kh, u_a, v_a, h_a, Kh_h_a, grid_Re_Kh_a, div_xx_h_a, sh_xx_h_a, &
                             h_u_a, h_v_a, sh_xx_a, str_xx_a, grad_vort_mag_h_a, &
                             grad_vort_mag_h_2d_a, grad_div_mag_h_a, dudx_a, dvdy_a, sh_xy_a, &
                             Del2vort_q_a, Ah_a, Kh_a, Shear_mag_a, vert_vort_mag_a, hrat_min_a, &
                             visc_bound_rem_a)

    if (CS%biharmonic) then
      if ((CS%Smagorinsky_Ah) .or. (CS%Leith_Ah) .or. (CS%use_Leithy)) then
        if (CS%use_Leithy) then
          call hor_visc_Leithy_m(Gcore, CS, metrics, nkblock, kmax, is_Kh, ie_Kh, js_Kh, je_Kh, &
                                 inv_PI6, Del2vort_q_a, vert_vort_mag_a, vort_xy_smooth_a, &
                                 m_leithy_a)
          if (CS%smooth_Ah) then
            call hor_visc_Leithy_m_smooth_host(host, Gcore, bxH0, nkblock, kmax, m_leithy_a)
          endif
          call hor_visc_Leithy_Ah(CS, metrics, nkblock, kmax, is_Kh, ie_Kh, js_Kh, je_Kh, inv_PI6, &
                                  Del2vort_q_a, vert_vort_mag_smooth_a, Ah_a, m_leithy_a, Ah_sq_a)
          if (CS%smooth_Ah) then
            call hor_visc_Leithy_Ah_sq_smooth_host(host, Gcore, bxH0, nkblock, kmax, Ah_sq_a)
          endif
          call hor_visc_Leithy_Ah_h(CS, metrics, nkblock, kstart, kmax, is_Kh, ie_Kh, js_Kh, je_Kh, &
                                    Ah_a, Ah_h_a, Ah_sq_a)
        endif
      endif ! Smagorinsky_Ah or Leith_Ah or Leith+E
    endif

    call hor_visc_h_coeffs_b(MEKE_a, VarMix_a, Gcore, CS, metrics, OBC, bxH, bxQ, nkblock, &
                             h_neglect, h_neglect3, apply_OBC, use_kh_struct, use_MEKE_Au, is_Kh, &
                             ie_Kh, js_Kh, je_Kh, u_a, v_a, h_a, Ah_h_a, Kh_h_a, grid_Re_Ah_a, &
                             visc_limit_h_a, visc_limit_h_flag_a, visc_limit_h_frac_a, &
                             BS_coeff_h_a, Del2u_a, h_u_a, Del2v_a, h_v_a, sh_xx_a, sh_xx_smooth_a, &
                             str_xx_a, bhstr_xx_a, m_leithy_a, dDel2vdx_a, dDel2udy_a, sh_xy_a, &
                             hq_a, Ah_a, Kh_a, Shear_mag_a, hrat_min_a, visc_bound_rem_a)

    call hor_visc_q_coeffs(MEKE_a, VarMix_a, Gcore, CS, metrics, bxH, bxQ, nkblock, inv_PI3, &
                           inv_PI6, find_FrictWork, use_kh_struct, use_MEKE_Ku, use_MEKE_Au, &
                           rescale_Kh, u_a, v_a, sh_xx_bt_a, sh_xy_bt_a, GME_effic_h_a, &
                           GME_effic_q_a, KH_u_GME_a, KH_v_GME_a, GME_coeff_h_a, GME_coeff_q_a, &
                           str_xx_GME_a, str_xy_GME_a, Ah_h_a, Ah_q_a, Kh_h_a, Kh_q_a, vort_xy_q_a, &
                           sh_xy_q_a, visc_limit_q_a, visc_limit_q_flag_a, visc_limit_q_frac_a, &
                           BS_coeff_q_a, sh_xx_a, str_xx_a, dDel2vdx_a, dDel2udy_a, sh_xy_a, &
                           sh_xy_smooth_a, str_xy_a, bhstr_xy_a, vort_xy_a, grad_vort_mag_q_a, &
                           grad_vort_mag_q_2d_a, Del2vort_q_a, grad_div_mag_q_a, hq_a, Ah_a, Kh_a, &
                           Shear_mag_a, vert_vort_mag_a, hrat_min_a, visc_bound_rem_a)

    if (CS%use_GME) then
      ! Applying GME diagonal term.  This is linear and the arguments can be rescaled.
      ! smooth_GME is an impure host routine, so str_xx_GME/str_xy_GME must be brought back
      ! from the device before this loop and pushed back afterward.
      call hor_visc_GME_smooth_host(host, Gcore, bxH0, bxQ0, nkblock, kmax, str_xx_GME_a, &
                                    str_xy_GME_a)
    endif

    call hor_visc_stress_div(MEKE_a, Gcore, CS, metrics, OBC, STOCH_a, bxH, bxQ, nkblock, &
                             h_neglect, apply_OBC, find_FrictWork, find_FrictWork_bh, &
                             backscat_subround, skeb_use_frict, u_a, v_a, h_a, uh_a, vh_a, diffu_a, &
                             diffv_a, str_xx_GME_a, str_xy_GME_a, FrictWork_a, FrictWork_bh_a, &
                             FrictWork_GME_a, visc_limit_h_flag_a, h_u_a, h_v_a, sh_xx_a, str_xx_a, &
                             bhstr_xx_a, sh_xy_a, str_xy_a, bhstr_xy_a, hq_a, H_to_RZ, s_to_T)

    ! Pass the velocity gradients and thickness to ZB2020
    if (CS%use_ZB2020) then
      call hor_visc_ZB2020_copy_host(host, nkblock, kstart, kmax, sh_xx_a, sh_xy_a, vort_xy_a, hq_a)
    endif
  enddo ! end of kstart block loop
  call bxH%free() ; call bxQ%free()
  call h_u_a%free() ; call Del2u_a%free() ; call h_v_a%free() ; call Del2v_a%free()
  call sh_xx_smooth_a%free() ; call sh_xx_a%free() ; call dvdy_a%free() ; call dudx_a%free()
  call grad_div_mag_h_a%free() ; call grad_vort_mag_h_2d_a%free() ; call grad_vort_mag_h_a%free()
  call bhstr_xx_a%free() ; call str_xx_a%free() ; call m_leithy_a%free() ; call dDel2udy_a%free()
  call dDel2vdx_a%free() ; call dudy_smooth_a%free() ; call dvdx_smooth_a%free()
  call dudy_a%free() ; call dvdx_a%free() ; call sh_xy_smooth_a%free() ; call sh_xy_a%free()
  call hq_a%free() ; call grad_div_mag_q_a%free() ; call Del2vort_q_a%free()
  call grad_vort_mag_q_2d_a%free() ; call grad_vort_mag_q_a%free() ; call vort_xy_smooth_a%free()
  call vort_xy_a%free() ; call bhstr_xy_a%free() ; call str_xy_a%free()
  call visc_bound_rem_a%free() ; call hrat_min_a%free() ; call vert_vort_mag_smooth_a%free()
  call vert_vort_mag_a%free() ; call Shear_mag_a%free() ; call Kh_a%free() ; call Ah_a%free()
  call div_xx_dy_a%free() ; call vort_xy_dx_smooth_a%free() ; call vort_xy_dx_a%free()
  call div_xx_dx_a%free() ; call vort_xy_dy_smooth_a%free() ; call vort_xy_dy_a%free()
  call Ah_sq_a%free()







end subroutine horizontal_viscosity_TR

!> Calculates the horizontal tension and shearing strain, the thicknesses at velocity points,
!! the Laplacian of velocity, and the vorticity for one tile of horizontal_viscosity.
subroutine hor_visc_strain(Gcore, CS, metrics, OBC, bxH, bxQ, nkblock, use_vort_xy, apply_OBC, &
                           apply_OBC_strain, use_cont_huv, is_Kh, ie_Kh, js_Kh, je_Kh, is_vort, &
                           ie_vort, js_vort, je_vort, u_a, v_a, h_a, hu_cont_a, hv_cont_a, &
                           u_smooth_a, v_smooth_a, NoSt_a, ShSt_a, Del2u_a, h_u_a, Del2v_a, h_v_a, &
                           sh_xx_a, sh_xx_smooth_a, dudx_a, dvdy_a, dvdx_a, dudy_a, dvdx_smooth_a, &
                           dudy_smooth_a, sh_xy_a, sh_xy_smooth_a, vort_xy_a, vort_xy_smooth_a)
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  type(hv_metrics_CS), intent(in)    :: metrics       !< Persistent copies of the horizontal
                                                     !! viscosity metric and coefficient arrays
  type(ocean_OBC_type), optional, pointer      :: OBC  !< Pointer to an open boundary condition type
  type(Box_t),                   intent(in)    :: bxH  !< The h-point iteration box of this tile,
                                                       !! [isc:iec, jsc:jec, kstart:kend]
  type(Box_t),                   intent(in)    :: bxQ  !< The B-grid-index iteration box of this
                                                       !! tile, [IscB:IecB, JscB:JecB, kstart:kend]
  integer,                       intent(in)    :: nkblock !< The effective k-block size [nondim]
  logical, intent(in)    :: use_vort_xy !< True if vort_xy must be computed
  logical, intent(in)    :: apply_OBC !< If true, open boundary conditions are applied on this PE
  logical, intent(in)    :: apply_OBC_strain !< If true, the specified open boundary strain configuration is applied
  logical, intent(in)    :: use_cont_huv !< If true, use the continuity solver thicknesses at velocity points
  integer, intent(in)    :: is_Kh !< Loop ranges for thickness point viscosities
  integer, intent(in)    :: ie_Kh !< Loop ranges for thickness point viscosities
  integer, intent(in)    :: js_Kh !< Loop ranges for thickness point viscosities
  integer, intent(in)    :: je_Kh !< Loop ranges for thickness point viscosities
  integer, intent(in)    :: is_vort !< Loop ranges for vorticity terms
  integer, intent(in)    :: ie_vort !< Loop ranges for vorticity terms
  integer, intent(in)    :: js_vort !< Loop ranges for vorticity terms
  integer, intent(in)    :: je_vort !< Loop ranges for vorticity terms
  type(RealArray_t), intent(in)    :: u_a !< The zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: v_a !< The meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: h_a !< Layer thicknesses [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: hu_cont_a !< Layer thickness at u-points [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: hv_cont_a !< Layer thickness at v-points [H ~> m or kg m-2].
  type(RealArray_t), intent(in)    :: u_smooth_a !< Zonal velocity, smoothed with a spatial low-pass filter
                                                 !! [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: v_smooth_a !< Meridional velocity, smoothed with a spatial low-pass
                                                 !! filter [L T-1 ~> m s-1]
  type(RealArray_t), intent(inout) :: NoSt_a !< A diagnostic array of normal stress [T-1 ~> s-1].
  type(RealArray_t), intent(inout) :: ShSt_a !< A diagnostic array of shear stress [T-1 ~> s-1].
  type(RealArray_t), intent(inout) :: Del2u_a !< The u-component of the Laplacian of velocity
                                              !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: h_u_a !< Thickness interpolated to u points [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: Del2v_a !< The v-component of the Laplacian of velocity
                                              !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: h_v_a !< Thickness interpolated to v points [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: sh_xx_a !< horizontal tension (du/dx - dv/dy) including metric terms
                                              !! [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xx_smooth_a !< horizontal tension from smoothed velocity including
                                                     !! metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: dudx_a !< components in the horizontal tension [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: dvdy_a !< components in the horizontal tension [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: dvdx_a !< components in the shearing strain [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: dudy_a !< components in the shearing strain [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: dvdx_smooth_a !< components in the shearing strain from smoothed
                                                    !! velocity [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: dudy_smooth_a !< components in the shearing strain from smoothed
                                                    !! velocity [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xy_a !< horizontal shearing strain (du/dy + dv/dx) including metric
                                              !! terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xy_smooth_a !< horizontal shearing strain from smoothed velocity
                                                     !! including metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: vort_xy_a !< Vertical vorticity (dv/dx - du/dy) including metric terms
                                                !! [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: vort_xy_smooth_a !< Vertical vorticity including metric terms, smoothed
                                                       !! [T-1 ~> s-1]

  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: vort_xy_smooth, vort_xy, sh_xy_smooth, sh_xy, &
    dudy_smooth, dvdx_smooth, dudy, dvdx, dvdy, dudx, sh_xx_smooth, sh_xx, h_v, Del2v, h_u, Del2u, &
    ShSt, NoSt, v_smooth, u_smooth, hv_cont, hu_cont, h, v, u
  real, dimension(:,:), contiguous, pointer :: IdyCu, IdxCv, IdyCv, IdxCu, mask2dT, dxBu, &
    mask2dBu, IareaBu, dyCv, dxCu
  real, dimension(:,:), contiguous, pointer :: dx2h, dy2h, dx_dyT, dy_dxT, dx2q, dy2q, dx_dyBu, dy_dxBu, &
    Idx2dyCu, Idxdy2u, Idx2dyCv, Idxdy2v
  real :: dudx_smooth, dvdy_smooth ! components in the horizontal tension from smoothed velocity [T-1 ~> s-1]
  integer :: i, j, k, kk, n, is, ie, js, je, Isq, Ieq, Jsq, Jeq, kstart, kmax

  nullify(hv_cont, hu_cont)
  if (hv_cont_a%associated()) call hv_cont_a%view(hv_cont)
  if (hu_cont_a%associated()) call hu_cont_a%view(hu_cont)
  call vort_xy_smooth_a%view(vort_xy_smooth) ; call vort_xy_a%view(vort_xy)
  call sh_xy_smooth_a%view(sh_xy_smooth) ; call sh_xy_a%view(sh_xy)
  call dudy_smooth_a%view(dudy_smooth) ; call dvdx_smooth_a%view(dvdx_smooth)
  call dudy_a%view(dudy) ; call dvdx_a%view(dvdx) ; call dvdy_a%view(dvdy)
  call dudx_a%view(dudx) ; call sh_xx_smooth_a%view(sh_xx_smooth) ; call sh_xx_a%view(sh_xx)
  call h_v_a%view(h_v) ; call Del2v_a%view(Del2v) ; call h_u_a%view(h_u)
  call Del2u_a%view(Del2u) ; call ShSt_a%view(ShSt) ; call NoSt_a%view(NoSt)
  call v_smooth_a%view(v_smooth) ; call u_smooth_a%view(u_smooth) ; call h_a%view(h)
  call v_a%view(v) ; call u_a%view(u)
  call Gcore%IdyCu%view(IdyCu) ; call Gcore%IdxCv%view(IdxCv) ; call Gcore%IdyCv%view(IdyCv)
  call Gcore%IdxCu%view(IdxCu) ; call Gcore%mask2dT%view(mask2dT) ; call Gcore%dxBu%view(dxBu)
  call Gcore%mask2dBu%view(mask2dBu) ; call Gcore%IareaBu%view(IareaBu)
  call Gcore%dyCv%view(dyCv) ; call Gcore%dxCu%view(dxCu)

  nullify(dx2h, dy2h, dx_dyT, dy_dxT, dx2q, dy2q, dx_dyBu, dy_dxBu, Idx2dyCu, Idxdy2u, Idx2dyCv, Idxdy2v)
  if (metrics%dx2h%associated()) call metrics%dx2h%view(dx2h)
  if (metrics%dy2h%associated()) call metrics%dy2h%view(dy2h)
  if (metrics%dx_dyT%associated()) call metrics%dx_dyT%view(dx_dyT)
  if (metrics%dy_dxT%associated()) call metrics%dy_dxT%view(dy_dxT)
  if (metrics%dx2q%associated()) call metrics%dx2q%view(dx2q)
  if (metrics%dy2q%associated()) call metrics%dy2q%view(dy2q)
  if (metrics%dx_dyBu%associated()) call metrics%dx_dyBu%view(dx_dyBu)
  if (metrics%dy_dxBu%associated()) call metrics%dy_dxBu%view(dy_dxBu)
  if (metrics%Idx2dyCu%associated()) call metrics%Idx2dyCu%view(Idx2dyCu)
  if (metrics%Idxdy2u%associated()) call metrics%Idxdy2u%view(Idxdy2u)
  if (metrics%Idx2dyCv%associated()) call metrics%Idx2dyCv%view(Idx2dyCv)
  if (metrics%Idxdy2v%associated()) call metrics%Idxdy2v%view(Idxdy2v)

  is  = bxH%idxS(1) ; ie  = bxH%idxE(1) ; js  = bxH%idxS(2) ; je  = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  kstart = bxH%idxS(3) ; kmax = bxH%idxE(3) - bxH%idxS(3) + 1

  k=kstart

  ! The following are the forms of the horizontal tension and horizontal
  ! shearing strain advocated by Smagorinsky (1993) and discussed in
  ! Griffies and Hallberg (2000).

  ! Calculate horizontal tension
  do concurrent (kk=1:kmax, j=Jsq-1:Jeq+2, i=Isq-1:Ieq+2) DO_LOCALITY(local(k))
    k = kstart + kk - 1
    dudx(i,j,kk) = DY_dxT(i,j)*((IdyCu(I,j) * u(I,j,k)) - &
                                  (IdyCu(I-1,j) * u(I-1,j,k)))
    dvdy(i,j,kk) = DX_dyT(i,j)*((IdxCv(i,J) * v(i,J,k)) - &
                                  (IdxCv(i,J-1) * v(i,J-1,k)))
    sh_xx(i,j,kk) = dudx(i,j,kk) - dvdy(i,j,kk)
  enddo

  ! Components for the shearing strain
  do concurrent (kk=1:kmax, J=js_vort:je_vort, I=is_vort:ie_vort) DO_LOCALITY(local(k))
    k = kstart + kk - 1
    dvdx(I,J,kk) = DY_dxBu(I,J)*((v(i+1,J,k)*IdyCv(i+1,J)) - (v(i,J,k)*IdyCv(i,J)))
    dudy(I,J,kk) = DX_dyBu(I,J)*((u(I,j+1,k)*IdxCu(I,j+1)) - (u(I,j,k)*IdxCu(I,j)))
  enddo

  if (CS%use_Leithy) then
    ! Calculate horizontal tension from smoothed velocity
    do concurrent (kk=1:kmax, j=Jsq:Jeq+1, i=Isq:Ieq+1) DO_LOCALITY(local(k, dudx_smooth, dvdy_smooth))
      k = kstart + kk - 1
      dudx_smooth = DY_dxT(i,j)*((IdyCu(I,j) * u_smooth(I,j,k)) - &
                                         (IdyCu(I-1,j) * u_smooth(I-1,j,k)))
      dvdy_smooth = DX_dyT(i,j)*((IdxCv(i,J) * v_smooth(i,J,k)) - &
                                         (IdxCv(i,J-1) * v_smooth(i,J-1,k)))
      sh_xx_smooth(i,j,kk) = dudx_smooth - dvdy_smooth
    enddo

    ! Components for the shearing strain from smoothed velocity
    do concurrent (kk=1:kmax, J=js_Kh-1:je_Kh, I=is_Kh-1:ie_Kh) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      dvdx_smooth(I,J,kk) = DY_dxBu(I,J) * &
                       ((v_smooth(i+1,J,k)*IdyCv(i+1,J)) - (v_smooth(i,J,k)*IdyCv(i,J)))
      dudy_smooth(I,J,kk) = DX_dyBu(I,J) * &
                       ((u_smooth(I,j+1,k)*IdxCu(I,j+1)) - (u_smooth(I,j,k)*IdxCu(I,j)))
    enddo
  endif ! use Leith+E

  if (CS%do_normstress) then
    do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      NoSt(i,j,k) = sh_xx(i,j,kk)
    enddo
  endif

  ! Interpolate the thicknesses to velocity points.
  ! The extra wide halos are to accommodate the cross-corner-point projections
  ! in OBCs, which are not ordinarily be necessary, and might not be necessary
  ! even with OBCs if the accelerations are zeroed at OBC points, in which
  ! case the j-loop for h_u could collapse to j=js=1,je+1. -RWH
  if (use_cont_huv) then
    do concurrent (kk=1:kmax, j=js-2:je+2, I=Isq-1:Ieq+1) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      h_u(I,j,kk) = hu_cont(I,j,k)
    enddo
    do concurrent (kk=1:kmax, J=Jsq-1:Jeq+1, i=is-2:ie+2) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      h_v(i,J,kk) = hv_cont(i,J,k)
    enddo
  elseif (CS%use_land_mask) then
    do concurrent (kk=1:kmax, j=js-2:je+2, I=is-2:Ieq+1) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      h_u(I,j,kk) = 0.5 * (mask2dT(i,j)*h(i,j,k) + mask2dT(i+1,j)*h(i+1,j,k))
    enddo
    do concurrent (kk=1:kmax, J=js-2:Jeq+1, i=is-2:ie+2) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      h_v(i,J,kk) = 0.5 * (mask2dT(i,j)*h(i,j,k) + mask2dT(i,j+1)*h(i,j+1,k))
    enddo
  else
    do concurrent (kk=1:kmax, j=js-2:je+2, I=is-2:Ieq+1) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      h_u(I,j,kk) = 0.5 * (h(i,j,k) + h(i+1,j,k))
    enddo
    do concurrent (kk=1:kmax, J=js-2:Jeq+1, i=is-2:ie+2) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      h_v(i,J,kk) = 0.5 * (h(i,j,k) + h(i,j+1,k))
    enddo
  endif

  ! Adjust contributions to shearing strain and interpolated values of
  ! thicknesses on open boundaries.
  if (apply_OBC) then
    ! TODO: Reindent this later
    do n=1,OBC%number_of_segments

    J = OBC%segment(n)%HI%JsdB ; I = OBC%segment(n)%HI%IsdB
    if (apply_OBC_strain) then
      if (OBC%segment(n)%is_N_or_S .and. (J >= Js_vort) .and. (J <= Je_vort)) then
        do concurrent (kk=1:kmax, &
                       I=max(OBC%segment(n)%HI%IsdB,Is_vort):min(OBC%segment(n)%HI%IedB,Ie_vort)) &
                       DO_LOCALITY(local(k))
          k = kstart + kk - 1
          select case (OBC%strain_config)
            case (OBC_STRAIN_ZERO)
              dvdx(I,J,kk) = 0. ; dudy(I,J,kk) = 0.
            case (OBC_STRAIN_FREESLIP)
              dudy(I,J,kk) = 0.
            case (OBC_STRAIN_COMPUTED)
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                dudy(I,J,kk) = 2.0*DX_dyBu(I,J)* &
                            (OBC%segment(n)%tangential_vel(I,J,k) - u(I,j,k))*IdxCu(I,j)
              else
                dudy(I,J,kk) = 2.0*DX_dyBu(I,J)* &
                            (u(I,j+1,k) - OBC%segment(n)%tangential_vel(I,J,k))*IdxCu(I,j+1)
              endif
            case (OBC_STRAIN_SPECIFIED)
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                dudy(I,J,kk) = DX_dyBu(I,J)*OBC%segment(n)%tangential_grad(I,J,k)*IdxCu(I,j)*dxBu(I,J)
              else
                dudy(I,J,kk) = DX_dyBu(I,J)*OBC%segment(n)%tangential_grad(I,J,k)*IdxCu(I,j+1)*dxBu(I,J)
              endif
          end select
          if (CS%use_Leithy) then
            dvdx_smooth(I,J,kk) = dvdx(I,J,kk)
            dudy_smooth(I,J,kk) = dudy(I,J,kk)
          endif
        enddo
      elseif (OBC%segment(n)%is_E_or_W .and. (I >= is_vort) .and. (I <= ie_vort)) then
        do concurrent (kk=1:kmax, &
                       J=max(OBC%segment(n)%HI%JsdB,js_vort):min(OBC%segment(n)%HI%JedB,je_vort)) &
                       DO_LOCALITY(local(k))
          k = kstart + kk - 1
          select case (OBC%strain_config)
            case (OBC_STRAIN_ZERO)
              dvdx(I,J,kk) = 0. ; dudy(I,J,kk) = 0.
            case (OBC_STRAIN_FREESLIP)
              dvdx(I,J,kk) = 0.
            case (OBC_STRAIN_COMPUTED)
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                dvdx(I,J,kk) = 2.0*DY_dxBu(I,J)* &
                            (OBC%segment(n)%tangential_vel(I,J,k) - v(i,J,k))*IdyCv(i,J)
              else
                dvdx(I,J,kk) = 2.0*DY_dxBu(I,J)* &
                            (v(i+1,J,k) - OBC%segment(n)%tangential_vel(I,J,k))*IdyCv(i+1,J)
              endif
            case (OBC_STRAIN_SPECIFIED)
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                dvdx(I,J,kk) = DY_dxBu(I,J)*OBC%segment(n)%tangential_grad(I,J,k)*IdyCv(i,J)*dxBu(I,J)
              else
                dvdx(I,J,kk) = DY_dxBu(I,J)*OBC%segment(n)%tangential_grad(I,J,k)*IdyCv(i+1,J)*dxBu(I,J)
              endif
          end select
          if (CS%use_Leithy) then
            dvdx_smooth(I,J,kk) = dvdx(I,J,kk)
            dudy_smooth(I,J,kk) = dudy(I,J,kk)
          endif
        enddo
      endif
    endif

    if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
      ! There are extra wide halos here to accommodate the cross-corner-point
      ! OBC projections, but they might not be necessary if the accelerations
      ! are always zeroed out at OBC points, in which case the i-loop below
      ! becomes do i=is-1,ie+1. -RWH
      if ((J >= js-2) .and. (J <= Jeq+1)) then
        do concurrent (kk=1:kmax, &
                       i=max(is-2,OBC%segment(n)%HI%isd):min(ie+2,OBC%segment(n)%HI%ied)) DO_LOCALITY(local(k))
          k = kstart + kk - 1
          h_v(i,J,kk) = h(i,j,k)
        enddo
      endif
    elseif (OBC%segment(n)%direction == OBC_DIRECTION_S) then
      if ((J >= js-2) .and. (J <= Jeq+1)) then
        do concurrent (kk=1:kmax, &
                       i=max(is-2,OBC%segment(n)%HI%isd):min(ie+2,OBC%segment(n)%HI%ied)) DO_LOCALITY(local(k))
          k = kstart + kk - 1
          h_v(i,J,kk) = h(i,j+1,k)
        enddo
      endif
    elseif (OBC%segment(n)%direction == OBC_DIRECTION_E) then
      if ((I >= is-2) .and. (I <= Ieq+1)) then
        do concurrent (kk=1:kmax, &
                       j=max(js-2,OBC%segment(n)%HI%jsd):min(je+2,OBC%segment(n)%HI%jed)) DO_LOCALITY(local(k))
          k = kstart + kk - 1
          h_u(I,j,kk) = h(i,j,k)
        enddo
      endif
    elseif (OBC%segment(n)%direction == OBC_DIRECTION_W) then
      if ((I >= is-2) .and. (I <= Ieq+1)) then
        do concurrent (kk=1:kmax, &
                       j=max(js-2,OBC%segment(n)%HI%jsd):min(je+2,OBC%segment(n)%HI%jed)) DO_LOCALITY(local(k))
          k = kstart + kk - 1
          h_u(I,j,kk) = h(i+1,j,k)
        enddo
      endif
    endif
  enddo ; endif
  ! Now project thicknesses across corner points on OBCs.
  if (apply_OBC) then ; do n=1,OBC%number_of_segments
    J = OBC%segment(n)%HI%JsdB ; I = OBC%segment(n)%HI%IsdB
    if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
      if ((J >= js-2) .and. (J <= je)) then
        do concurrent (kk=1:kmax, I=max(is-2,OBC%segment(n)%HI%IsdB):min(Ieq+1,OBC%segment(n)%HI%IedB))
          h_u(I,j+1,kk) = h_u(I,j,kk)
        enddo
      endif
    elseif (OBC%segment(n)%direction == OBC_DIRECTION_S) then
      if ((J >= js-1) .and. (J <= je+1)) then
        do concurrent (kk=1:kmax, I=max(is-2,OBC%segment(n)%HI%isd):min(Ieq+1,OBC%segment(n)%HI%ied))
          h_u(I,j,kk) = h_u(I,j+1,kk)
        enddo
      endif
    elseif (OBC%segment(n)%direction == OBC_DIRECTION_E) then
      if ((I >= is-2) .and. (I <= ie)) then
        do concurrent (kk=1:kmax, J=max(js-2,OBC%segment(n)%HI%jsd):min(Jeq+1,OBC%segment(n)%HI%jed))
          h_v(i+1,J,kk) = h_v(i,J,kk)
        enddo
      endif
    elseif (OBC%segment(n)%direction == OBC_DIRECTION_W) then
      if ((I >= is-1) .and. (I <= ie+1)) then
        do concurrent (kk=1:kmax, J=max(js-2,OBC%segment(n)%HI%jsd):min(Jeq+1,OBC%segment(n)%HI%jed))
          h_v(i,J,kk) = h_v(i+1,J,kk)
        enddo
      endif
    endif
  enddo
  ! TODO: Fix indentation
  endif

  ! Shearing strain (including no-slip boundary conditions at the 2-D land-sea mask).
  ! dudy and dvdx include modifications at OBCs from above.
  if (CS%no_slip) then
    do concurrent (kk=1:kmax, J=js-2:Jeq+1, I=is-2:Ieq+1) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      sh_xy(I,J,kk) = (2.0-mask2dBu(I,J)) * ( dvdx(I,J,kk) + dudy(I,J,kk) )
      if (CS%do_shearstress) ShSt(I,J,k) = sh_xy(I,J,kk)
    enddo
  else
    do concurrent (kk=1:kmax, J=js-2:Jeq+1, I=is-2:Ieq+1) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      sh_xy(I,J,kk) = mask2dBu(I,J) * ( dvdx(I,J,kk) + dudy(I,J,kk) )
      if (CS%do_shearstress) ShSt(I,J,k) = sh_xy(I,J,kk)
    enddo
  endif

  if (CS%use_Leithy) then
    ! Shearing strain (including no-slip boundary conditions at the 2-D land-sea mask).
    ! dudy_smooth and dvdx_smooth do not (yet) include modifications at OBCs from above.
    if (CS%no_slip) then
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
        sh_xy_smooth(I,J,kk) = (2.0-mask2dBu(I,J)) * ( dvdx_smooth(I,J,kk) + dudy_smooth(I,J,kk) )
      enddo
    else
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
        sh_xy_smooth(I,J,kk) = mask2dBu(I,J) * ( dvdx_smooth(I,J,kk) + dudy_smooth(I,J,kk) )
      enddo
    endif
  endif ! use Leith+E

  !  Evaluate Del2u = x.Div(Grad u) and Del2v = y.Div( Grad u)
  if (CS%biharmonic) then
    do concurrent (kk=1:kmax, j=js-1:Jeq+1, I=Isq-1:Ieq+1)
      Del2u(I,j,kk) = Idx2dyCu(I,j) * ((dx2q(I,J)*sh_xy(I,J,kk)) - (dx2q(I,J-1)*sh_xy(I,J-1,kk))) + &
                   Idxdy2u(I,j) * ((dy2h(i+1,j)*sh_xx(i+1,j,kk)) - (dy2h(i,j)*sh_xx(i,j,kk)))
    enddo
    do concurrent (kk=1:kmax, J=Jsq-1:Jeq+1, i=is-1:Ieq+1)
      Del2v(i,J,kk) = Idxdy2v(i,J) * ((dy2q(I,J)*sh_xy(I,J,kk)) - (dy2q(I-1,J)*sh_xy(I-1,J,kk))) - &
                   Idx2dyCv(i,J) * ((dx2h(i,j+1)*sh_xx(i,j+1,kk)) - (dx2h(i,j)*sh_xx(i,j,kk)))
    enddo
    if (apply_OBC) then ; if (OBC%zero_biharmonic) then
      do n=1,OBC%number_of_segments
        I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%is_N_or_S .and. (J >= Jsq-1) .and. (J <= Jeq+1)) then
          do concurrent (kk=1:kmax, I=OBC%segment(n)%HI%isd:OBC%segment(n)%HI%ied)
            Del2v(i,J,kk) = 0.
          enddo
        elseif (OBC%segment(n)%is_E_or_W .and. (I >= Isq-1) .and. (I <= Ieq+1)) then
          do concurrent (kk=1:kmax, j=OBC%segment(n)%HI%jsd:OBC%segment(n)%HI%jed)
            Del2u(I,j,kk) = 0.
          enddo
        endif
      enddo
    endif ; endif
  endif

  ! Vorticity
  if (use_vort_xy) then
    if (CS%no_slip) then
      do concurrent (kk=1:kmax, J=js_vort:je_vort, I=is_vort:ie_vort)
        vort_xy(I,J,kk) = (2.0-mask2dBu(I,J)) * ( dvdx(I,J,kk) - dudy(I,J,kk) )
      enddo
    else
      if (CS%use_circulation) then
        do concurrent (kk=1:kmax, J=js_vort:je_vort, I=is_vort:ie_vort) DO_LOCALITY(local(k))
          k = kstart + kk - 1
          vort_xy(I,J,kk) = mask2dBu(I,J) * IareaBu(I,J) * (  &
            ((v(i+1,J,k)*dyCv(i+1,J)) - (v(i,J,k)*dyCv(i,J)))  &
          - ((u(I,j+1,k)*dxCu(I,j+1)) - (u(I,j,k)*dxCu(I,j)))  &
           )
        enddo
      else
        do concurrent (kk=1:kmax, J=js_vort:je_vort, I=is_vort:ie_vort)
          vort_xy(I,J,kk) = mask2dBu(I,J) * ( dvdx(I,J,kk) - dudy(I,J,kk) )
        enddo
      endif
    endif
  endif

  if (CS%use_Leithy) then
    if (CS%no_slip) then
      do concurrent (kk=1:kmax, J=js_Kh-1:je_Kh, I=is_Kh-1:ie_Kh)
        vort_xy_smooth(I,J,kk) = (2.0-mask2dBu(I,J)) * ( dvdx_smooth(I,J,kk) - dudy_smooth(I,J,kk) )
      enddo
    else
      do concurrent (kk=1:kmax, J=js_Kh-1:je_Kh, I=is_Kh-1:ie_Kh)
        vort_xy_smooth(I,J,kk) = mask2dBu(I,J) * ( dvdx_smooth(I,J,kk) - dudy_smooth(I,J,kk) )
      enddo
    endif
  endif

end subroutine hor_visc_strain

!> Calculates the Laplacian viscosity and stress at h points, and starts the biharmonic
!! viscosity at h points (background, Smagorinsky and Leith terms) for one tile of
!! horizontal_viscosity.
subroutine hor_visc_h_coeffs_a(MEKE_a, VarMix_a, CS, metrics, bxH, bxQ, nkblock, h_neglect, &
                               inv_PI3, inv_PI6, use_kh_struct, use_MEKE_Ku, rescale_Kh, is_Kh, &
                               ie_Kh, js_Kh, je_Kh, u_a, v_a, h_a, Kh_h_a, grid_Re_Kh_a, &
                               div_xx_h_a, sh_xx_h_a, h_u_a, h_v_a, sh_xx_a, str_xx_a, &
                               grad_vort_mag_h_a, grad_vort_mag_h_2d_a, grad_div_mag_h_a, dudx_a, &
                               dvdy_a, sh_xy_a, Del2vort_q_a, Ah_a, Kh_a, Shear_mag_a, &
                               vert_vort_mag_a, hrat_min_a, visc_bound_rem_a)
  type(hor_visc_MEKE_type),      intent(inout) :: MEKE_a !< MEKE fields
                                                       !! related to Mesoscale Eddy Kinetic Energy.
  type(hor_visc_VarMix_type),    intent(inout) :: VarMix_a !< Variable mixing control structure
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  type(hv_metrics_CS), intent(in)    :: metrics       !< Persistent copies of the horizontal
                                                     !! viscosity metric and coefficient arrays
  type(Box_t),                   intent(in)    :: bxH  !< The h-point iteration box of this tile,
                                                       !! [isc:iec, jsc:jec, kstart:kend]
  type(Box_t),                   intent(in)    :: bxQ  !< The B-grid-index iteration box of this
                                                       !! tile, [IscB:IecB, JscB:JecB, kstart:kend]
  integer,                       intent(in)    :: nkblock !< The effective k-block size [nondim]
  real,    intent(in)    :: h_neglect !< thickness so small it can be lost in roundoff and so neglected
                                      !! [H ~> m or kg m-2]
  real,    intent(in)    :: inv_PI3 !< Powers of the inverse of pi [nondim]
  real,    intent(in)    :: inv_PI6 !< Powers of the inverse of pi [nondim]
  logical, intent(in)    :: use_kh_struct !< If true, the MEKE viscosity contributions are shaped with VarMix%BS_struct
  logical, intent(in)    :: use_MEKE_Ku !< If true, add the Laplacian viscosity derived from MEKE
  logical, intent(in)    :: rescale_Kh !< If true, the Laplacian viscosity is scaled by the resolution function
  integer, intent(in)    :: is_Kh !< Loop ranges for thickness point viscosities
  integer, intent(in)    :: ie_Kh !< Loop ranges for thickness point viscosities
  integer, intent(in)    :: js_Kh !< Loop ranges for thickness point viscosities
  integer, intent(in)    :: je_Kh !< Loop ranges for thickness point viscosities
  type(RealArray_t), intent(in)    :: u_a !< The zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: v_a !< The meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: h_a !< Layer thicknesses [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: Kh_h_a !< Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: grid_Re_Kh_a !< Grid Reynolds number for Laplacian horizontal viscosity
                                                   !! at h points [nondim]
  type(RealArray_t), intent(inout) :: div_xx_h_a !< horizontal divergence [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xx_h_a !< horizontal tension (du/dx - dv/dy) including metric terms
                                                !! [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: h_u_a !< Thickness interpolated to u points [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: h_v_a !< Thickness interpolated to v points [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: sh_xx_a !< horizontal tension (du/dx - dv/dy) including metric terms
                                              !! [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: str_xx_a !< str_xx is the diagonal term in the stress tensor
                                               !! [H L2 T-2 ~> m3 s-2 or kg s-2], but at some points in the
                                               !! code it is not yet layer integrated, so is in
                                               !! [L2 T-2 ~> m2 s-2].
  type(RealArray_t), intent(inout) :: grad_vort_mag_h_a !< Magnitude of vorticity gradient at h-points
                                                        !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: grad_vort_mag_h_2d_a !< Magnitude of 2d vorticity gradient at h-points
                                                           !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: grad_div_mag_h_a !< Magnitude of divergence gradient at h-points
                                                       !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: dudx_a !< components in the horizontal tension [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: dvdy_a !< components in the horizontal tension [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xy_a !< horizontal shearing strain (du/dy + dv/dx) including metric
                                              !! terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: Del2vort_q_a !< Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
  type(RealArray_t), intent(inout) :: Ah_a !< biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
  type(RealArray_t), intent(inout) :: Kh_a !< Laplacian  viscosity (h or q) [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: Shear_mag_a !< magnitude of the shear (h or q) [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: vert_vort_mag_a !< magnitude of the vertical vorticity gradient (h or q)
                                                      !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: hrat_min_a !< h_min divided by the thickness at the stress point (h or
                                                 !! q) [nondim]
  type(RealArray_t), intent(inout) :: visc_bound_rem_a !< fraction of overall viscous bounds that remain to be
                                                       !! applied (h or q) [nondim]

  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: visc_bound_rem, hrat_min, vert_vort_mag, &
    Shear_mag, Kh, Ah, Del2vort_q, sh_xy, dvdy, dudx, grad_div_mag_h, grad_vort_mag_h_2d, &
    grad_vort_mag_h, str_xx, sh_xx, h_v, h_u, sh_xx_h, div_xx_h, grid_Re_Kh, Kh_h, h, v, u
  real, dimension(:,:), contiguous, pointer :: Kh_bg_xx, Ah_bg_xx, Kh_Max_xx, n1n2_h, n1n1_m_n2n2_h, &
    grid_sp_h2, Laplac2_const_xx, Biharm6_const_xx, Laplac3_const_xx, Biharm_const_xx, Biharm_const2_xx
  real, dimension(:,:), contiguous, pointer :: MEKE_Ku, Res_fn_h
  real, dimension(:,:,:), contiguous, pointer :: BS_struct
  real :: Del2vort_h ! Laplacian of vorticity at h-points [L-2 T-1 ~> m-2 s-1]
  real :: AhSm       ! Smagorinsky biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: AhLth      ! 2D Leith biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: sh_xx_sq   ! Square of tension (sh_xx) [T-2 ~> s-2]
  real :: sh_xy_sq   ! Square of shearing strain (sh_xy) [T-2 ~> s-2]
  real :: h_min      ! Minimum h at the 4 neighboring velocity points [H ~> m]
  real :: Kh_max_here ! The local maximum Laplacian viscosity for stability [L2 T-1 ~> m2 s-1]
  real :: local_strain ! Local variable for interpolating computed strain rates [T-1 ~> s-1].
  real :: KE        ! Local kinetic energy [L2 T-2 ~> m2 s-2]
  real :: grad_vort ! Vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
  real :: grad_vort_qg ! QG-based vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
  real :: grid_Kh   ! Laplacian viscosity bound by grid [L2 T-1 ~> m2 s-1]
  integer :: i, j, k, kk, is, ie, js, je, Isq, Ieq, Jsq, Jeq, kstart, kmax

  call visc_bound_rem_a%view(visc_bound_rem) ; call hrat_min_a%view(hrat_min)
  call vert_vort_mag_a%view(vert_vort_mag) ; call Shear_mag_a%view(Shear_mag)
  call Kh_a%view(Kh) ; call Ah_a%view(Ah) ; call Del2vort_q_a%view(Del2vort_q)
  call sh_xy_a%view(sh_xy) ; call dvdy_a%view(dvdy) ; call dudx_a%view(dudx)
  call grad_div_mag_h_a%view(grad_div_mag_h) ; call grad_vort_mag_h_2d_a%view(grad_vort_mag_h_2d)
  call grad_vort_mag_h_a%view(grad_vort_mag_h) ; call str_xx_a%view(str_xx)
  call sh_xx_a%view(sh_xx) ; call h_v_a%view(h_v) ; call h_u_a%view(h_u)
  call sh_xx_h_a%view(sh_xx_h) ; call div_xx_h_a%view(div_xx_h)
  call grid_Re_Kh_a%view(grid_Re_Kh) ; call Kh_h_a%view(Kh_h) ; call h_a%view(h)
  call v_a%view(v) ; call u_a%view(u)

  nullify(Kh_bg_xx, Ah_bg_xx, Kh_Max_xx, n1n2_h, n1n1_m_n2n2_h, grid_sp_h2, Laplac2_const_xx, &
          Biharm6_const_xx, Laplac3_const_xx, Biharm_const_xx, Biharm_const2_xx)
  if (metrics%Kh_bg_xx%associated()) call metrics%Kh_bg_xx%view(Kh_bg_xx)
  if (metrics%Ah_bg_xx%associated()) call metrics%Ah_bg_xx%view(Ah_bg_xx)
  if (metrics%Kh_Max_xx%associated()) call metrics%Kh_Max_xx%view(Kh_Max_xx)
  if (metrics%n1n2_h%associated()) call metrics%n1n2_h%view(n1n2_h)
  if (metrics%n1n1_m_n2n2_h%associated()) call metrics%n1n1_m_n2n2_h%view(n1n1_m_n2n2_h)
  if (metrics%grid_sp_h2%associated()) call metrics%grid_sp_h2%view(grid_sp_h2)
  if (metrics%Laplac2_const_xx%associated()) call metrics%Laplac2_const_xx%view(Laplac2_const_xx)
  if (metrics%Biharm6_const_xx%associated()) call metrics%Biharm6_const_xx%view(Biharm6_const_xx)
  if (metrics%Laplac3_const_xx%associated()) call metrics%Laplac3_const_xx%view(Laplac3_const_xx)
  if (metrics%Biharm_const_xx%associated()) call metrics%Biharm_const_xx%view(Biharm_const_xx)
  if (metrics%Biharm_const2_xx%associated()) call metrics%Biharm_const2_xx%view(Biharm_const2_xx)

  nullify(MEKE_Ku, Res_fn_h, BS_struct)
  if (MEKE_a%Ku_a%associated()) call MEKE_a%Ku_a%view(MEKE_Ku)
  if (VarMix_a%Res_fn_h_a%associated()) call VarMix_a%Res_fn_h_a%view(Res_fn_h)
  if (VarMix_a%BS_struct_a%associated()) call VarMix_a%BS_struct_a%view(BS_struct)

  is  = bxH%idxS(1) ; ie  = bxH%idxE(1) ; js  = bxH%idxS(2) ; je  = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  kstart = bxH%idxS(3) ; kmax = bxH%idxE(3) - bxH%idxS(3) + 1

  if ((CS%Smagorinsky_Kh) .or. (CS%Smagorinsky_Ah)) then
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(sh_xx_sq, sh_xy_sq))
      sh_xx_sq = sh_xx(i,j,kk)**2
      sh_xy_sq = 0.25 * ( ((sh_xy(I-1,J-1,kk)**2) + (sh_xy(I,J,kk)**2)) &
                        + ((sh_xy(I-1,J,kk)**2) + (sh_xy(I,J-1,kk)**2)) )
      Shear_mag(i,j,kk) = sqrt(sh_xx_sq + sh_xy_sq)
    enddo
  endif

  if (CS%bound_Ah .or. CS%bound_Kh) then
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(k, h_min))
      k = kstart + kk - 1
      h_min = min(h_u(I,j,kk), h_u(I-1,j,kk), h_v(i,J,kk), h_v(i,J-1,kk))
      hrat_min(i,j,kk) = min(1.0, h_min / (h(i,j,k) + h_neglect))
    enddo
  endif

  if (CS%Laplacian) then
    ! Determine the Laplacian viscosity at h points, using the
    ! largest value from several parameterizations. Also get
    ! the Laplacian component of str_xx.

    if ((CS%Leith_Kh) .or. (CS%Leith_Ah) .or. (CS%use_Leithy)) then
      if (CS%use_QG_Leith_visc) then
        do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(grad_vort, grad_vort_qg))
          grad_vort = grad_vort_mag_h(i,j,kk) + grad_div_mag_h(i,j,kk)
          grad_vort_qg = 3. * grad_vort_mag_h_2d(i,j,kk)
          vert_vort_mag(i,j,kk) = min(grad_vort, grad_vort_qg)
        enddo
      else
        do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
          vert_vort_mag(i,j,kk) = grad_vort_mag_h(i,j,kk) + grad_div_mag_h(i,j,kk)
        enddo
      endif
    endif

    ! Static (pre-computed) background viscosity
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
      Kh(i,j,kk) = Kh_bg_xx(i,j)
    enddo

    ! NOTE: The following do-block can be decomposed and vectorized after the
    !   stack size has been reduced.
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
      if (CS%add_LES_viscosity) then
        if (CS%Smagorinsky_Kh) &
          Kh(i,j,kk) = Kh(i,j,kk) + Laplac2_const_xx(i,j) * Shear_mag(i,j,kk)
        if (CS%Leith_Kh) &
          Kh(i,j,kk) = Kh(i,j,kk) + Laplac3_const_xx(i,j) * vert_vort_mag(i,j,kk) * inv_PI3
      else
        if (CS%Smagorinsky_Kh) &
          Kh(i,j,kk) = max(Kh(i,j,kk), Laplac2_const_xx(i,j) * Shear_mag(i,j,kk))
        if (CS%Leith_Kh) &
          Kh(i,j,kk) = max(Kh(i,j,kk), Laplac3_const_xx(i,j) * vert_vort_mag(i,j,kk) * inv_PI3)
      endif
    enddo

    ! All viscosity contributions above are subject to resolution scaling

    if (rescale_Kh) then
      do kk=1,kmax ; do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
        Kh(i,j,kk) = Res_fn_h(i,j) * Kh(i,j,kk)
      enddo ; enddo ; enddo
    endif

    ! Place a floor on the viscosity, if desired.
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
      Kh(i,j,kk) = max(Kh(i,j,kk), CS%Kh_bg_min)
    enddo

    if (use_MEKE_Ku .and. .not. CS%EY24_EBT_BS) then
      ! *Add* the MEKE contribution (which might be negative)
      if (use_kh_struct) then
        if (CS%res_scale_MEKE) then
          do kk=1,kmax ; do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            k = kstart + kk - 1
            Kh(i,j,kk) = Kh(i,j,kk) + MEKE_Ku(i,j) * Res_fn_h(i,j) * BS_struct(i,j,k)
          enddo ; enddo ; enddo
        else
          do kk=1,kmax ; do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            k = kstart + kk - 1
            Kh(i,j,kk) = Kh(i,j,kk) + MEKE_Ku(i,j) * BS_struct(i,j,k)
          enddo ; enddo ; enddo
        endif
      else
        if (CS%res_scale_MEKE) then
          do kk=1,kmax ; do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            Kh(i,j,kk) = Kh(i,j,kk) + MEKE_Ku(i,j) * Res_fn_h(i,j)
          enddo ; enddo ; enddo
        else
          do kk=1,kmax ; do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            Kh(i,j,kk) = Kh(i,j,kk) + MEKE_Ku(i,j)
          enddo ; enddo ; enddo
        endif
      endif
    endif

    if (CS%anisotropic) then
      do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
        ! *Add* the tension component of anisotropic viscosity
        Kh(i,j,kk) = Kh(i,j,kk) + CS%Kh_aniso * (1. - n1n2_h(i,j)**2)
      enddo
    endif


    ! Newer method of bounding for stability
    if ((CS%bound_Kh) .and. (CS%bound_Ah)) then
      do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(Kh_max_here))
        visc_bound_rem(i,j,kk) = 1.0
        Kh_max_here = hrat_min(i,j,kk) * Kh_Max_xx(i,j)
        if (Kh(i,j,kk) >= Kh_max_here) then
          visc_bound_rem(i,j,kk) = 0.0
          Kh(i,j,kk) = Kh_max_here
        elseif ((Kh(i,j,kk) > 0.0) .or. (CS%backscatter_underbound .and. (Kh_max_here > 0.0))) then
          visc_bound_rem(i,j,kk) = 1.0 - Kh(i,j,kk) / Kh_max_here
        endif
      enddo
    elseif (CS%bound_Kh) then
      do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
        Kh(i,j,kk) = min(Kh(i,j,kk), hrat_min(i,j,kk) * Kh_Max_xx(i,j))
      enddo
    endif

    ! In Leith+E parameterization Kh is computed after Ah in the biharmonic loop.
    ! The harmonic component of str_xx is added in the biharmonic loop.
    if (CS%use_Leithy) then
      do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
        Kh(i,j,kk) = 0.
      enddo
    endif

    if (CS%do_Kh_h .or. CS%debug) then
      do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        Kh_h(i,j,k) = Kh(i,j,kk)
      enddo
    endif

    if (CS%do_grid_Re_Kh) then
      do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k, KE, grid_Kh))
        k = kstart + kk - 1
        KE = 0.125*(((u(I,j,k)+u(I-1,j,k))**2) + ((v(i,J,k)+v(i,J-1,k))**2))
        grid_Kh = max(Kh(i,j,kk), CS%min_grid_Kh)
        grid_Re_Kh(i,j,k) = (sqrt(KE) * sqrt(grid_sp_h2(i,j))) / grid_Kh
      enddo
    endif

    if (CS%do_div_xx_h) then
      do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        div_xx_h(i,j,k) = dudx(i,j,kk) + dvdy(i,j,kk)
      enddo
    endif

    if (CS%do_sh_xx_h) then
      do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        sh_xx_h(i,j,k) = sh_xx(i,j,kk)
      enddo
    endif

    do concurrent (kk=1:kmax, j=Jsq:Jeq+1, i=Isq:Ieq+1)
      str_xx(i,j,kk) = -Kh(i,j,kk) * sh_xx(i,j,kk)
    enddo
  else
    do concurrent (kk=1:kmax, j=Jsq:Jeq+1, i=Isq:Ieq+1)
      str_xx(i,j,kk) = 0.0
    enddo
  endif ! Get Kh at h points and get Laplacian component of str_xx

  if (CS%anisotropic) then
    do concurrent (kk=1:kmax, j=Jsq:Jeq+1, i=Isq:Ieq+1) DO_LOCALITY(local(local_strain))
      ! Shearing-strain averaged to h-points
      local_strain = 0.25 * ( (sh_xy(I,J,kk) + sh_xy(I-1,J-1,kk)) + (sh_xy(I-1,J,kk) + sh_xy(I,J-1,kk)) )
      ! *Add* the shear-strain contribution to the xx-component of stress
      str_xx(i,j,kk) = str_xx(i,j,kk) - CS%Kh_aniso * n1n2_h(i,j) * n1n1_m_n2n2_h(i,j) * local_strain
    enddo
  endif

  if (CS%biharmonic) then
    ! Determine the biharmonic viscosity at h points, using the
    ! largest value from several parameterizations. Also get the
    ! biharmonic component of str_xx.
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
      Ah(i,j,kk) = Ah_bg_xx(i,j)
    enddo

    if ((CS%Smagorinsky_Ah) .or. (CS%Leith_Ah) .or. (CS%use_Leithy)) then
      if (CS%Smagorinsky_Ah) then
        if (CS%bound_Coriolis) then
          do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(AhSm))
            AhSm = Shear_mag(i,j,kk) * (Biharm_const_xx(i,j) &
                + Biharm_const2_xx(i,j) * Shear_mag(i,j,kk))
            Ah(i,j,kk) = max(Ah(i,j,kk), AhSm)
          enddo
        else
          do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(AhSm))
            AhSm = Biharm_const_xx(i,j) * Shear_mag(i,j,kk)
            Ah(i,j,kk) = max(Ah(i,j,kk), AhSm)
          enddo
        endif
      endif

      if (CS%Leith_Ah) then
        do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(Del2vort_h, AhLth))
          Del2vort_h = 0.25 * ((Del2vort_q(I,J,kk) + Del2vort_q(I-1,J-1,kk)) + &
                               (Del2vort_q(I-1,J,kk) + Del2vort_q(I,J-1,kk)))
          AhLth = Biharm6_const_xx(i,j) * abs(Del2vort_h) * inv_PI6
          Ah(i,j,kk) = max(Ah(i,j,kk), AhLth)
        enddo
      endif
    endif ! Smagorinsky_Ah or Leith_Ah or Leith+E
  endif

end subroutine hor_visc_h_coeffs_a

!> Completes the biharmonic viscosity and stress at h points, adds the h-point backscatter,
!! and calculates the gradients of the Laplacian and the thicknesses at q points for one tile
!! of horizontal_viscosity.
subroutine hor_visc_h_coeffs_b(MEKE_a, VarMix_a, Gcore, CS, metrics, OBC, bxH, bxQ, nkblock, &
                               h_neglect, h_neglect3, apply_OBC, use_kh_struct, use_MEKE_Au, is_Kh, &
                               ie_Kh, js_Kh, je_Kh, u_a, v_a, h_a, Ah_h_a, Kh_h_a, grid_Re_Ah_a, &
                               visc_limit_h_a, visc_limit_h_flag_a, visc_limit_h_frac_a, &
                               BS_coeff_h_a, Del2u_a, h_u_a, Del2v_a, h_v_a, sh_xx_a, &
                               sh_xx_smooth_a, str_xx_a, bhstr_xx_a, m_leithy_a, dDel2vdx_a, &
                               dDel2udy_a, sh_xy_a, hq_a, Ah_a, Kh_a, Shear_mag_a, hrat_min_a, &
                               visc_bound_rem_a)
  type(hor_visc_MEKE_type),      intent(inout) :: MEKE_a !< MEKE fields
                                                       !! related to Mesoscale Eddy Kinetic Energy.
  type(hor_visc_VarMix_type),    intent(inout) :: VarMix_a !< Variable mixing control structure
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  type(hv_metrics_CS), intent(in)    :: metrics       !< Persistent copies of the horizontal
                                                     !! viscosity metric and coefficient arrays
  type(ocean_OBC_type), optional, pointer      :: OBC  !< Pointer to an open boundary condition type
  type(Box_t),                   intent(in)    :: bxH  !< The h-point iteration box of this tile,
                                                       !! [isc:iec, jsc:jec, kstart:kend]
  type(Box_t),                   intent(in)    :: bxQ  !< The B-grid-index iteration box of this
                                                       !! tile, [IscB:IecB, JscB:JecB, kstart:kend]
  integer,                       intent(in)    :: nkblock !< The effective k-block size [nondim]
  real,    intent(in)    :: h_neglect !< thickness so small it can be lost in roundoff and so neglected
                                      !! [H ~> m or kg m-2]
  real,    intent(in)    :: h_neglect3 !< h_neglect^3 [H3 ~> m3 or kg3 m-6]
  logical, intent(in)    :: apply_OBC !< If true, open boundary conditions are applied on this PE
  logical, intent(in)    :: use_kh_struct !< If true, the MEKE viscosity contributions are shaped with VarMix%BS_struct
  logical, intent(in)    :: use_MEKE_Au !< If true, add the biharmonic viscosity derived from MEKE
  integer, intent(in)    :: is_Kh !< Loop ranges for thickness point viscosities
  integer, intent(in)    :: ie_Kh !< Loop ranges for thickness point viscosities
  integer, intent(in)    :: js_Kh !< Loop ranges for thickness point viscosities
  integer, intent(in)    :: je_Kh !< Loop ranges for thickness point viscosities
  type(RealArray_t), intent(in)    :: u_a !< The zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: v_a !< The meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: h_a !< Layer thicknesses [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: Ah_h_a !< biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
  type(RealArray_t), intent(inout) :: Kh_h_a !< Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: grid_Re_Ah_a !< Grid Reynolds number for Biharmonic horizontal viscosity
                                                   !! at h points [nondim]
  type(RealArray_t), intent(inout) :: visc_limit_h_a !< Used to stabilize the EY24_EBT_BS backscatter [nondim]
  type(RealArray_t), intent(inout) :: visc_limit_h_flag_a !< determines whether backscatter is shut off
                                                          !! [nondim]
  type(RealArray_t), intent(inout) :: visc_limit_h_frac_a !< determines how close backscatter is to shutting
                                                          !! off [nondim]
  type(RealArray_t), intent(inout) :: BS_coeff_h_a !< A diagnostic array of the backscatter coefficient
                                                   !! [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: Del2u_a !< The u-component of the Laplacian of velocity
                                              !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: h_u_a !< Thickness interpolated to u points [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: Del2v_a !< The v-component of the Laplacian of velocity
                                              !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: h_v_a !< Thickness interpolated to v points [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: sh_xx_a !< horizontal tension (du/dx - dv/dy) including metric terms
                                              !! [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xx_smooth_a !< horizontal tension from smoothed velocity including
                                                     !! metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: str_xx_a !< str_xx is the diagonal term in the stress tensor
                                               !! [H L2 T-2 ~> m3 s-2 or kg s-2], but at some points in the
                                               !! code it is not yet layer integrated, so is in
                                               !! [L2 T-2 ~> m2 s-2].
  type(RealArray_t), intent(inout) :: bhstr_xx_a !< A copy of str_xx that only contains the biharmonic
                                                 !! contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
  type(RealArray_t), intent(inout) :: m_leithy_a !< Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
  type(RealArray_t), intent(inout) :: dDel2vdx_a !< Components in the biharmonic equivalent of the shearing
                                                 !! strain [L-2 T-1 ~> m-2 s-1]
  type(RealArray_t), intent(inout) :: dDel2udy_a !< Components in the biharmonic equivalent of the shearing
                                                 !! strain [L-2 T-1 ~> m-2 s-1]
  type(RealArray_t), intent(inout) :: sh_xy_a !< horizontal shearing strain (du/dy + dv/dx) including metric
                                              !! terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: hq_a !< harmonic mean of the harmonic means of the u- & v point
                                           !! thicknesses [H ~> m or kg m-2] This form guarantees that hq/hu <
                                           !! 4.
  type(RealArray_t), intent(inout) :: Ah_a !< biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
  type(RealArray_t), intent(inout) :: Kh_a !< Laplacian  viscosity (h or q) [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: Shear_mag_a !< magnitude of the shear (h or q) [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: hrat_min_a !< h_min divided by the thickness at the stress point (h or
                                                 !! q) [nondim]
  type(RealArray_t), intent(inout) :: visc_bound_rem_a !< fraction of overall viscous bounds that remain to be
                                                       !! applied (h or q) [nondim]

  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: visc_bound_rem, hrat_min, Shear_mag, Kh, Ah, hq, &
    sh_xy, dDel2udy, dDel2vdx, m_leithy, bhstr_xx, str_xx, sh_xx_smooth, sh_xx, h_v, Del2v, h_u, &
    Del2u, BS_coeff_h, visc_limit_h_frac, visc_limit_h_flag, visc_limit_h, grid_Re_Ah, Kh_h, Ah_h, &
    h, v, u
  real, dimension(:,:), contiguous, pointer :: IdyCu, IdxCv, IdyCv, IdxCu, mask2dBu, mask2dCu, &
    mask2dCv
  real, dimension(:,:), contiguous, pointer :: reduction_xx, Ah_Max_xx, Ah_Max_xx_KS, grid_sp_h3, dx_dyT, &
    dy_dxT, dx_dyBu, dy_dxBu, Re_Ah_const_xx
  real, dimension(:,:), contiguous, pointer :: MEKE_Au
  real :: sh_xx_sq   ! Square of tension (sh_xx) [T-2 ~> s-2]
  real :: sh_xy_sq   ! Square of shearing strain (sh_xy) [T-2 ~> s-2]
  real :: h2uq, h2vq ! temporary variables [H2 ~> m2 or kg2 m-4].
  real :: hu, hv     ! Thicknesses interpolated by arithmetic means to corner
                     ! points; these are first interpolated to u or v velocity
                     ! points where masks are applied [H ~> m or kg m-2].
  real :: h_min      ! Minimum h at the 4 neighboring velocity points [H ~> m]
  real :: KE        ! Local kinetic energy [L2 T-2 ~> m2 s-2]
  real :: d_del2u   ! dy-weighted Laplacian(u) diff in x [L-2 T-1 ~> m-2 s-1]
  real :: d_del2v   ! dx-weighted Laplacian(v) diff in y [L-2 T-1 ~> m-2 s-1]
  real :: d_str     ! Stress tensor update [L2 T-2 ~> m2 s-2]
  real :: grid_Ah   ! Biharmonic viscosity bound by grid [L4 T-1 ~> m4 s-1]
  real :: tmp
  integer :: i, j, k, kk, n, is, ie, js, je, Isq, Ieq, Jsq, Jeq, kstart, kmax

  call visc_bound_rem_a%view(visc_bound_rem) ; call hrat_min_a%view(hrat_min)
  call Shear_mag_a%view(Shear_mag) ; call Kh_a%view(Kh) ; call Ah_a%view(Ah) ; call hq_a%view(hq)
  call sh_xy_a%view(sh_xy) ; call dDel2udy_a%view(dDel2udy) ; call dDel2vdx_a%view(dDel2vdx)
  call m_leithy_a%view(m_leithy) ; call bhstr_xx_a%view(bhstr_xx) ; call str_xx_a%view(str_xx)
  call sh_xx_smooth_a%view(sh_xx_smooth) ; call sh_xx_a%view(sh_xx) ; call h_v_a%view(h_v)
  call Del2v_a%view(Del2v) ; call h_u_a%view(h_u) ; call Del2u_a%view(Del2u)
  call BS_coeff_h_a%view(BS_coeff_h) ; call visc_limit_h_frac_a%view(visc_limit_h_frac)
  call visc_limit_h_flag_a%view(visc_limit_h_flag) ; call visc_limit_h_a%view(visc_limit_h)
  call grid_Re_Ah_a%view(grid_Re_Ah) ; call Kh_h_a%view(Kh_h) ; call Ah_h_a%view(Ah_h)
  call h_a%view(h) ; call v_a%view(v) ; call u_a%view(u)
  call Gcore%IdyCu%view(IdyCu) ; call Gcore%IdxCv%view(IdxCv) ; call Gcore%IdyCv%view(IdyCv)
  call Gcore%IdxCu%view(IdxCu) ; call Gcore%mask2dBu%view(mask2dBu)
  call Gcore%mask2dCu%view(mask2dCu) ; call Gcore%mask2dCv%view(mask2dCv)

  nullify(reduction_xx, Ah_Max_xx, Ah_Max_xx_KS, grid_sp_h3, dx_dyT, dy_dxT, dx_dyBu, dy_dxBu, Re_Ah_const_xx)
  if (metrics%reduction_xx%associated()) call metrics%reduction_xx%view(reduction_xx)
  if (metrics%Ah_Max_xx%associated()) call metrics%Ah_Max_xx%view(Ah_Max_xx)
  if (metrics%Ah_Max_xx_KS%associated()) call metrics%Ah_Max_xx_KS%view(Ah_Max_xx_KS)
  if (metrics%grid_sp_h3%associated()) call metrics%grid_sp_h3%view(grid_sp_h3)
  if (metrics%dx_dyT%associated()) call metrics%dx_dyT%view(dx_dyT)
  if (metrics%dy_dxT%associated()) call metrics%dy_dxT%view(dy_dxT)
  if (metrics%dx_dyBu%associated()) call metrics%dx_dyBu%view(dx_dyBu)
  if (metrics%dy_dxBu%associated()) call metrics%dy_dxBu%view(dy_dxBu)
  if (metrics%Re_Ah_const_xx%associated()) call metrics%Re_Ah_const_xx%view(Re_Ah_const_xx)

  nullify(MEKE_Au)
  if (MEKE_a%Au_a%associated()) call MEKE_a%Au_a%view(MEKE_Au)

  is  = bxH%idxS(1) ; ie  = bxH%idxE(1) ; js  = bxH%idxS(2) ; je  = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  kstart = bxH%idxS(3) ; kmax = bxH%idxE(3) - bxH%idxS(3) + 1

  if (CS%biharmonic) then

    if (use_MEKE_Au) then
      ! *Add* the MEKE contribution
      do kk=1,kmax ; do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
        Ah(i,j,kk) = Ah(i,j,kk) + MEKE_Au(i,j)
      enddo ; enddo ; enddo
    endif

    if (CS%Re_Ah > 0.0) then
      do kk=1,kmax ; do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
        k = kstart + kk - 1
        KE = 0.125*(((u(I,j,k)+u(I-1,j,k))**2) + ((v(i,J,k)+v(i,J-1,k))**2))
        Ah(i,j,kk) = sqrt(KE) * Re_Ah_const_xx(i,j)
      enddo ; enddo ; enddo
    endif

    if (CS%bound_Ah) then
      if (CS%bound_Kh) then
        do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
          Ah(i,j,kk) = min(Ah(i,j,kk), visc_bound_rem(i,j,kk) * hrat_min(i,j,kk) * Ah_Max_xx(i,j))
        enddo
      else
        do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
          Ah(i,j,kk) = min(Ah(i,j,kk), hrat_min(i,j,kk) * Ah_Max_xx(i,j))
        enddo
      endif
    endif

    if (CS%EY24_EBT_BS) then
        do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(k, tmp))
          k = kstart + kk - 1
          tmp = CS%KS_coef * hrat_min(i,j,kk) * Ah_Max_xx_KS(i,j)
          visc_limit_h(i,j,k) = tmp
          visc_limit_h_frac(i,j,k) = Ah(i,j,kk) / (CS%KS_coef * hrat_min(i,j,kk) * Ah_Max_xx_KS(i,j))
          if (Ah(i,j,kk) >= tmp) then
            visc_limit_h_flag(i,j,k) = 1.
          endif
        enddo
    endif

    if ((CS%do_Ah_h) .or. CS%debug .or. CS%use_Leithy) then
      do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        Ah_h(i,j,k) = Ah(i,j,kk)
      enddo
    endif

    if (CS%use_Leithy) then
      ! Compute Leith+E Kh after bounds have been applied to Ah
      ! and after it has been smoothed. Kh = -m_leithy * Ah
      do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        Kh(i,j,kk) = -m_leithy(i,j,kk) * Ah(i,j,kk)
        Kh_h(i,j,k) = Kh(i,j,kk)
      enddo
    endif

    if (CS%do_grid_Re_Ah) then
      do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k, KE, grid_Ah))
        k = kstart + kk - 1
        KE = 0.125 * (((u(I,j,k) + u(I-1,j,k))**2) + ((v(i,J,k) + v(i,J-1,k))**2))
        grid_Ah = max(Ah(i,j,kk), CS%min_grid_Ah)
        grid_Re_Ah(i,j,k) = (sqrt(KE) * grid_sp_h3(i,j)) / grid_Ah
      enddo
    endif

    do concurrent (kk=1:kmax, j=Jsq:Jeq+1, i=Isq:Ieq+1) DO_LOCALITY(local(k, d_del2u, d_del2v, d_str))
      k = kstart + kk - 1
      d_del2u = (IdyCu(I,j) * Del2u(I,j,kk)) - (IdyCu(I-1,j) * Del2u(I-1,j,kk))
      d_del2v = (IdxCv(i,J) * Del2v(i,J,kk)) - (IdxCv(i,J-1) * Del2v(i,J-1,kk))
      d_str = Ah(i,j,kk) * ((DY_dxT(i,j) * d_del2u) - (DX_dyT(i,j) * d_del2v))

      str_xx(i,j,kk) = str_xx(i,j,kk) + d_str

      if (CS%use_Leithy) str_xx(i,j,kk) = str_xx(i,j,kk) - Kh(i,j,kk) * sh_xx_smooth(i,j,kk)

      ! Keep a copy of the biharmonic contribution for backscatter parameterization
      bhstr_xx(i,j,kk) = d_str * (h(i,j,k) * reduction_xx(i,j))
    enddo
  endif ! Get biharmonic coefficient at h points and biharmonic part of str_xx

  ! Backscatter using MEKE
  if (CS%EY24_EBT_BS) then
    call hor_visc_backscatter_h(CS, MEKE_a, VarMix_a, use_kh_struct, nkblock, kstart, kmax, Isq, &
                                Ieq, Jsq, Jeq, visc_limit_h_flag_a, sh_xx_a, str_xx_a, BS_coeff_h_a)
  endif ! Backscatter

  if (CS%biharmonic) then
    ! Gradient of Laplacian, for use in bi-harmonic term
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
      dDel2vdx(I,J,kk) = DY_dxBu(I,J)*((Del2v(i+1,J,kk)*IdyCv(i+1,J)) - (Del2v(i,J,kk)*IdyCv(i,J)))
      dDel2udy(I,J,kk) = DX_dyBu(I,J)*((Del2u(I,j+1,kk)*IdxCu(I,j+1)) - (Del2u(I,j,kk)*IdxCu(I,j)))
    enddo
    ! Adjust contributions to shearing strain on open boundaries.
    if (apply_OBC) then ; if ((OBC%strain_config == OBC_STRAIN_ZERO) .or. &
                              (OBC%strain_config == OBC_STRAIN_FREESLIP)) then
      do n=1,OBC%number_of_segments
        J = OBC%segment(n)%HI%JsdB ; I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%is_N_or_S .and. (J >= js-1) .and. (J <= Jeq)) then
          do concurrent (kk=1:kmax, I=OBC%segment(n)%HI%IsdB:OBC%segment(n)%HI%IedB)
            if (OBC%strain_config == OBC_STRAIN_ZERO) then
              dDel2vdx(I,J,kk) = 0. ; dDel2udy(I,J,kk) = 0.
            elseif (OBC%strain_config == OBC_STRAIN_FREESLIP) then
              dDel2udy(I,J,kk) = 0.
            endif
          enddo
        elseif (OBC%segment(n)%is_E_or_W .and. (I >= is-1) .and. (I <= Ieq)) then
          do concurrent (kk=1:kmax, J=OBC%segment(n)%HI%JsdB:OBC%segment(n)%HI%JedB)
            if (OBC%strain_config == OBC_STRAIN_ZERO) then
              dDel2vdx(I,J,kk) = 0. ; dDel2udy(I,J,kk) = 0.
            elseif (OBC%strain_config == OBC_STRAIN_FREESLIP) then
              dDel2vdx(I,J,kk) = 0.
            endif
          enddo
        endif
      enddo
    endif ; endif
  endif

  if ((CS%Smagorinsky_Kh) .or. (CS%Smagorinsky_Ah)) then
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(sh_xy_sq, sh_xx_sq))
      sh_xy_sq = sh_xy(I,J,kk)**2
      sh_xx_sq = 0.25 * ( ((sh_xx(i,j,kk)**2) + (sh_xx(i+1,j+1,kk)**2)) &
                        + ((sh_xx(i,j+1,kk)**2) + (sh_xx(i+1,j,kk)**2)) )
      Shear_mag(I,J,kk) = sqrt(sh_xy_sq + sh_xx_sq)
    enddo
  endif

  do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(h2uq, h2vq))
    h2uq = 4.0 * (h_u(I,j,kk) * h_u(I,j+1,kk))
    h2vq = 4.0 * (h_v(i,J,kk) * h_v(i+1,J,kk))
    hq(I,J,kk) = (2.0 * (h2uq * h2vq)) &
        / (h_neglect3 + (h2uq + h2vq) * ((h_u(I,j,kk) + h_u(I,j+1,kk)) + (h_v(i,J,kk) + h_v(i+1,J,kk))))
  enddo

  if (CS%bound_Ah .or. CS%bound_Kh) then
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(h_min))
      h_min = min(h_u(I,j,kk), h_u(I,j+1,kk), h_v(i,J,kk), h_v(i+1,J,kk))
      hrat_min(I,J,kk) = min(1.0, h_min / (hq(I,J,kk) + h_neglect))
    enddo

  endif

  ! TODO: GPU??  Are h_[uv] on CPU?  update to hrat_min?
  if (CS%no_slip) then
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(hu, hv))
      if (CS%no_slip .and. (mask2dBu(I,J) < 0.5)) then
        if ((mask2dCu(I,j) + mask2dCu(I,j+1)) + &
            (mask2dCv(i,J) + mask2dCv(i+1,J)) > 0.0) then
          ! This is a coastal vorticity point, so modify hq and hrat_min.

          hu = mask2dCu(I,j) * h_u(I,j,kk) + mask2dCu(I,j+1) * h_u(I,j+1,kk)
          hv = mask2dCv(i,J) * h_v(i,J,kk) + mask2dCv(i+1,J) * h_v(i+1,J,kk)
          if ((mask2dCu(I,j) + mask2dCu(I,j+1)) * &
              (mask2dCv(i,J) + mask2dCv(i+1,J)) == 0.0) then
            ! Only one of hu and hv is nonzero, so just add them.
            hq(I,J,kk) = hu + hv
            hrat_min(I,J,kk) = 1.0
          else
            ! Both hu and hv are nonzero, so take the harmonic mean.
            hq(I,J,kk) = 2.0 * (hu * hv) / ((hu + hv) + h_neglect)
            hrat_min(I,J,kk) = min(1.0, min(hu, hv) / (hq(I,J,kk) + h_neglect) )
          endif
        endif
      endif
    enddo
  endif

end subroutine hor_visc_h_coeffs_b

!> Calculates the Laplacian and biharmonic viscosities and stresses at q points, adds the
!! q-point backscatter, and starts the GME stress terms for one tile of horizontal_viscosity.
subroutine hor_visc_q_coeffs(MEKE_a, VarMix_a, Gcore, CS, metrics, bxH, bxQ, nkblock, inv_PI3, &
                             inv_PI6, find_FrictWork, use_kh_struct, use_MEKE_Ku, use_MEKE_Au, &
                             rescale_Kh, u_a, v_a, sh_xx_bt_a, sh_xy_bt_a, GME_effic_h_a, &
                             GME_effic_q_a, KH_u_GME_a, KH_v_GME_a, GME_coeff_h_a, GME_coeff_q_a, &
                             str_xx_GME_a, str_xy_GME_a, Ah_h_a, Ah_q_a, Kh_h_a, Kh_q_a, &
                             vort_xy_q_a, sh_xy_q_a, visc_limit_q_a, visc_limit_q_flag_a, &
                             visc_limit_q_frac_a, BS_coeff_q_a, sh_xx_a, str_xx_a, dDel2vdx_a, &
                             dDel2udy_a, sh_xy_a, sh_xy_smooth_a, str_xy_a, bhstr_xy_a, vort_xy_a, &
                             grad_vort_mag_q_a, grad_vort_mag_q_2d_a, Del2vort_q_a, &
                             grad_div_mag_q_a, hq_a, Ah_a, Kh_a, Shear_mag_a, vert_vort_mag_a, &
                             hrat_min_a, visc_bound_rem_a)
  type(hor_visc_MEKE_type),      intent(inout) :: MEKE_a !< MEKE fields
                                                       !! related to Mesoscale Eddy Kinetic Energy.
  type(hor_visc_VarMix_type),    intent(inout) :: VarMix_a !< Variable mixing control structure
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  type(hv_metrics_CS), intent(in)    :: metrics       !< Persistent copies of the horizontal
                                                     !! viscosity metric and coefficient arrays
  type(Box_t),                   intent(in)    :: bxH  !< The h-point iteration box of this tile,
                                                       !! [isc:iec, jsc:jec, kstart:kend]
  type(Box_t),                   intent(in)    :: bxQ  !< The B-grid-index iteration box of this
                                                       !! tile, [IscB:IecB, JscB:JecB, kstart:kend]
  integer,                       intent(in)    :: nkblock !< The effective k-block size [nondim]
  real,    intent(in)    :: inv_PI3 !< Powers of the inverse of pi [nondim]
  real,    intent(in)    :: inv_PI6 !< Powers of the inverse of pi [nondim]
  logical, intent(in)    :: find_FrictWork !< If true, the work done by lateral friction is calculated
  logical, intent(in)    :: use_kh_struct !< If true, the MEKE viscosity contributions are shaped with VarMix%BS_struct
  logical, intent(in)    :: use_MEKE_Ku !< If true, add the Laplacian viscosity derived from MEKE
  logical, intent(in)    :: use_MEKE_Au !< If true, add the biharmonic viscosity derived from MEKE
  logical, intent(in)    :: rescale_Kh !< If true, the Laplacian viscosity is scaled by the resolution function
  type(RealArray_t), intent(in)    :: u_a !< The zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: v_a !< The meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: sh_xx_bt_a !< barotropic horizontal tension (du/dx - dv/dy) including
                                                 !! metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: sh_xy_bt_a !< barotropic horizontal shearing strain (du/dy + dv/dx) inc.
                                                 !! metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: GME_effic_h_a !< The filtered efficiency of the GME terms at h points
                                                    !! [nondim]
  type(RealArray_t), intent(in)    :: GME_effic_q_a !< The filtered efficiency of the GME terms at q points
                                                    !! [nondim]
  type(RealArray_t), intent(in)    :: KH_u_GME_a !< Isopycnal height diffusivities in u-columns
                                                 !! [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(in)    :: KH_v_GME_a !< Isopycnal height diffusivities in v-columns
                                                 !! [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: GME_coeff_h_a !< GME coefficient at h-points [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: GME_coeff_q_a !< GME coeff. at q-points [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: str_xx_GME_a !< smoothed diagonal term in the stress tensor from GME
                                                   !! [L2 T-2 ~> m2 s-2]
  type(RealArray_t), intent(inout) :: str_xy_GME_a !< smoothed cross term in the stress tensor from GME
                                                   !! [L2 T-2 ~> m2 s-2]
  type(RealArray_t), intent(inout) :: Ah_h_a !< biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
  type(RealArray_t), intent(inout) :: Ah_q_a !< biharmonic viscosity at corner points [L4 T-1 ~> m4 s-1]
  type(RealArray_t), intent(inout) :: Kh_h_a !< Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: Kh_q_a !< Laplacian viscosity at corner points [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: vort_xy_q_a !< vertical vorticity at corner points [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xy_q_a !< horizontal shearing strain at corner points [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: visc_limit_q_a !< used to stabilize the EY24_EBT_BS backscatter [nondim]
  type(RealArray_t), intent(inout) :: visc_limit_q_flag_a !< determines whether backscatter is shut off
                                                          !! [nondim]
  type(RealArray_t), intent(inout) :: visc_limit_q_frac_a !< determines how close backscatter is to shutting
                                                          !! off [nondim]
  type(RealArray_t), intent(inout) :: BS_coeff_q_a !< A diagnostic array of the backscatter coefficient
                                                   !! [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: sh_xx_a !< horizontal tension (du/dx - dv/dy) including metric terms
                                              !! [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: str_xx_a !< str_xx is the diagonal term in the stress tensor
                                               !! [H L2 T-2 ~> m3 s-2 or kg s-2], but at some points in the
                                               !! code it is not yet layer integrated, so is in
                                               !! [L2 T-2 ~> m2 s-2].
  type(RealArray_t), intent(inout) :: dDel2vdx_a !< Components in the biharmonic equivalent of the shearing
                                                 !! strain [L-2 T-1 ~> m-2 s-1]
  type(RealArray_t), intent(inout) :: dDel2udy_a !< Components in the biharmonic equivalent of the shearing
                                                 !! strain [L-2 T-1 ~> m-2 s-1]
  type(RealArray_t), intent(inout) :: sh_xy_a !< horizontal shearing strain (du/dy + dv/dx) including metric
                                              !! terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xy_smooth_a !< horizontal shearing strain from smoothed velocity
                                                     !! including metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: str_xy_a !< str_xy is the cross term in the stress tensor
                                               !! [H L2 T-2 ~> m3 s-2 or kg s-2], but at some points in the
                                               !! code it is not yet layer integrated, so is in
                                               !! [L2 T-2 ~> m2 s-2].
  type(RealArray_t), intent(inout) :: bhstr_xy_a !< A copy of str_xy that only contains the biharmonic
                                                 !! contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
  type(RealArray_t), intent(inout) :: vort_xy_a !< Vertical vorticity (dv/dx - du/dy) including metric terms
                                                !! [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: grad_vort_mag_q_a !< Magnitude of vorticity gradient at q-points
                                                        !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: grad_vort_mag_q_2d_a !< Magnitude of 2d vorticity gradient at q-points
                                                           !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: Del2vort_q_a !< Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
  type(RealArray_t), intent(inout) :: grad_div_mag_q_a !< Magnitude of divergence gradient at q-points
                                                       !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: hq_a !< harmonic mean of the harmonic means of the u- & v point
                                           !! thicknesses [H ~> m or kg m-2] This form guarantees that hq/hu <
                                           !! 4.
  type(RealArray_t), intent(inout) :: Ah_a !< biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
  type(RealArray_t), intent(inout) :: Kh_a !< Laplacian  viscosity (h or q) [L2 T-1 ~> m2 s-1]
  type(RealArray_t), intent(inout) :: Shear_mag_a !< magnitude of the shear (h or q) [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: vert_vort_mag_a !< magnitude of the vertical vorticity gradient (h or q)
                                                      !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: hrat_min_a !< h_min divided by the thickness at the stress point (h or
                                                 !! q) [nondim]
  type(RealArray_t), intent(inout) :: visc_bound_rem_a !< fraction of overall viscous bounds that remain to be
                                                       !! applied (h or q) [nondim]

  ! Local variables
  real, dimension(:,:), contiguous, pointer :: GME_effic_q, GME_effic_h, sh_xy_bt, sh_xx_bt
  real, dimension(:,:), contiguous, pointer :: mask2dBu
  real, dimension(:,:,:), contiguous, pointer :: visc_bound_rem, hrat_min, vert_vort_mag, &
    Shear_mag, Kh, Ah, hq, grad_div_mag_q, Del2vort_q, grad_vort_mag_q_2d, grad_vort_mag_q, &
    vort_xy, bhstr_xy, str_xy, sh_xy_smooth, sh_xy, dDel2udy, dDel2vdx, str_xx, sh_xx, BS_coeff_q, &
    visc_limit_q_frac, visc_limit_q_flag, visc_limit_q, sh_xy_q, vort_xy_q, Kh_q, Kh_h, Ah_q, Ah_h, &
    str_xy_GME, str_xx_GME, GME_coeff_q, GME_coeff_h, KH_v_GME, KH_u_GME, v, u
  real, dimension(:,:), contiguous, pointer :: Kh_bg_xy, Ah_bg_xy, reduction_xy, Kh_Max_xy, Ah_Max_xy, &
    Ah_Max_xy_KS, n1n2_q, n1n1_m_n2n2_q, Laplac2_const_xy, Biharm6_const_xy, Laplac3_const_xy, &
    Biharm_const_xy, Biharm_const2_xy, Re_Ah_const_xy
  real, dimension(:,:), contiguous, pointer :: MEKE_Ku, MEKE_Au, Res_fn_q
  real, dimension(:,:,:), contiguous, pointer :: BS_struct
  real :: AhSm       ! Smagorinsky biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: AhLth      ! 2D Leith biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: Kh_max_here ! The local maximum Laplacian viscosity for stability [L2 T-1 ~> m2 s-1]
  real :: local_strain ! Local variable for interpolating computed strain rates [T-1 ~> s-1].
  real :: meke_res_fn ! A copy of the resolution scaling factor if being applied to MEKE [nondim]. Otherwise = 1.
  real :: GME_coeff ! The GME (negative) viscosity coefficient [L2 T-1 ~> m2 s-1]
  real :: KE        ! Local kinetic energy [L2 T-2 ~> m2 s-2]
  real :: d_str     ! Stress tensor update [L2 T-2 ~> m2 s-2]
  real :: grad_vort ! Vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
  real :: grad_vort_qg ! QG-based vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
  real :: tmp
  integer :: i, j, k, kk, is, ie, js, je, Isq, Ieq, Jsq, Jeq, kstart, kmax

  call visc_bound_rem_a%view(visc_bound_rem) ; call hrat_min_a%view(hrat_min)
  call vert_vort_mag_a%view(vert_vort_mag) ; call Shear_mag_a%view(Shear_mag)
  call Kh_a%view(Kh) ; call Ah_a%view(Ah) ; call hq_a%view(hq)
  call grad_div_mag_q_a%view(grad_div_mag_q) ; call Del2vort_q_a%view(Del2vort_q)
  call grad_vort_mag_q_2d_a%view(grad_vort_mag_q_2d)
  call grad_vort_mag_q_a%view(grad_vort_mag_q) ; call vort_xy_a%view(vort_xy)
  call bhstr_xy_a%view(bhstr_xy) ; call str_xy_a%view(str_xy)
  call sh_xy_smooth_a%view(sh_xy_smooth) ; call sh_xy_a%view(sh_xy)
  call dDel2udy_a%view(dDel2udy) ; call dDel2vdx_a%view(dDel2vdx) ; call str_xx_a%view(str_xx)
  call sh_xx_a%view(sh_xx) ; call BS_coeff_q_a%view(BS_coeff_q)
  call visc_limit_q_frac_a%view(visc_limit_q_frac)
  call visc_limit_q_flag_a%view(visc_limit_q_flag) ; call visc_limit_q_a%view(visc_limit_q)
  call sh_xy_q_a%view(sh_xy_q) ; call vort_xy_q_a%view(vort_xy_q) ; call Kh_q_a%view(Kh_q)
  call Kh_h_a%view(Kh_h) ; call Ah_q_a%view(Ah_q) ; call Ah_h_a%view(Ah_h)
  call str_xy_GME_a%view(str_xy_GME) ; call str_xx_GME_a%view(str_xx_GME)
  call GME_coeff_q_a%view(GME_coeff_q) ; call GME_coeff_h_a%view(GME_coeff_h)
  call KH_v_GME_a%view(KH_v_GME) ; call KH_u_GME_a%view(KH_u_GME)
  call GME_effic_q_a%view(GME_effic_q) ; call GME_effic_h_a%view(GME_effic_h)
  call sh_xy_bt_a%view(sh_xy_bt) ; call sh_xx_bt_a%view(sh_xx_bt) ; call v_a%view(v)
  call u_a%view(u)
  call Gcore%mask2dBu%view(mask2dBu)

  nullify(Kh_bg_xy, Ah_bg_xy, reduction_xy, Kh_Max_xy, Ah_Max_xy, Ah_Max_xy_KS, n1n2_q, n1n1_m_n2n2_q, &
          Laplac2_const_xy, Biharm6_const_xy, Laplac3_const_xy, Biharm_const_xy, Biharm_const2_xy, &
          Re_Ah_const_xy)
  if (metrics%Kh_bg_xy%associated()) call metrics%Kh_bg_xy%view(Kh_bg_xy)
  if (metrics%Ah_bg_xy%associated()) call metrics%Ah_bg_xy%view(Ah_bg_xy)
  if (metrics%reduction_xy%associated()) call metrics%reduction_xy%view(reduction_xy)
  if (metrics%Kh_Max_xy%associated()) call metrics%Kh_Max_xy%view(Kh_Max_xy)
  if (metrics%Ah_Max_xy%associated()) call metrics%Ah_Max_xy%view(Ah_Max_xy)
  if (metrics%Ah_Max_xy_KS%associated()) call metrics%Ah_Max_xy_KS%view(Ah_Max_xy_KS)
  if (metrics%n1n2_q%associated()) call metrics%n1n2_q%view(n1n2_q)
  if (metrics%n1n1_m_n2n2_q%associated()) call metrics%n1n1_m_n2n2_q%view(n1n1_m_n2n2_q)
  if (metrics%Laplac2_const_xy%associated()) call metrics%Laplac2_const_xy%view(Laplac2_const_xy)
  if (metrics%Biharm6_const_xy%associated()) call metrics%Biharm6_const_xy%view(Biharm6_const_xy)
  if (metrics%Laplac3_const_xy%associated()) call metrics%Laplac3_const_xy%view(Laplac3_const_xy)
  if (metrics%Biharm_const_xy%associated()) call metrics%Biharm_const_xy%view(Biharm_const_xy)
  if (metrics%Biharm_const2_xy%associated()) call metrics%Biharm_const2_xy%view(Biharm_const2_xy)
  if (metrics%Re_Ah_const_xy%associated()) call metrics%Re_Ah_const_xy%view(Re_Ah_const_xy)

  nullify(MEKE_Ku, MEKE_Au, Res_fn_q, BS_struct)
  if (MEKE_a%Ku_a%associated()) call MEKE_a%Ku_a%view(MEKE_Ku)
  if (MEKE_a%Au_a%associated()) call MEKE_a%Au_a%view(MEKE_Au)
  if (VarMix_a%Res_fn_q_a%associated()) call VarMix_a%Res_fn_q_a%view(Res_fn_q)
  if (VarMix_a%BS_struct_a%associated()) call VarMix_a%BS_struct_a%view(BS_struct)

  is  = bxH%idxS(1) ; ie  = bxH%idxE(1) ; js  = bxH%idxS(2) ; je  = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  kstart = bxH%idxS(3) ; kmax = bxH%idxE(3) - bxH%idxS(3) + 1

  !!$omp target update from(sh_xx, sh_xy)
  !!$omp target update from(h_u, h_v, hq)
  !!$omp target update from(str_xx)
  !!$omp target update from(Shear_mag) if (use_Smag)
  !!$omp target update from(Del2u, Del2v) if (CS%biharmonic)
  !!$omp target update from(dDel2vdx, dDel2udy) if (CS%biharmonic)

  if (CS%Laplacian) then
    ! Determine the Laplacian viscosity at q points, using the
    ! largest value from several parameterizations. Also get the
    ! Laplacian component of str_xy.

    if ((CS%Leith_Kh) .or. (CS%Leith_Ah)) then
      if (CS%use_QG_Leith_visc) then
        do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(grad_vort, grad_vort_qg))
          grad_vort = grad_vort_mag_q(I,J,kk) + grad_div_mag_q(I,J,kk)
          grad_vort_qg = 3. * grad_vort_mag_q_2d(I,J,kk)
          vert_vort_mag(I,J,kk) = min(grad_vort, grad_vort_qg)
        enddo
      else
        do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
          vert_vort_mag(I,J,kk) = grad_vort_mag_q(I,J,kk) + grad_div_mag_q(I,J,kk)
        enddo
      endif
    endif

    ! Static (pre-computed) background viscosity
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
      Kh(I,J,kk) = Kh_bg_xy(I,J)
    enddo

    if (CS%Smagorinsky_Kh) then
      if (CS%add_LES_viscosity) then
        do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
          Kh(I,J,kk) = Kh(I,J,kk) + Laplac2_const_xy(I,J) * Shear_mag(I,J,kk)
        enddo
      else
        do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
          Kh(I,J,kk) = max(Kh(I,J,kk), Laplac2_const_xy(I,J) * Shear_mag(I,J,kk) )
        enddo
      endif
    endif

    if (CS%Leith_Kh) then
      if (CS%add_LES_viscosity) then
        do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
          Kh(I,J,kk) = Kh(I,J,kk) + Laplac3_const_xy(I,J) * vert_vort_mag(I,J,kk) * inv_PI3 ! Is this right? -AJA
        enddo
      else
        do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
          Kh(I,J,kk) = max(Kh(I,J,kk), Laplac3_const_xy(I,J) * vert_vort_mag(I,J,kk) * inv_PI3)
        enddo
      endif
    endif

    ! All viscosity contributions above are subject to resolution scaling

    if (rescale_Kh) then
      do kk=1,kmax ; do J=js-1,Jeq ; do I=is-1,Ieq
        Kh(I,J,kk) = Res_fn_q(I,J) * Kh(I,J,kk)
      enddo ; enddo ; enddo
    endif

    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
      Kh(I,J,kk) = max(Kh(I,J,kk), CS%Kh_bg_min) ! Place a floor on the viscosity, if desired.
    enddo

    if (use_MEKE_Ku .and. .not. CS%EY24_EBT_BS) then
      if (use_kh_struct) then
        do kk=1,kmax ; do J=js-1,Jeq ; do I=is-1,Ieq
          k = kstart + kk - 1
          meke_res_fn = 1.
          if (CS%res_scale_MEKE) meke_res_fn = Res_fn_q(I,J)

          Kh(I,J,kk) = Kh(I,J,kk) + 0.25*( ((MEKE_Ku(i,j)*BS_struct(i,j,k)) + &
                                     (MEKE_Ku(i+1,j+1)*BS_struct(i+1,j+1,k))) + &
                                     ((MEKE_Ku(i+1,j)*BS_struct(i+1,j,k)) + &
                                     (MEKE_Ku(i,j+1)*BS_struct(i,j+1,k))) ) * meke_res_fn
        enddo ; enddo ; enddo
      else
        do kk=1,kmax ; do J=js-1,Jeq ; do I=is-1,Ieq
          meke_res_fn = 1.
          if (CS%res_scale_MEKE) meke_res_fn = Res_fn_q(I,J)

          Kh(I,J,kk) = Kh(I,J,kk) + 0.25 * ( &
              (MEKE_Ku(i,j) + MEKE_Ku(i+1,j+1)) + &
                                     (MEKE_Ku(i+1,j) + &
                                      MEKE_Ku(i,j+1)) ) * meke_res_fn
        enddo ; enddo ; enddo
      endif
    endif

    if (CS%anisotropic) then
      ! *Add* the shear component of anisotropic viscosity
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
          Kh(I,J,kk) = Kh(I,J,kk) + CS%Kh_aniso * n1n2_q(I,J)**2
      enddo
    endif

    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(Kh_max_here))
      ! Newer method of bounding for stability
      if ((CS%bound_Kh) .and. (CS%bound_Ah)) then
        visc_bound_rem(I,J,kk) = 1.0
        Kh_max_here = hrat_min(I,J,kk) * Kh_Max_xy(I,J)
        if (Kh(I,J,kk) >= Kh_max_here) then
          visc_bound_rem(I,J,kk) = 0.0
          Kh(I,J,kk) = Kh_max_here
        elseif ((Kh(I,J,kk) > 0.0) .or. (CS%backscatter_underbound .and. (Kh_max_here > 0.0))) then
          visc_bound_rem(I,J,kk) = 1.0 - Kh(I,J,kk) / Kh_max_here
        endif
      elseif (CS%bound_Kh) then
        Kh(I,J,kk) = min(Kh(I,J,kk), hrat_min(I,J,kk) * Kh_Max_xy(I,J))
      endif
    enddo

    if (CS%use_Leithy) then
      ! Leith+E doesn't recompute Kh at q points, it just interpolates it from h to q points
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        Kh(I,J,kk) = 0.25 * ((Kh_h(i,j,k) + Kh_h(i+1,j+1,k)) + (Kh_h(i,j+1,k) + Kh_h(i+1,j,k)))
      enddo
    endif

    if (CS%do_Kh_q .or. CS%debug) then
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        Kh_q(I,J,k) = Kh(I,J,kk)
      enddo
    endif

    if (CS%do_vort_xy_q) then
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        vort_xy_q(I,J,k) = vort_xy(I,J,kk)
      enddo
    endif

    if (CS%do_sh_xy_q) then
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        sh_xy_q(I,J,k) = sh_xy(I,J,kk)
      enddo
    endif

    if (.not. CS%use_Leithy) then
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
        str_xy(I,J,kk) = -Kh(I,J,kk) * sh_xy(I,J,kk)
      enddo
    else
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
        str_xy(I,J,kk) = -Kh(I,J,kk) * sh_xy_smooth(I,J,kk)
      enddo
    endif
  else
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
      str_xy(I,J,kk) = 0.
    enddo
  endif ! get harmonic coefficient Kh at q points and harmonic part of str_xy

  if (CS%anisotropic) then
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(local_strain))
      ! Horizontal-tension averaged to q-points
      local_strain = 0.25 * ( (sh_xx(i,j,kk) + sh_xx(i+1,j+1,kk)) + (sh_xx(i+1,j,kk) + sh_xx(i,j+1,kk)) )
      ! *Add* the tension contribution to the xy-component of stress
      str_xy(I,J,kk) = str_xy(I,J,kk) - CS%Kh_aniso * n1n2_q(I,J) * n1n1_m_n2n2_q(I,J) * local_strain
    enddo
  endif

  if (CS%biharmonic) then
    ! Determine the biharmonic viscosity at q points, using the
    ! largest value from several parameterizations. Also get the
    ! biharmonic component of str_xy.
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
      Ah(I,J,kk) = Ah_bg_xy(I,J)
    enddo

    if (CS%Smagorinsky_Ah .or. CS%Leith_Ah) then
      if (CS%Smagorinsky_Ah) then
        if (CS%bound_Coriolis) then
          do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(AhSm))
            AhSm = Shear_mag(I,J,kk) * (Biharm_const_xy(I,J) &
                + Biharm_const2_xy(I,J) * Shear_mag(I,J,kk))
            Ah(I,J,kk) = max(Ah(I,J,kk), AhSm)
          enddo
        else
          do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(AhSm))
            AhSm = Biharm_const_xy(I,J) * Shear_mag(I,J,kk)
            Ah(I,J,kk) = max(Ah(I,J,kk), AhSm)
          enddo
        endif
      endif

      if (CS%Leith_Ah) then
        do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(AhLth))
          AhLth = Biharm6_const_xy(I,J) * abs(Del2vort_q(I,J,kk)) * inv_PI6
          Ah(I,J,kk) = max(Ah(I,J,kk), AhLth)
        enddo
      endif
    endif ! Smagorinsky_Ah or Leith_Ah

    if (use_MEKE_Au) then
      ! *Add* the MEKE contribution
      do kk=1,kmax ; do J=js-1,Jeq ; do I=is-1,Ieq
        Ah(I,J,kk) = Ah(I,J,kk) + 0.25 * ( &
            (MEKE_Au(i,j) + MEKE_Au(i+1,j+1)) + (MEKE_Au(i+1,j) + MEKE_Au(i,j+1)) )
      enddo ; enddo ; enddo
    endif

    ! XXX: It is just overwrites the values!
    if (CS%Re_Ah > 0.0) then
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(k, KE))
        k = kstart + kk - 1
        KE = 0.125 * (((u(I,j,k) + u(I,j+1,k))**2) + ((v(i,J,k) + v(i+1,J,k))**2))
        Ah(I,J,kk) = sqrt(KE) * Re_Ah_const_xy(I,J)
      enddo
    endif

    if (CS%bound_Ah) then
      if (CS%bound_Kh) then
        do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
          Ah(I,J,kk) = min(Ah(I,J,kk), visc_bound_rem(I,J,kk) * hrat_min(I,J,kk) * Ah_Max_xy(I,J))
        enddo
      else
        do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
          Ah(I,J,kk) = min(Ah(I,J,kk), hrat_min(I,J,kk) * Ah_Max_xy(I,J))
        enddo
      endif
    endif

    if (CS%EY24_EBT_BS) then
        do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(k, tmp))
          k = kstart + kk - 1
          tmp = CS%KS_coef *hrat_min(I,J,kk) * Ah_Max_xy_KS(I,J)
          visc_limit_q(I,J,k) = tmp
          visc_limit_q_frac(i,j,k) = Ah(i,j,kk) / (CS%KS_coef * hrat_min(i,j,kk) * Ah_Max_xy_KS(i,j))
          if (Ah(I,J,kk) >= tmp) then
            visc_limit_q_flag(I,J,k) = 1.
          endif
        enddo
    endif

    ! Leith+E doesn't recompute Ah at q points, it just interpolates it from h to q points
    if (CS%use_Leithy) then
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        Ah(I,J,kk) = 0.25 * ((Ah_h(i,j,k) + Ah_h(i+1,j+1,k)) + (Ah_h(i,j+1,k) + Ah_h(i+1,j,k)))
      enddo
    endif

    if (CS%do_Ah_q .or. CS%debug) then
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        Ah_q(I,J,k) = Ah(I,J,kk)
      enddo
    endif

    ! Again, need to initialize str_xy as if its biharmonic
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(d_str))
      d_str = Ah(I,J,kk) * (dDel2vdx(I,J,kk) + dDel2udy(I,J,kk))

      str_xy(I,J,kk) = str_xy(I,J,kk) + d_str

      ! Keep a copy of the biharmonic contribution for backscatter parameterization
      bhstr_xy(I,J,kk) = d_str * (hq(I,J,kk) * mask2dBu(I,J) * reduction_xy(I,J))
    enddo
  endif ! Get Ah at q points and biharmonic part of str_xy

  ! Backscatter using MEKE
  if (CS%EY24_EBT_BS) then
    call hor_visc_backscatter_q(CS, MEKE_a, VarMix_a, use_kh_struct, nkblock, kstart, kmax, is, js, &
                                Ieq, Jeq, visc_limit_q_flag_a, sh_xy_a, str_xy_a, BS_coeff_q_a)
  endif ! Backscatter

  if (CS%use_GME) then

    ! The wider halo here is to permit one pass of smoothing without a halo update.
    do concurrent (kk=1:kmax, j=Jsq-1:Jeq+2, i=Isq-1:Ieq+2) DO_LOCALITY(local(k, GME_coeff))
      k = kstart + kk - 1
      GME_coeff = GME_effic_h(i,j) * 0.25 * &
          ((KH_u_GME(I,j,k)+KH_u_GME(I-1,j,k)) + (KH_v_GME(i,J,k)+KH_v_GME(i,J-1,k)))
      GME_coeff = MIN(GME_coeff, CS%GME_limiter)

      if ((CS%do_GME_coeff_h) .or. find_FrictWork) GME_coeff_h(i,j,k) = GME_coeff
      str_xx_GME(i,j,kk) = GME_coeff * sh_xx_bt(i,j)
    enddo

    ! The wider halo here is to permit one pass of smoothing without a halo update.
    do concurrent (kk=1:kmax, J=js-2:je+1, I=is-2:ie+1) DO_LOCALITY(local(k, GME_coeff))
      k = kstart + kk - 1
      GME_coeff = GME_effic_q(I,J) * 0.25 * &
          ((KH_u_GME(I,j,k)+KH_u_GME(I,j+1,k)) + (KH_v_GME(i,J,k)+KH_v_GME(i+1,J,k)))
      GME_coeff = MIN(GME_coeff, CS%GME_limiter)

      if (CS%do_GME_coeff_q) GME_coeff_q(I,J,k) = GME_coeff
      str_xy_GME(I,J,kk) = GME_coeff * sh_xy_bt(I,J)
    enddo
  endif

end subroutine hor_visc_q_coeffs

!> Converts the stresses to layer-integrated form, takes their divergence to get the
!! viscous accelerations, and accumulates the frictional work for one tile of
!! horizontal_viscosity.
subroutine hor_visc_stress_div(MEKE_a, Gcore, CS, metrics, OBC, STOCH_a, bxH, bxQ, nkblock, &
                               h_neglect, apply_OBC, find_FrictWork, find_FrictWork_bh, &
                               backscat_subround, skeb_use_frict, u_a, v_a, h_a, uh_a, vh_a, &
                               diffu_a, diffv_a, str_xx_GME_a, str_xy_GME_a, FrictWork_a, &
                               FrictWork_bh_a, FrictWork_GME_a, visc_limit_h_flag_a, h_u_a, h_v_a, &
                               sh_xx_a, str_xx_a, bhstr_xx_a, sh_xy_a, str_xy_a, bhstr_xy_a, hq_a, &
                               H_to_RZ, s_to_T)
  type(hor_visc_MEKE_type),      intent(inout) :: MEKE_a !< MEKE fields
                                                       !! related to Mesoscale Eddy Kinetic Energy.
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  real,                 intent(in)    :: H_to_RZ !< A constant that translates thickness units to the units
                                                 !! of mass per unit area [R Z H-1 ~> kg m-3 or 1]
  real,                 intent(in)    :: s_to_T  !< A constant that translates time intervals in seconds to
                                                 !! the units of time [T s-1 ~> 1]
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  type(hv_metrics_CS), intent(in)    :: metrics       !< Persistent copies of the horizontal
                                                     !! viscosity metric and coefficient arrays
  type(ocean_OBC_type), optional, pointer      :: OBC  !< Pointer to an open boundary condition type
  type(hor_visc_STOCH_type), intent(inout)           :: STOCH_a !< Stochastic control structure
  type(Box_t),                   intent(in)    :: bxH  !< The h-point iteration box of this tile,
                                                       !! [isc:iec, jsc:jec, kstart:kend]
  type(Box_t),                   intent(in)    :: bxQ  !< The B-grid-index iteration box of this
                                                       !! tile, [IscB:IecB, JscB:JecB, kstart:kend]
  integer,                       intent(in)    :: nkblock !< The effective k-block size [nondim]
  real,    intent(in)    :: h_neglect !< thickness so small it can be lost in roundoff and so neglected
                                      !! [H ~> m or kg m-2]
  logical, intent(in)    :: apply_OBC !< If true, open boundary conditions are applied on this PE
  logical, intent(in)    :: find_FrictWork !< If true, the work done by lateral friction is calculated
  logical, intent(in)    :: find_FrictWork_bh !< If true, the work done by biharmonic lateral friction is calculated
  real,    intent(in)    :: backscat_subround !< The ratio of f over Shear_mag that is so small that the backscatter
                                              !! calculation gives the same value as if f were 0 [nondim].
  logical, intent(in)    :: skeb_use_frict !< If true, accumulate the frictional dissipation for stochastic backscatter
  type(RealArray_t), intent(in)    :: u_a !< The zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: v_a !< The meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t), intent(in)    :: h_a !< Layer thicknesses [H ~> m or kg m-2].
  type(RealArray_t), intent(in)    :: uh_a !< The zonal volume transport [H L2 T-1 ~> m3 s-1].
  type(RealArray_t), intent(in)    :: vh_a !< The meridional volume transport [H L2 T-1 ~> m3 s-1].
  type(RealArray_t), intent(inout) :: diffu_a !< Zonal acceleration due to convergence of along-coordinate
                                              !! stress tensor [L T-2 ~> m s-2]
  type(RealArray_t), intent(inout) :: diffv_a !< Meridional acceleration due to convergence of
                                              !! along-coordinate stress tensor [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: str_xx_GME_a !< smoothed diagonal term in the stress tensor from GME
                                                   !! [L2 T-2 ~> m2 s-2]
  type(RealArray_t), intent(inout) :: str_xy_GME_a !< smoothed cross term in the stress tensor from GME
                                                   !! [L2 T-2 ~> m2 s-2]
  type(RealArray_t), intent(inout) :: FrictWork_a !< work done by MKE dissipation mechanisms
                                                  !! [R Z L2 T-3 ~> W m-2]
  type(RealArray_t), intent(inout) :: FrictWork_bh_a !< work done by the biharmonic MKE dissipation mechanisms
                                                     !! [R Z L2 T-3 ~> W m-2]
  type(RealArray_t), intent(inout) :: FrictWork_GME_a !< work done by GME [R Z L2 T-3 ~> W m-2]
  type(RealArray_t), intent(inout) :: visc_limit_h_flag_a !< determines whether backscatter is shut off
                                                          !! [nondim]
  type(RealArray_t), intent(inout) :: h_u_a !< Thickness interpolated to u points [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: h_v_a !< Thickness interpolated to v points [H ~> m or kg m-2].
  type(RealArray_t), intent(inout) :: sh_xx_a !< horizontal tension (du/dx - dv/dy) including metric terms
                                              !! [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: str_xx_a !< str_xx is the diagonal term in the stress tensor
                                               !! [H L2 T-2 ~> m3 s-2 or kg s-2], but at some points in the
                                               !! code it is not yet layer integrated, so is in
                                               !! [L2 T-2 ~> m2 s-2].
  type(RealArray_t), intent(inout) :: bhstr_xx_a !< A copy of str_xx that only contains the biharmonic
                                                 !! contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
  type(RealArray_t), intent(inout) :: sh_xy_a !< horizontal shearing strain (du/dy + dv/dx) including metric
                                              !! terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: str_xy_a !< str_xy is the cross term in the stress tensor
                                               !! [H L2 T-2 ~> m3 s-2 or kg s-2], but at some points in the
                                               !! code it is not yet layer integrated, so is in
                                               !! [L2 T-2 ~> m2 s-2].
  type(RealArray_t), intent(inout) :: bhstr_xy_a !< A copy of str_xy that only contains the biharmonic
                                                 !! contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
  type(RealArray_t), intent(inout) :: hq_a !< harmonic mean of the harmonic means of the u- & v point
                                           !! thicknesses [H ~> m or kg m-2] This form guarantees that hq/hu <
                                           !! 4.

  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: hq, bhstr_xy, str_xy, sh_xy, bhstr_xx, str_xx, &
    sh_xx, h_v, h_u, visc_limit_h_flag, FrictWork_GME, FrictWork_bh, FrictWork, str_xy_GME, &
    str_xx_GME, diffv, diffu, vh, uh, h, v, u
  real, dimension(:,:), contiguous, pointer :: mask2dBu, IdxCu, IdyCu, IareaCu, IdyCv, IdxCv, &
    IareaCv, IdxT, IdyT, IdyBu, IdxBu, IareaT, dxCu, dyCv, CoriolisBu
  real, dimension(:,:), contiguous, pointer :: reduction_xx, reduction_xy, dx2h, dy2h, dx2q, dy2q
  real, dimension(:,:), contiguous, pointer :: MEKE_mom_src, MEKE_mom_src_bh, MEKE_GME_snk
  real, dimension(:,:,:), contiguous, pointer :: skeb_diss
  real :: Shear_mag_bc  ! Shear_mag value in backscatter [T-1 ~> s-1]
  real :: RoScl     ! The scaling function for MEKE source term [nondim]
  real :: FatH      ! abs(f) at h-point for MEKE source term [T-1 ~> s-1]
  real :: Sh_F_pow  ! The ratio of shear over the absolute value of f raised to some power and rescaled [nondim]
  integer :: i, j, k, kk, n, is, ie, js, je, Isq, Ieq, Jsq, Jeq, kstart, kmax

  call hq_a%view(hq) ; call bhstr_xy_a%view(bhstr_xy) ; call str_xy_a%view(str_xy)
  call sh_xy_a%view(sh_xy) ; call bhstr_xx_a%view(bhstr_xx) ; call str_xx_a%view(str_xx)
  call sh_xx_a%view(sh_xx) ; call h_v_a%view(h_v) ; call h_u_a%view(h_u)
  call visc_limit_h_flag_a%view(visc_limit_h_flag) ; call FrictWork_GME_a%view(FrictWork_GME)
  call FrictWork_bh_a%view(FrictWork_bh) ; call FrictWork_a%view(FrictWork)
  call str_xy_GME_a%view(str_xy_GME) ; call str_xx_GME_a%view(str_xx_GME)
  call diffv_a%view(diffv) ; call diffu_a%view(diffu) ; call vh_a%view(vh) ; call uh_a%view(uh)
  call h_a%view(h) ; call v_a%view(v) ; call u_a%view(u)
  call Gcore%mask2dBu%view(mask2dBu) ; call Gcore%IdxCu%view(IdxCu) ; call Gcore%IdyCu%view(IdyCu)
  call Gcore%IareaCu%view(IareaCu) ; call Gcore%IdyCv%view(IdyCv) ; call Gcore%IdxCv%view(IdxCv)
  call Gcore%IareaCv%view(IareaCv) ; call Gcore%IdxT%view(IdxT) ; call Gcore%IdyT%view(IdyT)
  call Gcore%IdyBu%view(IdyBu) ; call Gcore%IdxBu%view(IdxBu) ; call Gcore%IareaT%view(IareaT)
  call Gcore%dxCu%view(dxCu) ; call Gcore%dyCv%view(dyCv) ; call Gcore%CoriolisBu%view(CoriolisBu)

  nullify(reduction_xx, reduction_xy, dx2h, dy2h, dx2q, dy2q)
  if (metrics%reduction_xx%associated()) call metrics%reduction_xx%view(reduction_xx)
  if (metrics%reduction_xy%associated()) call metrics%reduction_xy%view(reduction_xy)
  if (metrics%dx2h%associated()) call metrics%dx2h%view(dx2h)
  if (metrics%dy2h%associated()) call metrics%dy2h%view(dy2h)
  if (metrics%dx2q%associated()) call metrics%dx2q%view(dx2q)
  if (metrics%dy2q%associated()) call metrics%dy2q%view(dy2q)

  nullify(MEKE_mom_src, MEKE_mom_src_bh, MEKE_GME_snk, skeb_diss)
  if (MEKE_a%mom_src_a%associated()) call MEKE_a%mom_src_a%view(MEKE_mom_src)
  if (MEKE_a%mom_src_bh_a%associated()) call MEKE_a%mom_src_bh_a%view(MEKE_mom_src_bh)
  if (MEKE_a%GME_snk_a%associated()) call MEKE_a%GME_snk_a%view(MEKE_GME_snk)
  if (STOCH_a%skeb_diss_a%associated()) call STOCH_a%skeb_diss_a%view(skeb_diss)

  is  = bxH%idxS(1) ; ie  = bxH%idxE(1) ; js  = bxH%idxS(2) ; je  = bxH%idxE(2)
  Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
  kstart = bxH%idxS(3) ; kmax = bxH%idxE(3) - bxH%idxS(3) + 1

  if (CS%use_GME) then

    ! This changes the units of str_xx from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
    do concurrent (kk=1:kmax, j=Jsq:Jeq+1, i=Isq:Ieq+1) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      str_xx(i,j,kk) = (str_xx(i,j,kk) + str_xx_GME(i,j,kk)) * (h(i,j,k) * reduction_xx(i,j))
    enddo

    ! This adds in GME and changes the units of str_xx from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
    if (CS%no_slip) then
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
        str_xy(I,J,kk) = (str_xy(I,J,kk) + str_xy_GME(I,J,kk)) * (hq(I,J,kk) * reduction_xy(I,J))
      enddo
    else
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
        str_xy(I,J,kk) = (str_xy(I,J,kk) + str_xy_GME(I,J,kk)) * (hq(I,J,kk) * mask2dBu(I,J) * reduction_xy(I,J))
      enddo
    endif
  else ! .not. use_GME
    ! This changes the units of str_xx from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
    do concurrent (kk=1:kmax, j=Jsq:Jeq+1, i=Isq:Ieq+1) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      str_xx(i,j,kk) = str_xx(i,j,kk) * (h(i,j,k) * reduction_xx(i,j))
    enddo

    ! This changes the units of str_xy from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
    if (CS%no_slip) then
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
        str_xy(I,J,kk) = str_xy(I,J,kk) * (hq(I,J,kk) * reduction_xy(I,J))
      enddo
    else
      do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
        str_xy(I,J,kk) = str_xy(I,J,kk) * (hq(I,J,kk) * mask2dBu(I,J) * reduction_xy(I,J))
      enddo
    endif
  endif ! use_GME

  ! Evaluate 1/h x.Div(h Grad u) or the biharmonic equivalent.
  do concurrent (kk=1:kmax, j=js:je, I=Isq:Ieq) DO_LOCALITY(local(k))
    k = kstart + kk - 1
    diffu(I,j,k) = ((IdxCu(I,j)*((dx2q(I,J-1)*str_xy(I,J-1,kk)) - (dx2q(I,J)*str_xy(I,J,kk))) + &
                     IdyCu(I,j)*((dy2h(i,j)*str_xx(i,j,kk)) - (dy2h(i+1,j)*str_xx(i+1,j,kk)))) * &
                   IareaCu(I,j)) / (h_u(I,j,kk) + h_neglect)
  enddo

  if (apply_OBC) then
    ! This is not the right boundary condition. If all the masking of tendencies are done
    ! correctly later then eliminating this block should not change answers.
    do n=1,OBC%number_of_segments
      if (OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        do concurrent (kk=1:kmax, j=OBC%segment(n)%HI%jsd:OBC%segment(n)%HI%jed) DO_LOCALITY(local(k))
          k = kstart + kk - 1
          diffu(I,j,k) = 0.
        enddo
      endif
    enddo
  endif

  ! Evaluate 1/h y.Div(h Grad u) or the biharmonic equivalent.
  do concurrent (kk=1:kmax, J=Jsq:Jeq, i=is:ie) DO_LOCALITY(local(k))
    k = kstart + kk - 1
    diffv(i,J,k) = ((IdyCv(i,J)*((dy2q(I-1,J)*str_xy(I-1,J,kk)) - (dy2q(I,J)*str_xy(I,J,kk))) - &
                     IdxCv(i,J)*((dx2h(i,j)*str_xx(i,j,kk)) - (dx2h(i,j+1)*str_xx(i,j+1,kk)))) * &
                   IareaCv(i,J)) / (h_v(i,J,kk) + h_neglect)
  enddo

  if (apply_OBC) then
    ! This is not the right boundary condition. If all the masking of tendencies are done
    ! correctly later then eliminating this block should not change answers.
    do n=1,OBC%number_of_segments
      if (OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        do concurrent (kk=1:kmax, i=OBC%segment(n)%HI%isd:OBC%segment(n)%HI%ied) DO_LOCALITY(local(k))
          k = kstart + kk - 1
          diffv(i,J,k) = 0.
        enddo
      endif
    enddo
  endif


  if (find_FrictWork) then

    if (CS%FrictWork_bug) then
      ! Diagnose   str_xx*d_x u - str_yy*d_y v + str_xy*(d_y u + d_x v)
      ! This is the old formulation that includes energy diffusion
      do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        FrictWork(i,j,k) = H_to_RZ * ( &
                ((str_xx(i,j,kk) * (u(I,j,k)-u(I-1,j,k))*IdxT(i,j))    &
               - (str_xx(i,j,kk) * (v(i,J,k)-v(i,J-1,k))*IdyT(i,j)))   &
            + 0.25*(( (str_xy(I,J,kk) *                                  &
                       (((u(I,j+1,k)-u(I,j,k))*IdyBu(I,J))          &
                      + ((v(i+1,J,k)-v(i,J,k))*IdxBu(I,J))))        &
                    + (str_xy(I-1,J-1,kk) *                              &
                       (((u(I-1,j,k)-u(I-1,j-1,k))*IdyBu(I-1,J-1))  &
                      + ((v(i,J-1,k)-v(i-1,J-1,k))*IdxBu(I-1,J-1)))) ) &
                  + ( (str_xy(I-1,J,kk) *                                &
                       (((u(I-1,j+1,k)-u(I-1,j,k))*IdyBu(I-1,J))    &
                      + ((v(i,J,k)-v(i-1,J,k))*IdxBu(I-1,J))))      &
                    + (str_xy(I,J-1,kk) *                                &
                       (((u(I,j,k)-u(I,j-1,k))*IdyBu(I,J-1))        &
                      + ((v(i+1,J-1,k)-v(i,J-1,k))*IdxBu(I,J-1)))) ) ) )
      enddo
    else
      do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        FrictWork(i,j,k) = H_to_RZ * IareaT(i,j) * ( &
          ((str_xx(i,j,kk)*dy2h(i,j) * ( &
                (uh(I,j,k)*dxCu(I,j)*IdyCu(I,j)*IareaCu(I,j)/(h_u(I,j,kk)+h_neglect)) &
              - (uh(I-1,j,k)*dxCu(I-1,j)*IdyCu(I-1,j)*IareaCu(I-1,j)/(h_u(I-1,j,kk)+h_neglect)) ) ) &
         - (str_xx(i,j,kk)*dx2h(i,j) * ( &
                (vh(i,J,k)*dyCv(i,J)*IdxCv(i,J)*IareaCv(i,J)/(h_v(i,J,kk)+h_neglect)) &
              - (vh(i,J-1,k)*dyCv(i,J-1)*IdxCv(i,J-1)*IareaCv(i,J-1)/(h_v(i,J-1,kk)+h_neglect)) ) )) &
        + (0.25*(((str_xy(I,J,kk)*(                                     &
                   (dx2q(I,J)*((uh(I,j+1,k)*IareaCu(I,j+1)/(h_u(I,j+1,kk)+h_neglect)) &
                                - (uh(I,j,k)*IareaCu(I,j)/(h_u(I,j,kk)+h_neglect))))            &
                 + (dy2q(I,J)*((vh(i+1,J,k)*IareaCv(i+1,J)/(h_v(i+1,J,kk)+h_neglect)) &
                                - (vh(i,J,k)*IareaCv(i,J)/(h_v(i,J,kk)+h_neglect)))) ))          &
              +(str_xy(I-1,J-1,kk)*(                                 &
                   (dx2q(I-1,J-1)*((uh(I-1,j,k)*IareaCu(I-1,j)/(h_u(I-1,j,kk)+h_neglect)) &
                                    - (uh(I-1,j-1,k)*IareaCu(I-1,j-1)/(h_u(I-1,j-1,kk)+h_neglect))))    &
                 + (dy2q(I-1,J-1)*((vh(i,J-1,k)*IareaCv(i,J-1)/(h_v(i,J-1,kk)+h_neglect)) &
                                    - (vh(i-1,J-1,k)*IareaCv(i-1,J-1)/(h_v(i-1,J-1,kk)+h_neglect)))) )) ) &
             +((str_xy(I-1,J,kk)*(                                   &
                   (dx2q(I-1,J)*((uh(I-1,j+1,k)*IareaCu(I-1,j+1)/(h_u(I-1,j+1,kk)+h_neglect)) &
                                  - (uh(I-1,j,k)*IareaCu(I-1,j)/(h_u(I-1,j,kk)+h_neglect))))      &
                 + (dy2q(I-1,J)*((vh(i,J,k)*IareaCv(i,J)/(h_v(i,J,kk)+h_neglect)) &
                                  - (vh(i-1,J,k)*IareaCv(i-1,J)/(h_v(i-1,J,kk)+h_neglect)))) ))        &
              +(str_xy(I,J-1,kk)*(                                   &
                   (dx2q(I,J-1)*((uh(I,j,k)*IareaCu(I,j)/(h_u(I,j,kk)+h_neglect)) &
                                  - (uh(I,j-1,k)*IareaCu(I,j-1)/(h_u(I,j-1,kk)+h_neglect))))          &
                 + (dy2q(I,J-1)*((vh(i+1,J-1,k)*IareaCv(i+1,J-1)/(h_v(i+1,J-1,kk)+h_neglect)) &
                                  - (vh(i,J-1,k)*IareaCv(i,J-1)/(h_v(i,J-1,kk)+h_neglect)))) )) ) )) )

      enddo
    endif

    if (CS%EY24_EBT_BS) then
      do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        FrictWork(i,j,k) = (1. - visc_limit_h_flag(i,j,k)) * FrictWork(i,j,k)
      enddo
    endif
  endif

  if (find_FrictWork_bh) then
    if (CS%FrictWork_bug) then
      ! Diagnose   bhstr_xx*d_x u - bhstr_yy*d_y v + bhstr_xy*(d_y u + d_x v)
      ! This is the old formulation that includes energy diffusion !cyc
      do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        FrictWork_bh(i,j,k) = H_to_RZ * ( &
                ((bhstr_xx(i,j,kk) * (u(I,j,k)-u(I-1,j,k))*IdxT(i,j))  &
               - (bhstr_xx(i,j,kk) * (v(i,J,k)-v(i,J-1,k))*IdyT(i,j))) &
            + 0.25*(( (bhstr_xy(I,J,kk) *                              &
                     (((u(I,j+1,k)-u(I,j,k))*IdyBu(I,J))          &
                    + ((v(i+1,J,k)-v(i,J,k))*IdxBu(I,J))))        &
                  + (bhstr_xy(I-1,J-1,kk) *                            &
                     (((u(I-1,j,k)-u(I-1,j-1,k))*IdyBu(I-1,J-1))  &
                    + ((v(i,J-1,k)-v(i-1,J-1,k))*IdxBu(I-1,J-1)))) ) &
                + ( (bhstr_xy(I-1,J,kk) *                              &
                     (((u(I-1,j+1,k)-u(I-1,j,k))*IdyBu(I-1,J))    &
                    + ((v(i,J,k)-v(i-1,J,k))*IdxBu(I-1,J))))      &
                  + (bhstr_xy(I,J-1,kk) *                              &
                     (((u(I,j,k)-u(I,j-1,k))*IdyBu(I,J-1))        &
                    + ((v(i+1,J-1,k)-v(i,J-1,k))*IdxBu(I,J-1)))) ) ) )
      enddo
    else
      do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        ! Diagnose   bhstr_xx*d_x u - bhstr_yy*d_y v + bhstr_xy*(d_y u + d_x v)
        FrictWork_bh(i,j,k) = H_to_RZ * IareaT(i,j) * ( &
          ((bhstr_xx(i,j,kk)*dy2h(i,j) * ( &
                (uh(I,j,k)*dxCu(I,j)*IdyCu(I,j)*IareaCu(I,j)/(h_u(I,j,kk)+h_neglect)) &
              - (uh(I-1,j,k)*dxCu(I-1,j)*IdyCu(I-1,j)*IareaCu(I-1,j)/(h_u(I-1,j,kk)+h_neglect)) ) ) &
         - (bhstr_xx(i,j,kk)*dx2h(i,j) * ( &
                (vh(i,J,k)*dyCv(i,J)*IdxCv(i,J)*IareaCv(i,J)/(h_v(i,J,kk)+h_neglect)) &
              - (vh(i,J-1,k)*dyCv(i,J-1)*IdxCv(i,J-1)*IareaCv(i,J-1)/(h_v(i,J-1,kk)+h_neglect)) ) )) &
        + (0.25*(((bhstr_xy(I,J,kk)*(                                     &
                   (dx2q(I,J)*((uh(I,j+1,k)*IareaCu(I,j+1)/(h_u(I,j+1,kk)+h_neglect)) &
                                - (uh(I,j,k)*IareaCu(I,j)/(h_u(I,j,kk)+h_neglect))))            &
                 + (dy2q(I,J)*((vh(i+1,J,k)*IareaCv(i+1,J)/(h_v(i+1,J,kk)+h_neglect)) &
                                - (vh(i,J,k)*IareaCv(i,J)/(h_v(i,J,kk)+h_neglect)))) ))          &
              +(bhstr_xy(I-1,J-1,kk)*(                                 &
                   (dx2q(I-1,J-1)*((uh(I-1,j,k)*IareaCu(I-1,j)/(h_u(I-1,j,kk)+h_neglect)) &
                                    - (uh(I-1,j-1,k)*IareaCu(I-1,j-1)/(h_u(I-1,j-1,kk)+h_neglect))))    &
                 + (dy2q(I-1,J-1)*((vh(i,J-1,k)*IareaCv(i,J-1)/(h_v(i,J-1,kk)+h_neglect)) &
                                    - (vh(i-1,J-1,k)*IareaCv(i-1,J-1)/(h_v(i-1,J-1,kk)+h_neglect)))) )) ) &
             +((bhstr_xy(I-1,J,kk)*(                                   &
                   (dx2q(I-1,J)*((uh(I-1,j+1,k)*IareaCu(I-1,j+1)/(h_u(I-1,j+1,kk)+h_neglect)) &
                                  - (uh(I-1,j,k)*IareaCu(I-1,j)/(h_u(I-1,j,kk)+h_neglect))))      &
                 + (dy2q(I-1,J)*((vh(i,J,k)*IareaCv(i,J)/(h_v(i,J,kk)+h_neglect)) &
                                  - (vh(i-1,J,k)*IareaCv(i-1,J)/(h_v(i-1,J,kk)+h_neglect)))) ))        &
              +(bhstr_xy(I,J-1,kk)*(                                   &
                   (dx2q(I,J-1)*((uh(I,j,k)*IareaCu(I,j)/(h_u(I,j,kk)+h_neglect)) &
                                  - (uh(I,j-1,k)*IareaCu(I,j-1)/(h_u(I,j-1,kk)+h_neglect))))          &
                 + (dy2q(I,J-1)*((vh(i+1,J-1,k)*IareaCv(i+1,J-1)/(h_v(i+1,J-1,kk)+h_neglect)) &
                                  - (vh(i,J-1,k)*IareaCv(i,J-1)/(h_v(i,J-1,kk)+h_neglect)))) )) ) )) )
      enddo
    endif

    if (CS%EY24_EBT_BS) then
      do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
        k = kstart + kk - 1
        FrictWork_bh(i,j,k) = (1. - visc_limit_h_flag(i,j,k)) * FrictWork_bh(i,j,k)
      enddo
    endif
  endif

  if (CS%use_GME) then
    if (CS%FrictWork_bug) then ; do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
      k = kstart + kk - 1
    ! Diagnose   str_xx_GME*d_x u - str_yy_GME*d_y v + str_xy_GME*(d_y u + d_x v)
    ! This is the old formulation that includes energy diffusion
      FrictWork_GME(i,j,k) = H_to_RZ * ( &
              ((str_xx_GME(i,j,kk)*(u(I,j,k)-u(I-1,j,k))*IdxT(i,j))    &
             - (str_xx_GME(i,j,kk)*(v(i,J,k)-v(i,J-1,k))*IdyT(i,j)))   &
            + 0.25*(( (str_xy_GME(I,J,kk) *                              &
                       (((u(I,j+1,k)-u(I,j,k))*IdyBu(I,J))          &
                      + ((v(i+1,J,k)-v(i,J,k))*IdxBu(I,J))))        &
                    + (str_xy_GME(I-1,J-1,kk) *                          &
                       (((u(I-1,j,k)-u(I-1,j-1,k))*IdyBu(I-1,J-1))  &
                      + ((v(i,J-1,k)-v(i-1,J-1,k))*IdxBu(I-1,J-1)))) ) &
                  + ( (str_xy_GME(I-1,J,kk) *                            &
                       (((u(I-1,j+1,k)-u(I-1,j,k))*IdyBu(I-1,J))    &
                      + ((v(i,J,k)-v(i-1,J,k))*IdxBu(I-1,J))))      &
                    + (str_xy_GME(I,J-1,kk) *                            &
                       (((u(I,j,k)-u(I,j-1,k))*IdyBu(I,J-1))        &
                      + ((v(i+1,J-1,k)-v(i,J-1,k))*IdxBu(I,J-1)))) ) ) )
      enddo
    else ; do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      FrictWork_GME(i,j,k) = H_to_RZ * IareaT(i,j) * ( &
          ((str_xx_GME(i,j,kk)*dy2h(i,j) * ( &
                (uh(I,j,k)*dxCu(I,j)*IdyCu(I,j)*IareaCu(I,j)/(h_u(I,j,kk)+h_neglect)) &
              - (uh(I-1,j,k)*dxCu(I-1,j)*IdyCu(I-1,j)*IareaCu(I-1,j)/(h_u(I-1,j,kk)+h_neglect)) ) ) &
         - (str_xx_GME(i,j,kk)*dx2h(i,j) * ( &
                (vh(i,J,k)*dyCv(i,J)*IdxCv(i,J)*IareaCv(i,J)/(h_v(i,J,kk)+h_neglect)) &
              - (vh(i,J-1,k)*dyCv(i,J-1)*IdxCv(i,J-1)*IareaCv(i,J-1)/(h_v(i,J-1,kk)+h_neglect)) ) )) &
     + (0.25*(((str_xy_GME(I,J,kk)*(                                     &
                   (dx2q(I,J)*((uh(I,j+1,k)*IareaCu(I,j+1)/(h_u(I,j+1,kk)+h_neglect)) &
                                - (uh(I,j,k)*IareaCu(I,j)/(h_u(I,j,kk)+h_neglect))))            &
                 + (dy2q(I,J)*((vh(i+1,J,k)*IareaCv(i+1,J)/(h_v(i+1,J,kk)+h_neglect)) &
                                - (vh(i,J,k)*IareaCv(i,J)/(h_v(i,J,kk)+h_neglect)))) ))          &
              +(str_xy_GME(I-1,J-1,kk)*(                                 &
                   (dx2q(I-1,J-1)*((uh(I-1,j,k)*IareaCu(I-1,j)/(h_u(I-1,j,kk)+h_neglect)) &
                                    - (uh(I-1,j-1,k)*IareaCu(I-1,j-1)/(h_u(I-1,j-1,kk)+h_neglect))))    &
                 + (dy2q(I-1,J-1)*((vh(i,J-1,k)*IareaCv(i,J-1)/(h_v(i,J-1,kk)+h_neglect)) &
                                    - (vh(i-1,J-1,k)*IareaCv(i-1,J-1)/(h_v(i-1,J-1,kk)+h_neglect)))) )) ) &
             +((str_xy_GME(I-1,J,kk)*(                                   &
                   (dx2q(I-1,J)*((uh(I-1,j+1,k)*IareaCu(I-1,j+1)/(h_u(I-1,j+1,kk)+h_neglect)) &
                                  - (uh(I-1,j,k)*IareaCu(I-1,j)/(h_u(I-1,j,kk)+h_neglect))))      &
                 + (dy2q(I-1,J)*((vh(i,J,k)*IareaCv(i,J)/(h_v(i,J,kk)+h_neglect)) &
                                  - (vh(i-1,J,k)*IareaCv(i-1,J)/(h_v(i-1,J,kk)+h_neglect)))) ))        &
              +(str_xy_GME(I,J-1,kk)*(                                   &
                   (dx2q(I,J-1)*((uh(I,j,k)*IareaCu(I,j)/(h_u(I,j,kk)+h_neglect)) &
                                  - (uh(I,j-1,k)*IareaCu(I,j-1)/(h_u(I,j-1,kk)+h_neglect))))          &
                 + (dy2q(I,J-1)*((vh(i+1,J-1,k)*IareaCv(i+1,J-1)/(h_v(i+1,J-1,kk)+h_neglect)) &
                                  - (vh(i,J-1,k)*IareaCv(i,J-1)/(h_v(i,J-1,kk)+h_neglect)))) )) ) )) )
    enddo ; endif
    if (find_FrictWork .and. MEKE_a%mom_src_a%associated() .and. MEKE_a%GME_snk_a%associated()) then
      do kk=1,kmax ; do j=js,je ; do i=is,ie
        k = kstart+kk-1
        MEKE_GME_snk(i,j) = MEKE_GME_snk(i,j) + FrictWork_GME(i,j,k)
      enddo ; enddo ; enddo
    endif
  endif

  if (skeb_use_frict) then ; do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))
    k = kstart + kk - 1
    ! Note that the sign convention is FrictWork < 0 means energy dissipation.
    skeb_diss(i,j,k) = skeb_diss(i,j,k) - STOCH_a%skeb_frict_coef * &
                             FrictWork(i,j,k) / (H_to_RZ * (h(i,j,k) + h_neglect))
  enddo ; endif

  ! Make a similar calculation as for FrictWork above but accumulating into
  ! the vertically integrated MEKE source term, and adjusting for any
  ! energy loss seen as a reduction in the (biharmonic) frictional source term.
  if (find_FrictWork .and. MEKE_a%mom_src_a%associated()) then
    if (MEKE_a%backscatter_Ro_c /= 0.) then
      do kk=1,kmax ; do j=js,je ; do i=is,ie
        k = kstart + kk - 1
        FatH = 0.25*( (abs(CoriolisBu(I-1,J-1)) + abs(CoriolisBu(I,J))) + &
                      (abs(CoriolisBu(I-1,J)) + abs(CoriolisBu(I,J-1))) )
        Shear_mag_bc = sqrt(sh_xx(i,j,kk) * sh_xx(i,j,kk) + &
          0.25*(((sh_xy(I-1,J-1,kk)*sh_xy(I-1,J-1,kk)) + (sh_xy(I,J,kk)*sh_xy(I,J,kk))) + &
                ((sh_xy(I-1,J,kk)*sh_xy(I-1,J,kk)) + (sh_xy(I,J-1,kk)*sh_xy(I,J-1,kk)))))
        if ((CS%answer_date > 20190101) .and. (CS%answer_date < 20241201)) then
          FatH = (s_to_T*FatH)**MEKE_a%backscatter_Ro_pow ! f^n
          ! Note the hard-coded dimensional constant in the following line that can not
          ! be rescaled for dimensional consistency.
          Shear_mag_bc = (((s_to_T * Shear_mag_bc)**MEKE_a%backscatter_Ro_pow) + 1.e-30) &
                      * MEKE_a%backscatter_Ro_c ! c * D^n
          ! The Rossby number function is g(Ro) = 1/(1+c.Ro^n)
          ! RoScl = 1 - g(Ro)
          RoScl = Shear_mag_bc / (FatH + Shear_mag_bc) ! = 1 - f^n/(f^n+c*D^n)
        else
          if (FatH <= backscat_subround*Shear_mag_bc) then
            RoScl = 1.0
          else
            Sh_F_pow = MEKE_a%backscatter_Ro_c * (Shear_mag_bc / FatH)**MEKE_a%backscatter_Ro_pow
            RoScl = Sh_F_pow / (1.0 + Sh_F_pow) ! = 1 - f^n/(f^n+c*D^n)
          endif
        endif

        MEKE_mom_src(i,j) = MEKE_mom_src(i,j) + (FrictWork(i,j,k) - RoScl*FrictWork_bh(i,j,k))

        if (MEKE_a%mom_src_bh_a%associated()) &
          MEKE_mom_src_bh(i,j) = MEKE_mom_src_bh(i,j) &
              + (FrictWork_bh(i,j,k) - RoScl * FrictWork_bh(i,j,k))
      enddo ; enddo ; enddo
    else
      do kk=1,kmax ; do j=js,je ; do i=is,ie
        k = kstart + kk - 1
        MEKE_mom_src(i,j) = MEKE_mom_src(i,j) + FrictWork(i,j,k)
      enddo ; enddo ; enddo

      if (MEKE_a%mom_src_bh_a%associated()) then
        do kk=1,kmax ; do j=js,je ; do i=is,ie
          k = kstart + kk - 1
          MEKE_mom_src_bh(i,j) = MEKE_mom_src_bh(i,j) + FrictWork_bh(i,j,k)
        enddo ; enddo ; enddo
      endif
    endif ! MEKE_a%backscatter_Ro_c

  endif ! find_FrictWork and associated(mom_src)

end subroutine hor_visc_stress_div


!> Calculates the magnitude of the vorticity and divergence gradients, and the
!! Laplacian of vorticity, that are used within horizontal_viscosity by the Leith,
!! modified Leith, Leith+E, and QG Leith viscosity schemes. The caller must only
!! invoke this routine when CS%Leith_Kh, CS%Leith_Ah, or CS%use_Leithy is true.
subroutine hor_visc_Leith_grad_a(Gcore, CS, nkblock, kmax, is, ie, js, je, is_Kh, ie_Kh, js_Kh, &
                                 je_Kh, Ieq, Jeq, dudx_a, dvdy_a, vort_xy_a, vort_xy_smooth_a, &
                                 grad_vort_mag_h_2d_a, grad_div_mag_h_a, grad_vort_mag_q_2d_a, &
                                 grad_div_mag_q_a, Del2vort_q_a, vort_xy_dx_a, vort_xy_dx_smooth_a, &
                                 div_xx_dy_a, vort_xy_dy_a, vort_xy_dy_smooth_a, div_xx_dx_a)
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  integer,                       intent(in)    :: nkblock !< The k-block size used to size the following arrays [nondim]
  integer,                       intent(in)    :: kmax   !< The number of active k-layers in the current k-block
  integer,                       intent(in)    :: is     !< Start i-loop index for the h-point viscosities
  integer,                       intent(in)    :: ie     !< End i-loop index for the h-point viscosities
  integer,                       intent(in)    :: js     !< Start j-loop index for the h-point viscosities
  integer,                       intent(in)    :: je     !< End j-loop index for the h-point viscosities
  integer,                       intent(in)    :: is_Kh  !< Start i-loop index for the thickness point viscosities
  integer,                       intent(in)    :: ie_Kh  !< End i-loop index for the thickness point viscosities
  integer,                       intent(in)    :: js_Kh  !< Start j-loop index for the thickness point viscosities
  integer,                       intent(in)    :: je_Kh  !< End j-loop index for the thickness point viscosities
  integer,                       intent(in)    :: Ieq    !< The last i-index at q-points
  integer,                       intent(in)    :: Jeq    !< The last j-index at q-points
  type(RealArray_t), intent(in)    :: dudx_a !< x-derivative of the horizontal tension [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: dvdy_a !< y-derivative of the horizontal tension [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: vort_xy_a !< Vertical vorticity (dv/dx - du/dy) including
                                                !! metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: vort_xy_smooth_a !< Vertical vorticity including metric
                                                       !! terms, smoothed [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: grad_vort_mag_h_2d_a !< Magnitude of 2d vorticity gradient
                                                           !! at h-points [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: grad_div_mag_h_a !< Magnitude of divergence gradient at
                                                       !! h-points [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: grad_vort_mag_q_2d_a !< Magnitude of 2d vorticity gradient
                                                           !! at q-points [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: grad_div_mag_q_a !< Magnitude of divergence gradient at
                                                       !! q-points [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: Del2vort_q_a !< Laplacian of vorticity at q-points
                                                   !! [L-2 T-1 ~> m-2 s-1]
  type(RealArray_t), intent(inout) :: vort_xy_dx_a !< x-derivative of vertical vorticity (d/dx(dv/dx - du/dy))
                                                   !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: vort_xy_dx_smooth_a !< x-derivative of smoothed vertical vorticity
                                                          !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: div_xx_dy_a !< y-derivative of horizontal divergence (d/dy(du/dx +
                                                  !! dv/dy)) [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: vort_xy_dy_a !< y-derivative of vertical vorticity (d/dy(dv/dx - du/dy))
                                                   !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: vort_xy_dy_smooth_a !< y-derivative of smoothed vertical vorticity
                                                          !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: div_xx_dx_a !< x-derivative of horizontal divergence (d/dx(du/dx +
                                                  !! dv/dy)) [L-1 T-1 ~> m-1 s-1]
  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: div_xx_dx, vort_xy_dy_smooth, vort_xy_dy, &
    div_xx_dy, vort_xy_dx_smooth, vort_xy_dx, Del2vort_q, grad_div_mag_q, grad_vort_mag_q_2d, &
    grad_div_mag_h, grad_vort_mag_h_2d, vort_xy_smooth, vort_xy, dvdy, dudx
  real, dimension(:,:), contiguous, pointer :: dyBu, IdxBu, IdyCu, dxBu, IdyBu, IdxCv, IdyCv, &
    IdxCu, dF_dx, dF_dy
  real, dimension(dudx_a%lb(1):dudx_a%ub(1),dudx_a%lb(2):dudx_a%ub(2),nkblock) :: &
    div_xx               ! Estimate of horizontal divergence at h-points [T-1 ~> s-1]
  real :: DY_dxBu ! Ratio of meridional over zonal grid spacing at vertices [nondim]
  real :: DX_dyBu ! Ratio of zonal over meridional grid spacing at vertices [nondim]
  integer :: i, j, kk

  call div_xx_dx_a%view(div_xx_dx) ; call vort_xy_dy_smooth_a%view(vort_xy_dy_smooth)
  call vort_xy_dy_a%view(vort_xy_dy) ; call div_xx_dy_a%view(div_xx_dy)
  call vort_xy_dx_smooth_a%view(vort_xy_dx_smooth) ; call vort_xy_dx_a%view(vort_xy_dx)
  call Del2vort_q_a%view(Del2vort_q) ; call grad_div_mag_q_a%view(grad_div_mag_q)
  call grad_vort_mag_q_2d_a%view(grad_vort_mag_q_2d) ; call grad_div_mag_h_a%view(grad_div_mag_h)
  call grad_vort_mag_h_2d_a%view(grad_vort_mag_h_2d) ; call vort_xy_smooth_a%view(vort_xy_smooth)
  call vort_xy_a%view(vort_xy) ; call dvdy_a%view(dvdy) ; call dudx_a%view(dudx)
  call Gcore%dyBu%view(dyBu) ; call Gcore%IdxBu%view(IdxBu) ; call Gcore%IdyCu%view(IdyCu)
  call Gcore%dxBu%view(dxBu) ; call Gcore%IdyBu%view(IdyBu) ; call Gcore%IdxCv%view(IdxCv)
  call Gcore%IdyCv%view(IdyCv) ; call Gcore%IdxCu%view(IdxCu) ; call Gcore%dF_dx%view(dF_dx)
  call Gcore%dF_dy%view(dF_dy)

  ! Vorticity gradient
  do concurrent (kk=1:kmax, J=js-2:je_Kh, i=is_Kh-1:ie_Kh+1) DO_LOCALITY(local(DY_dxBu))
    DY_dxBu = dyBu(I,J) * IdxBu(I,J)
    vort_xy_dx(i,J,kk) = DY_dxBu * ((vort_xy(I,J,kk) * IdyCu(I,j)) - (vort_xy(I-1,J,kk) * IdyCu(I-1,j)))
  enddo

  do concurrent (kk=1:kmax, j=js_Kh-1:je_Kh+1, I=is-2:ie_Kh) DO_LOCALITY(local(DX_dyBu))
    DX_dyBu = dxBu(I,J) * IdyBu(I,J)
    vort_xy_dy(I,j,kk) = DX_dyBu * ((vort_xy(I,J,kk) * IdxCv(i,J)) - (vort_xy(I,J-1,kk) * IdxCv(i,J-1)))
  enddo

  if (CS%use_Leithy) then
    ! Gradient of smoothed vorticity
    do concurrent (kk=1:kmax, J=js_Kh-1:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(DY_dxBu))
      DY_dxBu = dyBu(I,J) * IdxBu(I,J)
      vort_xy_dx_smooth(i,J,kk) = DY_dxBu * &
                  ((vort_xy_smooth(I,J,kk) * IdyCu(I,j)) - (vort_xy_smooth(I-1,J,kk) * IdyCu(I-1,j)))
    enddo

    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, I=is_Kh-1:ie_Kh) DO_LOCALITY(local(DX_dyBu))
      DX_dyBu = dxBu(I,J) * IdyBu(I,J)
      vort_xy_dy_smooth(I,j,kk) = DX_dyBu * &
                  ((vort_xy_smooth(I,J,kk) * IdxCv(i,J)) - (vort_xy_smooth(I,J-1,kk) * IdxCv(i,J-1)))
    enddo
  endif ! If Leithy

  ! Laplacian of vorticity
  ! if (CS%Leith_Ah .or. CS%use_Leithy) then
  do concurrent (kk=1:kmax, J=js_Kh-1:je_Kh, I=is_Kh-1:ie_Kh) DO_LOCALITY(local(DY_dxBu, DX_dyBu))
    DY_dxBu = dyBu(I,J) * IdxBu(I,J)
    DX_dyBu = dxBu(I,J) * IdyBu(I,J)

    Del2vort_q(I,J,kk) = DY_dxBu * ((vort_xy_dx(i+1,J,kk) * IdyCv(i+1,J)) - (vort_xy_dx(i,J,kk) * IdyCv(i,J))) + &
                      DX_dyBu * ((vort_xy_dy(I,j+1,kk) * IdyCu(I,j+1)) - (vort_xy_dy(I,j,kk) * IdyCu(I,j)))
  enddo
  ! endif

  if (CS%modified_Leith) then

    ! Divergence
    do concurrent (kk=1:kmax, j=js_Kh-1:je_Kh+1, i=is_Kh-1:ie_Kh+1)
      div_xx(i,j,kk) = dudx(i,j,kk) + dvdy(i,j,kk)
    enddo

    ! Divergence gradient
    do concurrent (kk=1:kmax, j=js-1:je+1, I=is_Kh-1:ie_Kh)
      div_xx_dx(I,j,kk) = IdxCu(I,j)*(div_xx(i+1,j,kk) - div_xx(i,j,kk))
    enddo
    do concurrent (kk=1:kmax, J=js_Kh-1:je_Kh, i=is-1:ie+1)
      div_xx_dy(i,J,kk) = IdyCv(i,J)*(div_xx(i,j+1,kk) - div_xx(i,j,kk))
    enddo

    ! Magnitude of divergence gradient
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
      grad_div_mag_h(i,j,kk) = sqrt(((0.5*(div_xx_dx(I,j,kk) + div_xx_dx(I-1,j,kk)))**2) + &
                                 ((0.5*(div_xx_dy(i,J,kk) + div_xx_dy(i,J-1,kk)))**2))
    enddo
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
      grad_div_mag_q(I,J,kk) = sqrt(((0.5*(div_xx_dx(I,j,kk) + div_xx_dx(I,j+1,kk)))**2) + &
                                 ((0.5*(div_xx_dy(i,J,kk) + div_xx_dy(i+1,J,kk)))**2))
    enddo

  else

    do concurrent (kk=1:kmax, j=js-1:je+1, I=is_Kh-1:ie_Kh)
      div_xx_dx(I,j,kk) = 0.0
    enddo
    do concurrent (kk=1:kmax, J=js_Kh-1:je_Kh, i=is-1:ie+1)
      div_xx_dy(i,J,kk) = 0.0
    enddo
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
      grad_div_mag_h(i,j,kk) = 0.0
    enddo
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
      grad_div_mag_q(I,J,kk) = 0.0
    enddo

  endif ! CS%modified_Leith

  ! Add in beta for the Leith viscosity
  if (CS%use_beta_in_Leith) then
    do concurrent (kk=1:kmax, J=js-2:Jeq+1, i=is-1:ie+1)
      vort_xy_dx(i,J,kk) = vort_xy_dx(i,J,kk) + 0.5 * ( dF_dx(i,j) + dF_dx(i,j+1))
    enddo
    do concurrent (kk=1:kmax, j=js-1:je+1, I=is-2:Ieq+1)
      vort_xy_dy(I,j,kk) = vort_xy_dy(I,j,kk) + 0.5 * ( dF_dy(i,j) + dF_dy(i+1,j))
    enddo
  endif ! CS%use_beta_in_Leith

  if (CS%use_QG_Leith_visc) then

    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
      grad_vort_mag_h_2d(i,j,kk) = SQRT(((0.5*(vort_xy_dx(i,J,kk) + vort_xy_dx(i,J-1,kk)))**2) + &
                                     ((0.5*(vort_xy_dy(I,j,kk) + vort_xy_dy(I-1,j,kk)))**2) )
    enddo
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
      grad_vort_mag_q_2d(I,J,kk) = SQRT(((0.5*(vort_xy_dx(i,J,kk) + vort_xy_dx(i+1,J,kk)))**2) + &
                                     ((0.5*(vort_xy_dy(I,j,kk) + vort_xy_dy(I,j+1,kk)))**2) )
    enddo
  endif

end subroutine hor_visc_Leith_grad_a

!> Calculates the magnitudes of the vorticity gradients used within horizontal_viscosity by the
!! Leith viscosity schemes, after any QG Leith contribution has been added.
subroutine hor_visc_Leith_grad_b(CS, nkblock, kmax, is, js, is_Kh, ie_Kh, js_Kh, je_Kh, Ieq, Jeq, &
                                 grad_vort_mag_h_a, grad_vort_mag_q_a, vert_vort_mag_smooth_a, &
                                 vort_xy_dx_a, vort_xy_dx_smooth_a, vort_xy_dy_a, &
                                 vort_xy_dy_smooth_a)
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  integer,                       intent(in)    :: nkblock !< The k-block size used to size the following arrays [nondim]
  integer,                       intent(in)    :: kmax   !< The number of active k-layers in the current k-block
  integer,                       intent(in)    :: is     !< Start i-loop index for the h-point viscosities
  integer,                       intent(in)    :: js     !< Start j-loop index for the h-point viscosities
  integer,                       intent(in)    :: is_Kh  !< Start i-loop index for the thickness point viscosities
  integer,                       intent(in)    :: ie_Kh  !< End i-loop index for the thickness point viscosities
  integer,                       intent(in)    :: js_Kh  !< Start j-loop index for the thickness point viscosities
  integer,                       intent(in)    :: je_Kh  !< End j-loop index for the thickness point viscosities
  integer,                       intent(in)    :: Ieq    !< The last i-index at q-points
  integer,                       intent(in)    :: Jeq    !< The last j-index at q-points
  type(RealArray_t), intent(inout) :: grad_vort_mag_h_a !< Magnitude of vorticity gradient at
                                                        !! h-points [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: grad_vort_mag_q_a !< Magnitude of vorticity gradient at
                                                        !! q-points [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: vert_vort_mag_smooth_a !< Magnitude of gradient of smoothed
                                                             !! vertical vorticity (h or q) [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: vort_xy_dx_a !< x-derivative of vertical vorticity (d/dx(dv/dx - du/dy))
                                                   !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: vort_xy_dx_smooth_a !< x-derivative of smoothed vertical vorticity
                                                          !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: vort_xy_dy_a !< y-derivative of vertical vorticity (d/dy(dv/dx - du/dy))
                                                   !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: vort_xy_dy_smooth_a !< y-derivative of smoothed vertical vorticity
                                                          !! [L-1 T-1 ~> m-1 s-1]
  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: vort_xy_dy_smooth, vort_xy_dy, vort_xy_dx_smooth, &
    vort_xy_dx, vert_vort_mag_smooth, grad_vort_mag_q, grad_vort_mag_h
  integer :: i, j, kk

  call vort_xy_dy_smooth_a%view(vort_xy_dy_smooth) ; call vort_xy_dy_a%view(vort_xy_dy)
  call vort_xy_dx_smooth_a%view(vort_xy_dx_smooth) ; call vort_xy_dx_a%view(vort_xy_dx)
  call vert_vort_mag_smooth_a%view(vert_vort_mag_smooth)
  call grad_vort_mag_q_a%view(grad_vort_mag_q) ; call grad_vort_mag_h_a%view(grad_vort_mag_h)

  do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
    grad_vort_mag_h(i,j,kk) = SQRT(((0.5*(vort_xy_dx(i,J,kk) + vort_xy_dx(i,J-1,kk)))**2) + &
                                ((0.5*(vort_xy_dy(I,j,kk) + vort_xy_dy(I-1,j,kk)))**2) )
  enddo
  do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
    grad_vort_mag_q(I,J,kk) = SQRT(((0.5*(vort_xy_dx(i,J,kk) + vort_xy_dx(i+1,J,kk)))**2) + &
                                ((0.5*(vort_xy_dy(I,j,kk) + vort_xy_dy(I,j+1,kk)))**2) )
  enddo

  if (CS%use_Leithy) then
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
      vert_vort_mag_smooth(i,j,kk) = SQRT(((0.5*(vort_xy_dx_smooth(i,J,kk) + &
                                              vort_xy_dx_smooth(i,J-1,kk)))**2) + &
                                       ((0.5*(vort_xy_dy_smooth(I,j,kk) + &
                                              vort_xy_dy_smooth(I-1,j,kk)))**2) )
    enddo
  endif ! Leithy


end subroutine hor_visc_Leith_grad_b

!> Adds the QG Leith contribution to the vorticity gradients for each layer of one tile.
!! This calls into MOM_lateral_mixing_coeffs, so it always runs in Fortran.
subroutine hor_visc_QG_Leith_host(host, nkblock, kstart, kmax, h_a, dz_a, slope_x_a, slope_y_a, &
                                  div_xx_dx_a, div_xx_dy_a, vort_xy_dx_a, vort_xy_dy_a)
  type(c_ptr), value,   intent(in)    :: host  !< Opaque pointer to the hv_host_ctx that holds the Fortran
                                               !! objects used by the host helpers
  integer,                       intent(in)    :: nkblock !< The k-block size used to size the following arrays [nondim]
  integer,                       intent(in)    :: kstart !< The first absolute k-layer of the current k-block
  integer,                       intent(in)    :: kmax   !< The number of active k-layers in the current k-block
  type(RealArray_t), intent(in)    :: h_a !< Layer thicknesses [H ~> m or kg m-2].
  type(RealArray_t), intent(in)    :: dz_a !< Height change across layers [Z ~> m]
  type(RealArray_t), intent(inout) :: slope_x_a !< Isopycnal slope in i-direction [Z L-1 ~> nondim]
  type(RealArray_t), intent(inout) :: slope_y_a !< Isopycnal slope in j-direction [Z L-1 ~> nondim]
  type(RealArray_t), intent(in)    :: div_xx_dx_a !< x-derivative of horizontal divergence (d/dx(du/dx +
                                                  !! dv/dy)) [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(in)    :: div_xx_dy_a !< y-derivative of horizontal divergence (d/dy(du/dx +
                                                  !! dv/dy)) [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: vort_xy_dx_a !< x-derivative of vertical vorticity (d/dx(dv/dx - du/dy))
                                                   !! [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: vort_xy_dy_a !< y-derivative of vertical vorticity (d/dy(dv/dx - du/dy))
                                                   !! [L-1 T-1 ~> m-1 s-1]
  ! Local variables
  type(hv_host_ctx),       pointer :: ctx    ! The objects behind host
  type(ocean_grid_type),   pointer :: G      ! The ocean's grid structure.
  type(verticalGrid_type), pointer :: GV     ! The ocean's vertical grid structure.
  type(unit_scale_type),   pointer :: US     ! A dimensional unit scaling type
  type(VarMix_CS),         pointer :: VarMix ! Variable mixing control structure
  real, dimension(:,:,:), contiguous, pointer :: vort_xy_dy, vort_xy_dx, div_xx_dy, div_xx_dx, &
    slope_y, slope_x, dz, h
  integer :: k, kk

  call c_f_pointer(host, ctx)
  G => ctx%G ; GV => ctx%GV ; US => ctx%US ; VarMix => ctx%VarMix

  call vort_xy_dy_a%view(vort_xy_dy) ; call vort_xy_dx_a%view(vort_xy_dx)
  call div_xx_dy_a%view(div_xx_dy) ; call div_xx_dx_a%view(div_xx_dx)
  call slope_y_a%view(slope_y) ; call slope_x_a%view(slope_x) ; call dz_a%view(dz)
  call h_a%view(h)

  ! This accumulates terms, some of which are in VarMix.
  do kk=1,kmax
    k = kstart + kk - 1
    call calc_QG_Leith_viscosity(VarMix, G, GV, US, h, dz, k, div_xx_dx(:,:,kk), div_xx_dy(:,:,kk), &
                                 slope_x, slope_y, vort_xy_dx(:,:,kk), vort_xy_dy(:,:,kk))
  enddo

end subroutine hor_visc_QG_Leith_host

!> Adds the MEKE-based (EY24_EBT_BS) backscatter contribution to the diagonal term of
!! the stress tensor at h-points. The caller must only invoke this routine when
!! CS%EY24_EBT_BS is true.
subroutine hor_visc_backscatter_h(CS, MEKE_a, VarMix_a, use_kh_struct, nkblock, kstart, kmax, Isq, &
                                  Ieq, Jsq, Jeq, visc_limit_h_flag_a, sh_xx_a, str_xx_a, &
                                  BS_coeff_h_a)
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  type(hor_visc_MEKE_type), intent(in)    :: MEKE_a !< MEKE fields related to Mesoscale Eddy Kinetic Energy.
  type(hor_visc_VarMix_type), intent(in)    :: VarMix_a !< Variable mixing control structure
  logical,               intent(in)    :: use_kh_struct !< If true, shape the backscatter coefficient
                                                 !! with VarMix%BS_struct
  integer,               intent(in)    :: nkblock !< The k-block size used to size the following arrays [nondim]
  integer,               intent(in)    :: kstart  !< The first absolute k-layer of the current k-block
  integer,               intent(in)    :: kmax    !< The number of active k-layers in the current k-block
  integer,               intent(in)    :: Isq !< Start i-loop index for the stress tensor at h-points
  integer,               intent(in)    :: Ieq !< End i-loop index for the stress tensor at h-points
  integer,               intent(in)    :: Jsq !< Start j-loop index for the stress tensor at h-points
  integer,               intent(in)    :: Jeq !< End j-loop index for the stress tensor at h-points
  type(RealArray_t), intent(in)    :: visc_limit_h_flag_a !< determines whether backscatter is shut off [nondim]
  type(RealArray_t), intent(in)    :: sh_xx_a !< horizontal tension (du/dx - dv/dy) including
                                              !! metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: str_xx_a !< The diagonal term in the stress tensor
                                               !! [H L2 T-2 ~> m3 s-2 or kg s-2]
  type(RealArray_t), intent(inout) :: BS_coeff_h_a !< A diagnostic array of the backscatter
                                                   !! coefficient [L2 T-1 ~> m2 s-1]
  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: BS_coeff_h, str_xx, sh_xx, visc_limit_h_flag
  real, dimension(:,:), contiguous, pointer :: MEKE_Ku
  real, dimension(:,:,:), contiguous, pointer :: BS_struct
  real, dimension(sh_xx_a%lb(1):sh_xx_a%ub(1),sh_xx_a%lb(2):sh_xx_a%ub(2),nkblock) :: &
    Kh_BS, &     ! Laplacian antiviscosity [L2 T-1 ~> m2 s-1]
    str_xx_BS    ! The diagonal term in the stress tensor due to backscatter [H L2 T-2 ~> m3 s-2 or kg s-2]
  integer :: i, j, k, kk

  call BS_coeff_h_a%view(BS_coeff_h) ; call str_xx_a%view(str_xx) ; call sh_xx_a%view(sh_xx)
  call visc_limit_h_flag_a%view(visc_limit_h_flag)

  nullify(MEKE_Ku, BS_struct)
  if (MEKE_a%Ku_a%associated()) call MEKE_a%Ku_a%view(MEKE_Ku)
  if (VarMix_a%BS_struct_a%associated()) call VarMix_a%BS_struct_a%view(BS_struct)

  do kk=1,kmax ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    k = kstart + kk - 1
    if (visc_limit_h_flag(i,j,k) > 0) then
      Kh_BS(i,j,kk) = 0.
    else
      if (use_kh_struct) then
        Kh_BS(i,j,kk) = MEKE_Ku(i,j) * BS_struct(i,j,k)
      else
        Kh_BS(i,j,kk) = MEKE_Ku(i,j)
      endif
    endif
  enddo ; enddo ; enddo

  do concurrent (kk=1:kmax, j=Jsq:Jeq+1, i=Isq:Ieq+1)
    str_xx_BS(i,j,kk) = -Kh_BS(i,j,kk) * sh_xx(i,j,kk)
  enddo

  if (CS%do_BS_coeff_h) then
    do concurrent (kk=1:kmax, j=Jsq:Jeq+1, i=Isq:Ieq+1) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      BS_coeff_h(i,j,k) = Kh_BS(i,j,kk)
    enddo
  endif

  do concurrent (kk=1:kmax, j=Jsq:Jeq+1, i=Isq:Ieq+1)
    str_xx(i,j,kk) = str_xx(i,j,kk) + str_xx_BS(i,j,kk)
  enddo

end subroutine hor_visc_backscatter_h

!> Adds the MEKE-based (EY24_EBT_BS) backscatter contribution to the cross term of
!! the stress tensor at q-points. The caller must only invoke this routine when
!! CS%EY24_EBT_BS is true.
subroutine hor_visc_backscatter_q(CS, MEKE_a, VarMix_a, use_kh_struct, nkblock, kstart, kmax, is, &
                                  js, Ieq, Jeq, visc_limit_q_flag_a, sh_xy_a, str_xy_a, &
                                  BS_coeff_q_a)
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  type(hor_visc_MEKE_type), intent(in)    :: MEKE_a !< MEKE fields related to Mesoscale Eddy Kinetic Energy.
  type(hor_visc_VarMix_type), intent(in)    :: VarMix_a !< Variable mixing control structure
  logical,                 intent(in)    :: use_kh_struct !< If true, shape the backscatter coefficient
                                                 !! with VarMix%BS_struct
  integer,                 intent(in)    :: nkblock !< The k-block size used to size the following arrays [nondim]
  integer,                 intent(in)    :: kstart  !< The first absolute k-layer of the current k-block
  integer,                 intent(in)    :: kmax    !< The number of active k-layers in the current k-block
  integer,                 intent(in)    :: is  !< Start i-loop index for the q-point stress tensor
  integer,                 intent(in)    :: js  !< Start j-loop index for the q-point stress tensor
  integer,                 intent(in)    :: Ieq !< End i-loop index for the q-point stress tensor
  integer,                 intent(in)    :: Jeq !< End j-loop index for the q-point stress tensor
  type(RealArray_t), intent(in)    :: visc_limit_q_flag_a !< determines whether backscatter is shut off [nondim]
  type(RealArray_t), intent(in)    :: sh_xy_a !< horizontal shearing strain (du/dy + dv/dx) including
                                              !! metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: str_xy_a !< The cross term in the stress tensor
                                               !! [H L2 T-2 ~> m3 s-2 or kg s-2]
  type(RealArray_t), intent(inout) :: BS_coeff_q_a !< A diagnostic array of the backscatter
                                                   !! coefficient [L2 T-1 ~> m2 s-1]
  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: BS_coeff_q, str_xy, sh_xy, visc_limit_q_flag
  real, dimension(:,:), contiguous, pointer :: MEKE_Ku
  real, dimension(:,:,:), contiguous, pointer :: BS_struct
  real, dimension(sh_xy_a%lb(1):sh_xy_a%ub(1),sh_xy_a%lb(2):sh_xy_a%ub(2),nkblock) :: &
    Kh_BS, &     ! Laplacian antiviscosity [L2 T-1 ~> m2 s-1]
    str_xy_BS    ! The cross term in the stress tensor due to backscatter [H L2 T-2 ~> m3 s-2 or kg s-2]
  integer :: i, j, k, kk

  call BS_coeff_q_a%view(BS_coeff_q) ; call str_xy_a%view(str_xy) ; call sh_xy_a%view(sh_xy)
  call visc_limit_q_flag_a%view(visc_limit_q_flag)

  nullify(MEKE_Ku, BS_struct)
  if (MEKE_a%Ku_a%associated()) call MEKE_a%Ku_a%view(MEKE_Ku)
  if (VarMix_a%BS_struct_a%associated()) call VarMix_a%BS_struct_a%view(BS_struct)

  do kk=1,kmax ; do J=js-1,Jeq ; do I=is-1,Ieq
    k = kstart + kk - 1
    if (visc_limit_q_flag(I,J,k) > 0) then
      Kh_BS(I,J,kk) = 0.
    else
      if (use_kh_struct) then
        Kh_BS(I,J,kk) = 0.25*( ((MEKE_Ku(i,j)*BS_struct(i,j,k)) + &
                             (MEKE_Ku(i+1,j+1)*BS_struct(i+1,j+1,k))) + &
                            ((MEKE_Ku(i+1,j)*BS_struct(i+1,j,k)) + &
                             (MEKE_Ku(i,j+1)*BS_struct(i,j+1,k))) )
      else
        Kh_BS(I,J,kk) = 0.25*( (MEKE_Ku(i,j) + MEKE_Ku(i+1,j+1)) + &
                            (MEKE_Ku(i+1,j) + MEKE_Ku(i,j+1)) )
      endif
    endif
  enddo ; enddo ; enddo

  do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
    str_xy_BS(I,J,kk) = -Kh_BS(I,J,kk) * (sh_xy(I,J,kk))
  enddo

  if (CS%do_BS_coeff_q) then
    do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      BS_coeff_q(I,J,k) = Kh_BS(I,J,kk)
    enddo
  endif

  do concurrent (kk=1:kmax, J=js-1:Jeq, I=is-1:Ieq)
    str_xy(I,J,kk) = str_xy(I,J,kk) + str_xy_BS(I,J,kk)
  enddo

end subroutine hor_visc_backscatter_q

!> Computes the Leith+E antisymmetric-viscosity ratio m_leithy, before any smoothing. The
!! caller must only invoke this routine when CS%use_Leithy is true.
subroutine hor_visc_Leithy_m(Gcore, CS, metrics, nkblock, kmax, is_Kh, ie_Kh, js_Kh, je_Kh, &
                             inv_PI6, Del2vort_q_a, vert_vort_mag_a, vort_xy_smooth_a, m_leithy_a)
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  type(hv_metrics_CS), intent(in)    :: metrics       !< Persistent copies of the horizontal
                                                     !! viscosity metric and coefficient arrays
  integer,                 intent(in)    :: nkblock !< The k-block size used to size the following arrays [nondim]
  integer,                 intent(in)    :: kmax    !< The number of active k-layers in the current k-block
  integer,                 intent(in)    :: is_Kh !< Start i-loop index for the thickness point viscosities
  integer,                 intent(in)    :: ie_Kh !< End i-loop index for the thickness point viscosities
  integer,                 intent(in)    :: js_Kh !< Start j-loop index for the thickness point viscosities
  integer,                 intent(in)    :: je_Kh !< End j-loop index for the thickness point viscosities
  real,                    intent(in)    :: inv_PI6 !< The inverse of pi to the sixth power [nondim]
  type(RealArray_t), intent(in)    :: Del2vort_q_a !< Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
  type(RealArray_t), intent(in)    :: vert_vort_mag_a !< Magnitude of the vertical vorticity
                                                      !! gradient (h or q) [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(in)    :: vort_xy_smooth_a !< Vertical vorticity including metric
                                                       !! terms, smoothed [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: m_leithy_a !< Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: m_leithy, vort_xy_smooth, vert_vort_mag, Del2vort_q
  real, dimension(:,:), contiguous, pointer :: mask2dBu
  real, dimension(:,:), contiguous, pointer :: Ah_bg_xx, m_const_leithy, m_leithy_max, Biharm6_const_xx
  real :: Del2vort_h ! Laplacian of vorticity at h-points [L-2 T-1 ~> m-2 s-1]
  real :: AhLth      ! 2D Leith biharmonic viscosity [L4 T-1 ~> m4 s-1]
  integer :: i, j, kk

  call m_leithy_a%view(m_leithy) ; call vort_xy_smooth_a%view(vort_xy_smooth)
  call vert_vort_mag_a%view(vert_vort_mag) ; call Del2vort_q_a%view(Del2vort_q)
  call Gcore%mask2dBu%view(mask2dBu)

  nullify(Ah_bg_xx, m_const_leithy, m_leithy_max, Biharm6_const_xx)
  if (metrics%Ah_bg_xx%associated()) call metrics%Ah_bg_xx%view(Ah_bg_xx)
  if (metrics%m_const_leithy%associated()) call metrics%m_const_leithy%view(m_const_leithy)
  if (metrics%m_leithy_max%associated()) call metrics%m_leithy_max%view(m_leithy_max)
  if (metrics%Biharm6_const_xx%associated()) call metrics%Biharm6_const_xx%view(Biharm6_const_xx)

  ! Get m_leithy
  if (CS%smooth_Ah) then
    do kk=1,kmax
      m_leithy(:,:,kk) = 0.0 ! This is here to initialize domain edge halo values.
    enddo
  endif
  do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(Del2vort_h, AhLth))
    Del2vort_h = 0.25 * ((Del2vort_q(I,J,kk) + Del2vort_q(I-1,J-1,kk)) + &
                         (Del2vort_q(I-1,J,kk) + Del2vort_q(I,J-1,kk)))
    AhLth  = Biharm6_const_xx(i,j) * inv_PI6 * abs(Del2vort_h)
    if (AhLth <= Ah_bg_xx(i,j)) then
      m_leithy(i,j,kk) = 0.0
    else
      if ((m_const_leithy(i,j)*vert_vort_mag(i,j,kk)) < abs(vort_xy_smooth(i,j,kk))) then
        m_leithy(i,j,kk) = CS%c_K * (vert_vort_mag(i,j,kk) / vort_xy_smooth(i,j,kk))**2
      else
        m_leithy(i,j,kk) = m_leithy_max(i,j)
      endif
      m_leithy(i,j,kk) = mask2dBu(i,j) * m_leithy(i,j,kk)
    endif
  enddo

end subroutine hor_visc_Leithy_m

!> Computes the Leith+E antisymmetric-viscosity ratio m_leithy, and updates the
!! biharmonic viscosity Ah (and Ah_h) to include the Leith+E contribution, optionally
!! smoothing m_leithy and Ah in the process. The caller must only invoke this
!! routine when CS%use_Leithy is true.
subroutine hor_visc_Leithy_Ah(CS, metrics, nkblock, kmax, is_Kh, ie_Kh, js_Kh, je_Kh, inv_PI6, &
                              Del2vort_q_a, vert_vort_mag_smooth_a, Ah_a, m_leithy_a, Ah_sq_a)
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  type(hv_metrics_CS), intent(in)    :: metrics       !< Persistent copies of the horizontal
                                                     !! viscosity metric and coefficient arrays
  integer,                 intent(in)    :: nkblock !< The k-block size used to size the following arrays [nondim]
  integer,                 intent(in)    :: kmax    !< The number of active k-layers in the current k-block
  integer,                 intent(in)    :: is_Kh !< Start i-loop index for the thickness point viscosities
  integer,                 intent(in)    :: ie_Kh !< End i-loop index for the thickness point viscosities
  integer,                 intent(in)    :: js_Kh !< Start j-loop index for the thickness point viscosities
  integer,                 intent(in)    :: je_Kh !< End j-loop index for the thickness point viscosities
  real,                    intent(in)    :: inv_PI6 !< The inverse of pi to the sixth power [nondim]
  type(RealArray_t), intent(in)    :: Del2vort_q_a !< Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
  type(RealArray_t), intent(in)    :: vert_vort_mag_smooth_a !< Magnitude of gradient of smoothed
                                                             !! vertical vorticity (h or q) [L-1 T-1 ~> m-1 s-1]
  type(RealArray_t), intent(inout) :: Ah_a !< biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
  type(RealArray_t), intent(in)    :: m_leithy_a !< Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
  type(RealArray_t), intent(inout) :: Ah_sq_a !< The square of the biharmonic viscosity [L8 T-2 ~> m8 s-2]
  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: Ah_sq, m_leithy, Ah, vert_vort_mag_smooth, &
    Del2vort_q
  real, dimension(:,:), contiguous, pointer :: Ah_bg_xx, Biharm6_const_xx
  real :: Del2vort_h ! Laplacian of vorticity at h-points [L-2 T-1 ~> m-2 s-1]
  real :: AhLthy     ! 2D Leith+E biharmonic viscosity [L4 T-1 ~> m4 s-1]
  integer :: i, j, kk

  call Ah_sq_a%view(Ah_sq) ; call m_leithy_a%view(m_leithy) ; call Ah_a%view(Ah)
  call vert_vort_mag_smooth_a%view(vert_vort_mag_smooth) ; call Del2vort_q_a%view(Del2vort_q)

  nullify(Ah_bg_xx, Biharm6_const_xx)
  if (metrics%Ah_bg_xx%associated()) call metrics%Ah_bg_xx%view(Ah_bg_xx)
  if (metrics%Biharm6_const_xx%associated()) call metrics%Biharm6_const_xx%view(Biharm6_const_xx)

  ! Get Ah
  do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(Del2vort_h, AhLthy))
    Del2vort_h = 0.25 * ((Del2vort_q(I,J,kk) + Del2vort_q(I-1,J-1,kk)) + &
                         (Del2vort_q(I-1,J,kk) + Del2vort_q(I,J-1,kk)))
    AhLthy = Biharm6_const_xx(i,j) * inv_PI6 * &
            sqrt(max(0.,Del2vort_h**2 - m_leithy(i,j,kk)*vert_vort_mag_smooth(i,j,kk)**2))
    Ah(i,j,kk) = max(Ah_bg_xx(i,j), AhLthy)
  enddo
  if (CS%smooth_Ah) then
    ! Smooth Ah before applying upper bound.  Square Ah, then smooth, then take its square root.
    do kk=1,kmax
      Ah_sq(:,:,kk) = 0.0 ! This is here to initialize domain edge halo values.
    enddo
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
      Ah_sq(i,j,kk) = Ah(i,j,kk)**2
    enddo
  endif

end subroutine hor_visc_Leithy_Ah

!> Sets the Leith+E biharmonic viscosity Ah and Ah_h, from the smoothed square of Ah if
!! CS%smooth_Ah is true. The caller must only invoke this routine when CS%use_Leithy is true.
subroutine hor_visc_Leithy_Ah_h(CS, metrics, nkblock, kstart, kmax, is_Kh, ie_Kh, js_Kh, je_Kh, &
                                Ah_a, Ah_h_a, Ah_sq_a)
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  type(hv_metrics_CS), intent(in)    :: metrics       !< Persistent copies of the horizontal
                                                     !! viscosity metric and coefficient arrays
  integer,                 intent(in)    :: nkblock !< The k-block size used to size the following arrays [nondim]
  integer,                 intent(in)    :: kstart  !< The first absolute k-layer of the current k-block
  integer,                 intent(in)    :: kmax    !< The number of active k-layers in the current k-block
  integer,                 intent(in)    :: is_Kh !< Start i-loop index for the thickness point viscosities
  integer,                 intent(in)    :: ie_Kh !< End i-loop index for the thickness point viscosities
  integer,                 intent(in)    :: js_Kh !< Start j-loop index for the thickness point viscosities
  integer,                 intent(in)    :: je_Kh !< End j-loop index for the thickness point viscosities
  type(RealArray_t), intent(inout) :: Ah_a !< biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
  type(RealArray_t), intent(inout) :: Ah_h_a !< biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
  type(RealArray_t), intent(in)    :: Ah_sq_a !< The square of the biharmonic viscosity [L8 T-2 ~> m8 s-2]
  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: Ah_sq, Ah_h, Ah
  real, dimension(:,:), contiguous, pointer :: Ah_bg_xx
  integer :: i, j, k, kk

  call Ah_sq_a%view(Ah_sq) ; call Ah_h_a%view(Ah_h) ; call Ah_a%view(Ah)

  nullify(Ah_bg_xx)
  if (metrics%Ah_bg_xx%associated()) call metrics%Ah_bg_xx%view(Ah_bg_xx)

  if (CS%smooth_Ah) then
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      Ah_h(i,j,k) = max(Ah_bg_xx(i,j), sqrt(max(0., Ah_sq(i,j,kk))))
      Ah(i,j,kk)    = Ah_h(i,j,k)
    enddo
  else
    do concurrent (kk=1:kmax, j=js_Kh:je_Kh, i=is_Kh:ie_Kh) DO_LOCALITY(local(k))
      k = kstart + kk - 1
      Ah_h(i,j,k) = Ah(i,j,kk)
    enddo
  endif

end subroutine hor_visc_Leithy_Ah_h

!> Smooths m_leithy for each layer of one tile.
!! This includes halo updates, so it always runs in Fortran.
subroutine hor_visc_Leithy_m_smooth_host(host, Gcore, bxH0, nkblock, kmax, m_leithy_a)
  type(c_ptr), value,   intent(in)    :: host  !< Opaque pointer to the hv_host_ctx that holds the Fortran
                                               !! objects used by the host helpers
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),          intent(in)    :: bxH0  !< The h-point box of the whole column,
                                               !! [isc:iec, jsc:jec, 1:nz]
  integer,                 intent(in)    :: nkblock !< The k-block size used to size the following arrays [nondim]
  integer,                 intent(in)    :: kmax    !< The number of active k-layers in the current k-block
  type(RealArray_t), intent(inout) :: m_leithy_a !< Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
  ! Local variables
  type(hv_host_ctx),       pointer :: ctx    ! The objects behind host
  type(ocean_grid_type),   pointer :: G      ! The ocean's grid structure.
  type(RealArray_t) :: layer_a ! One layer of the field being smoothed [arbitrary]
  real, dimension(:,:,:), contiguous, pointer :: m_leithy
  integer :: kk

  call c_f_pointer(host, ctx)
  G => ctx%G

  call m_leithy_a%view(m_leithy)

  do kk=1,kmax
    ! Smooth m_leithy.  A single call smoothes twice.
    call pass_var(m_leithy(:,:,kk), G%Domain, halo=2)
    call layer_a%alloc(lb=m_leithy_a%lb(1:2), ub=m_leithy_a%ub(1:2), &
                       source=m_leithy(:,:,kk))
    call smooth_x9_h(Gcore, bxH0, layer_a, zero_land=.true.)
    call layer_a%copy2F(m_leithy(:,:,kk)) ; call layer_a%free()
    call pass_var(m_leithy(:,:,kk), G%Domain)
  enddo

end subroutine hor_visc_Leithy_m_smooth_host

!> Smooths the square of Ah for each layer of one tile.
!! This includes halo updates, so it always runs in Fortran.
subroutine hor_visc_Leithy_Ah_sq_smooth_host(host, Gcore, bxH0, nkblock, kmax, Ah_sq_a)
  type(c_ptr), value,   intent(in)    :: host  !< Opaque pointer to the hv_host_ctx that holds the Fortran
                                               !! objects used by the host helpers
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),          intent(in)    :: bxH0  !< The h-point box of the whole column,
                                               !! [isc:iec, jsc:jec, 1:nz]
  integer,                 intent(in)    :: nkblock !< The k-block size used to size the following arrays [nondim]
  integer,                 intent(in)    :: kmax    !< The number of active k-layers in the current k-block
  type(RealArray_t), intent(inout) :: Ah_sq_a !< The square of the biharmonic viscosity [L8 T-2 ~> m8 s-2]
  ! Local variables
  type(hv_host_ctx),       pointer :: ctx    ! The objects behind host
  type(ocean_grid_type),   pointer :: G      ! The ocean's grid structure.
  type(RealArray_t) :: layer_a ! One layer of the field being smoothed [arbitrary]
  real, dimension(:,:,:), contiguous, pointer :: Ah_sq
  integer :: kk

  call c_f_pointer(host, ctx)
  G => ctx%G

  call Ah_sq_a%view(Ah_sq)

  do kk=1,kmax
    call pass_var(Ah_sq(:,:,kk), G%Domain, halo=2)
    ! A single call smoothes twice.
    call layer_a%alloc(lb=Ah_sq_a%lb(1:2), ub=Ah_sq_a%ub(1:2), &
                       source=Ah_sq(:,:,kk))
    call smooth_x9_h(Gcore, bxH0, layer_a, zero_land=.false.)
    call layer_a%copy2F(Ah_sq(:,:,kk)) ; call layer_a%free()
    call pass_var(Ah_sq(:,:,kk), G%Domain)
  enddo

end subroutine hor_visc_Leithy_Ah_sq_smooth_host

!> Returns the effective k-block size for horizontal_viscosity calls.
integer function hor_visc_nkblock(CS, GV)
  type(hor_visc_CS),       intent(in) :: CS !< Horizontal viscosity control structure
  type(verticalGrid_type), intent(in) :: GV !< The ocean's vertical grid structure.
  hor_visc_nkblock = merge(GV%ke, CS%nkblock, CS%nkblock==0)
end function hor_visc_nkblock

!> Calculates the barotropic tension and shearing strain fields and the GME
!! efficiency and isopycnal height diffusivity fields that are used within
!! horizontal_viscosity when GME (CS%use_GME) is active. The caller must
!! only invoke this routine when CS%use_GME is true.
subroutine hor_visc_GME_setup(G, GV, US, CS, Gcore, bxH0, bxQ0, h, BT, TD, nkblock, dudx_bt, &
                              dvdy_bt, dvdx_bt, dudy_bt, sh_xx_bt, sh_xy_bt, GME_effic_h, &
                              GME_effic_q, KH_u_GME, KH_v_GME, GME_coeff_h, GME_coeff_q, &
                              str_xx_GME, str_xy_GME)
  type(ocean_grid_type),         intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),         intent(in)    :: US     !< A dimensional unit scaling type
  type(hor_visc_CS),             intent(inout) :: CS     !< Horizontal viscosity control structure
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),          intent(in)    :: bxH0  !< The h-point box of the whole column,
                                               !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),          intent(in)    :: bxQ0  !< The B-grid-index box of the whole column,
                                               !! [IscB:IecB, JscB:JecB, 1:nz]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                  intent(inout) :: h      !< Layer thicknesses [H ~> m or kg m-2].
  type(barotropic_CS),           optional, intent(in) :: BT !< Barotropic control structure
  type(thickness_diffuse_CS),    optional, intent(in) :: TD !< Thickness diffusion control structure
  integer,                       intent(in)    :: nkblock !< The k-block size used to size the following arrays [nondim]
  real, dimension(SZI_(G),SZJB_(G)), &
                                  intent(out)   :: dudx_bt !< x-component in the barotropic
                                                       !! horizontal tension [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJB_(G)), &
                                  intent(out)   :: dvdy_bt !< y-component in the barotropic
                                                       !! horizontal tension [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G)), &
                                  intent(out)   :: dvdx_bt !< x-component in the barotropic
                                                       !! shearing strain [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G)), &
                                  intent(out)   :: dudy_bt !< y-component in the barotropic
                                                       !! shearing strain [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G)), &
                                  intent(out)   :: sh_xx_bt !< Barotropic horizontal tension
                                                       !! (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G)), &
                                  intent(out)   :: sh_xy_bt !< Barotropic horizontal shearing strain
                                                       !! (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G)), &
                                  intent(out)   :: GME_effic_h !< The filtered efficiency of the
                                                       !! GME terms at h points [nondim]
  real, dimension(SZIB_(G),SZJB_(G)), &
                                  intent(out)   :: GME_effic_q !< The filtered efficiency of the
                                                       !! GME terms at q points [nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                  intent(out)   :: KH_u_GME !< Isopycnal height diffusivities in
                                                       !! u-columns [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                  intent(out)   :: KH_v_GME !< Isopycnal height diffusivities in
                                                       !! v-columns [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                  intent(out)   :: GME_coeff_h !< GME coefficient at h-points
                                                       !! [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)), &
                                  intent(out)   :: GME_coeff_q !< GME coeff. at q-points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),nkblock), &
                                  intent(out)   :: str_xx_GME !< Smoothed diagonal term in the
                                                       !! stress tensor from GME [L2 T-2 ~> m2 s-2]
  real, dimension(SZIB_(G),SZJB_(G),nkblock), &
                                  intent(out)   :: str_xy_GME !< Smoothed cross term in the
                                                       !! stress tensor from GME [L2 T-2 ~> m2 s-2]
  ! Local variables
  type(RealArray_t) :: h_a, ubtav_a, vbtav_a, dudx_bt_a, dvdy_bt_a, dvdx_bt_a, dudy_bt_a, &
    sh_xx_bt_a, sh_xy_bt_a, GME_effic_h_a, GME_effic_q_a
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    ubtav         ! zonal barotropic velocity averaged over a baroclinic time-step [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    vbtav         ! meridional barotropic velocity averaged over a baroclinic time-step [L T-1 ~> m s-1]

  ! Initialize diagnostic arrays with zeros
  GME_coeff_h(:,:,:) = 0.0
  GME_coeff_q(:,:,:) = 0.0
  str_xx_GME(:,:,:) = 0.0
  str_xy_GME(:,:,:) = 0.0

  ! Get barotropic velocities and their gradients
  call barotropic_get_tav(BT, ubtav, vbtav, G, US)

  call pass_vector(ubtav, vbtav, G%Domain)
  call pass_var(h, G%domain, halo=2)

  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), &
                 source=h)
  call ubtav_a%alloc(lb=LBOUND(ubtav), ub=UBOUND(ubtav), &
                     source=ubtav)
  call vbtav_a%alloc(lb=LBOUND(vbtav), ub=UBOUND(vbtav), &
                     source=vbtav)
  call dudx_bt_a%alloc(lb=LBOUND(dudx_bt), ub=UBOUND(dudx_bt), &
                       source=dudx_bt)
  call dvdy_bt_a%alloc(lb=LBOUND(dvdy_bt), ub=UBOUND(dvdy_bt), &
                       source=dvdy_bt)
  call dvdx_bt_a%alloc(lb=LBOUND(dvdx_bt), ub=UBOUND(dvdx_bt), &
                       source=dvdx_bt)
  call dudy_bt_a%alloc(lb=LBOUND(dudy_bt), ub=UBOUND(dudy_bt), &
                       source=dudy_bt)
  call sh_xx_bt_a%alloc(lb=LBOUND(sh_xx_bt), ub=UBOUND(sh_xx_bt), &
                        source=sh_xx_bt)
  call sh_xy_bt_a%alloc(lb=LBOUND(sh_xy_bt), ub=UBOUND(sh_xy_bt), &
                        source=sh_xy_bt)
  call GME_effic_h_a%alloc(lb=LBOUND(GME_effic_h), ub=UBOUND(GME_effic_h), &
                           source=GME_effic_h)
  call GME_effic_q_a%alloc(lb=LBOUND(GME_effic_q), ub=UBOUND(GME_effic_q), &
                           source=GME_effic_q)
  call hor_visc_GME_bt(Gcore, bxH0, bxQ0, CS%hv_opts_CS, CS%hv_metrics_CS, h_a, ubtav_a, vbtav_a, &
                       dudx_bt_a, dvdy_bt_a, dvdx_bt_a, dudy_bt_a, sh_xx_bt_a, sh_xy_bt_a, &
                       GME_effic_h_a, GME_effic_q_a)
  call h_a%copy2F(h) ; call dudx_bt_a%copy2F(dudx_bt) ; call dvdy_bt_a%copy2F(dvdy_bt)
  call dvdx_bt_a%copy2F(dvdx_bt) ; call dudy_bt_a%copy2F(dudy_bt)
  call sh_xx_bt_a%copy2F(sh_xx_bt) ; call sh_xy_bt_a%copy2F(sh_xy_bt)
  call GME_effic_h_a%copy2F(GME_effic_h) ; call GME_effic_q_a%copy2F(GME_effic_q)
  call h_a%free() ; call ubtav_a%free() ; call vbtav_a%free() ; call dudx_bt_a%free()
  call dvdy_bt_a%free() ; call dvdx_bt_a%free() ; call dudy_bt_a%free() ; call sh_xx_bt_a%free()
  call sh_xy_bt_a%free() ; call GME_effic_h_a%free() ; call GME_effic_q_a%free()

  call thickness_diffuse_get_KH(TD, KH_u_GME, KH_v_GME, G, GV)

  call pass_vector(KH_u_GME, KH_v_GME, G%domain, To_All+Scalar_Pair)

  if (CS%hv_opts_CS%debug) &
    call uvchksum("GME KH[u,v]_GME", KH_u_GME, KH_v_GME, G%HI, haloshift=2, unscale=US%L_to_m**2*US%s_to_T)

end subroutine hor_visc_GME_setup

!> Calculates the barotropic tension and shearing strain fields and the GME efficiency
!! fields from the barotropic velocities, for use within horizontal_viscosity when GME
!! (CS%use_GME) is active.
subroutine hor_visc_GME_bt(Gcore, bxH0, bxQ0, CS, metrics, h_a, ubtav_a, vbtav_a, dudx_bt_a, &
                           dvdy_bt_a, dvdx_bt_a, dudy_bt_a, sh_xx_bt_a, sh_xy_bt_a, GME_effic_h_a, &
                           GME_effic_q_a)
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),          intent(in)    :: bxH0  !< The h-point box of the whole column,
                                               !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),          intent(in)    :: bxQ0  !< The B-grid-index box of the whole column,
                                               !! [IscB:IecB, JscB:JecB, 1:nz]
  type(hv_opts_CS),    intent(in)    :: CS    !< Options used by the horizontal viscosity kernels
  type(hv_metrics_CS), intent(in)    :: metrics       !< Persistent copies of the horizontal
                                                     !! viscosity metric and coefficient arrays
  type(RealArray_t), intent(inout) :: h_a !< Layer thicknesses [H ~> m or kg m-2].
  type(RealArray_t), intent(in)    :: ubtav_a !< zonal barotropic velocity averaged over a baroclinic
                                              !! time-step [L T-1 ~> m s-1]
  type(RealArray_t), intent(in)    :: vbtav_a !< meridional barotropic velocity averaged over a baroclinic
                                              !! time-step [L T-1 ~> m s-1]
  type(RealArray_t), intent(inout) :: dudx_bt_a !< x-component in the barotropic
                                                !! horizontal tension [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: dvdy_bt_a !< y-component in the barotropic
                                                !! horizontal tension [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: dvdx_bt_a !< x-component in the barotropic
                                                !! shearing strain [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: dudy_bt_a !< y-component in the barotropic
                                                !! shearing strain [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xx_bt_a !< Barotropic horizontal tension
                                                 !! (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: sh_xy_bt_a !< Barotropic horizontal shearing strain
                                                 !! (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
  type(RealArray_t), intent(inout) :: GME_effic_h_a !< The filtered efficiency of the
                                                    !! GME terms at h points [nondim]
  type(RealArray_t), intent(inout) :: GME_effic_q_a !< The filtered efficiency of the
                                                    !! GME terms at q points [nondim]
  ! Local variables
  real, dimension(:,:), contiguous, pointer :: GME_effic_q, GME_effic_h, sh_xy_bt, sh_xx_bt, &
    dudy_bt, dvdx_bt, dvdy_bt, dudx_bt, vbtav, ubtav
  real, dimension(:,:), contiguous, pointer :: IdyCu, IdxCv, IdyCv, IdxCu, mask2dBu, mask2dCu, &
    mask2dCv, mask2dT
  real, dimension(:,:,:), contiguous, pointer :: h
  real, dimension(:,:), contiguous, pointer :: dx_dyT, dy_dxT, dx_dyBu, dy_dxBu
  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2)) :: &
    htot          ! The total thickness of all layers [H ~> m or kg m-2]
  real :: grad_vel_mag_bt_h ! Magnitude of the barotropic velocity gradient tensor squared at h-points [T-2 ~> s-2]
  real :: boundary_mask_h ! A mask that zeroes out cells with at least one land edge [nondim]
  real :: grad_vel_mag_bt_q ! Magnitude of the barotropic velocity gradient tensor squared at q-points [T-2 ~> s-2]
  real :: boundary_mask_q ! A mask that zeroes out cells with at least one land edge [nondim]
  real :: h_arith_q  ! The arithmetic mean total thickness at q points [H ~> m or kg m-2]
  real :: I_GME_h0   ! The inverse of GME tapering scale [H-1 ~> m-1 or m2 kg-1]
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz

  call GME_effic_q_a%view(GME_effic_q) ; call GME_effic_h_a%view(GME_effic_h)
  call sh_xy_bt_a%view(sh_xy_bt) ; call sh_xx_bt_a%view(sh_xx_bt) ; call dudy_bt_a%view(dudy_bt)
  call dvdx_bt_a%view(dvdx_bt) ; call dvdy_bt_a%view(dvdy_bt) ; call dudx_bt_a%view(dudx_bt)
  call vbtav_a%view(vbtav) ; call ubtav_a%view(ubtav) ; call h_a%view(h)
  call Gcore%IdyCu%view(IdyCu) ; call Gcore%IdxCv%view(IdxCv) ; call Gcore%IdyCv%view(IdyCv)
  call Gcore%IdxCu%view(IdxCu) ; call Gcore%mask2dBu%view(mask2dBu)
  call Gcore%mask2dCu%view(mask2dCu) ; call Gcore%mask2dCv%view(mask2dCv)
  call Gcore%mask2dT%view(mask2dT)

  nullify(dx_dyT, dy_dxT, dx_dyBu, dy_dxBu)
  if (metrics%dx_dyT%associated()) call metrics%dx_dyT%view(dx_dyT)
  if (metrics%dy_dxT%associated()) call metrics%dy_dxT%view(dy_dxT)
  if (metrics%dx_dyBu%associated()) call metrics%dx_dyBu%view(dx_dyBu)
  if (metrics%dy_dxBu%associated()) call metrics%dy_dxBu%view(dy_dxBu)

  is  = bxH0%idxS(1)  ; ie  = bxH0%idxE(1)  ; js  = bxH0%idxS(2)  ; je  = bxH0%idxE(2) ; nz = bxH0%idxE(3)
  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2)

  ! Calculate the barotropic horizontal tension
  do j=js-2,je+2 ; do i=is-2,ie+2
    dudx_bt(i,j) = DY_dxT(i,j)*((IdyCu(I,j) * ubtav(I,j)) - &
                                   (IdyCu(I-1,j) * ubtav(I-1,j)))
    dvdy_bt(i,j) = DX_dyT(i,j)*((IdxCv(i,J) * vbtav(i,J)) - &
                                   (IdxCv(i,J-1) * vbtav(i,J-1)))
  enddo ; enddo
  do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
    sh_xx_bt(i,j) = dudx_bt(i,j) - dvdy_bt(i,j)
  enddo ; enddo

  ! Components for the barotropic shearing strain
  do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
    dvdx_bt(I,J) = DY_dxBu(I,J)*((vbtav(i+1,J)*IdyCv(i+1,J)) &
                                  - (vbtav(i,J)*IdyCv(i,J)))
    dudy_bt(I,J) = DX_dyBu(I,J)*((ubtav(I,j+1)*IdxCu(I,j+1)) &
                                  - (ubtav(I,j)*IdxCu(I,j)))
  enddo ; enddo

  if (CS%no_slip) then
    do J=js-2,je+1 ; do I=is-2,ie+1
      sh_xy_bt(I,J) = (2.0-mask2dBu(I,J)) * ( dvdx_bt(I,J) + dudy_bt(I,J) )
    enddo ; enddo
  else
    do J=js-2,je+1 ; do I=is-2,ie+1
      sh_xy_bt(I,J) = mask2dBu(I,J) * ( dvdx_bt(I,J) + dudy_bt(I,J) )
    enddo ; enddo
  endif

  do j=js-2,je+2 ; do i=is-2,ie+2
    htot(i,j) = 0.0
  enddo ; enddo
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    htot(i,j) = htot(i,j) + h(i,j,k)
  enddo ; enddo ; enddo

  I_GME_h0 = 1.0 / CS%GME_h0
  do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
    boundary_mask_h = (mask2dCu(I,j) * mask2dCu(I-1,j)) * (mask2dCv(i,J) * mask2dCv(i,J-1))
    grad_vel_mag_bt_h = mask2dT(I,J) * boundary_mask_h * (dudx_bt(i,j)**2 + dvdy_bt(i,j)**2 + &
          (0.25*((dvdx_bt(I,J)+dvdx_bt(I-1,J-1)) + (dvdx_bt(I,J-1)+dvdx_bt(I-1,J))))**2 + &
          (0.25*((dudy_bt(I,J)+dudy_bt(I-1,J-1)) + (dudy_bt(I,J-1)+dudy_bt(I-1,J))))**2)
    ! Probably the following test could be simplified to
    ! if (boundary_mask_h * G%mask2dT(I,J) > 0.0) then
    if (grad_vel_mag_bt_h > 0.0) then
      GME_effic_h(i,j) = CS%GME_efficiency * mask2dT(I,J) * (MIN(htot(i,j) * I_GME_h0, 1.0)**2)
    else
      GME_effic_h(i,j) = 0.0
    endif
  enddo ; enddo

  do J=js-2,je+1 ; do I=is-2,ie+1
    boundary_mask_q = (mask2dCv(i,J) * mask2dCv(i+1,J)) * (mask2dCu(I,j) * mask2dCu(I,j+1))
    grad_vel_mag_bt_q = mask2dBu(I,J) * boundary_mask_q * (dvdx_bt(I,J)**2 + dudy_bt(I,J)**2 + &
          (0.25*((dudx_bt(i,j)+dudx_bt(i+1,j+1)) + (dudx_bt(i,j+1)+dudx_bt(i+1,j))))**2 + &
          (0.25*((dvdy_bt(i,j)+dvdy_bt(i+1,j+1)) + (dvdy_bt(i,j+1)+dvdy_bt(i+1,j))))**2)
    ! Probably the following test could be simplified to
    ! if (boundary_mask_q * G%mask2dBu(I,J) > 0.0) then
    if (grad_vel_mag_bt_q > 0.0) then
      h_arith_q = 0.25 * ((htot(i,j) + htot(i+1,j+1)) + (htot(i+1,j) + htot(i,j+1)))
      GME_effic_q(I,J) = CS%GME_efficiency * mask2dBu(I,J) * (MIN(h_arith_q * I_GME_h0, 1.0)**2)
    else
      GME_effic_q(I,J) = 0.0
    endif
  enddo ; enddo

end subroutine hor_visc_GME_bt

!> Allocates space for and calculates static variables used by horizontal_viscosity.
!! hor_visc_init calculates and stores the values of a number of metric functions that
!! are used in horizontal_viscosity.
subroutine hor_visc_init(Time, G, GV, US, param_file, diag, CS, ADp)
  type(time_type),         intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< Structure to regulate diagnostic output.
  type(hor_visc_CS),       intent(inout) :: CS   !< Horizontal viscosity control structure
  type(accel_diag_ptrs), intent(in), optional :: ADp !< Acceleration diagnostics

  ! u0v is the Laplacian sensitivities to the v velocities at u points, with u0u, v0u, and v0v defined analogously.
  real, dimension(SZIB_(G),SZJ_(G)) :: u0u, u0v ! Laplacian sensitivities at u points [L-2 ~> m-2]
  real, dimension(SZI_(G),SZJB_(G)) :: v0u, v0v ! Laplacian sensitivities at v points [L-2 ~> m-2]
  real :: grid_sp_h2       ! Harmonic mean of the squares of the grid [L2 ~> m2]
  real :: grid_sp_h3       ! Harmonic mean of the squares of the grid^(3/2) [L3 ~> m3]
  real :: grid_sp_q2       ! spacings at h and q points [L2 ~> m2]
  real :: grid_sp_q3       ! spacings at h and q points^(3/2) [L3 ~> m3]
  real :: min_grid_sp_h2   ! Minimum value of grid_sp_h2 [L2 ~> m2]
  real :: min_grid_sp_h4   ! Minimum value of grid_sp_h2**2 [L4 ~> m4]
  real :: Kh_Limit         ! A coefficient [T-1 ~> s-1] used, along with the
                           ! grid spacing, to limit Laplacian viscosity.
  real :: fmax             ! maximum absolute value of f at the four
                           ! vorticity points around a thickness point [T-1 ~> s-1]
  real :: BoundCorConst    ! A constant used when using viscosity to bound the Coriolis accelerations
                           ! [T2 L-2 ~> s2 m-2]
  real :: Ah_Limit         ! coefficient [T-1 ~> s-1] used, along with the
                           ! grid spacing, to limit biharmonic viscosity
  real :: Kh               ! Lapacian horizontal viscosity [L2 T-1 ~> m2 s-1]
  real :: Ah               ! biharmonic horizontal viscosity [L4 T-1 ~> m4 s-1]
  real :: Kh_vel_scale     ! this speed [L T-1 ~> m s-1] times grid spacing gives Laplacian viscosity
  real :: Ah_vel_scale     ! this speed [L T-1 ~> m s-1] times grid spacing cubed gives biharmonic viscosity
  real :: Ah_time_scale    ! damping time-scale for biharmonic visc [T ~> s]
  real :: Smag_Lap_const   ! nondimensional Laplacian Smagorinsky constant [nondim]
  real :: Smag_bi_const    ! nondimensional biharmonic Smagorinsky constant [nondim]
  real :: Leith_Lap_const  ! nondimensional Laplacian Leith constant [nondim]
  real :: Leith_bi_const   ! nondimensional biharmonic Leith constant [nondim]
  real :: dt               ! The dynamics time step [T ~> s]
  real :: Idt              ! The inverse of dt [T-1 ~> s-1]
  real :: denom            ! work variable; the denominator of a fraction [L-2 ~> m-2] or [L-4 ~> m-4]
  real :: maxvel           ! largest permitted velocity components [L T-1 ~> m s-1]
  real :: bound_Cor_vel    ! grid-scale velocity variations at which value
                           ! the quadratically varying biharmonic viscosity
                           ! balances Coriolis acceleration [L T-1 ~> m s-1]
  real :: Kh_sin_lat       ! Amplitude of latitudinally dependent viscosity [L2 T-1 ~> m2 s-1]
  real :: Kh_pwr_of_sine   ! Power used to raise sin(lat) when using Kh_sin_lat [nondim]
  logical :: bound_Cor_def ! parameter setting of BOUND_CORIOLIS
  logical :: split         ! If true, use the split time stepping scheme.
                           ! If false and USE_GME = True, issue a FATAL error.
  logical :: use_MEKE      ! If true, the MEKE parameterization is in use.
  logical :: enable_bugs   ! If true, the defaults for recently added bug-fix flags are set to
                           ! recreate the bugs, or if false bugs are only used if actively selected.
  real    :: backscatter_Ro_c ! Coefficient in Rossby number function for backscatter [nondim]
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags
  character(len=200) :: inputdir, filename ! Input file names and paths
  character(len=80) ::  Kh_var ! Input variable names
  real    :: deg2rad       ! Converts degrees to radians [radians degree-1]
  real    :: slat_fn       ! sin(lat)**Kh_pwr_of_sine [nondim]
  real    :: aniso_grid_dir(2) ! Vector (n1,n2) for anisotropic direction [nondim]
  integer :: aniso_mode    ! Selects the mode for setting the anisotropic direction
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: i, j
#ifdef __NVCOMPILER_OPENMP_GPU
  integer, parameter :: default_nkblock = 0
#else
  integer, parameter :: default_nkblock = 1
#endif
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_hor_visc"  ! module name
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = GV%ke
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! init control structure
  call ZB2020_init(Time, G, GV, US, param_file, diag, CS%ZB2020, CS%hv_opts_CS%use_ZB2020)

  CS%initialized = .true.

  CS%diag => diag
  ! Read parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "HORVISC_NKBLOCK", CS%nkblock, &
                 "The k-direction block size used in horizontal viscosity calculations. "//&
                 "The default 0 setting dynamically uses the full vertical column.", &
                 default=default_nkblock, layoutParam=.true.)
  if (CS%nkblock < 0) call MOM_error(FATAL, "HORVISC_NKBLOCK must be >= 0.")

  call get_param(param_file, mdl, "USE_CIRCULATION_IN_HORVISC", CS%hv_opts_CS%use_circulation, &
                 "Use circulation theorem to compute vorticity in horvisc module (for ZB20 or Leith)", &
                 default=.False.)

  ! All parameters are read in all cases to enable parameter spelling checks.
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)

  ! Determine whether HOR_VISC_ANSWER_DATE is used, and avoid logging it if it is not used.
  call get_param(param_file, mdl, "USE_MEKE", use_MEKE, &
                 default=.false., do_not_log=.true.)
  backscatter_Ro_c = 0.0
  if (use_MEKE) call get_param(param_file, mdl, "MEKE_BACKSCAT_RO_C", backscatter_Ro_c, &
                 "The coefficient in the Rossby number function for scaling the biharmonic "//&
                 "frictional energy source. Setting to non-zero enables the Rossby number function.", &
                 units="nondim", default=0.0, do_not_log=.true.)

  call get_param(param_file, mdl, "HOR_VISC_ANSWER_DATE", CS%hv_opts_CS%answer_date, &
                 "The vintage of the order of arithmetic and expressions in the horizontal "//&
                 "viscosity calculations.  Values between 20190102 and 20241201 recover the "//&
                 "answers from the end of 2018, while higher values use updated and more robust "//&
                 "forms of the same expressions.", &
                 default=default_answer_date, do_not_log=(.not.GV%Boussinesq).or.(backscatter_Ro_c==0.0))
  if (.not.GV%Boussinesq) CS%hv_opts_CS%answer_date = max(CS%hv_opts_CS%answer_date, 20241201)

  call get_param(param_file, mdl, "DEBUG", CS%hv_opts_CS%debug, default=.false.)
  call get_param(param_file, mdl, "USE_CONT_THICKNESS", CS%hv_opts_CS%use_cont_thick, &
                 "If true, use thickness at velocity points from continuity solver. This option "//&
                 "currently only works with split mode.", default=.false.)
  call get_param(param_file, mdl, "USE_CONT_THICKNESS_BUG", CS%use_cont_thick_bug, &
                   "If true, retain an answer-changing halo update bug when "//&
                   "USE_CONT_THICKNESS=True.  This is not recommended.", &
                   default=.false., do_not_log=.not.CS%hv_opts_CS%use_cont_thick)

  call get_param(param_file, mdl, "LAPLACIAN", CS%hv_opts_CS%Laplacian, &
                 "If true, use a Laplacian horizontal viscosity.", &
                 default=.false.)

  call get_param(param_file, mdl, "KH", Kh,                      &
                 "The background Laplacian horizontal viscosity.", &
                 units="m2 s-1", default=0.0, scale=US%m_to_L**2*US%T_to_s, &
                 do_not_log=.not.CS%hv_opts_CS%Laplacian)
  call get_param(param_file, mdl, "KH_BG_MIN", CS%hv_opts_CS%Kh_bg_min, &
                 "The minimum value allowed for Laplacian horizontal viscosity, KH.", &
                 units="m2 s-1", default=0.0, scale=US%m_to_L**2*US%T_to_s, &
                 do_not_log=.not.CS%hv_opts_CS%Laplacian)
  call get_param(param_file, mdl, "KH_VEL_SCALE", Kh_vel_scale, &
                 "The velocity scale which is multiplied by the grid "//&
                 "spacing to calculate the Laplacian viscosity. "//&
                 "The final viscosity is the largest of this scaled "//&
                 "viscosity, the Smagorinsky and Leith viscosities, and KH.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T, &
                 do_not_log=.not.CS%hv_opts_CS%Laplacian)
  call get_param(param_file, mdl, "KH_SIN_LAT", Kh_sin_lat, &
                 "The amplitude of a latitudinally-dependent background "//&
                 "viscosity of the form KH_SIN_LAT*(SIN(LAT)**KH_PWR_OF_SINE).", &
                 units="m2 s-1", default=0.0, scale=US%m_to_L**2*US%T_to_s, &
                 do_not_log=.not.CS%hv_opts_CS%Laplacian)
  call get_param(param_file, mdl, "KH_PWR_OF_SINE", Kh_pwr_of_sine, &
                 "The power used to raise SIN(LAT) when using a latitudinally "//&
                 "dependent background viscosity.", &
                 units="nondim", default=4.0, &
                 do_not_log=.not.(CS%hv_opts_CS%Laplacian .and. (Kh_sin_lat>0.)) )
  call get_param(param_file, mdl, "SMAGORINSKY_KH", CS%hv_opts_CS%Smagorinsky_Kh, &
                 "If true, use a Smagorinsky nonlinear eddy viscosity.", &
                 default=.false., do_not_log=.not.CS%hv_opts_CS%Laplacian)
  if (.not.CS%hv_opts_CS%Laplacian) CS%hv_opts_CS%Smagorinsky_Kh = .false.
  call get_param(param_file, mdl, "SMAG_LAP_CONST", Smag_Lap_const, &
                 "The nondimensional Laplacian Smagorinsky constant, "//&
                 "often 0.15.", units="nondim", default=0.0, &
                 fail_if_missing=CS%hv_opts_CS%Smagorinsky_Kh, do_not_log=.not.CS%hv_opts_CS%Smagorinsky_Kh)
  call get_param(param_file, mdl, "LEITH_KH", CS%hv_opts_CS%Leith_Kh, &
                 "If true, use a Leith nonlinear eddy viscosity.", &
                 default=.false., do_not_log=.not.CS%hv_opts_CS%Laplacian)
  if (.not.CS%hv_opts_CS%Laplacian) CS%hv_opts_CS%Leith_Kh = .false.
  call get_param(param_file, mdl, "LEITH_LAP_CONST", Leith_Lap_const, &
                 "The nondimensional Laplacian Leith constant, "//&
                 "often set to 1.0", units="nondim", default=0.0, &
                  fail_if_missing=CS%hv_opts_CS%Leith_Kh, do_not_log=.not.CS%hv_opts_CS%Leith_Kh)
  call get_param(param_file, mdl, "RES_SCALE_MEKE_VISC", CS%hv_opts_CS%res_scale_MEKE, &
                 "If true, the viscosity contribution from MEKE is scaled by "//&
                 "the resolution function.", default=.false., &
                 do_not_log=.not.(CS%hv_opts_CS%Laplacian.and.use_MEKE))
  if (.not.(CS%hv_opts_CS%Laplacian.and.use_MEKE)) CS%hv_opts_CS%res_scale_MEKE = .false.

  call get_param(param_file, mdl, "BOUND_KH", CS%hv_opts_CS%bound_Kh, &
                 "If true, the Laplacian coefficient is locally limited "//&
                 "to be stable.", default=.true., do_not_log=.not.CS%hv_opts_CS%Laplacian)
  call get_param(param_file, mdl, "EY24_EBT_BS", CS%hv_opts_CS%EY24_EBT_BS, &
                 "If true, use the backscatter scheme (EBT mode with kill switch) "//&
                 "developed by Yankovsky et al. (2024). ", &
                 default=.false., do_not_log=.not.CS%hv_opts_CS%Laplacian)
  if (.not.CS%hv_opts_CS%Laplacian) CS%hv_opts_CS%bound_Kh = .false.
  if (.not.(CS%hv_opts_CS%Laplacian.and.use_MEKE)) CS%hv_opts_CS%EY24_EBT_BS = .false.
  call get_param(param_file, mdl, "ANISOTROPIC_VISCOSITY", CS%hv_opts_CS%anisotropic, &
                 "If true, allow anistropic viscosity in the Laplacian "//&
                 "horizontal viscosity.", default=.false., &
                 do_not_log=.not.CS%hv_opts_CS%Laplacian)
  if (.not.CS%hv_opts_CS%Laplacian) &
    CS%hv_opts_CS%anisotropic = .false. ! This replicates the prior code, but is it intended?
  call get_param(param_file, mdl, "ADD_LES_VISCOSITY", CS%hv_opts_CS%add_LES_viscosity, &
                 "If true, adds the viscosity from Smagorinsky and Leith to the "//&
                 "background viscosity instead of taking the maximum.", default=.false., &
                 do_not_log=.not.CS%hv_opts_CS%Laplacian)

  call get_param(param_file, mdl, "KH_ANISO", CS%hv_opts_CS%Kh_aniso, &
                 "The background Laplacian anisotropic horizontal viscosity.", &
                 units="m2 s-1", default=0.0, scale=US%m_to_L**2*US%T_to_s, &
                 do_not_log=.not.CS%hv_opts_CS%anisotropic)
  call get_param(param_file, mdl, "ANISOTROPIC_MODE", aniso_mode, &
                 "Selects the mode for setting the direction of anisotropy.\n"//&
                 "\t 0 - Points along the grid i-direction.\n"//&
                 "\t 1 - Points towards East.\n"//&
                 "\t 2 - Points along the flow direction, U/|U|.", &
                 default=0, do_not_log=.not.CS%hv_opts_CS%anisotropic)
  if (aniso_mode == 0) then
    call get_param(param_file, mdl, "ANISO_GRID_DIR", aniso_grid_dir, &
                 "The vector pointing in the direction of anisotropy for horizontal viscosity. "//&
                 "n1,n2 are the i,j components relative to the grid.", &
                 units="nondim", fail_if_missing=CS%hv_opts_CS%anisotropic, do_not_log=.not.CS%hv_opts_CS%anisotropic)
  elseif (aniso_mode == 1) then
    call get_param(param_file, mdl, "ANISO_GRID_DIR", aniso_grid_dir, &
                 "The vector pointing in the direction of anisotropy for horizontal viscosity. "//&
                 "n1,n2 are the i,j components relative to the spherical coordinates.", &
                 units="nondim", fail_if_missing=CS%hv_opts_CS%anisotropic, do_not_log=.not.CS%hv_opts_CS%anisotropic)
  else
    call get_param(param_file, mdl, "ANISO_GRID_DIR", aniso_grid_dir, &
                 "The vector pointing in the direction of anisotropy for horizontal viscosity.", &
                 units="nondim", fail_if_missing=.false., do_not_log=.true.)
  endif

  call get_param(param_file, mdl, "BIHARMONIC", CS%hv_opts_CS%biharmonic, &
                 "If true, use a biharmonic horizontal viscosity. "//&
                 "BIHARMONIC may be used with LAPLACIAN.", &
                 default=.true.)
  call get_param(param_file, mdl, "AH", Ah, &
                 "The background biharmonic horizontal viscosity.", &
                 units="m4 s-1", default=0.0, scale=US%m_to_L**4*US%T_to_s, &
                 do_not_log=.not.CS%hv_opts_CS%biharmonic)
  call get_param(param_file, mdl, "AH_VEL_SCALE", Ah_vel_scale, &
                 "The velocity scale which is multiplied by the cube of "//&
                 "the grid spacing to calculate the biharmonic viscosity. "//&
                 "The final viscosity is the largest of this scaled "//&
                 "viscosity, the Smagorinsky and Leith viscosities, and AH.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T, do_not_log=.not.CS%hv_opts_CS%biharmonic)
  call get_param(param_file, mdl, "AH_TIME_SCALE", Ah_time_scale, &
                 "A time scale whose inverse is multiplied by the fourth "//&
                 "power of the grid spacing to calculate biharmonic viscosity. "//&
                 "The final viscosity is the largest of all viscosity "//&
                 "formulations in use. 0.0 means that it's not used.", &
                 units="s", default=0.0, scale=US%s_to_T, do_not_log=.not.CS%hv_opts_CS%biharmonic)
  call get_param(param_file, mdl, "SMAGORINSKY_AH", CS%hv_opts_CS%Smagorinsky_Ah, &
                 "If true, use a biharmonic Smagorinsky nonlinear eddy "//&
                 "viscosity.", default=.false., do_not_log=.not.CS%hv_opts_CS%biharmonic)
  if (.not.CS%hv_opts_CS%biharmonic) CS%hv_opts_CS%Smagorinsky_Ah = .false.
  call get_param(param_file, mdl, "LEITH_AH", CS%hv_opts_CS%Leith_Ah, &
                 "If true, use a biharmonic Leith nonlinear eddy "//&
                 "viscosity.", default=.false., do_not_log=.not.CS%hv_opts_CS%biharmonic)
  if (.not.CS%hv_opts_CS%biharmonic) CS%hv_opts_CS%Leith_Ah = .false.
  call get_param(param_file, mdl, "USE_LEITHY", CS%hv_opts_CS%use_Leithy, &
                 "If true, use a biharmonic Leith nonlinear eddy "//&
                 "viscosity together with a harmonic backscatter.", &
                 default=.false.)
  call get_param(param_file, mdl, "BOUND_AH", CS%hv_opts_CS%bound_Ah, &
                 "If true, the biharmonic coefficient is locally limited "//&
                 "to be stable.", default=.true., do_not_log=.not.CS%hv_opts_CS%biharmonic)
  if (.not.CS%hv_opts_CS%biharmonic) CS%hv_opts_CS%bound_Ah = .false.
  call get_param(param_file, mdl, "RE_AH", CS%hv_opts_CS%Re_Ah, &
                 "If nonzero, the biharmonic coefficient is scaled "//&
                 "so that the biharmonic Reynolds number is equal to this.", &
                 units="nondim", default=0.0, do_not_log=.not.CS%hv_opts_CS%biharmonic)

  call get_param(param_file, mdl, "BACKSCATTER_UNDERBOUND", CS%hv_opts_CS%backscatter_underbound, &
                 "If true, the bounds on the biharmonic viscosity are allowed to "//&
                 "increase where the Laplacian viscosity is negative (due to backscatter "//&
                 "parameterizations) beyond the largest timestep-dependent stable values of "//&
                 "biharmonic viscosity when no Laplacian viscosity is applied.  The default "//&
                 "is true for historical reasons, but this option probably should not be used "//&
                 "because it can contribute to numerical instabilities.", &
                 default=.false., do_not_log=.not.((CS%hv_opts_CS%bound_Kh).and.(CS%hv_opts_CS%bound_Ah)))

  call get_param(param_file, mdl, "SMAG_BI_CONST",Smag_bi_const, &
                 "The nondimensional biharmonic Smagorinsky constant, "//&
                 "typically 0.015 - 0.06.", units="nondim", default=0.0, &
                 fail_if_missing=CS%hv_opts_CS%Smagorinsky_Ah, do_not_log=.not.CS%hv_opts_CS%Smagorinsky_Ah)

  call get_param(param_file, mdl, "USE_BETA_IN_LEITH", CS%hv_opts_CS%use_beta_in_Leith, &
                 "If true, include the beta term in the Leith nonlinear eddy viscosity.", &
                 default=CS%hv_opts_CS%Leith_Kh, do_not_log=.not.(CS%hv_opts_CS%Leith_Kh .or. CS%hv_opts_CS%Leith_Ah) )
  call get_param(param_file, mdl, "MODIFIED_LEITH", CS%hv_opts_CS%modified_Leith, &
                 "If true, add a term to Leith viscosity which is "//&
                 "proportional to the gradient of divergence.", &
                 default=.false., do_not_log=.not.(CS%hv_opts_CS%Leith_Kh .or. CS%hv_opts_CS%Leith_Ah) )
  call get_param(param_file, mdl, "USE_QG_LEITH_VISC", CS%hv_opts_CS%use_QG_Leith_visc, &
                 "If true, use QG Leith nonlinear eddy viscosity.", &
                 default=.false., do_not_log=.not.(CS%hv_opts_CS%Leith_Kh .or. CS%hv_opts_CS%Leith_Ah) )
!  if (CS%use_QG_Leith_visc) then
!    call MOM_error(FATAL, "USE_QG_LEITH_VISC=True activates code that is a work-in-progress and "//&
!          "should not be used until a number of bugs are fixed.  Specifically it does not "//&
!          "reproduce across PE count or layout, and may use arrays that have not been properly "//&
!          "set or allocated.  See github.com/mom-ocean/MOM6/issues/1590 for a discussion.")
!  endif
  if (CS%hv_opts_CS%use_QG_Leith_visc .and. .not. (CS%hv_opts_CS%Leith_Kh .or. CS%hv_opts_CS%Leith_Ah) ) then
    call MOM_error(FATAL, "MOM_hor_visc.F90, hor_visc_init: "//&
                 "LEITH_KH or LEITH_AH must be True when USE_QG_LEITH_VISC=True.")
  endif

  call get_param(param_file, mdl, "BOUND_CORIOLIS", bound_Cor_def, default=.false.)
  call get_param(param_file, mdl, "BOUND_CORIOLIS_BIHARM", CS%hv_opts_CS%bound_Coriolis, &
                 "If true use a viscosity that increases with the square "//&
                 "of the velocity shears, so that the resulting viscous "//&
                 "drag is of comparable magnitude to the Coriolis terms "//&
                 "when the velocity differences between adjacent grid "//&
                 "points is 0.5*BOUND_CORIOLIS_VEL.  The default is the "//&
                 "value of BOUND_CORIOLIS (or false).", default=bound_Cor_def, &
                 do_not_log=.not.CS%hv_opts_CS%Smagorinsky_Ah)
  if (.not.CS%hv_opts_CS%Smagorinsky_Ah) CS%hv_opts_CS%bound_Coriolis = .false.
  call get_param(param_file, mdl, "MAXVEL", maxvel, &
                 units="m s-1", default=3.0e8, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "BOUND_CORIOLIS_VEL", bound_Cor_vel, &
                 "The velocity scale at which BOUND_CORIOLIS_BIHARM causes "//&
                 "the biharmonic drag to have comparable magnitude to the "//&
                 "Coriolis acceleration.  The default is set by MAXVEL.", &
                 units="m s-1", default=maxvel*US%L_T_to_m_s, scale=US%m_s_to_L_T, &
                 do_not_log=.not.(CS%hv_opts_CS%Smagorinsky_Ah .and. CS%hv_opts_CS%bound_Coriolis))
  call get_param(param_file, mdl, "LEITH_BI_CONST", Leith_bi_const, &
                 "The nondimensional biharmonic Leith constant, "//&
                 "typical values are thus far undetermined.", units="nondim", default=0.0, &
                 fail_if_missing=(CS%hv_opts_CS%Leith_Ah .or. CS%hv_opts_CS%use_Leithy), &
                 do_not_log=.not.(CS%hv_opts_CS%Leith_Ah .or. CS%hv_opts_CS%use_Leithy))
  call get_param(param_file, mdl, "USE_LAND_MASK_FOR_HVISC", CS%hv_opts_CS%use_land_mask, &
                 "If true, use the land mask for the computation of thicknesses "//&
                 "at velocity locations. This eliminates the dependence on arbitrary "//&
                 "values over land or outside of the domain.", default=.true.)
  call get_param(param_file, mdl, "HORVISC_BOUND_COEF", CS%bound_coef, &
                 "The nondimensional coefficient of the ratio of the "//&
                 "viscosity bounds to the theoretical maximum for "//&
                 "stability without considering other terms.", units="nondim", &
                 default=0.8, do_not_log=.not.(CS%hv_opts_CS%bound_Ah .or. CS%hv_opts_CS%bound_Kh))
  call get_param(param_file, mdl, "KILL_SWITCH_COEF", CS%hv_opts_CS%KS_coef, &
                 "A nondimensional coefficient on the biharmonic viscosity that "// &
                 "sets the kill switch for backscatter. Default is 1.0.", units="nondim", &
                 default=1.0, do_not_log=.not.(CS%hv_opts_CS%EY24_EBT_BS))
  call get_param(param_file, mdl, "NOSLIP", CS%hv_opts_CS%no_slip, &
                 "If true, no slip boundary conditions are used; otherwise "//&
                 "free slip boundary conditions are assumed. The "//&
                 "implementation of the free slip BCs on a C-grid is much "//&
                 "cleaner than the no slip BCs. The use of free slip BCs "//&
                 "is strongly encouraged, and no slip BCs are not used with "//&
                 "the biharmonic viscosity.", default=.false.)
  call get_param(param_file, mdl, "USE_KH_BG_2D", CS%use_Kh_bg_2d, &
                 "If true, read a file containing 2-d background harmonic "//&
                 "viscosities. The final viscosity is the maximum of the other "//&
                 "terms and this background value.", default=.false., do_not_log=.not.CS%hv_opts_CS%Laplacian)
  if (.not.CS%hv_opts_CS%Laplacian) CS%use_Kh_bg_2d = .false.
  call get_param(param_file, mdl, "KH_BG_2D_BUG", CS%Kh_bg_2d_bug, &
                 "If true, retain an answer-changing horizontal indexing bug in setting "//&
                 "the corner-point viscosities when USE_KH_BG_2D=True.  This is "//&
                 "not recommended.", default=.false., do_not_log=.not.CS%use_Kh_bg_2d)
  call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
  call get_param(param_file, mdl, "FRICTWORK_BUG", CS%hv_opts_CS%FrictWork_bug, &
                 "If true, retain an answer-changing bug in calculating the FrictWork, "//&
                 "which cancels the h in thickness flux and the h at velocity point. This is "//&
                 "not recommended.", default=.false.)
  call get_param(param_file, mdl, "OBC_SPECIFIED_STRAIN_BUG", CS%hv_opts_CS%OBC_strain_bug, &
                 "If true, recover a bug that specified shear strain option at open boundaries "//&
                 "cannot be applied.", default=.true.)
  call get_param(param_file, mdl, "USE_GME", CS%hv_opts_CS%use_GME, &
                 "If true, use the GM+E backscatter scheme in association \n"//&
                 "with the Gent and McWilliams parameterization.", default=.false.)
  call get_param(param_file, mdl, "SPLIT", split, &
                 "Use the split time stepping if true.", default=.true., do_not_log=.true.)
  if (CS%hv_opts_CS%use_Leithy) then
    if (.not.(CS%hv_opts_CS%biharmonic .and. CS%hv_opts_CS%Laplacian)) then
                   call MOM_error(FATAL, "MOM_hor_visc.F90, hor_visc_init: "//&
                   "LAPLACIAN and BIHARMONIC must both be True when USE_LEITHY=True.")
    endif
  endif
  call get_param(param_file, mdl, "LEITHY_CK", CS%hv_opts_CS%c_K, &
                 "Fraction of biharmonic dissipation that gets backscattered, "//&
                 "in Leith+E.", units="nondim", default=1.0, do_not_log=.not.CS%hv_opts_CS%use_Leithy)
  call get_param(param_file, mdl, "SMOOTH_AH", CS%hv_opts_CS%smooth_Ah, &
                 "If true, Ah and m_leithy are smoothed within Leith+E.  This requires "//&
                 "lots of blocking communications, which can be expensive", &
                 default=.true., do_not_log=.not.CS%hv_opts_CS%use_Leithy)

  if (CS%hv_opts_CS%use_GME .and. .not.split) call MOM_error(FATAL,"ERROR: Currently, USE_GME = True "// &
                                           "cannot be used with SPLIT=False.")

  if (CS%hv_opts_CS%use_GME) then
    call get_param(param_file, mdl, "GME_NUM_SMOOTHINGS", CS%num_smooth_gme, &
                   "Number of smoothing passes for the GME fluxes.", &
                   default=1)
    call get_param(param_file, mdl, "GME_H0", CS%hv_opts_CS%GME_h0, &
                   "The strength of GME tapers quadratically to zero when the bathymetric "//&
                   "depth is shallower than GME_H0.", &
                   units="m", scale=GV%m_to_H, default=1000.0)
    call get_param(param_file, mdl, "GME_EFFICIENCY", CS%hv_opts_CS%GME_efficiency, &
                   "The nondimensional prefactor multiplying the GME coefficient.", &
                   units="nondim", default=1.0)
    call get_param(param_file, mdl, "GME_LIMITER", CS%hv_opts_CS%GME_limiter, &
                   "The absolute maximum value the GME coefficient is allowed to take.", &
                   units="m2 s-1", scale=US%m_to_L**2*US%T_to_s, default=1.0e7)
  endif

  if (CS%hv_opts_CS%Laplacian .or. CS%hv_opts_CS%biharmonic) then
    call get_param(param_file, mdl, "DT", dt, &
                 "The (baroclinic) dynamics time step.", units="s", scale=US%s_to_T, &
                 fail_if_missing=.true.)
    Idt = 1.0 / dt
  endif
  call get_param(param_file, mdl, "KILL_SWITCH_TIMESCALE", CS%KS_timescale, &
                 "A timescale for computing the CFL limit for viscosity "// &
                 "that determines when backscatter is shut off. Default is DT.", &
                 default= dt ,  units="s", scale=US%s_to_T, do_not_log=.not.(CS%hv_opts_CS%EY24_EBT_BS))

  if (CS%hv_opts_CS%no_slip .and. CS%hv_opts_CS%biharmonic) &
    call MOM_error(FATAL,"ERROR: NOSLIP and BIHARMONIC cannot be defined "// &
                         "at the same time in MOM.")
  if (.not.(CS%hv_opts_CS%Laplacian .or. CS%hv_opts_CS%biharmonic)) then
    ! Only issue inviscid warning if not in single column mode (usually 2x2 domain)
    if ( max(G%domain%niglobal, G%domain%njglobal)>2 ) call MOM_error(WARNING, &
      "hor_visc_init:  It is usually a very bad idea not to use either "//&
      "LAPLACIAN or BIHARMONIC viscosity.")
    return ! We are not using either Laplacian or Bi-harmonic lateral viscosity
  endif

  deg2rad = atan(1.0) / 45.
  ALLOC_(CS%dx2h(isd:ied,jsd:jed))        ; CS%dx2h(:,:)    = 0.0
  ALLOC_(CS%dy2h(isd:ied,jsd:jed))        ; CS%dy2h(:,:)    = 0.0
  ALLOC_(CS%dx2q(IsdB:IedB,JsdB:JedB))    ; CS%dx2q(:,:)    = 0.0
  ALLOC_(CS%dy2q(IsdB:IedB,JsdB:JedB))    ; CS%dy2q(:,:)    = 0.0
  ALLOC_(CS%dx_dyT(isd:ied,jsd:jed))      ; CS%dx_dyT(:,:)  = 0.0
  ALLOC_(CS%dy_dxT(isd:ied,jsd:jed))      ; CS%dy_dxT(:,:)  = 0.0
  ALLOC_(CS%dx_dyBu(IsdB:IedB,JsdB:JedB)) ; CS%dx_dyBu(:,:) = 0.0
  ALLOC_(CS%dy_dxBu(IsdB:IedB,JsdB:JedB)) ; CS%dy_dxBu(:,:) = 0.0

  if (CS%hv_opts_CS%Laplacian) then
    ALLOC_(CS%grid_sp_h2(isd:ied,jsd:jed))   ; CS%grid_sp_h2(:,:) = 0.0
    ALLOC_(CS%Kh_bg_xx(isd:ied,jsd:jed))     ; CS%Kh_bg_xx(:,:) = 0.0
    ALLOC_(CS%Kh_bg_xy(IsdB:IedB,JsdB:JedB)) ; CS%Kh_bg_xy(:,:) = 0.0
    if (CS%hv_opts_CS%bound_Kh .or. CS%hv_opts_CS%EY24_EBT_BS) then
      allocate(CS%Kh_Max_xx(Isd:Ied,Jsd:Jed), source=0.0)
      allocate(CS%Kh_Max_xy(IsdB:IedB,JsdB:JedB), source=0.0)
    endif
    if (CS%hv_opts_CS%Smagorinsky_Kh .or. CS%hv_opts_CS%EY24_EBT_BS) then
      allocate(CS%Laplac2_const_xx(isd:ied,jsd:jed), source=0.0)
      allocate(CS%Laplac2_const_xy(IsdB:IedB,JsdB:JedB), source=0.0)
    endif
    if (CS%hv_opts_CS%Leith_Kh) then
      allocate(CS%Laplac3_const_xx(isd:ied,jsd:jed), source=0.0)
      allocate(CS%Laplac3_const_xy(IsdB:IedB,JsdB:JedB), source=0.0)
    endif
  endif
  ALLOC_(CS%reduction_xx(isd:ied,jsd:jed))     ; CS%reduction_xx(:,:) = 0.0
  ALLOC_(CS%reduction_xy(IsdB:IedB,JsdB:JedB)) ; CS%reduction_xy(:,:) = 0.0

  CS%dynamic_aniso = .false.
  if (CS%hv_opts_CS%anisotropic) then
    allocate(CS%n1n2_h(isd:ied,jsd:jed), source=0.0)
    allocate(CS%n1n1_m_n2n2_h(isd:ied,jsd:jed), source=0.0)
    allocate(CS%n1n2_q(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%n1n1_m_n2n2_q(IsdB:IedB,JsdB:JedB), source=0.0)
    select case (aniso_mode)
      case (0)
        call align_aniso_tensor_to_grid(CS, aniso_grid_dir(1), aniso_grid_dir(2))
      case (1)
      ! call align_aniso_tensor_to_grid(CS, aniso_grid_dir(1), aniso_grid_dir(2))
      case (2)
        CS%dynamic_aniso = .true.
      case default
        call MOM_error(FATAL, "MOM_hor_visc: "//&
             "Runtime parameter ANISOTROPIC_MODE is out of range.")
    end select
  endif

  call get_param(param_file, mdl, "KH_BG_2D_FILENAME", filename, &
                 'The filename containing a 2d map of "Kh".', &
                 default='KH_background_2d.nc', do_not_log=.not.CS%use_Kh_bg_2d)
  call get_param(param_file, mdl, "KH_BG_2D_VARNAME", Kh_var, &
                 'The name in the input file of the horizontal viscosity variable.', &
                 default='Kh', do_not_log=.not.CS%use_Kh_bg_2d)

  if (CS%use_Kh_bg_2d) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    allocate(CS%Kh_bg_2d(isd:ied,jsd:jed), source=0.0)
    call MOM_read_data(trim(inputdir)//trim(filename), Kh_var, CS%Kh_bg_2d, &
                       G%domain, timelevel=1, scale=US%m_to_L**2*US%T_to_s)
    call pass_var(CS%Kh_bg_2d, G%domain)
  endif
  if (CS%hv_opts_CS%biharmonic) then
    ALLOC_(CS%Idx2dyCu(IsdB:IedB,jsd:jed)) ; CS%Idx2dyCu(:,:) = 0.0
    ALLOC_(CS%Idx2dyCv(isd:ied,JsdB:JedB)) ; CS%Idx2dyCv(:,:) = 0.0
    ALLOC_(CS%Idxdy2u(IsdB:IedB,jsd:jed))  ; CS%Idxdy2u(:,:)  = 0.0
    ALLOC_(CS%Idxdy2v(isd:ied,JsdB:JedB))  ; CS%Idxdy2v(:,:)  = 0.0
    ALLOC_(CS%Ah_bg_xx(isd:ied,jsd:jed))     ; CS%Ah_bg_xx(:,:) = 0.0
    ALLOC_(CS%Ah_bg_xy(IsdB:IedB,JsdB:JedB)) ; CS%Ah_bg_xy(:,:) = 0.0
    ALLOC_(CS%grid_sp_h3(isd:ied,jsd:jed))   ; CS%grid_sp_h3(:,:) = 0.0
    if (CS%hv_opts_CS%bound_Ah) then
      allocate(CS%Ah_Max_xx(isd:ied,jsd:jed), source=0.0)
      allocate(CS%Ah_Max_xy(IsdB:IedB,JsdB:JedB), source=0.0)
    endif
    if (CS%hv_opts_CS%EY24_EBT_BS) then
      allocate(CS%Ah_Max_xx_KS(isd:ied,jsd:jed), source=0.0)
      allocate(CS%Ah_Max_xy_KS(IsdB:IedB,JsdB:JedB), source=0.0)
    endif
    if (CS%hv_opts_CS%Smagorinsky_Ah) then
      allocate(CS%Biharm_const_xx(isd:ied,jsd:jed), source=0.0)
      allocate(CS%Biharm_const_xy(IsdB:IedB,JsdB:JedB), source=0.0)
      if (CS%hv_opts_CS%bound_Coriolis) then
        allocate(CS%Biharm_const2_xx(isd:ied,jsd:jed), source=0.0)
        allocate(CS%Biharm_const2_xy(IsdB:IedB,JsdB:JedB), source=0.0)
      endif
    endif
    if ((CS%hv_opts_CS%Leith_Ah) .or. (CS%hv_opts_CS%use_Leithy)) then
      allocate(CS%biharm6_const_xx(isd:ied,jsd:jed), source=0.0)
      allocate(CS%biharm6_const_xy(IsdB:IedB,JsdB:JedB), source=0.0)
    endif
    if (CS%hv_opts_CS%use_Leithy) then
      allocate(CS%m_const_leithy(isd:ied,jsd:jed), source=0.0)
      allocate(CS%m_leithy_max(isd:ied,jsd:jed), source=0.0)
    endif
    if (CS%hv_opts_CS%Re_Ah > 0.0) then
      allocate(CS%Re_Ah_const_xx(isd:ied,jsd:jed), source=0.0)
      allocate(CS%Re_Ah_const_xy(IsdB:IedB,JsdB:JedB), source=0.0)
    endif
  endif

  do concurrent (J=js-2:Jeq+1, I=is-2:Ieq+1)
    CS%dx2q(I,J) = G%dxBu(I,J)*G%dxBu(I,J) ; CS%dy2q(I,J) = G%dyBu(I,J)*G%dyBu(I,J)
  enddo

  if (((CS%hv_opts_CS%Leith_Kh) .or. (CS%hv_opts_CS%Leith_Ah) .or. (CS%hv_opts_CS%use_Leithy)) .and. &
      ((G%isc-G%isd < 3) .or. (G%isc-G%isd < 3))) call MOM_error(FATAL, &
          "The minimum halo size is 3 when a Leith viscosity is being used.")
  if (CS%hv_opts_CS%use_Leithy) then
    do concurrent (J=js-3:Jeq+2, I=is-3:Ieq+2)
      CS%DX_dyBu(I,J) = G%dxBu(I,J)*G%IdyBu(I,J) ; CS%DY_dxBu(I,J) = G%dyBu(I,J)*G%IdxBu(I,J)
    enddo
  elseif ((CS%hv_opts_CS%Leith_Kh) .or. (CS%hv_opts_CS%Leith_Ah)) then
    do concurrent (J=Jsq-2:Jeq+2, I=Isq-2:Ieq+2)
      CS%DX_dyBu(I,J) = G%dxBu(I,J)*G%IdyBu(I,J) ; CS%DY_dxBu(I,J) = G%dyBu(I,J)*G%IdxBu(I,J)
    enddo
  else
    do concurrent (J=js-2:Jeq+1, I=is-2:Ieq+1)
      CS%DX_dyBu(I,J) = G%dxBu(I,J)*G%IdyBu(I,J) ; CS%DY_dxBu(I,J) = G%dyBu(I,J)*G%IdxBu(I,J)
    enddo
  endif

  do concurrent (j=js-2:Jeq+2, i=is-2:Ieq+2)
    CS%dx2h(i,j) = G%dxT(i,j)*G%dxT(i,j) ; CS%dy2h(i,j) = G%dyT(i,j)*G%dyT(i,j)
    CS%DX_dyT(i,j) = G%dxT(i,j)*G%IdyT(i,j) ; CS%DY_dxT(i,j) = G%dyT(i,j)*G%IdxT(i,j)
  enddo

  ! TODO: Remove this after every instance has been moved to GPU

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    CS%reduction_xx(i,j) = 1.0
    if ((G%dy_Cu(I,j) > 0.0) .and. (G%dy_Cu(I,j) < G%dyCu(I,j)) .and. &
        (G%dy_Cu(I,j) < G%dyCu(I,j) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dy_Cu(I,j) / (G%dyCu(I,j))
    if ((G%dy_Cu(I-1,j) > 0.0) .and. (G%dy_Cu(I-1,j) < G%dyCu(I-1,j)) .and. &
        (G%dy_Cu(I-1,j) < G%dyCu(I-1,j) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dy_Cu(I-1,j) / (G%dyCu(I-1,j))
    if ((G%dx_Cv(i,J) > 0.0) .and. (G%dx_Cv(i,J) < G%dxCv(i,J)) .and. &
        (G%dx_Cv(i,J) < G%dxCv(i,J) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dx_Cv(i,J) / (G%dxCv(i,J))
    if ((G%dx_Cv(i,J-1) > 0.0) .and. (G%dx_Cv(i,J-1) < G%dxCv(i,J-1)) .and. &
        (G%dx_Cv(i,J-1) < G%dxCv(i,J-1) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dx_Cv(i,J-1) / (G%dxCv(i,J-1))
  enddo ; enddo

  do J=js-1,Jeq ; do I=is-1,Ieq
    CS%reduction_xy(I,J) = 1.0
    if ((G%dy_Cu(I,j) > 0.0) .and. (G%dy_Cu(I,j) < G%dyCu(I,j)) .and. &
        (G%dy_Cu(I,j) < G%dyCu(I,j) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dy_Cu(I,j) / (G%dyCu(I,j))
    if ((G%dy_Cu(I,j+1) > 0.0) .and. (G%dy_Cu(I,j+1) < G%dyCu(I,j+1)) .and. &
        (G%dy_Cu(I,j+1) < G%dyCu(I,j+1) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dy_Cu(I,j+1) / (G%dyCu(I,j+1))
    if ((G%dx_Cv(i,J) > 0.0) .and. (G%dx_Cv(i,J) < G%dxCv(i,J)) .and. &
        (G%dx_Cv(i,J) < G%dxCv(i,J) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dx_Cv(i,J) / (G%dxCv(i,J))
    if ((G%dx_Cv(i+1,J) > 0.0) .and. (G%dx_Cv(i+1,J) < G%dxCv(i+1,J)) .and. &
        (G%dx_Cv(i+1,J) < G%dxCv(i+1,J) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dx_Cv(i+1,J) / (G%dxCv(i+1,J))
  enddo ; enddo
  if (CS%hv_opts_CS%Laplacian) then
   ! The 0.3 below was 0.4 in MOM1.10.  The change in hq requires
   ! this to be less than 1/3, rather than 1/2 as before.
    if (CS%hv_opts_CS%bound_Kh .or. CS%hv_opts_CS%bound_Ah) Kh_Limit = 0.3 / (dt*4.0)
    ! Calculate and store the background viscosity at h-points

    min_grid_sp_h2 = huge(1.)
    do j=js-1,Jeq+1 ; do i=is-1,Ieq+1
      ! Static factors in the Smagorinsky and Leith schemes
      grid_sp_h2 = (2.0*CS%dx2h(i,j)*CS%dy2h(i,j)) / (CS%dx2h(i,j) + CS%dy2h(i,j))
      CS%grid_sp_h2(i,j) = grid_sp_h2
      grid_sp_h3 = grid_sp_h2*sqrt(grid_sp_h2)
      if (CS%hv_opts_CS%Smagorinsky_Kh) CS%Laplac2_const_xx(i,j) = Smag_Lap_const * grid_sp_h2
      if (CS%hv_opts_CS%Leith_Kh)       CS%Laplac3_const_xx(i,j) = Leith_Lap_const * grid_sp_h3
      ! Maximum of constant background and MICOM viscosity
      CS%Kh_bg_xx(i,j) = MAX(Kh, Kh_vel_scale * sqrt(grid_sp_h2))
      ! Use the larger of the above and values read from a file
      if (CS%use_Kh_bg_2d) CS%Kh_bg_xx(i,j) = MAX(CS%Kh_bg_2d(i,j), CS%Kh_bg_xx(i,j))
      ! Use the larger of the above and a function of sin(latitude)
      if (Kh_sin_lat>0.) then
        slat_fn = abs( sin( deg2rad * G%geoLatT(i,j) ) ) ** Kh_pwr_of_sine
        CS%Kh_bg_xx(i,j) = MAX(Kh_sin_lat * slat_fn, CS%Kh_bg_xx(i,j))
      endif
      min_grid_sp_h2 = min(grid_sp_h2, min_grid_sp_h2)
    enddo ; enddo
    call min_across_PEs(min_grid_sp_h2)

    ! Calculate and store the background viscosity at q-points
    do J=js-1,Jeq ; do I=is-1,Ieq
      ! Static factors in the Smagorinsky and Leith schemes
      grid_sp_q2 = (2.0*CS%dx2q(I,J)*CS%dy2q(I,J)) / (CS%dx2q(I,J) + CS%dy2q(I,J))
      grid_sp_q3 = grid_sp_q2*sqrt(grid_sp_q2)
      if (CS%hv_opts_CS%Smagorinsky_Kh) CS%Laplac2_const_xy(I,J) = Smag_Lap_const * grid_sp_q2
      if (CS%hv_opts_CS%Leith_Kh)       CS%Laplac3_const_xy(I,J) = Leith_Lap_const * grid_sp_q3
      ! Maximum of constant background and MICOM viscosity
      CS%Kh_bg_xy(I,J) = MAX(Kh, Kh_vel_scale * sqrt(grid_sp_q2))
      ! Use the larger of the above and values read from a file
      if (CS%use_Kh_bg_2d) then
        if (CS%Kh_bg_2d_bug) then
          ! This option is unambiguously wrong but is needed to recover old answers
          CS%Kh_bg_xy(I,J) = MAX(CS%Kh_bg_2d(i,j), CS%Kh_bg_xy(I,J))
        else
          CS%Kh_bg_xy(I,J) = MAX(CS%Kh_bg_xy(I,J), &
              0.25*((CS%Kh_bg_2d(i,j) + CS%Kh_bg_2d(i+1,j+1)) + &
                    (CS%Kh_bg_2d(i+1,j) + CS%Kh_bg_2d(i,j+1))) )
        endif
      endif

      ! Use the larger of the above and a function of sin(latitude)
      if (Kh_sin_lat>0.) then
        slat_fn = abs( sin( deg2rad * G%geoLatBu(I,J) ) ) ** Kh_pwr_of_sine
        CS%Kh_bg_xy(I,J) = MAX(Kh_sin_lat * slat_fn, CS%Kh_bg_xy(I,J))
      endif
    enddo ; enddo
  endif
  if (CS%hv_opts_CS%biharmonic) then
    do j=js-1,Jeq+1 ; do I=is-2,Ieq+1
      CS%Idx2dyCu(I,j) = (G%IdxCu(I,j)*G%IdxCu(I,j)) * G%IdyCu(I,j)
      CS%Idxdy2u(I,j) = G%IdxCu(I,j) * (G%IdyCu(I,j)*G%IdyCu(I,j))
    enddo ; enddo
    do J=js-2,Jeq+1 ; do i=is-1,Ieq+1
      CS%Idx2dyCv(i,J) = (G%IdxCv(i,J)*G%IdxCv(i,J)) * G%IdyCv(i,J)
      CS%Idxdy2v(i,J) = G%IdxCv(i,J) * (G%IdyCv(i,J)*G%IdyCv(i,J))
    enddo ; enddo
    CS%Ah_bg_xy(:,:) = 0.0
    ! The 0.3 below was 0.4 in HIM 1.10.  The change in hq requires
    ! this to be less than 1/3, rather than 1/2 as before.
    if (CS%hv_opts_CS%bound_Ah) Ah_Limit = 0.3 / (dt*64.0)
    if (CS%hv_opts_CS%Smagorinsky_Ah .and. CS%hv_opts_CS%bound_Coriolis) &
      BoundCorConst = 1.0 / (5.0*(bound_Cor_vel*bound_Cor_vel))

    min_grid_sp_h4 = huge(1.)
    do j=js-1,Jeq+1 ; do i=is-1,Ieq+1
      grid_sp_h2 = (2.0*CS%dx2h(i,j)*CS%dy2h(i,j)) / (CS%dx2h(i,j)+CS%dy2h(i,j))
      grid_sp_h3 = grid_sp_h2*sqrt(grid_sp_h2)
      CS%grid_sp_h3(i,j) = grid_sp_h3
      if (CS%hv_opts_CS%Smagorinsky_Ah) then
        CS%Biharm_const_xx(i,j) = Smag_bi_const * (grid_sp_h2 * grid_sp_h2)
        if (CS%hv_opts_CS%bound_Coriolis) then
          fmax = MAX(abs(G%CoriolisBu(I-1,J-1)), abs(G%CoriolisBu(I,J-1)), &
                     abs(G%CoriolisBu(I-1,J)),   abs(G%CoriolisBu(I,J)))
          CS%Biharm_const2_xx(i,j) = (grid_sp_h2 * grid_sp_h2 * grid_sp_h2) * &
                                     (fmax * BoundCorConst)
        endif
      endif
      if (CS%hv_opts_CS%Leith_Ah) then
        CS%biharm6_const_xx(i,j) = Leith_bi_const * (grid_sp_h3 * grid_sp_h3)
      endif
      if (CS%hv_opts_CS%use_Leithy) then
        CS%biharm6_const_xx(i,j) = Leith_bi_const * max(G%dxT(i,j),G%dyT(i,j))**6
        CS%m_const_leithy(i,j) = 0.5 * sqrt(CS%hv_opts_CS%c_K) * max(G%dxT(i,j),G%dyT(i,j))
        CS%m_leithy_max(i,j) = 4. / max(G%dxT(i,j),G%dyT(i,j))**2
      endif
      CS%Ah_bg_xx(i,j) = MAX(Ah, Ah_vel_scale * grid_sp_h2 * sqrt(grid_sp_h2))
      if (CS%hv_opts_CS%Re_Ah > 0.0) CS%Re_Ah_const_xx(i,j) = grid_sp_h3 / CS%hv_opts_CS%Re_Ah
      if (Ah_time_scale > 0.) CS%Ah_bg_xx(i,j) = &
            MAX(CS%Ah_bg_xx(i,j), (grid_sp_h2 * grid_sp_h2) / Ah_time_scale)
      min_grid_sp_h4 = min(grid_sp_h2**2, min_grid_sp_h4)
    enddo ; enddo
    call min_across_PEs(min_grid_sp_h4)

    do J=js-1,Jeq ; do I=is-1,Ieq
      grid_sp_q2 = (2.0*CS%dx2q(I,J)*CS%dy2q(I,J)) / (CS%dx2q(I,J)+CS%dy2q(I,J))
      grid_sp_q3 = grid_sp_q2*sqrt(grid_sp_q2)
      if (CS%hv_opts_CS%Smagorinsky_Ah) then
        CS%Biharm_const_xy(I,J) = Smag_bi_const * (grid_sp_q2 * grid_sp_q2)
        if (CS%hv_opts_CS%bound_Coriolis) then
          CS%Biharm_const2_xy(I,J) = (grid_sp_q2 * grid_sp_q2 * grid_sp_q2) * &
                                     (abs(G%CoriolisBu(I,J)) * BoundCorConst)
        endif
      endif
      if ((CS%hv_opts_CS%Leith_Ah) .or. (CS%hv_opts_CS%use_Leithy))then
        CS%biharm6_const_xy(I,J) = Leith_bi_const * (grid_sp_q3 * grid_sp_q3)
      endif
      CS%Ah_bg_xy(I,J) = MAX(Ah, Ah_vel_scale * grid_sp_q2 * sqrt(grid_sp_q2))
      if (CS%hv_opts_CS%Re_Ah > 0.0) CS%Re_Ah_const_xy(i,j) = grid_sp_q3 / CS%hv_opts_CS%Re_Ah
      if (Ah_time_scale > 0.) CS%Ah_bg_xy(i,j) = &
            MAX(CS%Ah_bg_xy(i,j), (grid_sp_q2 * grid_sp_q2) / Ah_time_scale)
    enddo ; enddo
  endif
  ! The Laplacian bounds should avoid overshoots when CS%bound_coef < 1.
  if (CS%hv_opts_CS%Laplacian .and. CS%hv_opts_CS%bound_Kh) then
    do j=js-1,Jeq+1 ; do i=is-1,Ieq+1
      denom = max( &
         (CS%dy2h(i,j) * CS%DY_dxT(i,j) * (G%IdyCu(I,j) + G%IdyCu(I-1,j)) * &
          max(G%IdyCu(I,j)*G%IareaCu(I,j), G%IdyCu(I-1,j)*G%IareaCu(I-1,j)) ), &
         (CS%dx2h(i,j) * CS%DX_dyT(i,j) * (G%IdxCv(i,J) + G%IdxCv(i,J-1)) * &
          max(G%IdxCv(i,J)*G%IareaCv(i,J), G%IdxCv(i,J-1)*G%IareaCv(i,J-1)) ) )
      CS%Kh_Max_xx(i,j) = 0.0
      if (denom > 0.0) &
        CS%Kh_Max_xx(i,j) = CS%bound_coef * 0.25 * Idt / denom
    enddo ; enddo
    do J=js-1,Jeq ; do I=is-1,Ieq
      denom = max( &
         (CS%dx2q(I,J) * CS%DX_dyBu(I,J) * (G%IdxCu(I,j+1) + G%IdxCu(I,j)) * &
          max(G%IdxCu(I,j)*G%IareaCu(I,j), G%IdxCu(I,j+1)*G%IareaCu(I,j+1)) ), &
         (CS%dy2q(I,J) * CS%DY_dxBu(I,J) * (G%IdyCv(i+1,J) + G%IdyCv(i,J)) * &
          max(G%IdyCv(i,J)*G%IareaCv(i,J), G%IdyCv(i+1,J)*G%IareaCv(i+1,J)) ) )
      CS%Kh_Max_xy(I,J) = 0.0
      if (denom > 0.0) &
        CS%Kh_Max_xy(I,J) = CS%bound_coef * 0.25 * Idt / denom
    enddo ; enddo
    if (CS%hv_opts_CS%debug) then
      call hchksum(CS%Kh_Max_xx, "Kh_Max_xx", G%HI, haloshift=0, unscale=US%L_to_m**2*US%s_to_T)
      call Bchksum(CS%Kh_Max_xy, "Kh_Max_xy", G%HI, haloshift=0, unscale=US%L_to_m**2*US%s_to_T)
    endif
  endif
  ! The biharmonic bounds should avoid overshoots when CS%bound_coef < 0.5, but
  ! empirically work for CS%bound_coef <~ 1.0
  if (CS%hv_opts_CS%biharmonic .and. CS%hv_opts_CS%bound_Ah) then
    do j=js-1,Jeq+1 ; do I=is-2,Ieq+1
      u0u(I,j) = ((CS%Idxdy2u(I,j)*((CS%dy2h(i+1,j)*CS%DY_dxT(i+1,j)*(G%IdyCu(I+1,j) + G%IdyCu(I,j))) + &
                                   (CS%dy2h(i,j) * CS%DY_dxT(i,j) * (G%IdyCu(I,j) + G%IdyCu(I-1,j))) )) + &
                  (CS%Idx2dyCu(I,j)*((CS%dx2q(I,J) * CS%DX_dyBu(I,J) * (G%IdxCu(I,j+1) + G%IdxCu(I,j))) + &
                                   (CS%dx2q(I,J-1)*CS%DX_dyBu(I,J-1)*(G%IdxCu(I,j) + G%IdxCu(I,j-1))) )) )
      u0v(I,j) = ((CS%Idxdy2u(I,j)*((CS%dy2h(i+1,j)*CS%DX_dyT(i+1,j)*(G%IdxCv(i+1,J) + G%IdxCv(i+1,J-1))) + &
                                   (CS%dy2h(i,j) * CS%DX_dyT(i,j) * (G%IdxCv(i,J) + G%IdxCv(i,J-1))) )) + &
                  (CS%Idx2dyCu(I,j)*((CS%dx2q(I,J) * CS%DY_dxBu(I,J) * (G%IdyCv(i+1,J) + G%IdyCv(i,J))) + &
                                   (CS%dx2q(I,J-1)*CS%DY_dxBu(I,J-1)*(G%IdyCv(i+1,J-1) + G%IdyCv(i,J-1))) )) )
    enddo ; enddo
    do J=js-2,Jeq+1 ; do i=is-1,Ieq+1
      v0u(i,J) = ((CS%Idxdy2v(i,J)*((CS%dy2q(I,J) * CS%DX_dyBu(I,J) * (G%IdxCu(I,j+1) + G%IdxCu(I,j))) + &
                                    (CS%dy2q(I-1,J)*CS%DX_dyBu(I-1,J)*(G%IdxCu(I-1,j+1) + G%IdxCu(I-1,j))) )) + &
                  (CS%Idx2dyCv(i,J)*((CS%dx2h(i,j+1)*CS%DY_dxT(i,j+1)*(G%IdyCu(I,j+1) + G%IdyCu(I-1,j+1))) + &
                                    (CS%dx2h(i,j) * CS%DY_dxT(i,j) * (G%IdyCu(I,j) + G%IdyCu(I-1,j))) ) ))
      v0v(i,J) = ((CS%Idxdy2v(i,J)*((CS%dy2q(I,J) * CS%DY_dxBu(I,J) * (G%IdyCv(i+1,J) + G%IdyCv(i,J))) + &
                                    (CS%dy2q(I-1,J)*CS%DY_dxBu(I-1,J)*(G%IdyCv(i,J) + G%IdyCv(i-1,J))) )) + &
                  (CS%Idx2dyCv(i,J)*((CS%dx2h(i,j+1)*CS%DX_dyT(i,j+1)*(G%IdxCv(i,J+1) + G%IdxCv(i,J))) + &
                                    (CS%dx2h(i,j) * CS%DX_dyT(i,j) * (G%IdxCv(i,J) + G%IdxCv(i,J-1))) )) )
    enddo ; enddo
    do j=js-1,Jeq+1 ; do i=is-1,Ieq+1
      denom = max( &
         (CS%dy2h(i,j) * &
          ((CS%DY_dxT(i,j)*((G%IdyCu(I,j)*u0u(I,j)) + (G%IdyCu(I-1,j)*u0u(I-1,j))))  + &
           (CS%DX_dyT(i,j)*((G%IdxCv(i,J)*v0u(i,J)) + (G%IdxCv(i,J-1)*v0u(i,J-1))))) * &
          max(G%IdyCu(I,j)*G%IareaCu(I,j), G%IdyCu(I-1,j)*G%IareaCu(I-1,j)) ), &
         (CS%dx2h(i,j) * &
          ((CS%DY_dxT(i,j)*((G%IdyCu(I,j)*u0v(I,j)) + (G%IdyCu(I-1,j)*u0v(I-1,j))))  + &
           (CS%DX_dyT(i,j)*((G%IdxCv(i,J)*v0v(i,J)) + (G%IdxCv(i,J-1)*v0v(i,J-1))))) * &
          max(G%IdxCv(i,J)*G%IareaCv(i,J), G%IdxCv(i,J-1)*G%IareaCv(i,J-1)) ) )
      CS%Ah_Max_xx(I,J) = 0.0
      if (denom > 0.0) then
        CS%Ah_Max_xx(I,J) = CS%bound_coef * 0.5 * Idt / denom
        if (CS%hv_opts_CS%EY24_EBT_BS) then
          CS%Ah_Max_xx_KS(i,j) = CS%bound_coef * 0.5 / (CS%KS_timescale * denom)
        endif
      endif

    enddo ; enddo
    do J=js-1,Jeq ; do I=is-1,Ieq
      denom = max( &
         (CS%dx2q(I,J) * &
          ((CS%DX_dyBu(I,J)*((u0u(I,j+1)*G%IdxCu(I,j+1)) + (u0u(I,j)*G%IdxCu(I,j))))  + &
           (CS%DY_dxBu(I,J)*((v0u(i+1,J)*G%IdyCv(i+1,J)) + (v0u(i,J)*G%IdyCv(i,J))))) * &
          max(G%IdxCu(I,j)*G%IareaCu(I,j), G%IdxCu(I,j+1)*G%IareaCu(I,j+1)) ), &
         (CS%dy2q(I,J) * &
          ((CS%DX_dyBu(I,J)*((u0v(I,j+1)*G%IdxCu(I,j+1)) + (u0v(I,j)*G%IdxCu(I,j))))  + &
           (CS%DY_dxBu(I,J)*((v0v(i+1,J)*G%IdyCv(i+1,J)) + (v0v(i,J)*G%IdyCv(i,J))))) * &
          max(G%IdyCv(i,J)*G%IareaCv(i,J), G%IdyCv(i+1,J)*G%IareaCv(i+1,J)) ) )
      CS%Ah_Max_xy(I,J) = 0.0
      if (denom > 0.0) then
        CS%Ah_Max_xy(I,J) = CS%bound_coef * 0.5 * Idt / denom
        if (CS%hv_opts_CS%EY24_EBT_BS) then
          CS%Ah_Max_xy_KS(i,j) = CS%bound_coef * 0.5 / (CS%KS_timescale * denom)
        endif
      endif

    enddo ; enddo
    if (CS%hv_opts_CS%debug) then
      call hchksum(CS%Ah_Max_xx, "Ah_Max_xx", G%HI, haloshift=0, unscale=US%L_to_m**4*US%s_to_T)
      call Bchksum(CS%Ah_Max_xy, "Ah_Max_xy", G%HI, haloshift=0, unscale=US%L_to_m**4*US%s_to_T)
    endif
  endif
  ! Register fields for output from this module.
  CS%id_normstress = register_diag_field('ocean_model', 'NoSt', diag%axesTL, Time, &
      'Normal Stress', 's-1', conversion=US%s_to_T)
  CS%id_shearstress = register_diag_field('ocean_model', 'ShSt', diag%axesBL, Time, &
      'Shear Stress', 's-1', conversion=US%s_to_T)
  CS%id_diffu = register_diag_field('ocean_model', 'diffu', diag%axesCuL, Time, &
      'Zonal Acceleration from Horizontal Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_diffv = register_diag_field('ocean_model', 'diffv', diag%axesCvL, Time, &
      'Meridional Acceleration from Horizontal Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)

  !CS%id_hf_diffu = register_diag_field('ocean_model', 'hf_diffu', diag%axesCuL, Time, &
  !    'Fractional Thickness-weighted Zonal Acceleration from Horizontal Viscosity', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if ((CS%id_hf_diffu > 0) .and. (present(ADp))) then
  !  call safe_alloc_alloc(CS%hf_diffu,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  !  call safe_alloc_ptr(ADp%diag_hfrac_u,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  !endif

  !CS%id_hf_diffv = register_diag_field('ocean_model', 'hf_diffv', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Meridional Acceleration from Horizontal Viscosity', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if ((CS%id_hf_diffv > 0) .and. (present(ADp))) then
  !  call safe_alloc_alloc(CS%hf_diffv,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  !  call safe_alloc_ptr(ADp%diag_hfrac_v,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  !endif

  CS%id_hf_diffu_2d = register_diag_field('ocean_model', 'hf_diffu_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Horizontal Viscosity', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_hf_diffu_2d > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hfrac_u,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  endif

  CS%id_hf_diffv_2d = register_diag_field('ocean_model', 'hf_diffv_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Horizontal Viscosity', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_hf_diffv_2d > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hfrac_v,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  endif

  CS%id_h_diffu = register_diag_field('ocean_model', 'h_diffu', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if ((CS%id_h_diffu > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hu,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  endif

  CS%id_h_diffv = register_diag_field('ocean_model', 'h_diffv', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if ((CS%id_h_diffv > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hv,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  endif

  CS%id_intz_diffu_2d = register_diag_field('ocean_model', 'intz_diffu_2d', diag%axesCu1, Time, &
      'Depth-integral of Zonal Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if ((CS%id_intz_diffu_2d > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hu,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  endif

  CS%id_intz_diffv_2d = register_diag_field('ocean_model', 'intz_diffv_2d', diag%axesCv1, Time, &
      'Depth-integral of Meridional Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if ((CS%id_intz_diffv_2d > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hv,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  endif

  CS%id_diffu_visc_rem = register_diag_field('ocean_model', 'diffu_visc_rem', diag%axesCuL, Time, &
      'Zonal Acceleration from Horizontal Viscosity multiplied by viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_diffu_visc_rem > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%visc_rem_u,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  endif

  CS%id_diffv_visc_rem = register_diag_field('ocean_model', 'diffv_visc_rem', diag%axesCvL, Time, &
      'Meridional Acceleration from Horizontal Viscosity multiplied by viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_diffv_visc_rem > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%visc_rem_v,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  endif

  if (CS%hv_opts_CS%biharmonic) then
    CS%id_Ah_h = register_diag_field('ocean_model', 'Ahh', diag%axesTL, Time,    &
        'Biharmonic Horizontal Viscosity at h Points', 'm4 s-1', conversion=US%L_to_m**4*US%s_to_T, &
        cmor_field_name='difmxybo',                                             &
        cmor_long_name='Ocean lateral biharmonic viscosity',                     &
        cmor_standard_name='ocean_momentum_xy_biharmonic_diffusivity')
    CS%id_Ah_q = register_diag_field('ocean_model', 'Ahq', diag%axesBL, Time, &
        'Biharmonic Horizontal Viscosity at q Points', 'm4 s-1', conversion=US%L_to_m**4*US%s_to_T)
    CS%id_grid_Re_Ah = register_diag_field('ocean_model', 'grid_Re_Ah', diag%axesTL, Time, &
        'Grid Reynolds number for the Biharmonic horizontal viscosity at h points', 'nondim')
    if (CS%hv_opts_CS%EY24_EBT_BS) then
      CS%id_visc_limit_h_flag = register_diag_field('ocean_model', 'visc_limit_h_flag', diag%axesTL, Time, &
          'Locations where the biharmonic viscosity reached the better_bound limiter at h points', 'nondim')
      CS%id_visc_limit_q_flag = register_diag_field('ocean_model', 'visc_limit_q_flag', diag%axesBL, Time, &
          'Locations where the biharmonic viscosity reached the better_bound limiter at q points', 'nondim')
      CS%id_visc_limit_h = register_diag_field('ocean_model', 'visc_limit_h', diag%axesTL, Time, &
          'Value of the biharmonic viscosity limiter at h points', 'nondim')
      CS%id_visc_limit_q = register_diag_field('ocean_model', 'visc_limit_q', diag%axesBL, Time, &
          'Value of the biharmonic viscosity limiter at q points', 'nondim')
      CS%id_visc_limit_h_frac = register_diag_field('ocean_model', 'visc_limit_h_frac', diag%axesTL, Time, &
          'Value of the biharmonic viscosity limiter at h points', 'nondim')
      CS%id_visc_limit_q_frac = register_diag_field('ocean_model', 'visc_limit_q_frac', diag%axesBL, Time, &
          'Value of the biharmonic viscosity limiter at q points', 'nondim')
    endif

    if (CS%id_grid_Re_Ah > 0) &
      ! Compute the smallest biharmonic viscosity capable of modifying the
      ! velocity at floating point precision.
      CS%hv_opts_CS%min_grid_Ah = spacing(1.) * min_grid_sp_h4 * Idt
  endif
  if (CS%hv_opts_CS%Laplacian) then
    CS%id_Kh_h = register_diag_field('ocean_model', 'Khh', diag%axesTL, Time,   &
        'Laplacian Horizontal Viscosity at h Points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T, &
        cmor_field_name='difmxylo',                                             &
        cmor_long_name='Ocean lateral Laplacian viscosity',                     &
        cmor_standard_name='ocean_momentum_xy_laplacian_diffusivity')
    CS%id_Kh_q = register_diag_field('ocean_model', 'Khq', diag%axesBL, Time, &
        'Laplacian Horizontal Viscosity at q Points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
    CS%id_grid_Re_Kh = register_diag_field('ocean_model', 'grid_Re_Kh', diag%axesTL, Time, &
        'Grid Reynolds number for the Laplacian horizontal viscosity at h points', 'nondim')
    CS%id_vort_xy_q = register_diag_field('ocean_model', 'vort_xy_q', diag%axesBL, Time, &
      'Vertical vorticity at q Points', 's-1', conversion=US%s_to_T)
    CS%id_div_xx_h = register_diag_field('ocean_model', 'div_xx_h', diag%axesTL, Time, &
      'Horizontal divergence at h Points', 's-1', conversion=US%s_to_T)
    CS%id_sh_xy_q = register_diag_field('ocean_model', 'sh_xy_q', diag%axesBL, Time, &
      'Shearing strain at q Points', 's-1', conversion=US%s_to_T)
    CS%id_sh_xx_h = register_diag_field('ocean_model', 'sh_xx_h', diag%axesTL, Time, &
      'Horizontal tension at h Points', 's-1', conversion=US%s_to_T)

    if (CS%id_grid_Re_Kh > 0) &
      ! Compute a smallest Laplacian viscosity capable of modifying the
      ! velocity at floating point precision.
      CS%hv_opts_CS%min_grid_Kh = spacing(1.) * min_grid_sp_h2 * Idt
  endif
  if (CS%hv_opts_CS%use_GME) then
    CS%id_dudx_bt = register_diag_field('ocean_model', 'dudx_bt', diag%axesT1, Time, &
        'Zonal component of the barotropic shearing strain at h points', 's-1', &
        conversion=US%s_to_T)
    CS%id_dudy_bt = register_diag_field('ocean_model', 'dudy_bt', diag%axesB1, Time, &
        'Zonal component of the barotropic shearing strain at q points', 's-1', &
        conversion=US%s_to_T)
    CS%id_dvdy_bt = register_diag_field('ocean_model', 'dvdy_bt', diag%axesT1, Time, &
        'Meridional component of the barotropic shearing strain at h points', 's-1', &
        conversion=US%s_to_T)
    CS%id_dvdx_bt = register_diag_field('ocean_model', 'dvdx_bt', diag%axesB1, Time, &
        'Meridional component of the barotropic shearing strain at q points', 's-1', &
        conversion=US%s_to_T)
    CS%id_GME_coeff_h = register_diag_field('ocean_model', 'GME_coeff_h', diag%axesTL, Time, &
        'GME coefficient at h Points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
    CS%id_GME_coeff_q = register_diag_field('ocean_model', 'GME_coeff_q', diag%axesBL, Time, &
        'GME coefficient at q Points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
    CS%id_FrictWork_GME = register_diag_field('ocean_model','FrictWork_GME',diag%axesTL,Time,&
        'Integral work done by lateral friction terms in GME (excluding diffusion of energy)', &
        'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  endif

  if (CS%hv_opts_CS%EY24_EBT_BS) then
    CS%id_BS_coeff_h = register_diag_field('ocean_model', 'BS_coeff_h', diag%axesTL, Time, &
        'Backscatter coefficient at h points', units='m2 s-1', conversion=US%L_to_m**2*US%s_to_T)
    CS%id_BS_coeff_q = register_diag_field('ocean_model', 'BS_coeff_q', diag%axesBL, Time, &
        'Backscatter coefficient at q points', units='m2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  endif

  CS%id_FrictWork = register_diag_field('ocean_model','FrictWork',diag%axesTL,Time,&
      'Integral work done by lateral friction terms. If GME is turned on, this '//&
      'includes the GME contribution.', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  CS%id_FrictWorkIntz = register_diag_field('ocean_model','FrictWorkIntz',diag%axesT1,Time,      &
      'Depth integrated work done by lateral friction', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2, &
      cmor_field_name='dispkexyfo',                                                              &
      cmor_long_name='Depth integrated ocean kinetic energy dissipation due to lateral friction',&
      cmor_standard_name='ocean_kinetic_energy_dissipation_per_unit_area_due_to_xy_friction')
  CS%id_FrictWork_bh = register_diag_field('ocean_model','FrictWork_bh',diag%axesTL,Time,&
      'Integral work done by the biharmonic lateral friction terms.', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  CS%id_FrictWorkIntz_bh = register_diag_field('ocean_model','FrictWorkIntz_bh',diag%axesT1,Time,&
      'Depth integrated work done by the biharmonic lateral friction', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)


  ! TODO: Position these after their respective loops (and run loops on device)


  ! Persistent copies of the metric and coefficient arrays used by the horizontal viscosity kernels.
  if (CS%hv_opts_CS%Laplacian) &
    call CS%hv_metrics_CS%Kh_bg_xx%alloc(lb=LBOUND(CS%Kh_bg_xx), ub=UBOUND(CS%Kh_bg_xx), &
                                         source=CS%Kh_bg_xx)
  if (CS%hv_opts_CS%biharmonic) &
    call CS%hv_metrics_CS%Ah_bg_xx%alloc(lb=LBOUND(CS%Ah_bg_xx), ub=UBOUND(CS%Ah_bg_xx), &
                                         source=CS%Ah_bg_xx)
  call CS%hv_metrics_CS%reduction_xx%alloc(lb=LBOUND(CS%reduction_xx), ub=UBOUND(CS%reduction_xx), &
                                           source=CS%reduction_xx)
  if (allocated(CS%Kh_Max_xx)) &
    call CS%hv_metrics_CS%Kh_Max_xx%alloc(lb=LBOUND(CS%Kh_Max_xx), ub=UBOUND(CS%Kh_Max_xx), &
                                          source=CS%Kh_Max_xx)
  if (allocated(CS%Ah_Max_xx)) &
    call CS%hv_metrics_CS%Ah_Max_xx%alloc(lb=LBOUND(CS%Ah_Max_xx), ub=UBOUND(CS%Ah_Max_xx), &
                                          source=CS%Ah_Max_xx)
  if (allocated(CS%Ah_Max_xx_KS)) &
    call CS%hv_metrics_CS%Ah_Max_xx_KS%alloc(lb=LBOUND(CS%Ah_Max_xx_KS), ub=UBOUND(CS%Ah_Max_xx_KS), &
                                             source=CS%Ah_Max_xx_KS)
  if (allocated(CS%n1n2_h)) &
    call CS%hv_metrics_CS%n1n2_h%alloc(lb=LBOUND(CS%n1n2_h), ub=UBOUND(CS%n1n2_h), &
                                       source=CS%n1n2_h)
  if (allocated(CS%n1n1_m_n2n2_h)) &
    call CS%hv_metrics_CS%n1n1_m_n2n2_h%alloc(lb=LBOUND(CS%n1n1_m_n2n2_h), ub=UBOUND(CS%n1n1_m_n2n2_h), &
                                              source=CS%n1n1_m_n2n2_h)
  if (CS%hv_opts_CS%Laplacian) &
    call CS%hv_metrics_CS%grid_sp_h2%alloc(lb=LBOUND(CS%grid_sp_h2), ub=UBOUND(CS%grid_sp_h2), &
                                           source=CS%grid_sp_h2)
  if (CS%hv_opts_CS%biharmonic) &
    call CS%hv_metrics_CS%grid_sp_h3%alloc(lb=LBOUND(CS%grid_sp_h3), ub=UBOUND(CS%grid_sp_h3), &
                                           source=CS%grid_sp_h3)
  if (CS%hv_opts_CS%Laplacian) &
    call CS%hv_metrics_CS%Kh_bg_xy%alloc(lb=LBOUND(CS%Kh_bg_xy), ub=UBOUND(CS%Kh_bg_xy), &
                                         source=CS%Kh_bg_xy)
  if (CS%hv_opts_CS%biharmonic) &
    call CS%hv_metrics_CS%Ah_bg_xy%alloc(lb=LBOUND(CS%Ah_bg_xy), ub=UBOUND(CS%Ah_bg_xy), &
                                         source=CS%Ah_bg_xy)
  call CS%hv_metrics_CS%reduction_xy%alloc(lb=LBOUND(CS%reduction_xy), ub=UBOUND(CS%reduction_xy), &
                                           source=CS%reduction_xy)
  if (allocated(CS%Kh_Max_xy)) &
    call CS%hv_metrics_CS%Kh_Max_xy%alloc(lb=LBOUND(CS%Kh_Max_xy), ub=UBOUND(CS%Kh_Max_xy), &
                                          source=CS%Kh_Max_xy)
  if (allocated(CS%Ah_Max_xy)) &
    call CS%hv_metrics_CS%Ah_Max_xy%alloc(lb=LBOUND(CS%Ah_Max_xy), ub=UBOUND(CS%Ah_Max_xy), &
                                          source=CS%Ah_Max_xy)
  if (allocated(CS%Ah_Max_xy_KS)) &
    call CS%hv_metrics_CS%Ah_Max_xy_KS%alloc(lb=LBOUND(CS%Ah_Max_xy_KS), ub=UBOUND(CS%Ah_Max_xy_KS), &
                                             source=CS%Ah_Max_xy_KS)
  if (allocated(CS%n1n2_q)) &
    call CS%hv_metrics_CS%n1n2_q%alloc(lb=LBOUND(CS%n1n2_q), ub=UBOUND(CS%n1n2_q), &
                                       source=CS%n1n2_q)
  if (allocated(CS%n1n1_m_n2n2_q)) &
    call CS%hv_metrics_CS%n1n1_m_n2n2_q%alloc(lb=LBOUND(CS%n1n1_m_n2n2_q), ub=UBOUND(CS%n1n1_m_n2n2_q), &
                                              source=CS%n1n1_m_n2n2_q)
  call CS%hv_metrics_CS%dx2h%alloc(lb=LBOUND(CS%dx2h), ub=UBOUND(CS%dx2h), &
                                   source=CS%dx2h)
  call CS%hv_metrics_CS%dy2h%alloc(lb=LBOUND(CS%dy2h), ub=UBOUND(CS%dy2h), &
                                   source=CS%dy2h)
  call CS%hv_metrics_CS%dx_dyT%alloc(lb=LBOUND(CS%dx_dyT), ub=UBOUND(CS%dx_dyT), &
                                     source=CS%dx_dyT)
  call CS%hv_metrics_CS%dy_dxT%alloc(lb=LBOUND(CS%dy_dxT), ub=UBOUND(CS%dy_dxT), &
                                     source=CS%dy_dxT)
  if (allocated(CS%m_const_leithy)) &
    call CS%hv_metrics_CS%m_const_leithy%alloc(lb=LBOUND(CS%m_const_leithy), ub=UBOUND(CS%m_const_leithy), &
                                               source=CS%m_const_leithy)
  if (allocated(CS%m_leithy_max)) &
    call CS%hv_metrics_CS%m_leithy_max%alloc(lb=LBOUND(CS%m_leithy_max), ub=UBOUND(CS%m_leithy_max), &
                                             source=CS%m_leithy_max)
  call CS%hv_metrics_CS%dx2q%alloc(lb=LBOUND(CS%dx2q), ub=UBOUND(CS%dx2q), &
                                   source=CS%dx2q)
  call CS%hv_metrics_CS%dy2q%alloc(lb=LBOUND(CS%dy2q), ub=UBOUND(CS%dy2q), &
                                   source=CS%dy2q)
  call CS%hv_metrics_CS%dx_dyBu%alloc(lb=LBOUND(CS%dx_dyBu), ub=UBOUND(CS%dx_dyBu), &
                                      source=CS%dx_dyBu)
  call CS%hv_metrics_CS%dy_dxBu%alloc(lb=LBOUND(CS%dy_dxBu), ub=UBOUND(CS%dy_dxBu), &
                                      source=CS%dy_dxBu)
  if (CS%hv_opts_CS%biharmonic) &
    call CS%hv_metrics_CS%Idx2dyCu%alloc(lb=LBOUND(CS%Idx2dyCu), ub=UBOUND(CS%Idx2dyCu), &
                                         source=CS%Idx2dyCu)
  if (CS%hv_opts_CS%biharmonic) &
    call CS%hv_metrics_CS%Idxdy2u%alloc(lb=LBOUND(CS%Idxdy2u), ub=UBOUND(CS%Idxdy2u), &
                                        source=CS%Idxdy2u)
  if (CS%hv_opts_CS%biharmonic) &
    call CS%hv_metrics_CS%Idx2dyCv%alloc(lb=LBOUND(CS%Idx2dyCv), ub=UBOUND(CS%Idx2dyCv), &
                                         source=CS%Idx2dyCv)
  if (CS%hv_opts_CS%biharmonic) &
    call CS%hv_metrics_CS%Idxdy2v%alloc(lb=LBOUND(CS%Idxdy2v), ub=UBOUND(CS%Idxdy2v), &
                                        source=CS%Idxdy2v)
  if (allocated(CS%Laplac2_const_xx)) &
    call CS%hv_metrics_CS%Laplac2_const_xx%alloc(lb=LBOUND(CS%Laplac2_const_xx), ub=UBOUND(CS%Laplac2_const_xx), &
                                                 source=CS%Laplac2_const_xx)
  if (allocated(CS%Biharm6_const_xx)) &
    call CS%hv_metrics_CS%Biharm6_const_xx%alloc(lb=LBOUND(CS%Biharm6_const_xx), ub=UBOUND(CS%Biharm6_const_xx), &
                                                 source=CS%Biharm6_const_xx)
  if (allocated(CS%Laplac3_const_xx)) &
    call CS%hv_metrics_CS%Laplac3_const_xx%alloc(lb=LBOUND(CS%Laplac3_const_xx), ub=UBOUND(CS%Laplac3_const_xx), &
                                                 source=CS%Laplac3_const_xx)
  if (allocated(CS%Biharm_const_xx)) &
    call CS%hv_metrics_CS%Biharm_const_xx%alloc(lb=LBOUND(CS%Biharm_const_xx), ub=UBOUND(CS%Biharm_const_xx), &
                                                source=CS%Biharm_const_xx)
  if (allocated(CS%Biharm_const2_xx)) &
    call CS%hv_metrics_CS%Biharm_const2_xx%alloc(lb=LBOUND(CS%Biharm_const2_xx), ub=UBOUND(CS%Biharm_const2_xx), &
                                                 source=CS%Biharm_const2_xx)
  if (allocated(CS%Re_Ah_const_xx)) &
    call CS%hv_metrics_CS%Re_Ah_const_xx%alloc(lb=LBOUND(CS%Re_Ah_const_xx), ub=UBOUND(CS%Re_Ah_const_xx), &
                                               source=CS%Re_Ah_const_xx)
  if (allocated(CS%Laplac2_const_xy)) &
    call CS%hv_metrics_CS%Laplac2_const_xy%alloc(lb=LBOUND(CS%Laplac2_const_xy), ub=UBOUND(CS%Laplac2_const_xy), &
                                                 source=CS%Laplac2_const_xy)
  if (allocated(CS%Biharm6_const_xy)) &
    call CS%hv_metrics_CS%Biharm6_const_xy%alloc(lb=LBOUND(CS%Biharm6_const_xy), ub=UBOUND(CS%Biharm6_const_xy), &
                                                 source=CS%Biharm6_const_xy)
  if (allocated(CS%Laplac3_const_xy)) &
    call CS%hv_metrics_CS%Laplac3_const_xy%alloc(lb=LBOUND(CS%Laplac3_const_xy), ub=UBOUND(CS%Laplac3_const_xy), &
                                                 source=CS%Laplac3_const_xy)
  if (allocated(CS%Biharm_const_xy)) &
    call CS%hv_metrics_CS%Biharm_const_xy%alloc(lb=LBOUND(CS%Biharm_const_xy), ub=UBOUND(CS%Biharm_const_xy), &
                                                source=CS%Biharm_const_xy)
  if (allocated(CS%Biharm_const2_xy)) &
    call CS%hv_metrics_CS%Biharm_const2_xy%alloc(lb=LBOUND(CS%Biharm_const2_xy), ub=UBOUND(CS%Biharm_const2_xy), &
                                                 source=CS%Biharm_const2_xy)
  if (allocated(CS%Re_Ah_const_xy)) &
    call CS%hv_metrics_CS%Re_Ah_const_xy%alloc(lb=LBOUND(CS%Re_Ah_const_xy), ub=UBOUND(CS%Re_Ah_const_xy), &
                                               source=CS%Re_Ah_const_xy)

  ! Diagnostics that gate calculations within the horizontal viscosity kernels.
  CS%hv_opts_CS%do_grid_Re_Ah = (CS%id_grid_Re_Ah > 0)
  CS%hv_opts_CS%do_grid_Re_Kh = (CS%id_grid_Re_Kh > 0)
  CS%hv_opts_CS%do_Ah_h = (CS%id_Ah_h > 0)
  CS%hv_opts_CS%do_Ah_q = (CS%id_Ah_q > 0)
  CS%hv_opts_CS%do_Kh_h = (CS%id_Kh_h > 0)
  CS%hv_opts_CS%do_Kh_q = (CS%id_Kh_q > 0)
  CS%hv_opts_CS%do_GME_coeff_h = (CS%id_GME_coeff_h > 0)
  CS%hv_opts_CS%do_GME_coeff_q = (CS%id_GME_coeff_q > 0)
  CS%hv_opts_CS%do_vort_xy_q = (CS%id_vort_xy_q > 0)
  CS%hv_opts_CS%do_div_xx_h = (CS%id_div_xx_h > 0)
  CS%hv_opts_CS%do_sh_xy_q = (CS%id_sh_xy_q > 0)
  CS%hv_opts_CS%do_sh_xx_h = (CS%id_sh_xx_h > 0)
  CS%hv_opts_CS%do_FrictWork = (CS%id_FrictWork > 0)
  CS%hv_opts_CS%do_FrictWorkIntz = (CS%id_FrictWorkIntz > 0)
  CS%hv_opts_CS%do_FrictWork_bh = (CS%id_FrictWork_bh > 0)
  CS%hv_opts_CS%do_FrictWorkIntz_bh = (CS%id_FrictWorkIntz_bh > 0)
  CS%hv_opts_CS%do_normstress = (CS%id_normstress > 0)
  CS%hv_opts_CS%do_shearstress = (CS%id_shearstress > 0)
  CS%hv_opts_CS%do_BS_coeff_h = (CS%id_BS_coeff_h > 0)
  CS%hv_opts_CS%do_BS_coeff_q = (CS%id_BS_coeff_q > 0)

end subroutine hor_visc_init

!> hor_visc_vel_stencil returns the horizontal viscosity input velocity stencil size
function hor_visc_vel_stencil(CS) result(stencil)
  type(hor_visc_CS), intent(in) :: CS !< Control structure for horizontal viscosity
  integer ::  stencil !< The horizontal viscosity velocity stencil size with the current settings.

  stencil = 2

  if ((CS%hv_opts_CS%Leith_Kh) .or. (CS%hv_opts_CS%Leith_Ah) .or. (CS%hv_opts_CS%use_Leithy)) then
    stencil = 3
  endif
end function hor_visc_vel_stencil

!> Calculates factors in the anisotropic orientation tensor to be align with the grid.
!! With n1=1 and n2=0, this recovers the approach of Large et al, 2001.
subroutine align_aniso_tensor_to_grid(CS, n1, n2)
  type(hor_visc_CS), intent(inout) :: CS !< Control structure for horizontal viscosity
  real,              intent(in) :: n1 !< i-component of direction vector [nondim]
  real,              intent(in) :: n2 !< j-component of direction vector [nondim]
  ! Local variables
  real :: recip_n2_norm ! The inverse of the squared magnitude of n1 and n2 [nondim]
  ! For normalizing n=(n1,n2) in case arguments are not a unit vector
  recip_n2_norm = (n1**2) + (n2**2)
  if (recip_n2_norm > 0.) recip_n2_norm = 1. / recip_n2_norm
  CS%n1n2_h(:,:) = 2. * ( n1 * n2 ) * recip_n2_norm
  CS%n1n2_q(:,:) = 2. * ( n1 * n2 ) * recip_n2_norm
  CS%n1n1_m_n2n2_h(:,:) = ( (n1 * n1) - (n2 * n2) ) * recip_n2_norm
  CS%n1n1_m_n2n2_q(:,:) = ( (n1 * n1) - (n2 * n2) ) * recip_n2_norm
end subroutine align_aniso_tensor_to_grid

!> Apply a 1-1-4-1-1 Laplacian filter one time on GME diffusive flux to reduce any
!! horizontal two-grid-point noise
subroutine smooth_GME(CS, G, Gcore, bxH0, bxQ0, GME_flux_h, GME_flux_q)
  type(hor_visc_CS),                            intent(in)    :: CS        !< Control structure
  type(ocean_grid_type),                        intent(in)    :: G         !< Ocean grid
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),          intent(in)    :: bxH0  !< The h-point box of the whole column,
                                               !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),          intent(in)    :: bxQ0  !< The B-grid-index box of the whole column,
                                               !! [IscB:IecB, JscB:JecB, 1:nz]
  real, dimension(SZI_(G),SZJ_(G)),   optional, intent(inout) :: GME_flux_h!< GME diffusive flux
                                                              !! at h points [L2 T-2 ~> m2 s-2]
  real, dimension(SZIB_(G),SZJB_(G)), optional, intent(inout) :: GME_flux_q!< GME diffusive flux
                                                              !! at q points [L2 T-2 ~> m2 s-2]
  ! local variables
  type(RealArray_t) :: GME_flux_h_a ! The GME diffusive flux at h points [L2 T-2 ~> m2 s-2]
  type(RealArray_t) :: GME_flux_q_a ! The GME diffusive flux at q points [L2 T-2 ~> m2 s-2]
  integer :: s, halosz
  integer :: xh, xq  ! The number of valid extra halo points for h and q points.

  xh = 0 ; xq = 0

  do s=1,CS%num_smooth_gme
    if (present(GME_flux_h)) then
      if (xh < 0) then
        ! Update halos if needed, but avoid doing so more often than is needed.
        halosz = min(G%isc-G%isd, G%jsc-G%jsd, 2+CS%num_smooth_gme-s)
        call pass_var(GME_flux_h, G%Domain, halo=halosz)
        xh = halosz - 2
      endif
      call GME_flux_h_a%alloc(lb=LBOUND(GME_flux_h), ub=UBOUND(GME_flux_h), &
                              source=GME_flux_h)
      call smooth_GME_h_pass(Gcore, bxH0, bxQ0, GME_flux_h_a, xh)
      call GME_flux_h_a%copy2F(GME_flux_h) ; call GME_flux_h_a%free()
      xh = xh - 1
    endif
    if (present(GME_flux_q)) then
      if (xq < 0) then
        ! Update halos if needed, but avoid doing so more often than is needed.
        halosz = min(G%isc-G%isd, G%jsc-G%jsd, 2+CS%num_smooth_gme-s)
        call pass_var(GME_flux_q, G%Domain, position=CORNER, complete=.true., halo=halosz)
        xq = halosz - 2
      endif
      call GME_flux_q_a%alloc(lb=LBOUND(GME_flux_q), ub=UBOUND(GME_flux_q), &
                              source=GME_flux_q)
      call smooth_GME_q_pass(Gcore, bxH0, bxQ0, GME_flux_q_a, xq)
      call GME_flux_q_a%copy2F(GME_flux_q) ; call GME_flux_q_a%free()
      xq = xq - 1
    endif
  enddo ! s-loop
end subroutine smooth_GME

!> Applies one pass of the 1-1-4-1-1 Laplacian filter to a GME diffusive flux at h points,
!! over the computational domain extended by xh points.
subroutine smooth_GME_h_pass(Gcore, bxH0, bxQ0, GME_flux_h_a, xh)
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),          intent(in)    :: bxH0  !< The h-point box of the whole column,
                                               !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),          intent(in)    :: bxQ0  !< The B-grid-index box of the whole column,
                                               !! [IscB:IecB, JscB:JecB, 1:nz]
  type(RealArray_t), intent(inout) :: GME_flux_h_a !< GME diffusive flux
                                                   !! at h points [L2 T-2 ~> m2 s-2]
  integer,                                      intent(in)    :: xh  !< The number of valid extra halo
                                                              !! points for h points.
  ! local variables
  real, dimension(:,:), contiguous, pointer :: GME_flux_h
  real, dimension(:,:), contiguous, pointer :: mask2dT
  real, dimension(GME_flux_h_a%lb(1):GME_flux_h_a%ub(1),GME_flux_h_a%lb(2):GME_flux_h_a%ub(2)) :: &
    GME_flux_h_original ! The previous value of GME_flux_h [L2 T-2 ~> m2 s-2]
  real :: wc, ww, we, wn, ws ! averaging weights for smoothing [nondim]
  integer :: i, j
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq

  call GME_flux_h_a%view(GME_flux_h)
  call Gcore%mask2dT%view(mask2dT)

  is  = bxH0%idxS(1)  ; ie  = bxH0%idxE(1)  ; js  = bxH0%idxS(2)  ; je  = bxH0%idxE(2)
  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2)

    GME_flux_h_original(:,:) = GME_flux_h(:,:)
    ! apply smoothing on GME
    do j=Jsq-xh,Jeq+1+xh ; do i=Isq-xh,Ieq+1+xh
      ! skip land points
      if (mask2dT(i,j)==0.) cycle
      ! compute weights
      ww = 0.125 * mask2dT(i-1,j)
      we = 0.125 * mask2dT(i+1,j)
      ws = 0.125 * mask2dT(i,j-1)
      wn = 0.125 * mask2dT(i,j+1)
      wc = 1.0 - ((ww+we)+(wn+ws))
      GME_flux_h(i,j) =  wc * GME_flux_h_original(i,j)   &
                       + ((ww * GME_flux_h_original(i-1,j) + we * GME_flux_h_original(i+1,j)) &
                        + (ws * GME_flux_h_original(i,j-1) + wn * GME_flux_h_original(i,j+1)))
    enddo ; enddo

end subroutine smooth_GME_h_pass

!> Applies one pass of the 1-1-4-1-1 Laplacian filter to a GME diffusive flux at q points,
!! over the computational domain extended by xq points.
subroutine smooth_GME_q_pass(Gcore, bxH0, bxQ0, GME_flux_q_a, xq)
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),          intent(in)    :: bxH0  !< The h-point box of the whole column,
                                               !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),          intent(in)    :: bxQ0  !< The B-grid-index box of the whole column,
                                               !! [IscB:IecB, JscB:JecB, 1:nz]
  type(RealArray_t), intent(inout) :: GME_flux_q_a !< GME diffusive flux
                                                   !! at q points [L2 T-2 ~> m2 s-2]
  integer,                                      intent(in)    :: xq  !< The number of valid extra halo
                                                              !! points for q points.
  ! local variables
  real, dimension(:,:), contiguous, pointer :: GME_flux_q
  real, dimension(:,:), contiguous, pointer :: mask2dBu
  real, dimension(GME_flux_q_a%lb(1):GME_flux_q_a%ub(1),GME_flux_q_a%lb(2):GME_flux_q_a%ub(2)) :: &
    GME_flux_q_original ! The previous value of GME_flux_q [L2 T-2 ~> m2 s-2]
  real :: wc, ww, we, wn, ws ! averaging weights for smoothing [nondim]
  integer :: i, j
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq

  call GME_flux_q_a%view(GME_flux_q)
  call Gcore%mask2dBu%view(mask2dBu)

  is  = bxH0%idxS(1)  ; ie  = bxH0%idxE(1)  ; js  = bxH0%idxS(2)  ; je  = bxH0%idxE(2)
  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2)

    GME_flux_q_original(:,:) = GME_flux_q(:,:)
    ! apply smoothing on GME
    do J=js-1-xq,je+xq ; do I=is-1-xq,ie+xq
      ! skip land points
      if (mask2dBu(I,J)==0.) cycle
      ! compute weights
      ww = 0.125 * mask2dBu(I-1,J)
      we = 0.125 * mask2dBu(I+1,J)
      ws = 0.125 * mask2dBu(I,J-1)
      wn = 0.125 * mask2dBu(I,J+1)
      wc = 1.0 - ((ww+we)+(wn+ws))
      GME_flux_q(I,J) =  wc * GME_flux_q_original(I,J)   &
                       + ((ww * GME_flux_q_original(I-1,J) + we * GME_flux_q_original(I+1,J)) &
                        + (ws * GME_flux_q_original(I,J-1) + wn * GME_flux_q_original(I,J+1)))
    enddo ; enddo

end subroutine smooth_GME_q_pass

!> Smooths the GME stress terms for each layer of one tile.
!! This includes halo updates, so it always runs in Fortran.
subroutine hor_visc_GME_smooth_host(host, Gcore, bxH0, bxQ0, nkblock, kmax, str_xx_GME_a, &
                                    str_xy_GME_a)
  type(c_ptr), value,   intent(in)    :: host  !< Opaque pointer to the hv_host_ctx that holds the Fortran
                                               !! objects used by the host helpers
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),          intent(in)    :: bxH0  !< The h-point box of the whole column,
                                               !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),          intent(in)    :: bxQ0  !< The B-grid-index box of the whole column,
                                               !! [IscB:IecB, JscB:JecB, 1:nz]
  integer,                       intent(in)    :: nkblock !< The effective k-block size [nondim]
  integer,                       intent(in)    :: kmax   !< The number of active k-layers in the current k-block
  type(RealArray_t), intent(inout) :: str_xx_GME_a !< smoothed diagonal term in the stress tensor from GME
                                                   !! [L2 T-2 ~> m2 s-2]
  type(RealArray_t), intent(inout) :: str_xy_GME_a !< smoothed cross term in the stress tensor from GME
                                                   !! [L2 T-2 ~> m2 s-2]
  ! Local variables
  type(hv_host_ctx),       pointer :: ctx    ! The objects behind host
  type(ocean_grid_type),   pointer :: G      ! The ocean's grid structure.
  type(hor_visc_CS),       pointer :: CS     ! Horizontal viscosity control structure
  real, dimension(:,:,:), contiguous, pointer :: str_xy_GME, str_xx_GME
  integer :: kk

  call c_f_pointer(host, ctx)
  G => ctx%G ; CS => ctx%CS

  call str_xy_GME_a%view(str_xy_GME) ; call str_xx_GME_a%view(str_xx_GME)

    do kk=1,kmax
      call smooth_GME(CS, G, Gcore, bxH0, bxQ0, GME_flux_h=str_xx_GME(:,:,kk))
      call smooth_GME(CS, G, Gcore, bxH0, bxQ0, GME_flux_q=str_xy_GME(:,:,kk))
    enddo

end subroutine hor_visc_GME_smooth_host

!> Passes the velocity gradients and thickness of each layer of one tile to ZB2020.
!! This calls into MOM_Zanna_Bolton, so it always runs in Fortran.
subroutine hor_visc_ZB2020_copy_host(host, nkblock, kstart, kmax, sh_xx_a, sh_xy_a, vort_xy_a, hq_a)
  type(c_ptr), value,   intent(in)    :: host  !< Opaque pointer to the hv_host_ctx that holds the Fortran
                                               !! objects used by the host helpers
  integer,                       intent(in)    :: nkblock !< The effective k-block size [nondim]
  integer,                       intent(in)    :: kstart !< The first absolute k-layer of the current k-block
  integer,                       intent(in)    :: kmax   !< The number of active k-layers in the current k-block
  type(RealArray_t), intent(in)    :: sh_xx_a !< horizontal tension (du/dx - dv/dy) including metric terms
                                              !! [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: sh_xy_a !< horizontal shearing strain (du/dy + dv/dx) including metric
                                              !! terms [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: vort_xy_a !< Vertical vorticity (dv/dx - du/dy) including metric terms
                                                !! [T-1 ~> s-1]
  type(RealArray_t), intent(in)    :: hq_a !< harmonic mean of the harmonic means of the u- & v point
                                           !! thicknesses [H ~> m or kg m-2] This form guarantees that hq/hu <
                                           !! 4.
  ! Local variables
  type(hv_host_ctx),       pointer :: ctx    ! The objects behind host
  type(ocean_grid_type),   pointer :: G      ! The ocean's grid structure.
  type(verticalGrid_type), pointer :: GV     ! The ocean's vertical grid structure.
  type(hor_visc_CS),       pointer :: CS     ! Horizontal viscosity control structure
  real, dimension(:,:,:), contiguous, pointer :: hq, vort_xy, sh_xy, sh_xx
  integer :: k, kk

  call c_f_pointer(host, ctx)
  G => ctx%G ; GV => ctx%GV ; CS => ctx%CS

  call hq_a%view(hq) ; call vort_xy_a%view(vort_xy) ; call sh_xy_a%view(sh_xy)
  call sh_xx_a%view(sh_xx)

    do kk=1,kmax
      k = kstart + kk - 1
      call ZB2020_copy_gradient_and_thickness(sh_xx(:,:,kk), sh_xy(:,:,kk), vort_xy(:,:,kk), hq(:,:,kk), &
                                               G, GV, CS%ZB2020, k)
    enddo

end subroutine hor_visc_ZB2020_copy_host

!> Apply a 9-point smoothing filter twice to a field staggered at a thickness point to reduce
!! horizontal two-grid-point noise.
!! Note that this subroutine does not conserve mass, so don't use it in situations where you
!! need conservation.  Also note that it assumes that the input field has valid values in the
!! first two halo points upon entry.
subroutine smooth_x9_h(Gcore, bxH0, field_h_a, zero_land)
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),          intent(in)    :: bxH0  !< The h-point box of the whole column,
                                               !! [isc:iec, jsc:jec, 1:nz]
  type(RealArray_t), intent(inout) :: field_h_a !< h-point field to be smoothed [arbitrary]
  logical,                intent(in)              :: zero_land !< If present and false, return the average
                                                               !! of the surrounding ocean points when
                                                               !! smoothing, otherwise use a value of 0 for
                                                               !! land points and include them in the averages.

  ! Local variables
  real, dimension(:,:), contiguous, pointer :: field_h
  real, dimension(:,:), contiguous, pointer :: mask2dT
  real, dimension(field_h_a%lb(1):field_h_a%ub(1),field_h_a%lb(2):field_h_a%ub(2)) :: &
    fh_prev  ! The value of the h-point field at the previous iteration [arbitrary]
  real :: Iwts             ! The inverse of the sum of the weights [nondim]
  logical :: zero_land_val ! The value of the zero_land optional argument or .true. if it is absent.
  integer :: i, j, s, is, ie, js, je

  call field_h_a%view(field_h)
  call Gcore%mask2dT%view(mask2dT)

  is  = bxH0%idxS(1)  ; ie  = bxH0%idxE(1)  ; js  = bxH0%idxS(2)  ; je  = bxH0%idxE(2)

  zero_land_val = zero_land

  do s=1,0,-1
    fh_prev(:,:) = field_h(:,:)
    ! apply smoothing on field_h using rotationally symmetric expressions.
    do j=js-s,je+s ; do i=is-s,ie+s ; if (mask2dT(i,j) > 0.0) then
      Iwts = 0.0625
      if (.not. zero_land_val) &
        Iwts = 1.0 / ( (4.0*mask2dT(i,j) + &
                        ( 2.0*((mask2dT(i-1,j) + mask2dT(i+1,j)) + &
                               (mask2dT(i,j-1) + mask2dT(i,j+1))) + &
                         ((mask2dT(i-1,j-1) + mask2dT(i+1,j+1)) + &
                          (mask2dT(i-1,j+1) + mask2dT(i+1,j-1))) ) ) + 1.0e-16 )
      field_h(i,j) = Iwts * ( 4.0*mask2dT(i,j) * fh_prev(i,j) &
                            + (2.0*((mask2dT(i-1,j) * fh_prev(i-1,j) + mask2dT(i+1,j) * fh_prev(i+1,j)) + &
                                    (mask2dT(i,j-1) * fh_prev(i,j-1) + mask2dT(i,j+1) * fh_prev(i,j+1))) &
                              + ((mask2dT(i-1,j-1) * fh_prev(i-1,j-1) + mask2dT(i+1,j+1) * fh_prev(i+1,j+1)) + &
                                 (mask2dT(i-1,j+1) * fh_prev(i-1,j+1) + mask2dT(i+1,j-1) * fh_prev(i-1,j-1))) ))
    endif ; enddo ; enddo
  enddo

end subroutine smooth_x9_h

!> Apply a 9-point smoothing filter twice to a pair of velocity components to reduce
!! horizontal two-grid-point noise.
!! Note that this subroutine does not conserve angular momentum, so don't use it
!! in situations where you need conservation.  Also note that it assumes that the
!! input fields have valid values in the first two halo points upon entry.
subroutine smooth_x9_uv(Gcore, bxH0, bxQ0, field_u_a, field_v_a, zero_land)
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),          intent(in)    :: bxH0  !< The h-point box of the whole column,
                                               !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),          intent(in)    :: bxQ0  !< The B-grid-index box of the whole column,
                                               !! [IscB:IecB, JscB:JecB, 1:nz]
  type(RealArray_t), intent(inout) :: field_u_a !< u-point field to be smoothed [arbitrary]
  type(RealArray_t), intent(inout) :: field_v_a !< v-point field to be smoothed [arbitrary]
  logical,                 intent(in)              :: zero_land !< If present and false, return the average
                                                                !! of the surrounding ocean points when
                                                                !! smoothing, otherwise use a value of 0 for
                                                                !! land points and include them in the averages.

  ! Local variables.
  real, dimension(:,:), contiguous, pointer :: field_v, field_u
  real, dimension(:,:), contiguous, pointer :: mask2dCu, mask2dCv
  real, dimension(field_u_a%lb(1):field_u_a%ub(1),field_u_a%lb(2):field_u_a%ub(2)) :: &
    fu_prev  ! The value of the u-point field at the previous iteration [arbitrary]
  real, dimension(field_v_a%lb(1):field_v_a%ub(1),field_v_a%lb(2):field_v_a%ub(2)) :: &
    fv_prev  ! The value of the v-point field at the previous iteration [arbitrary]
  real :: Iwts             ! The inverse of the sum of the weights [nondim]
  logical :: zero_land_val ! The value of the zero_land optional argument or .true. if it is absent.
  integer :: i, j, s, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  call field_v_a%view(field_v) ; call field_u_a%view(field_u)
  call Gcore%mask2dCu%view(mask2dCu) ; call Gcore%mask2dCv%view(mask2dCv)

  is  = bxH0%idxS(1)  ; ie  = bxH0%idxE(1)  ; js  = bxH0%idxS(2)  ; je  = bxH0%idxE(2)
  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2)

  zero_land_val = zero_land

  do s=1,0,-1
    fu_prev(:,:) = field_u(:,:)
    ! apply smoothing on field_u using the original non-rotationally symmetric expressions.
    do j=js-s,je+s ; do I=Isq-s,Ieq+s ; if (mask2dCu(I,j) > 0.0) then
      Iwts = 0.0625
      if (.not. zero_land_val) &
        Iwts = 1.0 / ( (4.0*mask2dCu(I,j) + &
                        ( 2.0*((mask2dCu(I-1,j) + mask2dCu(I+1,j)) + &
                               (mask2dCu(I,j-1) + mask2dCu(I,j+1))) + &
                         ((mask2dCu(I-1,j-1) + mask2dCu(I+1,j+1)) + &
                          (mask2dCu(I-1,j+1) + mask2dCu(I+1,j-1))) ) ) + 1.0e-16 )
      field_u(I,j) = Iwts * ( 4.0*mask2dCu(I,j) * fu_prev(I,j) &
                            + (2.0*((mask2dCu(I-1,j) * fu_prev(I-1,j) + mask2dCu(I+1,j) * fu_prev(I+1,j)) + &
                                    (mask2dCu(I,j-1) * fu_prev(I,j-1) + mask2dCu(I,j+1) * fu_prev(I,j+1))) &
                              + ((mask2dCu(I-1,j-1) * fu_prev(I-1,j-1) + mask2dCu(I+1,j+1) * fu_prev(I+1,j+1)) + &
                                 (mask2dCu(I-1,j+1) * fu_prev(I-1,j+1) + mask2dCu(I+1,j-1) * fu_prev(I-1,j-1))) ))
    endif ; enddo ; enddo

    fv_prev(:,:) = field_v(:,:)
    ! apply smoothing on field_v using the original non-rotationally symmetric expressions.
    do J=Jsq-s,Jeq+s ; do i=is-s,ie+s ; if (mask2dCv(i,J) > 0.0) then
      Iwts = 0.0625
      if (.not. zero_land_val) &
        Iwts = 1.0 / ( (4.0*mask2dCv(i,J) + &
                        ( 2.0*((mask2dCv(i-1,J) + mask2dCv(i+1,J)) + &
                               (mask2dCv(i,J-1) + mask2dCv(i,J+1))) + &
                         ((mask2dCv(i-1,J-1) + mask2dCv(i+1,J+1)) + &
                          (mask2dCv(i-1,J+1) + mask2dCv(i+1,J-1))) ) ) + 1.0e-16 )
      field_v(i,J) = Iwts * ( 4.0*mask2dCv(i,J) * fv_prev(i,J) &
                            + (2.0*((mask2dCv(i-1,J) * fv_prev(i-1,J) + mask2dCv(i+1,J) * fv_prev(i+1,J)) + &
                                    (mask2dCv(i,J-1) * fv_prev(i,J-1) + mask2dCv(i,J+1) * fv_prev(i,J+1))) &
                              + ((mask2dCv(i-1,J-1) * fv_prev(i-1,J-1) + mask2dCv(i+1,J+1) * fv_prev(i+1,J+1)) + &
                                 (mask2dCv(i-1,J+1) * fv_prev(i-1,J+1) + mask2dCv(i+1,J-1) * fv_prev(i-1,J-1))) ))
    endif ; enddo ; enddo
  enddo

end subroutine smooth_x9_uv

!> Deallocates any variables allocated in hor_visc_init.
subroutine hor_visc_end(CS)
  type(hor_visc_CS), intent(inout) :: CS !< Horizontal viscosity control structure

  if (CS%hv_opts_CS%Laplacian .or. CS%hv_opts_CS%biharmonic) then
    DEALLOC_(CS%dx2h) ; DEALLOC_(CS%dx2q) ; DEALLOC_(CS%dy2h) ; DEALLOC_(CS%dy2q)
    DEALLOC_(CS%dx_dyT) ; DEALLOC_(CS%dy_dxT) ; DEALLOC_(CS%dx_dyBu) ; DEALLOC_(CS%dy_dxBu)
    DEALLOC_(CS%reduction_xx) ; DEALLOC_(CS%reduction_xy)
  endif
  if (CS%hv_opts_CS%Laplacian) then
    DEALLOC_(CS%Kh_bg_xx) ; DEALLOC_(CS%Kh_bg_xy)
    DEALLOC_(CS%grid_sp_h2)
    if (allocated(CS%Kh_bg_2d)) deallocate(CS%Kh_bg_2d)

    if (allocated(CS%Kh_Max_xx)) deallocate(CS%Kh_Max_xx)
    if (allocated(CS%Kh_Max_xy)) deallocate(CS%Kh_Max_xy)
    if (allocated(CS%Laplac2_const_xx)) deallocate(CS%Laplac2_const_xx)
    if (allocated(CS%Laplac2_const_xy)) deallocate(CS%Laplac2_const_xy)
    if (allocated(CS%Laplac3_const_xx)) deallocate(CS%Laplac3_const_xx)
    if (allocated(CS%Laplac3_const_xy)) deallocate(CS%Laplac3_const_xy)
  endif
  if (CS%hv_opts_CS%biharmonic) then
    DEALLOC_(CS%grid_sp_h3)
    DEALLOC_(CS%Idx2dyCu) ; DEALLOC_(CS%Idx2dyCv)
    DEALLOC_(CS%Idxdy2u) ; DEALLOC_(CS%Idxdy2v)
    DEALLOC_(CS%Ah_bg_xx) ; DEALLOC_(CS%Ah_bg_xy)

    if (allocated(CS%Ah_Max_xx)) deallocate(CS%Ah_Max_xx)
    if (allocated(CS%Ah_Max_xy)) deallocate(CS%Ah_Max_xy)
    if (allocated(CS%Ah_Max_xx_KS)) deallocate(CS%Ah_Max_xx_KS)
    if (allocated(CS%Ah_Max_xy_KS)) deallocate(CS%Ah_Max_xy_KS)
    if (allocated(CS%Biharm_const_xx)) deallocate(CS%Biharm_const_xx)
    if (allocated(CS%Biharm_const_xy)) deallocate(CS%Biharm_const_xy)
    if (allocated(CS%Biharm_const2_xx)) deallocate(CS%Biharm_const2_xx)
    if (allocated(CS%Biharm_const2_xy)) deallocate(CS%Biharm_const2_xy)
    if (allocated(CS%Biharm6_const_xx)) deallocate(CS%Biharm6_const_xx)
    if (allocated(CS%Biharm6_const_xy)) deallocate(CS%Biharm6_const_xy)
    if (allocated(CS%m_const_leithy)) deallocate(CS%m_const_leithy)
    if (allocated(CS%m_leithy_max)) deallocate(CS%m_leithy_max)
    if (allocated(CS%Re_Ah_const_xx)) deallocate(CS%Re_Ah_const_xx)
    if (allocated(CS%Re_Ah_const_xy)) deallocate(CS%Re_Ah_const_xy)
  endif

  if (allocated(CS%n1n2_h)) deallocate(CS%n1n2_h)
  if (allocated(CS%n1n2_q)) deallocate(CS%n1n2_q)
  if (allocated(CS%n1n1_m_n2n2_h)) deallocate(CS%n1n1_m_n2n2_h)
  if (allocated(CS%n1n1_m_n2n2_q)) deallocate(CS%n1n1_m_n2n2_q)

  if (CS%hv_opts_CS%use_ZB2020) then
    call ZB2020_end(CS%ZB2020)
  endif


  call CS%hv_metrics_CS%Kh_bg_xx%free() ; call CS%hv_metrics_CS%Ah_bg_xx%free()
  call CS%hv_metrics_CS%reduction_xx%free() ; call CS%hv_metrics_CS%Kh_Max_xx%free()
  call CS%hv_metrics_CS%Ah_Max_xx%free() ; call CS%hv_metrics_CS%Ah_Max_xx_KS%free()
  call CS%hv_metrics_CS%n1n2_h%free() ; call CS%hv_metrics_CS%n1n1_m_n2n2_h%free()
  call CS%hv_metrics_CS%grid_sp_h2%free() ; call CS%hv_metrics_CS%grid_sp_h3%free()
  call CS%hv_metrics_CS%Kh_bg_xy%free() ; call CS%hv_metrics_CS%Ah_bg_xy%free()
  call CS%hv_metrics_CS%reduction_xy%free() ; call CS%hv_metrics_CS%Kh_Max_xy%free()
  call CS%hv_metrics_CS%Ah_Max_xy%free() ; call CS%hv_metrics_CS%Ah_Max_xy_KS%free()
  call CS%hv_metrics_CS%n1n2_q%free() ; call CS%hv_metrics_CS%n1n1_m_n2n2_q%free()
  call CS%hv_metrics_CS%dx2h%free() ; call CS%hv_metrics_CS%dy2h%free()
  call CS%hv_metrics_CS%dx_dyT%free() ; call CS%hv_metrics_CS%dy_dxT%free()
  call CS%hv_metrics_CS%m_const_leithy%free() ; call CS%hv_metrics_CS%m_leithy_max%free()
  call CS%hv_metrics_CS%dx2q%free() ; call CS%hv_metrics_CS%dy2q%free()
  call CS%hv_metrics_CS%dx_dyBu%free() ; call CS%hv_metrics_CS%dy_dxBu%free()
  call CS%hv_metrics_CS%Idx2dyCu%free() ; call CS%hv_metrics_CS%Idxdy2u%free()
  call CS%hv_metrics_CS%Idx2dyCv%free() ; call CS%hv_metrics_CS%Idxdy2v%free()
  call CS%hv_metrics_CS%Laplac2_const_xx%free() ; call CS%hv_metrics_CS%Biharm6_const_xx%free()
  call CS%hv_metrics_CS%Laplac3_const_xx%free() ; call CS%hv_metrics_CS%Biharm_const_xx%free()
  call CS%hv_metrics_CS%Biharm_const2_xx%free() ; call CS%hv_metrics_CS%Re_Ah_const_xx%free()
  call CS%hv_metrics_CS%Laplac2_const_xy%free() ; call CS%hv_metrics_CS%Biharm6_const_xy%free()
  call CS%hv_metrics_CS%Laplac3_const_xy%free() ; call CS%hv_metrics_CS%Biharm_const_xy%free()
  call CS%hv_metrics_CS%Biharm_const2_xy%free() ; call CS%hv_metrics_CS%Re_Ah_const_xy%free()

end subroutine hor_visc_end

!> Allocates a container for, and copies in, each MEKE field used by horizontal_viscosity that is
!! allocated, along with the backscatter parameters.
subroutine hor_visc_MEKE_build_from(this, MEKE)
  class(hor_visc_MEKE_type), intent(inout) :: this !< Copies of the MEKE fields used by horizontal_viscosity
  type(MEKE_type),           intent(in)    :: MEKE !< MEKE fields related to Mesoscale Eddy Kinetic Energy.

  if (allocated(MEKE%Ku)) call this%Ku_a%alloc(lb=LBOUND(MEKE%Ku), ub=UBOUND(MEKE%Ku), source=MEKE%Ku)
  if (allocated(MEKE%Au)) call this%Au_a%alloc(lb=LBOUND(MEKE%Au), ub=UBOUND(MEKE%Au), source=MEKE%Au)
  if (allocated(MEKE%mom_src)) &
    call this%mom_src_a%alloc(lb=LBOUND(MEKE%mom_src), ub=UBOUND(MEKE%mom_src), source=MEKE%mom_src)
  if (allocated(MEKE%mom_src_bh)) &
    call this%mom_src_bh_a%alloc(lb=LBOUND(MEKE%mom_src_bh), ub=UBOUND(MEKE%mom_src_bh), &
                                 source=MEKE%mom_src_bh)
  if (allocated(MEKE%GME_snk)) &
    call this%GME_snk_a%alloc(lb=LBOUND(MEKE%GME_snk), ub=UBOUND(MEKE%GME_snk), source=MEKE%GME_snk)
  this%backscatter_Ro_pow = MEKE%backscatter_Ro_pow
  this%backscatter_Ro_c = MEKE%backscatter_Ro_c
end subroutine hor_visc_MEKE_build_from

!> Copies the MEKE source and sink fields written by horizontal_viscosity back into MEKE, then frees
!! every container.
subroutine hor_visc_MEKE_copy_back(this, MEKE)
  class(hor_visc_MEKE_type), intent(inout) :: this !< Copies of the MEKE fields used by horizontal_viscosity
  type(MEKE_type),           intent(inout) :: MEKE !< MEKE fields related to Mesoscale Eddy Kinetic Energy.

  if (this%mom_src_a%associated()) call this%mom_src_a%copy2F(MEKE%mom_src)
  if (this%mom_src_bh_a%associated()) call this%mom_src_bh_a%copy2F(MEKE%mom_src_bh)
  if (this%GME_snk_a%associated()) call this%GME_snk_a%copy2F(MEKE%GME_snk)
  call this%Ku_a%free() ; call this%Au_a%free()
  call this%mom_src_a%free() ; call this%mom_src_bh_a%free() ; call this%GME_snk_a%free()
end subroutine hor_visc_MEKE_copy_back

!> Allocates a container for, and copies in, each VarMix field used by horizontal_viscosity that is
!! allocated, along with the variable mixing switches.
subroutine hor_visc_VarMix_build_from(this, VarMix)
  class(hor_visc_VarMix_type), intent(inout) :: this   !< Copies of the variable mixing fields used by
                                                       !! horizontal_viscosity
  type(VarMix_CS),             intent(in)    :: VarMix !< Variable mixing control structure

  if (allocated(VarMix%Res_fn_h)) &
    call this%Res_fn_h_a%alloc(lb=LBOUND(VarMix%Res_fn_h), ub=UBOUND(VarMix%Res_fn_h), source=VarMix%Res_fn_h)
  if (allocated(VarMix%Res_fn_q)) &
    call this%Res_fn_q_a%alloc(lb=LBOUND(VarMix%Res_fn_q), ub=UBOUND(VarMix%Res_fn_q), source=VarMix%Res_fn_q)
  if (allocated(VarMix%BS_struct)) &
    call this%BS_struct_a%alloc(lb=LBOUND(VarMix%BS_struct), ub=UBOUND(VarMix%BS_struct), &
                                source=VarMix%BS_struct)
  this%use_variable_mixing = VarMix%use_variable_mixing
  this%Resoln_scaled_Kh = VarMix%Resoln_scaled_Kh
end subroutine hor_visc_VarMix_build_from

!> Frees every container of the variable mixing fields used by horizontal_viscosity.
subroutine hor_visc_VarMix_free(this)
  class(hor_visc_VarMix_type), intent(inout) :: this !< Copies of the variable mixing fields used by
                                                     !! horizontal_viscosity

  call this%Res_fn_h_a%free() ; call this%Res_fn_q_a%free() ; call this%BS_struct_a%free()
end subroutine hor_visc_VarMix_free

!> Allocates a container for, and copies in, the stochastic dissipation rate used by
!! horizontal_viscosity if it is allocated, along with the frictional SKEB switch and coefficient.
subroutine hor_visc_STOCH_build_from(this, STOCH)
  class(hor_visc_STOCH_type), intent(inout) :: this  !< Copies of the stochastic fields used by
                                                     !! horizontal_viscosity
  type(stochastic_CS),        intent(in)    :: STOCH !< Stochastic control structure

  if (allocated(STOCH%skeb_diss)) &
    call this%skeb_diss_a%alloc(lb=LBOUND(STOCH%skeb_diss), ub=UBOUND(STOCH%skeb_diss), source=STOCH%skeb_diss)
  this%skeb_use_frict = STOCH%skeb_use_frict
  this%skeb_frict_coef = STOCH%skeb_frict_coef
end subroutine hor_visc_STOCH_build_from

!> Copies the stochastic dissipation rate written by horizontal_viscosity back into STOCH, then frees
!! its container.
subroutine hor_visc_STOCH_copy_back(this, STOCH)
  class(hor_visc_STOCH_type), intent(inout) :: this  !< Copies of the stochastic fields used by
                                                     !! horizontal_viscosity
  type(stochastic_CS),        intent(inout) :: STOCH !< Stochastic control structure

  if (this%skeb_diss_a%associated()) call this%skeb_diss_a%copy2F(STOCH%skeb_diss)
  call this%skeb_diss_a%free()
end subroutine hor_visc_STOCH_copy_back

!> Converts a hor_visc_MEKE_type to its bind(C) mirror.
#ifdef _TIM
function hor_visc_MEKE_to_c(MEKE_a) result(cdesc)
  type(hor_visc_MEKE_type), target, intent(in) :: MEKE_a !< Copies of the MEKE fields
                                                         !! used by horizontal_viscosity
  type(hor_visc_MEKE_C) :: cdesc !< bind(C) mirror of MEKE_a
  cdesc%Ku = MEKE_a%Ku_a%to_c()
  cdesc%Au = MEKE_a%Au_a%to_c()
  cdesc%mom_src = MEKE_a%mom_src_a%to_c()
  cdesc%mom_src_bh = MEKE_a%mom_src_bh_a%to_c()
  cdesc%GME_snk = MEKE_a%GME_snk_a%to_c()
  cdesc%backscatter_Ro_pow = MEKE_a%backscatter_Ro_pow
  cdesc%backscatter_Ro_c = MEKE_a%backscatter_Ro_c
end function hor_visc_MEKE_to_c
#endif

!> Records the containers and scalars of a hor_visc_MEKE_type in a capture file.
subroutine record_hor_visc_MEKE(rec, prefix, MEKE_a, suffix)
  type(io_recorder), intent(inout) :: rec    !< The capture file being written
  character(len=*),  intent(in)    :: prefix !< The start of each record name
  type(hor_visc_MEKE_type), intent(in) :: MEKE_a !< Copies of the MEKE fields
                                                 !! used by horizontal_viscosity
  character(len=*),  intent(in)    :: suffix !< The end of each value record name, or empty to
                                             !! record which arrays are allocated and the scalars

  if (len(suffix) == 0) then
    call rec%add(prefix//"_Ku_allocated", MEKE_a%Ku_a%associated())
  elseif (MEKE_a%Ku_a%associated()) then
    call rec%add(prefix//"_Ku"//suffix, MEKE_a%Ku_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_Au_allocated", MEKE_a%Au_a%associated())
  elseif (MEKE_a%Au_a%associated()) then
    call rec%add(prefix//"_Au"//suffix, MEKE_a%Au_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_mom_src_allocated", MEKE_a%mom_src_a%associated())
  elseif (MEKE_a%mom_src_a%associated()) then
    call rec%add(prefix//"_mom_src"//suffix, MEKE_a%mom_src_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_mom_src_bh_allocated", MEKE_a%mom_src_bh_a%associated())
  elseif (MEKE_a%mom_src_bh_a%associated()) then
    call rec%add(prefix//"_mom_src_bh"//suffix, MEKE_a%mom_src_bh_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_GME_snk_allocated", MEKE_a%GME_snk_a%associated())
  elseif (MEKE_a%GME_snk_a%associated()) then
    call rec%add(prefix//"_GME_snk"//suffix, MEKE_a%GME_snk_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_backscatter_Ro_pow", MEKE_a%backscatter_Ro_pow)
    call rec%add(prefix//"_backscatter_Ro_c", MEKE_a%backscatter_Ro_c)
  endif
end subroutine record_hor_visc_MEKE

!> Converts a hor_visc_VarMix_type to its bind(C) mirror.
#ifdef _TIM
function hor_visc_VarMix_to_c(VarMix_a) result(cdesc)
  type(hor_visc_VarMix_type), target, intent(in) :: VarMix_a !< Copies of the variable mixing fields
                                                             !! used by horizontal_viscosity
  type(hor_visc_VarMix_C) :: cdesc !< bind(C) mirror of VarMix_a
  cdesc%Res_fn_h = VarMix_a%Res_fn_h_a%to_c()
  cdesc%Res_fn_q = VarMix_a%Res_fn_q_a%to_c()
  cdesc%BS_struct = VarMix_a%BS_struct_a%to_c()
  cdesc%use_variable_mixing = VarMix_a%use_variable_mixing
  cdesc%Resoln_scaled_Kh = VarMix_a%Resoln_scaled_Kh
end function hor_visc_VarMix_to_c
#endif

!> Records the containers and scalars of a hor_visc_VarMix_type in a capture file.
subroutine record_hor_visc_VarMix(rec, prefix, VarMix_a, suffix)
  type(io_recorder), intent(inout) :: rec    !< The capture file being written
  character(len=*),  intent(in)    :: prefix !< The start of each record name
  type(hor_visc_VarMix_type), intent(in) :: VarMix_a !< Copies of the variable mixing fields
                                                     !! used by horizontal_viscosity
  character(len=*),  intent(in)    :: suffix !< The end of each value record name, or empty to
                                             !! record which arrays are allocated and the scalars

  if (len(suffix) == 0) then
    call rec%add(prefix//"_Res_fn_h_allocated", VarMix_a%Res_fn_h_a%associated())
  elseif (VarMix_a%Res_fn_h_a%associated()) then
    call rec%add(prefix//"_Res_fn_h"//suffix, VarMix_a%Res_fn_h_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_Res_fn_q_allocated", VarMix_a%Res_fn_q_a%associated())
  elseif (VarMix_a%Res_fn_q_a%associated()) then
    call rec%add(prefix//"_Res_fn_q"//suffix, VarMix_a%Res_fn_q_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_BS_struct_allocated", VarMix_a%BS_struct_a%associated())
  elseif (VarMix_a%BS_struct_a%associated()) then
    call rec%add(prefix//"_BS_struct"//suffix, VarMix_a%BS_struct_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_use_variable_mixing", VarMix_a%use_variable_mixing)
    call rec%add(prefix//"_Resoln_scaled_Kh", VarMix_a%Resoln_scaled_Kh)
  endif
end subroutine record_hor_visc_VarMix

!> Converts a hor_visc_STOCH_type to its bind(C) mirror.
#ifdef _TIM
function hor_visc_STOCH_to_c(STOCH_a) result(cdesc)
  type(hor_visc_STOCH_type), target, intent(in) :: STOCH_a !< Copies of the stochastic fields
                                                           !! used by horizontal_viscosity
  type(hor_visc_STOCH_C) :: cdesc !< bind(C) mirror of STOCH_a
  cdesc%skeb_diss = STOCH_a%skeb_diss_a%to_c()
  cdesc%skeb_use_frict = STOCH_a%skeb_use_frict
  cdesc%skeb_frict_coef = STOCH_a%skeb_frict_coef
end function hor_visc_STOCH_to_c
#endif

!> Records the containers and scalars of a hor_visc_STOCH_type in a capture file.
subroutine record_hor_visc_STOCH(rec, prefix, STOCH_a, suffix)
  type(io_recorder), intent(inout) :: rec    !< The capture file being written
  character(len=*),  intent(in)    :: prefix !< The start of each record name
  type(hor_visc_STOCH_type), intent(in) :: STOCH_a !< Copies of the stochastic fields
                                                   !! used by horizontal_viscosity
  character(len=*),  intent(in)    :: suffix !< The end of each value record name, or empty to
                                             !! record which arrays are allocated and the scalars

  if (len(suffix) == 0) then
    call rec%add(prefix//"_skeb_diss_allocated", STOCH_a%skeb_diss_a%associated())
  elseif (STOCH_a%skeb_diss_a%associated()) then
    call rec%add(prefix//"_skeb_diss"//suffix, STOCH_a%skeb_diss_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_skeb_use_frict", STOCH_a%skeb_use_frict)
    call rec%add(prefix//"_skeb_frict_coef", STOCH_a%skeb_frict_coef)
  endif
end subroutine record_hor_visc_STOCH

!> Converts a hv_opts_CS to its bind(C) mirror.
function hv_opts_CS_to_c(opts) result(cdesc)
  type(hv_opts_CS), intent(in) :: opts  !< Options to convert
  type(hv_opts_CS_C)           :: cdesc !< bind(C) mirror of opts
  cdesc%Laplacian = opts%Laplacian
  cdesc%biharmonic = opts%biharmonic
  cdesc%debug = opts%debug
  cdesc%no_slip = opts%no_slip
  cdesc%bound_Kh = opts%bound_Kh
  cdesc%EY24_EBT_BS = opts%EY24_EBT_BS
  cdesc%bound_Ah = opts%bound_Ah
  cdesc%Re_Ah = opts%Re_Ah
  cdesc%KS_coef = opts%KS_coef
  cdesc%backscatter_underbound = opts%backscatter_underbound
  cdesc%Smagorinsky_Kh = opts%Smagorinsky_Kh
  cdesc%Smagorinsky_Ah = opts%Smagorinsky_Ah
  cdesc%Leith_Kh = opts%Leith_Kh
  cdesc%Modified_Leith = opts%Modified_Leith
  cdesc%use_beta_in_Leith = opts%use_beta_in_Leith
  cdesc%Leith_Ah = opts%Leith_Ah
  cdesc%use_Leithy = opts%use_Leithy
  cdesc%c_K = opts%c_K
  cdesc%smooth_Ah = opts%smooth_Ah
  cdesc%use_QG_Leith_visc = opts%use_QG_Leith_visc
  cdesc%bound_Coriolis = opts%bound_Coriolis
  cdesc%Kh_bg_min = opts%Kh_bg_min
  cdesc%FrictWork_bug = opts%FrictWork_bug
  cdesc%OBC_strain_bug = opts%OBC_strain_bug
  cdesc%use_land_mask = opts%use_land_mask
  cdesc%anisotropic = opts%anisotropic
  cdesc%add_LES_viscosity = opts%add_LES_viscosity
  cdesc%Kh_aniso = opts%Kh_aniso
  cdesc%res_scale_MEKE = opts%res_scale_MEKE
  cdesc%use_GME = opts%use_GME
  cdesc%answer_date = opts%answer_date
  cdesc%GME_h0 = opts%GME_h0
  cdesc%GME_efficiency = opts%GME_efficiency
  cdesc%GME_limiter = opts%GME_limiter
  cdesc%min_grid_Kh = opts%min_grid_Kh
  cdesc%min_grid_Ah = opts%min_grid_Ah
  cdesc%use_cont_thick = opts%use_cont_thick
  cdesc%use_ZB2020 = opts%use_ZB2020
  cdesc%use_circulation = opts%use_circulation
  cdesc%do_grid_Re_Ah = opts%do_grid_Re_Ah
  cdesc%do_grid_Re_Kh = opts%do_grid_Re_Kh
  cdesc%do_Ah_h = opts%do_Ah_h
  cdesc%do_Ah_q = opts%do_Ah_q
  cdesc%do_Kh_h = opts%do_Kh_h
  cdesc%do_Kh_q = opts%do_Kh_q
  cdesc%do_GME_coeff_h = opts%do_GME_coeff_h
  cdesc%do_GME_coeff_q = opts%do_GME_coeff_q
  cdesc%do_vort_xy_q = opts%do_vort_xy_q
  cdesc%do_div_xx_h = opts%do_div_xx_h
  cdesc%do_sh_xy_q = opts%do_sh_xy_q
  cdesc%do_sh_xx_h = opts%do_sh_xx_h
  cdesc%do_FrictWork = opts%do_FrictWork
  cdesc%do_FrictWorkIntz = opts%do_FrictWorkIntz
  cdesc%do_FrictWork_bh = opts%do_FrictWork_bh
  cdesc%do_FrictWorkIntz_bh = opts%do_FrictWorkIntz_bh
  cdesc%do_normstress = opts%do_normstress
  cdesc%do_shearstress = opts%do_shearstress
  cdesc%do_BS_coeff_h = opts%do_BS_coeff_h
  cdesc%do_BS_coeff_q = opts%do_BS_coeff_q
end function hv_opts_CS_to_c

!> Converts a hv_metrics_CS to its bind(C) mirror.
#ifdef _TIM
function hv_metrics_CS_to_c(metrics) result(cdesc)
  type(hv_metrics_CS), target, intent(in) :: metrics !< Metric and coefficient arrays to convert
  type(hv_metrics_CS_C)                   :: cdesc   !< bind(C) mirror of metrics
  cdesc%Kh_bg_xx = metrics%Kh_bg_xx%to_c()
  cdesc%Ah_bg_xx = metrics%Ah_bg_xx%to_c()
  cdesc%reduction_xx = metrics%reduction_xx%to_c()
  cdesc%Kh_Max_xx = metrics%Kh_Max_xx%to_c()
  cdesc%Ah_Max_xx = metrics%Ah_Max_xx%to_c()
  cdesc%Ah_Max_xx_KS = metrics%Ah_Max_xx_KS%to_c()
  cdesc%n1n2_h = metrics%n1n2_h%to_c()
  cdesc%n1n1_m_n2n2_h = metrics%n1n1_m_n2n2_h%to_c()
  cdesc%grid_sp_h2 = metrics%grid_sp_h2%to_c()
  cdesc%grid_sp_h3 = metrics%grid_sp_h3%to_c()
  cdesc%Kh_bg_xy = metrics%Kh_bg_xy%to_c()
  cdesc%Ah_bg_xy = metrics%Ah_bg_xy%to_c()
  cdesc%reduction_xy = metrics%reduction_xy%to_c()
  cdesc%Kh_Max_xy = metrics%Kh_Max_xy%to_c()
  cdesc%Ah_Max_xy = metrics%Ah_Max_xy%to_c()
  cdesc%Ah_Max_xy_KS = metrics%Ah_Max_xy_KS%to_c()
  cdesc%n1n2_q = metrics%n1n2_q%to_c()
  cdesc%n1n1_m_n2n2_q = metrics%n1n1_m_n2n2_q%to_c()
  cdesc%dx2h = metrics%dx2h%to_c()
  cdesc%dy2h = metrics%dy2h%to_c()
  cdesc%dx_dyT = metrics%dx_dyT%to_c()
  cdesc%dy_dxT = metrics%dy_dxT%to_c()
  cdesc%m_const_leithy = metrics%m_const_leithy%to_c()
  cdesc%m_leithy_max = metrics%m_leithy_max%to_c()
  cdesc%dx2q = metrics%dx2q%to_c()
  cdesc%dy2q = metrics%dy2q%to_c()
  cdesc%dx_dyBu = metrics%dx_dyBu%to_c()
  cdesc%dy_dxBu = metrics%dy_dxBu%to_c()
  cdesc%Idx2dyCu = metrics%Idx2dyCu%to_c()
  cdesc%Idxdy2u = metrics%Idxdy2u%to_c()
  cdesc%Idx2dyCv = metrics%Idx2dyCv%to_c()
  cdesc%Idxdy2v = metrics%Idxdy2v%to_c()
  cdesc%Laplac2_const_xx = metrics%Laplac2_const_xx%to_c()
  cdesc%Biharm6_const_xx = metrics%Biharm6_const_xx%to_c()
  cdesc%Laplac3_const_xx = metrics%Laplac3_const_xx%to_c()
  cdesc%Biharm_const_xx = metrics%Biharm_const_xx%to_c()
  cdesc%Biharm_const2_xx = metrics%Biharm_const2_xx%to_c()
  cdesc%Re_Ah_const_xx = metrics%Re_Ah_const_xx%to_c()
  cdesc%Laplac2_const_xy = metrics%Laplac2_const_xy%to_c()
  cdesc%Biharm6_const_xy = metrics%Biharm6_const_xy%to_c()
  cdesc%Laplac3_const_xy = metrics%Laplac3_const_xy%to_c()
  cdesc%Biharm_const_xy = metrics%Biharm_const_xy%to_c()
  cdesc%Biharm_const2_xy = metrics%Biharm_const2_xy%to_c()
  cdesc%Re_Ah_const_xy = metrics%Re_Ah_const_xy%to_c()
end function hv_metrics_CS_to_c
#endif
!> \namespace mom_hor_visc
!!
!! \section section_horizontal_viscosity Horizontal viscosity in MOM
!!
!! This module contains the subroutine horizontal_viscosity that calculates the
!! effects of horizontal viscosity, including parameterizations of the value of
!! the viscosity itself. Subroutine horizontal_viscosity calculates the acceleration due to
!! some combination of a biharmonic viscosity and a Laplacian viscosity. Either or
!! both may use a coefficient that depends on the shear and strain of the flow.
!! All metric terms are retained. The Laplacian is calculated as the divergence of
!! a stress tensor, using the form suggested by \cite Smagorinsky1993. The biharmonic
!! is calculated by twice applying the divergence of the stress tensor that is
!! used to calculate the Laplacian, but without the dependence on thickness in the
!! first pass. This form permits a variable viscosity, and indicates no
!! acceleration for either resting fluid or solid body rotation.
!!
!! The form of the viscous accelerations is discussed extensively in \cite griffies2000,
!! and the implementation here follows that discussion closely.
!! We use the notation of \cite Smith2003 with the exception that the
!! isotropic viscosity is \f$\kappa_h\f$.
!!
!! In general, the horizontal stress tensor can be written as
!! \f[
!! {\bf \sigma} =
!! \begin{pmatrix}
!! \frac{1}{2} \left( \sigma_D + \sigma_T \right) & \frac{1}{2} \sigma_S \\\\
!! \frac{1}{2} \sigma_S & \frac{1}{2} \left( \sigma_D - \sigma_T \right)
!! \end{pmatrix}
!! \f]
!! where \f$\sigma_D\f$, \f$\sigma_T\f$ and \f$\sigma_S\f$ are stresses associated with
!! invariant factors in the strain-rate tensor. For a Newtonian fluid, the stress
!! tensor is usually linearly related to the strain-rate tensor. The horizontal
!! strain-rate tensor is
!! \f[
!! \dot{\bf e} =
!! \begin{pmatrix}
!! \frac{1}{2} \left( \dot{e}_D + \dot{e}_T \right) & \frac{1}{2} \dot{e}_S \\\\
!! \frac{1}{2} \dot{e}_S & \frac{1}{2} \left( \dot{e}_D - \dot{e}_T \right)
!! \end{pmatrix}
!! \f]
!! where \f$\dot{e}_D = \partial_x u + \partial_y v\f$ is the horizontal divergence,
!! \f$\dot{e}_T = \partial_x u - \partial_y v\f$ is the horizontal tension, and
!! \f$\dot{e}_S = \partial_y u + \partial_x v\f$ is the horizontal shear strain.
!!
!! The trace of the stress tensor, \f$tr(\bf \sigma) = \sigma_D\f$, is usually
!! absorbed into the pressure and only the deviatoric stress tensor considered.
!! From here on, we drop \f$\sigma_D\f$. The trace of the strain tensor, \f$tr(\bf e) =
!! \dot{e}_D\f$ is non-zero for horizontally divergent flow but only enters the
!! stress tensor through \f$\sigma_D\f$ and so we will drop \f$\sigma_D\f$ from
!! calculations of the strain tensor in the code. Therefore the horizontal stress
!! tensor can be considered to be
!! \f[
!! {\bf \sigma} =
!! \begin{pmatrix}
!! \frac{1}{2} \sigma_T & \frac{1}{2} \sigma_S \\\\
!! \frac{1}{2} \sigma_S & - \frac{1}{2} \sigma_T
!! \end{pmatrix}
!! .\f]
!!
!! The stresses above are linearly related to the strain through a viscosity
!! coefficient, \f$\kappa_h\f$:
!! \f{eqnarray*}{
!! \sigma_T & = & 2 \kappa_h \dot{e}_T \\\\
!! \sigma_S & = & 2 \kappa_h \dot{e}_S
!! .
!! \f}
!!
!! The viscosity \f$\kappa_h\f$ may either be a constant or variable. For example,
!! \f$\kappa_h\f$ may vary with the shear, as proposed by \cite Smagorinsky1993.
!!
!! The accelerations resulting form the divergence of the stress tensor are
!! \f{eqnarray*}{
!! \hat{\bf x} \cdot \left( \nabla \cdot {\bf \sigma} \right)
!! & = &
!! \partial_x \left( \frac{1}{2} \sigma_T \right)
!! + \partial_y \left( \frac{1}{2} \sigma_S \right)
!! \\\\
!! & = &
!! \partial_x \left( \kappa_h \dot{e}_T \right)
!! + \partial_y \left( \kappa_h \dot{e}_S \right)
!! \\\\
!! \hat{\bf y} \cdot \left( \nabla \cdot {\bf \sigma} \right)
!! & = &
!! \partial_x \left( \frac{1}{2} \sigma_S \right)
!! + \partial_y \left( - \frac{1}{2} \sigma_T \right)
!! \\\\
!! & = &
!! \partial_x \left( \kappa_h \dot{e}_S \right)
!! + \partial_y \left( - \kappa_h \dot{e}_T \right)
!! .
!! \f}
!!
!! The form of the Laplacian viscosity in general coordinates is:
!! \f{eqnarray*}{
!! \hat{\bf x} \cdot \left( \nabla \cdot \sigma \right)
!! & = &
!! \frac{1}{h} \left[ \partial_x \left( \kappa_h h \dot{e}_T \right)
!! + \partial_y \left( \kappa_h h \dot{e}_S \right) \right]
!! \\\\
!! \hat{\bf y} \cdot \left( \nabla \cdot \sigma \right)
!! & = &
!! \frac{1}{h} \left[ \partial_x \left( \kappa_h h \dot{e}_S \right)
!! - \partial_y \left( \kappa_h h \dot{e}_T \right) \right]
!! .
!! \f}
!!
!! \subsection section_laplacian_viscosity_coefficient Laplacian viscosity coefficient
!!
!! The horizontal viscosity coefficient, \f$\kappa_h\f$, can have multiple components.
!! The isotropic components are:
!!   - A uniform background component, \f$\kappa_{bg}\f$.
!!   - A constant but spatially variable 2D map, \f$\kappa_{2d}(x,y)\f$.
!!   - A ''MICOM'' viscosity, \f$U_\nu \Delta(x,y)\f$, which uses a constant
!! velocity scale, \f$U_\nu\f$ and a measure of the grid-spacing \f$\Delta(x,y)^2 =
!! \frac{2 \Delta x^2 \Delta y^2}{\Delta x^2 + \Delta y^2}\f$.
!!   - A function of
!! latitude, \f$\kappa_{\phi}(x,y) = \kappa_{\pi/2} |\sin(\phi)|^n\f$.
!!   - A dynamic Smagorinsky viscosity, \f$\kappa_{Sm}(x,y,t) = C_{Sm} \Delta^2 \sqrt{\dot{e}_T^2 + \dot{e}_S^2}\f$.
!!   - A dynamic Leith viscosity, \f$\kappa_{Lth}(x,y,t) =
!!                                    C_{Lth} \Delta^3 \sqrt{|\nabla \zeta|^2 + |\nabla \dot{e}_D|^2}\f$.
!!
!! A maximum stable viscosity, \f$\kappa_{max}(x,y)\f$ is calculated based on the
!! grid-spacing and time-step and used to clip calculated viscosities.
!!
!! The static components of \f$\kappa_h\f$ are first combined as follows:
!! \f[
!! \kappa_{static} = \min \left[ \max\left(
!! \kappa_{bg},
!! U_\nu \Delta(x,y),
!! \kappa_{2d}(x,y),
!! \kappa_\phi(x,y)
!! \right)
!! , \kappa_{max}(x,y) \right]
!! \f]
!! and stored in the module control structure as variables <code>Kh_bg_xx</code> and
!! <code>Kh_bg_xy</code> for the tension (h-points) and shear (q-points) components
!! respectively.
!!
!! The full viscosity includes the dynamic components as follows:
!! \f[
!! \kappa_h(x,y,t) = r(\Delta,L_d)
!! \max \left( \kappa_{static}, \kappa_{Sm}, \kappa_{Lth} \right)
!! \f]
!! where \f$r(\Delta,L_d)\f$ is a resolution function.
!!
!! The dynamic Smagorinsky and Leith viscosity schemes are exclusive with each
!! other.
!!
!! \subsection section_viscous_boundary_conditions Viscous boundary conditions
!!
!! Free slip boundary conditions have been coded, although no slip boundary
!! conditions can be used with the Laplacian viscosity based on the 2D land-sea
!! mask. For a western boundary, for example, the boundary conditions with the
!! biharmonic operator would be written as:
!! \f[
!!   \partial_x v = 0 ; \partial_x^3 v = 0 ; u = 0 ; \partial_x^2 u = 0 ,
!! \f]
!! while for a Laplacian operator, they are simply
!! \f[
!!   \partial_x v = 0 ; u = 0 .
!! \f]
!! These boundary conditions are largely dictated by the use of an Arakawa
!! C-grid and by the varying layer thickness.
!!
!! \subsection section_anisotropic_viscosity Anisotropic viscosity
!!
!! \cite Large2001 proposed enhancing viscosity in a particular direction and the
!! approach was generalized in \cite Smith2003. We use the second form of their
!! two coefficient anisotropic viscosity (section 4.3). We also replace their
!! \f$A^\prime\f$ and $D$ such that \f$2A^\prime = 2 \kappa_h + D\f$ and
!! \f$\kappa_a = D\f$ so that \f$\kappa_h\f$ can be considered the isotropic
!! viscosity and \f$\kappa_a=D\f$ can be consider the anisotropic viscosity. The
!! direction of anisotropy is defined by a unit vector \f$\hat{\bf
!! n}=(n_1,n_2)\f$.
!!
!! The contributions to the stress tensor are
!! \f[
!! \begin{pmatrix}
!! \sigma_T \\\\ \sigma_S
!! \end{pmatrix}
!! =
!! \left[
!! \begin{pmatrix}
!! 2 \kappa_h + \kappa_a & 0 \\\\
!! 0 & 2 \kappa_h
!! \end{pmatrix}
!! + 2 \kappa_a n_1 n_2
!! \begin{pmatrix}
!! - 2 n_1 n_2 & n_1^2 - n_2^2 \\\\
!! n_1^2 - n_2^2 & 2 n_1 n_2
!! \end{pmatrix}
!! \right]
!! \begin{pmatrix}
!! \dot{e}_T \\\\ \dot{e}_S
!! \end{pmatrix}
!! \f]
!! Dissipation of kinetic energy requires \f$\kappa_h \geq 0\f$ and \f$2 \kappa_h + \kappa_a \geq 0\f$.
!! Note that when anisotropy is aligned with the x-direction, \f$n_1 = \pm 1\f$, then
!! \f$n_2 = 0\f$ and the cross terms vanish. The accelerations in this aligned limit
!! with constant coefficients become
!! \f{eqnarray*}{
!! \hat{\bf x} \cdot \nabla \cdot {\bf \sigma}
!! & = &
!! \partial_x \left( \left( \kappa_h + \frac{1}{2} \kappa_a \right) \dot{e}_T \right)
!! + \partial_y \left( \kappa_h \dot{e}_S \right)
!! \\\\
!! & = &
!! \left( \kappa_h + \kappa_a \right) \partial_{xx} u
!! + \kappa_h \partial_{yy} u
!! - \frac{1}{2} \kappa_a \partial_x \left( \partial_x u + \partial_y v \right)
!! \\\\
!! \hat{\bf y} \cdot \nabla \cdot {\bf \sigma}
!! & = &
!! \partial_x \left( \kappa_h \dot{e}_S \right)
!! - \partial_y \left( \left( \kappa_h + \frac{1}{2} \kappa_a \right) \dot{e}_T \right)
!! \\\\
!! & = &
!! \kappa_h \partial_{xx} v
!! + \left( \kappa_h + \kappa_a \right) \partial_{yy} v
!! - \frac{1}{2} \kappa_a \partial_y \left( \partial_x u + \partial_y v \right)
!! \f}
!! which has contributions akin to a negative divergence damping (a divergence
!! enhancement?) but which is weaker than the enhanced tension terms by half.
!!
!! \subsection section_viscous_discretization Discretization
!!
!! The horizontal tension, \f$\dot{e}_T\f$, is stored in variable <code>sh_xx</code> and
!! discretized as
!! \f[
!! \dot{e}_T
!! = \frac{\Delta y}{\Delta x} \delta_i \left( \frac{1}{\Delta y} u \right)
!! - \frac{\Delta x}{\Delta y} \delta_j \left( \frac{1}{\Delta x} v \right)
!! .
!! \f]
!! The horizontal divergent strain, \f$\dot{e}_D\f$, is stored in variable
!! <code>div_xx</code> and discretized as
!! \f[
!! \dot{e}_D
!! = \frac{1}{h A} \left( \delta_i \left( \overline{h}^i \Delta y \, u \right)
!! + \delta_j \left( \overline{h}^j\Delta x \, v \right) \right)
!! .
!! \f]
!! Note that for expediency this is the exact discretization used in the
!! continuity equation.
!!
!! The horizontal shear strain, \f$\dot{e}_S\f$, is stored in variable <code>sh_xy</code>
!! and discretized as
!! \f[
!! \dot{e}_S = v_x + u_y
!! \f]
!! where
!! \f{align*}{
!! v_x &= \frac{\Delta y}{\Delta x} \delta_i \left( \frac{1}{\Delta y} v \right) \\\\
!! u_y &= \frac{\Delta x}{\Delta y} \delta_j \left( \frac{1}{\Delta x} u \right)
!! \f}
!! which are calculated separately so that no-slip or free-slip boundary
!! conditions can be applied to \f$v_x\f$ and \f$u_y\f$ where appropriate.
!!
!! The tendency for the x-component of the divergence of stress is stored in
!! variable <code>diffu</code> and discretized as
!! \f[
!! \hat{\bf x} \cdot \left( \nabla \cdot {\bf \sigma} \right) =
!! \frac{1}{A \overline{h}^i} \left(
!! \frac{1}{\Delta y} \delta_i \left( h \Delta y^2 \kappa_h \dot{e}_T \right) +
!! \frac{1}{\Delta x} \delta_j \left( \tilde{h}^{ij} \Delta x^2 \kappa_h \dot{e}_S \right)
!! \right)
!! .
!! \f]
!!
!! The tendency for the y-component of the divergence of stress is stored in
!! variable <code>diffv</code> and discretized as
!! \f[
!! \hat{\bf y} \cdot \left( \nabla \cdot {\bf \sigma} \right) =
!! \frac{1}{A \overline{h}^j} \left(
!! \frac{1}{\Delta y} \delta_i \left( \tilde{h}^{ij} \Delta y^2 A_M \dot{e}_S \right)
!! - \frac{1}{\Delta x} \delta_j \left( h \Delta x^2 A_M \dot{e}_T \right)
!! \right)
!! .
!! \f]
!!
!! \subsection section_viscous_refs References
!!
!! Griffies, S.M., and Hallberg, R.W., 2000: Biharmonic friction with a
!! Smagorinsky-like viscosity for use in large-scale eddy-permitting ocean models.
!! Monthly Weather Review, 128(8), 2935-2946.
!! https://doi.org/10.1175/1520-0493(2000)128%3C2935:BFWASL%3E2.0.CO;2
!!
!! Large, W.G., Danabasoglu, G., McWilliams, J.C., Gent, P.R. and Bryan, F.O.,
!! 2001: Equatorial circulation of a global ocean climate model with
!! anisotropic horizontal viscosity.
!! Journal of Physical Oceanography, 31(2), pp.518-536.
!! https://doi.org/10.1175/1520-0485(2001)031%3C0518:ECOAGO%3E2.0.CO;2
!!
!! Smagorinsky, J., 1993: Some historical remarks on the use of nonlinear
!! viscosities. Large eddy simulation of complex engineering and geophysical
!! flows, 1, 69-106.
!!
!! Smith, R.D., and McWilliams, J.C., 2003: Anisotropic horizontal viscosity for
!! ocean models. Ocean Modelling, 5(2), 129-156.
!! https://doi.org/10.1016/S1463-5003(02)00016-1
end module MOM_hor_visc
