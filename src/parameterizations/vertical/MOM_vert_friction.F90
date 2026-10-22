! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!!SKILLS: 0.3.2

#include "do_concurrent_compat.h"

!> Implements vertical viscosity (vertvisc)
module MOM_vert_friction

use MOM_domains,       only : pass_var, To_All, Omit_corners
use MOM_domains,       only : pass_vector, Scalar_Pair
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : post_product_u, post_product_sum_u
use MOM_diag_mediator, only : post_product_v, post_product_sum_v
use MOM_diag_mediator, only : diag_ctrl, query_averaging_enabled
use MOM_domains,       only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,       only : To_North, To_East
use MOM_debugging,     only : uvchksum, hchksum
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,  only : mech_forcing, find_ustar
use MOM_get_input,     only : directories
use MOM_grid,          only : ocean_grid_type
use MOM_grid_containers, only : grid_core_type, grid_core
use MOM_io,            only : MOM_read_data, slasher
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary, only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_PointAccel,    only : write_u_accel, write_v_accel, PointAccel_init
use MOM_PointAccel,    only : PointAccel_CS
use MOM_time_manager,  only : time_type, time_minus_signed
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, vertvisc_type
use MOM_variables,     only : cont_diag_ptrs, accel_diag_ptrs
use MOM_variables,     only : ocean_internal_state
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS
use MOM_set_visc,      only : set_v_at_u, set_u_at_v
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use array_mod,         only : RealArray_t, RealArray_C, LogicalArray_t
use box_mod,           only : Box_t
use iso_c_binding,     only : c_int, c_bool, c_double

use CVMix_kpp,         only : cvmix_kpp_composite_Gshape

implicit none ; private

#include <MOM_memory.h>

public vertvisc, vertvisc_remnant, vertvisc_coef
public vertvisc_limit_vel, vertvisc_init, vertvisc_end
public updateCFLtruncationValue
public vertFPmix

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Options used in the calculation of the coupling coefficients across interfaces.
type :: coupling_CS
  real    :: Hmix            !< The mixed layer thickness [Z ~> m].
  real    :: Kvml_invZ2      !< The extra vertical viscosity scale in [H Z T-1 ~> m2 s-1 or Pa s] in
                             !! a surface mixed layer with a characteristic thickness given by Hmix,
                             !! and scaling proportional to (Hmix/z)^2, where z is the distance from
                             !! the surface; this can get very large with thin layers.
  real    :: Kv              !< The interior vertical viscosity [H Z T-1 ~> m2 s-1 or Pa s].
  real    :: Kv_extra_bbl    !< An extra vertical viscosity in the bottom boundary layer of
                             !! thickness Hbbl when there is not a bottom drag law in use
                             !! [H Z T-1 ~> m2 s-1 or Pa s].
  real    :: vonKar          !< The von Karman constant as used for mixed layer viscosity [nondim]
  logical :: bottomdraglaw  !< If true, the  bottom stress is calculated with a
                            !! drag law c_drag*|u|*u. The velocity magnitude
                            !! may be an assumed value or it may be based on the
                            !! actual velocity in the bottommost HBBL, depending
                            !! on whether linear_drag is true.
  logical :: dynamic_viscous_ML  !< If true, use the results from a dynamic
                            !! calculation, perhaps based on a bulk Richardson
                            !! number criterion, to determine the mixed layer
                            !! thickness for viscosity.
  logical :: fixed_LOTW_ML  !< If true, use a Law-of-the-wall prescription for the mixed layer
                            !! viscosity within a boundary layer that is the lesser of Hmix and the
                            !! total depth of the ocean in a column.
  logical :: apply_LOTW_floor !< If true, use a Law-of-the-wall prescription to set a lower bound
                            !! on the viscous coupling between layers within the surface boundary
                            !! layer, based the distance of interfaces from the surface.  This only
                            !! acts when there are large changes in the thicknesses of successive
                            !! layers or when the viscosity is set externally and the wind stress
                            !! has subsequently increased.
  integer :: answer_date    !< The vintage of the order of arithmetic and expressions in the viscous
                            !! calculations.  Values below 20190101 recover the answers from the end
                            !! of 2018, while higher values use expressions that do not use an
                            !! arbitrary and hard-coded maximum viscous coupling coefficient between
                            !! layers.  In non-Boussinesq cases, values below 20230601 recover a
                            !! form of the viscosity within  the mixed layer that breaks up the
                            !! magnitude of the wind stress with BULKMIXEDLAYER, DYNAMIC_VISCOUS_ML
                            !! or FIXED_DEPTH_LOTW_ML, but not LOTW_VISCOUS_ML_FLOOR.
end type coupling_CS

!> bind(C) mirror of coupling_CS, field-for-field, same order.
type, bind(C) :: coupling_CS_C
  real(c_double) :: Hmix !< The mixed layer thickness [Z ~> m].
  real(c_double) :: Kvml_invZ2 !< The extra vertical viscosity scale in [H Z T-1 ~> m2 s-1 or Pa s]
                               !! in a surface mixed layer with a characteristic thickness given by
                               !! Hmix, and scaling proportional to (Hmix/z)^2, where z is the
                               !! distance from the surface; this can get very large with thin
                               !! layers.
  real(c_double) :: Kv !< The interior vertical viscosity [H Z T-1 ~> m2 s-1 or Pa s].
  real(c_double) :: Kv_extra_bbl !< An extra vertical viscosity in the bottom boundary layer of
                                 !! thickness Hbbl when there is not a bottom drag law in use
                                 !! [H Z T-1 ~> m2 s-1 or Pa s].
  real(c_double) :: vonKar !< The von Karman constant as used for mixed layer viscosity [nondim]
  logical(c_bool) :: bottomdraglaw !< If true, the  bottom stress is calculated with a
                                   !! drag law c_drag*|u|*u. The velocity magnitude
                                   !! may be an assumed value or it may be based on the
                                   !! actual velocity in the bottommost HBBL, depending
                                   !! on whether linear_drag is true.
  logical(c_bool) :: dynamic_viscous_ML !< If true, use the results from a dynamic
                                        !! calculation, perhaps based on a bulk Richardson
                                        !! number criterion, to determine the mixed layer
                                        !! thickness for viscosity.
  logical(c_bool) :: fixed_LOTW_ML !< If true, use a Law-of-the-wall prescription for the mixed
                                   !! layer viscosity within a boundary layer that is the lesser of
                                   !! Hmix and the total depth of the ocean in a column.
  logical(c_bool) :: apply_LOTW_floor !< If true, use a Law-of-the-wall prescription to set a lower
                                      !! bound on the viscous coupling between layers within the
                                      !! surface boundary layer, based the distance of interfaces
                                      !! from the surface.  This only acts when there are large
                                      !! changes in the thicknesses of successive layers or when the
                                      !! viscosity is set externally and the wind stress has
                                      !! subsequently increased.
  integer(c_int) :: answer_date !< The vintage of the order of arithmetic and expressions in the
                                !! viscous calculations.  Values below 20190101 recover the answers
                                !! from the end of 2018, while higher values use expressions that do
                                !! not use an arbitrary and hard-coded maximum viscous coupling
                                !! coefficient between layers.  In non-Boussinesq cases, values
                                !! below 20230601 recover a form of the viscosity within  the mixed
                                !! layer that breaks up the magnitude of the wind stress with
                                !! BULKMIXEDLAYER, DYNAMIC_VISCOUS_ML or FIXED_DEPTH_LOTW_ML, but
                                !! not LOTW_VISCOUS_ML_FLOOR.
end type coupling_CS_C

!> Options used in the GL90 vertical viscosity parameterization.
type :: GL90_CS
  logical :: use_GL90_N2     !< If true, use GL90 vertical viscosity coefficient that is
                             !! depth-independent; this corresponds to a kappa_GM that scales as N^2
                             !! with depth.
  real    :: kappa_gl90      !< The scalar diffusivity used in the GL90 vertical viscosity scheme
                             !! [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
  logical :: read_kappa_gl90 !< If true, read a file containing the spatially varying kappa_gl90
  real    :: alpha_gl90      !< Coefficient used to compute a depth-independent GL90 vertical
                             !! viscosity via Kv_gl90 = alpha_gl90 * f^2. Note that the implied
                             !! Kv_gl90 corresponds to a kappa_gl90 that scales as N^2 with depth.
                             !! [H Z T ~> m2 s or kg s m-1]
  type(RealArray_t) :: kappa_gl90_2d !< 2D kappa_gl90 at h-points
                                     !! [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
end type GL90_CS

!> bind(C) mirror of GL90_CS, field-for-field, same order.
type, bind(C) :: GL90_CS_C
  logical(c_bool) :: use_GL90_N2 !< If true, use GL90 vertical viscosity coefficient that is
                                 !! depth-independent; this corresponds to a kappa_GM that scales as
                                 !! N^2 with depth.
  real(c_double) :: kappa_gl90 !< The scalar diffusivity used in the GL90 vertical viscosity scheme
                               !! [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
  logical(c_bool) :: read_kappa_gl90 !< If true, read a file containing the spatially varying
                                     !! kappa_gl90
  real(c_double) :: alpha_gl90 !< Coefficient used to compute a depth-independent GL90 vertical
                               !! viscosity via Kv_gl90 = alpha_gl90 * f^2. Note that the implied
                               !! Kv_gl90 corresponds to a kappa_gl90 that scales as N^2 with depth.
                               !! [H Z T ~> m2 s or kg s m-1]
  type(RealArray_C) :: kappa_gl90_2d !< 2D kappa_gl90 at h-points
                                     !! [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
end type GL90_CS_C

!> Options used in the vertical viscosity coupling coefficient calculations.
type :: coef_opts_CS
  type(coupling_CS) :: coupling_CS !< Options used in the calculation of the coupling coefficients
                                   !! across interfaces.
  real    :: Hbbl            !< The static bottom boundary layer thickness [Z ~> m].
  real    :: Hbbl_gl90       !< The static bottom boundary layer thickness used for GL90 [Z ~> m].
  logical :: use_GL90_in_SSW !< If true, use the GL90 parameterization in stacked shallow water mode
                             !! (SSW). The calculation of the GL90 viscosity coefficient uses the
                             !! fact that in SSW we simply have 1/N^2 = h/g^prime, where g^prime is
                             !! the reduced gravity. This identity does not generalize to non-SSW
                             !! setups.
  type(GL90_CS) :: GL90_CS !< Options used in the GL90 vertical viscosity parameterization.
  logical :: harmonic_visc  !< If true, the harmonic mean thicknesses are used
                            !! to calculate the viscous coupling between layers
                            !! except near the bottom.  Otherwise the arithmetic
                            !! mean thickness is used except near the bottom.
  real    :: harm_BL_val    !< A scale to determine when water is in the boundary
                            !! layers based solely on harmonic mean thicknesses
                            !! for the purpose of determining the extent to which
                            !! the thicknesses used in the viscosities are upwinded [nondim].
end type coef_opts_CS

!> bind(C) mirror of coef_opts_CS, field-for-field, same order.
type, bind(C) :: coef_opts_CS_C
  type(coupling_CS_C) :: coupling_CS !< Options used in the calculation of the coupling coefficients
                                     !! across interfaces.
  real(c_double) :: Hbbl !< The static bottom boundary layer thickness [Z ~> m].
  real(c_double) :: Hbbl_gl90 !< The static bottom boundary layer thickness used for GL90 [Z ~> m].
  logical(c_bool) :: use_GL90_in_SSW !< If true, use the GL90 parameterization in stacked shallow
                                     !! water mode (SSW). The calculation of the GL90 viscosity
                                     !! coefficient uses the fact that in SSW we simply have 1/N^2 =
                                     !! h/g^prime, where g^prime is the reduced gravity. This
                                     !! identity does not generalize to non-SSW setups.
  type(GL90_CS_C) :: GL90_CS !< Options used in the GL90 vertical viscosity parameterization.
  logical(c_bool) :: harmonic_visc !< If true, the harmonic mean thicknesses are used
                                   !! to calculate the viscous coupling between layers
                                   !! except near the bottom.  Otherwise the arithmetic
                                   !! mean thickness is used except near the bottom.
  real(c_double) :: harm_BL_val !< A scale to determine when water is in the boundary
                                !! layers based solely on harmonic mean thicknesses
                                !! for the purpose of determining the extent to which
                                !! the thicknesses used in the viscosities are upwinded [nondim].
end type coef_opts_CS_C

!> The vertical grid and unit-scaling constants used in the vertical viscosity coupling
!! coefficient calculations.
type :: vgrid_consts_CS
  integer :: ke !< The number of layers/levels in the vertical
  logical :: Boussinesq !< If true, make the Boussinesq approximation.
  integer :: nkml !< The number of layers at the top that should be treated
                  !! as parts of a homogeneous region.
  real :: dZ_subroundoff !< A thickness in height units that is so small that it can be added to a
                         !! vertical distance of Angstrom_Z or 1e-17 m without changing it at the
                         !! bit level [Z ~> m].  This is the height equivalent of H_subroundoff.
  real :: H_subroundoff !< A thickness that is so small that it can be added to a thickness of
                        !! Angstrom or larger without changing it at the bit level
                        !! [H ~> m or kg m-2]. If Angstrom is 0 or exceedingly small, this is
                        !! negligible compared to 1e-17 m.
  real :: H_to_m !< A constant that translates distances in the units of thickness
                 !! to m [m H-1 ~> 1 or m3 kg-1].
  real :: m_to_H !< A constant that translates distances in m to the units of
                 !! thickness [H m-1 ~> 1 or kg m-3].
  real :: RZ_to_H !< A constant that translates mass per unit area units to
                  !! thickness units [H R-1 Z-1 ~> m3 kg-2 or 1].
  real :: Z_to_H !< A constant that translates depth units to thickness units
                 !! depth [H Z-1 ~> 1 or kg m-3].
  real :: T_to_s !< A constant that translates the units of time to seconds
                 !! [s T-1 ~> 1]
  real :: Angstrom_H !< A one-Angstrom thickness in the model thickness units [H ~> m or kg m-2].
  real :: H_to_RZ !< A constant that translates thickness units to the units of
                  !! mass per unit area [R Z H-1 ~> kg m-3 or 1].
  real :: m_s_to_L_T !< Convert lateral velocities from m s-1 to L T-1
                     !! [L s T-1 m-1 ~> 1]
  type(RealArray_t) :: g_prime_a !< The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
end type vgrid_consts_CS

!> bind(C) mirror of vgrid_consts_CS, field-for-field, same order.
type, bind(C) :: vgrid_consts_C
  integer(c_int) :: ke !< The number of layers/levels in the vertical
  logical(c_bool) :: Boussinesq !< If true, make the Boussinesq approximation.
  integer(c_int) :: nkml !< The number of layers at the top that should be treated
                         !! as parts of a homogeneous region.
  real(c_double) :: dZ_subroundoff !< A thickness in height units that is so small that it can be
                                   !! added to a vertical distance of Angstrom_Z or 1e-17 m without
                                   !! changing it at the bit level [Z ~> m].  This is the height
                                   !! equivalent of H_subroundoff.
  real(c_double) :: H_subroundoff !< A thickness that is so small that it can be added to a
                                  !! thickness of Angstrom or larger without changing it at the bit
                                  !! level [H ~> m or kg m-2]. If Angstrom is 0 or exceedingly
                                  !! small, this is negligible compared to 1e-17 m.
  real(c_double) :: H_to_m !< A constant that translates distances in the units of thickness
                           !! to m [m H-1 ~> 1 or m3 kg-1].
  real(c_double) :: m_to_H !< A constant that translates distances in m to the units of
                           !! thickness [H m-1 ~> 1 or kg m-3].
  real(c_double) :: RZ_to_H !< A constant that translates mass per unit area units to
                            !! thickness units [H R-1 Z-1 ~> m3 kg-2 or 1].
  real(c_double) :: Z_to_H !< A constant that translates depth units to thickness units
                           !! depth [H Z-1 ~> 1 or kg m-3].
  real(c_double) :: T_to_s !< A constant that translates the units of time to seconds
                           !! [s T-1 ~> 1]
  real(c_double) :: Angstrom_H !< A one-Angstrom thickness in the model thickness units
                               !! [H ~> m or kg m-2].
  real(c_double) :: H_to_RZ !< A constant that translates thickness units to the units of
                            !! mass per unit area [R Z H-1 ~> kg m-3 or 1].
  real(c_double) :: m_s_to_L_T !< Convert lateral velocities from m s-1 to L T-1
                               !! [L s T-1 m-1 ~> 1]
  type(RealArray_C) :: g_prime_a !< The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
end type vgrid_consts_C

!> Options used in the vertical viscosity solver and the velocity truncations.
type :: vertvisc_opts_CS
  real    :: Hmix_stress     !< The mixed layer thickness over which the wind
                             !! stress is applied with direct_stress [H ~> m or kg m-2].
  logical :: direct_stress  !< If true, the wind stress is distributed over the topmost Hmix_stress
                            !! of fluid, and an added mixed layer viscosity or a physically based
                            !! boundary layer turbulence parameterization is not needed for
                            !! stability.
  logical :: StokesMixing   !< If true, do Stokes drift mixing via the Lagrangian current
                            !! (Eulerian plus Stokes drift).  False by default and set
                            !! via STOKES_MIXING_COMBINED.
  real    :: vel_underflow   !< Velocity components smaller than vel_underflow
                             !! are set to 0 [L T-1 ~> m s-1].
  real    :: CFL_trunc       !< Velocity components will be truncated when they
                             !! are large enough that the corresponding CFL number
                             !! exceeds this value [nondim].
  real    :: CFL_report      !< The value of the CFL number that will cause the
                             !! accelerations to be reported [nondim].  CFL_report
                             !! will often equal CFL_trunc.
  logical :: do_report_u = .false. !< If true, the columns with large u-velocities are reported.
  logical :: do_report_v = .false. !< If true, the columns with large v-velocities are reported.
  logical :: do_GLwork = .false. !< If true, the KE source from GL90 vertical viscosity is needed
                                 !! for diagnostics.
  logical :: do_du_dt_visc_gl90 = .false. !< If true, the zonal acceleration from GL90 vertical
                                          !! viscosity is needed for diagnostics.
  logical :: do_dv_dt_visc_gl90 = .false. !< If true, the meridional acceleration from GL90 vertical
                                          !! viscosity is needed for diagnostics.
end type vertvisc_opts_CS

!> bind(C) mirror of vertvisc_opts_CS, field-for-field, same order.
type, bind(C) :: vertvisc_opts_CS_C
  real(c_double) :: Hmix_stress !< The mixed layer thickness over which the wind
                                !! stress is applied with direct_stress [H ~> m or kg m-2].
  logical(c_bool) :: direct_stress !< If true, the wind stress is distributed over the topmost
                                   !! Hmix_stress of fluid, and an added mixed layer viscosity or a
                                   !! physically based boundary layer turbulence parameterization is
                                   !! not needed for stability.
  logical(c_bool) :: StokesMixing !< If true, do Stokes drift mixing via the Lagrangian current
                                  !! (Eulerian plus Stokes drift).  False by default and set
                                  !! via STOKES_MIXING_COMBINED.
  real(c_double) :: vel_underflow !< Velocity components smaller than vel_underflow
                                  !! are set to 0 [L T-1 ~> m s-1].
  real(c_double) :: CFL_trunc !< Velocity components will be truncated when they
                              !! are large enough that the corresponding CFL number
                              !! exceeds this value [nondim].
  real(c_double) :: CFL_report !< The value of the CFL number that will cause the
                               !! accelerations to be reported [nondim].  CFL_report
                               !! will often equal CFL_trunc.
  logical(c_bool) :: do_report_u !< If true, the columns with large u-velocities are reported.
  logical(c_bool) :: do_report_v !< If true, the columns with large v-velocities are reported.
  logical(c_bool) :: do_GLwork !< If true, the KE source from GL90 vertical viscosity is needed for
                               !! diagnostics.
  logical(c_bool) :: do_du_dt_visc_gl90 !< If true, the zonal acceleration from GL90 vertical
                                        !! viscosity is needed for diagnostics.
  logical(c_bool) :: do_dv_dt_visc_gl90 !< If true, the meridional acceleration from GL90 vertical
                                        !! viscosity is needed for diagnostics.
end type vertvisc_opts_CS_C

!> The velocities and flags used to report the columns in which velocity truncations occur,
!! each allocated only if those reports are being written.
type :: limit_vel_report_type
  type(RealArray_t)    :: u_old_a        !< The previous u-velocity [L T-1 ~> m s-1]
  type(RealArray_t)    :: v_old_a        !< The previous v-velocity [L T-1 ~> m s-1]
  type(RealArray_t)    :: vel_report_u_a !< The u-velocity to report [L T-1 ~> m s-1]
  type(RealArray_t)    :: vel_report_v_a !< The v-velocity to report [L T-1 ~> m s-1]
  type(LogicalArray_t) :: dowrite_u_a    !< True for the u-point columns to report
  type(LogicalArray_t) :: dowrite_v_a    !< True for the v-point columns to report
end type limit_vel_report_type

!> Copies of the acceleration diagnostics that the vertical viscosity calculations fill, each
!! allocated only if its accel_diag_ptrs counterpart is associated.
type :: vertvisc_AD_type
  type(RealArray_t) :: du_dt_visc_a      !< Zonal acceleration due to vertical viscosity
                                         !! [L T-2 ~> m s-2]
  type(RealArray_t) :: dv_dt_visc_a      !< Meridional acceleration due to vertical viscosity
                                         !! [L T-2 ~> m s-2]
  type(RealArray_t) :: du_dt_visc_gl90_a !< Zonal acceleration due to GL90 vertical viscosity
                                         !! (is included in du_dt_visc) [L T-2 ~> m s-2]
  type(RealArray_t) :: dv_dt_visc_gl90_a !< Meridional acceleration due to GL90 vertical viscosity
                                         !! (is included in dv_dt_visc) [L T-2 ~> m s-2]
  type(RealArray_t) :: du_dt_str_a       !< Zonal acceleration due to the surface stress (included
                                         !! in du_dt_visc) [L T-2 ~> m s-2]
  type(RealArray_t) :: dv_dt_str_a       !< Meridional acceleration due to the surface stress
                                         !! (included in dv_dt_visc) [L T-2 ~> m s-2]
contains
  procedure :: build_from => vertvisc_AD_build_from !< Copy in the associated fields
  procedure :: copy_back  => vertvisc_AD_copy_back  !< Copy out the allocated fields and free
end type vertvisc_AD_type

!> The viscous coupling coefficients and effective layer thicknesses at velocity points that
!! the coupling coefficient calculations store for the vertical viscosity solvers.
type :: coupling_state_type
  type(RealArray_t) :: a_u !< The u-drag coefficient across an interface
                           !! [H T-1 ~> m s-1 or Pa s m-1]
  type(RealArray_t) :: a_u_gl90 !< The u-drag coefficient associated with GL90 across an interface
                                !! [H T-1 ~> m s-1 or Pa s m-1]
  type(RealArray_t) :: h_u !< The effective layer thickness at u-points [H ~> m or kg m-2].
  type(RealArray_t) :: a_v !< The v-drag coefficient across an interface
                           !! [H T-1 ~> m s-1 or Pa s m-1]
  type(RealArray_t) :: a_v_gl90 !< The v-drag coefficient associated with GL90 across an interface
                                !! [H T-1 ~> m s-1 or Pa s m-1]
  type(RealArray_t) :: h_v !< The effective layer thickness at v-points [H ~> m or kg m-2].
  type(RealArray_t) :: a1_shelf_u !< The u-momentum coupling coefficient under ice shelves
                                  !! [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress
                                  !! under shelves.
  type(RealArray_t) :: a1_shelf_v !< The v-momentum coupling coefficient under ice shelves
                                  !! [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress
                                  !! under shelves.
end type coupling_state_type

!> The control structure with parameters and memory for the MOM_vert_friction module
type, public :: vertvisc_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  type(coef_opts_CS) :: coef_opts_CS !< Options used in the vertical viscosity coupling
                                     !! coefficient calculations.
  type(vertvisc_opts_CS) :: vertvisc_opts_CS !< Options used in the vertical viscosity solver
                                             !! and the velocity truncations.

  real    :: truncRampTime   !< The time-scale over which to ramp up the value of
                             !! CFL_trunc from CFL_truncS to CFL_truncE [T ~> s]
  real    :: CFL_truncS      !< The start value of CFL_trunc [nondim]
  real    :: CFL_truncE      !< The end/target value of CFL_trunc [nondim]
  logical :: CFLrampingIsActivated = .false. !< True if the ramping has been initialized
  type(time_type) :: rampStartTime !< The time at which the ramping of CFL_trunc starts

  type(coupling_state_type) :: coupling_state !< The viscous coupling coefficients and
                                              !! effective layer thicknesses at velocity points

  logical :: split          !< If true, use the split time stepping scheme.
  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  integer :: nkml           !< The number of layers in the mixed layer.
  integer, pointer :: ntrunc !< The number of times the velocity has been
                            !! truncated since the last call to write_energy.
  character(len=200) :: u_trunc_file  !< The complete path to a file in which a column of
                            !! u-accelerations are written if velocity truncations occur.
  character(len=200) :: v_trunc_file !< The complete path to a file in which a column of
                            !! v-accelerations are written if velocity truncations occur.

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  !>@{ Diagnostic identifiers
  integer :: id_du_dt_visc = -1, id_dv_dt_visc = -1, id_du_dt_visc_gl90 = -1, id_dv_dt_visc_gl90 = -1
  integer :: id_GLwork = -1
  integer :: id_au_vv = -1, id_av_vv = -1, id_au_gl90_vv = -1, id_av_gl90_vv = -1
  integer :: id_du_dt_str = -1, id_dv_dt_str = -1
  integer :: id_h_u = -1, id_h_v = -1, id_hML_u = -1 , id_hML_v = -1
  integer :: id_Omega_w2x = -1, id_FPtau2s  = -1 , id_FPtau2w = -1
  integer :: id_uE_h  = -1, id_vE_h  = -1
  integer :: id_uStk  = -1, id_vStk  = -1
  integer :: id_uStk0 = -1, id_vStk0 = -1
  integer :: id_uInc_h= -1, id_vInc_h= -1
  integer :: id_taux_bot = -1, id_tauy_bot = -1
  integer :: id_Kv_slow = -1, id_Kv_u = -1, id_Kv_v = -1
  integer :: id_Kv_gl90_u = -1, id_Kv_gl90_v = -1
  ! integer :: id_hf_du_dt_visc    = -1, id_hf_dv_dt_visc    = -1
  integer :: id_h_du_dt_visc    = -1, id_h_dv_dt_visc    = -1
  integer :: id_hf_du_dt_visc_2d = -1, id_hf_dv_dt_visc_2d = -1
  integer :: id_h_du_dt_str    = -1, id_h_dv_dt_str    = -1
  integer :: id_du_dt_str_visc_rem = -1, id_dv_dt_str_visc_rem = -1
  !>@}

  type(PointAccel_CS), pointer :: PointAccel_CSp => NULL() !< A pointer to the control structure
                              !! for recording accelerations leading to velocity truncations

  type(group_pass_type) :: pass_KE_uv !< A handle used for group halo passes
end type vertvisc_CS

!> Copies of the vertvisc_type fields used by the vertvisc, vertvisc_remnant and vertvisc_coef
!! calculations, each allocated only if its vertvisc_type counterpart is allocated.
type, public :: vertvisc_container_type
  type(RealArray_t) :: Ray_u_a             !< The Rayleigh drag velocity to be applied to each layer
                                           !! at u-points [H T-1 ~> m s-1 or Pa s m-1].
  type(RealArray_t) :: Ray_v_a             !< The Rayleigh drag velocity to be applied to each layer
                                           !! at v-points [H T-1 ~> m s-1 or Pa s m-1].
  type(RealArray_t) :: taux_shelf_a        !< The zonal stresses on the ocean under shelves
                                           !! [R Z L T-2 ~> Pa].
  type(RealArray_t) :: tauy_shelf_a        !< The meridional stresses on the ocean under shelves
                                           !! [R Z L T-2 ~> Pa].
  type(RealArray_t) :: kv_bbl_u_a          !< The bottom boundary layer viscosity at the u-points
                                           !! [H Z T-1 ~> m2 s-1 or Pa s]
  type(RealArray_t) :: kv_bbl_v_a          !< The bottom boundary layer viscosity at the v-points
                                           !! [H Z T-1 ~> m2 s-1 or Pa s]
  type(RealArray_t) :: bbl_thick_u_a       !< The bottom boundary layer thickness at the u-points
                                           !! [Z ~> m].
  type(RealArray_t) :: bbl_thick_v_a       !< The bottom boundary layer thickness at the v-points
                                           !! [Z ~> m].
  type(RealArray_t) :: tbl_thick_shelf_u_a !< Thickness of the viscous top boundary layer under ice
                                           !! shelves at u-points [Z ~> m].
  type(RealArray_t) :: tbl_thick_shelf_v_a !< Thickness of the viscous top boundary layer under ice
                                           !! shelves at v-points [Z ~> m].
  type(RealArray_t) :: kv_tbl_shelf_u_a    !< Viscosity in the viscous top boundary layer under ice
                                           !! shelves at u-points [H Z T-1 ~> m2 s-1 or Pa s]
  type(RealArray_t) :: kv_tbl_shelf_v_a    !< Viscosity in the viscous top boundary layer under ice
                                           !! shelves at v-points [H Z T-1 ~> m2 s-1 or Pa s]
  type(RealArray_t) :: nkml_visc_u_a       !< The number of layers in the viscous surface mixed
                                           !! layer at u-points [nondim]. This is not an integer
                                           !! because there may be fractional layers, and it is
                                           !! stored in terms of layers, not depth, to facilitate
                                           !! the movement of the viscous boundary layer with the
                                           !! flow.
  type(RealArray_t) :: nkml_visc_v_a       !< The number of layers in the viscous surface mixed
                                           !! layer at v-points [nondim].
  type(RealArray_t) :: Kv_shear_a          !< The shear-driven turbulent vertical viscosity at the
                                           !! interfaces between layers in tracer columns
                                           !! [H Z T-1 ~> m2 s-1 or Pa s]
  type(RealArray_t) :: Kv_shear_Bu_a       !< The shear-driven turbulent vertical viscosity at the
                                           !! interfaces between layers in corner columns
                                           !! [H Z T-1 ~> m2 s-1 or Pa s]
contains
  procedure :: build_from => vertvisc_container_build_from !< Copy in the allocated fields
  procedure :: copy_back  => vertvisc_container_copy_back !< Copy out the shelf stresses and free
  procedure :: free       => vertvisc_container_free !< Free every allocated field
end type vertvisc_container_type

!> Copies of the mech_forcing fields used by the vertvisc and vertvisc_coef calculations,
!! each allocated only if its mech_forcing counterpart is associated.
type, public :: mech_forcing_container_type
  type(RealArray_t) :: taux_a         !< zonal wind stress [R L Z T-2 ~> Pa]
  type(RealArray_t) :: tauy_a         !< meridional wind stress [R L Z T-2 ~> Pa]
  type(RealArray_t) :: frac_shelf_u_a !< Fractional ice shelf coverage of u-cells, nondimensional
                                      !! from 0 to 1 [nondim]. This is only associated if ice
                                      !! shelves are enabled, and is exactly 0 away from shelves or
                                      !! on land.
  type(RealArray_t) :: frac_shelf_v_a !< Fractional ice shelf coverage of v-cells, nondimensional
                                      !! from 0 to 1 [nondim]. This is only associated if ice
                                      !! shelves are enabled, and is exactly 0 away from shelves or
                                      !! on land.
contains
  procedure :: build_from => mech_forcing_container_build_from !< Copy in the associated fields
  procedure :: free       => mech_forcing_container_free !< Free every allocated field
end type mech_forcing_container_type

!> The diagnostic arrays that the vertical viscosity coupling coefficient calculations fill,
!! each allocated only if it is needed for a diagnostic or for debugging checksums.
type :: vertvisc_coef_diag_type
  type(RealArray_t) :: hML_u_a     !< Diagnostic of the mixed layer depth at u points [Z ~> m].
  type(RealArray_t) :: hML_v_a     !< Diagnostic of the mixed layer depth at v points [Z ~> m].
  type(RealArray_t) :: Kv_u_a      !< Total vertical viscosity at u-points in
                                   !! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  type(RealArray_t) :: Kv_v_a      !< Total vertical viscosity at v-points in
                                   !! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  type(RealArray_t) :: Kv_gl90_u_a !< GL90 vertical viscosity at u-points in
                                   !! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  type(RealArray_t) :: Kv_gl90_v_a !< GL90 vertical viscosity at v-points in
                                   !! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
end type vertvisc_coef_diag_type

contains

!> Add nonlocal stress increments to ui^n and vi^n.
subroutine vertFPmix(ui, vi, uold, vold, hbl_h, h, forces, dt, lpost, Cemp_NL, G, GV, US, CS, OBC, Waves)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: ui     !< Zonal velocity after vertvisc [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vi     !< Meridional velocity after vertvisc [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: uold   !< Old Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vold   !< Old Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: hbl_h !<  boundary layer depth [H ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h       !< Layer thicknesses [H ~> m or kg m-2]
  type(mech_forcing),      intent(in) :: forces  !< A structure with the driving mechanical forces
  real,                    intent(in) :: dt      !< Time increment [T ~> s]
  real,                    intent(in) :: Cemp_NL !< empirical coefficient of non-local momentum mixing [nondim]
  logical,                 intent(in) :: lpost   !< Compute and make available FPMix diagnostics
  type(unit_scale_type),   intent(in) :: US      !< A dimensional unit scaling type
  type(vertvisc_CS),       pointer    :: CS      !< Vertical viscosity control structure
  type(ocean_OBC_type),    pointer    :: OBC     !< Open boundary condition structure
  type(wave_parameters_CS), &
                   optional, pointer  :: Waves   !< Container for wave/Stokes information

  ! local variables
  real, dimension(SZIB_(G),SZJ_(G))  :: hbl_u   !< boundary layer depth (u-pts) [H ~> m]
  real, dimension(SZI_(G),SZJB_(G))  :: hbl_v   !< boundary layer depth (v-pts) [H ~> m]
  real, dimension(SZIB_(G),SZJ_(G))  :: taux_u  !< kinematic zonal wind stress (u-pts) [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJB_(G))  :: tauy_v  !< kinematic merid wind stress (v-pts) [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJ_(G))   :: uS0     !< surface zonal Stokes drift h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G))   :: vS0     !< surface zonal Stokes drift h-pts [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uE_u    !< zonal Eulerian u-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: uE_h    !< zonal Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vE_v    !< merid Eulerian v-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: vE_h    !< merid Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uInc_u  !< zonal Eulerian u-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: uInc_h  !< zonal Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vInc_v  !< merid Eulerian v-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: vInc_h  !< merid Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: uStk    !< zonal Stokes Drift (h-pts) [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: vStk    !< merid Stokes Drift (h-pts) [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)+1) :: omega_tau2s !< angle stress to shear (h-pts) [rad]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)+1) :: omega_tau2w !< angle stress to wind  (h-pts) [rad]
  real :: omega_tmp, omega_s2x, omega_tau2x                    !< temporary angle wrt the x axis [rad]
  real :: Irho0        !< Inverse of the mean density rescaled to [Z L-1 R-1 ~> m3 kg-1]
  real :: pi           !< ! The ratio of the circumference of a circle to its diameter [nondim]
  real :: tmp_u, tmp_v !< temporary ocean mask weights on u and v points [nondim]
  real :: fexp         !< temporary exponential function [nondim]
  real :: sigma        !< temporary normalize boundary layer coordinate [nondim]
  real :: Gat1, Gsig, dGdsig !< Shape parameters [nondim]
  real :: du, dv       !< Intermediate velocity differences [L T-1 ~> m s-1]
  real :: depth        !< Cumulative of thicknesses [H ~> m]
  real, dimension(:,:,:), contiguous, pointer :: h_u, h_v ! The effective layer thicknesses at
                                                    ! u- and v-points [H ~> m or kg m-2]
  integer :: b, kp1, k, nz !< band and vertical indices
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq !< horizontal indices

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = GV%ke

  call CS%coupling_state%h_u%view(h_u)
  call CS%coupling_state%h_v%view(h_v)

  pi = 4. * atan2(1.,1.)
  Irho0 = 1.0 / GV%Rho0

  call pass_var(hbl_h , G%Domain, halo=1)

  ! u-points
  do j = js,je
    do I = Isq,Ieq
      taux_u(I,j)  = forces%taux(I,j) * Irho0
      if ( (G%mask2dCu(I,j) > 0.5) ) then
        ! h to u-pts
        tmp_u  = MAX (1.0 ,(G%mask2dT(i,j) + G%mask2dT(i+1,j) ) )
        hbl_u(I,j) = ((G%mask2dT(i,j) * hbl_h(i,j)) + (G%mask2dT(i+1,j) * hbl_h(i+1,j))) / tmp_u
        depth   = 0.
        Gat1  = 0.
        do k=1, nz
          ! cell center
          depth = depth + 0.5*h_u(I,j,k)
          uE_u(I,j,k) = ui(I,j,k) - waves%Us_x(I,j,k)
          if ( depth < hbl_u(I,j) )     then
            sigma = depth / hbl_u(i,j)
            ! cell bottom
            depth = depth + 0.5*h_u(I,j,k)
            call cvmix_kpp_composite_Gshape(sigma,Gat1,Gsig,dGdsig)
            ! nonlocal boundary-layer increment
            uInc_u(I,j,k)  = dt * Cemp_NL * taux_u(I,j) * dGdsig / hbl_u(I,j)
            ui(I,j,k) = ui(I,j,k) + uInc_u(I,j,k)
          else
            uInc_u(I,j,k) = 0.0
          endif
        enddo
      else
        do k=1, nz
          uInc_u(I,j,k) = 0.0
        enddo
      endif
    enddo
  enddo

  ! v-points
  do J = Jsq,Jeq
    do i = is,ie
      tauy_v(i,J)  = forces%tauy(i,J) * Irho0
      if ( (G%mask2dCv(i,J) > 0.5) ) then
        ! h to v-pts
        tmp_v  = max( 1.0 ,(G%mask2dT(i,j) + G%mask2dT(i,j+1)))
        hbl_v(i,J) = (G%mask2dT(i,j) * hbl_h(i,J) + G%mask2dT(i,j+1) * hbl_h(i,j+1)) / tmp_v
        depth = 0.
        Gat1  = 0.
        do k=1, nz
          ! cell center
          depth = depth + 0.5* h_v(i,J,k)
          vE_v(i,J,k) = vi(i,J,k) - waves%Us_y(i,J,k)
          if ( depth < hbl_v(i,J) )    then
            sigma = depth / hbl_v(i,J)
            ! cell bottom
            depth = depth + 0.5* h_v(i,J,k)
            call cvmix_kpp_composite_Gshape(sigma,Gat1,Gsig,dGdsig)
            ! nonlocal boundary-layer increment
            vInc_v(i,J,k) = dt * Cemp_NL * tauy_v(i,J) * dGdsig / hbl_v(i,J)
            vi(i,J,k) = vi(i,J,k) + vInc_v(i,J,k)
          else
            vInc_v(i,J,k)  = 0.0
          endif
        enddo
      else
        do k=1, nz
          vInc_v(i,J,k)  = 0.0
        enddo
      endif
    enddo
  enddo

  ! Compute and store diagnostics, only during the corrector step.
  if (lpost)  then
    call pass_vector(uE_u  ,  vE_v  , G%Domain, To_All)
    call pass_vector(uInc_u, vInc_v , G%Domain, To_All)
    uStk = 0.0
    vStk = 0.0
    uS0  = 0.0
    vS0  = 0.0

    do j = js,je
      do i = is,ie
        if (G%mask2dT(i,j) > 0.5)  then
          ! u to h-pts
          tmp_u  = max( 1.0 ,(G%mask2dCu(i,j) + G%mask2dCu(i-1,j)))
          ! v to h-pts
          tmp_v  = max( 1.0 ,(G%mask2dCv(i,j) + G%mask2dCv(i,j-1)))
          do k = 1,nz
            uE_h(i,j,k)   = (G%mask2dCu(i,j) *   uE_u(i,j,k) + G%mask2dCu(i-1,j) *   uE_u(i-1,j,k)) / tmp_u
            uInc_h(i,j,k) = (G%mask2dCu(i,j) * uInc_u(i,j,k) + G%mask2dCu(i-1,j) * uInc_u(i-1,j,k)) / tmp_u
            vE_h(i,j,k)   = (G%mask2dCv(i,j) *   vE_v(i,j,k) + G%mask2dCv(i,j-1) *   vE_v(i,j-1,k)) / tmp_v
            vInc_h(i,j,k) = (G%mask2dCv(i,j) * vInc_v(i,j,k) + G%mask2dCv(i,j-1) * vInc_v(i,j-1,k)) / tmp_v
          enddo
          ! Wind, Stress and Shear align at surface
          Omega_tau2w(i,j,1) = 0.0
          Omega_tau2s(i,j,1) = 0.0
          do k = 1,nz
            kp1 = min( nz , k+1)
            du = uE_h(i,j,k) - uE_h(i,j,kp1)
            dv = vE_h(i,j,k) - vE_h(i,j,kp1)
            omega_s2x = atan2( dv , du )

            du = du + uInc_h(i,j,k) - uInc_h(i,j,kp1)
            dv = dv + vInc_h(i,j,k) - vInc_h(i,j,kp1)
            omega_tau2x = atan2( dv , du )

            omega_tmp = omega_tau2x - forces%omega_w2x(i,j)
            if ( (omega_tmp  >   pi   ) )  omega_tmp = omega_tmp - 2.*pi
            if ( (omega_tmp  < (0.-pi)) )  omega_tmp = omega_tmp + 2.*pi
            Omega_tau2w(i,j,kp1) = omega_tmp

            omega_tmp = omega_tau2x - omega_s2x
            if ( (omega_tmp  >   pi   ) )  omega_tmp = omega_tmp - 2.*pi
            if ( (omega_tmp  < (0.-pi)) )  omega_tmp = omega_tmp + 2.*pi
            Omega_tau2s(i,j,kp1) = omega_tmp

          enddo
        endif

        ! Stokes drift
        do b=1,waves%NumBands
          uS0(i,j)  = uS0(i,j) + waves%UStk_Hb(i,j,b)    ! or forces%UStkb(i,j,b)
          vS0(i,j)  = vS0(i,j) + waves%VStk_Hb(i,j,b)    ! or forces%VStkb(i,j,b)
        enddo
        depth = 0.0
        do k = 1,nz
          do b  = 1, waves%NumBands
            ! cell center
            fexp = exp(-2. * waves%WaveNum_Cen(b) * (depth+0.5*h(i,j,k)) )
            uStk(i,j,k) = uStk(i,j,k) + waves%UStk_Hb(i,j,b) * fexp
            vStk(i,j,k) = vStk(i,j,k) + waves%VStk_Hb(i,j,b) * fexp
          enddo
          ! cell bottom
          depth = depth + h(i,j,k)
        enddo
      enddo
    enddo

    ! post FPmix diagnostics
    if (CS%id_uE_h    > 0) call post_data(CS%id_uE_h     , uE_h   , CS%diag)
    if (CS%id_vE_h    > 0) call post_data(CS%id_vE_h   , vE_h   , CS%diag)
    if (CS%id_uInc_h  > 0) call post_data(CS%id_uInc_h , uInc_h , CS%diag)
    if (CS%id_vInc_h  > 0) call post_data(CS%id_vInc_h , vInc_h , CS%diag)
    if (CS%id_FPtau2s > 0) call post_data(CS%id_FPtau2s, Omega_tau2s, CS%diag)
    if (CS%id_FPtau2w > 0) call post_data(CS%id_FPtau2w, Omega_tau2w, CS%diag)
    if (CS%id_uStk0   > 0) call post_data(CS%id_uStk0  , uS0 , CS%diag)
    if (CS%id_vStk0   > 0) call post_data(CS%id_vStk0  , vS0    , CS%diag)
    if (CS%id_uStk    > 0) call post_data(CS%id_uStk   , uStk   , CS%diag)
    if (CS%id_vStk    > 0) call post_data(CS%id_vStk   , vStk   , CS%diag)
    if (CS%id_Omega_w2x > 0) call post_data(CS%id_Omega_w2x, forces%omega_w2x, CS%diag)

  endif

end subroutine vertFPmix


!> Compute coupling coefficient associated with vertical viscosity parameterization as in Greatbatch and Lamb
!! (1990), hereafter referred to as the GL90 vertical viscosity parameterization. This vertical viscosity scheme
!! redistributes momentum in the vertical, and is the equivalent of the Gent & McWilliams (1990) parameterization,
!! but in a TWA (thickness-weighted averaged) set of equations. The vertical viscosity coefficient nu is computed
!! from kappa_GM via thermal wind balance, and the following relation:
!! nu = kappa_GM * f^2 / N^2.
!! In the following subroutine kappa_GM is assumed either (a) constant or (b) horizontally varying. In both cases,
!! (a) and  (b), one can additionally impose an EBT structure in the vertical for kappa_GM.
!! A third possible formulation of nu is depth-independent:
!! nu = f^2 * alpha
!! The latter formulation would be equivalent to a kappa_GM that varies as N^2 with depth.
!! The vertical viscosity del_z ( nu del_z u) is applied to the momentum equation with stress-free boundary
!! conditions at the top and bottom.
!!
!! In SSW mode, we have 1/N^2 = h/g'. The coupling coefficient is therefore equal to
!! a_cpl_gl90 = nu / h = kappa_GM * f^2 / g'
!! or
!! a_cpl_gl90 = nu / h = f^2 * alpha / h

pure subroutine find_coupling_coef_gl90(a_cpl_gl90, hvel, i, j, z_i, ke, dZ_subroundoff, g_prime, &
                                        CoriolisBu, CS, kappa_gl90_2d, kdgl90_struct, work_on_u)
  !$omp declare target
  integer, intent(in) :: ke !< The number of layers/levels in the vertical
  real, intent(in) :: dZ_subroundoff !< A thickness in height units that is so small that it can be
                                     !! added to a vertical distance of Angstrom_Z or 1e-17 m
                                     !! without changing it at the bit level [Z ~> m].  This is the
                                     !! height equivalent of H_subroundoff.
  real, dimension(ke+1), intent(in) :: g_prime !< The reduced gravity at each interface
                                               !! [L2 Z-1 T-2 ~> m s-2].
  real, dimension(:,:), pointer, &
                       intent(in) :: CoriolisBu !< The Coriolis parameter at corner points
                                                !! [T-1 ~> s-1].
  real, dimension(ke), intent(in) :: hvel !< Distance between interfaces
                                                !! at velocity points [Z ~> m]
  integer, intent(in) :: i                      !< Column i-index
  integer, intent(in) :: j                      !< Column j-index
  real, dimension(ke+1), intent(in) :: z_i  !< Estimate of interface heights above the
                                                !! bottom, normalized by the GL90 bottom
                                                !! boundary layer thickness [nondim]
  real, dimension(ke+1),intent(out) :: a_cpl_gl90   !< Coupling coefficient associated
                                                !! with GL90 across interfaces; is not
                                                !! included in a_cpl [H T-1 ~> m s-1 or Pa s m-1].
  type(GL90_CS),     intent(in) :: CS           !< Vertical viscosity control structure
  real, dimension(:,:), pointer, &
                    intent(in) :: kappa_gl90_2d !< 2D kappa_gl90 at h-points
                                                !! [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:,:), pointer, &
                    intent(in) :: kdgl90_struct !< Vertical structure function used in GL90
                                                !! diffusivity [nondim]
  logical, intent(in) :: work_on_u              !< If true, u-points are being calculated,
                                                !! otherwise they are v-points.

  ! local variables
  logical :: kdgl90_use_vert_struct  ! use vertical structure for GL90 coefficient
  integer :: k, nz
  real    :: f2         !< Squared Coriolis parameter at a velocity grid point [T-2 ~> s-2].
  real    :: h_neglect  ! A vertical distance that is so small it is usually lost in roundoff error
                        ! and can be neglected [Z ~> m].
  real    :: botfn      ! A function that is 1 at the bottom and small far from it [nondim]
  real    :: z2         ! The distance from the bottom, normalized by Hbbl_gl90 [nondim]

  nz = ke
  h_neglect = dZ_subroundoff
  kdgl90_use_vert_struct = .false.

  kdgl90_use_vert_struct = associated(kdgl90_struct)

  a_cpl_gl90(:) = 0.

  do K=2,nz
    if (work_on_u) then
      ! compute coupling coefficient at u-points
      f2 = 0.25 * (CoriolisBu(I,J-1) + CoriolisBu(I,J))**2
      if (CS%use_GL90_N2) then
        a_cpl_gl90(K) = 2. * f2 * CS%alpha_gl90 / (hvel(k) + hvel(k-1) + h_neglect)
      else
        if (CS%read_kappa_gl90) then
          a_cpl_gl90(K) = f2 * 0.5 * (kappa_gl90_2d(i,j) + kappa_gl90_2d(i+1,j)) / g_prime(K)
        else
          a_cpl_gl90(K) = f2 * CS%kappa_gl90 / g_prime(K)
        endif
        if (kdgl90_use_vert_struct) then
          a_cpl_gl90(K) = a_cpl_gl90(K) * 0.5 &
              * (kdgl90_struct(i,j,k-1) + kdgl90_struct(i+1,j,k-1))
        endif
      endif
      ! botfn determines when a point is within the influence of the GL90 bottom boundary layer,
      ! going from 1 at the bottom to 0 in the interior.
      z2 = z_i(k)
      botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

      a_cpl_gl90(K) = a_cpl_gl90(K) * (1. - botfn)
    else
      ! compute viscosities at v-points
      f2 = 0.25 * (CoriolisBu(I-1,J) + CoriolisBu(I,J))**2

      if (CS%use_GL90_N2) then
        a_cpl_gl90(K) = 2. * f2 * CS%alpha_gl90 / (hvel(k) + hvel(k-1) + h_neglect)
      else
        if (CS%read_kappa_gl90) then
          a_cpl_gl90(K) = f2 * 0.5 * (kappa_gl90_2d(i,j) + kappa_gl90_2d(i,j+1)) / g_prime(K)
        else
          a_cpl_gl90(K) = f2 * CS%kappa_gl90 / g_prime(K)
        endif
        if (kdgl90_use_vert_struct) then
          a_cpl_gl90(K) = a_cpl_gl90(K) * 0.5 &
              * (kdgl90_struct(i,j,k-1) + kdgl90_struct(i,j+1,k-1))
        endif
      endif
      ! botfn determines when a point is within the influence of the GL90 bottom boundary layer,
      ! going from 1 at the bottom to 0 in the interior.
      z2 = z_i(k)
      botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

      a_cpl_gl90(K) = a_cpl_gl90(K) * (1. - botfn)
    endif
  enddo
end subroutine find_coupling_coef_gl90


!> Copies each allocated vertvisc_type field used by these calculations into this container set.
subroutine vertvisc_container_build_from(this, visc)
  class(vertvisc_container_type), intent(inout) :: this !< Copies of the viscosities and bottom drag
  type(vertvisc_type),            intent(in)    :: visc !< Viscosities and bottom drag

  if (allocated(visc%Ray_u)) &
    call this%Ray_u_a%alloc(lb=LBOUND(visc%Ray_u), ub=UBOUND(visc%Ray_u), source=visc%Ray_u)
  if (allocated(visc%Ray_v)) &
    call this%Ray_v_a%alloc(lb=LBOUND(visc%Ray_v), ub=UBOUND(visc%Ray_v), source=visc%Ray_v)
  if (allocated(visc%taux_shelf)) &
    call this%taux_shelf_a%alloc(lb=LBOUND(visc%taux_shelf), ub=UBOUND(visc%taux_shelf), &
                                 source=visc%taux_shelf)
  if (allocated(visc%tauy_shelf)) &
    call this%tauy_shelf_a%alloc(lb=LBOUND(visc%tauy_shelf), ub=UBOUND(visc%tauy_shelf), &
                                 source=visc%tauy_shelf)
  if (allocated(visc%kv_bbl_u)) &
    call this%kv_bbl_u_a%alloc(lb=LBOUND(visc%kv_bbl_u), ub=UBOUND(visc%kv_bbl_u), &
                               source=visc%kv_bbl_u)
  if (allocated(visc%kv_bbl_v)) &
    call this%kv_bbl_v_a%alloc(lb=LBOUND(visc%kv_bbl_v), ub=UBOUND(visc%kv_bbl_v), &
                               source=visc%kv_bbl_v)
  if (allocated(visc%bbl_thick_u)) &
    call this%bbl_thick_u_a%alloc(lb=LBOUND(visc%bbl_thick_u), ub=UBOUND(visc%bbl_thick_u), &
                                  source=visc%bbl_thick_u)
  if (allocated(visc%bbl_thick_v)) &
    call this%bbl_thick_v_a%alloc(lb=LBOUND(visc%bbl_thick_v), ub=UBOUND(visc%bbl_thick_v), &
                                  source=visc%bbl_thick_v)
  if (allocated(visc%tbl_thick_shelf_u)) &
    call this%tbl_thick_shelf_u_a%alloc(lb=LBOUND(visc%tbl_thick_shelf_u), &
                                        ub=UBOUND(visc%tbl_thick_shelf_u), &
                                        source=visc%tbl_thick_shelf_u)
  if (allocated(visc%tbl_thick_shelf_v)) &
    call this%tbl_thick_shelf_v_a%alloc(lb=LBOUND(visc%tbl_thick_shelf_v), &
                                        ub=UBOUND(visc%tbl_thick_shelf_v), &
                                        source=visc%tbl_thick_shelf_v)
  if (allocated(visc%kv_tbl_shelf_u)) &
    call this%kv_tbl_shelf_u_a%alloc(lb=LBOUND(visc%kv_tbl_shelf_u), &
                                     ub=UBOUND(visc%kv_tbl_shelf_u), &
                                     source=visc%kv_tbl_shelf_u)
  if (allocated(visc%kv_tbl_shelf_v)) &
    call this%kv_tbl_shelf_v_a%alloc(lb=LBOUND(visc%kv_tbl_shelf_v), &
                                     ub=UBOUND(visc%kv_tbl_shelf_v), &
                                     source=visc%kv_tbl_shelf_v)
  if (allocated(visc%nkml_visc_u)) &
    call this%nkml_visc_u_a%alloc(lb=LBOUND(visc%nkml_visc_u), ub=UBOUND(visc%nkml_visc_u), &
                                  source=visc%nkml_visc_u)
  if (allocated(visc%nkml_visc_v)) &
    call this%nkml_visc_v_a%alloc(lb=LBOUND(visc%nkml_visc_v), ub=UBOUND(visc%nkml_visc_v), &
                                  source=visc%nkml_visc_v)
  if (associated(visc%Kv_shear)) &
    call this%Kv_shear_a%alloc(lb=LBOUND(visc%Kv_shear), ub=UBOUND(visc%Kv_shear), &
                               source=visc%Kv_shear)
  if (associated(visc%Kv_shear_Bu)) &
    call this%Kv_shear_Bu_a%alloc(lb=LBOUND(visc%Kv_shear_Bu), ub=UBOUND(visc%Kv_shear_Bu), &
                                  source=visc%Kv_shear_Bu)
end subroutine vertvisc_container_build_from

!> Copies the shelf stresses, the only fields these calculations change, back into visc, then
!! frees every field in this container set.
subroutine vertvisc_container_copy_back(this, visc)
  class(vertvisc_container_type), intent(inout) :: this !< Copies of the viscosities and bottom drag
  type(vertvisc_type),            intent(inout) :: visc !< Viscosities and bottom drag

  if (this%taux_shelf_a%associated()) call this%taux_shelf_a%copy2F(visc%taux_shelf)
  if (this%tauy_shelf_a%associated()) call this%tauy_shelf_a%copy2F(visc%tauy_shelf)
  call this%free()
end subroutine vertvisc_container_copy_back

!> Frees every allocated field in this container set.
subroutine vertvisc_container_free(this)
  class(vertvisc_container_type), intent(inout) :: this !< Copies of the viscosities and bottom drag

  if (this%Ray_u_a%associated()) call this%Ray_u_a%free()
  if (this%Ray_v_a%associated()) call this%Ray_v_a%free()
  if (this%taux_shelf_a%associated()) call this%taux_shelf_a%free()
  if (this%tauy_shelf_a%associated()) call this%tauy_shelf_a%free()
  if (this%kv_bbl_u_a%associated()) call this%kv_bbl_u_a%free()
  if (this%kv_bbl_v_a%associated()) call this%kv_bbl_v_a%free()
  if (this%bbl_thick_u_a%associated()) call this%bbl_thick_u_a%free()
  if (this%bbl_thick_v_a%associated()) call this%bbl_thick_v_a%free()
  if (this%tbl_thick_shelf_u_a%associated()) call this%tbl_thick_shelf_u_a%free()
  if (this%tbl_thick_shelf_v_a%associated()) call this%tbl_thick_shelf_v_a%free()
  if (this%kv_tbl_shelf_u_a%associated()) call this%kv_tbl_shelf_u_a%free()
  if (this%kv_tbl_shelf_v_a%associated()) call this%kv_tbl_shelf_v_a%free()
  if (this%nkml_visc_u_a%associated()) call this%nkml_visc_u_a%free()
  if (this%nkml_visc_v_a%associated()) call this%nkml_visc_v_a%free()
  if (this%Kv_shear_a%associated()) call this%Kv_shear_a%free()
  if (this%Kv_shear_Bu_a%associated()) call this%Kv_shear_Bu_a%free()
end subroutine vertvisc_container_free

!> Copies each associated mech_forcing field used by these calculations into this container set.
subroutine mech_forcing_container_build_from(this, forces)
  class(mech_forcing_container_type), intent(inout) :: this   !< Copies of the mechanical forcing
  type(mech_forcing),                 intent(in)    :: forces !< A structure with the driving
                                                              !! mechanical forces

  if (associated(forces%taux)) &
    call this%taux_a%alloc(lb=LBOUND(forces%taux), ub=UBOUND(forces%taux), source=forces%taux)
  if (associated(forces%tauy)) &
    call this%tauy_a%alloc(lb=LBOUND(forces%tauy), ub=UBOUND(forces%tauy), source=forces%tauy)
  if (associated(forces%frac_shelf_u)) &
    call this%frac_shelf_u_a%alloc(lb=LBOUND(forces%frac_shelf_u), ub=UBOUND(forces%frac_shelf_u), &
                                   source=forces%frac_shelf_u)
  if (associated(forces%frac_shelf_v)) &
    call this%frac_shelf_v_a%alloc(lb=LBOUND(forces%frac_shelf_v), ub=UBOUND(forces%frac_shelf_v), &
                                   source=forces%frac_shelf_v)
end subroutine mech_forcing_container_build_from

!> Frees every allocated field in this container set.
subroutine mech_forcing_container_free(this)
  class(mech_forcing_container_type), intent(inout) :: this !< Copies of the mechanical forcing

  if (this%taux_a%associated()) call this%taux_a%free()
  if (this%tauy_a%associated()) call this%tauy_a%free()
  if (this%frac_shelf_u_a%associated()) call this%frac_shelf_u_a%free()
  if (this%frac_shelf_v_a%associated()) call this%frac_shelf_v_a%free()
end subroutine mech_forcing_container_free

!> Fills the vertical grid and unit-scaling constants used by the vertical viscosity
!! calculations.
subroutine vgrid_consts_build(vgrid, GV, US)
  type(vgrid_consts_CS),   intent(inout) :: vgrid !< The vertical grid and unit-scaling constants
  type(verticalGrid_type), intent(in)    :: GV    !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type

  vgrid%ke = GV%ke ; vgrid%Boussinesq = GV%Boussinesq ; vgrid%nkml = GV%nkml
  vgrid%dZ_subroundoff = GV%dZ_subroundoff ; vgrid%H_subroundoff = GV%H_subroundoff
  vgrid%H_to_m = GV%H_to_m ; vgrid%m_to_H = GV%m_to_H
  vgrid%RZ_to_H = GV%RZ_to_H ; vgrid%Z_to_H = GV%Z_to_H ; vgrid%T_to_s = US%T_to_s
  call vgrid%g_prime_a%alloc(lb=LBOUND(GV%g_prime), ub=UBOUND(GV%g_prime), source=GV%g_prime)
  vgrid%Angstrom_H = GV%Angstrom_H ; vgrid%H_to_RZ = GV%H_to_RZ
  vgrid%m_s_to_L_T = US%m_s_to_L_T
end subroutine vgrid_consts_build

!> Copies each associated acceleration diagnostic in AD into its container.
subroutine vertvisc_AD_build_from(this, AD)
  class(vertvisc_AD_type), intent(inout) :: this !< Copies of the acceleration diagnostics
  type(accel_diag_ptrs),   intent(in)    :: AD   !< Accelerations in the momentum
                                                 !! equations for diagnostics

  if (associated(AD%du_dt_visc)) &
    call this%du_dt_visc_a%alloc(lb=LBOUND(AD%du_dt_visc), ub=UBOUND(AD%du_dt_visc), &
                                 source=AD%du_dt_visc)
  if (associated(AD%dv_dt_visc)) &
    call this%dv_dt_visc_a%alloc(lb=LBOUND(AD%dv_dt_visc), ub=UBOUND(AD%dv_dt_visc), &
                                 source=AD%dv_dt_visc)
  if (associated(AD%du_dt_visc_gl90)) &
    call this%du_dt_visc_gl90_a%alloc(lb=LBOUND(AD%du_dt_visc_gl90), &
                                      ub=UBOUND(AD%du_dt_visc_gl90), source=AD%du_dt_visc_gl90)
  if (associated(AD%dv_dt_visc_gl90)) &
    call this%dv_dt_visc_gl90_a%alloc(lb=LBOUND(AD%dv_dt_visc_gl90), &
                                      ub=UBOUND(AD%dv_dt_visc_gl90), source=AD%dv_dt_visc_gl90)
  if (associated(AD%du_dt_str)) &
    call this%du_dt_str_a%alloc(lb=LBOUND(AD%du_dt_str), ub=UBOUND(AD%du_dt_str), &
                                source=AD%du_dt_str)
  if (associated(AD%dv_dt_str)) &
    call this%dv_dt_str_a%alloc(lb=LBOUND(AD%dv_dt_str), ub=UBOUND(AD%dv_dt_str), &
                                source=AD%dv_dt_str)
end subroutine vertvisc_AD_build_from

!> Copies each allocated container back into AD, then frees it.
subroutine vertvisc_AD_copy_back(this, AD)
  class(vertvisc_AD_type), intent(inout) :: this !< Copies of the acceleration diagnostics
  type(accel_diag_ptrs),   intent(inout) :: AD   !< Accelerations in the momentum
                                                 !! equations for diagnostics

  if (this%du_dt_visc_a%associated()) then
    call this%du_dt_visc_a%copy2F(AD%du_dt_visc) ; call this%du_dt_visc_a%free()
  endif
  if (this%dv_dt_visc_a%associated()) then
    call this%dv_dt_visc_a%copy2F(AD%dv_dt_visc) ; call this%dv_dt_visc_a%free()
  endif
  if (this%du_dt_visc_gl90_a%associated()) then
    call this%du_dt_visc_gl90_a%copy2F(AD%du_dt_visc_gl90) ; call this%du_dt_visc_gl90_a%free()
  endif
  if (this%dv_dt_visc_gl90_a%associated()) then
    call this%dv_dt_visc_gl90_a%copy2F(AD%dv_dt_visc_gl90) ; call this%dv_dt_visc_gl90_a%free()
  endif
  if (this%du_dt_str_a%associated()) then
    call this%du_dt_str_a%copy2F(AD%du_dt_str) ; call this%du_dt_str_a%free()
  endif
  if (this%dv_dt_str_a%associated()) then
    call this%dv_dt_str_a%copy2F(AD%dv_dt_str) ; call this%dv_dt_str_a%free()
  endif
end subroutine vertvisc_AD_copy_back

!> Converts a coupling_CS to its bind(C) mirror.
function coupling_CS_to_c(opts) result(cdesc)
  type(coupling_CS), intent(in) :: opts !< Options to convert
  type(coupling_CS_C) :: cdesc          !< bind(C) mirror of opts
  cdesc%Hmix = opts%Hmix
  cdesc%Kvml_invZ2 = opts%Kvml_invZ2
  cdesc%Kv = opts%Kv
  cdesc%Kv_extra_bbl = opts%Kv_extra_bbl
  cdesc%vonKar = opts%vonKar
  cdesc%bottomdraglaw = opts%bottomdraglaw
  cdesc%dynamic_viscous_ML = opts%dynamic_viscous_ML
  cdesc%fixed_LOTW_ML = opts%fixed_LOTW_ML
  cdesc%apply_LOTW_floor = opts%apply_LOTW_floor
  cdesc%answer_date = opts%answer_date
end function coupling_CS_to_c

!> Converts a vertvisc_opts_CS to its bind(C) mirror.
function vertvisc_opts_CS_to_c(opts) result(cdesc)
  type(vertvisc_opts_CS), intent(in) :: opts !< Options to convert
  type(vertvisc_opts_CS_C) :: cdesc          !< bind(C) mirror of opts
  cdesc%Hmix_stress = opts%Hmix_stress
  cdesc%direct_stress = opts%direct_stress
  cdesc%StokesMixing = opts%StokesMixing
  cdesc%vel_underflow = opts%vel_underflow
  cdesc%CFL_trunc = opts%CFL_trunc
  cdesc%CFL_report = opts%CFL_report
  cdesc%do_report_u = opts%do_report_u
  cdesc%do_report_v = opts%do_report_v
  cdesc%do_GLwork = opts%do_GLwork
  cdesc%do_du_dt_visc_gl90 = opts%do_du_dt_visc_gl90
  cdesc%do_dv_dt_visc_gl90 = opts%do_dv_dt_visc_gl90
end function vertvisc_opts_CS_to_c

#ifdef _TIM
!> Converts a GL90_CS to its bind(C) mirror.
function GL90_CS_to_c(opts) result(cdesc)
  type(GL90_CS), target, intent(in) :: opts !< Options to convert
  type(GL90_CS_C) :: cdesc          !< bind(C) mirror of opts
  cdesc%use_GL90_N2 = opts%use_GL90_N2
  cdesc%kappa_gl90 = opts%kappa_gl90
  cdesc%read_kappa_gl90 = opts%read_kappa_gl90
  cdesc%alpha_gl90 = opts%alpha_gl90
  cdesc%kappa_gl90_2d = opts%kappa_gl90_2d%to_c()
end function GL90_CS_to_c
#endif

#ifdef _TIM
!> Converts a coef_opts_CS to its bind(C) mirror.
function coef_opts_CS_to_c(opts) result(cdesc)
  type(coef_opts_CS), target, intent(in) :: opts !< Options to convert
  type(coef_opts_CS_C) :: cdesc          !< bind(C) mirror of opts
  cdesc%coupling_CS = coupling_CS_to_c(opts%coupling_CS)
  cdesc%Hbbl = opts%Hbbl
  cdesc%Hbbl_gl90 = opts%Hbbl_gl90
  cdesc%use_GL90_in_SSW = opts%use_GL90_in_SSW
  cdesc%GL90_CS = GL90_CS_to_c(opts%GL90_CS)
  cdesc%harmonic_visc = opts%harmonic_visc
  cdesc%harm_BL_val = opts%harm_BL_val
end function coef_opts_CS_to_c

!> Converts a vgrid_consts_CS to its bind(C) mirror.
function vgrid_consts_to_c(vgrid) result(cdesc)
  type(vgrid_consts_CS), target, intent(in) :: vgrid !< Constants to convert
  type(vgrid_consts_C) :: cdesc          !< bind(C) mirror of vgrid
  cdesc%ke = vgrid%ke
  cdesc%Boussinesq = vgrid%Boussinesq
  cdesc%nkml = vgrid%nkml
  cdesc%dZ_subroundoff = vgrid%dZ_subroundoff
  cdesc%H_subroundoff = vgrid%H_subroundoff
  cdesc%H_to_m = vgrid%H_to_m
  cdesc%m_to_H = vgrid%m_to_H
  cdesc%RZ_to_H = vgrid%RZ_to_H
  cdesc%Z_to_H = vgrid%Z_to_H
  cdesc%T_to_s = vgrid%T_to_s
  cdesc%Angstrom_H = vgrid%Angstrom_H
  cdesc%H_to_RZ = vgrid%H_to_RZ
  cdesc%m_s_to_L_T = vgrid%m_s_to_L_T
  cdesc%g_prime_a = vgrid%g_prime_a%to_c()
end function vgrid_consts_to_c
#endif


!> Perform a fully implicit vertical diffusion
!! of momentum.  Stress top and bottom boundary conditions are used.
!!
!! This is solving the tridiagonal system
!! \f[ \left(h_k + a_{k + 1/2} + a_{k - 1/2} + r_k\right) u_k^{n+1}
!!     = h_k u_k^n + a_{k + 1/2} u_{k+1}^{n+1} + a_{k - 1/2} u_{k-1}^{n+1} \f]
!! where \f$a_{k + 1/2} = \Delta t \nu_{k + 1/2} / h_{k + 1/2}\f$
!! is the <em>interfacial coupling thickness per time step</em>,
!! encompassing background viscosity as well as contributions from
!! enhanced mixed and bottom layer viscosities.
!! $r_k$ is a Rayleigh drag term due to channel drag.
!! There is an additional stress term on the right-hand side
!! if DIRECT_STRESS is true, applied to the surface layer.
subroutine vertvisc(u, v, h, forces, visc, dt, OBC, ADp, CDp, G, GV, US, CS, &
                    taux_bot, tauy_bot, fpmix, Waves)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  type(mech_forcing),    intent(in)      :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),   intent(inout)   :: visc   !< Viscosities and bottom drag
  real,                  intent(in)      :: dt     !< Time increment [T ~> s]
  type(ocean_OBC_type),  pointer         :: OBC    !< Open boundary condition structure
  type(accel_diag_ptrs), intent(inout)   :: ADp    !< Accelerations in the momentum
                                                   !! equations for diagnostics
  type(cont_diag_ptrs),  intent(inout)   :: CDp    !< Continuity equation terms
  type(vertvisc_CS)                      :: CS     !< Vertical viscosity control structure
  real, dimension(SZIB_(G),SZJ_(G)), &
                   optional, intent(out) :: taux_bot !< Zonal bottom stress from ocean to
                                                     !! rock [R L Z T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)), &
                   optional, intent(out) :: tauy_bot !< Meridional bottom stress from ocean to
                                                     !! rock [R L Z T-2 ~> Pa]
  logical,         optional, intent(in)  :: fpmix !< fpmix along Eulerian shear
  type(wave_parameters_CS), &
                   optional, pointer     :: Waves !< Container for wave/Stokes information

  type(vertvisc_container_type)     :: visc_a   ! Copies of the viscosities and bottom drag in visc
  type(grid_core_type), pointer :: Gcore ! Persistent copies of the core ocean grid fields
  type(Box_t) :: bxU0 ! The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
  type(Box_t) :: bxV0 ! The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
  type(vgrid_consts_CS) :: vgrid ! The vertical grid and unit-scaling constants
  type(vertvisc_AD_type) :: AD_a ! Copies of the accelerations in the momentum equations for
                                 ! diagnostics
  type(limit_vel_report_type) :: report_a ! The velocities and flags used to report the columns in
                                          ! which velocity truncations occur
  integer :: ntrunc_sum ! The number of velocity truncations in this call [nondim]
  type(mech_forcing_container_type) :: forces_a ! Copies of the mechanical forcing fields in forces
  type(RealArray_t) :: u_a  ! Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t) :: v_a  ! Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t) :: h_a  ! Layer thickness [H ~> m or kg m-2]
  type(RealArray_t) :: taux_bot_a  ! Zonal bottom stress from ocean to rock [R L Z T-2 ~> Pa]
  type(RealArray_t) :: tauy_bot_a  ! Meridional bottom stress from ocean to rock [R L Z T-2 ~> Pa]
  type(RealArray_t) :: KE_u_a ! The area integral of a KE term in a layer at u-points
                              ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  type(RealArray_t) :: KE_v_a ! The area integral of a KE term in a layer at v-points
                              ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  logical :: lfpmix                ! fpmix along Eulerian shear, false if fpmix is absent
  logical :: DoStokesMixing
  type(RealArray_t) :: Us_x_a ! 3d zonal Stokes drift profile [L T-1 ~> m s-1]
  type(RealArray_t) :: Us_y_a ! 3d meridional Stokes drift profile [L T-1 ~> m s-1]

  if (.not.CS%initialized) call MOM_error(FATAL,"MOM_vert_friction(visc): "// &
         "Module must be initialized before it is used.")

  lfpmix = .false. ; if (present(fpmix)) lfpmix = fpmix
  !Check if Stokes mixing allowed if requested (present and associated)
  DoStokesMixing=.false.
  if (CS%vertvisc_opts_CS%StokesMixing) then
    if (present(Waves)) DoStokesMixing = associated(Waves)
    if (.not. DoStokesMixing) &
      call MOM_error(FATAL, "Stokes Mixing called without associated Waves Control Structure")
  endif
  if (DoStokesMixing) then
    call Us_x_a%alloc(lb=LBOUND(Waves%Us_x), ub=UBOUND(Waves%Us_x), source=Waves%Us_x)
    call Us_y_a%alloc(lb=LBOUND(Waves%Us_y), ub=UBOUND(Waves%Us_y), source=Waves%Us_y)
  endif
  Gcore => grid_core(G)
  call bxU0%safe_alloc(ndims=3) ; call bxV0%safe_alloc(ndims=3)
  call bxU0%set(idxS=[G%IscB,G%jsc,1], idxE=[G%IecB,G%jec,GV%ke])
  call bxV0%set(idxS=[G%isc,G%JscB,1], idxE=[G%iec,G%JecB,GV%ke])
  call vgrid_consts_build(vgrid, GV, US)
  call visc_a%build_from(visc)
  call forces_a%build_from(forces)
  call AD_a%build_from(ADp)
  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  if (present(taux_bot)) &
    call taux_bot_a%alloc(lb=LBOUND(taux_bot), ub=UBOUND(taux_bot), source=taux_bot)
  if (present(tauy_bot)) &
    call tauy_bot_a%alloc(lb=LBOUND(tauy_bot), ub=UBOUND(tauy_bot), source=tauy_bot)
  if (CS%vertvisc_opts_CS%do_report_u) then
    call report_a%u_old_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=0.0)
    call report_a%vel_report_u_a%alloc(lb=[G%IsdB,G%JsdB], ub=[G%IedB,G%JedB], source=0.0)
    call report_a%dowrite_u_a%alloc(lb=[G%IsdB,G%JsdB], ub=[G%IedB,G%JedB], source=.false.)
  endif
  if (CS%vertvisc_opts_CS%do_report_v) then
    call report_a%v_old_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=0.0)
    call report_a%vel_report_v_a%alloc(lb=[G%IsdB,G%JsdB], ub=[G%IedB,G%JedB], source=0.0)
    call report_a%dowrite_v_a%alloc(lb=[G%IsdB,G%JsdB], ub=[G%IedB,G%JedB], source=.false.)
  endif
  ntrunc_sum = 0
  if (CS%id_GLwork > 0) then
    call KE_u_a%alloc(lb=[G%IsdB,G%jsd,1], ub=[G%IedB,G%jed,GV%ke], source=0.0)
    call KE_v_a%alloc(lb=[G%isd,G%JsdB,1], ub=[G%ied,G%JedB,GV%ke], source=0.0)
  endif

  call vertvisc_TR(bxU0, bxV0, u_a, v_a, h_a, forces_a, visc_a, dt, OBC, AD_a, Gcore, vgrid, &
                   CS%vertvisc_opts_CS, CS%coupling_state, taux_bot_a, tauy_bot_a, KE_u_a, &
                   KE_v_a, report_a, ntrunc_sum, lfpmix, Us_x_a, Us_y_a)

  call AD_a%copy_back(ADp)
  CS%ntrunc = CS%ntrunc + ntrunc_sum
  if (CS%vertvisc_opts_CS%do_report_u .or. CS%vertvisc_opts_CS%do_report_v) &
    call vertvisc_limit_vel_report(report_a, h, ADp, CDp, dt, forces%taux, forces%tauy, G, GV, &
                                   US, CS)

  call vertvisc_diagnostics(KE_u_a, KE_v_a, taux_bot_a, tauy_bot_a, ADp, G, Gcore, GV, CS)

  call u_a%copy2F(u)
  call v_a%copy2F(v)
  call visc_a%copy_back(visc)
  call forces_a%free()
  call u_a%free()
  call v_a%free()
  call h_a%free()
  call KE_u_a%free()
  call KE_v_a%free()
  call Us_x_a%free() ; call Us_y_a%free()
  call report_a%u_old_a%free() ; call report_a%v_old_a%free()
  call report_a%vel_report_u_a%free() ; call report_a%vel_report_v_a%free()
  call report_a%dowrite_u_a%free() ; call report_a%dowrite_v_a%free()
  call vgrid%g_prime_a%free()
  call bxU0%free() ; call bxV0%free()
  if (present(taux_bot)) then
    call taux_bot_a%copy2F(taux_bot)
    call taux_bot_a%free()
  endif
  if (present(tauy_bot)) then
    call tauy_bot_a%copy2F(tauy_bot)
    call tauy_bot_a%free()
  endif

end subroutine vertvisc

!> The implementation of vertvisc (the root of its call tree).
subroutine vertvisc_TR(bxU0, bxV0, u_a, v_a, h_a, forces_a, visc_a, dt, OBC, AD_a, Gcore, vgrid, &
                       CS, coupling_state, taux_bot_a, tauy_bot_a, KE_u_a, KE_v_a, report_a, &
                       ntrunc_sum, fpmix, Us_x_a, Us_y_a)
  type(Box_t),             intent(in)    :: bxU0   !< The u-point box of the whole column,
                                                   !! [IscB:IecB, jsc:jec, 1:nz]
  type(Box_t),             intent(in)    :: bxV0   !< The v-point box of the whole column,
                                                   !! [isc:iec, JscB:JecB, 1:nz]
  type(grid_core_type),   intent(in)    :: Gcore  !< Persistent copies of the core ocean grid fields
  type(vgrid_consts_CS),   intent(in)    :: vgrid  !< The vertical grid and unit-scaling constants
  type(RealArray_t),       intent(inout) :: u_a    !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t),       intent(inout) :: v_a    !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t),       intent(in)    :: h_a    !< Layer thickness [H ~> m or kg m-2]
  type(mech_forcing_container_type), &
                         intent(in)      :: forces_a !< Copies of the driving mechanical forces
  type(vertvisc_container_type), &
                         intent(inout)   :: visc_a !< Copies of the viscosities and bottom drag
  real,                  intent(in)      :: dt     !< Time increment [T ~> s]
  type(ocean_OBC_type),  pointer         :: OBC    !< Open boundary condition structure
  type(vertvisc_AD_type), intent(inout)  :: AD_a   !< Copies of the accelerations in the momentum
                                                   !! equations for diagnostics
  type(vertvisc_opts_CS),  intent(in)    :: CS     !< Vertical viscosity control structure
  type(coupling_state_type), &
                           intent(in)    :: coupling_state !< The viscous coupling coefficients and
                                                           !! effective layer thicknesses at
                                                           !! velocity points
  type(RealArray_t),       intent(inout) :: taux_bot_a !< Zonal bottom stress from ocean to
                                                       !! rock [R L Z T-2 ~> Pa]
  type(RealArray_t),       intent(inout) :: tauy_bot_a !< Meridional bottom stress from ocean to
                                                       !! rock [R L Z T-2 ~> Pa]
  type(RealArray_t),       intent(inout) :: KE_u_a !< The area integral of a KE term in a layer at
                                                   !! u-points [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  type(RealArray_t),       intent(inout) :: KE_v_a !< The area integral of a KE term in a layer at
                                                   !! v-points [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  type(limit_vel_report_type), &
                           intent(inout) :: report_a !< The velocities and flags used to report the
                                                     !! columns in which velocity truncations occur
  integer,                 intent(inout) :: ntrunc_sum !< The number of velocity truncations, added
                                                       !! to the value on entry [nondim]
  logical,                   intent(in)  :: fpmix !< fpmix along Eulerian shear
  type(RealArray_t),       intent(in)    :: Us_x_a !< 3d zonal Stokes drift profile [L T-1 ~> m s-1]
                                                   !! Horizontal -> U points
                                                   !! Vertical -> Mid-points
  type(RealArray_t),       intent(in)    :: Us_y_a !< 3d meridional Stokes drift profile
                                                   !! [L T-1 ~> m s-1]
                                                   !! Horizontal -> V points
                                                   !! Vertical -> Mid-points

  ! Fields from forces used in this subroutine:
  !   taux: Zonal wind stress [R L Z T-2 ~> Pa].
  !   tauy: Meridional wind stress [R L Z T-2 ~> Pa].

  ! Local variables

  real :: b1
    ! A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
  real :: c1(vgrid%ke)
    ! A variable used by the tridiagonal solver [nondim].
  real :: d1
    ! d1=1-c1 is used by the tridiagonal solver [nondim].
  real :: Ray
    ! Ray is the Rayleigh-drag velocity [H T-1 ~> m s-1 or Pa s m-1]
  real :: b_denom_1
    ! The first term in the denominator of b1 [H ~> m or kg m-2].

  real :: Hmix             ! The mixed layer thickness over which stress
                           ! is applied with direct_stress [H ~> m or kg m-2].
  real :: I_Hmix           ! The inverse of Hmix [H-1 ~> m-1 or m2 kg-1].
  real :: Idt              ! The inverse of the time step [T-1 ~> s-1].
  real :: dt_Rho0          ! The time step divided by the mean density [T H Z-1 R-1 ~> s m3 kg-1 or s].
  real :: h_neglect        ! A thickness that is so small it is usually lost
                           ! in roundoff and can be neglected [H ~> m or kg m-2].

  real :: stress           !   The surface stress times the time step, divided
                           ! by the density [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real :: accel_underflow  ! An acceleration magnitude that is so small that values that are less
                           ! than this are diagnosed as 0 [L T-2 ~> m s-2].
  real :: zDS, h_avg       ! Temporary thickness variables used with direct_stress
                           ! [H ~> m or kg m-2]
  real :: hfr              ! Temporary ratio of thicknesses used with direct_stress [nondim]
  real :: surface_stress(u_a%lb(1):u_a%ub(1), v_a%lb(2):v_a%ub(2))
    ! The same as stress, unless the wind stress is applied as a body force
    ! [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(:,:,:), contiguous, pointer :: KE_u ! The area integral of a KE term in a layer at
                                                      ! u-points [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  real, dimension(:,:,:), contiguous, pointer :: KE_v ! The area integral of a KE term in a layer at
                                                      ! v-points [H L4 T-3 ~> m5 s-3 or kg m2 s-3]

  logical :: DoStokesMixing
  logical :: lfpmix

  real, dimension(:,:), contiguous, pointer :: taux, tauy ! The zonal and meridional wind stresses
                                                         ! [R L Z T-2 ~> Pa]
  real, dimension(:,:), contiguous, pointer :: taux_shelf, tauy_shelf ! The zonal and meridional
                                                         ! stresses on the ocean under shelves
                                                         ! [R Z L T-2 ~> Pa]
  real, dimension(:,:,:), contiguous, pointer :: Ray_u, Ray_v ! The Rayleigh drag velocities at
                                                         ! u- and v-points
                                                         ! [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(:,:), contiguous, pointer :: taux_bot, tauy_bot ! The zonal and meridional
                                                         ! bottom stresses from ocean to rock
                                                         ! [R L Z T-2 ~> Pa]

  real, dimension(:,:,:), contiguous, pointer :: u, v ! The zonal and meridional velocities
                                                    ! [L T-1 ~> m s-1]
  real, dimension(:,:,:), contiguous, pointer :: h    ! Layer thickness [H ~> m or kg m-2]

  real, dimension(:,:,:), contiguous, pointer :: a_u, a_v ! The u- and v-drag coefficients
                                                    ! across interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(:,:,:), contiguous, pointer :: h_u, h_v ! The effective layer thicknesses at
                                                    ! u- and v-points [H ~> m or kg m-2]
  real, dimension(:,:,:), contiguous, pointer :: a_u_gl90, a_v_gl90 ! The u- and v-drag
                                                    ! coefficients associated with GL90 across
                                                    ! interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(:,:), contiguous, pointer :: a1_shelf_u, a1_shelf_v ! The u- and v-momentum
                                                    ! coupling coefficients under ice shelves
                                                    ! [H T-1 ~> m s-1 or Pa s m-1]

  real, dimension(:,:), contiguous, pointer :: mask2dCu ! 0 for boundary points and 1 for ocean
                                                        ! points on the u grid [nondim]
  real, dimension(:,:), contiguous, pointer :: mask2dCv ! 0 for boundary points and 1 for ocean
                                                        ! points on the v grid [nondim]
  real, dimension(:,:), contiguous, pointer :: areaCu ! The areas of the u-grid cells [L2 ~> m2]
  real, dimension(:,:), contiguous, pointer :: areaCv ! The areas of the v-grid cells [L2 ~> m2]
  real, dimension(:,:,:), contiguous, pointer :: du_dt_visc, dv_dt_visc ! The zonal and meridional
                                                         ! accelerations due to vertical viscosity
                                                         ! [L T-2 ~> m s-2]
  real, dimension(:,:,:), contiguous, pointer :: du_dt_visc_gl90, dv_dt_visc_gl90 ! The zonal and
                                                         ! meridional accelerations due to GL90
                                                         ! vertical viscosity [L T-2 ~> m s-2]
  real, dimension(:,:,:), contiguous, pointer :: du_dt_str, dv_dt_str ! The zonal and meridional
                                                         ! accelerations due to the surface stress
                                                         ! [L T-2 ~> m s-2]
  real, dimension(:,:,:), contiguous, pointer :: Us_x ! 3d zonal Stokes drift profile
                                                    ! [L T-1 ~> m s-1]
  real, dimension(:,:,:), contiguous, pointer :: Us_y ! 3d meridional Stokes drift profile
                                                    ! [L T-1 ~> m s-1]
  integer :: i, j, k, nz, n
  nz = vgrid%ke

  call u_a%view(u)
  call v_a%view(v)
  call h_a%view(h)
  call coupling_state%a_u%view(a_u)
  call coupling_state%a_v%view(a_v)
  call coupling_state%h_u%view(h_u)
  call coupling_state%h_v%view(h_v)
  call coupling_state%a_u_gl90%view(a_u_gl90)
  call coupling_state%a_v_gl90%view(a_v_gl90)
  nullify(a1_shelf_u, a1_shelf_v)
  if (coupling_state%a1_shelf_u%associated()) call coupling_state%a1_shelf_u%view(a1_shelf_u)
  if (coupling_state%a1_shelf_v%associated()) call coupling_state%a1_shelf_v%view(a1_shelf_v)
  call Gcore%mask2dCu%view(mask2dCu)
  call Gcore%mask2dCv%view(mask2dCv)
  call Gcore%areaCu%view(areaCu)
  call Gcore%areaCv%view(areaCv)

  nullify(taux, tauy)
  if (forces_a%taux_a%associated()) call forces_a%taux_a%view(taux)
  if (forces_a%tauy_a%associated()) call forces_a%tauy_a%view(tauy)
  nullify(taux_shelf, tauy_shelf, Ray_u, Ray_v)
  if (visc_a%taux_shelf_a%associated()) call visc_a%taux_shelf_a%view(taux_shelf)
  if (visc_a%tauy_shelf_a%associated()) call visc_a%tauy_shelf_a%view(tauy_shelf)
  if (visc_a%Ray_u_a%associated()) call visc_a%Ray_u_a%view(Ray_u)
  if (visc_a%Ray_v_a%associated()) call visc_a%Ray_v_a%view(Ray_v)
  nullify(taux_bot, tauy_bot)
  if (taux_bot_a%associated()) call taux_bot_a%view(taux_bot)
  if (tauy_bot_a%associated()) call tauy_bot_a%view(tauy_bot)
  nullify(KE_u, KE_v)
  if (KE_u_a%associated()) call KE_u_a%view(KE_u)
  if (KE_v_a%associated()) call KE_v_a%view(KE_v)
  nullify(du_dt_visc, dv_dt_visc, du_dt_visc_gl90, &
          dv_dt_visc_gl90, du_dt_str, dv_dt_str)
  if (AD_a%du_dt_visc_a%associated()) call AD_a%du_dt_visc_a%view(du_dt_visc)
  if (AD_a%dv_dt_visc_a%associated()) call AD_a%dv_dt_visc_a%view(dv_dt_visc)
  if (AD_a%du_dt_visc_gl90_a%associated()) call AD_a%du_dt_visc_gl90_a%view(du_dt_visc_gl90)
  if (AD_a%dv_dt_visc_gl90_a%associated()) call AD_a%dv_dt_visc_gl90_a%view(dv_dt_visc_gl90)
  if (AD_a%du_dt_str_a%associated()) call AD_a%du_dt_str_a%view(du_dt_str)
  if (AD_a%dv_dt_str_a%associated()) call AD_a%dv_dt_str_a%view(dv_dt_str)

  if (CS%direct_stress) then
    Hmix = CS%Hmix_stress
    I_Hmix = 1.0 / Hmix
  endif
  dt_Rho0 = dt / vgrid%H_to_RZ
  h_neglect = vgrid%H_subroundoff
  Idt = 1.0 / dt

  accel_underflow = CS%vel_underflow * Idt

  nullify(Us_x, Us_y)
  if (Us_x_a%associated()) call Us_x_a%view(Us_x)
  if (Us_y_a%associated()) call Us_y_a%view(Us_y)
  DoStokesMixing = associated(Us_x)
  lfpmix = fpmix

  !   Update the zonal velocity component using a modification of a standard
  ! tridiagonal solver.

  ! WGL: Brandon Reichl says the following is obsolete. u(I,j,k) already
  ! includes Stokes.
  ! When mixing down Eulerian current + Stokes drift add before calling solver
  if (DoStokesMixing) then
    do concurrent (k=bxU0%idxS(3):bxU0%idxE(3), j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1), mask2dCu(I,j) > 0.)
      u(I,j,k) = u(I,j,k) + Us_x(I,j,k)
    enddo
  endif

  if (lfpmix) then
    do concurrent (k=bxU0%idxS(3):bxU0%idxE(3), j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1), mask2dCu(I,j) > 0.)
      u(I,j,k) = u(I,j,k) - Us_x(I,j,k)
    enddo
  endif

  if (associated(du_dt_visc)) then
    do concurrent (k=bxU0%idxS(3):bxU0%idxE(3), j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1))
      du_dt_visc(I,j,k) = u(I,j,k)
    enddo
  endif

  !$omp target enter data map(to: AD_a)
  !$omp target enter data map(alloc: surface_stress)
  !$omp target enter data map(alloc: dv_dt_str)
  !$omp target enter data map(alloc: du_dt_str)

  if (associated(du_dt_visc_gl90)) then
    do concurrent (k=bxU0%idxS(3):bxU0%idxE(3), j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1))
      du_dt_visc_gl90(I,j,k) = u(I,j,k)
    enddo
  endif

  if (associated(du_dt_str)) then
    do concurrent (k=bxU0%idxS(3):bxU0%idxE(3), j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1))
      du_dt_str(I,j,k) = 0.0
    enddo
  endif

  !   One option is to have the wind stress applied as a body force
  ! over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
  ! the wind stress is applied as a stress boundary condition.
  if (CS%direct_stress) then
    do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1), mask2dCu(i,j) > 0.0) &
        DO_LOCALITY(local(k, zDS, stress, h_avg, hfr))
      surface_stress(I,j) = 0.0
      zDS = 0.0
      stress = dt_Rho0 * taux(I,j)
      do k=1,nz
        h_avg = 0.5 * (h(i,j,k) + h(i+1,j,k)) + h_neglect
        hfr = 1.0 ; if ((zDS+h_avg) > Hmix) hfr = (Hmix - zDS) / h_avg
        u(I,j,k) = u(I,j,k) + I_Hmix * hfr * stress
        if (associated(du_dt_str)) du_dt_str(I,j,k) = (I_Hmix * hfr * stress) * Idt
        zDS = zDS + h_avg ; if (zDS >= Hmix) exit
      enddo
    enddo
  else
    do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1))
      surface_stress(I,j) = dt_Rho0 * (mask2dCu(I,j)*taux(I,j))
    enddo
  endif

  ! perform forward elimination on the tridiagonal system
  !
  ! denote the diagonal of the system as b_k, the subdiagonal as a_k
  ! and the superdiagonal as c_k. The right-hand side terms are d_k.
  !
  ! ignoring the Rayleigh drag contribution,
  ! we have a_k = -dt * a_u(k)
  !         b_k = h_u(k) + dt * (a_u(k) + a_u(k+1))
  !         c_k = -dt * a_u(k+1)
  !
  ! for forward elimination, we want to:
  ! calculate c'_k = - c_k                / (b_k + a_k c'_(k-1))
  ! and       d'_k = (d_k - a_k d'_(k-1)) / (b_k + a_k c'_(k-1))
  ! where c'_1 = c_1/b_1 and d'_1 = d_1/b_1
  !
  ! This form is mathematically equivalent to Thomas' tridiagonal matrix algorithm, but it
  ! does not suffer from the acute sensitivity to truncation errors of the Thomas algorithm
  ! because it involves no subtraction, as discussed by Schopf & Loughe, MWR, 1995.
  !
  ! b1 is the denominator term 1 / (b_k + a_k c'_(k-1))
  ! b_denom_1 is (b_k + a_k + c_k) - a_k(1 - c'_(k-1))
  !            = (b_k + c_k + c'_(k-1))
  ! this is done so that d1 = b1 * b_denom_1 = 1 - c'_(k-1)
  ! c1(k) is -c'_(k - 1)
  ! and the right-hand-side is destructively updated to be d'_k

  !$omp target enter data map(alloc: c1)

  do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                 I=bxU0%idxS(1):bxU0%idxE(1), mask2dCu(I,j) > 0.) &
      DO_LOCALITY(local(k, b1, d1, Ray, b_denom_1))
    Ray = 0.
    if (associated(Ray_u)) Ray = Ray_u(I,j,1)

    b_denom_1 = h_u(I,j,1) + dt * (Ray + a_u(I,j,1))
    b1 = 1. / (b_denom_1 + dt * a_u(I,j,2))
    d1 = b_denom_1 * b1
    u(I,j,1) = b1 * (h_u(I,j,1) * u(I,j,1) + surface_stress(I,j))

    if (associated(du_dt_str)) then
      du_dt_str(I,j,1) = b1 * (h_u(I,j,1) * du_dt_str(I,j,1) + surface_stress(I,j) * Idt)
    endif

    do k=2,nz
      if (associated(Ray_u)) Ray = Ray_u(I,j,k)

      c1(k) = dt * a_u(I,j,K) * b1
      b_denom_1 = h_u(I,j,k) + dt * (Ray + a_u(I,j,K) * d1)
      b1 = 1. / (b_denom_1 + dt * a_u(I,j,K+1))
      d1 = b_denom_1 * b1
      u(I,j,k) = (h_u(I,j,k) * u(I,j,k) + dt * a_u(I,j,K) * u(I,j,k-1)) * b1

      if (associated(du_dt_str)) then
        du_dt_str(I,j,k) = (h_u(I,j,k) * du_dt_str(I,j,k) &
            + dt * a_u(I,j,K) * du_dt_str(I,j,k-1)) * b1
      endif

      !### Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
      ! XXX: Check GPU behavior
      if (dt < 0) exit
    enddo

    if (associated(du_dt_str)) then
      if (abs(du_dt_str(I,j,nz)) < accel_underflow) &
        du_dt_str(I,j,nz) = 0.
    endif

    do k=nz-1,1,-1
      u(I,j,k) = u(I,j,k) + c1(k+1) * u(I,j,k+1)

      if (associated(du_dt_str)) then
        du_dt_str(I,j,k) = du_dt_str(I,j,k) + c1(k+1) * du_dt_str(I,j,k+1)

        if (abs(du_dt_str(I,j,k)) < accel_underflow) &
          du_dt_str(I,j,k) = 0.0
      endif
    enddo
  enddo

  ! compute vertical velocity tendency that arises from GL90 viscosity;
  ! follow tridiagonal solve method as above; to avoid corrupting u,
  ! use ADp%du_dt_visc_gl90 as a placeholder for updated u (due to GL90) until last do loop
  if ((CS%do_du_dt_visc_gl90) .or. (CS%do_GLwork)) then
    if (associated(du_dt_visc_gl90)) then
      do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                     I=bxU0%idxS(1):bxU0%idxE(1), mask2dCu(I,j) > 0.) &
          DO_LOCALITY(local(k, b1, d1, b_denom_1))
        b_denom_1 = h_u(I,j,1)  ! CS%a_u_gl90(I,j,1) is zero
        b1 = 1.0 / (b_denom_1 + dt * a_u_gl90(I,j,2))
        d1 = b_denom_1 * b1

        du_dt_visc_gl90(I,j,1) = b1 * (h_u(I,j,1) * du_dt_visc_gl90(I,j,1))

        do k=2,nz
          c1(k) = dt * a_u_gl90(I,j,K) * b1
          b_denom_1 = h_u(I,j,k) + dt * (a_u_gl90(I,j,K)*d1)
          b1 = 1.0 / (b_denom_1 + dt * a_u_gl90(I,j,K+1))
          d1 = b_denom_1 * b1

          du_dt_visc_gl90(I,j,k) = (h_u(I,j,k) * du_dt_visc_gl90(I,j,k) &
              + dt * a_u_gl90(I,j,K) * du_dt_visc_gl90(I,j,k-1)) * b1
        enddo

        ! back substitute to solve for new velocities, held by ADp%du_dt_visc_gl90
        do k=nz-1,1,-1
          du_dt_visc_gl90(I,j,k) = &
              du_dt_visc_gl90(I,j,k) + c1(k+1) * du_dt_visc_gl90(I,j,k+1)
        enddo

        do k=1,nz
          ! now fill ADp%du_dt_visc_gl90(I,j,k) with actual velocity tendency due to GL90;
          ! note that on RHS: ADp%du_dt_visc(I,j,k) holds the original velocity value u(I,j,k)
          ! and ADp%du_dt_visc_gl90(I,j,k) the updated velocity due to GL90
          du_dt_visc_gl90(I,j,k) = &
              (du_dt_visc_gl90(I,j,k) - du_dt_visc(I,j,k)) * Idt

          if (abs(du_dt_visc_gl90(I,j,k)) < accel_underflow) then
            du_dt_visc_gl90(I,j,k) = 0.0
          endif
        enddo

        ! to compute energetics, we need to multiply by u*h, where u is original velocity before
        ! velocity update; note that ADp%du_dt_visc(I,j,k) holds the original velocity value u(I,j,k)
        if (CS%do_GLwork) then
          do k=1,nz
            KE_u(I,j,k) = du_dt_visc(I,j,k) * h_u(I,j,k) * areaCu(I,j) * &
                          du_dt_visc_gl90(I,j,k)
          enddo
        endif
      enddo
    endif
  endif

  if (associated(du_dt_visc)) then
    do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1)) DO_LOCALITY(local(k))
      do k=1,nz
        du_dt_visc(I,j,k) = (u(I,j,k) - du_dt_visc(I,j,k)) * Idt

        if (abs(du_dt_visc(I,j,k)) < accel_underflow) &
          du_dt_visc(I,j,k) = 0.0
      enddo
    enddo
  endif

  if (associated(taux_shelf)) then
    do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1))
      taux_shelf(I,j) = -vgrid%H_to_RZ * a1_shelf_u(I,j) * u(I,j,1) ! - u_shelf?
    enddo
  endif

  if (taux_bot_a%associated()) then
    do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1))
      taux_bot(I,j) = vgrid%H_to_RZ * (u(I,j,nz) * a_u(I,j,nz+1))
    enddo

    if (associated(Ray_u)) then
      do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                     I=bxU0%idxS(1):bxU0%idxE(1)) DO_LOCALITY(local(k))
        do k=1,nz
          taux_bot(I,j) = taux_bot(I,j) + vgrid%H_to_RZ * (Ray_u(I,j,k) * u(I,j,k))
        enddo
      enddo
    endif
  endif

  ! When mixing down Eulerian current + Stokes drift subtract after calling solver
  if (DoStokesMixing) then
    do concurrent (k=bxU0%idxS(3):bxU0%idxE(3), j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1), mask2dCu(I,j) > 0.)
      u(I,j,k) = u(I,j,k) - Us_x(I,j,k)
    enddo
  endif

  if (lfpmix) then
    do concurrent (k=bxU0%idxS(3):bxU0%idxE(3), j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1), mask2dCu(I,j) > 0.)
      u(I,j,k) = u(I,j,k) + Us_x(I,j,k)
    enddo
  endif

  ! == Now work on the meridional velocity component.

  ! When mixing down Eulerian current + Stokes drift add before calling solver
  if (DoStokesMixing) then
    do concurrent (k=bxV0%idxS(3):bxV0%idxE(3), J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1), mask2dCv(i,J) > 0.)
      v(i,j,k) = v(i,j,k) + Us_y(i,j,k)
    enddo
  endif

  if (lfpmix) then
    do concurrent (k=bxV0%idxS(3):bxV0%idxE(3), J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1), mask2dCv(i,J) > 0.)
      v(i,j,k) = v(i,j,k) - Us_y(i,j,k)
    enddo
  endif

  if (associated(dv_dt_visc)) then
    do concurrent (k=bxV0%idxS(3):bxV0%idxE(3), J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1))
      dv_dt_visc(i,J,k) = v(i,J,k)
    enddo
  endif

  if (associated(dv_dt_visc_gl90)) then
    do concurrent (k=bxV0%idxS(3):bxV0%idxE(3), J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1))
      dv_dt_visc_gl90(i,J,k) = v(i,J,k)
    enddo
  endif

  if (associated(dv_dt_str)) then
    do concurrent (k=bxV0%idxS(3):bxV0%idxE(3), J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1))
      dv_dt_str(i,J,k) = 0.0
    enddo
  endif

  !   One option is to have the wind stress applied as a body force
  ! over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
  ! the wind stress is applied as a stress boundary condition.
  if (CS%direct_stress) then
    do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1), mask2dCv(i,J) > 0.0) &
        DO_LOCALITY(local(k, zDS, stress, h_avg, hfr))
      surface_stress(i,J) = 0.0
      zDS = 0.0
      stress = dt_Rho0 * tauy(i,J)
      do k=1,nz
        h_avg = 0.5 * (h(i,J,k) + h(i,J+1,k)) + h_neglect
        hfr = 1.0 ; if ((zDS+h_avg) > Hmix) hfr = (Hmix - zDS) / h_avg
        v(i,J,k) = v(i,J,k) + I_Hmix * hfr * stress
        if (associated(dv_dt_str)) dv_dt_str(i,J,k) = (I_Hmix * hfr * stress) * Idt
        zDS = zDS + h_avg ; if (zDS >= Hmix) exit
      enddo
    enddo
  else
    do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1))
      surface_stress(i,J) = dt_Rho0 * (mask2dCv(i,J) * tauy(i,J))
    enddo
  endif

  do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                 i=bxV0%idxS(1):bxV0%idxE(1), mask2dCv(i,J) > 0.) &
      DO_LOCALITY(local(k, b1, d1, Ray, b_denom_1))
    Ray = 0.
    if (associated(Ray_v)) Ray = Ray_v(i,J,1)

    b_denom_1 = h_v(i,J,1) + dt * (Ray + a_v(i,J,1))
    b1 = 1.0 / (b_denom_1 + dt*a_v(i,J,2))
    d1 = b_denom_1 * b1
    v(i,J,1) = b1 * (h_v(i,J,1) * v(i,J,1) + surface_stress(i,J))

    if (associated(dv_dt_str)) then
      dv_dt_str(i,J,1) = b1 * (h_v(i,J,1) * dv_dt_str(i,J,1) + surface_stress(i,J) * Idt)
    endif

    do k=2,nz
      if (associated(Ray_v)) Ray = Ray_v(i,J,k)

      c1(k) = dt * a_v(i,J,K) * b1
      b_denom_1 = h_v(i,J,k) + dt * (Ray + a_v(i,J,K) * d1)
      b1 = 1. / (b_denom_1 + dt * a_v(i,J,K+1))
      d1 = b_denom_1 * b1
      v(i,J,k) = (h_v(i,J,k) * v(i,J,k) + dt * a_v(i,J,K) * v(i,J,k-1)) * b1

      if (associated(dv_dt_str)) then
        dv_dt_str(i,J,k) = (h_v(i,J,k) * dv_dt_str(i,J,k) &
            + dt * a_v(i,J,K) * dv_dt_str(i,J,k-1)) * b1
      endif

      !### Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
      ! XXX: Check GPU behavior
      if (dt < 0) exit
    enddo

    if (associated(dv_dt_str)) then
      if (abs(dv_dt_str(i,J,nz)) < accel_underflow) &
        dv_dt_str(i,J,nz) = 0.0
    endif

    do k=nz-1,1,-1
      v(i,J,k) = v(i,J,k) + c1(k+1) * v(i,J,k+1)

      if (associated(dv_dt_str)) then
        dv_dt_str(i,J,k) = dv_dt_str(i,J,k) + c1(k+1) * dv_dt_str(i,J,k+1)

        if (abs(dv_dt_str(i,J,k)) < accel_underflow) &
          dv_dt_str(i,J,k) = 0.0
      endif
    enddo
  enddo

  ! compute vertical velocity tendency that arises from GL90 viscosity;
  ! follow tridiagonal solve method as above; to avoid corrupting v,
  ! use ADp%dv_dt_visc_gl90 as a placeholder for updated v (due to GL90) until last do loop
  if ((CS%do_dv_dt_visc_gl90) .or. (CS%do_GLwork)) then
    if (associated(dv_dt_visc_gl90)) then
      do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                     i=bxV0%idxS(1):bxV0%idxE(1), mask2dCv(i,J) > 0.) &
          DO_LOCALITY(local(k, b1, d1, b_denom_1))
        b_denom_1 = h_v(i,J,1)  ! CS%a_v_gl90(i,J,1) is zero
        b1 = 1.0 / (b_denom_1 + dt*a_v_gl90(i,J,2))
        d1 = b_denom_1 * b1
        dv_dt_visc_gl90(I,J,1) = b1 * (h_v(i,J,1) * dv_dt_visc_gl90(i,J,1))

        do k=2,nz
          c1(k) = dt * a_v_gl90(i,J,K) * b1
          b_denom_1 = h_v(i,J,k) + dt * (a_v_gl90(i,J,K) * d1)
          b1 = 1.0 / (b_denom_1 + dt * a_v_gl90(i,J,K+1))
          d1 = b_denom_1 * b1
          dv_dt_visc_gl90(i,J,k) = (h_v(i,J,k) * dv_dt_visc_gl90(i,J,k) &
              + dt * a_v_gl90(i,J,K) * dv_dt_visc_gl90(i,J,k-1)) * b1
        enddo

        ! back substitute to solve for new velocities, held by ADp%dv_dt_visc_gl90
        do k=nz-1,1,-1
          dv_dt_visc_gl90(i,J,k) = dv_dt_visc_gl90(i,J,k) + c1(k+1) * dv_dt_visc_gl90(i,J,k+1)
        enddo
      enddo

      do concurrent (k=bxV0%idxS(3):bxV0%idxE(3), J=bxV0%idxS(2):bxV0%idxE(2), &
                     i=bxV0%idxS(1):bxV0%idxE(1), mask2dCv(i,J) > 0.)
        ! now fill ADp%dv_dt_visc_gl90(i,J,k) with actual velocity tendency due to GL90;
        ! note that on RHS: ADp%dv_dt_visc(i,J,k) holds the original velocity value v(i,J,k)
        ! and ADp%dv_dt_visc_gl90(i,J,k) the updated velocity due to GL90
        dv_dt_visc_gl90(i,J,k) = (dv_dt_visc_gl90(i,J,k) - dv_dt_visc(i,J,k)) * Idt

        if (abs(dv_dt_visc_gl90(i,J,k)) < accel_underflow) &
          dv_dt_visc_gl90(i,J,k) = 0.0
      enddo

      ! to compute energetics, we need to multiply by v*h, where u is original velocity before
      ! velocity update; note that ADp%dv_dt_visc(I,j,k) holds the original velocity value v(i,J,k)
      if (CS%do_GLwork) then
        do concurrent (k=bxV0%idxS(3):bxV0%idxE(3), J=bxV0%idxS(2):bxV0%idxE(2), &
                       i=bxV0%idxS(1):bxV0%idxE(1), mask2dCv(i,J) > 0.)
          ! note that on RHS: ADp%dv_dt_visc(I,j,k) holds the original velocity value v(I,j,k)
          KE_v(I,j,k) = dv_dt_visc(i,J,k) * h_v(i,J,k) * areaCv(i,J) * &
                        dv_dt_visc_gl90(i,J,k)
        enddo
      endif
    endif
  endif

  if (associated(dv_dt_visc)) then
    do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1)) DO_LOCALITY(local(k))
      do k=1,nz
        dv_dt_visc(i,J,k) = (v(i,J,k) - dv_dt_visc(i,J,k))*Idt
        if (abs(dv_dt_visc(i,J,k)) < accel_underflow) dv_dt_visc(i,J,k) = 0.0
      enddo
    enddo
  endif

  if (associated(tauy_shelf)) then
    do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1))
      tauy_shelf(i,J) = -vgrid%H_to_RZ * a1_shelf_v(i,J) * v(i,J,1) ! - v_shelf?
    enddo
  endif

  ! JORGE TODO: this has to be malloced
  if (tauy_bot_a%associated()) then
    do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1))
      tauy_bot(i,J) = vgrid%H_to_RZ * (v(i,J,nz) * a_v(i,J,nz+1))
    enddo

    if (associated(Ray_v)) then
      do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                     i=bxV0%idxS(1):bxV0%idxE(1)) DO_LOCALITY(local(k))
        do k=1,nz
          tauy_bot(i,J) = tauy_bot(i,J) + vgrid%H_to_RZ * (Ray_v(i,J,k)*v(i,J,k))
        enddo
      enddo
    endif
  endif

  ! When mixing down Eulerian current + Stokes drift subtract after calling solver
  if (DoStokesMixing) then
    do concurrent (k=bxV0%idxS(3):bxV0%idxE(3), J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1), mask2dCv(i,J) > 0.)
      v(i,J,k) = v(i,J,k) - Us_y(i,J,k)
    enddo
  endif

  if (lfpmix) then
    do concurrent (k=bxV0%idxS(3):bxV0%idxE(3), J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1), mask2dCv(i,J) > 0.)
      v(i,J,k) = v(i,J,k) + Us_y(i,J,k)
    enddo
  endif

  call vertvisc_limit_vel(bxU0, bxV0, u_a, v_a, dt, Gcore, vgrid, CS, &
                          coupling_state, report_a, ntrunc_sum)

  !$omp target exit data map(delete: c1)

  !$omp target exit data map(from: du_dt_str, dv_dt_str)
  !$omp target exit data map(delete: AD_a)
  !$omp target exit data map(delete: surface_stress)

  ! Here the velocities associated with open boundary conditions are applied.
  if (associated(OBC)) then
    do n=1,OBC%number_of_segments
      if (OBC%segment(n)%specified) then
        if (OBC%segment(n)%is_N_or_S) then
          J = OBC%segment(n)%HI%JsdB
          do k=1,nz ; do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
            v(i,J,k) = OBC%segment(n)%normal_vel(i,J,k)
          enddo ; enddo
        elseif (OBC%segment(n)%is_E_or_W) then
          I = OBC%segment(n)%HI%IsdB
          do k=1,nz ; do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
            u(I,j,k) = OBC%segment(n)%normal_vel(I,j,k)
          enddo ; enddo
        endif
      endif
    enddo
  endif

end subroutine vertvisc_TR

!> Offers the diagnostics of the vertical viscosity calculations, including the KE source from
!! GL90 vertical viscosity.
subroutine vertvisc_diagnostics(KE_u_a, KE_v_a, taux_bot_a, tauy_bot_a, ADp, G, Gcore, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(grid_core_type),   intent(in)    :: Gcore  !< Persistent copies of the core ocean grid fields
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(RealArray_t),       intent(in)    :: KE_u_a !< The area integral of a KE term in a layer at
                                                   !! u-points [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  type(RealArray_t),       intent(in)    :: KE_v_a !< The area integral of a KE term in a layer at
                                                   !! v-points [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  type(RealArray_t),       intent(in)    :: taux_bot_a !< Zonal bottom stress from ocean to
                                                       !! rock [R L Z T-2 ~> Pa]
  type(RealArray_t),       intent(in)    :: tauy_bot_a !< Meridional bottom stress from ocean to
                                                       !! rock [R L Z T-2 ~> Pa]
  type(accel_diag_ptrs),   intent(in)    :: ADp    !< Accelerations in the momentum
                                                   !! equations for diagnostics
  type(vertvisc_CS),       intent(inout) :: CS     !< Vertical viscosity control structure

  real, allocatable, dimension(:,:,:) :: KE_term ! A term in the kinetic energy budget
                                                 ! [H L2 T-3 ~> m3 s-3 or W m-2]
  real, dimension(:,:,:), contiguous, pointer :: KE_u ! The area integral of a KE term in a layer at
                                                      ! u-points [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  real, dimension(:,:,:), contiguous, pointer :: KE_v ! The area integral of a KE term in a layer at
                                                      ! v-points [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  real, dimension(:,:), contiguous, pointer :: taux_bot, tauy_bot ! The zonal and meridional
                                                         ! bottom stresses from ocean to rock
                                                         ! [R L Z T-2 ~> Pa]
  real, dimension(:,:), contiguous, pointer :: IareaT ! 1/areaT [L-2 ~> m-2]
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  nullify(KE_u, KE_v, taux_bot, tauy_bot)
  if (KE_u_a%associated()) call KE_u_a%view(KE_u)
  if (KE_v_a%associated()) call KE_v_a%view(KE_v)
  if (taux_bot_a%associated()) call taux_bot_a%view(taux_bot)
  if (tauy_bot_a%associated()) call tauy_bot_a%view(tauy_bot)
  call Gcore%IareaT%view(IareaT)

  if (CS%id_GLwork > 0) then
    allocate(KE_term(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0)
    if (.not.G%symmetric) &
      call create_group_pass(CS%pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)
  endif

  ! Calculate the KE source from GL90 vertical viscosity [H L2 T-3 ~> m3 s-3].
  ! We do the KE-rate calculation here (rather than in MOM_diagnostics) to ensure
  ! a sign-definite term. MOM_diagnostics does not have access to the velocities
  ! and thicknesses used in the vertical solver, but rather uses a time-mean
  ! barotropic transport [uv]h.
  if (CS%id_GLwork > 0) then
    if (.not.G%symmetric) &
      call do_group_pass(CS%pass_KE_uv, G%domain)
    do concurrent (k=1:nz, j=js:je, i=is:ie)
      KE_term(i,j,k) = 0.5 * IareaT(i,j) &
          * (KE_u(I,j,k) + KE_u(I-1,j,k) + KE_v(i,J,k) + KE_v(i,J-1,k))
    enddo
    call post_data(CS%id_GLwork, KE_term, CS%diag)
  endif

  ! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_du_dt_visc > 0) &
      call post_data(CS%id_du_dt_visc, ADp%du_dt_visc, CS%diag)
    if (CS%id_du_dt_visc_gl90 > 0) &
      call post_data(CS%id_du_dt_visc_gl90, ADp%du_dt_visc_gl90, CS%diag)
    if (CS%id_dv_dt_visc > 0) &
      call post_data(CS%id_dv_dt_visc, ADp%dv_dt_visc, CS%diag)
    if (CS%id_dv_dt_visc_gl90 > 0) &
      call post_data(CS%id_dv_dt_visc_gl90, ADp%dv_dt_visc_gl90, CS%diag)
    if (taux_bot_a%associated() .and. (CS%id_taux_bot > 0)) &
      call post_data(CS%id_taux_bot, taux_bot, CS%diag)
    if (tauy_bot_a%associated() .and. (CS%id_tauy_bot > 0)) &
      call post_data(CS%id_tauy_bot, tauy_bot, CS%diag)
    if (CS%id_du_dt_str > 0) &
      call post_data(CS%id_du_dt_str, ADp%du_dt_str, CS%diag)
    if (CS%id_dv_dt_str > 0) &
      call post_data(CS%id_dv_dt_str, ADp%dv_dt_str, CS%diag)

    if (associated(ADp%du_dt_visc) .and. associated(ADp%dv_dt_visc)) then
      ! Diagnostics of the fractional thicknesses times momentum budget terms
      ! 3D diagnostics of hf_du(dv)_dt_visc are commented because there is no clarity on proper remapping grid option.
      ! The code is retained for debugging purposes in the future.
      !if (CS%id_hf_du_dt_visc > 0) &
      !  call post_product_u(CS%id_hf_du_dt_visc, ADp%du_dt_visc, ADp%diag_hfrac_u, G, nz, CS%diag)
      !if (CS%id_hf_dv_dt_visc > 0) &
      !  call post_product_v(CS%id_hf_dv_dt_visc, ADp%dv_dt_visc, ADp%diag_hfrac_v, G, nz, CS%diag)

      ! Diagnostics for thickness-weighted vertically averaged viscous accelerations
      if (CS%id_hf_du_dt_visc_2d > 0) &
        call post_product_sum_u(CS%id_hf_du_dt_visc_2d, ADp%du_dt_visc, ADp%diag_hfrac_u, G, nz, CS%diag)
      if (CS%id_hf_dv_dt_visc_2d > 0) &
        call post_product_sum_v(CS%id_hf_dv_dt_visc_2d, ADp%dv_dt_visc, ADp%diag_hfrac_v, G, nz, CS%diag)

      ! Diagnostics for thickness x viscous accelerations
      if (CS%id_h_du_dt_visc > 0) &
        call post_product_u(CS%id_h_du_dt_visc, ADp%du_dt_visc, ADp%diag_hu, G, nz, CS%diag)
      if (CS%id_h_dv_dt_visc > 0) &
        call post_product_v(CS%id_h_dv_dt_visc, ADp%dv_dt_visc, ADp%diag_hv, G, nz, CS%diag)
    endif

    if (associated(ADp%du_dt_str) .and.  associated(ADp%dv_dt_str)) then
      ! Diagnostics for thickness x wind stress accelerations
      if (CS%id_h_du_dt_str > 0) &
        call post_product_u(CS%id_h_du_dt_str, ADp%du_dt_str, ADp%diag_hu, G, nz, CS%diag)
      if (CS%id_h_dv_dt_str > 0) &
        call post_product_v(CS%id_h_dv_dt_str, ADp%dv_dt_str, ADp%diag_hv, G, nz, CS%diag)

      ! Diagnostics for wind stress accelerations multiplied by visc_rem_[uv],
      if (CS%id_du_dt_str_visc_rem > 0) &
        call post_product_u(CS%id_du_dt_str_visc_rem, ADp%du_dt_str, ADp%visc_rem_u, G, nz, CS%diag)
      if (CS%id_dv_dt_str_visc_rem > 0) &
        call post_product_v(CS%id_dv_dt_str_visc_rem, ADp%dv_dt_str, ADp%visc_rem_v, G, nz, CS%diag)
    endif
  endif

end subroutine vertvisc_diagnostics


!> Calculate the fraction of momentum originally in a layer that remains in the water column
!! after a time-step of viscosity, equivalently the fraction of a time-step's worth of
!! barotropic acceleration that a layer experiences after viscosity is applied.
subroutine vertvisc_remnant(visc, visc_rem_u, visc_rem_v, dt, G, GV, US, CS)
  type(ocean_grid_type), intent(in)   :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV   !< Ocean vertical grid structure
  type(vertvisc_type),   intent(in)   :: visc !< Viscosities and bottom drag
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                         intent(inout) :: visc_rem_u !< Fraction of a time-step's worth of a
                                              !! barotropic acceleration that a layer experiences after
                                              !! viscosity is applied in the zonal direction [nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                         intent(inout) :: visc_rem_v !< Fraction of a time-step's worth of a
                                              !! barotropic acceleration that a layer experiences after
                                              !! viscosity is applied in the meridional direction [nondim]
  real,                  intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  type(vertvisc_CS)                    :: CS  !< Vertical viscosity control structure

  type(vertvisc_container_type) :: visc_a ! Copies of the viscosities and bottom drag in visc
  type(grid_core_type), pointer :: Gcore ! Persistent copies of the core ocean grid fields
  type(Box_t) :: bxU0 ! The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
  type(Box_t) :: bxV0 ! The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
  type(RealArray_t) :: visc_rem_u_a ! Fraction of a time-step's worth of a barotropic
                                 ! acceleration that a layer experiences after viscosity
                                 ! is applied in the zonal direction [nondim]
  type(RealArray_t) :: visc_rem_v_a ! Fraction of a time-step's worth of a barotropic
                                 ! acceleration that a layer experiences after viscosity
                                 ! is applied in the meridional direction [nondim]

  if (.not.CS%initialized) call MOM_error(FATAL,"MOM_vert_friction(remnant): "// &
         "Module must be initialized before it is used.")

  Gcore => grid_core(G)
  call bxU0%safe_alloc(ndims=3) ; call bxV0%safe_alloc(ndims=3)
  call bxU0%set(idxS=[G%IscB,G%jsc,1], idxE=[G%IecB,G%jec,GV%ke])
  call bxV0%set(idxS=[G%isc,G%JscB,1], idxE=[G%iec,G%JecB,GV%ke])
  call visc_a%build_from(visc)
  call visc_rem_u_a%alloc(lb=LBOUND(visc_rem_u), ub=UBOUND(visc_rem_u), source=visc_rem_u)
  call visc_rem_v_a%alloc(lb=LBOUND(visc_rem_v), ub=UBOUND(visc_rem_v), source=visc_rem_v)

  call vertvisc_remnant_TR(bxU0, bxV0, visc_a, visc_rem_u_a, visc_rem_v_a, dt, Gcore, GV%ke, &
                           CS%coupling_state)

  call visc_rem_u_a%copy2F(visc_rem_u)
  call visc_rem_v_a%copy2F(visc_rem_v)
  call vertvisc_remnant_diagnostics(visc_rem_u, visc_rem_v, G, GV, CS)
  call visc_a%free()
  call visc_rem_u_a%free()
  call visc_rem_v_a%free()
  call bxU0%free() ; call bxV0%free()

end subroutine vertvisc_remnant

!> The implementation of vertvisc_remnant (the root of its call tree).
subroutine vertvisc_remnant_TR(bxU0, bxV0, visc_a, visc_rem_u_a, visc_rem_v_a, dt, Gcore, ke, &
                               coupling_state)
  type(Box_t),           intent(in)   :: bxU0 !< The u-point box of the whole column,
                                             !! [IscB:IecB, jsc:jec, 1:nz]
  type(Box_t),           intent(in)   :: bxV0 !< The v-point box of the whole column,
                                             !! [isc:iec, JscB:JecB, 1:nz]
  type(grid_core_type), intent(in)    :: Gcore !< Persistent copies of the core ocean grid fields
  integer,               intent(in)   :: ke   !< The number of layers/levels in the vertical
  type(vertvisc_container_type), &
                         intent(in)   :: visc_a !< Copies of the viscosities and bottom drag
  type(RealArray_t),     intent(inout) :: visc_rem_u_a !< Fraction of a time-step's worth of a
                                                       !! barotropic acceleration that a layer
                                                       !! experiences after viscosity is applied in
                                                       !! the zonal direction [nondim]
  type(RealArray_t),     intent(inout) :: visc_rem_v_a !< Fraction of a time-step's worth of a
                                                       !! barotropic acceleration that a layer
                                                       !! experiences after viscosity is applied in
                                                       !! the meridional direction [nondim]
  real,                  intent(in)    :: dt  !< Time increment [T ~> s]
  type(coupling_state_type), &
                         intent(in)    :: coupling_state !< The viscous coupling coefficients
                                                         !! and effective layer thicknesses at
                                                         !! velocity points

  ! Local variables

  real :: b1
    ! A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
  real :: c1(ke)
    ! A variable used by the tridiagonal solver [nondim].
  real :: d1
    ! d1=1-c1 is used by the tridiagonal solver [nondim].
  real :: Ray
    ! Ray is the Rayleigh-drag velocity [H T-1 ~> m s-1 or Pa s m-1]
  real :: b_denom_1
    ! The first term in the denominator of b1 [H ~> m or kg m-2].

  real, dimension(:,:,:), contiguous, pointer :: Ray_u, Ray_v ! The Rayleigh drag velocities at
                                                         ! u- and v-points
                                                         ! [H T-1 ~> m s-1 or Pa s m-1]

  real, dimension(:,:,:), contiguous, pointer :: visc_rem_u, visc_rem_v ! Fraction of a
                                     ! time-step's worth of a barotropic acceleration that a
                                     ! layer experiences after viscosity is applied in the
                                     ! zonal and meridional directions [nondim]

  real, dimension(:,:,:), contiguous, pointer :: a_u, a_v ! The u- and v-drag coefficients
                                                    ! across interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(:,:,:), contiguous, pointer :: h_u, h_v ! The effective layer thicknesses at
                                                    ! u- and v-points [H ~> m or kg m-2]

  real, dimension(:,:), contiguous, pointer :: mask2dCu ! 0 for boundary points and 1 for ocean
                                                        ! points on the u grid [nondim]
  real, dimension(:,:), contiguous, pointer :: mask2dCv ! 0 for boundary points and 1 for ocean
                                                        ! points on the v grid [nondim]
  integer :: i, j, k, nz
  nz = ke

  call visc_rem_u_a%view(visc_rem_u)
  call visc_rem_v_a%view(visc_rem_v)
  call coupling_state%a_u%view(a_u)
  call coupling_state%a_v%view(a_v)
  call coupling_state%h_u%view(h_u)
  call coupling_state%h_v%view(h_v)
  call Gcore%mask2dCu%view(mask2dCu)
  call Gcore%mask2dCv%view(mask2dCv)

  nullify(Ray_u, Ray_v)
  if (visc_a%Ray_u_a%associated()) call visc_a%Ray_u_a%view(Ray_u)
  if (visc_a%Ray_v_a%associated()) call visc_a%Ray_v_a%view(Ray_v)

  ! Find the zonal viscous remnant using a modification of a standard tridagonal solver.

  do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                 I=bxU0%idxS(1):bxU0%idxE(1), mask2dCu(I,j) > 0.) &
      DO_LOCALITY(local(k, b1, d1, Ray, b_denom_1))
    Ray = 0.
    if (associated(Ray_u)) Ray = Ray_u(I,j,1)

    b_denom_1 = h_u(I,j,1) + dt * (Ray + a_u(I,j,1))
    b1 = 1.0 / (b_denom_1 + dt * a_u(I,j,2))
    d1 = b_denom_1 * b1
    visc_rem_u(I,j,1) = b1 * h_u(I,j,1)

    do k=2,nz
      if (associated(Ray_u)) Ray = Ray_u(I,j,k)

      c1(k) = dt * a_u(I,j,K) * b1
      b_denom_1 = h_u(I,j,k) + dt * (Ray + a_u(I,j,K) * d1)
      b1 = 1.0 / (b_denom_1 + dt * a_u(I,j,K+1))
      d1 = b_denom_1 * b1
      visc_rem_u(I,j,k) = (h_u(I,j,k) + dt * a_u(I,j,K) * visc_rem_u(I,j,k-1)) * b1

      !### Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
      ! XXX: Check GPU behavior
      if (dt < 0) exit
    enddo

    do k=nz-1,1,-1
      visc_rem_u(I,j,k) = visc_rem_u(I,j,k) + c1(k+1) * visc_rem_u(I,j,k+1)
    enddo
  enddo

  ! Now find the meridional viscous remnant using the robust tridiagonal solver.

  do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                 i=bxV0%idxS(1):bxV0%idxE(1), mask2dCv(i,J) > 0.) &
      DO_LOCALITY(local(k, b1, d1, Ray, b_denom_1))
    Ray = 0.
    if (associated(Ray_v)) Ray = Ray_v(i,J,1)

    b_denom_1 = h_v(i,J,1) + dt * (Ray + a_v(i,J,1))
    b1 = 1.0 / (b_denom_1 + dt*a_v(i,J,2))
    d1 = b_denom_1 * b1
    visc_rem_v(i,J,1) = b1 * h_v(i,J,1)

    do k=2,nz
      if (associated(Ray_v)) Ray = Ray_v(i,J,k)

      c1(k) = dt * a_v(i,J,K) * b1
      b_denom_1 = h_v(i,J,k) + dt * (Ray + a_v(i,J,K) * d1)
      b1 = 1.0 / (b_denom_1 + dt * a_v(i,J,K+1))
      d1 = b_denom_1 * b1
      visc_rem_v(i,J,k) = (h_v(i,J,k) + dt * a_v(i,J,K) * visc_rem_v(i,J,k-1)) * b1

      !### Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
      ! XXX: Check GPU behavior
      if (dt < 0) exit
    enddo

    do k=nz-1,1,-1
      visc_rem_v(i,J,k) = visc_rem_v(i,J,k) + c1(k+1) * visc_rem_v(i,J,k+1)
    enddo
  enddo

end subroutine vertvisc_remnant_TR

!> Writes the debugging checksums of the fractions of a time-step's worth of a barotropic
!! acceleration that a layer experiences after viscosity is applied.
subroutine vertvisc_remnant_diagnostics(visc_rem_u, visc_rem_v, G, GV, CS)
  type(ocean_grid_type), intent(in)   :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV   !< Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                         intent(in)  :: visc_rem_u !< Fraction of a time-step's worth of a
                                                   !! barotropic acceleration that a layer
                                                   !! experiences after viscosity is applied in the
                                                   !! zonal direction [nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                         intent(in)  :: visc_rem_v !< Fraction of a time-step's worth of a
                                                   !! barotropic acceleration that a layer
                                                   !! experiences after viscosity is applied in the
                                                   !! meridional direction [nondim]
  type(vertvisc_CS),     intent(in)  :: CS  !< Vertical viscosity control structure

  if (CS%debug) then
    !$omp target update from(visc_rem_u, visc_rem_v)
    call uvchksum("visc_rem_[uv]", visc_rem_u, visc_rem_v, G%HI, haloshift=0, &
                  scalar_pair=.true.)
  endif
end subroutine vertvisc_remnant_diagnostics


!> Calculate the coupling coefficients (CS%a_u, CS%a_v, CS%a_u_gl90, CS%a_v_gl90)
!! and effective layer thicknesses (CS%h_u and CS%h_v) for later use in the
!! applying the implicit vertical viscosity via vertvisc().
subroutine vertvisc_coef(u, v, h, dz, forces, visc, tv, dt, G, GV, US, CS, OBC, VarMix)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dz     !< Vertical distance across layers [Z ~> m]
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),     intent(in)    :: visc   !< Viscosities and bottom drag
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any available
                                                   !! thermodynamic fields.
  real,                    intent(in)    :: dt     !< Time increment [T ~> s]
  type(vertvisc_CS),       intent(inout) :: CS     !< Vertical viscosity control structure
  type(ocean_OBC_type),    pointer       :: OBC    !< Open boundary condition structure
  type(VarMix_CS),         intent(in) :: VarMix !< Variable mixing coefficients

  type(RealArray_t) :: &
    Ustar_2d_a  ! The wind friction velocity, calculated using the Boussinesq reference density or
                ! the time-evolving surface density in non-Boussinesq mode [Z T-1 ~> m s-1]
  real, dimension(:,:), contiguous, pointer :: Ustar_2d ! The wind friction velocity
                                                        ! [Z T-1 ~> m s-1]
  type(vertvisc_container_type)     :: visc_a   ! Copies of the viscosities and bottom drag in visc
  type(grid_core_type), pointer :: Gcore ! Persistent copies of the core ocean grid fields
  type(Box_t) :: bxU0 ! The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
  type(Box_t) :: bxV0 ! The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
  type(vgrid_consts_CS) :: vgrid ! The vertical grid and unit-scaling constants
  type(mech_forcing_container_type) :: forces_a ! Copies of the mechanical forcing fields in forces
  type(RealArray_t) :: u_a  ! Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t) :: v_a  ! Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t) :: h_a  ! Layer thickness [H ~> m or kg m-2]
  type(RealArray_t) :: dz_a ! Vertical distance across layers [Z ~> m]
  type(vertvisc_coef_diag_type) :: coef_diag_a ! The diagnostic arrays filled by vertvisc_coef_TR
  type(RealArray_t) :: SpV_avg_a ! The layer averaged in situ specific volume [R-1 ~> m3 kg-1]
  type(RealArray_t) :: kdgl90_struct_a ! Vertical structure function used in GL90 diffusivity
                                       ! [nondim]

  if (.not.CS%initialized) call MOM_error(FATAL,"MOM_vert_friction(coef): "// &
         "Module must be initialized before it is used.")

  Gcore => grid_core(G)
  call bxU0%safe_alloc(ndims=3) ; call bxV0%safe_alloc(ndims=3)
  call bxU0%set(idxS=[G%IscB,G%jsc,1], idxE=[G%IecB,G%jec,GV%ke])
  call bxV0%set(idxS=[G%isc,G%JscB,1], idxE=[G%iec,G%JecB,GV%ke])
  call vgrid_consts_build(vgrid, GV, US)
  call visc_a%build_from(visc)
  call forces_a%build_from(forces)
  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call dz_a%alloc(lb=LBOUND(dz), ub=UBOUND(dz), source=dz)
  call Ustar_2d_a%alloc(lb=[G%isd,G%jsd], ub=[G%ied,G%jed])
  call Ustar_2d_a%view(Ustar_2d)
  if (allocated(tv%SpV_avg)) &
    call SpV_avg_a%alloc(lb=LBOUND(tv%SpV_avg), ub=UBOUND(tv%SpV_avg), source=tv%SpV_avg)
  if (VarMix%use_variable_mixing) then
    if (allocated(VarMix%kdgl90_struct)) &
      call kdgl90_struct_a%alloc(lb=LBOUND(VarMix%kdgl90_struct), ub=UBOUND(VarMix%kdgl90_struct), &
                                 source=VarMix%kdgl90_struct)
  endif
  if (CS%id_Kv_u > 0) call coef_diag_a%Kv_u_a%alloc(lb=[G%IsdB,G%jsd,1], &
                                               ub=[G%IedB,G%jed,GV%ke], source=0.0)
  if (CS%id_Kv_v > 0) call coef_diag_a%Kv_v_a%alloc(lb=[G%isd,G%JsdB,1], &
                                               ub=[G%ied,G%JedB,GV%ke], source=0.0)
  if (CS%id_Kv_gl90_u > 0) call coef_diag_a%Kv_gl90_u_a%alloc(lb=[G%IsdB,G%jsd,1], &
                                                         ub=[G%IedB,G%jed,GV%ke], source=0.0)
  if (CS%id_Kv_gl90_v > 0) call coef_diag_a%Kv_gl90_v_a%alloc(lb=[G%isd,G%JsdB,1], &
                                                         ub=[G%ied,G%JedB,GV%ke], source=0.0)
  if (CS%debug .or. (CS%id_hML_u > 0)) &
    call coef_diag_a%hML_u_a%alloc(lb=[G%IsdB,G%jsd], ub=[G%IedB,G%jed], source=0.0)
  if (CS%debug .or. (CS%id_hML_v > 0)) &
    call coef_diag_a%hML_v_a%alloc(lb=[G%isd,G%JsdB], ub=[G%ied,G%JedB], source=0.0)

  !$omp target enter data map(alloc: Ustar_2d)
  call find_ustar(forces, tv, Ustar_2d, G, GV, US, halo=1)

  call vertvisc_coef_TR(bxU0, bxV0, u_a, v_a, h_a, dz_a, forces_a, Ustar_2d_a, visc_a, &
                        SpV_avg_a, dt, Gcore, vgrid, CS%coef_opts_CS, CS%coupling_state, OBC, &
                        kdgl90_struct_a, coef_diag_a)
  !$omp target exit data map(delete: Ustar_2d)

  call vertvisc_coef_diagnostics(coef_diag_a, visc%Kv_slow, G, GV, US, CS)

  call visc_a%free()
  call forces_a%free()
  call u_a%free()
  call v_a%free()
  call h_a%free()
  call dz_a%free()
  call Ustar_2d_a%free()
  call SpV_avg_a%free()
  call kdgl90_struct_a%free()
  call coef_diag_a%hML_u_a%free()
  call coef_diag_a%hML_v_a%free()
  call coef_diag_a%Kv_u_a%free()
  call coef_diag_a%Kv_v_a%free()
  call coef_diag_a%Kv_gl90_u_a%free()
  call coef_diag_a%Kv_gl90_v_a%free()
  call vgrid%g_prime_a%free()
  call bxU0%free() ; call bxV0%free()

end subroutine vertvisc_coef

!> The implementation of vertvisc_coef (the root of its call tree).
subroutine vertvisc_coef_TR(bxU0, bxV0, u_a, v_a, h_a, dz_a, forces_a, Ustar_2d_a, visc_a, &
                            SpV_avg_a, dt, Gcore, vgrid, CS, coupling_state, OBC, kdgl90_struct_a, &
                            coef_diag_a)
  type(Box_t),             intent(in)    :: bxU0   !< The u-point box of the whole column,
                                                   !! [IscB:IecB, jsc:jec, 1:nz]
  type(Box_t),             intent(in)    :: bxV0   !< The v-point box of the whole column,
                                                   !! [isc:iec, JscB:JecB, 1:nz]
  type(grid_core_type),   intent(in)    :: Gcore  !< Persistent copies of the core ocean grid fields
  type(vgrid_consts_CS),   intent(in)    :: vgrid  !< The vertical grid and unit-scaling constants
  type(RealArray_t),       intent(in)    :: u_a    !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t),       intent(in)    :: v_a    !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t),       intent(in)    :: h_a    !< Layer thickness [H ~> m or kg m-2]
  type(RealArray_t),       intent(in)    :: dz_a   !< Vertical distance across layers [Z ~> m]
  type(mech_forcing_container_type), &
                           intent(in)    :: forces_a !< Copies of the driving mechanical forces
  type(RealArray_t),       intent(in)    :: Ustar_2d_a !< The wind friction velocity, calculated
                                                       !! using the Boussinesq reference density or
                                                       !! the time-evolving surface density in
                                                       !! non-Boussinesq mode [Z T-1 ~> m s-1]
  type(vertvisc_container_type), &
                           intent(in)    :: visc_a !< Copies of the viscosities and bottom drag
  type(RealArray_t),       intent(in)    :: SpV_avg_a !< The layer averaged in situ specific volume
                                                      !! [R-1 ~> m3 kg-1].
  real,                    intent(in)    :: dt     !< Time increment [T ~> s]
  type(coef_opts_CS),      intent(in)    :: CS     !< Vertical viscosity control structure
  type(coupling_state_type), &
                           intent(inout) :: coupling_state !< The viscous coupling coefficients and
                                                           !! effective layer thicknesses at
                                                           !! velocity points
  type(ocean_OBC_type),    pointer       :: OBC    !< Open boundary condition structure
  type(RealArray_t),       intent(in)    :: kdgl90_struct_a !< Vertical structure function used in
                                                            !! GL90 diffusivity [nondim]
  type(vertvisc_coef_diag_type), &
                           intent(inout) :: coef_diag_a !< The diagnostic arrays filled here, each
                                                        !! allocated only if it is needed
  ! Field from forces used in this subroutine:
  !   ustar: the friction velocity [Z T-1 ~> m s-1], used here as the mixing
  !     velocity in the mixed layer if NKML > 1 in a bulk mixed layer.

  ! Local variables

  real, dimension(vgrid%ke) :: &
    hvel, &     ! hvel is the thickness used at a velocity grid point [H ~> m or kg m-2].
    dz_harm, &  ! Harmonic mean of the vertical distances around a velocity grid point,
                ! given by 2*(h+ * h-)/(h+ + h-) [Z ~> m].
    dz_vel, &   ! The vertical distance between interfaces used at a velocity grid point [Z ~> m].
    hvel_shelf, & ! The equivalent of hvel under shelves [H ~> m or kg m-2].
    dz_vel_shelf ! The equivalent of dz_vel under shelves [Z ~> m].
  real :: &
    h_harm, &   ! Harmonic mean of the thicknesses around a velocity grid point,
                ! given by 2*(h+ * h-)/(h+ + h-) [H ~> m or kg m-2].
    h_arith, &  ! The arithmetic mean thickness [H ~> m or kg m-2].
    h_delta, &  ! The lateral difference of thickness [H ~> m or kg m-2].
    dz_arith    ! The arithmetic mean of the vertical distances around a velocity grid point [Z ~> m]
  real, dimension(vgrid%ke+1) :: &
    z_i, &      ! An estimate of each interface's height above the bottom,
                ! normalized by the bottom boundary layer thickness [nondim]
    z_i_gl90, & ! An estimate of each interface's height above the bottom,
                ! normalized by the GL90 bottom boundary layer thickness [nondim]
    a_cpl, &    ! The drag coefficients across interfaces [H T-1 ~> m s-1 or Pa s m-1].  a_cpl times
                ! the velocity difference gives the stress across an interface.
    a_cpl_gl90, & ! The drag coefficients across interfaces associated with GL90 [H T-1 ~> m s-1 or Pa s m-1].
                ! a_cpl_gl90 times the velocity difference gives the GL90 stress across an interface.
                ! a_cpl_gl90 is part of a_cpl.
    a_shelf     ! The drag coefficients across interfaces in water columns under
                ! ice shelves [H T-1 ~> m s-1 or Pa s m-1].
  real :: &
    kv_bbl, &     ! The bottom boundary layer viscosity [H Z T-1 ~> m2 s-1 or Pa s].
    bbl_thick, &  ! The bottom boundary layer thickness [Z ~> m].
    I_Hbbl, &     ! The inverse of the bottom boundary layer thickness [Z-1 ~> m-1].
    I_Hbbl_gl90, &! The inverse of the bottom boundary layer thickness used for the GL90 scheme
                  ! [Z-1 ~> m-1].
    I_HTbl, &     ! The inverse of the top boundary layer thickness [Z-1 ~> m-1].
    Ztop_min, &   ! The deeper of the two adjacent surface heights [Z ~> m].
    Dmin, &       ! The shallower of the two adjacent bottom depths [Z ~> m].
    zh, &         ! An estimate of the interface's distance from the bottom
                  ! based on harmonic mean thicknesses [Z ~> m].
    h_ml          ! The mixed layer depth [Z ~> m].
  real, dimension(:,:), contiguous, pointer :: hML_u ! Diagnostic of the mixed layer depth at u
                                                     ! points [Z ~> m].
  real, dimension(:,:), contiguous, pointer :: hML_v ! Diagnostic of the mixed layer depth at v
                                                     ! points [Z ~> m].
  real, dimension(:,:,:), contiguous, pointer :: Kv_u ! Total vertical viscosity at u-points in
                                                      ! thickness-based units
                                                      ! [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, dimension(:,:,:), contiguous, pointer :: Kv_v ! Total vertical viscosity at v-points in
                                                      ! thickness-based units
                                                      ! [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, dimension(:,:,:), contiguous, pointer :: Kv_gl90_u ! GL90 vertical viscosity at u-points in
                                                           ! thickness-based units
                                                           ! [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, dimension(:,:,:), contiguous, pointer :: Kv_gl90_v ! GL90 vertical viscosity at v-points in
                                                           ! thickness-based units
                                                           ! [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real :: zcol    ! The height of an interface at h-points [Z ~> m].
  real :: zcol_p1 ! An adjacent east/north h-point interface height [Z ~> m].
  real :: botfn   ! A function which goes from 1 at the bottom to 0 much more
                  ! than Hbbl into the interior [nondim].
  real :: topfn   ! A function which goes from 1 at the top to 0 much more
                  ! than Htbl into the interior [nondim].
  real :: z2      ! The distance from the bottom, normalized by Hbbl [nondim]
  real :: z2_wt   ! A nondimensional (0-1) weight used when calculating z2 [nondim].
  real :: z_clear ! The clearance of an interface above the surrounding topography [Z ~> m].
  real :: a_cpl_max  ! The maximum drag coefficient across interfaces, set so that it will be
                     ! representable as a 32-bit float in MKS units [H T-1 ~> m s-1 or Pa s m-1]
  real :: h_neglect  ! A thickness that is so small it is usually lost
                     ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: dz_neglect ! A vertical distance that is so small it is usually lost
                     ! in roundoff and can be neglected [Z ~> m].

  real :: I_valBL ! The inverse of a scaling factor determining when water is
                  ! still within the boundary layer, as determined by the sum
                  ! of the harmonic mean thicknesses [nondim].
  logical :: do_any_shelf
  integer :: zi_dir
    ! A ternary logical indicating which thickness to use for finding z_clear.
  real, dimension(:,:), contiguous, pointer :: frac_shelf_u, frac_shelf_v ! Fractional ice shelf
                                                         ! coverage of u- and v-cells [nondim]
  real, dimension(:,:), contiguous, pointer :: taux_shelf, tauy_shelf ! The zonal and meridional
                                                         ! stresses on the ocean under shelves
                                                         ! [R Z L T-2 ~> Pa]
  real, dimension(:,:), contiguous, pointer :: Kv_bbl_u, Kv_bbl_v ! The bottom boundary layer
                                                         ! viscosity at u- and v-points
                                                         ! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:), contiguous, pointer :: bbl_thick_u, bbl_thick_v ! The bottom boundary layer
                                                         ! thickness at u- and v-points [Z ~> m]
  real, dimension(:,:), contiguous, pointer :: tbl_thick_shelf_u, tbl_thick_shelf_v ! Thickness of
                                                         ! the viscous top boundary layer under ice
                                                         ! shelves at u- and v-points [Z ~> m]
  real, dimension(:,:), contiguous, pointer :: Kv_tbl_shelf_u, Kv_tbl_shelf_v ! Viscosity in the
                                                         ! viscous top boundary layer under ice
                                                         ! shelves at u- and v-points
                                                         ! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:), contiguous, pointer :: nkml_visc_u, nkml_visc_v ! The number of layers in
                                                         ! the viscous surface mixed layer at u- and
                                                         ! v-points [nondim]
  real, dimension(:,:,:), contiguous, pointer :: Kv_shear, Kv_shear_Bu ! The shear-driven turbulent
                                                         ! vertical viscosity at the interfaces in
                                                         ! tracer and corner columns
                                                         ! [H Z T-1 ~> m2 s-1 or Pa s]

  real, dimension(:,:,:), contiguous, pointer :: u, v ! The zonal and meridional velocities
                                                    ! [L T-1 ~> m s-1]
  real, dimension(:,:,:), contiguous, pointer :: h    ! Layer thickness [H ~> m or kg m-2]
  real, dimension(:,:,:), contiguous, pointer :: dz   ! Vertical distance across layers [Z ~> m]
  real, dimension(:,:), contiguous, pointer :: Ustar_2d ! The wind friction velocity
                                                       ! [Z T-1 ~> m s-1]

  real, dimension(:,:,:), contiguous, pointer :: a_u, a_v ! The u- and v-drag coefficients
                                                    ! across interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(:,:,:), contiguous, pointer :: h_u, h_v ! The effective layer thicknesses at
                                                    ! u- and v-points [H ~> m or kg m-2]
  real, dimension(:,:,:), contiguous, pointer :: a_u_gl90, a_v_gl90 ! The u- and v-drag
                                                    ! coefficients associated with GL90 across
                                                    ! interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(:,:), contiguous, pointer :: a1_shelf_u, a1_shelf_v ! The u- and v-momentum
                                                    ! coupling coefficients under ice shelves
                                                    ! [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(:,:), contiguous, pointer :: kappa_gl90_2d ! 2D kappa_gl90 at h-points
                                                    ! [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]

  real, dimension(:,:), contiguous, pointer :: mask2dCu ! 0 for boundary points and 1 for ocean
                                                        ! points on the u grid [nondim]
  real, dimension(:,:), contiguous, pointer :: mask2dCv ! 0 for boundary points and 1 for ocean
                                                        ! points on the v grid [nondim]
  real, dimension(:,:), contiguous, pointer :: bathyT ! Ocean bottom depth, referenced to Z_ref at
                                                      ! tracer points [Z ~> m]
  real, dimension(:,:), contiguous, pointer :: CoriolisBu ! The Coriolis parameter at corner points
                                                          ! [T-1 ~> s-1]
  real, dimension(:,:,:), contiguous, pointer :: SpV_avg ! The layer averaged in situ specific
                                                         ! volume [R-1 ~> m3 kg-1]
  real, dimension(:,:,:), contiguous, pointer :: kdgl90_struct ! Vertical structure function used
                                                         ! in GL90 diffusivity [nondim]
  real, dimension(:), contiguous, pointer :: g_prime ! The reduced gravity at each interface
                                                     ! [L2 Z-1 T-2 ~> m s-2].
  integer :: i, j, k, nz

  nz = vgrid%ke

  call u_a%view(u)
  call v_a%view(v)
  call h_a%view(h)
  call dz_a%view(dz)
  call Ustar_2d_a%view(Ustar_2d)
  call vgrid%g_prime_a%view(g_prime)
  nullify(SpV_avg, kdgl90_struct)
  if (SpV_avg_a%associated()) call SpV_avg_a%view(SpV_avg)
  if (kdgl90_struct_a%associated()) call kdgl90_struct_a%view(kdgl90_struct)
  call coupling_state%a_u%view(a_u)
  call coupling_state%a_v%view(a_v)
  call coupling_state%h_u%view(h_u)
  call coupling_state%h_v%view(h_v)
  call coupling_state%a_u_gl90%view(a_u_gl90)
  call coupling_state%a_v_gl90%view(a_v_gl90)
  nullify(kappa_gl90_2d)
  if (CS%GL90_CS%kappa_gl90_2d%associated()) call CS%GL90_CS%kappa_gl90_2d%view(kappa_gl90_2d)
  call Gcore%mask2dCu%view(mask2dCu)
  call Gcore%mask2dCv%view(mask2dCv)
  call Gcore%bathyT%view(bathyT)
  call Gcore%CoriolisBu%view(CoriolisBu)

  nullify(frac_shelf_u, frac_shelf_v)
  if (forces_a%frac_shelf_u_a%associated()) call forces_a%frac_shelf_u_a%view(frac_shelf_u)
  if (forces_a%frac_shelf_v_a%associated()) call forces_a%frac_shelf_v_a%view(frac_shelf_v)
  nullify(taux_shelf, tauy_shelf, Kv_bbl_u, Kv_bbl_v, bbl_thick_u, bbl_thick_v, tbl_thick_shelf_u, &
          tbl_thick_shelf_v, Kv_tbl_shelf_u, Kv_tbl_shelf_v, nkml_visc_u, nkml_visc_v, Kv_shear, &
          Kv_shear_Bu)
  if (visc_a%taux_shelf_a%associated()) call visc_a%taux_shelf_a%view(taux_shelf)
  if (visc_a%tauy_shelf_a%associated()) call visc_a%tauy_shelf_a%view(tauy_shelf)
  if (visc_a%Kv_bbl_u_a%associated()) call visc_a%Kv_bbl_u_a%view(Kv_bbl_u)
  if (visc_a%Kv_bbl_v_a%associated()) call visc_a%Kv_bbl_v_a%view(Kv_bbl_v)
  if (visc_a%bbl_thick_u_a%associated()) call visc_a%bbl_thick_u_a%view(bbl_thick_u)
  if (visc_a%bbl_thick_v_a%associated()) call visc_a%bbl_thick_v_a%view(bbl_thick_v)
  if (visc_a%tbl_thick_shelf_u_a%associated()) &
    call visc_a%tbl_thick_shelf_u_a%view(tbl_thick_shelf_u)
  if (visc_a%tbl_thick_shelf_v_a%associated()) &
    call visc_a%tbl_thick_shelf_v_a%view(tbl_thick_shelf_v)
  if (visc_a%Kv_tbl_shelf_u_a%associated()) call visc_a%Kv_tbl_shelf_u_a%view(Kv_tbl_shelf_u)
  if (visc_a%Kv_tbl_shelf_v_a%associated()) call visc_a%Kv_tbl_shelf_v_a%view(Kv_tbl_shelf_v)
  if (visc_a%nkml_visc_u_a%associated()) call visc_a%nkml_visc_u_a%view(nkml_visc_u)
  if (visc_a%nkml_visc_v_a%associated()) call visc_a%nkml_visc_v_a%view(nkml_visc_v)
  if (visc_a%Kv_shear_a%associated()) call visc_a%Kv_shear_a%view(Kv_shear)
  if (visc_a%Kv_shear_Bu_a%associated()) call visc_a%Kv_shear_Bu_a%view(Kv_shear_Bu)
  nullify(hML_u, hML_v, Kv_u, Kv_v, Kv_gl90_u, Kv_gl90_v)
  if (coef_diag_a%hML_u_a%associated()) call coef_diag_a%hML_u_a%view(hML_u)
  if (coef_diag_a%hML_v_a%associated()) call coef_diag_a%hML_v_a%view(hML_v)
  if (coef_diag_a%Kv_u_a%associated()) call coef_diag_a%Kv_u_a%view(Kv_u)
  if (coef_diag_a%Kv_v_a%associated()) call coef_diag_a%Kv_v_a%view(Kv_v)
  if (coef_diag_a%Kv_gl90_u_a%associated()) call coef_diag_a%Kv_gl90_u_a%view(Kv_gl90_u)
  if (coef_diag_a%Kv_gl90_v_a%associated()) call coef_diag_a%Kv_gl90_v_a%view(Kv_gl90_v)

  h_neglect = vgrid%H_subroundoff
  dz_neglect = vgrid%dZ_subroundoff
  a_cpl_max = 1.0e37 * vgrid%m_to_H * vgrid%T_to_s
  I_valBL = 0.0 ; if (CS%harm_BL_val > 0.0) I_valBL = 1.0 / CS%harm_BL_val

  if ((associated(taux_shelf) .or. associated(frac_shelf_u)) .and. &
      .not.coupling_state%a1_shelf_u%associated()) then
    call coupling_state%a1_shelf_u%alloc(lb=u_a%lb(1:2), ub=u_a%ub(1:2), source=0.0)
  endif
  if ((associated(tauy_shelf) .or. associated(frac_shelf_v)) .and. &
      .not.coupling_state%a1_shelf_v%associated()) then
    call coupling_state%a1_shelf_v%alloc(lb=v_a%lb(1:2), ub=v_a%ub(1:2), source=0.0)
  endif
  nullify(a1_shelf_u, a1_shelf_v)
  if (coupling_state%a1_shelf_u%associated()) call coupling_state%a1_shelf_u%view(a1_shelf_u)
  if (coupling_state%a1_shelf_v%associated()) call coupling_state%a1_shelf_v%view(a1_shelf_v)

  ! First do u-points

  ! TODO: tv and VarMix probably should be conditionally transferred (if at all)
  !$omp target enter data map(alloc: z_i, z_i_gl90, dz_harm, hvel, dz_vel, a_cpl, a_cpl_gl90, &
  !$omp& SpV_avg, kdgl90_struct, hvel_shelf, dz_vel_shelf, a_shelf)

  ! These are used in diagnostics, so they need to be mapped back and forth
  !$omp target enter data map(to: hML_u, kv_u, kv_gl90_u )
  !$omp target enter data map(to: hML_v, kv_v, kv_gl90_v)

  !$omp target update to(Kv_shear) if (associated(Kv_shear))
  !$omp target update to(Kv_shear_Bu) if (associated(Kv_shear_Bu))

  do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                 I=bxU0%idxS(1):bxU0%idxE(1), mask2dCu(I,j) > 0.) &
      DO_LOCALITY(local(k, I_Hbbl, I_Hbbl_gl90, kv_bbl, bbl_thick, Dmin, zi_dir, zh)) &
      DO_LOCALITY(local(zcol, zcol_p1, h_harm, h_arith, h_delta, dz_arith, z2, botfn)) &
      DO_LOCALITY(local(z_clear, z2_wt, h_ml, I_HTbl, Ztop_min, topfn, do_any_shelf))
    I_Hbbl = 1. / (CS%Hbbl + dz_neglect)
    if (CS%use_GL90_in_SSW) then
      I_Hbbl_gl90 = 1. / (CS%Hbbl_gl90 + dz_neglect)
    endif

    if (CS%coupling_CS%bottomdraglaw) then
      kv_bbl = Kv_bbl_u(I,j)
      bbl_thick = bbl_thick_u(I,j) + dz_neglect
      I_Hbbl = 1. / bbl_thick
    endif

    Dmin = min(bathyT(i,j), bathyT(i+1,j))
    zi_dir = 0

    ! Project thickness outward across OBCs using a zero-gradient condition.
    if (associated(OBC)) then
      if (OBC%u_E_OBCs_on_PE) then
        if (OBC%segnum_u(I,j) > 0) then
          Dmin = bathyT(i,j)
          zi_dir = -1
        endif
      endif

      if (OBC%u_W_OBCs_on_PE) then
        if (OBC%segnum_u(I,j) < 0) then
          Dmin = bathyT(i+1,j)
          zi_dir = 1
        endif
      endif
    endif

    ! The following block calculates the thicknesses at velocity grid points for
    ! the vertical viscosity (hvel and dz_vel).  Near the bottom an upwind biased
    ! thickness is used to control the effect of spurious Montgomery potential
    ! gradients at the bottom where nearly massless layers layers ride over the
    ! topography.

    z_i(nz+1) = 0.

    if (.not. CS%harmonic_visc) then
      zh = 0.
      zcol = -bathyT(i,j)
      zcol_p1 = -bathyT(i+1,j)
    endif

    if (CS%use_GL90_in_SSW) then
      z_i_gl90(nz+1) = 0.
    endif

    do k=nz,1,-1
      h_harm = 2. * h(i,j,k) * h(i+1,j,k) / (h(i,j,k) + h(i+1,j,k) + h_neglect)
      h_arith = 0.5 * (h(i+1,j,k) + h(i,j,k))
      h_delta = h(i+1,j,k) - h(i,j,k)
      dz_harm(k) = 2. * dz(i,j,k) * dz(i+1,j,k) / (dz(i,j,k) + dz(i+1,j,k) + dz_neglect)
      dz_arith = 0.5 * (dz(i+1,j,k) + dz(i,j,k))

      ! Project thickness outward across OBCs using a zero-gradient condition.
      if (associated(OBC)) then
        if (OBC%u_E_OBCs_on_PE) then
          if (OBC%segnum_u(I,j) > 0) then
            h_harm = h(i,j,k)
            h_arith = h(i,j,k)
            h_delta = 0.
            dz_harm(k) = dz(i,j,k)
            dz_arith = dz(i,j,k)
          endif
        endif

        if (OBC%u_W_OBCs_on_PE) then
          if (OBC%segnum_u(I,j) < 0) then
            h_harm = h(i+1,j,k)
            h_arith = h(i+1,j,k)
            h_delta = 0.
            dz_harm(k) = dz(i+1,j,k)
            dz_arith = dz(i+1,j,k)
          endif
        endif
      endif

      if (CS%harmonic_visc) then
        ! The following block calculates the thicknesses at velocity grid points
        ! for the vertical viscosity (hvel and dz_vel).  Near the bottom an
        ! upwind biased thickness is used to control the effect of spurious
        ! Montgomery potential gradients at the bottom where nearly massless
        ! layers ride over the topography.

        hvel(k) = h_harm
        dz_vel(k) = dz_harm(k)

        if (u(I,j,k) * h_delta < 0) then
          z2 = z_i(k+1)
          botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

          hvel(k) = (1. - botfn) * h_harm + botfn * h_arith
          dz_vel(k) = (1. - botfn) * dz_harm(k) + botfn * dz_arith
        endif

        z_i(k) =  z_i(k+1) + dz_harm(k) * I_Hbbl
      else
        zcol = zcol + dz(i,j,k)
        zcol_p1 = zcol_p1 + dz(i+1,j,k)

        zh = zh + dz_harm(k)

        z_clear = max(zcol, zcol_p1) + Dmin
        if (zi_dir < 0) z_clear = zcol + Dmin
        if (zi_dir > 0) z_clear = zcol_p1 + Dmin

        z_i(k) = max(zh, z_clear) * I_Hbbl

        hvel(k) = h_arith
        dz_vel(k) = dz_arith

        if (u(I,j,k) * h_delta > 0.) then
          if (zh * I_Hbbl < CS%harm_BL_val) then
            hvel(k) = h_harm
            dz_vel(k) = dz_harm(k)
          else
            z2_wt = 1.
            if (zh * I_Hbbl < 2. * CS%harm_BL_val) &
              z2_wt = max(0., min(1., zh * I_Hbbl * I_valBL - 1.))

            z2 = z2_wt * (max(zh, z_clear) * I_Hbbl)
            botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

            hvel(k) = (1. - botfn) * h_arith + botfn * h_harm
            dz_vel(k) = (1. - botfn) * dz_arith + botfn * dz_harm(k)
          endif
        endif
      endif

      if (CS%use_GL90_in_SSW) then
        ! The following block calculates the normalized height above the GL90 BBL
        ! (z_i_gl90), using a harmonic mean between layer thicknesses. For the
        ! GL90 BBL we use simply a constant (Hbbl_gl90). The purpose isthat the
        ! GL90 coupling coefficient is zeroed out within Hbbl_gl90, to ensure
        ! that no momentum gets fluxed into vanished layers. The scheme is not
        ! sensitive to the exact value of Hbbl_gl90, as long as it is in a
        ! reasonable range (~1-20 m): large enough to capture vanished layers
        ! over topography, small enough to not contaminate the interior.

        z_i_gl90(k) = z_i_gl90(k+1) + dz_harm(k) * I_Hbbl_gl90
      endif
    enddo

    call find_coupling_coef_k(a_cpl, dz_vel, i, j, dz_harm, bbl_thick, kv_bbl, z_i, h_ml, dt, &
        vgrid%ke, vgrid%Boussinesq, vgrid%nkml, vgrid%dZ_subroundoff, vgrid%H_to_m, vgrid%RZ_to_H, &
        vgrid%Z_to_H, CoriolisBu, CS%coupling_CS, Kv_shear, Kv_shear_Bu, Kv_tbl_shelf_u, &
        Kv_tbl_shelf_v, tbl_thick_shelf_u, tbl_thick_shelf_v, nkml_visc_u, nkml_visc_v, Ustar_2d, &
        SpV_avg, work_on_u=.true., OBC=OBC)

    if (associated(hML_u)) hML_u(I,j) = h_ml

    if (CS%use_GL90_in_SSW) then
      call find_coupling_coef_gl90(a_cpl_gl90, dz_vel, i, j, z_i_gl90, vgrid%ke, &
          vgrid%dZ_subroundoff, g_prime, CoriolisBu, CS%GL90_CS, kappa_gl90_2d, kdgl90_struct, &
          work_on_u=.true.)
    endif

    do_any_shelf = .false.
    if (associated(frac_shelf_u)) then
      a1_shelf_u(I,j) = 0.
      do_any_shelf = frac_shelf_u(I,j) > 0.

      if (do_any_shelf) then
        if (.not. CS%harmonic_visc) then
          zh = 0.
          Ztop_min = min(zcol, zcol_p1)
          I_HTbl = 1. / (tbl_thick_shelf_u(I,j) + dz_neglect)
        endif

        do k=1,nz
          if (CS%harmonic_visc) then
            hvel_shelf(k) = hvel(k)
            dz_vel_shelf(k) = dz_vel(k)
          else
            ! Find upwind-biased thickness near the surface.
            ! (Perhaps this needs to be done more carefully, via find_eta.)

            h_harm = 2. * h(i,j,k) * h(i+1,j,k) &
                / (h(i,j,k) + h(i+1,j,k) + h_neglect)
            h_arith = 0.5 * (h(i+1,j,k) + h(i,j,k))
            h_delta = h(i+1,j,k) - h(i,j,k)
            dz_arith = 0.5 * (dz(i+1,j,k) + dz(i,j,k))

            if (associated(OBC)) then
              if (OBC%u_E_OBCs_on_PE) then
                if (OBC%segnum_u(I,j) > 0) then
                  h_harm = h(i,j,k)
                  h_arith = h(i,j,k)
                  h_delta = 0.
                  dz_arith = dz(i,j,k)
                endif
              endif

              if (OBC%u_W_OBCs_on_PE) then
                if (OBC%segnum_u(I,j) < 0) then
                  h_harm = h(i+1,j,k)
                  h_arith = h(i+1,j,k)
                  h_delta = 0.
                  dz_arith = dz(i+1,j,k)
                endif
              endif
            endif

            zcol = zcol - dz(i,j,k)
            zcol_p1 = zcol_p1 - dz(i+1,j,k)

            zh = zh + dz_harm(k)

            hvel_shelf(k) = hvel(k)
            dz_vel_shelf(k) = dz_vel(k)

            if (u(I,j,k) * h_delta > 0.) then
              if (zh * I_HTbl < CS%harm_BL_val) then
                hvel_shelf(k) = min(hvel(k), h_harm)
                dz_vel_shelf(k) = min(dz_vel(k), dz_harm(k))
              else
                z2_wt = 1.
                if (zh * I_HTbl < 2. * CS%harm_BL_val) then
                  z2_wt = max(0., min(1., zh * I_HTbl * I_valBL - 1.))
                endif

                z2 = z2_wt * (max(zh, Ztop_min - min(zcol, zcol_p1)) * I_HTbl)
                ! TODO: replace **6 with multiply
                topfn = 1. / (1. + 0.09 * z2**6)

                hvel_shelf(k) = min(hvel(k), (1. - topfn) * h_arith + topfn * h_harm)
                dz_vel_shelf(k) = min(dz_vel(k), (1. - topfn) * dz_arith + topfn * dz_harm(k))
              endif
            endif
          endif
        enddo

        call find_coupling_coef(a_shelf, dz_vel_shelf, i, j, dz_harm, bbl_thick, kv_bbl, z_i, &
            h_ml, dt, vgrid%ke, vgrid%Boussinesq, vgrid%nkml, vgrid%dZ_subroundoff, vgrid%H_to_m, &
            vgrid%RZ_to_H, vgrid%Z_to_H, CoriolisBu, CS%coupling_CS, Kv_shear, Kv_shear_Bu, &
            Kv_tbl_shelf_u, Kv_tbl_shelf_v, tbl_thick_shelf_u, tbl_thick_shelf_v, nkml_visc_u, &
            nkml_visc_v, Ustar_2d, SpV_avg, work_on_u=.true., OBC=OBC, shelf=.true.)

        a1_shelf_u(I,j) = a_shelf(1)
      endif
    endif

    if (do_any_shelf) then
      if (CS%use_GL90_in_SSW) then
        do K=1,nz+1
          a_u(I,j,K) = min(a_cpl_max, (frac_shelf_u(I,j) * a_shelf(K) + &
                                         (1. - frac_shelf_u(I,j)) * a_cpl(K)) + a_cpl_gl90(K))

          ! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
          ! CS%a_u(I,j,K) = min(a_cpl_max, forces%frac_shelf_u(I,j) * max(a_shelf(K), a_cpl(K)) + &
          !                                (1. - forces%frac_shelf_u(I,j)) * a_cpl(K))

          a_u_gl90(I,j,K) = min(a_cpl_max, a_cpl_gl90(K))
        enddo
      else
        do K=1,nz+1
          a_u(I,j,K) = min(a_cpl_max, (frac_shelf_u(I,j) * a_shelf(K) + &
                                         (1. - frac_shelf_u(I,j)) * a_cpl(K)))

          ! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
          ! CS%a_u(I,j,K) = min(a_cpl_max, forces%frac_shelf_u(I,j) * max(a_shelf(K), a_cpl(K)) + &
          !                                (1. - forces%frac_shelf_u(I,j)) * a_cpl(K))
        enddo
      endif

      do k=1,nz
        ! Should we instead take the inverse of the average of the inverses?
        h_u(I,j,k) = frac_shelf_u(I,j) * hvel_shelf(k) &
            + (1. - frac_shelf_u(I,j)) * hvel(k) + h_neglect
      enddo
    else
      if (CS%use_GL90_in_SSW) then
        do K=1,nz+1
          a_cpl(K) = a_cpl(K) + a_cpl_gl90(K)
        enddo

        do K=1,nz+1
          a_u_gl90(I,j,K) = min(a_cpl_max, a_cpl_gl90(K))
        enddo
      endif

      do K=1,nz+1
        a_u(I,j,K) = min(a_cpl_max, a_cpl(K))
      enddo

      do k=1,nz
        h_u(I,j,k) = hvel(k) + h_neglect
      enddo
    endif

    ! Diagnose total Kv at u-points
    if (associated(Kv_u)) then
      do k=1,nz
        Kv_u(I,j,k) = 0.5 * (a_u(I,j,K) + a_u(I,j,K+1)) * h_u(I,j,k)
      enddo
    endif

    ! Diagnose GL90 Kv at u-points
    if (associated(Kv_gl90_u)) then
      do k=1,nz
        Kv_gl90_u(I,j,k) = 0.5 * (a_u_gl90(I,j,K) + a_u_gl90(I,j,K+1)) * h_u(I,j,k)
      enddo
    endif
  enddo

  ! Now work on v-points.

  do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                 i=bxV0%idxS(1):bxV0%idxE(1), mask2dCv(i,J) > 0.) &
      DO_LOCALITY(local(k, I_Hbbl, I_Hbbl_gl90, kv_bbl, bbl_thick, Dmin, zi_dir, zh)) &
      DO_LOCALITY(local(zcol, zcol_p1, h_harm, h_arith, h_delta, dz_arith, z2, botfn)) &
      DO_LOCALITY(local(z_clear, z2_wt, h_ml, I_HTbl, Ztop_min, topfn, do_any_shelf))
    I_Hbbl = 1. / (CS%Hbbl + dz_neglect)
    if (CS%use_GL90_in_SSW) then
      I_Hbbl_gl90 = 1. / (CS%Hbbl_gl90 + dz_neglect)
    endif

    if (CS%coupling_CS%bottomdraglaw) then
      kv_bbl = Kv_bbl_v(i,J)
      bbl_thick = bbl_thick_v(i,J) + dz_neglect
      I_Hbbl = 1. / bbl_thick
    endif

    Dmin = min(bathyT(i,j), bathyT(i,j+1))
    zi_dir = 0

    ! Project thickness outward across OBCs using a zero-gradient condition.
    if (associated(OBC)) then
      if (OBC%v_N_OBCs_on_PE) then
        if (OBC%segnum_v(i,J) > 0) then
          Dmin = bathyT(i,j)
          zi_dir = -1
        endif
      endif

      if (OBC%v_S_OBCs_on_PE) then
        if (OBC%segnum_v(i,J) < 0) then
          Dmin = bathyT(i,j+1)
          zi_dir = 1
        endif
      endif
    endif

    z_i(nz+1) = 0.

    if (.not. CS%harmonic_visc) then
      zh = 0.
      zcol = -bathyT(i,j)
      zcol_p1 = -bathyT(i,j+1)
    endif

    if (CS%use_GL90_in_SSW) then
      z_i_gl90(nz+1) = 0.
    endif

    do k=nz,1,-1
      h_harm = 2. * h(i,j,k) * h(i,j+1,k) / (h(i,j,k) + h(i,j+1,k) + h_neglect)
      h_arith = 0.5 * (h(i,j+1,k) + h(i,j,k))
      h_delta = h(i,j+1,k) - h(i,j,k)
      dz_harm(k) = 2. * dz(i,j,k) * dz(i,j+1,k) / (dz(i,j,k) + dz(i,j+1,k) + dz_neglect)
      dz_arith = 0.5 * (dz(i,j+1,k) + dz(i,j,k))

      ! Project thickness outward across OBCs using a zero-gradient condition.
      if (associated(OBC)) then
        if (OBC%v_N_OBCs_on_PE) then
          if (OBC%segnum_v(i,J) > 0) then
            h_harm = h(i,j,k)
            h_arith = h(i,j,k)
            h_delta = 0.
            dz_harm(k) = dz(i,j,k)
            dz_arith = dz(i,j,k)
          endif
        endif

        if (OBC%v_S_OBCs_on_PE) then
          if (OBC%segnum_v(i,J) < 0) then
            h_harm = h(i,j+1,k)
            h_arith = h(i,j+1,k)
            h_delta = 0.
            dz_harm(k) = dz(i,j+1,k)
            dz_arith = dz(i,j+1,k)
          endif
        endif
      endif

      if (CS%harmonic_visc) then
        ! The following block calculates the thicknesses at velocity grid points
        ! for the vertical viscosity (hvel and dz_vel).  Near the bottom an
        ! upwind biased thickness is used to control the effect of spurious
        ! Montgomery potential gradients at the bottom where nearly massless
        ! layers ride over the topography.

        hvel(k) = h_harm
        dz_vel(k) = dz_harm(k)

        if (v(i,J,k) * h_delta < 0) then
          z2 = z_i(k+1)
          botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

          hvel(k) = (1. - botfn) * h_harm + botfn * h_arith
          dz_vel(k) = (1. - botfn) * dz_harm(k) + botfn * dz_arith
        endif

        z_i(k) = z_i(k+1) + dz_harm(k) * I_Hbbl
      else
        zcol = zcol + dz(i,j,k)
        zcol_p1 = zcol_p1 + dz(i,j+1,k)

        zh = zh + dz_harm(k)

        z_clear = max(zcol, zcol_p1) + Dmin
        if (zi_dir < 0) z_clear = zcol + Dmin
        if (zi_dir > 0) z_clear = zcol_p1 + Dmin

        z_i(k) = max(zh, z_clear) * I_Hbbl

        hvel(k) = h_arith
        dz_vel(k) = dz_arith

        if (v(i,J,k) * h_delta > 0) then
          if (zh * I_Hbbl < CS%harm_BL_val) then
            hvel(k) = h_harm
            dz_vel(k) = dz_harm(k)
          else
            z2_wt = 1.
            if (zh * I_Hbbl < 2. * CS%harm_BL_val) &
              z2_wt = max(0., min(1., zh * I_Hbbl * I_valBL - 1.))

            ! TODO: should z_clear be used here?
            z2 = z2_wt * (max(zh, max(zcol, zcol_p1) + Dmin) * I_Hbbl)
            botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

            hvel(k) = (1. - botfn) * h_arith + botfn * h_harm
            dz_vel(k) = (1. - botfn) * dz_arith + botfn * dz_harm(k)
          endif
        endif
      endif

      if (CS%use_GL90_in_SSW) then
        ! The following block calculates the normalized height above the GL90 BBL
        ! (z_i_gl90), using a harmonic mean between layer thicknesses. For the
        ! GL90 BBL we use simply a constant (Hbbl_gl90). The purpose is that the
        ! GL90 coupling coefficient is zeroed out within Hbbl_gl90, to ensure
        ! that no momentum gets fluxed into vanished layers. The scheme is not
        ! sensitive to the exact value of Hbbl_gl90, as long as it is in a
        ! reasonable range (~1-20 m): large enough to capture vanished layers
        ! over topography, small enough to not contaminate the interior.

        z_i_gl90(k) = z_i_gl90(k+1) + dz_harm(k) * I_Hbbl_gl90
      endif
    enddo

    call find_coupling_coef_k(a_cpl, dz_vel, i, j, dz_harm, bbl_thick, kv_bbl, z_i, h_ml, dt, &
        vgrid%ke, vgrid%Boussinesq, vgrid%nkml, vgrid%dZ_subroundoff, vgrid%H_to_m, vgrid%RZ_to_H, &
        vgrid%Z_to_H, CoriolisBu, CS%coupling_CS, Kv_shear, Kv_shear_Bu, Kv_tbl_shelf_u, &
        Kv_tbl_shelf_v, tbl_thick_shelf_u, tbl_thick_shelf_v, nkml_visc_u, nkml_visc_v, Ustar_2d, &
        SpV_avg, work_on_u=.false., OBC=OBC)

    if (associated(hML_v)) hML_v(i,J) = h_ml

    if (CS%use_GL90_in_SSW) then
      call find_coupling_coef_gl90(a_cpl_gl90, dz_vel, i, j, z_i_gl90, vgrid%ke, &
          vgrid%dZ_subroundoff, g_prime, CoriolisBu, CS%GL90_CS, kappa_gl90_2d, kdgl90_struct, &
          work_on_u=.false.)
    endif

    do_any_shelf = .false.
    if (associated(frac_shelf_v)) then
      a1_shelf_v(i,J) = 0.
      do_any_shelf = frac_shelf_v(i,J) > 0.

      if (do_any_shelf) then
        ! Initialize non-harmonic depths
        if (.not. CS%harmonic_visc) then
          zh = 0.
          Ztop_min = min(zcol, zcol_p1)
          I_HTbl = 1. / (tbl_thick_shelf_v(i,J) + dz_neglect)
        endif

        do k=1,nz
          if (CS%harmonic_visc) then
            hvel_shelf(k) = hvel(k)
            dz_vel_shelf(k) = dz_vel(k)
          else
            ! Find upwind-biased thickness near the surface.
            ! Perhaps this needs to be done more carefully, via find_eta.
            h_harm = 2. * h(i,j,k) * h(i,j+1,k) &
                / (h(i,j,k) + h(i,j+1,k) + h_neglect)
            h_arith = 0.5 * (h(i,j+1,k) + h(i,j,k))
            h_delta = h(i,j+1,k) - h(i,j,k)
            dz_arith = 0.5 * (dz(i,j+1,k) + dz(i,j,k))

            ! Project thickness outward across OBCs using a zero-gradient condition.
            if (associated(OBC)) then
              if (OBC%v_N_OBCs_on_PE) then
                if (OBC%segnum_v(i,J) > 0) then
                  h_harm = h(i,j,k)
                  h_arith = h(i,j,k)
                  h_delta = 0.
                  dz_arith = dz(i,j,k)
                endif
              endif

              if (OBC%v_S_OBCs_on_PE) then
                if (OBC%segnum_v(i,J) < 0) then
                  h_harm = h(i,j+1,k)
                  h_arith = h(i,j+1,k)
                  h_delta = 0.
                  dz_arith = dz(i,j+1,k)
                endif
              endif
            endif

            zcol = zcol - dz(i,j,k)
            zcol_p1 = zcol_p1 - dz(i,j+1,k)

            zh = zh + dz_harm(k)

            hvel_shelf(k) = hvel(k)
            dz_vel_shelf(k) = dz_vel(k)

            if (v(i,J,k) * h_delta > 0.) then
              if (zh * I_HTbl < CS%harm_BL_val) then
                hvel_shelf(k) = min(hvel(k), h_harm)
                dz_vel_shelf(k) = min(dz_vel(k), dz_harm(k))
              else
                z2_wt = 1.
                if (zh * I_HTbl < 2. * CS%harm_BL_val) &
                  z2_wt = max(0., min(1., zh * I_HTbl * I_valBL - 1.))

                z2 = z2_wt * (max(zh, Ztop_min - min(zcol, zcol_p1)) * I_HTbl)
                ! TODO: Replace **6
                topfn = 1. / (1. + 0.09 * z2**6)

                hvel_shelf(k) = min(hvel(k), (1. - topfn) * h_arith + topfn * h_harm)
                dz_vel_shelf(k) = min(dz_vel(k), (1. - topfn) * dz_arith + topfn * dz_harm(k))
              endif
            endif
          endif
        enddo

        call find_coupling_coef(a_shelf, dz_vel_shelf, i, j, dz_harm, bbl_thick, kv_bbl, z_i, &
            h_ml, dt, vgrid%ke, vgrid%Boussinesq, vgrid%nkml, vgrid%dZ_subroundoff, vgrid%H_to_m, &
            vgrid%RZ_to_H, vgrid%Z_to_H, CoriolisBu, CS%coupling_CS, Kv_shear, Kv_shear_Bu, &
            Kv_tbl_shelf_u, Kv_tbl_shelf_v, tbl_thick_shelf_u, tbl_thick_shelf_v, nkml_visc_u, &
            nkml_visc_v, Ustar_2d, SpV_avg, work_on_u=.false., OBC=OBC, shelf=.true.)

        a1_shelf_v(i,J) = a_shelf(1)
      endif
    endif

    if (do_any_shelf) then
      if (CS%use_GL90_in_SSW) then
        do K=1,nz+1
          a_v(I,j,K) = min(a_cpl_max, (frac_shelf_v(I,j) * a_shelf(K) + &
                                         (1. - frac_shelf_v(I,j)) * a_cpl(K)) + a_cpl_gl90(K))

          ! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
          ! CS%a_v(I,j,K) = min(a_cpl_max, forces%frac_shelf_v(I,j) * max(a_shelf(K), a_cpl(K)) + &
          !                                (1. - forces%frac_shelf_v(I,j)) * a_cpl(K))

          a_v_gl90(I,j,K) = min(a_cpl_max, a_cpl_gl90(K))
        enddo
      else
        do K=1,nz+1
          a_v(I,j,K) = min(a_cpl_max, (frac_shelf_v(I,j) * a_shelf(K) + &
                                           (1. - frac_shelf_v(I,j)) * a_cpl(K)))
          ! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
          ! CS%a_v(I,j,K) = min(a_cpl_max, forces%frac_shelf_v(I,j) * max(a_shelf(K), a_cpl(K)) + &
          !                                (1. - forces%frac_shelf_v(I,j)) * a_cpl(K))
        enddo
      endif

      do k=1,nz
        ! Should we instead take the inverse of the average of the inverses?
        h_v(I,j,k) = frac_shelf_v(I,j) * hvel_shelf(k) &
            + (1. - frac_shelf_v(I,j)) * hvel(k) + h_neglect
      enddo
    else
      if (CS%use_GL90_in_SSW) then
        do K=1,nz+1
          a_cpl(K) = a_cpl(K) + a_cpl_gl90(K)
        enddo

        do K=1,nz+1
          a_v_gl90(i,J,K) = min(a_cpl_max, a_cpl_gl90(K))
        enddo
      endif

      do K=1,nz+1
        a_v(i,J,K) = min(a_cpl_max, a_cpl(K))
      enddo

      do k=1,nz
        h_v(i,J,k) = hvel(k) + h_neglect
      enddo
    endif

    ! Diagnose total Kv at v-points
    if (associated(Kv_v)) then
      do k=1,nz
        Kv_v(i,J,k) = 0.5 * (a_v(i,J,K) + a_v(i,J,K+1)) * h_v(i,J,k)
      enddo
    endif

    ! Diagnose GL90 Kv at v-points
    if (associated(Kv_gl90_v)) then
      do k=1,nz
        Kv_gl90_v(i,J,k) = 0.5 * (a_v_gl90(i,J,K) + a_v_gl90(i,J,K+1)) * h_v(i,J,k)
      enddo
    endif
  enddo

  !$omp target exit data map(delete: z_i, z_i_gl90, dz_harm, hvel, dz_vel, a_cpl, a_cpl_gl90, &
  !$omp& SpV_avg, kdgl90_struct, hvel_shelf, dz_vel_shelf, a_shelf, hml_u, kv_u, kv_gl90_u)

  ! These are used in diagnostics, so they need to be mapped back and forth
  !$omp target exit data map(from: hML_u, kv_u, kv_gl90_u )
  !$omp target exit data map(from: hML_v, kv_v, kv_gl90_v)

end subroutine vertvisc_coef_TR

!> Writes the debugging checksums and offers the diagnostics of the vertical viscosity coupling
!! coefficients and effective layer thicknesses.
subroutine vertvisc_coef_diagnostics(coef_diag_a, Kv_slow, G, GV, US, CS)
  type(ocean_grid_type),   intent(in) :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  type(vertvisc_coef_diag_type), &
                           intent(in) :: coef_diag_a !< The diagnostic arrays filled by the
                                                     !! coupling coefficient calculations
  real, dimension(:,:,:), &
                  pointer, intent(in) :: Kv_slow !< The turbulent vertical viscosity component due
                                                 !! to "slow" processes (e.g., tidal, background,
                                                 !! convection etc) [H Z T-1 ~> m2 s-1 or Pa s]
  type(vertvisc_CS),       intent(in) :: CS !< Vertical viscosity control structure

  real, dimension(:,:,:), contiguous, pointer :: a_u, a_v ! The u- and v-drag coefficients
                                                    ! across interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(:,:,:), contiguous, pointer :: h_u, h_v ! The effective layer thicknesses at
                                                    ! u- and v-points [H ~> m or kg m-2]
  real, dimension(:,:,:), contiguous, pointer :: a_u_gl90, a_v_gl90 ! The u- and v-drag
                                                    ! coefficients associated with GL90 across
                                                    ! interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(:,:), contiguous, pointer :: hML_u ! Diagnostic of the mixed layer depth at u
                                                     ! points [Z ~> m].
  real, dimension(:,:), contiguous, pointer :: hML_v ! Diagnostic of the mixed layer depth at v
                                                     ! points [Z ~> m].
  real, dimension(:,:,:), contiguous, pointer :: Kv_u ! Total vertical viscosity at u-points in
                                                      ! thickness-based units
                                                      ! [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, dimension(:,:,:), contiguous, pointer :: Kv_v ! Total vertical viscosity at v-points in
                                                      ! thickness-based units
                                                      ! [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, dimension(:,:,:), contiguous, pointer :: Kv_gl90_u ! GL90 vertical viscosity at u-points in
                                                           ! thickness-based units
                                                           ! [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, dimension(:,:,:), contiguous, pointer :: Kv_gl90_v ! GL90 vertical viscosity at v-points in
                                                           ! thickness-based units
                                                           ! [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].

  call CS%coupling_state%a_u%view(a_u)
  call CS%coupling_state%a_v%view(a_v)
  call CS%coupling_state%h_u%view(h_u)
  call CS%coupling_state%h_v%view(h_v)
  call CS%coupling_state%a_u_gl90%view(a_u_gl90)
  call CS%coupling_state%a_v_gl90%view(a_v_gl90)
  nullify(hML_u, hML_v, Kv_u, Kv_v, Kv_gl90_u, Kv_gl90_v)
  if (coef_diag_a%hML_u_a%associated()) call coef_diag_a%hML_u_a%view(hML_u)
  if (coef_diag_a%hML_v_a%associated()) call coef_diag_a%hML_v_a%view(hML_v)
  if (coef_diag_a%Kv_u_a%associated()) call coef_diag_a%Kv_u_a%view(Kv_u)
  if (coef_diag_a%Kv_v_a%associated()) call coef_diag_a%Kv_v_a%view(Kv_v)
  if (coef_diag_a%Kv_gl90_u_a%associated()) call coef_diag_a%Kv_gl90_u_a%view(Kv_gl90_u)
  if (coef_diag_a%Kv_gl90_v_a%associated()) call coef_diag_a%Kv_gl90_v_a%view(Kv_gl90_v)

  if (CS%debug) then
    !$omp target update from(h_u, h_v, a_u, a_v)
    call uvchksum("loop vertvisc_coef h_[uv]", h_u, h_v, G%HI, haloshift=0, &
                  unscale=GV%H_to_m, scalar_pair=.true.)
    call uvchksum("vertvisc_coef a_[uv]", a_u, a_v, G%HI, haloshift=0, &
                  unscale=GV%H_to_m*US%s_to_T, scalar_pair=.true.)
    if (associated(hML_u) .and. associated(hML_v)) &
      call uvchksum("vertvisc_coef hML_[uv]", hML_u, hML_v, G%HI, &
                    haloshift=0, unscale=US%Z_to_m, scalar_pair=.true.)
  endif

! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (associated(Kv_slow) .and. (CS%id_Kv_slow > 0)) &
        call post_data(CS%id_Kv_slow, Kv_slow, CS%diag)
    if (CS%id_Kv_u > 0) call post_data(CS%id_Kv_u, Kv_u, CS%diag)
    if (CS%id_Kv_v > 0) call post_data(CS%id_Kv_v, Kv_v, CS%diag)
    if (CS%id_Kv_gl90_u > 0) call post_data(CS%id_Kv_gl90_u, Kv_gl90_u, CS%diag)
    if (CS%id_Kv_gl90_v > 0) call post_data(CS%id_Kv_gl90_v, Kv_gl90_v, CS%diag)
    if (CS%id_au_vv > 0) call post_data(CS%id_au_vv, a_u, CS%diag)
    if (CS%id_av_vv > 0) call post_data(CS%id_av_vv, a_v, CS%diag)
    if (CS%id_au_gl90_vv > 0) call post_data(CS%id_au_gl90_vv, a_u_gl90, CS%diag)
    if (CS%id_av_gl90_vv > 0) call post_data(CS%id_av_gl90_vv, a_v_gl90, CS%diag)
    if (CS%id_h_u > 0) call post_data(CS%id_h_u, h_u, CS%diag)
    if (CS%id_h_v > 0) call post_data(CS%id_h_v, h_v, CS%diag)
    if (CS%id_hML_u > 0) call post_data(CS%id_hML_u, hML_u, CS%diag)
    if (CS%id_hML_v > 0) call post_data(CS%id_hML_v, hML_v, CS%diag)
  endif

end subroutine vertvisc_coef_diagnostics

!> Calculate the 'coupling coefficient' (a_cpl) at the interfaces.
!! If BOTTOMDRAGLAW is defined, the minimum of Hbbl and half the adjacent
!! layer thicknesses are used to calculate a_cpl near the bottom.
pure subroutine find_coupling_coef_k(a_cpl, hvel, i, j, h_harm, bbl_thick, kv_bbl, z_i, h_ml, dt, &
                                     ke, Boussinesq, nkml, dZ_subroundoff, H_to_m, RZ_to_H, &
                                     Z_to_H, CoriolisBu, CS, Kv_shear, Kv_shear_Bu, &
                                     Kv_tbl_shelf_u, Kv_tbl_shelf_v, tbl_thick_shelf_u, &
                                     tbl_thick_shelf_v, nkml_visc_u, nkml_visc_v, Ustar_2d, &
                                     SpV_avg, work_on_u, OBC, shelf)
  !$omp declare target
  integer,                   intent(in)  :: ke !< The number of layers/levels in the vertical
  logical,                   intent(in)  :: Boussinesq !< If true, make the Boussinesq
                                                       !! approximation.
  integer,                   intent(in)  :: nkml !< The number of layers at the top that should be
                                                 !! treated as parts of a homogeneous region.
  real,                      intent(in)  :: dZ_subroundoff !< A thickness in height units that is so
                                                           !! small that it can be added to a
                                                           !! vertical distance of Angstrom_Z or
                                                           !! 1e-17 m without changing it at the bit
                                                           !! level [Z ~> m].  This is the height
                                                           !! equivalent of H_subroundoff.
  real,                      intent(in)  :: H_to_m !< A constant that translates distances in the
                                                   !! units of thickness to m
                                                   !! [m H-1 ~> 1 or m3 kg-1].
  real,                      intent(in)  :: RZ_to_H !< A constant that translates mass per unit area
                                                    !! units to thickness units
                                                    !! [H R-1 Z-1 ~> m3 kg-2 or 1].
  real,                      intent(in)  :: Z_to_H !< A constant that translates depth units to
                                                   !! thickness units depth [H Z-1 ~> 1 or kg m-3].
  real, dimension(:,:), &
                   pointer, intent(in)  :: CoriolisBu !< The Coriolis parameter at corner points
                                                      !! [T-1 ~> s-1].
  real, dimension(ke+1), &
                             intent(out) :: a_cpl !< Coupling coefficient across interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(ke), &
                             intent(in)  :: hvel !< Distance between interfaces at velocity points [Z ~> m]
  integer,                   intent(in)  :: i    !< Column i-index
  integer,                   intent(in)  :: j    !< Column j-index
  real, dimension(ke), &
                             intent(in)  :: h_harm !< Harmonic mean of thicknesses around a velocity
                                                   !! grid point [Z ~> m]
  real, intent(in)  :: bbl_thick !< Bottom boundary layer thickness [Z ~> m]
  real, intent(in)  :: kv_bbl !< Bottom boundary layer viscosity, exclusive of
                                                   !! any depth-dependent contributions from
                                                   !! visc%Kv_shear [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(ke+1), intent(in) :: z_i !< Estimate of interface heights above the bottom,
                                                 !! normalized by the bottom boundary layer thickness [nondim]
  real,                      intent(out) :: h_ml !< Mixed layer depth [Z ~> m]
  real,                      intent(in)  :: dt   !< Time increment [T ~> s]
  type(coupling_CS),         intent(in)  :: CS   !< Vertical viscosity control structure
  real, dimension(:,:,:), &
                   pointer, intent(in)  :: Kv_shear !< The shear-driven turbulent vertical viscosity
                                                    !! at the interfaces between layers in tracer
                                                    !! columns [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:,:), &
                   pointer, intent(in)  :: Kv_shear_Bu !< The shear-driven turbulent vertical
                                                       !! viscosity at the interfaces between layers
                                                       !! in corner columns
                                                       !! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:), &
                   pointer, intent(in)  :: Kv_tbl_shelf_u !< Viscosity in the viscous top boundary
                                                          !! layer under ice shelves at u-points
                                                          !! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:), &
                   pointer, intent(in)  :: Kv_tbl_shelf_v !< Viscosity in the viscous top boundary
                                                          !! layer under ice shelves at v-points
                                                          !! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:), &
                   pointer, intent(in)  :: tbl_thick_shelf_u !< Thickness of the viscous top
                                                             !! boundary layer under ice shelves at
                                                             !! u-points [Z ~> m].
  real, dimension(:,:), &
                   pointer, intent(in)  :: tbl_thick_shelf_v !< Thickness of the viscous top
                                                             !! boundary layer under ice shelves at
                                                             !! v-points [Z ~> m].
  real, dimension(:,:), &
                   pointer, intent(in)  :: nkml_visc_u !< The number of layers in the viscous
                                                       !! surface mixed layer at u-points [nondim].
                                                       !! This is not an integer because there may
                                                       !! be fractional layers, and it is stored in
                                                       !! terms of layers, not depth, to facilitate
                                                       !! the movement of the viscous boundary layer
                                                       !! with the flow.
  real, dimension(:,:), &
                   pointer, intent(in)  :: nkml_visc_v !< The number of layers in the viscous
                                                       !! surface mixed layer at v-points [nondim].
  real, dimension(:,:), &
                   pointer, intent(in)  :: Ustar_2d !< The wind friction velocity, calculated using
                                                 !! the Boussinesq reference density or the
                                                 !! time-evolving surface density in non-Boussinesq
                                                 !! mode [Z T-1 ~> m s-1]
  real, dimension(:,:,:), &
                   pointer, intent(in)  :: SpV_avg !< The layer averaged in situ specific volume
                                                   !! [R-1 ~> m3 kg-1].
  logical,                   intent(in)  :: work_on_u !< If true, u-points are being calculated,
                                                  !! otherwise they are v-points
  type(ocean_OBC_type),      pointer     :: OBC   !< Open boundary condition structure
  logical,         optional, intent(in)  :: shelf !< If present and true, use a surface boundary
                                                  !! condition appropriate for an ice shelf.

  real :: u_star    ! ustar at a velocity point [Z T-1 ~> m s-1]
  real :: tau_mag   ! The magnitude of the wind stress at a velocity point including gustiness [H Z T-2 ~> m2 s-2 or Pa]
  real :: absf      ! The average of the neighboring absolute values of f [T-1 ~> s-1].
  real :: rho_av1   ! The harmonic mean surface layer density at velocity points [R ~> kg m-3]
  real :: z_t       ! The distance from the top, sometimes normalized
                    ! by Hmix, [Z ~> m] or [nondim].
  real :: kv_TBL    ! The viscosity in a top boundary layer under ice [H Z T-1 ~> m2 s-1 or Pa s]
  real :: tbl_thick ! The thickness of the top boundary layer [Z ~> m]
  real :: Kv_add    ! A viscosity to add [H Z T-1 ~> m2 s-1 or Pa s]
  real :: Kv_tot    ! The total viscosity at an interface [H Z T-1 ~> m2 s-1 or Pa s]
  integer :: nk_in_ml   ! The index of the deepest interface in the mixed layer.
  real :: h_shear ! The distance over which shears occur [Z ~> m].
  real :: dhc     ! The distance between the center of adjacent layers [Z ~> m].
  real :: visc_ml ! The mixed layer viscosity [H Z T-1 ~> m2 s-1 or Pa s].
  real :: I_Hmix  ! The inverse of the mixed layer thickness [Z-1 ~> m-1].
  real :: a_ml    ! The layer coupling coefficient across an interface in
                  ! the mixed layer [H T-1 ~> m s-1 or Pa s m-1].
  real :: a_floor ! A lower bound on the layer coupling coefficient across an interface in
                  ! the mixed layer [H T-1 ~> m s-1 or Pa s m-1].
  real :: I_amax  ! The inverse of the maximum coupling coefficient [T H-1 ~> s m-1 or s m2 kg-1].
  real :: temp1   ! A temporary variable [Z2 ~> m2]
  real :: ustar2_denom ! A temporary variable in the surface boundary layer turbulence
                  ! calculations [H Z-1 T-1 ~> s-1 or kg m-3 s-1]
  real :: h_neglect ! A vertical distance that is so small it is usually lost
                  ! in roundoff and can be neglected [Z ~> m].
  real :: z2      ! A copy of z_i [nondim]
  real :: botfn   ! A function that is 1 at the bottom and small far from it [nondim]
  real :: topfn   ! A function that is 1 at the top and small far from it [nondim]
  real :: kv_top  ! A viscosity associated with the top boundary layer [H Z T-1 ~> m2 s-1 or Pa s]
  logical :: do_shelf, do_OBCs, can_exit
  integer :: k
  integer :: nz, max_nk

  nz = ke

  h_neglect = dZ_subroundoff

  if (CS%answer_date < 20190101) then
    !   The maximum coupling coefficient was originally introduced to avoid
    ! truncation error problems in the tridiagonal solver. Effectively, the 1e-10
    ! sets the maximum coupling coefficient increment to 1e10 m per timestep.
    I_amax = (1.0e-10*H_to_m) * dt
  else
    I_amax = 0.0
  endif

  do_shelf = .false. ; if (present(shelf)) do_shelf = shelf

  do_OBCs = .false.
  !if (associated(OBC)) then
  !  if (work_on_u) then
  !    do_OBCS = OBC%u_E_OBCs_on_PE .or. OBC%u_W_OBCs_on_PE
  !  else
  !    do_OBCS = OBC%v_N_OBCs_on_PE .or. OBC%v_S_OBCs_on_PE
  !  endif
  !endif

  a_cpl(:) = 0.
  h_ml = 0.

  if (CS%Kvml_invZ2 > 0. .and. .not. do_shelf) then
    I_Hmix = 1. / (CS%Hmix + h_neglect)
    z_t = h_neglect * I_Hmix
  endif

  do K=2,nz
    Kv_tot = CS%Kv

    if (CS%Kvml_invZ2 > 0. .and. .not. do_shelf) then
      ! This is an older (vintage ~1997) way to prevent wind stresses from driving very
      ! large flows in nearly massless near-surface layers when there is not a physically-
      ! based surface boundary layer parameterization.  It does not have a plausible
      ! physical basis, and probably should not be used.
      z_t = z_t + h_harm(k-1) * I_Hmix
      Kv_tot = CS%Kv + CS%Kvml_invZ2 / ((z_t * z_t) *  &
               (1. + 0.09 * z_t * z_t * z_t * z_t * z_t * z_t))
    endif

    if (associated(Kv_shear)) then
      ! Add in viscosities that are determined by physical processes that are handled in
      ! other modules, and which do not respond immediately to the changing layer thicknesses.
      ! These processes may include shear-driven mixing or contributions from some boundary
      ! layer turbulence schemes.  Other viscosity contributions that respond to the evolving
      ! layer thicknesses or the surface wind stresses are added later.
      if (work_on_u) then
        Kv_add = 0.5 * (Kv_shear(i,j,k) + Kv_shear(i+1,j,k))

        if (do_OBCs) then
          if (OBC%u_E_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) > 0) then
              Kv_add = Kv_shear(i,j,k)
            endif
          endif

          if (OBC%u_W_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) < 0) then
              Kv_add = Kv_shear(i+1,j,k)
            endif
          endif
        endif

        Kv_tot = Kv_tot + Kv_add
      else
        Kv_add = 0.5 * (Kv_shear(i,j,k) + Kv_shear(i,j+1,k))

        if (do_OBCs) then
          if (OBC%v_N_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) > 0) then
              Kv_add = Kv_shear(i,j,k)
            endif
          endif

          if (OBC%v_S_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) < 0) then
              Kv_add = Kv_shear(i,j+1,k)
            endif
          endif
        endif

        Kv_tot = Kv_tot + Kv_add
      endif
    endif

    if (associated(Kv_shear_Bu)) then
      ! This is similar to what was done above, but for contributions coming from the corner
      ! (vorticity) points.  Because OBCs run through the faces and corners there is no need
      ! to further modify these viscosities here to take OBCs into account.
      if (work_on_u) then
        Kv_tot = Kv_tot + 0.5 * (Kv_shear_Bu(I,J-1,k) + Kv_shear_Bu(I,J,k))
      else
        Kv_tot = Kv_tot + 0.5 * (Kv_shear_Bu(I-1,J,k) + Kv_shear_Bu(I,J,k))
      endif
    endif

    ! Set the viscous coupling coefficients, excluding surface mixed layer contributions
    ! for now, but including viscous bottom drag, working up from the bottom.
    if (CS%bottomdraglaw) then
      !    botfn determines when a point is within the influence of the bottom
      !  boundary layer, going from 1 at the bottom to 0 in the interior.
      z2 = z_i(k)
      botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

      Kv_tot = Kv_tot + (kv_bbl - CS%Kv) * botfn
      dhc = 0.5 * (hvel(k) + hvel(k-1))
      if (dhc > bbl_thick) then
        h_shear = ((1. - botfn) * dhc + botfn * bbl_thick) + h_neglect
      else
        h_shear = dhc + h_neglect
      endif

      ! Calculate the coupling coefficients from the viscosities.
      a_cpl(K) = Kv_tot / (h_shear + (I_amax * Kv_tot))
    elseif (abs(CS%Kv_extra_bbl) > 0.0) then
      ! There is a simple enhancement of the near-bottom viscosities, but no
      ! adjustment of the viscous coupling length scales to give a particular
      ! bottom stress.

      !    botfn determines when a point is within the influence of the bottom
      !  boundary layer, going from 1 at the bottom to 0 in the interior.
      z2 = z_i(k)
      botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

      Kv_tot = Kv_tot + CS%Kv_extra_bbl * botfn
      h_shear = 0.5 * (hvel(k) + hvel(k-1) + h_neglect)

      ! Calculate the coupling coefficients from the viscosities.
      a_cpl(K) = Kv_tot / (h_shear + I_amax * Kv_tot)
    else
      ! Any near-bottom viscous enhancements were already incorporated into
      ! Kv_tot, and there is no adjustment of the viscous coupling length
      ! scales to give a particular bottom stress.

      h_shear = 0.5 * (hvel(k) + hvel(k-1) + h_neglect)
      ! Calculate the coupling coefficients from the viscosities.
      a_cpl(K) = Kv_tot / (h_shear + I_amax * Kv_tot)
    endif
  enddo

  ! Assign the bottom coupling coefficients
  if (CS%bottomdraglaw) then
    dhc = hvel(nz) * 0.5
    a_cpl(nz+1) = kv_bbl / ((min(dhc, bbl_thick) + h_neglect) + I_amax * kv_bbl)
  elseif (abs(CS%Kv_extra_bbl) > 0.0) then
    a_cpl(nz+1) = (CS%Kv + CS%Kv_extra_bbl) &
        / ((0.5 * hvel(nz) + h_neglect) + I_amax * (CS%Kv + CS%Kv_extra_bbl))
  else
    a_cpl(nz+1) = CS%Kv / ((0.5 * hvel(nz) + h_neglect) + I_amax * CS%Kv)
  endif

  ! Add surface intensified viscous coupling, either as a no-slip boundary condition under a
  ! rigid ice-shelf, or due to wind-stress driven surface boundary layer mixing that has not
  ! already been added via visc%Kv_shear.
  if (do_shelf) then
    ! Set the coefficients to include the no-slip surface stress.
    if (work_on_u) then
      kv_TBL = Kv_tbl_shelf_u(I,j)
      tbl_thick = tbl_thick_shelf_u(I,j) + h_neglect
    else
      kv_TBL = Kv_tbl_shelf_v(i,J)
      tbl_thick = tbl_thick_shelf_v(i,J) + h_neglect
    endif
    z_t = 0.0

    ! If a_cpl(1) were not already 0, it would be added here.
    if (0.5 * hvel(1) > tbl_thick) then
      a_cpl(1) = kv_TBL / (tbl_thick + I_amax * kv_TBL)
    else
      a_cpl(1) = kv_TBL / ((0.5 * hvel(1) + h_neglect) + I_amax * kv_TBL)
    endif

    do K=2,nz
      z_t = z_t + hvel(k-1) / tbl_thick
      topfn = 1. / (1. + 0.09 * z_t**6)

      dhc = 0.5 * (hvel(k) + hvel(k-1))
      if (dhc > tbl_thick) then
        h_shear = ((1. - topfn) * dhc + topfn * tbl_thick) + h_neglect
      else
        h_shear = dhc + h_neglect
      endif

      kv_top = topfn * kv_TBL
      a_cpl(K) = a_cpl(K) + kv_top / (h_shear + I_amax * kv_top)
    enddo
  elseif (CS%dynamic_viscous_ML .or. (nkml>0) .or. CS%fixed_LOTW_ML .or. CS%apply_LOTW_floor) then

    ! Find the friction velocity and the absolute value of the Coriolis parameter at this point.
    ! Zero out the friction velocity on land points.
    u_star = 0.
    tau_mag = 0.

    if (associated(SpV_avg)) then
      rho_av1 = 0.

      if (work_on_u) then
        u_star = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i+1,j))
        rho_av1 = 2. / (SpV_avg(i,j,1) + SpV_avg(i+1,j,1))
        absf = 0.5 * (abs(CoriolisBu(I,J-1)) + abs(CoriolisBu(I,J)))

        if (do_OBCs) then
          if (OBC%u_E_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) > 0) then
              u_star = Ustar_2d(i,j)
              rho_av1 = 1. / SpV_avg(i,j,1)
            endif
          endif

          if (OBC%u_W_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) < 0) then
              u_star = Ustar_2d(i+1,j)
              rho_av1 = 1. / SpV_avg(i+1,j,1)
            endif
          endif
        endif
      else
        u_star = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i,j+1))
        rho_av1 = 2. / (SpV_avg(i,j,1) + SpV_avg(i,j+1,1))
        absf = 0.5 * (abs(CoriolisBu(I-1,J)) + abs(CoriolisBu(I,J)))

        if (do_OBCs) then
          if (OBC%v_N_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) > 0) then
              u_star = Ustar_2d(i,j)
              rho_av1 = 1. / SpV_avg(i,j,1)
            endif
          endif

          if (OBC%v_S_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) < 0) then
              u_star = Ustar_2d(i,j+1)
              rho_av1 = 1. / SpV_avg(i,j+1,1)
            endif
          endif
        endif
      endif

      tau_mag = RZ_to_H * rho_av1 * u_star**2
    else ! (.not.allocated(tv%SpV_avg))
      if (work_on_u) then
        u_star = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i+1,j))
        absf = 0.5 * (abs(CoriolisBu(I,J-1)) + abs(CoriolisBu(I,J)))

        if (do_OBCs) then
          if (OBC%u_E_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) > 0) then
              u_star = Ustar_2d(i,j)
            endif
          endif

          if (OBC%u_W_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) < 0) then
              u_star = Ustar_2d(i+1,j)
            endif
          endif
        endif
      else
        u_star = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i,j+1))
        absf = 0.5 * (abs(CoriolisBu(I-1,J)) + abs(CoriolisBu(I,J)))

        if (do_OBCs) then
          if (OBC%v_N_OBCs_on_PE) then
              if (OBC%segnum_v(i,J) > 0) then
                u_star = Ustar_2d(i,j)
              endif
          endif

          if (OBC%v_S_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) < 0) then
              u_star = Ustar_2d(i,j+1)
            endif
          endif
        endif
      endif

      tau_mag = Z_to_H * u_star**2
    endif

    ! Determine the thickness of the surface ocean boundary layer and its extent in index space.
    nk_in_ml = 0

    if (CS%dynamic_viscous_ML) then
      ! The fractional number of layers that are within the viscous boundary layer were
      ! previously stored in visc%nkml_visc_[uv].
      h_ml = h_neglect
      max_nk = 0

      if (work_on_u) then
        nk_in_ml = ceiling(nkml_visc_u(I,j))
        max_nk = max(max_nk, nk_in_ml)

        do k=1,max_nk
          if (k <= nkml_visc_u(I,j)) then ! This layer is all in the ML.
            h_ml = h_ml + hvel(k)
          elseif (k < nkml_visc_u(I,j) + 1.) then ! Part of this layer is in the ML.
            h_ml = h_ml + ((nkml_visc_u(I,j) + 1.) - k) * hvel(k)
          endif
        enddo
      else
        nk_in_ml = ceiling(nkml_visc_v(i,J))
        max_nk = max(max_nk, nk_in_ml)

        do k=1,max_nk
          if (k <= nkml_visc_v(i,J)) then ! This layer is all in the ML.
            h_ml = h_ml + hvel(k)
          elseif (k < nkml_visc_v(i,J) + 1.) then ! Part of this layer is in the ML.
            h_ml = h_ml + ((nkml_visc_v(i,J) + 1.) - k) * hvel(k)
          endif
        enddo
      endif
    elseif (nkml>0) then
      ! This is a simple application of a refined-bulk mixed layer with GV%nkml sublayers.
      max_nk = nkml
      nk_in_ml = nkml

      h_ml = h_neglect

      do k=1,nkml
        h_ml = h_ml + hvel(k)
      enddo
    elseif (CS%fixed_LOTW_ML .or. CS%apply_LOTW_floor) then
      ! Determine which interfaces are within CS%Hmix of the surface, and set the viscous
      ! boundary layer thickness to the the smaller of CS%Hmix and the depth of the ocean.
      h_ml = 0.0
      do k=1,nz
        can_exit = .true.
        if (h_ml < CS%Hmix) then
          nk_in_ml = k

          if (h_ml + hvel(k) < CS%Hmix) then
            h_ml = h_ml + hvel(k)
            can_exit = .false.  ! Part of the next deeper layer is also in the mixed layer.
          else
            h_ml = CS%Hmix
          endif
        endif

        if (can_exit) exit  ! All remaining layers in this row are below the mixed layer depth.
      enddo

      ! This made more sense in the slab/layer solvers...
      max_nk = 0
      max_nk = max(max_nk, nk_in_ml)
    endif

    ! Avoid working on land or on columns where the viscous coupling could not be increased.
    if (u_star <= 0.) nk_in_ml = 0

    ! Set the viscous coupling at the interfaces as the larger of what was previously
    ! set and the contributions from the surface boundary layer.
    z_t = 0.0
    if (CS%apply_LOTW_floor .and. &
        (CS%dynamic_viscous_ML .or. nkml > 0 .or. CS%fixed_LOTW_ML)) then
      do K=2,max_nk
        if (k <= nk_in_ml) then
          z_t = z_t + hvel(k-1)

          !   The viscosity in visc_ml is set to go to 0 at the mixed layer top and bottom
          ! (in a log-layer) and be further limited by rotation to give the natural Ekman length.
          temp1 = (z_t * h_ml - z_t * z_t)
          if (Boussinesq) then
            ustar2_denom = (CS%vonKar * Z_to_H * u_star**2) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          else
            ustar2_denom = (CS%vonKar * tau_mag) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          endif

          visc_ml = temp1 * ustar2_denom
          ! Set the viscous coupling based on the model's vertical resolution.  The omission of
          ! the I_amax factor here is consistent with answer dates above 20190101.
          a_ml = visc_ml / (0.25 * (hvel(k) + hvel(k-1) + h_neglect))

          ! As a floor on the viscous coupling, assume that the length scale in the denominator can
          ! not be larger than the distance from the surface, consistent with a logarithmic velocity
          ! profile.  This is consistent with visc_ml, but cancels out common factors of z_t.
          a_floor = (h_ml - z_t) * ustar2_denom

          ! Choose the largest estimate of a_cpl.
          a_cpl(K) = max(a_cpl(K), a_ml, a_floor)
          ! An option could be added to change this to: a_cpl(i,K) = max(a_cpl(i,K) + a_ml, a_floor)
        endif
      enddo
    elseif (CS%apply_LOTW_floor) then
      do K=2,max_nk
        if (k <= nk_in_ml) then
          z_t = z_t + hvel(k-1)

          temp1 = (z_t * h_ml - z_t * z_t)
          if (Boussinesq) then
            ustar2_denom = (CS%vonKar * Z_to_H * u_star**2) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          else
            ustar2_denom = (CS%vonKar * tau_mag) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          endif

          ! As a floor on the viscous coupling, assume that the length scale in the denominator can not
          ! be larger than the distance from the surface, consistent with a logarithmic velocity profile.
          a_cpl(K) = max(a_cpl(K), (h_ml - z_t) * ustar2_denom)
        endif
      enddo
    else
      do K=2,max_nk
        if (k <= nk_in_ml) then
          z_t = z_t + hvel(k-1)

          temp1 = (z_t * h_ml - z_t * z_t)
          !   This viscosity is set to go to 0 at the mixed layer top and bottom (in a log-layer)
          ! and be further limited by rotation to give the natural Ekman length.
          ! The following expressions are mathematically equivalent.
          if (Boussinesq .or. (CS%answer_date < 20230601)) then
            visc_ml = u_star * CS%vonKar * (Z_to_H * temp1 * u_star) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          else
            visc_ml = CS%vonKar * (temp1 * tau_mag) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          endif
          a_ml = visc_ml / (0.25 * (hvel(k) + hvel(k-1) + h_neglect) + 0.5 * I_amax * visc_ml)

          ! Choose the largest estimate of a_cpl, but these could be changed to be additive.
          a_cpl(K) = max(a_cpl(K), a_ml)
          ! An option could be added to change this to: a_cpl(i,K) = a_cpl(i,K) + a_ml
        endif
      enddo
    endif
  endif
end subroutine find_coupling_coef_k


!> Calculate the 'coupling coefficient' (a_cpl) at the interfaces.
!! If BOTTOMDRAGLAW is defined, the minimum of Hbbl and half the adjacent
!! layer thicknesses are used to calculate a_cpl near the bottom.
pure subroutine find_coupling_coef(a_cpl, hvel, i, j, h_harm, bbl_thick, kv_bbl, z_i, h_ml, dt, &
                                   ke, Boussinesq, nkml, dZ_subroundoff, H_to_m, RZ_to_H, Z_to_H, &
                                   CoriolisBu, CS, Kv_shear, Kv_shear_Bu, Kv_tbl_shelf_u, &
                                   Kv_tbl_shelf_v, tbl_thick_shelf_u, tbl_thick_shelf_v, &
                                   nkml_visc_u, nkml_visc_v, Ustar_2d, SpV_avg, work_on_u, OBC, &
                                   shelf)
  !$omp declare target
  integer,                   intent(in)  :: ke !< The number of layers/levels in the vertical
  logical,                   intent(in)  :: Boussinesq !< If true, make the Boussinesq
                                                       !! approximation.
  integer,                   intent(in)  :: nkml !< The number of layers at the top that should be
                                                 !! treated as parts of a homogeneous region.
  real,                      intent(in)  :: dZ_subroundoff !< A thickness in height units that is so
                                                           !! small that it can be added to a
                                                           !! vertical distance of Angstrom_Z or
                                                           !! 1e-17 m without changing it at the bit
                                                           !! level [Z ~> m].  This is the height
                                                           !! equivalent of H_subroundoff.
  real,                      intent(in)  :: H_to_m !< A constant that translates distances in the
                                                   !! units of thickness to m
                                                   !! [m H-1 ~> 1 or m3 kg-1].
  real,                      intent(in)  :: RZ_to_H !< A constant that translates mass per unit area
                                                    !! units to thickness units
                                                    !! [H R-1 Z-1 ~> m3 kg-2 or 1].
  real,                      intent(in)  :: Z_to_H !< A constant that translates depth units to
                                                   !! thickness units depth [H Z-1 ~> 1 or kg m-3].
  real, dimension(:,:), &
                   pointer, intent(in)  :: CoriolisBu !< The Coriolis parameter at corner points
                                                      !! [T-1 ~> s-1].
  real, dimension(ke+1), &
                             intent(out) :: a_cpl !< Coupling coefficient across interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(ke), &
                             intent(in)  :: hvel !< Distance between interfaces at velocity points [Z ~> m]
  integer,                   intent(in)  :: i    !< Column i-index
  integer,                   intent(in)  :: j    !< Column j-index
  real, dimension(ke), &
                             intent(in)  :: h_harm !< Harmonic mean of thicknesses around a velocity
                                                   !! grid point [Z ~> m]
  real, intent(in)  :: bbl_thick !< Bottom boundary layer thickness [Z ~> m]
  real, intent(in)  :: kv_bbl !< Bottom boundary layer viscosity, exclusive of
                                                   !! any depth-dependent contributions from
                                                   !! visc%Kv_shear [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(ke+1), &
                             intent(in)  :: z_i  !< Estimate of interface heights above the bottom,
                                                 !! normalized by the bottom boundary layer thickness [nondim]
  real, intent(out) :: h_ml !< Mixed layer depth [Z ~> m]
  real,                      intent(in)  :: dt   !< Time increment [T ~> s]
  type(coupling_CS),         intent(in)  :: CS   !< Vertical viscosity control structure
  real, dimension(:,:,:), &
                   pointer, intent(in)  :: Kv_shear !< The shear-driven turbulent vertical viscosity
                                                    !! at the interfaces between layers in tracer
                                                    !! columns [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:,:), &
                   pointer, intent(in)  :: Kv_shear_Bu !< The shear-driven turbulent vertical
                                                       !! viscosity at the interfaces between layers
                                                       !! in corner columns
                                                       !! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:), &
                   pointer, intent(in)  :: Kv_tbl_shelf_u !< Viscosity in the viscous top boundary
                                                          !! layer under ice shelves at u-points
                                                          !! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:), &
                   pointer, intent(in)  :: Kv_tbl_shelf_v !< Viscosity in the viscous top boundary
                                                          !! layer under ice shelves at v-points
                                                          !! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:), &
                   pointer, intent(in)  :: tbl_thick_shelf_u !< Thickness of the viscous top
                                                             !! boundary layer under ice shelves at
                                                             !! u-points [Z ~> m].
  real, dimension(:,:), &
                   pointer, intent(in)  :: tbl_thick_shelf_v !< Thickness of the viscous top
                                                             !! boundary layer under ice shelves at
                                                             !! v-points [Z ~> m].
  real, dimension(:,:), &
                   pointer, intent(in)  :: nkml_visc_u !< The number of layers in the viscous
                                                       !! surface mixed layer at u-points [nondim].
                                                       !! This is not an integer because there may
                                                       !! be fractional layers, and it is stored in
                                                       !! terms of layers, not depth, to facilitate
                                                       !! the movement of the viscous boundary layer
                                                       !! with the flow.
  real, dimension(:,:), &
                   pointer, intent(in)  :: nkml_visc_v !< The number of layers in the viscous
                                                       !! surface mixed layer at v-points [nondim].
  real, dimension(:,:), &
                   pointer, intent(in)  :: Ustar_2d !< The wind friction velocity, calculated using
                                                 !! the Boussinesq reference density or the
                                                 !! time-evolving surface density in non-Boussinesq
                                                 !! mode [Z T-1 ~> m s-1]
  real, dimension(:,:,:), &
                   pointer, intent(in)  :: SpV_avg !< The layer averaged in situ specific volume
                                                   !! [R-1 ~> m3 kg-1].
  logical,                   intent(in)  :: work_on_u !< If true, u-points are being calculated,
                                                  !! otherwise they are v-points
  type(ocean_OBC_type),      pointer     :: OBC   !< Open boundary condition structure
  logical,         optional, intent(in)  :: shelf !< If present and true, use a surface boundary
                                                  !! condition appropriate for an ice shelf.

  ! Local variables

  real :: &
    u_star, &   ! ustar at a velocity point [Z T-1 ~> m s-1]
    tau_mag, &  ! The magnitude of the wind stress at a velocity point including gustiness [H Z T-2 ~> m2 s-2 or Pa]
    absf, &     ! The average of the neighboring absolute values of f [T-1 ~> s-1].
    rho_av1, &  ! The harmonic mean surface layer density at velocity points [R ~> kg m-3]
    z_t, &      ! The distance from the top, sometimes normalized
                ! by Hmix, [Z ~> m] or [nondim].
    kv_TBL, &   ! The viscosity in a top boundary layer under ice [H Z T-1 ~> m2 s-1 or Pa s]
    tbl_thick, &! The thickness of the top boundary layer [Z ~> m]
    Kv_add, &   ! A viscosity to add [H Z T-1 ~> m2 s-1 or Pa s]
    Kv_tot      ! The total viscosity at an interface [H Z T-1 ~> m2 s-1 or Pa s]
  integer :: &
    nk_in_ml      ! The index of the deepest interface in the mixed layer.
  real :: h_shear ! The distance over which shears occur [Z ~> m].
  real :: dhc     ! The distance between the center of adjacent layers [Z ~> m].
  real :: visc_ml ! The mixed layer viscosity [H Z T-1 ~> m2 s-1 or Pa s].
  real :: I_Hmix  ! The inverse of the mixed layer thickness [Z-1 ~> m-1].
  real :: a_ml    ! The layer coupling coefficient across an interface in
                  ! the mixed layer [H T-1 ~> m s-1 or Pa s m-1].
  real :: a_floor ! A lower bound on the layer coupling coefficient across an interface in
                  ! the mixed layer [H T-1 ~> m s-1 or Pa s m-1].
  real :: I_amax  ! The inverse of the maximum coupling coefficient [T H-1 ~> s m-1 or s m2 kg-1].
  real :: temp1   ! A temporary variable [Z2 ~> m2]
  real :: ustar2_denom ! A temporary variable in the surface boundary layer turbulence
                  ! calculations [H Z-1 T-1 ~> s-1 or kg m-3 s-1]
  real :: h_neglect ! A vertical distance that is so small it is usually lost
                  ! in roundoff and can be neglected [Z ~> m].
  real :: z2      ! A copy of z_i [nondim]
  real :: botfn   ! A function that is 1 at the bottom and small far from it [nondim]
  real :: topfn   ! A function that is 1 at the top and small far from it [nondim]
  real :: kv_top  ! A viscosity associated with the top boundary layer [H Z T-1 ~> m2 s-1 or Pa s]
  logical :: do_shelf, do_OBCs, can_exit
  integer :: k
  integer :: nz, max_nk

  nz = ke

  h_neglect = dZ_subroundoff

  if (CS%answer_date < 20190101) then
    !   The maximum coupling coefficient was originally introduced to avoid
    ! truncation error problems in the tridiagonal solver. Effectively, the 1e-10
    ! sets the maximum coupling coefficient increment to 1e10 m per timestep.
    I_amax = (1.0e-10*H_to_m) * dt
  else
    I_amax = 0.0
  endif

  do_shelf = .false. ; if (present(shelf)) do_shelf = shelf

  do_OBCs = .false.
  if (associated(OBC)) then
    if (work_on_u) then
      do_OBCS = OBC%u_E_OBCs_on_PE .or. OBC%u_W_OBCs_on_PE
    else
      do_OBCS = OBC%v_N_OBCs_on_PE .or. OBC%v_S_OBCs_on_PE
    endif
  endif

  a_cpl(:) = 0.
  h_ml = 0.

  if (CS%Kvml_invZ2 > 0. .and. .not. do_shelf) then
    I_Hmix = 1. / (CS%Hmix + h_neglect)
    z_t = h_neglect * I_Hmix
  endif

  do K=2,nz
    Kv_tot = CS%Kv

    if (CS%Kvml_invZ2 > 0. .and. .not. do_shelf) then
      ! This is an older (vintage ~1997) way to prevent wind stresses from driving very
      ! large flows in nearly massless near-surface layers when there is not a physically-
      ! based surface boundary layer parameterization.  It does not have a plausible
      ! physical basis, and probably should not be used.
      z_t = z_t + h_harm(k-1) * I_Hmix
      Kv_tot = CS%Kv + CS%Kvml_invZ2 / ((z_t * z_t) *  &
               (1. + 0.09 * z_t * z_t * z_t * z_t * z_t * z_t))
    endif

    if (associated(Kv_shear)) then
      ! Add in viscosities that are determined by physical processes that are handled in
      ! other modules, and which do not respond immediately to the changing layer thicknesses.
      ! These processes may include shear-driven mixing or contributions from some boundary
      ! layer turbulence schemes.  Other viscosity contributions that respond to the evolving
      ! layer thicknesses or the surface wind stresses are added later.
      if (work_on_u) then
        Kv_add = 0.5 * (Kv_shear(i,j,k) + Kv_shear(i+1,j,k))

        if (do_OBCs) then
          if (OBC%u_E_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) > 0) then
              Kv_add = Kv_shear(i,j,k)
            endif
          endif

          if (OBC%u_W_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) < 0) then
              Kv_add = Kv_shear(i+1,j,k)
            endif
          endif
        endif

        Kv_tot = Kv_tot + Kv_add
      else
        Kv_add = 0.5 * (Kv_shear(i,j,k) + Kv_shear(i,j+1,k))

        if (do_OBCs) then
          if (OBC%v_N_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) > 0) then
              Kv_add = Kv_shear(i,j,k)
            endif
          endif

          if (OBC%v_S_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) < 0) then
              Kv_add = Kv_shear(i,j+1,k)
            endif
          endif
        endif

        Kv_tot = Kv_tot + Kv_add
      endif
    endif

    if (associated(Kv_shear_Bu)) then
      ! This is similar to what was done above, but for contributions coming from the corner
      ! (vorticity) points.  Because OBCs run through the faces and corners there is no need
      ! to further modify these viscosities here to take OBCs into account.
      if (work_on_u) then
        Kv_tot = Kv_tot + 0.5 * (Kv_shear_Bu(I,J-1,k) + Kv_shear_Bu(I,J,k))
      else
        Kv_tot = Kv_tot + 0.5 * (Kv_shear_Bu(I-1,J,k) + Kv_shear_Bu(I,J,k))
      endif
    endif

    ! Set the viscous coupling coefficients, excluding surface mixed layer contributions
    ! for now, but including viscous bottom drag, working up from the bottom.
    if (CS%bottomdraglaw) then
      !    botfn determines when a point is within the influence of the bottom
      !  boundary layer, going from 1 at the bottom to 0 in the interior.
      z2 = z_i(k)
      botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

      Kv_tot = Kv_tot + (kv_bbl - CS%Kv) * botfn
      dhc = 0.5 * (hvel(k) + hvel(k-1))
      if (dhc > bbl_thick) then
        h_shear = ((1. - botfn) * dhc + botfn * bbl_thick) + h_neglect
      else
        h_shear = dhc + h_neglect
      endif

      ! Calculate the coupling coefficients from the viscosities.
      a_cpl(K) = Kv_tot / (h_shear + (I_amax * Kv_tot))
    elseif (abs(CS%Kv_extra_bbl) > 0.0) then
      ! There is a simple enhancement of the near-bottom viscosities, but no
      ! adjustment of the viscous coupling length scales to give a particular
      ! bottom stress.

      !    botfn determines when a point is within the influence of the bottom
      !  boundary layer, going from 1 at the bottom to 0 in the interior.
      z2 = z_i(k)
      botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

      Kv_tot = Kv_tot + CS%Kv_extra_bbl * botfn
      h_shear = 0.5 * (hvel(k) + hvel(k-1) + h_neglect)

      ! Calculate the coupling coefficients from the viscosities.
      a_cpl(K) = Kv_tot / (h_shear + I_amax * Kv_tot)
    else
      ! Any near-bottom viscous enhancements were already incorporated into
      ! Kv_tot, and there is no adjustment of the viscous coupling length
      ! scales to give a particular bottom stress.

      h_shear = 0.5 * (hvel(k) + hvel(k-1) + h_neglect)
      ! Calculate the coupling coefficients from the viscosities.
      a_cpl(K) = Kv_tot / (h_shear + I_amax * Kv_tot)
    endif
  enddo

  ! Assign the bottom coupling coefficients
  if (CS%bottomdraglaw) then
    dhc = hvel(nz) * 0.5
    a_cpl(nz+1) = kv_bbl / ((min(dhc, bbl_thick) + h_neglect) + I_amax * kv_bbl)
  elseif (abs(CS%Kv_extra_bbl) > 0.0) then
    a_cpl(nz+1) = (CS%Kv + CS%Kv_extra_bbl) &
        / ((0.5 * hvel(nz) + h_neglect) + I_amax * (CS%Kv + CS%Kv_extra_bbl))
  else
    a_cpl(nz+1) = CS%Kv / ((0.5 * hvel(nz) + h_neglect) + I_amax * CS%Kv)
  endif

  ! Add surface intensified viscous coupling, either as a no-slip boundary condition under a
  ! rigid ice-shelf, or due to wind-stress driven surface boundary layer mixing that has not
  ! already been added via visc%Kv_shear.
  if (do_shelf) then
    ! Set the coefficients to include the no-slip surface stress.
    if (work_on_u) then
      kv_TBL = Kv_tbl_shelf_u(I,j)
      tbl_thick = tbl_thick_shelf_u(I,j) + h_neglect
    else
      kv_TBL = Kv_tbl_shelf_v(i,J)
      tbl_thick = tbl_thick_shelf_v(i,J) + h_neglect
    endif

    z_t = 0.0

    ! If a_cpl(1) were not already 0, it would be added here.
    if (0.5 * hvel(1) > tbl_thick) then
      a_cpl(1) = kv_TBL / (tbl_thick + I_amax * kv_TBL)
    else
      a_cpl(1) = kv_TBL / ((0.5 * hvel(1) + h_neglect) + I_amax * kv_TBL)
    endif

    do K=2,nz
      z_t = z_t + hvel(k-1) / tbl_thick
      topfn = 1. / (1. + 0.09 * z_t**6)

      dhc = 0.5 * (hvel(k) + hvel(k-1))
      if (dhc > tbl_thick) then
        h_shear = ((1. - topfn) * dhc + topfn * tbl_thick) + h_neglect
      else
        h_shear = dhc + h_neglect
      endif

      kv_top = topfn * kv_TBL
      a_cpl(K) = a_cpl(K) + kv_top / (h_shear + I_amax * kv_top)
    enddo
  elseif (CS%dynamic_viscous_ML .or. (nkml>0) .or. CS%fixed_LOTW_ML .or. CS%apply_LOTW_floor) then

    ! Find the friction velocity and the absolute value of the Coriolis parameter at this point.
    u_star = 0.   ! Zero out the friction velocity on land points.
    tau_mag = 0.  ! Zero out the friction velocity on land points.

    if (associated(SpV_avg)) then
      rho_av1 = 0.

      if (work_on_u) then
        u_star = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i+1,j))
        rho_av1 = 2. / (SpV_avg(i,j,1) + SpV_avg(i+1,j,1))
        absf = 0.5 * (abs(CoriolisBu(I,J-1)) + abs(CoriolisBu(I,J)))

        if (do_OBCs) then
          if (OBC%u_E_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) > 0) then
              u_star = Ustar_2d(i,j)
              rho_av1 = 1. / SpV_avg(i,j,1)
            endif
          endif

          if (OBC%u_W_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) < 0) then
              u_star = Ustar_2d(i+1,j)
              rho_av1 = 1. / SpV_avg(i+1,j,1)
            endif
          endif
        endif
      else
        u_star = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i,j+1))
        rho_av1 = 2. / (SpV_avg(i,j,1) + SpV_avg(i,j+1,1))
        absf = 0.5 * (abs(CoriolisBu(I-1,J)) + abs(CoriolisBu(I,J)))

        if (do_OBCs) then
          if (OBC%v_N_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) > 0) then
              u_star = Ustar_2d(i,j)
              rho_av1 = 1. / SpV_avg(i,j,1)
            endif
          endif

          if (OBC%v_S_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) < 0) then
              u_star = Ustar_2d(i,j+1)
              rho_av1 = 1. / SpV_avg(i,j+1,1)
            endif
          endif
        endif
      endif

      tau_mag = RZ_to_H * rho_av1 * u_star**2
    else ! (.not.allocated(tv%SpV_avg))
      if (work_on_u) then
        u_star = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i+1,j))
        absf = 0.5 * (abs(CoriolisBu(I,J-1)) + abs(CoriolisBu(I,J)))

        if (do_OBCs) then
          if (OBC%u_E_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) > 0) then
              u_star = Ustar_2d(i,j)
            endif
          endif

          if (OBC%u_W_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) < 0) then
              u_star = Ustar_2d(i+1,j)
            endif
          endif
        endif
      else
        u_star = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i,j+1))
        absf = 0.5 * (abs(CoriolisBu(I-1,J)) + abs(CoriolisBu(I,J)))

        if (do_OBCs) then
          if (OBC%v_N_OBCs_on_PE) then
              if (OBC%segnum_v(i,J) > 0) then
                u_star = Ustar_2d(i,j)
              endif
          endif

          if (OBC%v_S_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) < 0) then
              u_star = Ustar_2d(i,j+1)
            endif
          endif
        endif
      endif

      tau_mag = Z_to_H * u_star**2
    endif

    ! Determine the thickness of the surface ocean boundary layer and its extent in index space.
    nk_in_ml = 0
    if (CS%dynamic_viscous_ML) then
      ! The fractional number of layers that are within the viscous boundary layer were
      ! previously stored in visc%nkml_visc_[uv].
      h_ml = h_neglect
      max_nk = 0

      if (work_on_u) then
        nk_in_ml = ceiling(nkml_visc_u(I,j))
        max_nk = max(max_nk, nk_in_ml)

        do k=1,max_nk
          if (k <= nkml_visc_u(I,j)) then ! This layer is all in the ML.
            h_ml = h_ml + hvel(k)
          elseif (k < nkml_visc_u(I,j) + 1.) then ! Part of this layer is in the ML.
            h_ml = h_ml + ((nkml_visc_u(I,j) + 1.) - k) * hvel(k)
          endif
        enddo
      else
        nk_in_ml = ceiling(nkml_visc_v(i,J))
        max_nk = max(max_nk, nk_in_ml)

        do k=1,max_nk
          if (k <= nkml_visc_v(i,J)) then ! This layer is all in the ML.
            h_ml = h_ml + hvel(k)
          elseif (k < nkml_visc_v(i,J) + 1.) then ! Part of this layer is in the ML.
            h_ml = h_ml + ((nkml_visc_v(i,J) + 1.) - k) * hvel(k)
          endif
        enddo
      endif
    elseif (nkml>0) then
      ! This is a simple application of a refined-bulk mixed layer with GV%nkml sublayers.
      max_nk = nkml
      nk_in_ml = nkml

      h_ml = h_neglect

      do k=1,nkml
        h_ml = h_ml + hvel(k)
      enddo
    elseif (CS%fixed_LOTW_ML .or. CS%apply_LOTW_floor) then
      ! Determine which interfaces are within CS%Hmix of the surface, and set the viscous
      ! boundary layer thickness to the smaller of CS%Hmix and the depth of the ocean.
      h_ml = 0.0
      do k=1,nz
        can_exit = .true.
        if (h_ml < CS%Hmix) then
          nk_in_ml = k

          if (h_ml + hvel(k) < CS%Hmix) then
            h_ml = h_ml + hvel(k)
            can_exit = .false.  ! Part of the next deeper layer is also in the mixed layer.
          else
            h_ml = CS%Hmix
          endif
        endif

        if (can_exit) exit  ! All remaining layers in this row are below the mixed layer depth.
      enddo

      max_nk = max(0, nk_in_ml)
    endif

    ! Avoid working on columns where the viscous coupling could not be increased.
    if (u_star <= 0.) nk_in_ml = 0

    ! Set the viscous coupling at the interfaces as the larger of what was previously
    ! set and the contributions from the surface boundary layer.
    z_t = 0.
    if (CS%apply_LOTW_floor .and. &
        (CS%dynamic_viscous_ML .or. nkml > 0 .or. CS%fixed_LOTW_ML)) then
      do K=2,max_nk
        if (k <= nk_in_ml) then
          z_t = z_t + hvel(k-1)

          !   The viscosity in visc_ml is set to go to 0 at the mixed layer top and bottom
          ! (in a log-layer) and be further limited by rotation to give the natural Ekman length.
          temp1 = (z_t * h_ml - z_t * z_t)
          if (Boussinesq) then
            ustar2_denom = (CS%vonKar * Z_to_H * u_star**2) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          else
            ustar2_denom = (CS%vonKar * tau_mag) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          endif

          visc_ml = temp1 * ustar2_denom
          ! Set the viscous coupling based on the model's vertical resolution.  The omission of
          ! the I_amax factor here is consistent with answer dates above 20190101.
          a_ml = visc_ml / (0.25 * (hvel(k) + hvel(k-1) + h_neglect))

          ! As a floor on the viscous coupling, assume that the length scale in the denominator can
          ! not be larger than the distance from the surface, consistent with a logarithmic velocity
          ! profile.  This is consistent with visc_ml, but cancels out common factors of z_t.
          a_floor = (h_ml - z_t) * ustar2_denom

          ! Choose the largest estimate of a_cpl.
          a_cpl(K) = max(a_cpl(K), a_ml, a_floor)
          ! An option could be added to change this to: a_cpl(i,K) = max(a_cpl(i,K) + a_ml, a_floor)
        endif
      enddo
    elseif (CS%apply_LOTW_floor) then
      do K=2,max_nk
        if (k <= nk_in_ml) then
          z_t = z_t + hvel(k-1)

          temp1 = (z_t * h_ml - z_t * z_t)
          if (Boussinesq) then
            ustar2_denom = (CS%vonKar * Z_to_H * u_star**2) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          else
            ustar2_denom = (CS%vonKar * tau_mag) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          endif

          ! As a floor on the viscous coupling, assume that the length scale in the denominator can not
          ! be larger than the distance from the surface, consistent with a logarithmic velocity profile.
          a_cpl(K) = max(a_cpl(K), (h_ml - z_t) * ustar2_denom)
        endif
      enddo
    else
      do K=2,max_nk
        if (k <= nk_in_ml) then
          z_t = z_t + hvel(k-1)

          temp1 = (z_t * h_ml - z_t * z_t)
          !   This viscosity is set to go to 0 at the mixed layer top and bottom (in a log-layer)
          ! and be further limited by rotation to give the natural Ekman length.
          ! The following expressions are mathematically equivalent.
          if (Boussinesq .or. (CS%answer_date < 20230601)) then
            visc_ml = u_star * CS%vonKar * (Z_to_H * temp1 * u_star) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          else
            visc_ml = CS%vonKar * (temp1 * tau_mag) &
                / (absf * temp1 + (h_ml + h_neglect) * u_star)
          endif
          a_ml = visc_ml / (0.25 * (hvel(k) + hvel(k-1) + h_neglect) + 0.5 * I_amax * visc_ml)

          ! Choose the largest estimate of a_cpl, but these could be changed to be additive.
          a_cpl(K) = max(a_cpl(K), a_ml)
          ! An option could be added to change this to: a_cpl(i,K) = a_cpl(i,K) + a_ml
        endif
      enddo
    endif
  endif
end subroutine find_coupling_coef


!> Velocity components which exceed a threshold for physically reasonable values are truncated,
!! and the running sum of the number of trunctionas within the non-symmetric memory computational
!! domain is incremented.  Optionally, any column with excessive velocities may be sent
!! to a diagnostic reporting subroutine.
subroutine vertvisc_limit_vel(bxU0, bxV0, u_a, v_a, dt, Gcore, vgrid, CS, coupling_state, &
                              report_a, ntrunc_sum)
  type(Box_t),             intent(in)    :: bxU0   !< The u-point box of the whole column,
                                                   !! [IscB:IecB, jsc:jec, 1:nz]
  type(Box_t),             intent(in)    :: bxV0   !< The v-point box of the whole column,
                                                   !! [isc:iec, JscB:JecB, 1:nz]
  type(grid_core_type),   intent(in)    :: Gcore  !< Persistent copies of the core ocean grid fields
  type(vgrid_consts_CS),   intent(in)    :: vgrid  !< The vertical grid and unit-scaling constants
  type(RealArray_t),       intent(inout) :: u_a    !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t),       intent(inout) :: v_a    !< Meridional velocity [L T-1 ~> m s-1]
  real,                    intent(in)    :: dt     !< Time increment [T ~> s]
  type(vertvisc_opts_CS),  intent(in)    :: CS     !< Vertical viscosity control structure
  type(coupling_state_type), &
                           intent(in)    :: coupling_state !< The viscous coupling coefficients and
                                                           !! effective layer thicknesses at
                                                           !! velocity points
  type(limit_vel_report_type), &
                           intent(inout) :: report_a !< The velocities and flags used to report the
                                                     !! columns in which velocity truncations occur
  integer,                 intent(inout) :: ntrunc_sum !< The number of velocity truncations, added
                                                       !! to the value on entry [nondim]

  ! Local variables
  real :: CFL              ! The local CFL number [nondim]
  real :: H_report         ! A thickness below which not to report truncations [H ~> m or kg m-2]
  real, dimension(:,:), contiguous, pointer :: vel_report ! The velocity to report [L T-1 ~> m s-1]
  real, dimension(:,:,:), contiguous, pointer :: u_old ! The previous u-velocity [L T-1 ~> m s-1]
  real, dimension(:,:,:), contiguous, pointer :: v_old ! The previous v-velocity [L T-1 ~> m s-1]
  logical, dimension(:,:), contiguous, pointer :: dowrite ! True for the columns to report
  logical :: trunc_any
  logical :: do_any_write
  integer :: ntrunc         ! Thread-safe truncation counter
  real, dimension(:,:,:), contiguous, pointer :: u, v ! The zonal and meridional velocities
                                                    ! [L T-1 ~> m s-1]
  real, dimension(:,:,:), contiguous, pointer :: h_u, h_v ! The effective layer thicknesses at
                                                    ! u- and v-points [H ~> m or kg m-2]
  real, dimension(:,:), contiguous, pointer :: areaT ! The area of an h-cell [L2 ~> m2]
  real, dimension(:,:), contiguous, pointer :: IareaT ! 1/areaT [L-2 ~> m-2]
  real, dimension(:,:), contiguous, pointer :: dx_Cv ! The unblocked lengths of the v-faces of the
                                                     ! h-cell [L ~> m]
  real, dimension(:,:), contiguous, pointer :: dy_Cu ! The unblocked lengths of the u-faces of the
                                                     ! h-cell [L ~> m]
  integer :: i, j, k, is, ie, js, je, nz
  is = bxV0%idxS(1) ; ie = bxV0%idxE(1) ; js = bxU0%idxS(2) ; je = bxU0%idxE(2) ; nz = vgrid%ke

  call u_a%view(u)
  call v_a%view(v)
  call coupling_state%h_u%view(h_u)
  call coupling_state%h_v%view(h_v)
  call Gcore%areaT%view(areaT)
  call Gcore%IareaT%view(IareaT)
  call Gcore%dx_Cv%view(dx_Cv)
  call Gcore%dy_Cu%view(dy_Cu)

  nullify(u_old, v_old, vel_report, dowrite)
  if (report_a%u_old_a%associated()) call report_a%u_old_a%view(u_old)
  if (report_a%v_old_a%associated()) call report_a%v_old_a%view(v_old)
  if (report_a%vel_report_u_a%associated()) call report_a%vel_report_u_a%view(vel_report)
  if (report_a%dowrite_u_a%associated()) call report_a%dowrite_u_a%view(dowrite)

  H_report = 3.0 * vgrid%Angstrom_H

  !$omp target enter data map(alloc: dowrite, vel_report)
  !$omp target enter data map(alloc: u_old, v_old)

  if (CS%do_report_u) then
    do_any_write = .false.
    trunc_any = .false.

    do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1))
      dowrite(I,j) = .false.
      vel_report(I,j) = 3.0e8 * vgrid%m_s_to_L_T
    enddo

    do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1)) &
        DO_LOCALITY(local(k, CFL) reduce(.or.: trunc_any, do_any_write))
      do k=1,nz
        if (abs(u(I,j,k)) < CS%vel_underflow) u(I,j,k) = 0.0
        if (u(I,j,k) < 0.0) then
          CFL = (-u(I,j,k) * dt) * (dy_Cu(I,j) * IareaT(i+1,j))
        else
          CFL = (u(I,j,k) * dt) * (dy_Cu(I,j) * IareaT(i,j))
        endif
        if (CFL > CS%CFL_trunc) trunc_any = .true.
        if (CFL > CS%CFL_report) then
          dowrite(I,j) = .true.
          do_any_write = .true.
          vel_report(I,j) = min(vel_report(I,j), abs(u(I,j,k)))
        endif
      enddo
    enddo

    do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1), dowrite(I,j))
      u_old(I,j,:) = u(I,j,:)
    enddo

    if (trunc_any) then
      ntrunc = 0
      do concurrent (k=bxU0%idxS(3):bxU0%idxE(3), j=bxU0%idxS(2):bxU0%idxE(2), &
                     I=bxU0%idxS(1):bxU0%idxE(1)) DO_LOCALITY(reduce(+: ntrunc))
        if ((u(I,j,k) * (dt * dy_Cu(I,j))) * IareaT(i+1,j) < -CS%CFL_trunc) then
          u(I,j,k) = (-0.9*CS%CFL_trunc) * (areaT(i+1,j) / (dt * dy_Cu(I,j)))
          if (((I >= is) .and. (I <= ie) .and. (j >= js) .and. (j <= je)) .and. &
              (h_u(I,j,k) > H_report)) ntrunc = ntrunc + 1
        elseif ((u(I,j,k) * (dt * dy_Cu(I,j))) * IareaT(i,j) > CS%CFL_trunc) then
          u(I,j,k) = (0.9*CS%CFL_trunc) * (areaT(i,j) / (dt * dy_Cu(I,j)))
          if (((I >= is) .and. (I <= ie) .and. (j >= js) .and. (j <= je)) .and. &
              (h_u(I,j,k) > H_report)) ntrunc = ntrunc + 1
        endif
      enddo
      ntrunc_sum = ntrunc_sum + ntrunc
    endif

  else  ! Do not report accelerations leading to large velocities.
    ntrunc = 0
    do concurrent (k=bxU0%idxS(3):bxU0%idxE(3), j=bxU0%idxS(2):bxU0%idxE(2), &
                   I=bxU0%idxS(1):bxU0%idxE(1)) DO_LOCALITY(reduce(+: ntrunc))
      if (abs(u(I,j,k)) < CS%vel_underflow) then ; u(I,j,k) = 0.0
      elseif ((u(I,j,k) * (dt * dy_Cu(I,j))) * IareaT(i+1,j) < -CS%CFL_trunc) then
        u(I,j,k) = (-0.9*CS%CFL_trunc) * (areaT(i+1,j) / (dt * dy_Cu(I,j)))
        if (((I >= is) .and. (I <= ie) .and. (j >= js) .and. (j <= je)) .and. &
            (h_u(I,j,k) > H_report)) ntrunc = ntrunc + 1
      elseif ((u(I,j,k) * (dt * dy_Cu(I,j))) * IareaT(i,j) > CS%CFL_trunc) then
        u(I,j,k) = (0.9*CS%CFL_trunc) * (areaT(i,j) / (dt * dy_Cu(I,j)))
        if (((I >= is) .and. (I <= ie) .and. (j >= js) .and. (j <= je)) .and. &
            (h_u(I,j,k) > H_report)) ntrunc = ntrunc + 1
      endif
    enddo
    ntrunc_sum = ntrunc_sum + ntrunc
  endif

  if (CS%do_report_v) then
    call report_a%vel_report_v_a%view(vel_report)
    call report_a%dowrite_v_a%view(dowrite)
    do_any_write =.false.
    trunc_any = .false.

    do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1))
      dowrite(i,J) = .false.
      vel_report(i,J) = 3.0e8 * vgrid%m_s_to_L_T
    enddo

    do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1)) &
        DO_LOCALITY(local(k, CFL) reduce(.or.: trunc_any, do_any_write))
      do k=1,nz
        if (abs(v(i,J,k)) < CS%vel_underflow) v(i,J,k) = 0.0
        if (v(i,J,k) < 0.0) then
          CFL = (-v(i,J,k) * dt) * (dx_Cv(i,J) * IareaT(i,j+1))
        else
          CFL = (v(i,J,k) * dt) * (dx_Cv(i,J) * IareaT(i,j))
        endif
        if (CFL > CS%CFL_trunc) trunc_any = .true.
        if (CFL > CS%CFL_report) then
          dowrite(i,J) = .true.
          do_any_write = .true.
          vel_report(i,J) = min(vel_report(i,J), abs(v(i,J,k)))
        endif
      enddo
    enddo

    do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1), dowrite(i,J))
      v_old(i,J,:) = v(i,J,:)
    enddo

    if (trunc_any) then
      ntrunc = 0
      do concurrent (k=bxV0%idxS(3):bxV0%idxE(3), J=bxV0%idxS(2):bxV0%idxE(2), &
                     i=bxV0%idxS(1):bxV0%idxE(1)) DO_LOCALITY(reduce(+: ntrunc))
        if ((v(i,J,k) * (dt * dx_Cv(i,J))) * IareaT(i,j+1) < -CS%CFL_trunc) then
          v(i,J,k) = (-0.9*CS%CFL_trunc) * (areaT(i,j+1) / (dt * dx_Cv(i,J)))
          if (((i >= is) .and. (i <= ie) .and. (J >= js) .and. (J <= je)) .and. &
              (h_v(i,J,k) > H_report)) ntrunc = ntrunc + 1
        elseif ((v(i,J,k) * (dt * dx_Cv(i,J))) * IareaT(i,j) > CS%CFL_trunc) then
          v(i,J,k) = (0.9*CS%CFL_trunc) * (areaT(i,j) / (dt * dx_Cv(i,J)))
          if (((i >= is) .and. (i <= ie) .and. (J >= js) .and. (J <= je)) .and. &
              (h_v(i,J,k) > H_report)) ntrunc = ntrunc + 1
        endif
      enddo
      ntrunc_sum = ntrunc_sum + ntrunc
    endif

  else  ! Do not report accelerations leading to large velocities.
    ntrunc = 0
    do concurrent (k=bxV0%idxS(3):bxV0%idxE(3), J=bxV0%idxS(2):bxV0%idxE(2), &
                   i=bxV0%idxS(1):bxV0%idxE(1)) DO_LOCALITY(reduce(+: ntrunc))
      if (abs(v(i,J,k)) < CS%vel_underflow) then ; v(i,J,k) = 0.0
      elseif ((v(i,J,k) * (dt * dx_Cv(i,J))) * IareaT(i,j+1) < -CS%CFL_trunc) then
        v(i,J,k) = (-0.9*CS%CFL_trunc) * (areaT(i,j+1) / (dt * dx_Cv(i,J)))
        if (((i >= is) .and. (i <= ie) .and. (J >= js) .and. (J <= je)) .and. &
            (h_v(i,J,k) > H_report)) ntrunc = ntrunc + 1
      elseif ((v(i,J,k) * (dt * dx_Cv(i,J))) * IareaT(i,j) > CS%CFL_trunc) then
        v(i,J,k) = (0.9*CS%CFL_trunc) * (areaT(i,j) / (dt * dx_Cv(i,J)))
        if (((i >= is) .and. (i <= ie) .and. (J >= js) .and. (J <= je)) .and. &
            (h_v(i,J,k) > H_report)) ntrunc = ntrunc + 1
      endif
    enddo
    ntrunc_sum = ntrunc_sum + ntrunc
  endif

  !$omp target exit data map(release: u_old, v_old, dowrite, vel_report)
end subroutine vertvisc_limit_vel

!> Writes the accelerations in the columns with large velocities to the truncation files.
subroutine vertvisc_limit_vel_report(report_a, h, ADp, CDp, dt, taux, tauy, G, GV, US, CS)
  type(ocean_grid_type),   intent(in) :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(limit_vel_report_type), &
                           intent(in) :: report_a !< The velocities and flags used to report the
                                                  !! columns in which velocity truncations occur
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h   !< Layer thickness [H ~> m or kg m-2]
  type(accel_diag_ptrs),   intent(in) :: ADp !< Acceleration diagnostic pointers
  type(cont_diag_ptrs),    intent(in) :: CDp !< Continuity diagnostic pointers
  real,                    intent(in) :: dt  !< Time increment [T ~> s]
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in) :: taux !< zonal wind stress [R L Z T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(in) :: tauy !< meridional wind stress [R L Z T-2 ~> Pa]
  type(vertvisc_CS),       intent(in) :: CS  !< Vertical viscosity control structure

  ! Local variables
  real, dimension(:,:), contiguous, pointer :: vel_report ! The velocity to report [L T-1 ~> m s-1]
  real, dimension(:,:,:), contiguous, pointer :: u_old ! The previous u-velocity [L T-1 ~> m s-1]
  real, dimension(:,:,:), contiguous, pointer :: v_old ! The previous v-velocity [L T-1 ~> m s-1]
  logical, dimension(:,:), contiguous, pointer :: dowrite ! True for the columns to report
  logical :: do_any_write
  real, dimension(:,:,:), contiguous, pointer :: a_u, a_v ! The u- and v-drag coefficients
                                                    ! across interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(:,:,:), contiguous, pointer :: h_u, h_v ! The effective layer thicknesses at
                                                    ! u- and v-points [H ~> m or kg m-2]
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  call CS%coupling_state%a_u%view(a_u)
  call CS%coupling_state%a_v%view(a_v)
  call CS%coupling_state%h_u%view(h_u)
  call CS%coupling_state%h_v%view(h_v)

  if (CS%vertvisc_opts_CS%do_report_u) then
    call report_a%u_old_a%view(u_old)
    call report_a%vel_report_u_a%view(vel_report)
    call report_a%dowrite_u_a%view(dowrite)
    do_any_write = any(dowrite)
    if (do_any_write) then
      !$omp target update from (u_old, vel_report)
      do j=js,je ; do I=Isq,Ieq ; if (dowrite(I,j)) then
        ! Call a diagnostic reporting subroutines are called if unphysically large values are found.
        call write_u_accel(I, j, u_old, h, ADp, CDp, dt, G, GV, US, CS%PointAccel_CSp, &
                           vel_report(I,j), taux(I,j), a=a_u, hv=h_u)
      endif ; enddo ; enddo
    endif
  endif

  if (CS%vertvisc_opts_CS%do_report_v) then
    call report_a%v_old_a%view(v_old)
    call report_a%vel_report_v_a%view(vel_report)
    call report_a%dowrite_v_a%view(dowrite)
    do_any_write = any(dowrite)
    if (do_any_write) then
      !$omp target update from(v_old, vel_report)
      do J=Jsq,Jeq ; do i=is,ie ; if (dowrite(i,J)) then
        ! Call a diagnostic reporting subroutines are called if unphysically large values are found.
        call write_v_accel(i, J, v_old, h, ADp, CDp, dt, G, GV, US, CS%PointAccel_CSp, &
                           vel_report(i,J), tauy(i,J), a=a_v, hv=h_v)
      endif ; enddo ; enddo
    endif
  endif
end subroutine vertvisc_limit_vel_report


!> Initialize the vertical friction module
subroutine vertvisc_init(MIS, Time, G, GV, US, param_file, diag, ADp, dirs, &
                          ntrunc, CS, fpmix)
  type(ocean_internal_state), &
                   target, intent(in)    :: MIS    !< The "MOM Internal State", a set of pointers
                                                   !! to the fields and accelerations that make
                                                   !! up the ocean's physical state
  type(time_type), target, intent(in)    :: Time   !< Current model time
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< File to parse for parameters
  type(diag_ctrl), target, intent(inout) :: diag   !< Diagnostic control structure
  type(accel_diag_ptrs),   intent(inout) :: ADp    !< Acceleration diagnostic pointers
  type(directories),       intent(in)    :: dirs   !< Relevant directory paths
  integer, target,         intent(inout) :: ntrunc !< Number of velocity truncations
  type(vertvisc_CS),       pointer       :: CS     !< Vertical viscosity control structure
  logical, optional,       intent(in)    :: fpmix  !< Nonlocal momentum mixing

  ! Local variables

  real :: Kv_BBL  ! A viscosity in the bottom boundary layer with a simple scheme [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(:,:), contiguous, pointer :: kappa_gl90_2d ! 2D kappa_gl90 at h-points
                                                    ! [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
  real :: Kv_back_z  ! A background kinematic viscosity [Z2 T-1 ~> m2 s-1]
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  logical :: lfpmix
  character(len=200) :: kappa_gl90_file, inputdir, kdgl90_varname
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_vert_friction" ! This module's name.
  character(len=40)  :: thickness_units
  real :: Kv_mks ! KVML in MKS [m2 s-1]

  CS%initialized = .true.

  if (GV%Boussinesq) then ; thickness_units = "m"
  else ; thickness_units = "kg m-2" ; endif

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%diag => diag ; CS%ntrunc => ntrunc ; ntrunc = 0

  lfpmix = .false.
  if (present(fpmix)) lfpmix = fpmix

! Default, read and log parameters
  call log_version(param_file, mdl, version, "", log_to_all=.true., debugging=.true.)
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "VERT_FRICTION_ANSWER_DATE", &
                 CS%coef_opts_CS%coupling_CS%answer_date, &
                 "The vintage of the order of arithmetic and expressions in the viscous "//&
                 "calculations.  Values below 20190101 recover the answers from the end of 2018, "//&
                 "while higher values use expressions that do not use an arbitrary hard-coded "//&
                 "maximum viscous coupling coefficient between layers.  Values below 20230601 "//&
                 "recover a form of the viscosity within the mixed layer that breaks up the "//&
                 "magnitude of the wind stress in some non-Boussinesq cases.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%coef_opts_CS%coupling_CS%answer_date = &
      max(CS%coef_opts_CS%coupling_CS%answer_date, 20230701)

  call get_param(param_file, mdl, "BOTTOMDRAGLAW", CS%coef_opts_CS%coupling_CS%bottomdraglaw, &
                 "If true, the bottom stress is calculated with a drag "//&
                 "law of the form c_drag*|u|*u. The velocity magnitude "//&
                 "may be an assumed value or it may be based on the "//&
                 "actual velocity in the bottommost HBBL, depending on "//&
                 "LINEAR_DRAG.", default=.true.)
  call get_param(param_file, mdl, "DIRECT_STRESS", CS%vertvisc_opts_CS%direct_stress, &
                 "If true, the wind stress is distributed over the topmost HMIX_STRESS of fluid "//&
                 "(like in HYCOM), and an added mixed layer viscosity or a physically based "//&
                 "boundary layer turbulence parameterization is not needed for stability.", &
                 default=.false.)
  call get_param(param_file, mdl, "DYNAMIC_VISCOUS_ML", &
                 CS%coef_opts_CS%coupling_CS%dynamic_viscous_ML, &
                 "If true, use a bulk Richardson number criterion to "//&
                 "determine the mixed layer thickness for viscosity.", &
                 default=.false.)
  call get_param(param_file, mdl, "FIXED_DEPTH_LOTW_ML", &
                 CS%coef_opts_CS%coupling_CS%fixed_LOTW_ML, &
                 "If true, use a Law-of-the-wall prescription for the mixed layer viscosity "//&
                 "within a boundary layer that is the lesser of HMIX_FIXED and the total "//&
                 "depth of the ocean in a column.", default=.false.)
  call get_param(param_file, mdl, "LOTW_VISCOUS_ML_FLOOR", &
                 CS%coef_opts_CS%coupling_CS%apply_LOTW_floor, &
                 "If true, use a Law-of-the-wall prescription to set a lower bound on the "//&
                 "viscous coupling between layers within the surface boundary layer, based "//&
                 "the distance of interfaces from the surface.  This only acts when there "//&
                 "are large changes in the thicknesses of successive layers or when the "//&
                 "viscosity is set externally and the wind stress has subsequently increased.", &
                 default=.false.)
  call get_param(param_file, mdl, 'VON_KARMAN_CONST', CS%coef_opts_CS%coupling_CS%vonKar, &
                 'The value the von Karman constant as used for mixed layer viscosity.', &
                 units='nondim', default=0.41)
  call get_param(param_file, mdl, "U_TRUNC_FILE", CS%u_trunc_file, &
                 "The absolute path to a file into which the accelerations "//&
                 "leading to zonal velocity truncations are written. "//&
                 "Undefine this for efficiency if this diagnostic is not needed.", &
                 default=" ", debuggingParam=.true.)
  call get_param(param_file, mdl, "V_TRUNC_FILE", CS%v_trunc_file, &
                 "The absolute path to a file into which the accelerations "//&
                 "leading to meridional velocity truncations are written. "//&
                 "Undefine this for efficiency if this diagnostic is not needed.", &
                 default=" ", debuggingParam=.true.)
  call get_param(param_file, mdl, "HARMONIC_VISC", CS%coef_opts_CS%harmonic_visc, &
                 "If true, use the harmonic mean thicknesses for "//&
                 "calculating the vertical viscosity.", default=.false.)
  call get_param(param_file, mdl, "HARMONIC_BL_SCALE", CS%coef_opts_CS%harm_BL_val, &
                 "A scale to determine when water is in the boundary "//&
                 "layers based solely on harmonic mean thicknesses for "//&
                 "the purpose of determining the extent to which the "//&
                 "thicknesses used in the viscosities are upwinded.", &
                 default=0.0, units="nondim")
  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)

  if (GV%nkml < 1) then
    call get_param(param_file, mdl, "HMIX_FIXED", CS%coef_opts_CS%coupling_CS%Hmix, &
                 "The prescribed depth over which the near-surface viscosity and "//&
                 "diffusivity are elevated when the bulk mixed layer is not used.", &
                 units="m", scale=US%m_to_Z, fail_if_missing=.true.)
  endif
  if (CS%vertvisc_opts_CS%direct_stress) then
    if (GV%nkml < 1) then
      call get_param(param_file, mdl, "HMIX_STRESS", CS%vertvisc_opts_CS%Hmix_stress, &
                 "The depth over which the wind stress is applied if DIRECT_STRESS is true.", &
                 units="m", default=US%Z_to_m*CS%coef_opts_CS%coupling_CS%Hmix, scale=GV%m_to_H)
    else
      call get_param(param_file, mdl, "HMIX_STRESS", CS%vertvisc_opts_CS%Hmix_stress, &
                 "The depth over which the wind stress is applied if DIRECT_STRESS is true.", &
                 units="m", fail_if_missing=.true., scale=GV%m_to_H)
    endif
    if (CS%vertvisc_opts_CS%Hmix_stress <= 0.0) call MOM_error(FATAL, "vertvisc_init: " // &
       "HMIX_STRESS must be set to a positive value if DIRECT_STRESS is true.")
  endif
  call get_param(param_file, mdl, "KV", Kv_back_z, &
                 "The background kinematic viscosity in the interior. "//&
                 "The molecular value, ~1e-6 m2 s-1, may be used.", &
                 units="m2 s-1", fail_if_missing=.true., scale=US%m2_s_to_Z2_T)
  ! Convert input kinematic viscosity to dynamic viscosity when non-Boussinesq.
  CS%coef_opts_CS%coupling_CS%Kv = (US%Z2_T_to_m2_s*GV%m2_s_to_HZ_T) * Kv_back_z

  call get_param(param_file, mdl, "USE_GL90_IN_SSW", CS%coef_opts_CS%use_GL90_in_SSW, &
                 "If true, use simpler method to calculate 1/N^2 in GL90 vertical "// &
                 "viscosity coefficient. This method is valid in stacked shallow water mode.", &
                 default=.false.)
  call get_param(param_file, mdl, "KD_GL90", CS%coef_opts_CS%GL90_CS%kappa_gl90, &
                 "The scalar diffusivity used in GL90 vertical viscosity scheme.", &
                 units="m2 s-1", default=0.0, scale=US%m_to_L*US%Z_to_L*GV%m_to_H*US%T_to_s, &
                 do_not_log=.not.CS%coef_opts_CS%use_GL90_in_SSW)
  call get_param(param_file, mdl, "READ_KD_GL90", CS%coef_opts_CS%GL90_CS%read_kappa_gl90, &
                 "If true, read a file (given by KD_GL90_FILE) containing the "//&
                 "spatially varying diffusivity KD_GL90 used in the GL90 scheme.", default=.false., &
                 do_not_log=.not.CS%coef_opts_CS%use_GL90_in_SSW)
  if (CS%coef_opts_CS%GL90_CS%read_kappa_gl90) then
    if (CS%coef_opts_CS%GL90_CS%kappa_gl90 > 0) then
        call MOM_error(FATAL, "MOM_vert_friction.F90, vertvisc_init: KD_GL90 > 0 "// &
              "is not compatible with READ_KD_GL90 = .TRUE. ")
    endif
    call get_param(param_file, mdl, "INPUTDIR", inputdir, &
                 "The directory in which all input files are found.", &
                 default=".", do_not_log=.true.)
    inputdir = slasher(inputdir)
    call get_param(param_file, mdl, "KD_GL90_FILE", kappa_gl90_file, &
                 "The file containing the spatially varying diffusivity used in the "// &
                 "GL90 scheme.", default="kd_gl90.nc", &
                 do_not_log=.not.CS%coef_opts_CS%use_GL90_in_SSW)
    call get_param(param_file, mdl, "KD_GL90_VARIABLE", kdgl90_varname, &
                 "The name of the GL90 diffusivity variable to read "//&
                 "from KD_GL90_FILE.", default="kd_gl90", &
                 do_not_log=.not.CS%coef_opts_CS%use_GL90_in_SSW)
    kappa_gl90_file = trim(inputdir) // trim(kappa_gl90_file)

    call CS%coef_opts_CS%GL90_CS%kappa_gl90_2d%alloc(lb=[G%isd,G%jsd], ub=[G%ied,G%jed], source=0.0)
    call CS%coef_opts_CS%GL90_CS%kappa_gl90_2d%view(kappa_gl90_2d)
    call MOM_read_data(kappa_gl90_file, kdgl90_varname, kappa_gl90_2d(:,:), G%domain, &
                       scale=US%m_to_L*US%Z_to_L*GV%m_to_H*US%T_to_s)
    call pass_var(kappa_gl90_2d, G%domain)
  endif
  call get_param(param_file, mdl, "USE_GL90_N2", CS%coef_opts_CS%GL90_CS%use_GL90_N2, &
                 "If true, use GL90 vertical viscosity coefficient that is depth-independent; "// &
                 "this corresponds to a kappa_GM that scales as N^2 with depth.", &
                 default=.false., do_not_log=.not.CS%coef_opts_CS%use_GL90_in_SSW)
  if (CS%coef_opts_CS%GL90_CS%use_GL90_N2) then
    if (.not. CS%coef_opts_CS%use_GL90_in_SSW) call MOM_error(FATAL, &
           "MOM_vert_friction.F90, vertvisc_init: "//&
           "When USE_GL90_N2=True, USE_GL90_in_SSW must also be True.")
    if (CS%coef_opts_CS%GL90_CS%kappa_gl90 > 0) then
        call MOM_error(FATAL, "MOM_vert_friction.F90, vertvisc_init: KD_GL90 > 0 "// &
              "is not compatible with USE_GL90_N2 = .TRUE. ")
    endif
    if (CS%coef_opts_CS%GL90_CS%read_kappa_gl90) call MOM_error(FATAL, &
           "MOM_vert_friction.F90, vertvisc_init: "//&
           "READ_KD_GL90 = .TRUE. is not compatible with USE_GL90_N2 = .TRUE.")
    call get_param(param_file, mdl, "alpha_GL90", CS%coef_opts_CS%GL90_CS%alpha_gl90, &
                   "Coefficient used to compute a depth-independent GL90 vertical "//&
                   "viscosity via Kv_GL90 = alpha_GL90 * f2. Is only used "// &
                   "if USE_GL90_N2 is true. Note that the implied Kv_GL90 "// &
                   "corresponds to a KD_GL90 that scales as N^2 with depth.", &
                   units="m2 s", default=0.0, scale=GV%m_to_H*US%m_to_Z*US%s_to_T, &
                   do_not_log=.not.CS%coef_opts_CS%use_GL90_in_SSW)
  endif
  call get_param(param_file, mdl, "HBBL_GL90", CS%coef_opts_CS%Hbbl_gl90, &
                 "The thickness of the GL90 bottom boundary layer, "//&
                 "which defines the range over which the GL90 coupling "//&
                 "coefficient is zeroed out, in order to avoid fluxing "//&
                 "momentum into vanished layers over steep topography.", &
                 units="m", default=5.0, scale=US%m_to_Z, &
                 do_not_log=.not.CS%coef_opts_CS%use_GL90_in_SSW)

  CS%coef_opts_CS%coupling_CS%Kvml_invZ2 = 0.0
  if (GV%nkml < 1) then
    call get_param(param_file, mdl, "KVML", Kv_mks, &
                 "The scale for an extra kinematic viscosity in the mixed layer", &
                 units="m2 s-1", default=-1.0, do_not_log=.true.)
    if (Kv_mks >= 0.0) then
      call MOM_error(WARNING, "KVML is a deprecated parameter. Use KV_ML_INVZ2 instead.")
    else
      Kv_mks = 0.0
    endif
    call get_param(param_file, mdl, "KV_ML_INVZ2", CS%coef_opts_CS%coupling_CS%Kvml_invZ2, &
                 "An extra kinematic viscosity in a mixed layer of thickness HMIX_FIXED, "//&
                 "with the actual viscosity scaling as 1/(z*HMIX_FIXED)^2, where z is the "//&
                 "distance from the surface, to allow for finite wind stresses to be "//&
                 "transmitted through infinitesimally thin surface layers.  This is an "//&
                 "older option for numerical convenience without a strong physical basis, "//&
                 "and its use is now discouraged.", &
                 units="m2 s-1", default=Kv_mks, scale=GV%m2_s_to_HZ_T)
  endif

  if (.not.CS%coef_opts_CS%coupling_CS%bottomdraglaw) then
    call get_param(param_file, mdl, "KV_EXTRA_BBL", CS%coef_opts_CS%coupling_CS%Kv_extra_bbl, &
                 "An extra kinematic viscosity in the benthic boundary layer. "//&
                 "KV_EXTRA_BBL is not used if BOTTOMDRAGLAW is true.", &
                 units="m2 s-1", default=0.0, scale=GV%m2_s_to_HZ_T, do_not_log=.true.)
    if (CS%coef_opts_CS%coupling_CS%Kv_extra_bbl == 0.0) then
      call get_param(param_file, mdl, "KVBBL", Kv_BBL, &
                 "An extra kinematic viscosity in the benthic boundary layer. "//&
                 "KV_EXTRA_BBL is not used if BOTTOMDRAGLAW is true.", &
                 units="m2 s-1", default=US%Z2_T_to_m2_s*Kv_back_z, scale=GV%m2_s_to_HZ_T, &
                 do_not_log=.true.)
      if (abs(Kv_BBL - CS%coef_opts_CS%coupling_CS%Kv) > &
          1.0e-15*abs(CS%coef_opts_CS%coupling_CS%Kv)) then
        call MOM_error(WARNING, "KVBBL is a deprecated parameter. Use KV_EXTRA_BBL instead.")
        CS%coef_opts_CS%coupling_CS%Kv_extra_bbl = Kv_BBL - CS%coef_opts_CS%coupling_CS%Kv
      endif
    endif
    call log_param(param_file, mdl, "KV_EXTRA_BBL", CS%coef_opts_CS%coupling_CS%Kv_extra_bbl, &
                 "An extra kinematic viscosity in the benthic boundary layer. "//&
                 "KV_EXTRA_BBL is not used if BOTTOMDRAGLAW is true.", &
                 units="m2 s-1", default=0.0, unscale=GV%HZ_T_to_m2_s)
  endif
  call get_param(param_file, mdl, "HBBL", CS%coef_opts_CS%Hbbl, &
                 "The thickness of a bottom boundary layer with a viscosity increased by "//&
                 "KV_EXTRA_BBL if BOTTOMDRAGLAW is not defined, or the thickness over which "//&
                 "near-bottom velocities are averaged for the drag law if BOTTOMDRAGLAW is "//&
                 "defined but LINEAR_DRAG is not.", &
                 units="m", fail_if_missing=.true., scale=US%m_to_Z)
  call get_param(param_file, mdl, "CFL_TRUNCATE", CS%vertvisc_opts_CS%CFL_trunc, &
                 "The value of the CFL number that will cause velocity "//&
                 "components to be truncated; instability can occur past 0.5.", &
                 units="nondim", default=0.5)
  call get_param(param_file, mdl, "CFL_REPORT", CS%vertvisc_opts_CS%CFL_report, &
                 "The value of the CFL number that causes accelerations "//&
                 "to be reported; the default is CFL_TRUNCATE.", &
                 units="nondim", default=CS%vertvisc_opts_CS%CFL_trunc)
  call get_param(param_file, mdl, "CFL_TRUNCATE_RAMP_TIME", CS%truncRampTime, &
                 "The time over which the CFL truncation value is ramped "//&
                 "up at the beginning of the run.", &
                 units="s", default=0., scale=US%s_to_T)
  CS%CFL_truncE = CS%vertvisc_opts_CS%CFL_trunc
  call get_param(param_file, mdl, "CFL_TRUNCATE_START", CS%CFL_truncS, &
                 "The start value of the truncation CFL number used when "//&
                 "ramping up CFL_TRUNC.", &
                 units="nondim", default=0.)
  call get_param(param_file, mdl, "STOKES_MIXING_COMBINED", CS%vertvisc_opts_CS%StokesMixing, &
                 "Flag to use Stokes drift Mixing via the Lagrangian "//&
                 " current (Eulerian plus Stokes drift). "//&
                 " Still needs work and testing, so not recommended for use.",&
                 default=.false.)
  !BGR 04/04/2018{
  ! StokesMixing is required for MOM6 for some Langmuir mixing parameterization.
  !   The code used here has not been developed for vanishing layers or in
  !   conjunction with any bottom friction.  Therefore, the following line is
  !   added so this functionality cannot be used without user intervention in
  !   the code.  This will prevent general use of this functionality until proper
  !   care is given to the previously mentioned issues.  Comment out the following
  !   MOM_error to use, but do so at your own risk and with these points in mind.
  !}
  if (CS%vertvisc_opts_CS%StokesMixing) then
    call MOM_error(FATAL, "Stokes mixing requires user intervention in the code.\n"//&
                          "  Model now exiting.  See MOM_vert_friction.F90 for \n"//&
                          "  details (search 'BGR 04/04/2018' to locate comment).")
  endif
  call get_param(param_file, mdl, "VEL_UNDERFLOW", CS%vertvisc_opts_CS%vel_underflow, &
                 "A negligibly small velocity magnitude below which velocity "//&
                 "components are set to 0.  A reasonable value might be "//&
                 "1e-30 m/s, which is less than an Angstrom divided by "//&
                 "the age of the universe.", units="m s-1", default=0.0, scale=US%m_s_to_L_T)

  call CS%coupling_state%a_u%alloc(lb=[IsdB,jsd,1], ub=[IedB,jed,nz+1], source=0.0)
  call CS%coupling_state%a_u_gl90%alloc(lb=[IsdB,jsd,1], ub=[IedB,jed,nz+1], source=0.0)
  call CS%coupling_state%h_u%alloc(lb=[IsdB,jsd,1], ub=[IedB,jed,nz], source=0.0)
  call CS%coupling_state%a_v%alloc(lb=[isd,JsdB,1], ub=[ied,JedB,nz+1], source=0.0)
  call CS%coupling_state%a_v_gl90%alloc(lb=[isd,JsdB,1], ub=[ied,JedB,nz+1], source=0.0)
  call CS%coupling_state%h_v%alloc(lb=[isd,JsdB,1], ub=[ied,JedB,nz], source=0.0)

  ! TODO: Conditional?

  CS%id_Kv_slow = register_diag_field('ocean_model', 'Kv_slow', diag%axesTi, Time, &
      'Slow varying vertical viscosity', 'm2 s-1', conversion=GV%HZ_T_to_m2_s)

  CS%id_Kv_u = register_diag_field('ocean_model', 'Kv_u', diag%axesCuL, Time, &
      'Total vertical viscosity at u-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_Kv_v = register_diag_field('ocean_model', 'Kv_v', diag%axesCvL, Time, &
      'Total vertical viscosity at v-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_Kv_gl90_u = register_diag_field('ocean_model', 'Kv_gl90_u', diag%axesCuL, Time, &
      'GL90 vertical viscosity at u-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_Kv_gl90_v = register_diag_field('ocean_model', 'Kv_gl90_v', diag%axesCvL, Time, &
      'GL90 vertical viscosity at v-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_au_vv = register_diag_field('ocean_model', 'au_visc', diag%axesCui, Time, &
      'Zonal Viscous Vertical Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_av_vv = register_diag_field('ocean_model', 'av_visc', diag%axesCvi, Time, &
      'Meridional Viscous Vertical Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_au_gl90_vv = register_diag_field('ocean_model', 'au_gl90_visc', diag%axesCui, Time, &
      'Zonal Viscous Vertical GL90 Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_av_gl90_vv = register_diag_field('ocean_model', 'av_gl90_visc', diag%axesCvi, Time, &
      'Meridional Viscous Vertical GL90 Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_h_u = register_diag_field('ocean_model', 'Hu_visc', diag%axesCuL, Time, &
      'Thickness at Zonal Velocity Points for Viscosity', &
      thickness_units, conversion=GV%H_to_MKS)
      ! Alternately, to always give this variable in 'm' use the following line instead:
      ! 'm', conversion=GV%H_to_m)

  CS%id_h_v = register_diag_field('ocean_model', 'Hv_visc', diag%axesCvL, Time, &
      'Thickness at Meridional Velocity Points for Viscosity', &
      thickness_units, conversion=GV%H_to_MKS)

  CS%id_hML_u = register_diag_field('ocean_model', 'HMLu_visc', diag%axesCu1, Time, &
      'Mixed Layer Thickness at Zonal Velocity Points for Viscosity', &
      thickness_units, conversion=US%Z_to_m)

  CS%id_hML_v = register_diag_field('ocean_model', 'HMLv_visc', diag%axesCv1, Time, &
      'Mixed Layer Thickness at Meridional Velocity Points for Viscosity', &
      thickness_units, conversion=US%Z_to_m)

 if (lfpmix) then
  CS%id_uE_h = register_diag_field('ocean_model', 'uE_h' , CS%diag%axesTL, &
      Time, 'x-zonal Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vE_h = register_diag_field('ocean_model', 'vE_h' , CS%diag%axesTL, &
      Time, 'y-merid Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_uInc_h = register_diag_field('ocean_model','uInc_h',CS%diag%axesTL, &
      Time, 'x-zonal Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vInc_h = register_diag_field('ocean_model','vInc_h',CS%diag%axesTL, &
      Time, 'x-zonal Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_uStk = register_diag_field('ocean_model', 'uStk' , CS%diag%axesTL, &
      Time, 'x-FP du increment' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vStk = register_diag_field('ocean_model', 'vStk' , CS%diag%axesTL, &
      Time, 'y-FP dv increment' , 'm s-1', conversion=US%L_T_to_m_s)

  CS%id_FPtau2s = register_diag_field('ocean_model','Omega_tau2s',CS%diag%axesTi, &
      Time, 'Stress direction from shear','radians')
  CS%id_FPtau2w = register_diag_field('ocean_model','Omega_tau2w',CS%diag%axesTi, &
      Time, 'Stress direction from wind','radians')
  CS%id_uStk0 = register_diag_field('ocean_model', 'uStk0' , diag%axesT1, &
      Time, 'Zonal Surface Stokes', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vStk0 = register_diag_field('ocean_model', 'vStk0' , diag%axesT1, &
      Time, 'Merid Surface Stokes', 'm s-1', conversion=US%L_T_to_m_s)
  endif

  CS%id_du_dt_visc = register_diag_field('ocean_model', 'du_dt_visc', diag%axesCuL, Time, &
      'Zonal Acceleration from Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_du_dt_visc > 0) call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  CS%id_dv_dt_visc = register_diag_field('ocean_model', 'dv_dt_visc', diag%axesCvL, Time, &
      'Meridional Acceleration from Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_dv_dt_visc > 0) call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  CS%id_GLwork = register_diag_field('ocean_model', 'GLwork', diag%axesTL, Time, &
      'Sign-definite Kinetic Energy Source from GL90 Vertical Viscosity', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_du_dt_visc_gl90 = register_diag_field('ocean_model', 'du_dt_visc_gl90', diag%axesCuL, Time, &
      'Zonal Acceleration from GL90 Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_du_dt_visc_gl90 > 0) .or. (CS%id_GLwork > 0)) then
    call safe_alloc_ptr(ADp%du_dt_visc_gl90,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  endif
  CS%id_dv_dt_visc_gl90 = register_diag_field('ocean_model', 'dv_dt_visc_gl90', diag%axesCvL, Time, &
      'Meridional Acceleration from GL90 Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_dv_dt_visc_gl90 > 0) .or. (CS%id_GLwork > 0)) then
    call safe_alloc_ptr(ADp%dv_dt_visc_gl90,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  endif
  CS%id_du_dt_str = register_diag_field('ocean_model', 'du_dt_str', diag%axesCuL, Time, &
      'Zonal Acceleration from Surface Wind Stresses', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_du_dt_str > 0) call safe_alloc_ptr(ADp%du_dt_str,IsdB,IedB,jsd,jed,nz)
  CS%id_dv_dt_str = register_diag_field('ocean_model', 'dv_dt_str', diag%axesCvL, Time, &
      'Meridional Acceleration from Surface Wind Stresses', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_dv_dt_str > 0) call safe_alloc_ptr(ADp%dv_dt_str,isd,ied,JsdB,JedB,nz)

  CS%id_taux_bot = register_diag_field('ocean_model', 'taux_bot', diag%axesCu1, &
      Time, 'Zonal Bottom Stress from Ocean to Earth', &
      'Pa', conversion=US%RZ_to_kg_m2*US%L_T2_to_m_s2)
  CS%id_tauy_bot = register_diag_field('ocean_model', 'tauy_bot', diag%axesCv1, &
      Time, 'Meridional Bottom Stress from Ocean to Earth', &
      'Pa', conversion=US%RZ_to_kg_m2*US%L_T2_to_m_s2)

  !CS%id_hf_du_dt_visc = register_diag_field('ocean_model', 'hf_du_dt_visc', diag%axesCuL, Time, &
  !    'Fractional Thickness-weighted Zonal Acceleration from Vertical Viscosity', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_du_dt_visc > 0) then
  !  call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  !  call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  !endif

  !CS%id_hf_dv_dt_visc = register_diag_field('ocean_model', 'hf_dv_dt_visc', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Meridional Acceleration from Vertical Viscosity', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_dv_dt_visc > 0) then
  !  call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  !  call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)
  !endif

  CS%id_hf_du_dt_visc_2d = register_diag_field('ocean_model', 'hf_du_dt_visc_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Vertical Viscosity', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_du_dt_visc_2d > 0) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_hf_dv_dt_visc_2d = register_diag_field('ocean_model', 'hf_dv_dt_visc_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Vertical Viscosity', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_dv_dt_visc_2d > 0) then
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)
  endif

  CS%id_h_du_dt_visc = register_diag_field('ocean_model', 'h_du_dt_visc', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_du_dt_visc > 0) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%diag_hu,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_h_dv_dt_visc = register_diag_field('ocean_model', 'h_dv_dt_visc', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_dv_dt_visc > 0) then
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%diag_hv,isd,ied,JsdB,JedB,nz)
  endif

  CS%id_h_du_dt_str = register_diag_field('ocean_model', 'h_du_dt_str', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Acceleration from Surface Wind Stresses', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_du_dt_str > 0) then
    call safe_alloc_ptr(ADp%du_dt_str,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%diag_hu,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_h_dv_dt_str = register_diag_field('ocean_model', 'h_dv_dt_str', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Acceleration from Surface Wind Stresses', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_dv_dt_str > 0) then
    call safe_alloc_ptr(ADp%dv_dt_str,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%diag_hv,isd,ied,JsdB,JedB,nz)
  endif

  CS%id_du_dt_str_visc_rem = register_diag_field('ocean_model', 'du_dt_str_visc_rem', diag%axesCuL, Time, &
      'Zonal Acceleration from Surface Wind Stresses multiplied by viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_du_dt_str_visc_rem > 0) then
    call safe_alloc_ptr(ADp%du_dt_str,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%visc_rem_u,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_dv_dt_str_visc_rem = register_diag_field('ocean_model', 'dv_dt_str_visc_rem', diag%axesCvL, Time, &
      'Meridional Acceleration from Surface Wind Stresses multiplied by viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_dv_dt_str_visc_rem > 0) then
    call safe_alloc_ptr(ADp%dv_dt_str,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%visc_rem_v,isd,ied,JsdB,JedB,nz)
  endif

  if ((len_trim(CS%u_trunc_file) > 0) .or. (len_trim(CS%v_trunc_file) > 0)) &
    call PointAccel_init(MIS, Time, G, param_file, diag, dirs, CS%PointAccel_CSp)

  CS%vertvisc_opts_CS%do_report_u = (len_trim(CS%u_trunc_file) > 0)
  CS%vertvisc_opts_CS%do_report_v = (len_trim(CS%v_trunc_file) > 0)
  CS%vertvisc_opts_CS%do_GLwork = (CS%id_GLwork > 0)
  CS%vertvisc_opts_CS%do_du_dt_visc_gl90 = (CS%id_du_dt_visc_gl90 > 0)
  CS%vertvisc_opts_CS%do_dv_dt_visc_gl90 = (CS%id_dv_dt_visc_gl90 > 0)

end subroutine vertvisc_init

!> Update the CFL truncation value as a function of time.
!! If called with the optional argument activate=.true., record the
!! value of Time as the beginning of the ramp period.
subroutine updateCFLtruncationValue(Time, CS, US, activate)
  type(time_type), target, intent(in)    :: Time     !< Current model time
  type(vertvisc_CS),       pointer       :: CS       !< Vertical viscosity control structure
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  logical, optional,       intent(in)    :: activate !< Specify whether to record the value of
                                                     !! Time as the beginning of the ramp period

  ! Local variables
  real :: deltaTime ! The time since CS%rampStartTime [T ~> s], which may be negative.
  real :: wghtA     ! The relative weight of the final value [nondim]
  character(len=12) :: msg

  if (CS%truncRampTime==0.) return ! This indicates to ramping is turned off

  ! We use the optional argument to indicate this Time should be recorded as the
  ! beginning of the ramp-up period.
  if (present(activate)) then
    if (activate) then
      CS%rampStartTime = Time ! Record the current time
      CS%CFLrampingIsActivated = .true.
    endif
  endif
  if (.not.CS%CFLrampingIsActivated) return
  deltaTime = max(0., time_minus_signed(Time, CS%rampStartTime, scale=US%s_to_T))
  if (deltaTime >= CS%truncRampTime) then
    CS%vertvisc_opts_CS%CFL_trunc = CS%CFL_truncE
    CS%truncRampTime = 0. ! This turns off ramping after this call
  else
    wghtA = min( 1., deltaTime / CS%truncRampTime ) ! Linear profile in time
    !wghtA = wghtA*wghtA ! Convert linear profile to parabolic profile in time
    !wghtA = wghtA*wghtA*(3. - 2.*wghtA) ! Convert linear profile to cosine profile
    wghtA = 1. - ( (1. - wghtA)**2 ) ! Convert linear profile to inverted parabolic profile
    CS%vertvisc_opts_CS%CFL_trunc = CS%CFL_truncS + wghtA * ( CS%CFL_truncE - CS%CFL_truncS )
  endif
  write(msg(1:12),'(es12.3)') CS%vertvisc_opts_CS%CFL_trunc
  call MOM_error(NOTE, "MOM_vert_friction: updateCFLtruncationValue set CFL limit to "//trim(msg))
end subroutine updateCFLtruncationValue

!> Clean up and deallocate the vertical friction module
subroutine vertvisc_end(CS)
  type(vertvisc_CS), intent(inout) :: CS  !< Vertical viscosity control structure that
                                          !! will be deallocated in this subroutine.

  if ((len_trim(CS%u_trunc_file) > 0) .or. (len_trim(CS%v_trunc_file) > 0)) &
    deallocate(CS%PointAccel_CSp)

  call CS%coupling_state%a_u%free() ; call CS%coupling_state%h_u%free()
  call CS%coupling_state%a_v%free() ; call CS%coupling_state%h_v%free()
  call CS%coupling_state%a1_shelf_u%free()
  call CS%coupling_state%a1_shelf_v%free()
  call CS%coef_opts_CS%GL90_CS%kappa_gl90_2d%free()
end subroutine vertvisc_end

!> \namespace mom_vert_friction
!! \author Robert Hallberg
!! \date April 1994 - October 2006
!!
!!  The vertical diffusion of momentum is fully implicit.  This is
!!  necessary to allow for vanishingly small layers.  The coupling
!!  is based on the distance between the centers of adjacent layers,
!!  except where a layer is close to the bottom compared with a
!!  bottom boundary layer thickness when a bottom drag law is used.
!!  A stress top b.c. and a no slip bottom  b.c. are used.  There
!!  is no limit on the time step for vertvisc.
!!
!!  Near the bottom, the horizontal thickness interpolation scheme
!!  changes to an upwind biased estimate to control the effect of
!!  spurious Montgomery potential gradients at the bottom where
!!  nearly massless layers layers ride over the topography.  Within a
!!  few boundary layer depths of the bottom, the harmonic mean
!!  thickness (i.e. (2 h+ h-) / (h+ + h-) ) is used if the velocity
!!  is from the thinner side and the arithmetic mean thickness
!!  (i.e. (h+ + h-)/2) is used if the velocity is from the thicker
!!  side.  Both of these thickness estimates are second order
!!  accurate.  Above this the arithmetic mean thickness is used.
!!
!!  In addition, vertvisc truncates any velocity component that exceeds a
!!  maximum CFL number to a fraction of this value.  This basically keeps
!!  instabilities spatially localized.  The number of times the velocity is
!!  truncated is reported each time the energies are saved, and if
!!  exceeds CS%Maxtrunc the model will stop itself and change the time
!!  to a large value.  This has proven very useful in (1) diagnosing
!!  model failures and (2) letting the model settle down to a
!!  meaningful integration from a poorly specified initial condition.
!!
!!  The same code is used for the two velocity components, by
!!  indirectly referencing the velocities and defining a handful of
!!  direction-specific defined variables.
!!
!!  Macros written all in capital letters are defined in MOM_memory.h.
!!
!!     A small fragment of the grid is shown below:
!! \verbatim
!!    j+1  x ^ x ^ x   At x:  q
!!    j+1  > o > o >   At ^:  v, frhatv, tauy
!!    j    x ^ x ^ x   At >:  u, frhatu, taux
!!    j    > o > o >   At o:  h
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!! \endverbatim
!!
!!  The boundaries always run through q grid points (x).
end module MOM_vert_friction
