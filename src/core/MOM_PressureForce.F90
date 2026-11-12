! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!!SKILLS: 0.4.0

!> A thin wrapper for Boussinesq/non-Boussinesq forms of the pressure force calculation.
module MOM_PressureForce

use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_PressureForce_FV, only : PressureForce_FV_Bouss, PressureForce_FV_nonBouss
use MOM_PressureForce_FV, only : PressureForce_FV_init
use MOM_PressureForce_FV, only : PressureForce_FV_CS
use MOM_PressureForce_FV, only : PF_FV_host_ctx
use MOM_PressureForce_FV, only : PressureForce_FV_Bouss_pre, PressureForce_FV_Bouss_post
use MOM_PressureForce_FV, only : PressureForce_FV_nonBouss_pre, PressureForce_FV_nonBouss_post
use MOM_PressureForce_Mont, only : PressureForce_Mont_Bouss, PressureForce_Mont_nonBouss
use MOM_PressureForce_Mont, only : PressureForce_Mont_init
use MOM_PressureForce_Mont, only : PressureForce_Mont_CS
use MOM_PressureForce_Mont, only : PF_Mont_host_ctx
use MOM_PressureForce_Mont, only : PressureForce_Mont_Bouss_pre, PressureForce_Mont_Bouss_post
use MOM_PressureForce_Mont, only : PressureForce_Mont_nonBouss_pre, PressureForce_Mont_nonBouss_post
use MOM_self_attr_load, only : SAL_CS
use MOM_tidal_forcing, only : tidal_forcing_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs, accel_diag_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_ALE, only: ALE_CS
use array_mod, only : RealArray_t
use iso_c_binding, only : c_ptr, c_loc
implicit none ; private

#include <MOM_memory.h>

public PressureForce, PressureForce_init

!> Pressure force control structure
type, public :: PressureForce_CS ; private
  logical :: Analytic_FV_PGF !< If true, use the analytic finite volume form
                             !! (Adcroft et al., Ocean Mod. 2008) of the PGF.
  !> Control structure for the analytically integrated finite volume pressure force
  type(PressureForce_FV_CS) :: PressureForce_FV
  !> Control structure for the Montgomery potential form of pressure force
  type(PressureForce_Mont_CS) :: PressureForce_Mont
end type PressureForce_CS

contains

!> A thin layer between the model and the Boussinesq and non-Boussinesq pressure force routines.
subroutine PressureForce(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, ADp, p_atm, pbce, eta)
  type(ocean_grid_type),   target, intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), target, intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   target, intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   target, intent(in)  :: tv   !< A structure pointing to various thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: PFu  !< Zonal pressure force acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out) :: PFv  !< Meridional pressure force acceleration [L T-2 ~> m s-2]
  type(PressureForce_CS),  target, intent(inout) :: CS !< Pressure force control structure
  type(ALE_CS),            pointer     :: ALE_CSp !< ALE control structure
  type(accel_diag_ptrs),   pointer     :: ADp !< Acceleration diagnostic pointers
  real, dimension(:,:),    pointer     :: p_atm !< The pressure at the ice-ocean or
                                               !! atmosphere-ocean interface [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out) :: pbce !< The baroclinic pressure anomaly in each layer
                                               !! due to eta anomalies [L2 T-2 H-1 ~> m s-2 or m4 s-2 kg-1].
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(out) :: eta  !< The bottom mass used to calculate PFu and PFv,
                                               !! [H ~> m or kg m-2], with any tidal contributions.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    e_sal, &          ! The bottom geopotential anomaly due to self-attraction and loading [Z ~> m].
    e_tidal_eq, &     ! The bottom geopotential anomaly due to tidal forces from astronomical sources [Z ~> m].
    e_tidal_sal, &    ! The bottom geopotential anomaly due to harmonic self-attraction and loading
                      ! specific to tides [Z ~> m].
    e_sal_and_tide    ! The summation of self-attraction and loading and tidal forcing [Z ~> m].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    MassWt_u, &       ! The fractional mass weighting at a u-point [nondim].
    intx_dpa          ! The change in intx_pa through a layer [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    MassWt_v, &       ! The fractional mass weighting at a v-point [nondim].
    inty_dpa          ! The change in inty_pa through a layer [R L2 T-2 ~> Pa].
  type(PF_FV_host_ctx),   target :: FV_ctx   ! The objects used by the finite volume form's host helpers
  type(PF_Mont_host_ctx), target :: Mont_ctx ! The objects used by the Montgomery potential form's host helpers
  type(c_ptr) :: host ! Opaque pointer to the host context of the form in use
  logical :: use_EOS ! If true, density is calculated from T & S using an equation of state.
  type(RealArray_t) :: h_a, PFu_a, PFv_a, p_atm_a, e_sal_a, e_tidal_eq_a, e_tidal_sal_a, e_sal_and_tide_a
  type(RealArray_t) :: MassWt_u_a, MassWt_v_a, intx_dpa_a, inty_dpa_a, pbce_a, eta_a

  if (CS%Analytic_FV_PGF) then
    FV_ctx%G => G ; FV_ctx%GV => GV ; FV_ctx%US => US ; FV_ctx%tv => tv
    FV_ctx%ALE_CSp => ALE_CSp ; FV_ctx%CS => CS%PressureForce_FV
    FV_ctx%pbce%G => G ; FV_ctx%pbce%GV => GV ; FV_ctx%pbce%US => US
    host = c_loc(FV_ctx)
  else
    Mont_ctx%G => G ; Mont_ctx%GV => GV ; Mont_ctx%US => US ; Mont_ctx%tv => tv
    Mont_ctx%CS => CS%PressureForce_Mont
    Mont_ctx%pbce%G => G ; Mont_ctx%pbce%GV => GV ; Mont_ctx%pbce%US => US
    host = c_loc(Mont_ctx)
  endif

  use_EOS = associated(tv%eqn_of_state)

  if (CS%Analytic_FV_PGF) then
    if (GV%Boussinesq) then
      call PressureForce_FV_Bouss_pre(G, US, CS%PressureForce_FV, e_tidal_eq, e_tidal_sal)
    else
      call PressureForce_FV_nonBouss_pre(G, GV, US, CS%PressureForce_FV, e_tidal_eq, e_tidal_sal)
    endif
  else
    if (GV%Boussinesq) then
      call PressureForce_Mont_Bouss_pre(G, US, CS%PressureForce_Mont, use_EOS, tv, e_tidal_eq, e_tidal_sal)
    else
      call PressureForce_Mont_nonBouss_pre(G, US, CS%PressureForce_Mont, use_EOS, tv, e_tidal_eq, e_tidal_sal)
    endif
  endif

  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call PFu_a%alloc(lb=LBOUND(PFu), ub=UBOUND(PFu), source=PFu)
  call PFv_a%alloc(lb=LBOUND(PFv), ub=UBOUND(PFv), source=PFv)
  if (associated(p_atm)) call p_atm_a%alloc(lb=LBOUND(p_atm), ub=UBOUND(p_atm), source=p_atm)
  call e_sal_a%alloc(lb=LBOUND(e_sal), ub=UBOUND(e_sal), source=e_sal)
  call e_tidal_eq_a%alloc(lb=LBOUND(e_tidal_eq), ub=UBOUND(e_tidal_eq), source=e_tidal_eq)
  call e_tidal_sal_a%alloc(lb=LBOUND(e_tidal_sal), ub=UBOUND(e_tidal_sal), source=e_tidal_sal)
  call e_sal_and_tide_a%alloc(lb=LBOUND(e_sal_and_tide), ub=UBOUND(e_sal_and_tide), source=e_sal_and_tide)
  call MassWt_u_a%alloc(lb=LBOUND(MassWt_u), ub=UBOUND(MassWt_u), source=MassWt_u)
  call MassWt_v_a%alloc(lb=LBOUND(MassWt_v), ub=UBOUND(MassWt_v), source=MassWt_v)
  call intx_dpa_a%alloc(lb=LBOUND(intx_dpa), ub=UBOUND(intx_dpa), source=intx_dpa)
  call inty_dpa_a%alloc(lb=LBOUND(inty_dpa), ub=UBOUND(inty_dpa), source=inty_dpa)
  if (present(pbce)) call pbce_a%alloc(lb=LBOUND(pbce), ub=UBOUND(pbce), source=pbce)
  if (present(eta)) call eta_a%alloc(lb=LBOUND(eta), ub=UBOUND(eta), source=eta)

  call PressureForce_TR(h_a, use_EOS, PFu_a, PFv_a, G, GV, US, CS, associated(ALE_CSp), p_atm_a, host, e_sal_a, &
                        e_tidal_eq_a, e_tidal_sal_a, e_sal_and_tide_a, MassWt_u_a, MassWt_v_a, intx_dpa_a, &
                        inty_dpa_a, pbce_a, eta_a)

  call PFu_a%copy2F(PFu)
  call PFv_a%copy2F(PFv)
  call e_sal_a%copy2F(e_sal)
  call e_tidal_eq_a%copy2F(e_tidal_eq)
  call e_tidal_sal_a%copy2F(e_tidal_sal)
  call e_sal_and_tide_a%copy2F(e_sal_and_tide)
  call MassWt_u_a%copy2F(MassWt_u)
  call MassWt_v_a%copy2F(MassWt_v)
  call intx_dpa_a%copy2F(intx_dpa)
  call inty_dpa_a%copy2F(inty_dpa)
  if (present(pbce)) call pbce_a%copy2F(pbce)
  if (present(eta)) call eta_a%copy2F(eta)
  call h_a%free() ; call PFu_a%free() ; call PFv_a%free() ; call p_atm_a%free()
  call e_sal_a%free() ; call e_tidal_eq_a%free() ; call e_tidal_sal_a%free() ; call e_sal_and_tide_a%free()
  call MassWt_u_a%free() ; call MassWt_v_a%free() ; call intx_dpa_a%free() ; call inty_dpa_a%free()
  call pbce_a%free() ; call eta_a%free()

  if (CS%Analytic_FV_PGF) then
    if (GV%Boussinesq) then
      call PressureForce_FV_Bouss_post(h, tv, G, GV, CS%PressureForce_FV, ADp, p_atm, e_sal, e_tidal_eq, &
                                       e_tidal_sal, e_sal_and_tide, MassWt_u, MassWt_v, intx_dpa, inty_dpa)
    else
      call PressureForce_FV_nonBouss_post(G, GV, CS%PressureForce_FV, ADp, e_sal, e_tidal_eq, e_tidal_sal, &
                                          e_sal_and_tide, MassWt_u, MassWt_v)
    endif
  else
    if (GV%Boussinesq) then
      call PressureForce_Mont_Bouss_post(G, CS%PressureForce_Mont, e_sal, e_tidal_eq, e_tidal_sal)
    else
      call PressureForce_Mont_nonBouss_post(G, CS%PressureForce_Mont, e_sal, e_tidal_eq, e_tidal_sal)
    endif
  endif

end subroutine PressureForce

!> The implementation of PressureForce (the root of its call tree): calls the finite volume or
!! Montgomery potential form of the pressure force, in its Boussinesq or non-Boussinesq version.
subroutine PressureForce_TR(h_a, use_EOS, PFu_a, PFv_a, G, GV, US, CS, have_ALE_CSp, p_atm_a, host, e_sal_a, &
                            e_tidal_eq_a, e_tidal_sal_a, e_sal_and_tide_a, MassWt_u_a, MassWt_v_a, intx_dpa_a, &
                            inty_dpa_a, pbce_a, eta_a)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  type(RealArray_t),       intent(in)  :: h_a !< Layer thicknesses [H ~> m or kg m-2]
  logical,                 intent(in)  :: use_EOS !< If true, density is calculated from T & S using an
                                                  !! equation of state.
  type(RealArray_t),       intent(inout) :: PFu_a !< Zonal pressure force acceleration [L T-2 ~> m s-2]
  type(RealArray_t),       intent(inout) :: PFv_a !< Meridional pressure force acceleration [L T-2 ~> m s-2]
  type(PressureForce_CS),  intent(inout) :: CS !< Pressure force control structure
  logical,                 intent(in)  :: have_ALE_CSp !< If true, the ALE control structure (ALE_CSp)
                                                       !! is associated.
  type(RealArray_t),       intent(in)  :: p_atm_a !< The pressure at the ice-ocean or
                                                  !! atmosphere-ocean interface [R L2 T-2 ~> Pa].
  type(c_ptr),      value, intent(in)  :: host !< Opaque pointer to the PF_FV_host_ctx or PF_Mont_host_ctx
                                               !! that holds the Fortran objects used by the host helpers
  type(RealArray_t),       intent(inout) :: e_sal_a !< The bottom geopotential anomaly due to self-attraction and
                                                    !! loading [Z ~> m].
  type(RealArray_t),       intent(inout) :: e_tidal_eq_a !< The bottom geopotential anomaly due to tidal forces
                                                         !! from astronomical sources [Z ~> m].
  type(RealArray_t),       intent(inout) :: e_tidal_sal_a !< The bottom geopotential anomaly due to harmonic
                                                          !! self-attraction and loading specific to tides
                                                          !! [Z ~> m].
  type(RealArray_t),       intent(inout) :: e_sal_and_tide_a !< The summation of self-attraction and loading and
                                                             !! tidal forcing [Z ~> m].
  type(RealArray_t),       intent(inout) :: MassWt_u_a !< The fractional mass weighting at a u-point [nondim].
  type(RealArray_t),       intent(inout) :: MassWt_v_a !< The fractional mass weighting at a v-point [nondim].
  type(RealArray_t),       intent(inout) :: intx_dpa_a !< The change in intx_pa through a layer [R L2 T-2 ~> Pa].
  type(RealArray_t),       intent(inout) :: inty_dpa_a !< The change in inty_pa through a layer [R L2 T-2 ~> Pa].
  type(RealArray_t),       intent(inout) :: pbce_a !< The baroclinic pressure anomaly in each layer
                                                   !! due to eta anomalies [L2 T-2 H-1 ~> m s-2 or m4 s-2 kg-1].
  type(RealArray_t),       intent(inout) :: eta_a  !< The bottom mass used to calculate PFu and PFv,
                                                   !! [H ~> m or kg m-2], with any tidal contributions.
  ! Local variables
  type(RealArray_t) :: PFu_bc3d_a, PFv_bc3d_a

  if (CS%Analytic_FV_PGF) then
    if (GV%Boussinesq) then
      call PressureForce_FV_Bouss(h_a, use_EOS, PFu_a, PFv_a, G, GV, US, CS%PressureForce_FV%PF_FV_opts_CS, &
                                  CS%PressureForce_FV%PF_vgrid_consts_CS, have_ALE_CSp, p_atm_a, host, e_sal_a, &
                                  e_tidal_eq_a, e_tidal_sal_a, e_sal_and_tide_a, MassWt_u_a, MassWt_v_a, &
                                  intx_dpa_a, inty_dpa_a, pbce_a, eta_a)
    else
      call PressureForce_FV_nonBouss(h_a, use_EOS, PFu_a, PFv_a, G, GV, US, CS%PressureForce_FV%PF_FV_opts_CS, &
                                     CS%PressureForce_FV%PF_vgrid_consts_CS, have_ALE_CSp, p_atm_a, host, &
                                     e_sal_a, e_tidal_eq_a, e_tidal_sal_a, e_sal_and_tide_a, MassWt_u_a, &
                                     MassWt_v_a, pbce_a, eta_a)
    endif
  else
    if (allocated(CS%PressureForce_Mont%PFu_bc)) &
      call PFu_bc3d_a%alloc(lb=LBOUND(CS%PressureForce_Mont%PFu_bc), ub=UBOUND(CS%PressureForce_Mont%PFu_bc), &
                            source=CS%PressureForce_Mont%PFu_bc)
    if (allocated(CS%PressureForce_Mont%PFv_bc)) &
      call PFv_bc3d_a%alloc(lb=LBOUND(CS%PressureForce_Mont%PFv_bc), ub=UBOUND(CS%PressureForce_Mont%PFv_bc), &
                            source=CS%PressureForce_Mont%PFv_bc)
    if (GV%Boussinesq) then
      call PressureForce_Mont_Bouss(h_a, use_EOS, PFu_a, PFv_a, G, GV, US, &
                                    CS%PressureForce_Mont%PF_Mont_opts_CS, &
                                    CS%PressureForce_Mont%PF_vgrid_consts_CS, p_atm_a, host, e_sal_a, &
                                    e_tidal_eq_a, e_tidal_sal_a, pbce_a, eta_a, PFu_bc3d_a, PFv_bc3d_a)
    else
      call PressureForce_Mont_nonBouss(h_a, use_EOS, PFu_a, PFv_a, G, GV, US, &
                                       CS%PressureForce_Mont%PF_Mont_opts_CS, &
                                       CS%PressureForce_Mont%PF_vgrid_consts_CS, p_atm_a, host, e_sal_a, &
                                       e_tidal_eq_a, e_tidal_sal_a, pbce_a, eta_a, PFu_bc3d_a, PFv_bc3d_a)
    endif
    if (allocated(CS%PressureForce_Mont%PFu_bc)) call PFu_bc3d_a%copy2F(CS%PressureForce_Mont%PFu_bc)
    if (allocated(CS%PressureForce_Mont%PFv_bc)) call PFv_bc3d_a%copy2F(CS%PressureForce_Mont%PFv_bc)
    call PFu_bc3d_a%free() ; call PFv_bc3d_a%free()
  endif

end subroutine PressureForce_TR

!> Initialize the pressure force control structure
subroutine PressureForce_init(Time, G, GV, US, param_file, diag, CS, ADp, SAL_CSp, tides_CSp)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(PressureForce_CS),  intent(inout) :: CS   !< Pressure force control structure
  type(accel_diag_ptrs),   pointer       :: ADp !< Acceleration diagnostic pointers
  type(SAL_CS),           intent(in), optional :: SAL_CSp !< SAL control structure
  type(tidal_forcing_CS), intent(in), optional :: tides_CSp !< Tide control structure
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_PressureForce" ! This module's name.

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ANALYTIC_FV_PGF", CS%Analytic_FV_PGF, &
                 "If true the pressure gradient forces are calculated "//&
                 "with a finite volume form that analytically integrates "//&
                 "the equations of state in pressure to avoid any "//&
                 "possibility of numerical thermobaric instability, as "//&
                 "described in Adcroft et al., O. Mod. (2008).", default=.true.)

  if (CS%Analytic_FV_PGF) then
    !$omp target enter data map(alloc: CS%PressureForce_FV)
    call PressureForce_FV_init(Time, G, GV, US, param_file, diag, &
             CS%PressureForce_FV, ADp, SAL_CSp, tides_CSp)
  else
    call PressureForce_Mont_init(Time, G, GV, US, param_file, diag, &
             CS%PressureForce_Mont, SAL_CSp, tides_CSp)
  endif
  !$omp target update to(CS)
end subroutine PressureForce_init

!> \namespace mom_pressureforce
!!
!! This thin module provides a branch to two forms of the horizontal accelerations
!! due to pressure gradients. The two options currently available are a
!! Montgomery potential form (used in traditional isopycnal layer models), and the
!! analytic finite volume form.

end module MOM_PressureForce
