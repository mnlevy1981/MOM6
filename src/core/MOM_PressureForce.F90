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
use MOM_PressureForce_Mont, only : PressureForce_Mont_Bouss, PressureForce_Mont_nonBouss
use MOM_PressureForce_Mont, only : PressureForce_Mont_init
use MOM_PressureForce_Mont, only : PressureForce_Mont_CS
use MOM_self_attr_load, only : SAL_CS
use MOM_tidal_forcing, only : tidal_forcing_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs, accel_diag_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_ALE, only: ALE_CS
use MOM_grid_containers, only : grid_core_type, grid_core
use array_mod, only : RealArray_t
use box_mod, only : Box_t
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
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)  :: tv   !< A structure pointing to various thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: PFu  !< Zonal pressure force acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out) :: PFv  !< Meridional pressure force acceleration [L T-2 ~> m s-2]
  type(PressureForce_CS),  intent(inout) :: CS !< Pressure force control structure
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

  type(RealArray_t) :: h_a, PFu_a, PFv_a, p_atm_a, pbce_a, eta_a, Rlay_a, g_prime_a
  type(grid_core_type), pointer :: Gcore ! Persistent copies of the core ocean grid fields
  type(Box_t) :: bxH0  ! The h-point box of the whole column, [isc:iec, jsc:jec, 1:nz]
                       ! with no halo. Set_pbce_Bouss and Set_pbce_nonBouss read its bounds,
                       ! PressureForce_Mont_nonBouss and Set_pbce_nonBouss iterate its k range,
                       ! and the forms pass it through.
                       ! bxHkint is grown from it.
  type(Box_t) :: bxQ0  ! The B-grid-index box of the whole column, [IscB:IecB, JscB:JecB, 1:nz]
                       ! with no halo. bxHpgf is grown from it; Set_pbce_Bouss and
                       ! Set_pbce_nonBouss read its bounds; the forms pass it through.
  type(Box_t) :: bxU0  ! The u-point box of the whole column, [IscB:IecB, jsc:jec, 1:nz]
                       ! with no halo, iterated by PressureForce_FV_Bouss and
                       ! PressureForce_FV_nonBouss.
  type(Box_t) :: bxV0  ! The v-point box of the whole column, [isc:iec, JscB:JecB, 1:nz]
                       ! with no halo, iterated by PressureForce_FV_Bouss and
                       ! PressureForce_FV_nonBouss.
  type(Box_t) :: bxHpgf ! The h points the u and v pressure gradients read,
                        ! [IscB:IecB+1, JscB:JecB+1, 1:nz]. A gradient at the last face, IecB
                        ! or JecB, reads h at i+1 or j+1, so the range ends one point beyond the
                        ! compute domain; it is built from bxQ0 because it starts at IscB and
                        ! JscB. Iterated by the four PressureForce forms and Set_pbce_Bouss.
  type(Box_t) :: bxHkint ! The interior interfaces, [isc:iec, jsc:jec, 2:nz]: a jump across
                         ! interface K takes layers k-1 and k, so k starts at 2. Grown from
                         ! bxH0. Iterated by PressureForce_Mont_nonBouss and Set_pbce_nonBouss;
                         ! PressureForce_FV_nonBouss passes it through.

  Gcore => grid_core(G)
  call bxH0%safe_alloc(ndims=3) ; call bxQ0%safe_alloc(ndims=3)
  call bxU0%safe_alloc(ndims=3) ; call bxV0%safe_alloc(ndims=3)
  call bxH0%set(idxS=[G%isc,G%jsc,1], idxE=[G%iec,G%jec,GV%ke])
  call bxQ0%set(idxS=[G%IscB,G%JscB,1], idxE=[G%IecB,G%JecB,GV%ke])
  call bxU0%set(idxS=[G%IscB,G%jsc,1], idxE=[G%IecB,G%jec,GV%ke])
  call bxV0%set(idxS=[G%isc,G%JscB,1], idxE=[G%iec,G%JecB,GV%ke])
  bxHpgf = bxQ0%growBy(lo=[0,0,0], hi=[1,1,0])
  bxHkint = bxH0%growLo(dim=3, n=-1)

  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call PFu_a%alloc(lb=LBOUND(PFu), ub=UBOUND(PFu), source=PFu)
  call PFv_a%alloc(lb=LBOUND(PFv), ub=UBOUND(PFv), source=PFv)
  if (associated(p_atm)) call p_atm_a%alloc(lb=LBOUND(p_atm), ub=UBOUND(p_atm), source=p_atm)
  if (present(pbce)) call pbce_a%alloc(lb=LBOUND(pbce), ub=UBOUND(pbce), source=pbce)
  if (present(eta)) call eta_a%alloc(lb=LBOUND(eta), ub=UBOUND(eta), source=eta)
  call Rlay_a%alloc(lb=LBOUND(GV%Rlay), ub=UBOUND(GV%Rlay), source=GV%Rlay)
  call g_prime_a%alloc(lb=LBOUND(GV%g_prime), ub=UBOUND(GV%g_prime), source=GV%g_prime)

  call PressureForce_TR(bxH0, bxQ0, bxU0, bxV0, bxHpgf, bxHkint, h_a, tv, PFu_a, PFv_a, G, Gcore, &
                        GV, US, CS, ALE_CSp, ADp, p_atm_a, pbce_a, eta_a, Rlay_a, g_prime_a, &
                        G%Z_ref, GV%g_Earth, GV%H_to_Z, GV%Z_to_H, GV%H_to_RZ, GV%Rho0, &
                        GV%H_subroundoff, GV%dZ_subroundoff, GV%nk_rho_varies, GV%semi_Boussinesq, &
                        GV%Boussinesq)
  call bxH0%free() ; call bxQ0%free()
  call bxU0%free() ; call bxV0%free()
  call bxHpgf%free()
  call bxHkint%free()

  call PFu_a%copy2F(PFu)
  call PFv_a%copy2F(PFv)
  if (present(pbce)) call pbce_a%copy2F(pbce)
  if (present(eta)) call eta_a%copy2F(eta)
  call h_a%free() ; call PFu_a%free() ; call PFv_a%free() ; call p_atm_a%free()
  call pbce_a%free() ; call eta_a%free() ; call Rlay_a%free() ; call g_prime_a%free()

end subroutine PressureForce

!> The implementation of PressureForce (the root of its call tree).
subroutine PressureForce_TR(bxH0, bxQ0, bxU0, bxV0, bxHpgf, bxHkint, h_a, tv, PFu_a, PFv_a, G, &
                            Gcore, GV, US, CS, ALE_CSp, ADp, p_atm_a, pbce_a, eta_a, Rlay_a, &
                            g_prime_a, Z_ref, g_Earth, H_to_Z, Z_to_H, H_to_RZ, Rho0, &
                            H_subroundoff, dZ_subroundoff, nk_rho_varies, semi_Boussinesq, &
                            Boussinesq)
  type(Box_t),             intent(in)  :: bxH0 !< The h-point box of the whole column,
                                               !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),             intent(in)  :: bxQ0 !< The B-grid-index box of the whole column,
                                               !! [IscB:IecB, JscB:JecB, 1:nz]
  type(Box_t),             intent(in)  :: bxU0 !< The u-point box of the whole column,
                                               !! [IscB:IecB, jsc:jec, 1:nz]
  type(Box_t),             intent(in)  :: bxV0 !< The v-point box of the whole column,
                                               !! [isc:iec, JscB:JecB, 1:nz]
  type(Box_t),             intent(in)  :: bxHpgf !< The h points the u and v pressure gradients
                                                 !! read, [IscB:IecB+1, JscB:JecB+1, 1:nz]
  type(Box_t),             intent(in)  :: bxHkint !< The interior interfaces of the whole column,
                                                  !! [isc:iec, jsc:jec, 2:nz]
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(grid_core_type),    intent(in)  :: Gcore !< Persistent copies of the core ocean grid fields
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  type(RealArray_t),       intent(in)  :: h_a  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)  :: tv   !< A structure pointing to various thermodynamic
                                               !! variables
  type(RealArray_t), intent(inout) :: PFu_a !< Zonal pressure force acceleration [L T-2 ~> m s-2]
  type(RealArray_t), intent(inout) :: PFv_a !< Meridional pressure force acceleration
                                            !! [L T-2 ~> m s-2]
  type(PressureForce_CS),  intent(inout) :: CS !< Pressure force control structure
  type(ALE_CS),            pointer     :: ALE_CSp !< ALE control structure
  type(accel_diag_ptrs),   pointer     :: ADp !< Acceleration diagnostic pointers
  type(RealArray_t),       intent(in)  :: p_atm_a !< The pressure at the ice-ocean or
                                               !! atmosphere-ocean interface [R L2 T-2 ~> Pa].
  type(RealArray_t), &
                 intent(inout) :: pbce_a !< The baroclinic pressure anomaly in each layer
                                               !! due to eta anomalies
                                               !! [L2 T-2 H-1 ~> m s-2 or m4 s-2 kg-1].
  type(RealArray_t), &
                 intent(inout) :: eta_a  !< The bottom mass used to calculate PFu and PFv,
                                               !! [H ~> m or kg m-2], with any tidal contributions.
  type(RealArray_t), intent(in) :: Rlay_a !< The target coordinate value (potential density) in each
                                          !! layer [R ~> kg m-3].
  type(RealArray_t), intent(in) :: g_prime_a !< The reduced gravity at each interface
                                             !! [L2 Z-1 T-2 ~> m s-2].
  real,              intent(in) :: Z_ref !< A reference value for all geometric height fields, such
                                         !! as bathyT [Z ~> m].
  real,              intent(in) :: g_Earth !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real,              intent(in) :: H_to_Z !< A constant that translates thickness units to the units
                                          !! of depth [Z H-1 ~> 1 or m3 kg-1].
  real,              intent(in) :: Z_to_H !< A constant that translates depth units to thickness
                                          !! units depth [H Z-1 ~> 1 or kg m-3].
  real,              intent(in) :: H_to_RZ !< A constant that translates thickness units to the
                                           !! units of mass per unit area [R Z H-1 ~> kg m-3 or 1].
  real,              intent(in) :: Rho0 !< The density used in the Boussinesq approximation or
                                        !! nominal density used to convert depths into mass units
                                        !! [R ~> kg m-3].
  real,              intent(in) :: H_subroundoff !< A thickness that is so small that it can be
                                                 !! added to a thickness of Angstrom or larger
                                                 !! without changing it at the bit level
                                                 !! [H ~> m or kg m-2]. If Angstrom is 0 or
                                                 !! exceedingly small, this is negligible compared
                                                 !! to 1e-17 m.
  real,              intent(in) :: dZ_subroundoff !< A thickness in height units that is so small
                                                  !! that it can be added to a vertical distance of
                                                  !! Angstrom_Z or 1e-17 m without changing it at
                                                  !! the bit level [Z ~> m].  This is the height
                                                  !! equivalent of H_subroundoff.
  integer,           intent(in) :: nk_rho_varies !< The number of layers at the top where the
                                                 !! density does not track any target density.
  logical,           intent(in) :: semi_Boussinesq !< If true, do non-Boussinesq pressure force
                                                   !! calculations and use mass-based "thicknesses,
                                                   !! but use Rho0 to convert layer thicknesses into
                                                   !! certain height changes.  This only applies if
                                                   !! BOUSSINESQ is false.
  logical,           intent(in) :: Boussinesq !< If true, make the Boussinesq approximation.

  if (CS%Analytic_FV_PGF) then
    if (Boussinesq) then
      call PressureForce_FV_Bouss(bxH0, bxQ0, bxU0, bxV0, bxHpgf, h_a, tv, PFu_a, PFv_a, G, Gcore, &
                                  GV, US, CS%PressureForce_FV, ALE_CSp, ADp, p_atm_a, pbce_a, &
                                  eta_a, Rlay_a, g_prime_a, Z_ref, g_Earth, H_to_Z, Z_to_H, &
                                  H_to_RZ, Rho0, H_subroundoff, dZ_subroundoff, nk_rho_varies)
    else
      call PressureForce_FV_nonBouss(bxH0, bxQ0, bxU0, bxV0, bxHpgf, bxHkint, h_a, tv, PFu_a, &
                                     PFv_a, G, Gcore, GV, US, CS%PressureForce_FV, ALE_CSp, ADp, &
                                     p_atm_a, pbce_a, eta_a, Rlay_a, Z_ref, g_Earth, H_to_RZ, &
                                     H_subroundoff, nk_rho_varies, semi_Boussinesq)
    endif
  else
    if (Boussinesq) then
      call PressureForce_Mont_Bouss(bxH0, bxQ0, bxHpgf, h_a, tv, PFu_a, PFv_a, G, Gcore, US, &
                                    CS%PressureForce_Mont, p_atm_a, pbce_a, eta_a, Rlay_a, &
                                    g_prime_a, g_Earth, H_to_Z, Z_to_H, Rho0, dZ_subroundoff, &
                                    nk_rho_varies)
    else
      call PressureForce_Mont_nonBouss(bxH0, bxQ0, bxHpgf, bxHkint, h_a, tv, PFu_a, PFv_a, G, &
                                       Gcore, US, CS%PressureForce_Mont, p_atm_a, pbce_a, &
                                       eta_a, Rlay_a, g_Earth, H_to_RZ, H_subroundoff, &
                                       nk_rho_varies)
    endif
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
