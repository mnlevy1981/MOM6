! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!!SKILLS: 0.4.0

!> Provides the Montgomery potential form of pressure gradient
module MOM_PressureForce_Mont

use MOM_density_integrals, only : int_specific_vol_dp
use MOM_diag_mediator, only : post_data, register_diag_field
use MOM_diag_mediator, only : safe_alloc_ptr, diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_grid_containers, only : grid_core_type
use MOM_self_attr_load, only : calc_SAL, SAL_CS
use MOM_tidal_forcing, only : calc_tidal_forcing, tidal_forcing_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs
use MOM_EOS, only : query_compressible
use array_mod, only : RealArray_t
use box_mod, only : Box_t

implicit none ; private

#include <MOM_memory.h>

public PressureForce_Mont_Bouss, PressureForce_Mont_nonBouss, Set_pbce_Bouss
public Set_pbce_nonBouss, PressureForce_Mont_init

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for the Montgomery potential form of pressure gradient
type, public :: PressureForce_Mont_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: calculate_SAL  !< If true, calculate self-attraction and loading.
  logical :: tides          !< If true, apply tidal momentum forcing.
  real    :: Rho0           !< The density used in the Boussinesq
                            !! approximation [R ~> kg m-3].
  real    :: GFS_scale      !< Ratio between gravity applied to top interface and the
                            !! gravitational acceleration of the planet [nondim].
                            !! Usually this ratio is 1.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate
                            !! the timing of diagnostic output.
  real, allocatable :: PFu_bc(:,:,:) !< Zonal accelerations due to pressure gradients
                            !! deriving from density gradients within layers [L T-2 ~> m s-2].
  real, allocatable :: PFv_bc(:,:,:) !< Meridional accelerations due to pressure gradients
                            !! deriving from density gradients within layers [L T-2 ~> m s-2].
  !>@{ Diagnostic IDs
  integer :: id_PFu_bc = -1, id_PFv_bc = -1, id_e_sal = -1
  integer :: id_e_tide = -1,  id_e_tide_eq = -1, id_e_tide_sal = -1
  !>@}
  type(SAL_CS), pointer :: SAL_CSp => NULL() !< SAL control structure
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL() !< The tidal forcing control structure
end type PressureForce_Mont_CS

contains

!> \brief Non-Boussinesq Montgomery-potential form of pressure gradient
!!
!! Determines the acceleration due to pressure forces in a
!! non-Boussinesq fluid using the compressibility compensated (if appropriate)
!! Montgomery-potential form described in Hallberg (Ocean Mod., 2005).
!!
!! To work, the following fields must be set outside of the usual (is:ie,js:je)
!! range before this subroutine is called:
!!   h(isB:ie+1,jsB:je+1), T(isB:ie+1,jsB:je+1), and S(isB:ie+1,jsB:je+1).
subroutine PressureForce_Mont_nonBouss(bxH0, bxQ0, bxHpgf, bxHkint, h_a, tv, PFu_a, PFv_a, G, &
                                       Gcore, US, CS, p_atm_a, pbce_a, eta_a, Rlay_a, g_Earth, &
                                       H_to_RZ, H_subroundoff, nk_rho_varies)
  type(Box_t),          intent(in)  :: bxH0  !< The h-point box of the whole column,
                                             !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),          intent(in)  :: bxQ0  !< The B-grid-index box of the whole column,
                                             !! [IscB:IecB, JscB:JecB, 1:nz]
  type(Box_t),          intent(in)  :: bxHpgf !< The h points the u and v pressure gradients
                                              !! read, [IscB:IecB+1, JscB:JecB+1, 1:nz]
  type(Box_t),          intent(in)  :: bxHkint !< The interior interfaces of the whole column,
                                               !! [isc:iec, jsc:jec, 2:nz]
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(grid_core_type), intent(in)  :: Gcore !< Persistent copies of the core ocean grid fields
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  type(RealArray_t), intent(in) :: h_a !< Layer thickness, [H ~> kg m-2].
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables.
  type(RealArray_t), intent(inout) :: PFu_a !< Zonal acceleration due to pressure gradients
                                                                 !! (equal to -dM/dx) [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: PFv_a !< Meridional acceleration due to pressure gradients
                                                                 !! (equal to -dM/dy) [L T-2 ~> m s-2].
  type(PressureForce_Mont_CS),                intent(inout) :: CS  !< Control structure for Montgomery potential PGF
  type(RealArray_t), intent(in) :: p_atm_a !< The pressure at the ice-ocean or
                                                                 !! atmosphere-ocean [R L2 T-2 ~> Pa].
  type(RealArray_t), intent(inout) :: pbce_a !< The baroclinic pressure anomaly in
                                                                 !! each layer due to free surface height anomalies,
                                                                 !! [L2 T-2 H-1 ~> m s-2 or m4 kg-1 s-2].
  type(RealArray_t), intent(inout) :: eta_a !< The total column mass used to calculate
                                                                 !! PFu and PFv [H ~> kg m-2].
  type(RealArray_t), intent(in) :: Rlay_a !< The target coordinate value (potential density) in each
                                          !! layer [R ~> kg m-3].
  real,              intent(in) :: g_Earth !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real,              intent(in) :: H_to_RZ !< A constant that translates thickness units to the
                                           !! units of mass per unit area [R Z H-1 ~> kg m-3 or 1].
  real,              intent(in) :: H_subroundoff !< A thickness that is so small that it can be
                                                 !! added to a thickness of Angstrom or larger
                                                 !! without changing it at the bit level
                                                 !! [H ~> m or kg m-2]. If Angstrom is 0 or
                                                 !! exceedingly small, this is negligible compared
                                                 !! to 1e-17 m.
  integer,           intent(in) :: nk_rho_varies !< The number of layers at the top where the
                                                 !! density does not track any target density.
  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: h, PFu, PFv, pbce
  real, dimension(:,:), contiguous, pointer :: p_atm, eta
  real, dimension(:,:), contiguous, pointer :: bathyT, meanSL, IdxCu, IdyCv
  real, dimension(:), contiguous, pointer :: Rlay
  type(RealArray_t) :: p_a, alpha_star_a
  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2),h_a%lb(3):h_a%ub(3)) :: &
    M, &          ! The Montgomery potential, M = (p/rho + gz)  [L2 T-2 ~> m2 s-2].
    dz_geo        !   The change in geopotential across a layer [L2 T-2 ~> m2 s-2].
  real, dimension(:,:,:), contiguous, pointer :: &
    alpha_star    ! Compression adjusted specific volume [R-1 ~> m3 kg-1].
  real, dimension(:,:,:), contiguous, pointer :: p ! Interface pressure [R L2 T-2 ~> Pa].
                ! p may be adjusted (with a nonlinear equation of state) so that
                ! its derivative compensates for the adiabatic compressibility
                ! in seawater, but p will still be close to the pressure.
  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2),h_a%lb(3):h_a%ub(3)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [C ~> degC].
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [S ~> ppt].

  real, dimension(h_a%lb(1):h_a%ub(1)) :: Rho_cv_BL  ! The coordinate potential density in the
                  ! deepest variable density near-surface layer [R ~> kg m-3].

  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2)) :: &
    dM, &         !   A barotropic correction to the Montgomery potentials to enable the use
                  ! of a reduced gravity form of the equations [L2 T-2 ~> m2 s-2].
    dp_star, &    ! Layer thickness after compensation for compressibility [R L2 T-2 ~> Pa].
    SSH, &        ! The sea surface height anomaly, in depth units [Z ~> m].
    e_sal, &      ! Bottom geopotential anomaly due to self-attraction and loading [Z ~> m].
    e_tide_eq,  & ! Bottom geopotential anomaly due to tidal forces from astronomical sources [Z ~> m].
    e_tide_sal, & ! Bottom geopotential anomaly due to harmonic self-attraction and loading
                  ! specific to tides [Z ~> m].
    geopot_bot    !   Bottom geopotential relative to a temporally fixed reference value,
                  ! including any tidal contributions [L2 T-2 ~> m2 s-2].
  real :: p_ref(h_a%lb(1):h_a%ub(1))     !   The pressure used to calculate the coordinate
                             ! density [R L2 T-2 ~> Pa] (usually 2e7 Pa = 2000 dbar).
  real :: rho_in_situ(h_a%lb(1):h_a%ub(1)) !In-situ density of a layer [R ~> kg m-3].
  real :: PFu_bc, PFv_bc     ! The pressure gradient force due to along-layer
                             ! compensated density gradients [L T-2 ~> m s-2]
  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [R L2 T-2 ~> Pa].
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_EOS         ! If true, density is calculated from T & S using an equation of state.
  logical :: is_split        ! A flag indicating whether the pressure gradient terms are to be
                             ! split into barotropic and baroclinic pieces.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.

  real :: I_gEarth           ! The inverse of g_Earth [T2 Z L-2 ~> s2 m-1]
!  real :: dalpha
  real :: Pa_to_H     ! A factor to convert from R L2 T-2 to the thickness units (H)
                      ! [H T2 R-1 L-2 ~> m2 s2 kg-1 or s2 m-1].
  real :: alpha_Lay(h_a%lb(3):h_a%ub(3)) ! The specific volume of each layer [R-1 ~> m3 kg-1].
  real :: dalpha_int(h_a%lb(3):h_a%ub(3)+1) ! The change in specific volume across each
                             ! interface [R-1 ~> m3 kg-1].
  integer, dimension(2) :: EOSdom ! The computational domain for the equation of state
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: i, j, k
  is = bxH0%idxS(1) ; ie = bxH0%idxE(1) ; js = bxH0%idxS(2) ; je = bxH0%idxE(2) ; nz = bxH0%idxE(3)
  nkmb=nk_rho_varies
  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2)
  EOSdom(1) = Isq - (h_a%lb(1)-1) ;  EOSdom(2) = bxH0%idxE(1)+1 - (h_a%lb(1)-1)

  call h_a%view(h) ; call PFu_a%view(PFu) ; call PFv_a%view(PFv)
  call Rlay_a%view(Rlay)
  call Gcore%bathyT%view(bathyT) ; call Gcore%meanSL%view(meanSL)
  call Gcore%IdxCu%view(IdxCu) ; call Gcore%IdyCv%view(IdyCv)
  call p_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),bxH0%idxE(3)+1]) ; call p_a%view(p)
  call alpha_star_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),bxH0%idxE(3)])
  call alpha_star_a%view(alpha_star)
  nullify(p_atm, pbce, eta)
  if (p_atm_a%associated()) call p_atm_a%view(p_atm)
  if (pbce_a%associated()) call pbce_a%view(pbce)
  if (eta_a%associated()) call eta_a%view(eta)

  use_p_atm = associated(p_atm)
  is_split = pbce_a%associated()
  use_EOS = associated(tv%eqn_of_state)

  if (.not.CS%initialized) call MOM_error(FATAL, &
      "MOM_PressureForce_Mont: Module must be initialized before it is used.")

  if (use_EOS) then
    if (query_compressible(tv%eqn_of_state)) call MOM_error(FATAL, &
      "PressureForce_Mont_nonBouss: The Montgomery form of the pressure force "//&
      "can no longer be used with a compressible EOS. Use #define ANALYTIC_FV_PGF.")
  endif

  I_gEarth = 1.0 / g_Earth
  dp_neglect = g_Earth * H_to_RZ * H_subroundoff
  do concurrent (k=bxH0%idxS(3):bxH0%idxE(3)) ; alpha_Lay(k) = 1.0 / (Rlay(k)) ; enddo
  do concurrent (k=bxHkint%idxS(3):bxHkint%idxE(3)) ; dalpha_int(K) = alpha_Lay(k-1) - alpha_Lay(k) ; enddo

  if (use_p_atm) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1 ; p(i,j,1) = p_atm(i,j) ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1 ; p(i,j,1) = 0.0 ; enddo ; enddo
  endif
  !$OMP parallel do default(shared)
  do j=Jsq,Jeq+1 ; do k=1,nz ; do i=Isq,Ieq+1
    p(i,j,K+1) = p(i,j,K) + g_Earth * H_to_RZ * h(i,j,k)
  enddo ; enddo ; enddo

  if (eta_a%associated()) then
    Pa_to_H = 1.0 / (g_Earth * H_to_RZ)
    if (use_p_atm) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = (p(i,j,nz+1) - p_atm(i,j)) * Pa_to_H ! eta has the same units as h.
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = p(i,j,nz+1) * Pa_to_H ! eta has the same units as h.
      enddo ; enddo
    endif
  endif

  !$OMP parallel do default(shared)
  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    geopot_bot(i,j) = -g_Earth * bathyT(i,j)
  enddo ; enddo

  ! Calculate and add the self-attraction and loading geopotential anomaly.
  if (CS%calculate_SAL) then
    !   Determine the sea surface height anomalies, to enable the calculation
    ! of self-attraction and loading.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      SSH(i,j) = min(-bathyT(i,j) - meanSL(i,j), 0.0)
    enddo ; enddo
    if (use_EOS) then
      !$OMP parallel do default(shared)
      do k=1,nz
        call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p(:,:,k), p(:,:,k+1), &
                                 0.0, G%HI, tv%eqn_of_state, US, dz_geo(:,:,k), halo_size=1)
      enddo
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do k=1,nz ; do i=Isq,Ieq+1
        SSH(i,j) = SSH(i,j) + I_gEarth * dz_geo(i,j,k)
      enddo ; enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do k=1,nz ; do i=Isq,Ieq+1
        SSH(i,j) = SSH(i,j) + H_to_RZ * h(i,j,k) * alpha_Lay(k)
      enddo ; enddo ; enddo
    endif

    call calc_SAL(SSH, e_sal, G, CS%SAL_CSp, tmp_scale=US%Z_to_m)
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      geopot_bot(i,j) = geopot_bot(i,j) - g_Earth*e_sal(i,j)
    enddo ; enddo
  endif

  ! Calculate and add the tidal geopotential anomaly.
  if (CS%tides) then
    call calc_tidal_forcing(CS%Time, e_tide_eq, e_tide_sal, G, US, CS%tides_CSp)
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      geopot_bot(i,j) = geopot_bot(i,j) - g_Earth*(e_tide_eq(i,j) + e_tide_sal(i,j))
    enddo ; enddo
  endif

  if (use_EOS) then
    !   Calculate in-situ specific volumes (alpha_star).

    !   With a bulk mixed layer, replace the T & S of any layers that are
    ! lighter than the buffer layer with the properties of the buffer
    ! layer.  These layers will be massless anyway, and it avoids any
    ! formal calculations with hydrostatically unstable profiles.
    if (nkmb>0) then
      tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
      tv_tmp%eqn_of_state => tv%eqn_of_state
      do concurrent (i=bxHpgf%idxS(1):bxHpgf%idxE(1)) ; p_ref(i) = tv%P_Ref ; enddo
      !$OMP parallel do default(shared) private(Rho_cv_BL)
      do j=Jsq,Jeq+1
        do k=1,nkmb ; do i=Isq,Ieq+1
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        enddo ; enddo
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, Rho_cv_BL(:), &
                               tv%eqn_of_state, EOSdom)
        do k=nkmb+1,nz ; do i=Isq,Ieq+1
          if (Rlay(k) < Rho_cv_BL(i)) then
            tv_tmp%T(i,j,k) = tv%T(i,j,nkmb) ; tv_tmp%S(i,j,k) = tv%S(i,j,nkmb)
          else
            tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
          endif
        enddo ; enddo
      enddo
    else
      tv_tmp%T => tv%T ; tv_tmp%S => tv%S
      tv_tmp%eqn_of_state => tv%eqn_of_state
      do concurrent (i=bxHpgf%idxS(1):bxHpgf%idxE(1)) ; p_ref(i) = 0 ; enddo
    endif
    !$OMP parallel do default(shared) private(rho_in_situ)
    do k=1,nz ; do j=Jsq,Jeq+1
      call calculate_density(tv_tmp%T(:,j,k), tv_tmp%S(:,j,k), p_ref, rho_in_situ, &
                             tv%eqn_of_state, EOSdom)
      do i=Isq,Ieq+1 ; alpha_star(i,j,k) = 1.0 / rho_in_situ(i) ; enddo
    enddo ; enddo
  endif                                               ! use_EOS

  if (use_EOS) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        M(i,j,nz) = geopot_bot(i,j) + p(i,j,nz+1) * alpha_star(i,j,nz)
      enddo
      do k=nz-1,1,-1 ; do i=Isq,Ieq+1
        M(i,j,k) = M(i,j,k+1) + p(i,j,K+1) * (alpha_star(i,j,k) - alpha_star(i,j,k+1))
      enddo ; enddo
    enddo
  else ! not use_EOS
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        M(i,j,nz) = geopot_bot(i,j) + p(i,j,nz+1) * alpha_Lay(nz)
      enddo
      do k=nz-1,1,-1 ; do i=Isq,Ieq+1
        M(i,j,k) = M(i,j,k+1) + p(i,j,K+1) * dalpha_int(K+1)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (CS%GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      dM(i,j) = (CS%GFS_scale - 1.0) * M(i,j,1)
    enddo ; enddo
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,k) = M(i,j,k) + dM(i,j)
    enddo ; enddo ; enddo

    !   Could instead do the following, to avoid taking small differences
    ! of large numbers...
!   do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
!     M(i,j,1) = CS%GFS_scale * M(i,j,1)
!   enddo ; enddo
!   if (use_EOS) then
!     do k=2,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
!       M(i,j,k) = M(i,j,k-1) - p(i,j,K) * (alpha_star(i,j,k-1) - alpha_star(i,j,k))
!     enddo ; enddo ; enddo
!   else ! not use_EOS
!     do k=2,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
!        M(i,j,k) = M(i,j,k-1) - p(i,j,K) * dalpha_int(K)
!     enddo ; enddo ; enddo
!   endif ! use_EOS

  endif

  ! Note that ddM/dPb = alpha_star(i,j,1)
  if (pbce_a%associated()) then
    call Set_pbce_nonBouss(bxH0, bxQ0, bxHkint, p_a, tv_tmp, h_a%lb(1), g_Earth, H_to_RZ, &
                           H_subroundoff, Rlay_a, CS%GFS_scale, pbce_a, alpha_star_a)
  endif

!    Calculate the pressure force. On a Cartesian grid,
!      PFu = - dM/dx   and  PFv = - dM/dy.
  if (use_EOS) then
    !$OMP parallel do default(shared) private(dp_star,PFu_bc,PFv_bc)
    do k=1,nz
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dp_star(i,j) = (p(i,j,K+1) - p(i,j,K)) + dp_neglect
      enddo ; enddo
      do j=js,je ; do I=Isq,Ieq
        ! PFu_bc = p* grad alpha*
        PFu_bc = (alpha_star(i+1,j,k) - alpha_star(i,j,k)) * (IdxCu(I,j) * &
            ((dp_star(i,j)*dp_star(i+1,j) + ((p(i,j,K)*dp_star(i+1,j)) + (p(i+1,j,K)*dp_star(i,j)))) / &
             (dp_star(i,j) + dp_star(i+1,j))))
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * IdxCu(I,j) + PFu_bc
        if (allocated(CS%PFu_bc)) CS%PFu_bc(i,j,k) = PFu_bc
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv_bc = (alpha_star(i,j+1,k) - alpha_star(i,j,k)) * (IdyCv(i,J) * &
            ((dp_star(i,j)*dp_star(i,j+1) + ((p(i,j,K)*dp_star(i,j+1)) + (p(i,j+1,K)*dp_star(i,j)))) / &
             (dp_star(i,j) + dp_star(i,j+1))))
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * IdyCv(i,J) + PFv_bc
        if (allocated(CS%PFv_bc)) CS%PFv_bc(i,j,k) = PFv_bc
      enddo ; enddo
    enddo ! k-loop
  else ! .not. use_EOS
    !$OMP parallel do default(shared)
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * IdxCu(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * IdyCv(i,J)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (CS%id_PFu_bc>0) call post_data(CS%id_PFu_bc, CS%PFu_bc, CS%diag)
  if (CS%id_PFv_bc>0) call post_data(CS%id_PFv_bc, CS%PFv_bc, CS%diag)
  ! To be consistent with old runs, tidal forcing diagnostic also includes total SAL.
  ! New diagnostics are given for each individual field.
  if (CS%id_e_tide>0) call post_data(CS%id_e_tide, e_sal+e_tide_eq+e_tide_sal, CS%diag)
  if (CS%id_e_sal>0) call post_data(CS%id_e_sal, e_sal, CS%diag)
  if (CS%id_e_tide_eq>0) call post_data(CS%id_e_tide_eq, e_tide_eq, CS%diag)
  if (CS%id_e_tide_sal>0) call post_data(CS%id_e_tide_sal, e_tide_sal, CS%diag)

  call p_a%free() ; call alpha_star_a%free()

end subroutine PressureForce_Mont_nonBouss

!> \brief Boussinesq Montgomery-potential form of pressure gradient
!!
!! Determines the acceleration due to pressure forces.
!!
!! To work, the following fields must be set outside of the usual (is:ie,js:je)
!! range before this subroutine is called:
!!   h(isB:ie+1,jsB:je+1), T(isB:ie+1,jsB:je+1), and S(isB:ie+1,jsB:je+1).
subroutine PressureForce_Mont_Bouss(bxH0, bxQ0, bxHpgf, h_a, tv, PFu_a, PFv_a, G, Gcore, US, &
                                    CS, p_atm_a, pbce_a, eta_a, Rlay_a, g_prime_a, g_Earth, &
                                    H_to_Z, Z_to_H, Rho0, dZ_subroundoff, nk_rho_varies)
  type(Box_t),          intent(in)  :: bxH0  !< The h-point box of the whole column,
                                             !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),          intent(in)  :: bxQ0  !< The B-grid-index box of the whole column,
                                             !! [IscB:IecB, JscB:JecB, 1:nz]
  type(Box_t),          intent(in)  :: bxHpgf !< The h points the u and v pressure gradients
                                              !! read, [IscB:IecB+1, JscB:JecB+1, 1:nz]
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(grid_core_type), intent(in)  :: Gcore !< Persistent copies of the core ocean grid fields
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  type(RealArray_t), intent(in) :: h_a !< Layer thickness [H ~> m].
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables.
  type(RealArray_t), intent(inout) :: PFu_a !< Zonal acceleration due to pressure gradients
                                                                 !! (equal to -dM/dx) [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: PFv_a !< Meridional acceleration due to pressure gradients
                                                                 !! (equal to -dM/dy) [L T-2 ~> m s-2].
  type(PressureForce_Mont_CS),                intent(inout) :: CS  !< Control structure for Montgomery potential PGF
  type(RealArray_t), intent(in) :: p_atm_a !< The pressure at the ice-ocean or
                                                                !! atmosphere-ocean [R L2 T-2 ~> Pa].
  type(RealArray_t), intent(inout) :: pbce_a !< The baroclinic pressure anomaly in
                                                                !! each layer due to free surface height anomalies
                                                                !! [L2 T-2 H-1 ~> m s-2].
  type(RealArray_t), intent(inout) :: eta_a !< Free surface height [H ~> m].
  type(RealArray_t), intent(in) :: Rlay_a !< The target coordinate value (potential density) in each
                                          !! layer [R ~> kg m-3].
  type(RealArray_t), intent(in) :: g_prime_a !< The reduced gravity at each interface
                                             !! [L2 Z-1 T-2 ~> m s-2].
  real,              intent(in) :: g_Earth !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real,              intent(in) :: H_to_Z !< A constant that translates thickness units to the units
                                          !! of depth [Z H-1 ~> 1 or m3 kg-1].
  real,              intent(in) :: Z_to_H !< A constant that translates depth units to thickness
                                          !! units depth [H Z-1 ~> 1 or kg m-3].
  real,              intent(in) :: Rho0 !< The density used in the Boussinesq approximation or
                                        !! nominal density used to convert depths into mass units
                                        !! [R ~> kg m-3].
  real,              intent(in) :: dZ_subroundoff !< A thickness in height units that is so small
                                                  !! that it can be added to a vertical distance of
                                                  !! Angstrom_Z or 1e-17 m without changing it at
                                                  !! the bit level [Z ~> m].  This is the height
                                                  !! equivalent of H_subroundoff.
  integer,           intent(in) :: nk_rho_varies !< The number of layers at the top where the
                                                 !! density does not track any target density.
  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: h, PFu, PFv, pbce
  real, dimension(:,:), contiguous, pointer :: p_atm, eta
  real, dimension(:,:), contiguous, pointer :: bathyT, meanSL, IdxCu, IdyCv
  real, dimension(:), contiguous, pointer :: Rlay, g_prime
  type(RealArray_t) :: e_a, rho_star_a
  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2),h_a%lb(3):h_a%ub(3)) :: &
    M           ! The Montgomery potential, M = (p/rho + gz) [L2 T-2 ~> m2 s-2].
  real, dimension(:,:,:), contiguous, pointer :: &
    rho_star    ! In-situ density divided by the derivative with depth of the
                ! corrected e times (G_Earth/Rho0) [L2 Z-1 T-2 ~> m s-2].
  real, dimension(:,:,:), contiguous, pointer :: e ! Interface height [Z ~> m].
                ! e may be adjusted (with a nonlinear equation of state) so that
                ! its derivative compensates for the adiabatic compressibility
                ! in seawater, but e will still be close to the interface depth.
  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2),h_a%lb(3):h_a%ub(3)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [C ~> degC].
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [S ~> ppt].

  real :: Rho_cv_BL(h_a%lb(1):h_a%ub(1)) !   The coordinate potential density in
                ! the deepest variable density near-surface layer [R ~> kg m-3].
  real :: h_star(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2)) ! Layer thickness after compensation
                             ! for compressibility [Z ~> m].
  real :: SSH(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2)) ! The sea surface height anomaly, in depth units [Z ~> m].
  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2)) :: &
    e_sal ! The bottom geopotential anomaly due to self-attraction and loading [Z ~> m].
  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2)) :: &
    e_tide_eq ! Bottom geopotential anomaly due to tidal forces from astronomical sources
                                     ! [Z ~> m].
  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2)) :: &
    e_tide_sal ! Bottom geopotential anomaly due to harmonic self-attraction and loading
                                      ! specific to tides, in depth units [Z ~> m].
  real :: p_ref(h_a%lb(1):h_a%ub(1))     !   The pressure used to calculate the coordinate
                             ! density [R L2 T-2 ~> Pa] (usually 2e7 Pa = 2000 dbar).
  real :: I_Rho0             ! 1/Rho0 [R-1 ~> m3 kg-1].
  real :: G_Rho0             ! G_Earth / Rho0 [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1].
  real :: PFu_bc, PFv_bc     ! The pressure gradient force due to along-layer
                             ! compensated density gradients [L T-2 ~> m s-2]
  real :: dz_neglect         ! A vertical distance that is so small it is usually lost
                             ! in roundoff and can be neglected [Z ~> m].
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_EOS         ! If true, density is calculated from T & S using
                             ! an equation of state.
  logical :: is_split        ! A flag indicating whether the pressure
                             ! gradient terms are to be split into
                             ! barotropic and baroclinic pieces.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.
  integer, dimension(2) :: EOSdom ! The computational domain for the equation of state
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: i, j, k

  is = bxH0%idxS(1) ; ie = bxH0%idxE(1) ; js = bxH0%idxS(2) ; je = bxH0%idxE(2) ; nz = bxH0%idxE(3)
  nkmb=nk_rho_varies
  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2)
  EOSdom(1) = Isq - (h_a%lb(1)-1) ;  EOSdom(2) = bxH0%idxE(1)+1 - (h_a%lb(1)-1)

  call h_a%view(h) ; call PFu_a%view(PFu) ; call PFv_a%view(PFv)
  call Rlay_a%view(Rlay) ; call g_prime_a%view(g_prime)
  call Gcore%bathyT%view(bathyT) ; call Gcore%meanSL%view(meanSL)
  call Gcore%IdxCu%view(IdxCu) ; call Gcore%IdyCv%view(IdyCv)
  call e_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),bxH0%idxE(3)+1]) ; call e_a%view(e)
  call rho_star_a%alloc(lb=[h_a%lb(1),h_a%lb(2),1], ub=[h_a%ub(1),h_a%ub(2),bxH0%idxE(3)])
  call rho_star_a%view(rho_star)
  nullify(p_atm, pbce, eta)
  if (p_atm_a%associated()) call p_atm_a%view(p_atm)
  if (pbce_a%associated()) call pbce_a%view(pbce)
  if (eta_a%associated()) call eta_a%view(eta)

  use_p_atm = associated(p_atm)
  is_split = pbce_a%associated()
  use_EOS = associated(tv%eqn_of_state)

  if (.not.CS%initialized) call MOM_error(FATAL, &
       "MOM_PressureForce_Mont: Module must be initialized before it is used.")

  if (use_EOS) then
    if (query_compressible(tv%eqn_of_state)) call MOM_error(FATAL, &
      "PressureForce_Mont_Bouss: The Montgomery form of the pressure force "//&
      "can no longer be used with a compressible EOS. Use #define ANALYTIC_FV_PGF.")
  endif

  dz_neglect = dZ_subroundoff
  I_Rho0 = 1.0/CS%Rho0
  G_Rho0 = g_Earth / Rho0

  !$OMP parallel do default(shared)
  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    e(i,j,nz+1) = -bathyT(i,j)
  enddo ; enddo

  ! Calculate and add the self-attraction and loading geopotential anomaly.
  if (CS%calculate_SAL) then
    !   Determine the surface height anomaly for calculating self attraction
    ! and loading.  This should really be based on bottom pressure anomalies,
    ! but that is not yet implemented, and the current form is correct for
    ! barotropic tides.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1 ; SSH(i,j) = min(-bathyT(i,j) - meanSL(i,j), 0.0) ; enddo
      do k=1,nz ; do i=Isq,Ieq+1
        SSH(i,j) = SSH(i,j) + h(i,j,k)*H_to_Z
      enddo ; enddo
    enddo
    call calc_SAL(SSH, e_sal, G, CS%SAL_CSp, tmp_scale=US%Z_to_m)
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e(i,j,nz+1) = e(i,j,nz+1) - e_sal(i,j)
    enddo ; enddo
  endif

  ! Calculate and add the tidal geopotential anomaly.
  if (CS%tides) then
    call calc_tidal_forcing(CS%Time, e_tide_eq, e_tide_sal, G, US, CS%tides_CSp)
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e(i,j,nz+1) = e(i,j,nz+1) - (e_tide_eq(i,j) + e_tide_sal(i,j))
    enddo ; enddo
  endif

  !$OMP parallel do default(shared)
  do j=Jsq,Jeq+1 ; do k=nz,1,-1 ; do i=Isq,Ieq+1
    e(i,j,K) = e(i,j,K+1) + h(i,j,k)*H_to_Z
  enddo ; enddo ; enddo

  if (use_EOS) then
!   Calculate in-situ densities (rho_star).

! With a bulk mixed layer, replace the T & S of any layers that are
! lighter than the buffer layer with the properties of the buffer
! layer.  These layers will be massless anyway, and it avoids any
! formal calculations with hydrostatically unstable profiles.

    if (nkmb>0) then
      tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
      tv_tmp%eqn_of_state => tv%eqn_of_state

      do concurrent (i=bxHpgf%idxS(1):bxHpgf%idxE(1)) ; p_ref(i) = tv%P_Ref ; enddo
      !$OMP parallel do default(shared) private(Rho_cv_BL)
      do j=Jsq,Jeq+1
        do k=1,nkmb ; do i=Isq,Ieq+1
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        enddo ; enddo
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, Rho_cv_BL(:), &
                               tv%eqn_of_state, EOSdom)

        do k=nkmb+1,nz ; do i=Isq,Ieq+1
          if (Rlay(k) < Rho_cv_BL(i)) then
            tv_tmp%T(i,j,k) = tv%T(i,j,nkmb) ; tv_tmp%S(i,j,k) = tv%S(i,j,nkmb)
          else
            tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
          endif
        enddo ; enddo
      enddo
    else
      tv_tmp%T => tv%T ; tv_tmp%S => tv%S
      tv_tmp%eqn_of_state => tv%eqn_of_state
      do concurrent (i=bxHpgf%idxS(1):bxHpgf%idxE(1)) ; p_ref(i) = 0.0 ; enddo
    endif

    ! This no longer includes any pressure dependency, since this routine
    ! will come down with a fatal error if there is any compressibility.
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=Jsq,Jeq+1
      call calculate_density(tv_tmp%T(:,j,k), tv_tmp%S(:,j,k), p_ref, rho_star(:,j,k), &
                             tv%eqn_of_state, EOSdom)
      do i=Isq,Ieq+1 ; rho_star(i,j,k) = G_Rho0*rho_star(i,j,k) ; enddo
    enddo ; enddo
  endif                                               ! use_EOS

!    Here the layer Montgomery potentials, M, are calculated.
  if (use_EOS) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        M(i,j,1) = CS%GFS_scale * (rho_star(i,j,1) * e(i,j,1))
        if (use_p_atm) M(i,j,1) = M(i,j,1) + p_atm(i,j) * I_Rho0
      enddo
      do k=2,nz ; do i=Isq,Ieq+1
        M(i,j,k) = M(i,j,k-1) + (rho_star(i,j,k) - rho_star(i,j,k-1)) * e(i,j,K)
      enddo ; enddo
    enddo
  else ! not use_EOS
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        M(i,j,1) = g_prime(1) * e(i,j,1)
        if (use_p_atm) M(i,j,1) = M(i,j,1) + p_atm(i,j) * I_Rho0
      enddo
      do k=2,nz ; do i=Isq,Ieq+1
        M(i,j,k) = M(i,j,k-1) + g_prime(K) * e(i,j,K)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (pbce_a%associated()) then
    call Set_pbce_Bouss(bxH0, bxQ0, bxHpgf, e_a, tv_tmp, Gcore, h_a%lb(1), h_a%lb(2), H_to_Z, &
                        g_Earth, Rho0, dZ_subroundoff, g_prime_a, CS%Rho0, &
                        CS%GFS_scale, pbce_a, rho_star_a, 1)
  endif

!    Calculate the pressure force. On a Cartesian grid,
!      PFu = - dM/dx   and  PFv = - dM/dy.
  if (use_EOS) then
    !$OMP parallel do default(shared) private(h_star,PFu_bc,PFv_bc)
    do k=1,nz
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        h_star(i,j) = (e(i,j,K) - e(i,j,K+1)) + dz_neglect
      enddo ; enddo
      do j=js,je ; do I=Isq,Ieq
        PFu_bc = -1.0*(rho_star(i+1,j,k) - rho_star(i,j,k)) * (IdxCu(I,j) * &
          ((h_star(i,j) * h_star(i+1,j) - ((e(i,j,K) * h_star(i+1,j)) + &
          (e(i+1,j,K) * h_star(i,j)))) / (h_star(i,j) + h_star(i+1,j))))
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * IdxCu(I,j) + PFu_bc
        if (allocated(CS%PFu_bc)) CS%PFu_bc(i,j,k) = PFu_bc
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv_bc = -1.0*(rho_star(i,j+1,k) - rho_star(i,j,k)) * (IdyCv(i,J) * &
          ((h_star(i,j) * h_star(i,j+1) - ((e(i,j,K) * h_star(i,j+1)) + &
          (e(i,j+1,K) * h_star(i,j)))) / (h_star(i,j) + h_star(i,j+1))))
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * IdyCv(i,J) + PFv_bc
        if (allocated(CS%PFv_bc)) CS%PFv_bc(i,j,k) = PFv_bc
      enddo ; enddo
    enddo ! k-loop
  else ! .not. use_EOS
    !$OMP parallel do default(shared)
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * IdxCu(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * IdyCv(i,J)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (eta_a%associated()) then
    ! eta is the sea surface height relative to a time-invariant geoid, for
    ! comparison with what is used for eta in btstep.  See how e was calculated
    ! about 200 lines above.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta(i,j) = e(i,j,1)*Z_to_H
    enddo ; enddo
    if (CS%tides) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = eta(i,j) + (e_tide_eq(i,j)+e_tide_sal(i,j))*Z_to_H
      enddo ; enddo
    endif
    if (CS%calculate_SAL) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = eta(i,j) + e_sal(i,j)*Z_to_H
      enddo ; enddo
    endif
  endif

  if (CS%id_PFu_bc>0) call post_data(CS%id_PFu_bc, CS%PFu_bc, CS%diag)
  if (CS%id_PFv_bc>0) call post_data(CS%id_PFv_bc, CS%PFv_bc, CS%diag)
  ! To be consistent with old runs, tidal forcing diagnostic also includes total SAL.
  ! New diagnostics are given for each individual field.
  if (CS%id_e_tide>0) call post_data(CS%id_e_tide, e_sal+e_tide_eq+e_tide_sal, CS%diag)
  if (CS%id_e_sal>0) call post_data(CS%id_e_sal, e_sal, CS%diag)
  if (CS%id_e_tide_eq>0) call post_data(CS%id_e_tide_eq, e_tide_eq, CS%diag)
  if (CS%id_e_tide_sal>0) call post_data(CS%id_e_tide_sal, e_tide_sal, CS%diag)

  call e_a%free() ; call rho_star_a%free()

end subroutine PressureForce_Mont_Bouss

!> Determines the partial derivative of the acceleration due
!! to pressure forces with the free surface height.
subroutine Set_pbce_Bouss(bxH0, bxQ0, bxHpgf, e_a, tv, Gcore, isd, jsd, H_to_Z, g_Earth, GV_Rho0, &
                          dZ_subroundoff, g_prime_a, Rho0, GFS_scale, pbce_a, rho_star_a, nkblock)
  type(Box_t),          intent(in)  :: bxH0 !< The h-point box of the whole column,
                                            !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),          intent(in)  :: bxQ0 !< The B-grid-index box of the whole column,
                                            !! [IscB:IecB, JscB:JecB, 1:nz]
  type(Box_t),          intent(in)  :: bxHpgf !< The h points the u and v pressure gradients
                                              !! read, [IscB:IecB+1, JscB:JecB+1, 1:nz]
  type(RealArray_t), intent(in) :: e_a !< Interface height [Z ~> m].
  type(thermo_var_ptrs),                intent(in)  :: tv   !< Thermodynamic variables
  type(grid_core_type), intent(in)  :: Gcore !< Persistent copies of the core ocean grid fields
  integer, intent(in) :: isd !< The start i-index of cell centers within the data domain
  integer, intent(in) :: jsd !< The start j-index of cell centers within the data domain
  real,    intent(in) :: H_to_Z !< A constant that translates thickness units to the units of
                                !! depth [Z H-1 ~> 1 or m3 kg-1].
  real,    intent(in) :: g_Earth !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real,    intent(in) :: GV_Rho0 !< The density used in the Boussinesq approximation or nominal
                                 !! density used to convert depths into mass units [R ~> kg m-3].
  real,    intent(in) :: dZ_subroundoff !< A thickness in height units that is so small that it
                     !! can be added to a vertical distance of Angstrom_Z or 1e-17 m without
                     !! changing it at the bit level [Z ~> m].  This is the height equivalent
                     !! of H_subroundoff.
  type(RealArray_t), intent(in) :: g_prime_a !< The reduced gravity at each interface
                                             !! [L2 Z-1 T-2 ~> m s-2].
  real,                                 intent(in)  :: Rho0 !< The "Boussinesq" ocean density [R ~> kg m-3].
  real,                                 intent(in)  :: GFS_scale !< Ratio between gravity applied to top
                                                            !! interface and the gravitational acceleration of
                                                            !! the planet [nondim]. Usually this ratio is 1.
  type(RealArray_t), intent(inout) :: pbce_a !< The baroclinic pressure anomaly in each layer due
                                                            !! to free surface height anomalies
                                                            !! [L2 T-2 H-1 ~> m s-2].
  type(RealArray_t), &
                    intent(in) :: rho_star_a !< The layer densities (maybe compressibility
                                                            !! compensated), times g/rho_0 [L2 Z-1 T-2 ~> m s-2].
  integer, intent(in) :: nkblock !< The number of vertical levels in each block for the variable density
                                 !! near-surface layer.  This is needed to determine whether rho_star is
                                 !! compressibility compensated or not.
  real, dimension(:,:,:), contiguous, pointer :: e, pbce, rho_star
  real, dimension(:,:), contiguous, pointer :: meanSL
  real, dimension(:), contiguous, pointer :: g_prime
  real :: Ihtot(e_a%lb(1):e_a%ub(1),e_a%lb(2):e_a%ub(2))  ! The inverse of the sum of the layer
                                  ! thicknesses [H-1 ~> m-1 or m2 kg-1].
  real :: press(e_a%lb(1):e_a%ub(1),e_a%lb(2):e_a%ub(2),nkblock) ! Interface pressure
                                  ! [R L2 T-2 ~> Pa].
  real :: T_int(e_a%lb(1):e_a%ub(1),e_a%lb(2):e_a%ub(2),nkblock) ! Interface temperature [C ~> degC]
  real :: S_int(e_a%lb(1):e_a%ub(1),e_a%lb(2):e_a%ub(2),nkblock)  ! Interface salinity [S ~> ppt]
  real :: dR_dT(e_a%lb(1):e_a%ub(1),e_a%lb(2):e_a%ub(2),nkblock)  ! Partial derivative of density
                                  ! with temperature [R C-1 ~> kg m-3 degC-1]
  real :: dR_dS(e_a%lb(1):e_a%ub(1),e_a%lb(2):e_a%ub(2),nkblock)  ! Partial derivative of density
                                  ! with salinity [R S-1 ~> kg m-3 ppt-1].
  real :: rho_in_situ(e_a%lb(1):e_a%ub(1),e_a%lb(2):e_a%ub(2))  ! In-situ density at the top of a
                                  ! layer [R ~> kg m-3].
  real :: G_Rho0             ! A scaled version of g_Earth / Rho0 [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]
  real :: Rho0xG             ! g_Earth * Rho0 [R L2 Z-1 T-2 ~> kg s-2 m-2]
  logical :: use_EOS         ! If true, density is calculated from T & S using
                             ! an equation of state.
  real :: dz_neglect         ! A vertical distance that is so small it is usually lost
                             ! in roundoff and can be neglected [Z ~> m].
  integer :: EOSdom(3,2)     ! The computational domain for the equation of state
  integer :: Isq, Ieq, Jsq, Jeq, nz, i, j, k, kstart, kend, kk
  integer :: iec, jec

  call e_a%view(e) ; call pbce_a%view(pbce)
  call Gcore%meanSL%view(meanSL) ; call g_prime_a%view(g_prime)
  nullify(rho_star)
  if (rho_star_a%associated()) call rho_star_a%view(rho_star)

  !$omp target data map(alloc: Ihtot)
  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2) ; nz = bxH0%idxE(3)
  iec = bxH0%idxE(1) ; jec = bxH0%idxE(2)

  use_EOS = associated(tv%eqn_of_state)
  dz_neglect = dZ_subroundoff

  if (use_EOS) then
    if (rho_star_a%associated()) then
      do concurrent (j=bxHpgf%idxS(2):bxHpgf%idxE(2))
        do concurrent (i=bxHpgf%idxS(1):bxHpgf%idxE(1))
          Ihtot(i,j) = H_to_Z / ((e(i,j,1) - e(i,j,nz+1)) + dz_neglect)
          pbce(i,j,1) = GFS_scale * rho_star(i,j,1) * H_to_Z
        enddo

        do k=2,nz
          do concurrent (i=bxHpgf%idxS(1):bxHpgf%idxE(1))
            pbce(i,j,k) = pbce(i,j,k-1) + (rho_star(i,j,k) - rho_star(i,j,k-1)) &
                * ((e(i,j,K) - e(i,j,nz+1)) * Ihtot(i,j))
          enddo
        enddo
      enddo
    else
      !$omp target data &
      !$omp   map(alloc: EOSdom, press, T_int, S_int, rho_in_situ) &
      !$omp   map(alloc: dR_dT, dR_dS)

      Rho0xG = Rho0 * g_Earth
      G_Rho0 = g_Earth / GV_Rho0

      EOSdom(1,:) = [Isq - (isd-1), iec+1 - (isd-1)]
      EOSdom(2,:) = [Jsq - (jsd-1), jec+1 - (jsd-1)]
      do concurrent (j=bxHpgf%idxS(2):bxHpgf%idxE(2), i=bxHpgf%idxS(1):bxHpgf%idxE(1))
        Ihtot(i,j) = H_to_Z / ((e(i,j,1) - e(i,j,nz+1)) + dz_neglect)
        press(i,j,1) = -Rho0xG * (e(i,j,1) - meanSL(i,j))
      enddo

      call calculate_density(tv%T(:,:,1), tv%S(:,:,1), press(:,:,1), rho_in_situ, &
                             tv%eqn_of_state, EOSdom(1:2,:))

      do concurrent (j=bxHpgf%idxS(2):bxHpgf%idxE(2), i=bxHpgf%idxS(1):bxHpgf%idxE(1))
        pbce(i,j,1) = G_Rho0 * (GFS_scale * rho_in_situ(i,j)) * H_to_Z
      enddo

      EOSdom(3,1) = 1
      do kstart=2,nz,nkblock
        kend = min(kstart+nkblock-1, nz)
        do concurrent (k=kstart:kend, j=Jsq:Jeq+1, i=Isq:Ieq+1)
          kk=k-kstart+1
          press(i,j,kk) = -Rho0xG * (e(i,j,K) - meanSL(i,j))
          T_int(i,j,kk) = 0.5 * (tv%T(i,j,k-1) + tv%T(i,j,k))
          S_int(i,j,kk) = 0.5 * (tv%S(i,j,k-1) + tv%S(i,j,k))
        enddo
        EOSdom(3,2) = kend-kstart+1
        call calculate_density_derivs(T_int, S_int, press, dR_dT, dR_dS, &
                                      tv%eqn_of_state, EOSdom)

        do concurrent(j=bxHpgf%idxS(2):bxHpgf%idxE(2))
          do k=kstart,kend ; do concurrent (i=bxHpgf%idxS(1):bxHpgf%idxE(1))
          kk=k-kstart+1
          pbce(i,j,k) = pbce(i,j,k-1) + G_Rho0 * &
             ((e(i,j,K) - e(i,j,nz+1)) * Ihtot(i,j)) * &
             (dR_dT(i,j,kk) * (tv%T(i,j,k) - tv%T(i,j,k-1)) + &
              dR_dS(i,j,kk) * (tv%S(i,j,k) - tv%S(i,j,k-1)))
          enddo ; enddo
        enddo
      enddo
      !$omp end target data
    endif
  else
    do concurrent (j=bxHpgf%idxS(2):bxHpgf%idxE(2))
      do concurrent (i=bxHpgf%idxS(1):bxHpgf%idxE(1))
        Ihtot(i,j) = 1.0 / ((e(i,j,1) - e(i,j,nz+1)) + dz_neglect)
        pbce(i,j,1) = g_prime(1) * H_to_Z
      enddo
      do k=2,nz
        do concurrent (i=bxHpgf%idxS(1):bxHpgf%idxE(1))
          pbce(i,j,k) = pbce(i,j,k-1) + (g_prime(K) * H_to_Z) &
              * ((e(i,j,K) - e(i,j,nz+1)) * Ihtot(i,j))
        enddo
      enddo
    enddo
  endif
  !$omp end target data
end subroutine Set_pbce_Bouss

!> Determines the partial derivative of the acceleration due
!! to pressure forces with the column mass.
subroutine Set_pbce_nonBouss(bxH0, bxQ0, bxHkint, p_a, tv, isd, g_Earth, H_to_RZ, H_subroundoff, &
                             Rlay_a, GFS_scale, pbce_a, alpha_star_a)
  type(Box_t),          intent(in)  :: bxH0 !< The h-point box of the whole column,
                                            !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),          intent(in)  :: bxQ0 !< The B-grid-index box of the whole column,
                                            !! [IscB:IecB, JscB:JecB, 1:nz]
  type(Box_t),          intent(in)  :: bxHkint !< The interior interfaces of the whole column,
                                               !! [isc:iec, jsc:jec, 2:nz]
  type(RealArray_t), intent(in) :: p_a !< Interface pressures [R L2 T-2 ~> Pa].
  type(thermo_var_ptrs),                intent(in)  :: tv !< Thermodynamic variables
  integer, intent(in) :: isd !< The start i-index of cell centers within the data domain
  real,    intent(in) :: g_Earth !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real,    intent(in) :: H_to_RZ !< A constant that translates thickness units to the units of
                                 !! mass per unit area [R Z H-1 ~> kg m-3 or 1].
  real,    intent(in) :: H_subroundoff !< A thickness that is so small that it can be added to
                     !! a thickness of Angstrom or larger without changing it at the bit level
                     !! [H ~> m or kg m-2]. If Angstrom is 0 or exceedingly small, this is
                     !! negligible compared to 1e-17 m.
  type(RealArray_t), intent(in) :: Rlay_a !< The target coordinate value (potential density)
                                          !! in each layer [R ~> kg m-3].
  real,                                 intent(in)  :: GFS_scale !< Ratio between gravity applied to top
                                                          !! interface and the gravitational acceleration of
                                                          !! the planet [nondim]. Usually this ratio is 1.
  type(RealArray_t), intent(inout) :: pbce_a !< The baroclinic pressure anomaly in each
                                                          !! layer due to free surface height anomalies
                                                          !! [L2 H-1 T-2 ~> m4 kg-1 s-2].
  type(RealArray_t), intent(in) :: alpha_star_a !< The layer specific volumes
                                                          !! (maybe compressibility compensated) [R-1 ~> m3 kg-1].
  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: p, pbce, alpha_star
  real, dimension(:), contiguous, pointer :: Rlay
  real, dimension(p_a%lb(1):p_a%ub(1),p_a%lb(2):p_a%ub(2)) :: &
    dpbce, &      !   A barotropic correction to the pbce to enable the use of
                  ! a reduced gravity form of the equations [L2 H-1 T-2 ~> m4 kg-1 s-2].
    C_htot        ! dP_dH divided by the total ocean pressure [H-1 ~> m2 kg-1].
  real :: T_int(p_a%lb(1):p_a%ub(1))     ! Interface temperature [C ~> degC]
  real :: S_int(p_a%lb(1):p_a%ub(1))     ! Interface salinity [S ~> ppt]
  real :: dR_dT(p_a%lb(1):p_a%ub(1))     ! Partial derivative of density with temperature
                             ! [R C-1 ~> kg m-3 degC-1]
  real :: dR_dS(p_a%lb(1):p_a%ub(1))     ! Partial derivative of density with salinity
                             ! [R S-1 ~> kg m-3 ppt-1].
  real :: rho_in_situ(p_a%lb(1):p_a%ub(1)) ! In-situ density at an interface [R ~> kg m-3].
  real :: alpha_Lay(pbce_a%lb(3):pbce_a%ub(3)) ! The specific volume of each layer [R-1 ~> m3 kg-1].
  real :: dalpha_int(p_a%lb(3):p_a%ub(3)) ! The change in specific volume across each interface
                             ! [R-1 ~> m3 kg-1].
  real :: dP_dH              ! A factor that converts from thickness to pressure times other dimensional
                             ! conversion factors [R L2 T-2 H-1 ~> Pa m2 kg-1].
  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [R L2 T-2 ~> Pa].
  logical :: use_EOS         ! If true, density is calculated from T & S using an equation of state.
  integer, dimension(2) :: EOSdom ! The computational domain for the equation of state
  integer :: Isq, Ieq, Jsq, Jeq, nz, i, j, k
  integer :: iec

  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2) ; nz = bxH0%idxE(3)
  iec = bxH0%idxE(1)
  EOSdom(1) = Isq - (isd-1) ;  EOSdom(2) = iec+1 - (isd-1)

  call p_a%view(p) ; call pbce_a%view(pbce) ; call Rlay_a%view(Rlay)
  nullify(alpha_star)
  if (alpha_star_a%associated()) call alpha_star_a%view(alpha_star)

  use_EOS = associated(tv%eqn_of_state)

  dP_dH = g_Earth * H_to_RZ
  dp_neglect = g_Earth * H_to_RZ * H_subroundoff

  if (use_EOS) then
    if (alpha_star_a%associated()) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1
        do i=Isq,Ieq+1
          C_htot(i,j) = dP_dH / ((p(i,j,nz+1)-p(i,j,1)) + dp_neglect)
          pbce(i,j,nz) = dP_dH * alpha_star(i,j,nz)
        enddo
        do k=nz-1,1,-1 ; do i=Isq,Ieq+1
          pbce(i,j,k) = pbce(i,j,k+1) + ((p(i,j,K+1)-p(i,j,1)) * C_htot(i,j)) * &
              (alpha_star(i,j,k) - alpha_star(i,j,k+1))
        enddo ; enddo
      enddo
    else
      !$OMP parallel do default(shared) private(T_int,S_int,dR_dT,dR_dS,rho_in_situ)
      do j=Jsq,Jeq+1
        call calculate_density(tv%T(:,j,nz), tv%S(:,j,nz), p(:,j,nz+1), rho_in_situ, &
                               tv%eqn_of_state, EOSdom)
        do i=Isq,Ieq+1
          C_htot(i,j) = dP_dH / ((p(i,j,nz+1)-p(i,j,1)) + dp_neglect)
          pbce(i,j,nz) = dP_dH / (rho_in_situ(i))
        enddo
        do k=nz-1,1,-1
          do i=Isq,Ieq+1
            T_int(i) = 0.5*(tv%T(i,j,k)+tv%T(i,j,k+1))
            S_int(i) = 0.5*(tv%S(i,j,k)+tv%S(i,j,k+1))
          enddo
          call calculate_density(T_int, S_int, p(:,j,k+1), rho_in_situ, tv%eqn_of_state, EOSdom)
          call calculate_density_derivs(T_int, S_int, p(:,j,k+1), dR_dT, dR_dS, &
                                        tv%eqn_of_state, EOSdom)
          do i=Isq,Ieq+1
            pbce(i,j,k) = pbce(i,j,k+1) + ((p(i,j,K+1)-p(i,j,1))*C_htot(i,j)) *  &
                ((dR_dT(i)*(tv%T(i,j,k+1)-tv%T(i,j,k)) + &
                  dR_dS(i)*(tv%S(i,j,k+1)-tv%S(i,j,k))) / (rho_in_situ(i)**2))
          enddo
        enddo
      enddo
    endif
  else ! not use_EOS

    do concurrent (k=bxH0%idxS(3):bxH0%idxE(3)) ; alpha_Lay(k) = 1.0 / (Rlay(k)) ; enddo
    do concurrent (k=bxHkint%idxS(3):bxHkint%idxE(3)) ; dalpha_int(K) = alpha_Lay(k-1) - alpha_Lay(k) ; enddo

    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        C_htot(i,j) = dP_dH / ((p(i,j,nz+1)-p(i,j,1)) + dp_neglect)
        pbce(i,j,nz) = dP_dH * alpha_Lay(nz)
      enddo
      do k=nz-1,1,-1 ; do i=Isq,Ieq+1
        pbce(i,j,k) = pbce(i,j,k+1) + ((p(i,j,K+1)-p(i,j,1))*C_htot(i,j)) * dalpha_int(K+1)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      dpbce(i,j) = (GFS_scale - 1.0) * pbce(i,j,1)
    enddo ; enddo
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pbce(i,j,k) = pbce(i,j,k) + dpbce(i,j)
    enddo ; enddo ; enddo
  endif

end subroutine Set_pbce_nonBouss

!> Initialize the Montgomery-potential form of PGF control structure
subroutine PressureForce_Mont_init(Time, G, GV, US, param_file, diag, CS, SAL_CSp, tides_CSp)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G  !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(PressureForce_Mont_CS), intent(inout) :: CS !< Montgomery PGF control structure
  type(SAL_CS), intent(in), target, optional :: SAL_CSp !< SAL control structure
  type(tidal_forcing_CS), intent(in), target, optional :: tides_CSp !< Tides control structure

  ! Local variables
  logical :: use_EOS
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl   ! This module's name.

  CS%initialized = .true.
  CS%diag => diag ; CS%Time => Time
  if (present(tides_CSp)) &
    CS%tides_CSp => tides_CSp
  if (present(SAL_CSp)) &
    CS%SAL_CSp => SAL_CSp

  mdl = "MOM_PressureForce_Mont"
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%R_to_kg_m3)
  call get_param(param_file, mdl, "TIDES", CS%tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  call get_param(param_file, mdl, "CALCULATE_SAL", CS%calculate_SAL, &
                 "If true, calculate self-attraction and loading.", default=CS%tides)
  call get_param(param_file, mdl, "USE_EOS", use_EOS, default=.true., &
                 do_not_log=.true.) ! Input for diagnostic use only.

  if (use_EOS) then
    CS%id_PFu_bc = register_diag_field('ocean_model', 'PFu_bc', diag%axesCuL, Time, &
         'Density Gradient Zonal Pressure Force Accel.', "meter second-2", conversion=US%L_T2_to_m_s2)
    CS%id_PFv_bc = register_diag_field('ocean_model', 'PFv_bc', diag%axesCvL, Time, &
         'Density Gradient Meridional Pressure Force Accel.', "meter second-2", conversion=US%L_T2_to_m_s2)
    if (CS%id_PFu_bc > 0) &
      allocate(CS%PFu_bc(G%IsdB:G%IedB,G%jsd:G%jed,GV%ke), source=0.)
    if (CS%id_PFv_bc > 0) &
      allocate(CS%PFv_bc(G%isd:G%ied,G%JsdB:G%JedB,GV%ke), source=0.)
  endif

  if (CS%calculate_SAL) then
    CS%id_e_sal = register_diag_field('ocean_model', 'e_sal', diag%axesT1, Time, &
        'Self-attraction and loading height anomaly', 'meter', conversion=US%Z_to_m)
  endif
  if (CS%tides) then
    CS%id_e_tide = register_diag_field('ocean_model', 'e_tidal', diag%axesT1, Time, &
        'Tidal Forcing Astronomical and SAL Height Anomaly', 'meter', conversion=US%Z_to_m)
    CS%id_e_tide_eq  = register_diag_field('ocean_model', 'e_tide_eq', diag%axesT1, Time, &
        'Equilibrium tides height anomaly', 'meter', conversion=US%Z_to_m)
    CS%id_e_tide_sal = register_diag_field('ocean_model', 'e_tide_sal', diag%axesT1, Time, &
        'Read-in tidal self-attraction and loading height anomaly', 'meter', conversion=US%Z_to_m)
  endif

  CS%GFS_scale = 1.0
  if (GV%g_prime(1) /= GV%g_Earth) CS%GFS_scale = GV%g_prime(1) / GV%g_Earth

  call log_param(param_file, mdl, "GFS / G_EARTH", CS%GFS_scale, units="nondim")

end subroutine PressureForce_Mont_init

!>\namespace mom_pressureforce_mont
!!
!! Provides the Boussunesq and non-Boussinesq forms of the horizontal
!! accelerations due to pressure gradients using the Montgomery potential. A
!! second-order accurate, centered scheme is used. If a split time stepping
!! scheme is used, the vertical decomposition into barotropic and baroclinic
!! contributions described by Hallberg (J Comp Phys 1997) is used.  With a
!! nonlinear equation of state, compressibility is added along the lines proposed
!! by Sun et al. (JPO 1999), but with compressibility coefficients based on a fit
!! to a user-provided reference profile.

end module MOM_PressureForce_Mont
