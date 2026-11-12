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
use array_mod, only : RealArray_t, RealArray_C
use box_mod, only : Box_t
use turbotmp_helperF, only : io_recorder
use iso_c_binding, only : c_ptr, c_loc, c_f_pointer, c_bool, c_double

implicit none ; private

#include <MOM_memory.h>
#include "do_concurrent_compat.h"

public PressureForce_Mont_Bouss, PressureForce_Mont_nonBouss, Set_pbce_Bouss
public Set_pbce_nonBouss, PressureForce_Mont_init
public PressureForce_Mont_Bouss_pre, PressureForce_Mont_Bouss_post
public PressureForce_Mont_nonBouss_pre, PressureForce_Mont_nonBouss_post
public PF_Mont_opts_CS_C, PF_Mont_opts_CS_to_c, record_PF_Mont_opts_CS
public PF_vgrid_consts_CS_C, PF_vgrid_consts_build, record_PF_vgrid_consts_CS
#ifdef _TIM
public PF_vgrid_consts_CS_to_c
#endif

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Persistent containers of the vertical grid arrays used by the pressure force calculations.
type, public :: PF_vgrid_consts_CS
  type(RealArray_t) :: Rlay_a    !< The target coordinate value (potential density) in each layer [R ~> kg m-3].
  type(RealArray_t) :: g_prime_a !< The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
end type PF_vgrid_consts_CS

!> bind(C) mirror of PF_vgrid_consts_CS, field-for-field, same order.
type, bind(C) :: PF_vgrid_consts_CS_C
  type(RealArray_C) :: Rlay_a    !< The target coordinate value (potential density) in each layer [R ~> kg m-3].
  type(RealArray_C) :: g_prime_a !< The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
end type PF_vgrid_consts_CS_C

!> Options used by PressureForce_Mont_Bouss and PressureForce_Mont_nonBouss.
type, public :: PF_Mont_opts_CS
  logical :: calculate_SAL  !< If true, calculate self-attraction and loading.
  logical :: tides          !< If true, apply tidal momentum forcing.
  real    :: Rho0           !< The density used in the Boussinesq
                            !! approximation [R ~> kg m-3].
  real    :: GFS_scale      !< Ratio between gravity applied to top interface and the
                            !! gravitational acceleration of the planet [nondim].
                            !! Usually this ratio is 1.
end type PF_Mont_opts_CS

!> bind(C) mirror of PF_Mont_opts_CS, field-for-field, same order.
type, bind(C) :: PF_Mont_opts_CS_C
  logical(c_bool) :: calculate_SAL !< If true, calculate self-attraction and loading.
  logical(c_bool) :: tides         !< If true, apply tidal momentum forcing.
  real(c_double)  :: Rho0          !< The density used in the Boussinesq
                                   !! approximation [R ~> kg m-3].
  real(c_double)  :: GFS_scale     !< Ratio between gravity applied to top interface and the
                                   !! gravitational acceleration of the planet [nondim].
                                   !! Usually this ratio is 1.
end type PF_Mont_opts_CS_C

!> Control structure for the Montgomery potential form of pressure gradient
type, public :: PressureForce_Mont_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  type(PF_Mont_opts_CS), public :: PF_Mont_opts_CS !< Options used by PressureForce_Mont_Bouss and
                            !! PressureForce_Mont_nonBouss
  type(PF_vgrid_consts_CS), public :: PF_vgrid_consts_CS !< Persistent containers of the vertical grid
                            !! arrays used by the pressure force calculations
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate
                            !! the timing of diagnostic output.
  real, allocatable, public :: PFu_bc(:,:,:) !< Zonal accelerations due to pressure gradients
                            !! deriving from density gradients within layers [L T-2 ~> m s-2].
  real, allocatable, public :: PFv_bc(:,:,:) !< Meridional accelerations due to pressure gradients
                            !! deriving from density gradients within layers [L T-2 ~> m s-2].
  !>@{ Diagnostic IDs
  integer :: id_PFu_bc = -1, id_PFv_bc = -1, id_e_sal = -1
  integer :: id_e_tide = -1,  id_e_tide_eq = -1, id_e_tide_sal = -1
  !>@}
  type(SAL_CS), pointer :: SAL_CSp => NULL() !< SAL control structure
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL() !< The tidal forcing control structure
end type PressureForce_Mont_CS

!> The Fortran objects that the Set_pbce_Bouss and Set_pbce_nonBouss host helpers need. Each form of
!! the pressure force passes one to its Set_pbce call as an opaque c_ptr.
type, public :: pbce_host_ctx
  type(ocean_grid_type),   pointer :: G => NULL()  !< The ocean's grid structure.
  type(verticalGrid_type), pointer :: GV => NULL() !< The ocean's vertical grid structure.
  type(unit_scale_type),   pointer :: US => NULL() !< A dimensional unit scaling type
  type(thermo_var_ptrs),   pointer :: tv => NULL() !< The thermodynamic variables that Set_pbce uses, in which
                                                   !! T and S may be temporary arrays
  type(thermo_var_ptrs)   :: tv_tmp !< A structure of temporary T & S.
  real, allocatable       :: T_tmp(:,:,:) !< Temporary array of temperatures where layers that are lighter
                                          !! than the mixed layer have the mixed layer's properties [C ~> degC].
  real, allocatable       :: S_tmp(:,:,:) !< Temporary array of salinities where layers that are lighter
                                          !! than the mixed layer have the mixed layer's properties [S ~> ppt].
end type pbce_host_ctx

!> The Fortran objects that the Montgomery potential form's host helpers need. PressureForce builds
!! one and passes it down the tree as an opaque c_ptr.
type, public :: PF_Mont_host_ctx
  type(ocean_grid_type),       pointer :: G => NULL()  !< The ocean's grid structure.
  type(verticalGrid_type),     pointer :: GV => NULL() !< The ocean's vertical grid structure.
  type(unit_scale_type),       pointer :: US => NULL() !< A dimensional unit scaling type
  type(thermo_var_ptrs),       pointer :: tv => NULL() !< Thermodynamic variables.
  type(PressureForce_Mont_CS), pointer :: CS => NULL() !< Control structure for Montgomery potential PGF
  type(pbce_host_ctx) :: pbce !< The objects that the Set_pbce host helpers need
end type PF_Mont_host_ctx

contains

!> Converts a PF_Mont_opts_CS to its bind(C) mirror.
function PF_Mont_opts_CS_to_c(opts) result(cdesc)
  type(PF_Mont_opts_CS), intent(in) :: opts !< Options to convert
  type(PF_Mont_opts_CS_C) :: cdesc          !< bind(C) mirror of opts
  cdesc%calculate_SAL = opts%calculate_SAL
  cdesc%tides = opts%tides
  cdesc%Rho0 = opts%Rho0
  cdesc%GFS_scale = opts%GFS_scale
end function PF_Mont_opts_CS_to_c

!> Records the members of a PF_Mont_opts_CS in a capture file, with names starting with prefix.
subroutine record_PF_Mont_opts_CS(rec, prefix, opts)
  type(io_recorder),     intent(inout) :: rec    !< The capture file being written
  character(len=*),      intent(in)    :: prefix !< The start of each record name
  type(PF_Mont_opts_CS), intent(in)    :: opts   !< The options to record

  call rec%add(prefix//"_calculate_SAL", opts%calculate_SAL)
  call rec%add(prefix//"_tides", opts%tides)
  call rec%add(prefix//"_Rho0", opts%Rho0)
  call rec%add(prefix//"_GFS_scale", opts%GFS_scale)
end subroutine record_PF_Mont_opts_CS

!> Converts a PF_vgrid_consts_CS to its bind(C) mirror.
#ifdef _TIM
function PF_vgrid_consts_CS_to_c(vgrid) result(cdesc)
  type(PF_vgrid_consts_CS), intent(in) :: vgrid !< Vertical grid constants to convert
  type(PF_vgrid_consts_CS_C) :: cdesc           !< bind(C) mirror of vgrid
  cdesc%Rlay_a = vgrid%Rlay_a%to_c()
  cdesc%g_prime_a = vgrid%g_prime_a%to_c()
end function PF_vgrid_consts_CS_to_c
#endif

!> Records which arrays in a PF_vgrid_consts_CS are allocated (with an empty suffix), or the values
!! of the allocated ones (with a non-empty suffix), in a capture file.
subroutine record_PF_vgrid_consts_CS(rec, prefix, vgrid, suffix)
  type(io_recorder),        intent(inout) :: rec    !< The capture file being written
  character(len=*),         intent(in)    :: prefix !< The start of each record name
  type(PF_vgrid_consts_CS), intent(in)    :: vgrid  !< The vertical grid constants to record
  character(len=*),         intent(in)    :: suffix !< The end of each value record name,
                                                    !! or empty to record which arrays are
                                                    !! allocated

  if (len(suffix) == 0) then
    call rec%add(prefix//"_Rlay_allocated", vgrid%Rlay_a%associated())
  elseif (vgrid%Rlay_a%associated()) then
    call rec%add(prefix//"_Rlay"//suffix, vgrid%Rlay_a)
  endif
  if (len(suffix) == 0) then
    call rec%add(prefix//"_g_prime_allocated", vgrid%g_prime_a%associated())
  elseif (vgrid%g_prime_a%associated()) then
    call rec%add(prefix//"_g_prime"//suffix, vgrid%g_prime_a)
  endif
end subroutine record_PF_vgrid_consts_CS

!> Copies the vertical grid arrays used by the pressure force calculations into their persistent
!! containers.
subroutine PF_vgrid_consts_build(vgrid, GV)
  type(PF_vgrid_consts_CS), intent(inout) :: vgrid !< Persistent containers of GV%Rlay and GV%g_prime
  type(verticalGrid_type),  intent(in)    :: GV    !< Vertical grid structure

  if (allocated(GV%Rlay)) &
    call vgrid%Rlay_a%alloc(lb=LBOUND(GV%Rlay), ub=UBOUND(GV%Rlay), source=GV%Rlay)
  if (allocated(GV%g_prime)) &
    call vgrid%g_prime_a%alloc(lb=LBOUND(GV%g_prime), ub=UBOUND(GV%g_prime), source=GV%g_prime)
end subroutine PF_vgrid_consts_build

!> \brief Non-Boussinesq Montgomery-potential form of pressure gradient
!!
!! Determines the acceleration due to pressure forces in a
!! non-Boussinesq fluid using the compressibility compensated (if appropriate)
!! Montgomery-potential form described in Hallberg (Ocean Mod., 2005).
!!
!! To work, the following fields must be set outside of the usual (is:ie,js:je)
!! range before this subroutine is called:
!!   h(isB:ie+1,jsB:je+1), T(isB:ie+1,jsB:je+1), and S(isB:ie+1,jsB:je+1).
subroutine PressureForce_Mont_nonBouss(bxH0, bxQ0, bxU0, bxV0, bxHgrad, h_a, use_EOS, PFu_a, PFv_a, Gcore, g_Earth, &
                                       H_to_RZ, H_subroundoff, CS, vgrid, p_atm_a, host, e_sal_a, e_tide_eq_a, &
                                       e_tide_sal_a, pbce_a, eta_a, PFu_bc3d_a, PFv_bc3d_a)
  type(grid_core_type),                       intent(in)  :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),                                intent(in)  :: bxH0 !< The h-point box of the whole column,
                                                                  !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),                                intent(in)  :: bxQ0 !< The B-grid-index box of the whole column,
                                                                  !! [IscB:IecB, JscB:JecB, 1:nz]
  type(Box_t),                                intent(in)  :: bxU0 !< The u-point box of the whole column,
                                                                  !! [IscB:IecB, jsc:jec, 1:nz]
  type(Box_t),                                intent(in)  :: bxV0 !< The v-point box of the whole column,
                                                                  !! [isc:iec, JscB:JecB, 1:nz]
  type(Box_t),                                intent(in)  :: bxHgrad !< The h points on both sides of every u and v
                                                                     !! face, [IscB:IecB+1, JscB:JecB+1, 1:nz]
  real,                                       intent(in)  :: g_Earth !< The gravitational acceleration
                                                                     !! [L2 Z-1 T-2 ~> m s-2].
  real,                                       intent(in)  :: H_to_RZ !< A constant that translates thickness units to
                                                                     !! the units of mass per unit area
                                                                     !! [R Z H-1 ~> kg m-3 or 1].
  real,                                       intent(in)  :: H_subroundoff !< A thickness that is so small that it can
                                                                           !! be added to a thickness of Angstrom or
                                                                           !! larger without changing it at the bit
                                                                           !! level [H ~> m or kg m-2]. If Angstrom is
                                                                           !! 0 or exceedingly small, this is
                                                                           !! negligible compared to 1e-17 m.
  type(RealArray_t),                          intent(in)  :: h_a !< Layer thickness, [H ~> kg m-2].
  logical,                                    intent(in)  :: use_EOS !< If true, density is calculated from T & S
                                                                     !! using an equation of state.
  type(RealArray_t),                          intent(inout) :: PFu_a !< Zonal acceleration due to pressure gradients
                                                                     !! (equal to -dM/dx) [L T-2 ~> m s-2].
  type(RealArray_t),                          intent(inout) :: PFv_a !< Meridional acceleration due to pressure
                                                                     !! gradients (equal to -dM/dy) [L T-2 ~> m s-2].
  type(PF_Mont_opts_CS),                      intent(in)  :: CS  !< Montgomery potential PGF options
  type(PF_vgrid_consts_CS),                   intent(in)  :: vgrid !< Persistent containers of GV%Rlay and GV%g_prime
  type(RealArray_t),                          intent(in)  :: p_atm_a !< The pressure at the ice-ocean or
                                                                     !! atmosphere-ocean [R L2 T-2 ~> Pa].
  type(c_ptr),                         value, intent(in)  :: host !< Opaque pointer to the PF_Mont_host_ctx that holds
                                                                  !! the Fortran objects used by the host helpers
  type(RealArray_t), intent(inout) :: &
    e_sal_a     !< Bottom geopotential anomaly due to self-attraction and loading [Z ~> m].
  type(RealArray_t), intent(inout) :: &
    e_tide_eq_a !< Bottom geopotential anomaly due to tidal forces from astronomical sources [Z ~> m].
  type(RealArray_t), intent(inout) :: &
    e_tide_sal_a !< Bottom geopotential anomaly due to harmonic self-attraction and loading
                 !! specific to tides [Z ~> m].
  type(RealArray_t),                          intent(inout) :: pbce_a !< The baroclinic pressure anomaly in
                                                                 !! each layer due to free surface height anomalies,
                                                                 !! [L2 T-2 H-1 ~> m s-2 or m4 kg-1 s-2].
  type(RealArray_t),                          intent(inout) :: eta_a !< The total column mass used to calculate
                                                                 !! PFu and PFv [H ~> kg m-2].
  type(RealArray_t), intent(inout) :: &
    PFu_bc3d_a  !< Zonal accelerations due to pressure gradients
                !! deriving from density gradients within layers [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: &
    PFv_bc3d_a  !< Meridional accelerations due to pressure gradients
                !! deriving from density gradients within layers [L T-2 ~> m s-2].
  ! Local variables
  type(Box_t) :: bxH0_sk ! bxH0 with its k start moved in, [isc:iec, jsc:jec, 2:nz]
  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2),h_a%lb(3):h_a%ub(3)) :: &
    M, &          ! The Montgomery potential, M = (p/rho + gz)  [L2 T-2 ~> m2 s-2].
    dz_geo        !   The change in geopotential across a layer [L2 T-2 ~> m2 s-2].
  type(RealArray_t) :: alpha_star_a ! Compression adjusted specific volume [R-1 ~> m3 kg-1].
  type(RealArray_t) :: p_a ! Interface pressure [R L2 T-2 ~> Pa].
                ! p may be adjusted (with a nonlinear equation of state) so that
                ! its derivative compensates for the adiabatic compressibility
                ! in seawater, but p will still be close to the pressure.

  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2)) :: &
    dM, &         !   A barotropic correction to the Montgomery potentials to enable the use
                  ! of a reduced gravity form of the equations [L2 T-2 ~> m2 s-2].
    dp_star, &    ! Layer thickness after compensation for compressibility [R L2 T-2 ~> Pa].
    SSH, &        ! The sea surface height anomaly, in depth units [Z ~> m].
    geopot_bot    !   Bottom geopotential relative to a temporally fixed reference value,
                  ! including any tidal contributions [L2 T-2 ~> m2 s-2].
  real :: PFu_bc, PFv_bc     ! The pressure gradient force due to along-layer
                             ! compensated density gradients [L T-2 ~> m s-2]
  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [R L2 T-2 ~> Pa].
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: is_split        ! A flag indicating whether the pressure gradient terms are to be
                             ! split into barotropic and baroclinic pieces.
  type(PF_Mont_host_ctx), pointer :: ctx ! The objects behind host

  real :: I_gEarth           ! The inverse of g_Earth [T2 Z L-2 ~> s2 m-1]
!  real :: dalpha
  real :: Pa_to_H     ! A factor to convert from R L2 T-2 to the thickness units (H)
                      ! [H T2 R-1 L-2 ~> m2 s2 kg-1 or s2 m-1].
  real :: alpha_Lay(h_a%lb(3):h_a%ub(3)) ! The specific volume of each layer [R-1 ~> m3 kg-1].
  real :: dalpha_int(bxH0%idxS(3):bxH0%idxE(3)+1) ! The change in specific volume across each
                             ! interface [R-1 ~> m3 kg-1].
  real, dimension(:), contiguous, pointer :: Rlay ! The target coordinate value (potential density)
                                                  ! in each layer [R ~> kg m-3].
  real, dimension(:,:,:), contiguous, pointer :: pbce, PFu_bc3d, PFv_bc3d
  real, dimension(:,:,:), contiguous, pointer :: h, PFu, PFv, p, alpha_star
  real, dimension(:,:),   contiguous, pointer :: p_atm, e_sal, e_tide_eq, e_tide_sal
  real, dimension(:,:),   contiguous, pointer :: eta
  real, dimension(:,:),   contiguous, pointer :: IdxCu, IdyCv, bathyT, meanSL
  logical :: use_PFu_bc3d, use_PFv_bc3d
  integer :: Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k
  nz = bxH0%idxE(3)
  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2)
  bxH0_sk = bxH0%growLo(3,-1)

  nullify(Rlay)
  if (vgrid%Rlay_a%associated()) call vgrid%Rlay_a%view(Rlay)
  nullify(pbce, eta, PFu_bc3d, PFv_bc3d)
  if (pbce_a%associated()) call pbce_a%view(pbce)
  if (eta_a%associated()) call eta_a%view(eta)
  if (PFu_bc3d_a%associated()) call PFu_bc3d_a%view(PFu_bc3d)
  if (PFv_bc3d_a%associated()) call PFv_bc3d_a%view(PFv_bc3d)
  use_PFu_bc3d = PFu_bc3d_a%associated()
  use_PFv_bc3d = PFv_bc3d_a%associated()
  call h_a%view(h)
  call PFu_a%view(PFu)
  call PFv_a%view(PFv)
  nullify(p_atm)
  if (p_atm_a%associated()) call p_atm_a%view(p_atm)
  call e_sal_a%view(e_sal)
  call e_tide_eq_a%view(e_tide_eq)
  call e_tide_sal_a%view(e_tide_sal)
  call Gcore%IdxCu%view(IdxCu)
  call Gcore%IdyCv%view(IdyCv)
  call Gcore%bathyT%view(bathyT)
  call Gcore%meanSL%view(meanSL)
  call p_a%alloc(lb=[h_a%lb(1), h_a%lb(2), bxH0%idxS(3)], ub=[h_a%ub(1), h_a%ub(2), bxH0%idxE(3)+1])
  call p_a%view(p)
  call alpha_star_a%alloc(lb=[h_a%lb(1), h_a%lb(2), h_a%lb(3)], ub=[h_a%ub(1), h_a%ub(2), h_a%ub(3)])
  call alpha_star_a%view(alpha_star)

  use_p_atm = associated(p_atm)
  is_split = pbce_a%associated()

  I_gEarth = 1.0 / g_Earth
  dp_neglect = g_Earth * H_to_RZ * H_subroundoff
  do concurrent (k=bxH0%idxS(3):bxH0%idxE(3)) ; alpha_Lay(k) = 1.0 / (Rlay(k)) ; enddo
  do concurrent (k=bxH0_sk%idxS(3):bxH0_sk%idxE(3)) ; dalpha_int(K) = alpha_Lay(k-1) - alpha_Lay(k) ; enddo

  if (use_p_atm) then
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), &
                   i=bxHgrad%idxS(1):bxHgrad%idxE(1)) ; p(i,j,1) = p_atm(i,j) ; enddo
  else
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1)) ; p(i,j,1) = 0.0 ; enddo
  endif
  do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2)) DO_LOCALITY(local(k))
  do k=1,nz ; do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
    p(i,j,K+1) = p(i,j,K) + g_Earth * H_to_RZ * h(i,j,k)
  enddo ; enddo ; enddo

  if (eta_a%associated()) then
    Pa_to_H = 1.0 / (g_Earth * H_to_RZ)
    if (use_p_atm) then
      do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        eta(i,j) = (p(i,j,nz+1) - p_atm(i,j)) * Pa_to_H ! eta has the same units as h.
      enddo
    else
      do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        eta(i,j) = p(i,j,nz+1) * Pa_to_H ! eta has the same units as h.
      enddo
    endif
  endif

  do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
    geopot_bot(i,j) = -g_Earth * bathyT(i,j)
  enddo

  ! Calculate and add the self-attraction and loading geopotential anomaly.
  if (CS%calculate_SAL) then
    !   Determine the sea surface height anomalies, to enable the calculation
    ! of self-attraction and loading.
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
      SSH(i,j) = min(-bathyT(i,j) - meanSL(i,j), 0.0)
    enddo
    if (use_EOS) then
      call Mont_nonBouss_dz_geo_host(host, nz, p, dz_geo)
      do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2)) DO_LOCALITY(local(k))
      do k=1,nz ; do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        SSH(i,j) = SSH(i,j) + I_gEarth * dz_geo(i,j,k)
      enddo ; enddo ; enddo
    else
      do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2)) DO_LOCALITY(local(k))
      do k=1,nz ; do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        SSH(i,j) = SSH(i,j) + H_to_RZ * h(i,j,k) * alpha_Lay(k)
      enddo ; enddo ; enddo
    endif

    call Mont_calc_SAL_host(host, SSH, e_sal)
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
      geopot_bot(i,j) = geopot_bot(i,j) - g_Earth*e_sal(i,j)
    enddo
  endif

  ! Calculate and add the tidal geopotential anomaly.
  if (CS%tides) then
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
      geopot_bot(i,j) = geopot_bot(i,j) - g_Earth*(e_tide_eq(i,j) + e_tide_sal(i,j))
    enddo
  endif

  if (use_EOS) then
    call Mont_nonBouss_alpha_star_host(host, Isq, Ieq, Jsq, Jeq, nz, alpha_star)
  endif                                               ! use_EOS

  if (use_EOS) then
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2)) DO_LOCALITY(local(k))
      do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        M(i,j,nz) = geopot_bot(i,j) + p(i,j,nz+1) * alpha_star(i,j,nz)
      enddo
      do k=nz-1,1,-1 ; do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        M(i,j,k) = M(i,j,k+1) + p(i,j,K+1) * (alpha_star(i,j,k) - alpha_star(i,j,k+1))
      enddo ; enddo
    enddo
  else ! not use_EOS
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2)) DO_LOCALITY(local(k))
      do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        M(i,j,nz) = geopot_bot(i,j) + p(i,j,nz+1) * alpha_Lay(nz)
      enddo
      do k=nz-1,1,-1 ; do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        M(i,j,k) = M(i,j,k+1) + p(i,j,K+1) * dalpha_int(K+1)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (CS%GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
      dM(i,j) = (CS%GFS_scale - 1.0) * M(i,j,1)
    enddo
    do concurrent (k=bxHgrad%idxS(3):bxHgrad%idxE(3), j=bxHgrad%idxS(2):bxHgrad%idxE(2), &
                   i=bxHgrad%idxS(1):bxHgrad%idxE(1))
      M(i,j,k) = M(i,j,k) + dM(i,j)
    enddo

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
    call c_f_pointer(host, ctx)
    call Set_pbce_nonBouss(bxQ0, bxHgrad, p_a, use_EOS, g_Earth, H_to_RZ, H_subroundoff, vgrid, CS%GFS_scale, pbce_a, &
                           c_loc(ctx%pbce), alpha_star_a)
  endif

!    Calculate the pressure force. On a Cartesian grid,
!      PFu = - dM/dx   and  PFv = - dM/dy.
  if (use_EOS) then
    do k=1,nz
      do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        dp_star(i,j) = (p(i,j,K+1) - p(i,j,K)) + dp_neglect
      enddo
      do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), I=bxU0%idxS(1):bxU0%idxE(1)) DO_LOCALITY(local(PFu_bc))
        ! PFu_bc = p* grad alpha*
        PFu_bc = (alpha_star(i+1,j,k) - alpha_star(i,j,k)) * (IdxCu(I,j) * &
            ((dp_star(i,j)*dp_star(i+1,j) + ((p(i,j,K)*dp_star(i+1,j)) + (p(i+1,j,K)*dp_star(i,j)))) / &
             (dp_star(i,j) + dp_star(i+1,j))))
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * IdxCu(I,j) + PFu_bc
        if (use_PFu_bc3d) PFu_bc3d(i,j,k) = PFu_bc
      enddo
      do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), i=bxV0%idxS(1):bxV0%idxE(1)) DO_LOCALITY(local(PFv_bc))
        PFv_bc = (alpha_star(i,j+1,k) - alpha_star(i,j,k)) * (IdyCv(i,J) * &
            ((dp_star(i,j)*dp_star(i,j+1) + ((p(i,j,K)*dp_star(i,j+1)) + (p(i,j+1,K)*dp_star(i,j)))) / &
             (dp_star(i,j) + dp_star(i,j+1))))
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * IdyCv(i,J) + PFv_bc
        if (use_PFv_bc3d) PFv_bc3d(i,j,k) = PFv_bc
      enddo
    enddo ! k-loop
  else ! .not. use_EOS
    do concurrent (k=bxH0%idxS(3):bxH0%idxE(3))
      do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), I=bxU0%idxS(1):bxU0%idxE(1))
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * IdxCu(I,j)
      enddo
      do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), i=bxV0%idxS(1):bxV0%idxE(1))
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * IdyCv(i,J)
      enddo
    enddo
  endif ! use_EOS

  call p_a%free() ; call alpha_star_a%free()
  call bxH0_sk%free()
end subroutine PressureForce_Mont_nonBouss

!> The Fortran-only work PressureForce_Mont_nonBouss needs done before it: checks that the module is initialized
!! and that the equation of state is incompressible, and finds the tidal geopotential anomalies.
subroutine PressureForce_Mont_nonBouss_pre(G, US, CS, use_EOS, tv, e_tide_eq, e_tide_sal)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  type(PressureForce_Mont_CS),                intent(in)  :: CS !< Control structure for Montgomery potential PGF
  logical,                                    intent(in)  :: use_EOS !< If true, density is calculated from T & S
                                                                     !! using an equation of state.
  type(thermo_var_ptrs),                      intent(in)  :: tv !< Thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: &
    e_tide_eq   !< Bottom geopotential anomaly due to tidal forces from astronomical sources [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: &
    e_tide_sal  !< Bottom geopotential anomaly due to harmonic self-attraction and loading
                !! specific to tides [Z ~> m].

  if (.not.CS%initialized) call MOM_error(FATAL, &
      "MOM_PressureForce_Mont: Module must be initialized before it is used.")

  if (use_EOS) then
    if (query_compressible(tv%eqn_of_state)) call MOM_error(FATAL, &
      "PressureForce_Mont_nonBouss: The Montgomery form of the pressure force "//&
      "can no longer be used with a compressible EOS. Use #define ANALYTIC_FV_PGF.")
  endif

  if (CS%PF_Mont_opts_CS%tides) then
    call calc_tidal_forcing(CS%Time, e_tide_eq, e_tide_sal, G, US, CS%tides_CSp)
  endif

end subroutine PressureForce_Mont_nonBouss_pre

!> Posts the diagnostics of PressureForce_Mont_nonBouss, after it.
subroutine PressureForce_Mont_nonBouss_post(G, CS, e_sal, e_tide_eq, e_tide_sal)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(PressureForce_Mont_CS),                intent(in)  :: CS !< Control structure for Montgomery potential PGF
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: &
    e_sal       !< Bottom geopotential anomaly due to self-attraction and loading [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: &
    e_tide_eq   !< Bottom geopotential anomaly due to tidal forces from astronomical sources [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: &
    e_tide_sal  !< Bottom geopotential anomaly due to harmonic self-attraction and loading
                !! specific to tides [Z ~> m].

  if (CS%id_PFu_bc>0) call post_data(CS%id_PFu_bc, CS%PFu_bc, CS%diag)
  if (CS%id_PFv_bc>0) call post_data(CS%id_PFv_bc, CS%PFv_bc, CS%diag)
  ! To be consistent with old runs, tidal forcing diagnostic also includes total SAL.
  ! New diagnostics are given for each individual field.
  if (CS%id_e_tide>0) call post_data(CS%id_e_tide, e_sal+e_tide_eq+e_tide_sal, CS%diag)
  if (CS%id_e_sal>0) call post_data(CS%id_e_sal, e_sal, CS%diag)
  if (CS%id_e_tide_eq>0) call post_data(CS%id_e_tide_eq, e_tide_eq, CS%diag)
  if (CS%id_e_tide_sal>0) call post_data(CS%id_e_tide_sal, e_tide_sal, CS%diag)

end subroutine PressureForce_Mont_nonBouss_post

!> \brief Boussinesq Montgomery-potential form of pressure gradient
!!
!! Determines the acceleration due to pressure forces.
!!
!! To work, the following fields must be set outside of the usual (is:ie,js:je)
!! range before this subroutine is called:
!!   h(isB:ie+1,jsB:je+1), T(isB:ie+1,jsB:je+1), and S(isB:ie+1,jsB:je+1).
subroutine PressureForce_Mont_Bouss(bxH0, bxQ0, bxU0, bxV0, bxHgrad, h_a, use_EOS, PFu_a, PFv_a, Gcore, H_to_Z, &
                                    Z_to_H, dZ_subroundoff, CS, vgrid, p_atm_a, host, e_sal_a, e_tide_eq_a, &
                                    e_tide_sal_a, pbce_a, eta_a, PFu_bc3d_a, PFv_bc3d_a)
  type(grid_core_type),                       intent(in)  :: Gcore !< Persistent copies of the core ocean grid fields
  type(Box_t),                                intent(in)  :: bxH0 !< The h-point box of the whole column,
                                                                  !! [isc:iec, jsc:jec, 1:nz]
  type(Box_t),                                intent(in)  :: bxQ0 !< The B-grid-index box of the whole column,
                                                                  !! [IscB:IecB, JscB:JecB, 1:nz]
  type(Box_t),                                intent(in)  :: bxU0 !< The u-point box of the whole column,
                                                                  !! [IscB:IecB, jsc:jec, 1:nz]
  type(Box_t),                                intent(in)  :: bxV0 !< The v-point box of the whole column,
                                                                  !! [isc:iec, JscB:JecB, 1:nz]
  type(Box_t),                                intent(in)  :: bxHgrad !< The h points on both sides of every u and v
                                                                     !! face, [IscB:IecB+1, JscB:JecB+1, 1:nz]
  real,                                       intent(in)  :: H_to_Z !< A constant that translates thickness units to
                                                                    !! the units of depth [Z H-1 ~> 1 or m3 kg-1].
  real,                                       intent(in)  :: Z_to_H !< A constant that translates depth units to
                                                                    !! thickness units depth [H Z-1 ~> 1 or kg m-3].
  real,                                       intent(in)  :: dZ_subroundoff !< A thickness in height units that is so
                                                                            !! small that it can be added to a
                                                                            !! vertical distance of Angstrom_Z or
                                                                            !! 1e-17 m without changing it at the bit
                                                                            !! level [Z ~> m].  This is the height
                                                                            !! equivalent of H_subroundoff.
  type(RealArray_t),                          intent(in)  :: h_a !< Layer thickness [H ~> m].
  logical,                                    intent(in)  :: use_EOS !< If true, density is calculated from T & S
                                                                     !! using an equation of state.
  type(RealArray_t),                          intent(inout) :: PFu_a !< Zonal acceleration due to pressure gradients
                                                                     !! (equal to -dM/dx) [L T-2 ~> m s-2].
  type(RealArray_t),                          intent(inout) :: PFv_a !< Meridional acceleration due to pressure
                                                                     !! gradients (equal to -dM/dy) [L T-2 ~> m s-2].
  type(PF_Mont_opts_CS),                      intent(in)  :: CS  !< Montgomery potential PGF options
  type(PF_vgrid_consts_CS),                   intent(in)  :: vgrid !< Persistent containers of GV%Rlay and GV%g_prime
  type(RealArray_t),                          intent(in)  :: p_atm_a !< The pressure at the ice-ocean or
                                                                     !! atmosphere-ocean [R L2 T-2 ~> Pa].
  type(c_ptr),                         value, intent(in)  :: host !< Opaque pointer to the PF_Mont_host_ctx that holds
                                                                  !! the Fortran objects used by the host helpers
  type(RealArray_t), intent(inout) :: &
    e_sal_a     !< The bottom geopotential anomaly due to self-attraction and loading [Z ~> m].
  type(RealArray_t), intent(inout) :: &
    e_tide_eq_a !< Bottom geopotential anomaly due to tidal forces from astronomical sources
                !! [Z ~> m].
  type(RealArray_t), intent(inout) :: &
    e_tide_sal_a !< Bottom geopotential anomaly due to harmonic self-attraction and loading
                 !! specific to tides, in depth units [Z ~> m].
  type(RealArray_t),                          intent(inout) :: pbce_a !< The baroclinic pressure anomaly in
                                                                !! each layer due to free surface height anomalies
                                                                !! [L2 T-2 H-1 ~> m s-2].
  type(RealArray_t),                          intent(inout) :: eta_a !< Free surface height [H ~> m].
  type(RealArray_t), intent(inout) :: &
    PFu_bc3d_a  !< Zonal accelerations due to pressure gradients
                !! deriving from density gradients within layers [L T-2 ~> m s-2].
  type(RealArray_t), intent(inout) :: &
    PFv_bc3d_a  !< Meridional accelerations due to pressure gradients
                !! deriving from density gradients within layers [L T-2 ~> m s-2].
  ! Local variables
  real, dimension(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2),h_a%lb(3):h_a%ub(3)) :: &
    M           ! The Montgomery potential, M = (p/rho + gz) [L2 T-2 ~> m2 s-2].
  type(RealArray_t) :: rho_star_a ! In-situ density divided by the derivative with depth of the
                                  ! corrected e times (G_Earth/Rho0) [L2 Z-1 T-2 ~> m s-2].
  type(RealArray_t) :: e_a ! Interface height [Z ~> m].
                ! e may be adjusted (with a nonlinear equation of state) so that
                ! its derivative compensates for the adiabatic compressibility
                ! in seawater, but e will still be close to the interface depth.
  real :: h_star(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2)) ! Layer thickness after compensation
                             ! for compressibility [Z ~> m].
  real :: SSH(h_a%lb(1):h_a%ub(1),h_a%lb(2):h_a%ub(2)) ! The sea surface height anomaly, in depth units [Z ~> m].
  real :: I_Rho0             ! 1/Rho0 [R-1 ~> m3 kg-1].
  real :: PFu_bc, PFv_bc     ! The pressure gradient force due to along-layer
                             ! compensated density gradients [L T-2 ~> m s-2]
  real :: dz_neglect         ! A vertical distance that is so small it is usually lost
                             ! in roundoff and can be neglected [Z ~> m].
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: is_split        ! A flag indicating whether the pressure
                             ! gradient terms are to be split into
                             ! barotropic and baroclinic pieces.
  type(PF_Mont_host_ctx), pointer :: ctx ! The objects behind host
  real, dimension(:), contiguous, pointer :: g_prime ! The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
  real, dimension(:,:,:), contiguous, pointer :: pbce, PFu_bc3d, PFv_bc3d
  real, dimension(:,:,:), contiguous, pointer :: h, PFu, PFv, e, rho_star
  real, dimension(:,:),   contiguous, pointer :: p_atm, e_sal, e_tide_eq, e_tide_sal
  real, dimension(:,:),   contiguous, pointer :: eta
  real, dimension(:,:),   contiguous, pointer :: IdxCu, IdyCv, bathyT, meanSL
  logical :: use_PFu_bc3d, use_PFv_bc3d
  integer :: Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k

  nz = bxH0%idxE(3)
  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2)

  nullify(g_prime)
  if (vgrid%g_prime_a%associated()) call vgrid%g_prime_a%view(g_prime)
  nullify(pbce, eta, PFu_bc3d, PFv_bc3d)
  if (pbce_a%associated()) call pbce_a%view(pbce)
  if (eta_a%associated()) call eta_a%view(eta)
  if (PFu_bc3d_a%associated()) call PFu_bc3d_a%view(PFu_bc3d)
  if (PFv_bc3d_a%associated()) call PFv_bc3d_a%view(PFv_bc3d)
  use_PFu_bc3d = PFu_bc3d_a%associated()
  use_PFv_bc3d = PFv_bc3d_a%associated()
  call h_a%view(h)
  call PFu_a%view(PFu)
  call PFv_a%view(PFv)
  nullify(p_atm)
  if (p_atm_a%associated()) call p_atm_a%view(p_atm)
  call e_sal_a%view(e_sal)
  call e_tide_eq_a%view(e_tide_eq)
  call e_tide_sal_a%view(e_tide_sal)
  call Gcore%IdxCu%view(IdxCu)
  call Gcore%IdyCv%view(IdyCv)
  call Gcore%bathyT%view(bathyT)
  call Gcore%meanSL%view(meanSL)
  call e_a%alloc(lb=[h_a%lb(1), h_a%lb(2), bxH0%idxS(3)], ub=[h_a%ub(1), h_a%ub(2), bxH0%idxE(3)+1])
  call e_a%view(e)
  call rho_star_a%alloc(lb=[h_a%lb(1), h_a%lb(2), h_a%lb(3)], ub=[h_a%ub(1), h_a%ub(2), h_a%ub(3)])
  call rho_star_a%view(rho_star)

  use_p_atm = associated(p_atm)
  is_split = pbce_a%associated()

  dz_neglect = dZ_subroundoff
  I_Rho0 = 1.0/CS%Rho0

  do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
    e(i,j,nz+1) = -bathyT(i,j)
  enddo

  ! Calculate and add the self-attraction and loading geopotential anomaly.
  if (CS%calculate_SAL) then
    !   Determine the surface height anomaly for calculating self attraction
    ! and loading.  This should really be based on bottom pressure anomalies,
    ! but that is not yet implemented, and the current form is correct for
    ! barotropic tides.
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2)) DO_LOCALITY(local(k))
      do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1)) ; SSH(i,j) = min(-bathyT(i,j) - meanSL(i,j), 0.0) ; enddo
      do k=1,nz ; do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        SSH(i,j) = SSH(i,j) + h(i,j,k)*H_to_Z
      enddo ; enddo
    enddo
    call Mont_calc_SAL_host(host, SSH, e_sal)
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
      e(i,j,nz+1) = e(i,j,nz+1) - e_sal(i,j)
    enddo
  endif

  ! Calculate and add the tidal geopotential anomaly.
  if (CS%tides) then
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
      e(i,j,nz+1) = e(i,j,nz+1) - (e_tide_eq(i,j) + e_tide_sal(i,j))
    enddo
  endif

  do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2)) DO_LOCALITY(local(k))
  do k=nz,1,-1 ; do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
    e(i,j,K) = e(i,j,K+1) + h(i,j,k)*H_to_Z
  enddo ; enddo ; enddo

  if (use_EOS) then
    call Mont_Bouss_rho_star_host(host, Isq, Ieq, Jsq, Jeq, nz, rho_star)
  endif                                               ! use_EOS

!    Here the layer Montgomery potentials, M, are calculated.
  if (use_EOS) then
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2)) DO_LOCALITY(local(k))
      do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        M(i,j,1) = CS%GFS_scale * (rho_star(i,j,1) * e(i,j,1))
        if (use_p_atm) M(i,j,1) = M(i,j,1) + p_atm(i,j) * I_Rho0
      enddo
      do k=2,nz ; do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        M(i,j,k) = M(i,j,k-1) + (rho_star(i,j,k) - rho_star(i,j,k-1)) * e(i,j,K)
      enddo ; enddo
    enddo
  else ! not use_EOS
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2)) DO_LOCALITY(local(k))
      do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        M(i,j,1) = g_prime(1) * e(i,j,1)
        if (use_p_atm) M(i,j,1) = M(i,j,1) + p_atm(i,j) * I_Rho0
      enddo
      do k=2,nz ; do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        M(i,j,k) = M(i,j,k-1) + g_prime(K) * e(i,j,K)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (pbce_a%associated()) then
    call c_f_pointer(host, ctx)
    call Set_pbce_Bouss(bxQ0, bxHgrad, e_a, use_EOS, H_to_Z, dZ_subroundoff, vgrid, CS%Rho0, CS%GFS_scale, pbce_a, &
                        c_loc(ctx%pbce), rho_star_a, 1)
  endif

!    Calculate the pressure force. On a Cartesian grid,
!      PFu = - dM/dx   and  PFv = - dM/dy.
  if (use_EOS) then
    do k=1,nz
      do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        h_star(i,j) = (e(i,j,K) - e(i,j,K+1)) + dz_neglect
      enddo
      do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), I=bxU0%idxS(1):bxU0%idxE(1)) DO_LOCALITY(local(PFu_bc))
        PFu_bc = -1.0*(rho_star(i+1,j,k) - rho_star(i,j,k)) * (IdxCu(I,j) * &
          ((h_star(i,j) * h_star(i+1,j) - ((e(i,j,K) * h_star(i+1,j)) + &
          (e(i+1,j,K) * h_star(i,j)))) / (h_star(i,j) + h_star(i+1,j))))
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * IdxCu(I,j) + PFu_bc
        if (use_PFu_bc3d) PFu_bc3d(i,j,k) = PFu_bc
      enddo
      do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), i=bxV0%idxS(1):bxV0%idxE(1)) DO_LOCALITY(local(PFv_bc))
        PFv_bc = -1.0*(rho_star(i,j+1,k) - rho_star(i,j,k)) * (IdyCv(i,J) * &
          ((h_star(i,j) * h_star(i,j+1) - ((e(i,j,K) * h_star(i,j+1)) + &
          (e(i,j+1,K) * h_star(i,j)))) / (h_star(i,j) + h_star(i,j+1))))
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * IdyCv(i,J) + PFv_bc
        if (use_PFv_bc3d) PFv_bc3d(i,j,k) = PFv_bc
      enddo
    enddo ! k-loop
  else ! .not. use_EOS
    do concurrent (k=bxH0%idxS(3):bxH0%idxE(3))
      do concurrent (j=bxU0%idxS(2):bxU0%idxE(2), I=bxU0%idxS(1):bxU0%idxE(1))
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * IdxCu(I,j)
      enddo
      do concurrent (J=bxV0%idxS(2):bxV0%idxE(2), i=bxV0%idxS(1):bxV0%idxE(1))
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * IdyCv(i,J)
      enddo
    enddo
  endif ! use_EOS

  if (eta_a%associated()) then
    ! eta is the sea surface height relative to a time-invariant geoid, for
    ! comparison with what is used for eta in btstep.  See how e was calculated
    ! about 200 lines above.
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
      eta(i,j) = e(i,j,1)*Z_to_H
    enddo
    if (CS%tides) then
      do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        eta(i,j) = eta(i,j) + (e_tide_eq(i,j)+e_tide_sal(i,j))*Z_to_H
      enddo
    endif
    if (CS%calculate_SAL) then
      do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        eta(i,j) = eta(i,j) + e_sal(i,j)*Z_to_H
      enddo
    endif
  endif

  call e_a%free() ; call rho_star_a%free()
end subroutine PressureForce_Mont_Bouss

!> The Fortran-only work PressureForce_Mont_Bouss needs done before it: checks that the module is initialized
!! and that the equation of state is incompressible, and finds the tidal geopotential anomalies.
subroutine PressureForce_Mont_Bouss_pre(G, US, CS, use_EOS, tv, e_tide_eq, e_tide_sal)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  type(PressureForce_Mont_CS),                intent(in)  :: CS !< Control structure for Montgomery potential PGF
  logical,                                    intent(in)  :: use_EOS !< If true, density is calculated from T & S
                                                                     !! using an equation of state.
  type(thermo_var_ptrs),                      intent(in)  :: tv !< Thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: &
    e_tide_eq   !< Bottom geopotential anomaly due to tidal forces from astronomical sources
                !! [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: &
    e_tide_sal  !< Bottom geopotential anomaly due to harmonic self-attraction and loading
                !! specific to tides, in depth units [Z ~> m].

  if (.not.CS%initialized) call MOM_error(FATAL, &
       "MOM_PressureForce_Mont: Module must be initialized before it is used.")

  if (use_EOS) then
    if (query_compressible(tv%eqn_of_state)) call MOM_error(FATAL, &
      "PressureForce_Mont_Bouss: The Montgomery form of the pressure force "//&
      "can no longer be used with a compressible EOS. Use #define ANALYTIC_FV_PGF.")
  endif

  if (CS%PF_Mont_opts_CS%tides) then
    call calc_tidal_forcing(CS%Time, e_tide_eq, e_tide_sal, G, US, CS%tides_CSp)
  endif

end subroutine PressureForce_Mont_Bouss_pre

!> Posts the diagnostics of PressureForce_Mont_Bouss, after it.
subroutine PressureForce_Mont_Bouss_post(G, CS, e_sal, e_tide_eq, e_tide_sal)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(PressureForce_Mont_CS),                intent(in)  :: CS !< Control structure for Montgomery potential PGF
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: &
    e_sal       !< The bottom geopotential anomaly due to self-attraction and loading [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: &
    e_tide_eq   !< Bottom geopotential anomaly due to tidal forces from astronomical sources
                !! [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: &
    e_tide_sal  !< Bottom geopotential anomaly due to harmonic self-attraction and loading
                !! specific to tides, in depth units [Z ~> m].

  if (CS%id_PFu_bc>0) call post_data(CS%id_PFu_bc, CS%PFu_bc, CS%diag)
  if (CS%id_PFv_bc>0) call post_data(CS%id_PFv_bc, CS%PFv_bc, CS%diag)
  ! To be consistent with old runs, tidal forcing diagnostic also includes total SAL.
  ! New diagnostics are given for each individual field.
  if (CS%id_e_tide>0) call post_data(CS%id_e_tide, e_sal+e_tide_eq+e_tide_sal, CS%diag)
  if (CS%id_e_sal>0) call post_data(CS%id_e_sal, e_sal, CS%diag)
  if (CS%id_e_tide_eq>0) call post_data(CS%id_e_tide_eq, e_tide_eq, CS%diag)
  if (CS%id_e_tide_sal>0) call post_data(CS%id_e_tide_sal, e_tide_sal, CS%diag)

end subroutine PressureForce_Mont_Bouss_post

!> Determines the partial derivative of the acceleration due
!! to pressure forces with the free surface height.
subroutine Set_pbce_Bouss(bxQ0, bxHgrad, e_a, use_EOS, H_to_Z, dZ_subroundoff, vgrid, Rho0, GFS_scale, pbce_a, host, &
                          rho_star_a, nkblock)
  type(Box_t),                          intent(in)  :: bxQ0 !< The B-grid-index box of the whole column,
                                                            !! [IscB:IecB, JscB:JecB, 1:nz]
  type(Box_t),                          intent(in)  :: bxHgrad !< The h points on both sides of every u and v face,
                                                               !! [IscB:IecB+1, JscB:JecB+1, 1:nz]
  real,                                 intent(in)  :: H_to_Z !< A constant that translates thickness units to the
                                                              !! units of depth [Z H-1 ~> 1 or m3 kg-1].
  real,                                 intent(in)  :: dZ_subroundoff !< A thickness in height units that is so small
                                                                      !! that it can be added to a vertical distance
                                                                      !! of Angstrom_Z or 1e-17 m without changing it
                                                                      !! at the bit level [Z ~> m].  This is the
                                                                      !! height equivalent of H_subroundoff.
  type(RealArray_t),                    intent(in)  :: e_a !< Interface height [Z ~> m].
  logical,                              intent(in)  :: use_EOS !< If true, density is calculated from T & S using
                                                               !! an equation of state.
  type(PF_vgrid_consts_CS),             intent(in)  :: vgrid !< Persistent containers of GV%Rlay and GV%g_prime
  real,                                 intent(in)  :: Rho0 !< The "Boussinesq" ocean density [R ~> kg m-3].
  real,                                 intent(in)  :: GFS_scale !< Ratio between gravity applied to top
                                                            !! interface and the gravitational acceleration of
                                                            !! the planet [nondim]. Usually this ratio is 1.
  type(RealArray_t),                    intent(inout) :: pbce_a !< The baroclinic pressure anomaly in each layer due
                                                                !! to free surface height anomalies
                                                                !! [L2 T-2 H-1 ~> m s-2].
  type(c_ptr),                   value, intent(in)  :: host !< Opaque pointer to the pbce_host_ctx that holds the
                                                            !! Fortran objects used by the host helpers
  type(RealArray_t),                    intent(in)  :: rho_star_a !< The layer densities (maybe compressibility
                                                            !! compensated), times g/rho_0 [L2 Z-1 T-2 ~> m s-2].
  integer, intent(in) :: nkblock !< The number of vertical levels in each block for the variable density
                                 !! near-surface layer.  This is needed to determine whether rho_star is
                                 !! compressibility compensated or not.
  real :: Ihtot(e_a%lb(1):e_a%ub(1), &
                e_a%lb(2):e_a%ub(2))  ! The inverse of the sum of the layer thicknesses [H-1 ~> m-1 or m2 kg-1].
  real :: dz_neglect         ! A vertical distance that is so small it is usually lost
                             ! in roundoff and can be neglected [Z ~> m].
  real, dimension(:), contiguous, pointer :: g_prime ! The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
  real, dimension(:,:,:), contiguous, pointer :: rho_star
  real, dimension(:,:,:), contiguous, pointer :: e, pbce
  logical :: use_rho_star
  integer :: Isq, Ieq, Jsq, Jeq, nz, i, j, k

  nullify(g_prime)
  if (vgrid%g_prime_a%associated()) call vgrid%g_prime_a%view(g_prime)
  nullify(rho_star)
  if (rho_star_a%associated()) call rho_star_a%view(rho_star)
  use_rho_star = rho_star_a%associated()
  call e_a%view(e)
  call pbce_a%view(pbce)

  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2) ; nz = bxQ0%idxE(3)

  dz_neglect = dZ_subroundoff

  if (use_EOS) then
    if (use_rho_star) then
      do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2))
        do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
          Ihtot(i,j) = H_to_Z / ((e(i,j,1) - e(i,j,nz+1)) + dz_neglect)
          pbce(i,j,1) = GFS_scale * rho_star(i,j,1) * H_to_Z
        enddo

        do k=2,nz
          do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
            pbce(i,j,k) = pbce(i,j,k-1) + (rho_star(i,j,k) - rho_star(i,j,k-1)) &
                * ((e(i,j,K) - e(i,j,nz+1)) * Ihtot(i,j))
          enddo
        enddo
      enddo
    else
      call Set_pbce_Bouss_EOS_host(host, Isq, Ieq, Jsq, Jeq, nz, nkblock, e, Rho0, GFS_scale, &
                                   dz_neglect, Ihtot, pbce)
    endif
  else
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2))
      do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        Ihtot(i,j) = 1.0 / ((e(i,j,1) - e(i,j,nz+1)) + dz_neglect)
        pbce(i,j,1) = g_prime(1) * H_to_Z
      enddo
      do k=2,nz
        do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
          pbce(i,j,k) = pbce(i,j,k-1) + (g_prime(K) * H_to_Z) &
              * ((e(i,j,K) - e(i,j,nz+1)) * Ihtot(i,j))
        enddo
      enddo
    enddo
  endif
end subroutine Set_pbce_Bouss

!> Determines the partial derivative of the acceleration due
!! to pressure forces with the column mass.
subroutine Set_pbce_nonBouss(bxQ0, bxHgrad, p_a, use_EOS, g_Earth, H_to_RZ, H_subroundoff, vgrid, GFS_scale, pbce_a, &
                             host, alpha_star_a)
  type(Box_t),                          intent(in)  :: bxQ0 !< The B-grid-index box of the whole column,
                                                            !! [IscB:IecB, JscB:JecB, 1:nz]
  type(Box_t),                          intent(in)  :: bxHgrad !< The h points on both sides of every u and v face,
                                                               !! [IscB:IecB+1, JscB:JecB+1, 1:nz]
  real,                                 intent(in)  :: g_Earth !< The gravitational acceleration
                                                               !! [L2 Z-1 T-2 ~> m s-2].
  real,                                 intent(in)  :: H_to_RZ !< A constant that translates thickness units to the
                                                               !! units of mass per unit area
                                                               !! [R Z H-1 ~> kg m-3 or 1].
  real,                                 intent(in)  :: H_subroundoff !< A thickness that is so small that it can be
                                                                     !! added to a thickness of Angstrom or larger
                                                                     !! without changing it at the bit level
                                                                     !! [H ~> m or kg m-2]. If Angstrom is 0 or
                                                                     !! exceedingly small, this is negligible compared
                                                                     !! to 1e-17 m.
  type(RealArray_t),                    intent(in)  :: p_a !< Interface pressures [R L2 T-2 ~> Pa].
  logical,                              intent(in)  :: use_EOS !< If true, density is calculated from T & S using an
                                                               !! equation of state.
  type(PF_vgrid_consts_CS),             intent(in)  :: vgrid !< Persistent containers of GV%Rlay and GV%g_prime
  real,                                 intent(in)  :: GFS_scale !< Ratio between gravity applied to top
                                                          !! interface and the gravitational acceleration of
                                                          !! the planet [nondim]. Usually this ratio is 1.
  type(RealArray_t),                    intent(inout) :: pbce_a !< The baroclinic pressure anomaly in each
                                                                !! layer due to free surface height anomalies
                                                                !! [L2 H-1 T-2 ~> m4 kg-1 s-2].
  type(c_ptr),                   value, intent(in)  :: host !< Opaque pointer to the pbce_host_ctx that holds the
                                                            !! Fortran objects used by the host helpers
  type(RealArray_t),                    intent(in)  :: alpha_star_a !< The layer specific volumes
                                                          !! (maybe compressibility compensated) [R-1 ~> m3 kg-1].
  ! Local variables
  type(Box_t) :: bxQ0_sk ! bxQ0 with its k start moved in, [IscB:IecB, JscB:JecB, 2:nz]
  real, dimension(p_a%lb(1):p_a%ub(1),p_a%lb(2):p_a%ub(2)) :: &
    dpbce, &      !   A barotropic correction to the pbce to enable the use of
                  ! a reduced gravity form of the equations [L2 H-1 T-2 ~> m4 kg-1 s-2].
    C_htot        ! dP_dH divided by the total ocean pressure [H-1 ~> m2 kg-1].
  real :: alpha_Lay(pbce_a%lb(3):pbce_a%ub(3)) ! The specific volume of each layer [R-1 ~> m3 kg-1].
  real :: dalpha_int(p_a%lb(3):p_a%ub(3)) ! The change in specific volume across each interface [R-1 ~> m3 kg-1].
  real :: dP_dH              ! A factor that converts from thickness to pressure times other dimensional
                             ! conversion factors [R L2 T-2 H-1 ~> Pa m2 kg-1].
  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [R L2 T-2 ~> Pa].
  real, dimension(:), contiguous, pointer :: Rlay ! The target coordinate value (potential density)
                                                  ! in each layer [R ~> kg m-3].
  real, dimension(:,:,:), contiguous, pointer :: alpha_star
  real, dimension(:,:,:), contiguous, pointer :: p, pbce
  integer :: Isq, Ieq, Jsq, Jeq, nz, i, j, k

  Isq = bxQ0%idxS(1) ; Ieq = bxQ0%idxE(1) ; Jsq = bxQ0%idxS(2) ; Jeq = bxQ0%idxE(2) ; nz = bxQ0%idxE(3)
  bxQ0_sk = bxQ0%growLo(3,-1)

  nullify(Rlay)
  if (vgrid%Rlay_a%associated()) call vgrid%Rlay_a%view(Rlay)
  nullify(alpha_star)
  if (alpha_star_a%associated()) call alpha_star_a%view(alpha_star)
  call p_a%view(p)
  call pbce_a%view(pbce)

  dP_dH = g_Earth * H_to_RZ
  dp_neglect = g_Earth * H_to_RZ * H_subroundoff

  if (use_EOS) then
    if (alpha_star_a%associated()) then
      do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2)) DO_LOCALITY(local(k))
        do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
          C_htot(i,j) = dP_dH / ((p(i,j,nz+1)-p(i,j,1)) + dp_neglect)
          pbce(i,j,nz) = dP_dH * alpha_star(i,j,nz)
        enddo
        do k=nz-1,1,-1 ; do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
          pbce(i,j,k) = pbce(i,j,k+1) + ((p(i,j,K+1)-p(i,j,1)) * C_htot(i,j)) * &
              (alpha_star(i,j,k) - alpha_star(i,j,k+1))
        enddo ; enddo
      enddo
    else
      call Set_pbce_nonBouss_EOS_host(host, Isq, Ieq, Jsq, Jeq, nz, p, dP_dH, dp_neglect, C_htot, pbce)
    endif
  else ! not use_EOS

    do concurrent (k=bxQ0%idxS(3):bxQ0%idxE(3)) ; alpha_Lay(k) = 1.0 / (Rlay(k)) ; enddo
    do concurrent (k=bxQ0_sk%idxS(3):bxQ0_sk%idxE(3)) ; dalpha_int(K) = alpha_Lay(k-1) - alpha_Lay(k) ; enddo

    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2)) DO_LOCALITY(local(k))
      do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        C_htot(i,j) = dP_dH / ((p(i,j,nz+1)-p(i,j,1)) + dp_neglect)
        pbce(i,j,nz) = dP_dH * alpha_Lay(nz)
      enddo
      do k=nz-1,1,-1 ; do concurrent (i=bxHgrad%idxS(1):bxHgrad%idxE(1))
        pbce(i,j,k) = pbce(i,j,k+1) + ((p(i,j,K+1)-p(i,j,1))*C_htot(i,j)) * dalpha_int(K+1)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    do concurrent (j=bxHgrad%idxS(2):bxHgrad%idxE(2), i=bxHgrad%idxS(1):bxHgrad%idxE(1))
      dpbce(i,j) = (GFS_scale - 1.0) * pbce(i,j,1)
    enddo
    do concurrent (k=bxHgrad%idxS(3):bxHgrad%idxE(3), j=bxHgrad%idxS(2):bxHgrad%idxE(2), &
                   i=bxHgrad%idxS(1):bxHgrad%idxE(1))
      pbce(i,j,k) = pbce(i,j,k) + dpbce(i,j)
    enddo
  endif
  call bxQ0_sk%free()

end subroutine Set_pbce_nonBouss

!> Calls Set_pbce_Bouss_EOS with the grid, vertical grid and thermodynamic variables in the
!! pbce_host_ctx that host points to.
subroutine Set_pbce_Bouss_EOS_host(host, Isq, Ieq, Jsq, Jeq, nz, nkblock, e, Rho0, GFS_scale, &
                                   dz_neglect, Ihtot, pbce)
  type(c_ptr),     value, intent(in)    :: host !< Opaque pointer to the pbce_host_ctx that holds the
                                                !! Fortran objects used by the host helpers
  integer,                intent(in)    :: Isq !< The start i-index of cell vertices within the computational domain
  integer,                intent(in)    :: Ieq !< The end i-index of cell vertices within the computational domain
  integer,                intent(in)    :: Jsq !< The start j-index of cell vertices within the computational domain
  integer,                intent(in)    :: Jeq !< The end j-index of cell vertices within the computational domain
  integer,                intent(in)    :: nz !< The number of layers/levels in the vertical
  integer,                intent(in)    :: nkblock !< The number of vertical levels in each block for the variable
                                                   !! density near-surface layer.  This is needed to determine whether
                                                   !! rho_star is compressibility compensated or not.
  real, dimension(:,:,:), intent(in)    :: e !< Interface height [Z ~> m].
  real,                   intent(in)    :: Rho0 !< The "Boussinesq" ocean density [R ~> kg m-3].
  real,                   intent(in)    :: GFS_scale !< Ratio between gravity applied to top
                                                     !! interface and the gravitational acceleration of
                                                     !! the planet [nondim]. Usually this ratio is 1.
  real,                   intent(in)    :: dz_neglect !< A vertical distance that is so small it is usually lost
                                                      !! in roundoff and can be neglected [Z ~> m].
  real, dimension(:,:),   intent(inout) :: Ihtot !< The inverse of the sum of the layer thicknesses
                                                 !! [H-1 ~> m-1 or m2 kg-1].
  real, dimension(:,:,:), intent(out)   :: pbce !< The baroclinic pressure anomaly in each layer due
                                                !! to free surface height anomalies
                                                !! [L2 T-2 H-1 ~> m s-2].
  ! Local variables
  type(pbce_host_ctx), pointer :: ctx ! The objects behind host

  call c_f_pointer(host, ctx)
  call Set_pbce_Bouss_EOS(ctx%G, ctx%GV, ctx%tv, Isq, Ieq, Jsq, Jeq, nz, nkblock, e, Rho0, GFS_scale, &
                          dz_neglect, Ihtot, pbce)

end subroutine Set_pbce_Bouss_EOS_host

!> Determines pbce for Set_pbce_Bouss from the equation of state, when the layer densities are
!! not provided.
subroutine Set_pbce_Bouss_EOS(G, GV, tv, Isq, Ieq, Jsq, Jeq, nz, nkblock, e, Rho0, GFS_scale, &
                              dz_neglect, Ihtot, pbce)
  type(ocean_grid_type),                intent(in)  :: G    !< Ocean grid structure
  type(verticalGrid_type),              intent(in)  :: GV   !< Vertical grid structure
  type(thermo_var_ptrs),                intent(in)  :: tv   !< Thermodynamic variables
  integer,                              intent(in)  :: Isq !< The start i-index of cell vertices within the
                                                           !! computational domain
  integer,                              intent(in)  :: Ieq !< The end i-index of cell vertices within the
                                                           !! computational domain
  integer,                              intent(in)  :: Jsq !< The start j-index of cell vertices within the
                                                           !! computational domain
  integer,                              intent(in)  :: Jeq !< The end j-index of cell vertices within the
                                                           !! computational domain
  integer,                              intent(in)  :: nz !< The number of layers/levels in the vertical
  integer, intent(in) :: nkblock !< The number of vertical levels in each block for the variable density
                                 !! near-surface layer.  This is needed to determine whether rho_star is
                                 !! compressibility compensated or not.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: e !< Interface height [Z ~> m].
  real,                                 intent(in)  :: Rho0 !< The "Boussinesq" ocean density [R ~> kg m-3].
  real,                                 intent(in)  :: GFS_scale !< Ratio between gravity applied to top
                                                            !! interface and the gravitational acceleration of
                                                            !! the planet [nondim]. Usually this ratio is 1.
  real,                                 intent(in)  :: dz_neglect !< A vertical distance that is so small it is
                                                                  !! usually lost in roundoff and can be neglected
                                                                  !! [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)),     intent(inout) :: Ihtot !< The inverse of the sum of the layer thicknesses
                                                               !! [H-1 ~> m-1 or m2 kg-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                        intent(out) :: pbce !< The baroclinic pressure anomaly in each layer due
                                                            !! to free surface height anomalies
                                                            !! [L2 T-2 H-1 ~> m s-2].
  ! Local variables
  real :: press(SZI_(G),SZJ_(G),nkblock)  ! Interface pressure [R L2 T-2 ~> Pa].
  real :: T_int(SZI_(G),SZJ_(G),nkblock)  ! Interface temperature [C ~> degC]
  real :: S_int(SZI_(G),SZJ_(G),nkblock)  ! Interface salinity [S ~> ppt]
  real :: dR_dT(SZI_(G),SZJ_(G),nkblock)  ! Partial derivative of density with temperature [R C-1 ~> kg m-3 degC-1]
  real :: dR_dS(SZI_(G),SZJ_(G),nkblock)  ! Partial derivative of density with salinity [R S-1 ~> kg m-3 ppt-1].
  real :: rho_in_situ(SZI_(G),SZJ_(G))  ! In-situ density at the top of a layer [R ~> kg m-3].
  real :: G_Rho0             ! A scaled version of g_Earth / Rho0 [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]
  real :: Rho0xG             ! g_Earth * Rho0 [R L2 Z-1 T-2 ~> kg s-2 m-2]
  integer :: EOSdom(3,2)     ! The computational domain for the equation of state
  integer :: i, j, k, kstart, kend, kk

  !$omp target data &
  !$omp   map(alloc: EOSdom, press, T_int, S_int, rho_in_situ) &
  !$omp   map(alloc: dR_dT, dR_dS)

  Rho0xG = Rho0 * GV%g_Earth
  G_Rho0 = GV%g_Earth / GV%Rho0

  EOSdom(1,:) = [Isq - (G%isd-1), G%iec+1 - (G%isd-1)]
  EOSdom(2,:) = [Jsq - (G%jsd-1), G%jec+1 - (G%jsd-1)]
  do concurrent (j=Jsq:Jeq+1, i=Isq:Ieq+1)
    Ihtot(i,j) = GV%H_to_Z / ((e(i,j,1) - e(i,j,nz+1)) + dz_neglect)
    press(i,j,1) = -Rho0xG * (e(i,j,1) - G%meanSL(i,j))
  enddo

  call calculate_density(tv%T(:,:,1), tv%S(:,:,1), press(:,:,1), rho_in_situ, &
                         tv%eqn_of_state, EOSdom(1:2,:))

  do concurrent (j=Jsq:Jeq+1, i=Isq:Ieq+1)
    pbce(i,j,1) = G_Rho0 * (GFS_scale * rho_in_situ(i,j)) * GV%H_to_Z
  enddo

  EOSdom(3,1) = 1
  do kstart=2,nz,nkblock
    kend = min(kstart+nkblock-1, nz)
    do concurrent (k=kstart:kend, j=Jsq:Jeq+1, i=Isq:Ieq+1)
      kk=k-kstart+1
      press(i,j,kk) = -Rho0xG * (e(i,j,K) - G%meanSL(i,j))
      T_int(i,j,kk) = 0.5 * (tv%T(i,j,k-1) + tv%T(i,j,k))
      S_int(i,j,kk) = 0.5 * (tv%S(i,j,k-1) + tv%S(i,j,k))
    enddo
    EOSdom(3,2) = kend-kstart+1
    call calculate_density_derivs(T_int, S_int, press, dR_dT, dR_dS, &
                                  tv%eqn_of_state, EOSdom)

    do concurrent(j=Jsq:Jeq+1)
      do k=kstart,kend ; do concurrent (i=Isq:Ieq+1)
      kk=k-kstart+1
      pbce(i,j,k) = pbce(i,j,k-1) + G_Rho0 * &
         ((e(i,j,K) - e(i,j,nz+1)) * Ihtot(i,j)) * &
         (dR_dT(i,j,kk) * (tv%T(i,j,k) - tv%T(i,j,k-1)) + &
          dR_dS(i,j,kk) * (tv%S(i,j,k) - tv%S(i,j,k-1)))
      enddo ; enddo
    enddo
  enddo
  !$omp end target data
end subroutine Set_pbce_Bouss_EOS

!> Calls Set_pbce_nonBouss_EOS with the grid, vertical grid and thermodynamic variables in the
!! pbce_host_ctx that host points to.
subroutine Set_pbce_nonBouss_EOS_host(host, Isq, Ieq, Jsq, Jeq, nz, p, dP_dH, dp_neglect, C_htot, pbce)
  type(c_ptr),     value, intent(in)    :: host !< Opaque pointer to the pbce_host_ctx that holds the
                                                !! Fortran objects used by the host helpers
  integer,                intent(in)    :: Isq !< The start i-index of cell vertices within the computational domain
  integer,                intent(in)    :: Ieq !< The end i-index of cell vertices within the computational domain
  integer,                intent(in)    :: Jsq !< The start j-index of cell vertices within the computational domain
  integer,                intent(in)    :: Jeq !< The end j-index of cell vertices within the computational domain
  integer,                intent(in)    :: nz !< The number of layers/levels in the vertical
  real, dimension(:,:,:), intent(in)    :: p !< Interface pressures [R L2 T-2 ~> Pa].
  real,                   intent(in)    :: dP_dH !< A factor that converts from thickness to pressure times other
                                                 !! dimensional conversion factors [R L2 T-2 H-1 ~> Pa m2 kg-1].
  real,                   intent(in)    :: dp_neglect !< A thickness that is so small it is usually lost
                                                      !! in roundoff and can be neglected [R L2 T-2 ~> Pa].
  real, dimension(:,:),   intent(inout) :: C_htot !< dP_dH divided by the total ocean pressure [H-1 ~> m2 kg-1].
  real, dimension(:,:,:), intent(out)   :: pbce !< The baroclinic pressure anomaly in each
                                                !! layer due to free surface height anomalies
                                                !! [L2 H-1 T-2 ~> m4 kg-1 s-2].
  ! Local variables
  type(pbce_host_ctx), pointer :: ctx ! The objects behind host

  call c_f_pointer(host, ctx)
  call Set_pbce_nonBouss_EOS(ctx%G, ctx%GV, ctx%tv, Isq, Ieq, Jsq, Jeq, nz, p, dP_dH, dp_neglect, &
                             C_htot, pbce)

end subroutine Set_pbce_nonBouss_EOS_host

!> Determines pbce for Set_pbce_nonBouss from the equation of state, when the layer specific
!! volumes are not provided.
subroutine Set_pbce_nonBouss_EOS(G, GV, tv, Isq, Ieq, Jsq, Jeq, nz, p, dP_dH, dp_neglect, C_htot, pbce)
  type(ocean_grid_type),                intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),              intent(in)  :: GV !< Vertical grid structure
  type(thermo_var_ptrs),                intent(in)  :: tv !< Thermodynamic variables
  integer,                              intent(in)  :: Isq !< The start i-index of cell vertices within the
                                                           !! computational domain
  integer,                              intent(in)  :: Ieq !< The end i-index of cell vertices within the
                                                           !! computational domain
  integer,                              intent(in)  :: Jsq !< The start j-index of cell vertices within the
                                                           !! computational domain
  integer,                              intent(in)  :: Jeq !< The end j-index of cell vertices within the
                                                           !! computational domain
  integer,                              intent(in)  :: nz !< The number of layers/levels in the vertical
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: p !< Interface pressures [R L2 T-2 ~> Pa].
  real,                                 intent(in)  :: dP_dH !< A factor that converts from thickness to pressure
                                                             !! times other dimensional conversion factors
                                                             !! [R L2 T-2 H-1 ~> Pa m2 kg-1].
  real,                                 intent(in)  :: dp_neglect !< A thickness that is so small it is usually lost
                                                                  !! in roundoff and can be neglected
                                                                  !! [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G)),     intent(inout) :: C_htot !< dP_dH divided by the total ocean pressure
                                                                !! [H-1 ~> m2 kg-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: pbce !< The baroclinic pressure anomaly in each
                                                          !! layer due to free surface height anomalies
                                                          !! [L2 H-1 T-2 ~> m4 kg-1 s-2].
  ! Local variables
  real :: T_int(SZI_(G))     ! Interface temperature [C ~> degC]
  real :: S_int(SZI_(G))     ! Interface salinity [S ~> ppt]
  real :: dR_dT(SZI_(G))     ! Partial derivative of density with temperature [R C-1 ~> kg m-3 degC-1]
  real :: dR_dS(SZI_(G))     ! Partial derivative of density with salinity [R S-1 ~> kg m-3 ppt-1].
  real :: rho_in_situ(SZI_(G)) ! In-situ density at an interface [R ~> kg m-3].
  integer, dimension(2) :: EOSdom ! The computational domain for the equation of state
  integer :: i, j, k

  EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)

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
end subroutine Set_pbce_nonBouss_EOS

!> Finds the change in geopotential across each layer for PressureForce_Mont_nonBouss with
!! int_specific_vol_dp, using the grid, unit scaling and thermodynamic variables in the
!! PF_Mont_host_ctx that host points to.
subroutine Mont_nonBouss_dz_geo_host(host, nz, p, dz_geo)
  type(c_ptr),     value, intent(in)    :: host !< Opaque pointer to the PF_Mont_host_ctx that holds
                                                !! the Fortran objects used by the host helpers
  integer,                intent(in)    :: nz !< The number of layers/levels in the vertical
  real, dimension(:,:,:), intent(in)    :: p !< Interface pressure [R L2 T-2 ~> Pa].
  real, dimension(:,:,:), intent(inout) :: dz_geo !< The change in geopotential across a layer [L2 T-2 ~> m2 s-2].
  ! Local variables
  type(PF_Mont_host_ctx), pointer :: ctx ! The objects behind host
  type(ocean_grid_type),  pointer :: G   ! The ocean's grid structure.
  type(unit_scale_type),  pointer :: US  ! A dimensional unit scaling type
  type(thermo_var_ptrs),  pointer :: tv  ! Thermodynamic variables.
  integer :: k

  call c_f_pointer(host, ctx)
  G => ctx%G ; US => ctx%US ; tv => ctx%tv

  !$OMP parallel do default(shared)
  do k=1,nz
    call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p(:,:,k), p(:,:,k+1), &
                             0.0, G%HI, tv%eqn_of_state, US, dz_geo(:,:,k), halo_size=1)
  enddo

end subroutine Mont_nonBouss_dz_geo_host

!> Calls Mont_nonBouss_alpha_star with the grid, vertical grid, thermodynamic variables and temporary
!! T and S storage in the PF_Mont_host_ctx that host points to, then points the thermodynamic
!! variables that the Set_pbce host helpers use at the temporary ones.
subroutine Mont_nonBouss_alpha_star_host(host, Isq, Ieq, Jsq, Jeq, nz, alpha_star)
  type(c_ptr),     value, intent(in)    :: host !< Opaque pointer to the PF_Mont_host_ctx that holds
                                                !! the Fortran objects used by the host helpers
  integer,                intent(in)    :: Isq !< The start i-index of cell vertices within the computational domain
  integer,                intent(in)    :: Ieq !< The end i-index of cell vertices within the computational domain
  integer,                intent(in)    :: Jsq !< The start j-index of cell vertices within the computational domain
  integer,                intent(in)    :: Jeq !< The end j-index of cell vertices within the computational domain
  integer,                intent(in)    :: nz !< The number of layers/levels in the vertical
  real, dimension(:,:,:), intent(inout) :: alpha_star !< Compression adjusted specific volume [R-1 ~> m3 kg-1].
  ! Local variables
  type(PF_Mont_host_ctx),  pointer :: ctx ! The objects behind host
  type(ocean_grid_type),   pointer :: G   ! The ocean's grid structure.
  type(verticalGrid_type), pointer :: GV  ! The ocean's vertical grid structure.

  call c_f_pointer(host, ctx)
  G => ctx%G ; GV => ctx%GV
  if (.not.allocated(ctx%pbce%T_tmp)) allocate(ctx%pbce%T_tmp(SZI_(G),SZJ_(G),SZK_(GV)))
  if (.not.allocated(ctx%pbce%S_tmp)) allocate(ctx%pbce%S_tmp(SZI_(G),SZJ_(G),SZK_(GV)))
  call Mont_nonBouss_alpha_star(G, GV, ctx%tv, ctx%pbce%tv_tmp, ctx%pbce%T_tmp, ctx%pbce%S_tmp, &
                                Isq, Ieq, Jsq, Jeq, nz, alpha_star)
  ctx%pbce%tv => ctx%pbce%tv_tmp

end subroutine Mont_nonBouss_alpha_star_host

!> Calculates the in-situ specific volumes (alpha_star) for PressureForce_Mont_nonBouss from the
!! equation of state, setting up tv_tmp with the temperatures and salinities it uses.
subroutine Mont_nonBouss_alpha_star(G, GV, tv, tv_tmp, T_tmp, S_tmp, Isq, Ieq, Jsq, Jeq, nz, alpha_star)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure.
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables.
  type(thermo_var_ptrs),                      intent(inout) :: tv_tmp !< A structure of temporary T & S.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target, intent(inout) :: &
    T_tmp       !< Temporary array of temperatures where layers that are lighter
                !! than the mixed layer have the mixed layer's properties [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target, intent(inout) :: &
    S_tmp       !< Temporary array of salinities where layers that are lighter
                !! than the mixed layer have the mixed layer's properties [S ~> ppt].
  integer,                                    intent(in)  :: Isq !< The start i-index of cell vertices within the
                                                                 !! computational domain
  integer,                                    intent(in)  :: Ieq !< The end i-index of cell vertices within the
                                                                 !! computational domain
  integer,                                    intent(in)  :: Jsq !< The start j-index of cell vertices within the
                                                                 !! computational domain
  integer,                                    intent(in)  :: Jeq !< The end j-index of cell vertices within the
                                                                 !! computational domain
  integer,                                    intent(in)  :: nz !< The number of layers/levels in the vertical
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: &
    alpha_star  !< Compression adjusted specific volume [R-1 ~> m3 kg-1].
  ! Local variables
  real, dimension(SZI_(G)) :: Rho_cv_BL  ! The coordinate potential density in the
                  ! deepest variable density near-surface layer [R ~> kg m-3].
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density [R L2 T-2 ~> Pa] (usually 2e7 Pa = 2000 dbar).
  real :: rho_in_situ(SZI_(G)) !In-situ density of a layer [R ~> kg m-3].
  integer, dimension(2) :: EOSdom ! The computational domain for the equation of state
  integer :: i, j, k, nkmb

  nkmb=GV%nk_rho_varies
  EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)

  !   Calculate in-situ specific volumes (alpha_star).

  !   With a bulk mixed layer, replace the T & S of any layers that are
  ! lighter than the buffer layer with the properties of the buffer
  ! layer.  These layers will be massless anyway, and it avoids any
  ! formal calculations with hydrostatically unstable profiles.
  if (nkmb>0) then
    tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
    tv_tmp%eqn_of_state => tv%eqn_of_state
    do i=Isq,Ieq+1 ; p_ref(i) = tv%P_Ref ; enddo
    !$OMP parallel do default(shared) private(Rho_cv_BL)
    do j=Jsq,Jeq+1
      do k=1,nkmb ; do i=Isq,Ieq+1
        tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
      enddo ; enddo
      call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, Rho_cv_BL(:), &
                             tv%eqn_of_state, EOSdom)
      do k=nkmb+1,nz ; do i=Isq,Ieq+1
        if (GV%Rlay(k) < Rho_cv_BL(i)) then
          tv_tmp%T(i,j,k) = tv%T(i,j,nkmb) ; tv_tmp%S(i,j,k) = tv%S(i,j,nkmb)
        else
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        endif
      enddo ; enddo
    enddo
  else
    tv_tmp%T => tv%T ; tv_tmp%S => tv%S
    tv_tmp%eqn_of_state => tv%eqn_of_state
    do i=Isq,Ieq+1 ; p_ref(i) = 0 ; enddo
  endif
  !$OMP parallel do default(shared) private(rho_in_situ)
  do k=1,nz ; do j=Jsq,Jeq+1
    call calculate_density(tv_tmp%T(:,j,k), tv_tmp%S(:,j,k), p_ref, rho_in_situ, &
                           tv%eqn_of_state, EOSdom)
    do i=Isq,Ieq+1 ; alpha_star(i,j,k) = 1.0 / rho_in_situ(i) ; enddo
  enddo ; enddo
end subroutine Mont_nonBouss_alpha_star

!> Calls Mont_Bouss_rho_star with the grid, vertical grid, thermodynamic variables and temporary
!! T and S storage in the PF_Mont_host_ctx that host points to, then points the thermodynamic
!! variables that the Set_pbce host helpers use at the temporary ones.
subroutine Mont_Bouss_rho_star_host(host, Isq, Ieq, Jsq, Jeq, nz, rho_star)
  type(c_ptr),     value, intent(in)    :: host !< Opaque pointer to the PF_Mont_host_ctx that holds
                                                !! the Fortran objects used by the host helpers
  integer,                intent(in)    :: Isq !< The start i-index of cell vertices within the computational domain
  integer,                intent(in)    :: Ieq !< The end i-index of cell vertices within the computational domain
  integer,                intent(in)    :: Jsq !< The start j-index of cell vertices within the computational domain
  integer,                intent(in)    :: Jeq !< The end j-index of cell vertices within the computational domain
  integer,                intent(in)    :: nz !< The number of layers/levels in the vertical
  real, dimension(:,:,:), intent(inout) :: rho_star !< In-situ density divided by the derivative with depth of the
                                                    !! corrected e times (G_Earth/Rho0) [L2 Z-1 T-2 ~> m s-2].
  ! Local variables
  type(PF_Mont_host_ctx),  pointer :: ctx ! The objects behind host
  type(ocean_grid_type),   pointer :: G   ! The ocean's grid structure.
  type(verticalGrid_type), pointer :: GV  ! The ocean's vertical grid structure.

  call c_f_pointer(host, ctx)
  G => ctx%G ; GV => ctx%GV
  if (.not.allocated(ctx%pbce%T_tmp)) allocate(ctx%pbce%T_tmp(SZI_(G),SZJ_(G),SZK_(GV)))
  if (.not.allocated(ctx%pbce%S_tmp)) allocate(ctx%pbce%S_tmp(SZI_(G),SZJ_(G),SZK_(GV)))
  call Mont_Bouss_rho_star(G, GV, ctx%tv, ctx%pbce%tv_tmp, ctx%pbce%T_tmp, ctx%pbce%S_tmp, &
                           Isq, Ieq, Jsq, Jeq, nz, rho_star)
  ctx%pbce%tv => ctx%pbce%tv_tmp

end subroutine Mont_Bouss_rho_star_host

!> Calculates the in-situ densities (rho_star) for PressureForce_Mont_Bouss from the equation of
!! state, setting up tv_tmp with the temperatures and salinities it uses.
subroutine Mont_Bouss_rho_star(G, GV, tv, tv_tmp, T_tmp, S_tmp, Isq, Ieq, Jsq, Jeq, nz, rho_star)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure.
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables.
  type(thermo_var_ptrs),                      intent(inout) :: tv_tmp !< A structure of temporary T & S.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target, intent(inout) :: &
    T_tmp       !< Temporary array of temperatures where layers that are lighter
                !! than the mixed layer have the mixed layer's properties [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target, intent(inout) :: &
    S_tmp       !< Temporary array of salinities where layers that are lighter
                !! than the mixed layer have the mixed layer's properties [S ~> ppt].
  integer,                                    intent(in)  :: Isq !< The start i-index of cell vertices within the
                                                                 !! computational domain
  integer,                                    intent(in)  :: Ieq !< The end i-index of cell vertices within the
                                                                 !! computational domain
  integer,                                    intent(in)  :: Jsq !< The start j-index of cell vertices within the
                                                                 !! computational domain
  integer,                                    intent(in)  :: Jeq !< The end j-index of cell vertices within the
                                                                 !! computational domain
  integer,                                    intent(in)  :: nz !< The number of layers/levels in the vertical
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: &
    rho_star    !< In-situ density divided by the derivative with depth of the
                !! corrected e times (G_Earth/Rho0) [L2 Z-1 T-2 ~> m s-2].
  ! Local variables
  real :: Rho_cv_BL(SZI_(G)) !   The coordinate potential density in
                ! the deepest variable density near-surface layer [R ~> kg m-3].
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density [R L2 T-2 ~> Pa] (usually 2e7 Pa = 2000 dbar).
  real :: G_Rho0             ! G_Earth / Rho0 [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1].
  integer, dimension(2) :: EOSdom ! The computational domain for the equation of state
  integer :: i, j, k, nkmb

  nkmb=GV%nk_rho_varies
  EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)
  G_Rho0 = GV%g_Earth / GV%Rho0

!   Calculate in-situ densities (rho_star).

! With a bulk mixed layer, replace the T & S of any layers that are
! lighter than the buffer layer with the properties of the buffer
! layer.  These layers will be massless anyway, and it avoids any
! formal calculations with hydrostatically unstable profiles.

  if (nkmb>0) then
    tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
    tv_tmp%eqn_of_state => tv%eqn_of_state

    do i=Isq,Ieq+1 ; p_ref(i) = tv%P_Ref ; enddo
    !$OMP parallel do default(shared) private(Rho_cv_BL)
    do j=Jsq,Jeq+1
      do k=1,nkmb ; do i=Isq,Ieq+1
        tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
      enddo ; enddo
      call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, Rho_cv_BL(:), &
                             tv%eqn_of_state, EOSdom)

      do k=nkmb+1,nz ; do i=Isq,Ieq+1
        if (GV%Rlay(k) < Rho_cv_BL(i)) then
          tv_tmp%T(i,j,k) = tv%T(i,j,nkmb) ; tv_tmp%S(i,j,k) = tv%S(i,j,nkmb)
        else
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        endif
      enddo ; enddo
    enddo
  else
    tv_tmp%T => tv%T ; tv_tmp%S => tv%S
    tv_tmp%eqn_of_state => tv%eqn_of_state
    do i=Isq,Ieq+1 ; p_ref(i) = 0.0 ; enddo
  endif

  ! This no longer includes any pressure dependency, since this routine
  ! will come down with a fatal error if there is any compressibility.
  !$OMP parallel do default(shared)
  do k=1,nz ; do j=Jsq,Jeq+1
    call calculate_density(tv_tmp%T(:,j,k), tv_tmp%S(:,j,k), p_ref, rho_star(:,j,k), &
                           tv%eqn_of_state, EOSdom)
    do i=Isq,Ieq+1 ; rho_star(i,j,k) = G_Rho0*rho_star(i,j,k) ; enddo
  enddo ; enddo
end subroutine Mont_Bouss_rho_star

!> Calls calc_SAL with the grid, unit scaling and self-attraction and loading control structure in
!! the PF_Mont_host_ctx that host points to.
subroutine Mont_calc_SAL_host(host, eta, e_sal)
  type(c_ptr),     value, intent(in)    :: host !< Opaque pointer to the PF_Mont_host_ctx that holds
                                                !! the Fortran objects used by the host helpers
  real, dimension(:,:),   intent(in)    :: eta !< The sea surface height anomaly from a time-mean geoid or total
                                               !! bottom pressure [Z ~> m] or [R L2 T-2 ~> Pa].
  real, dimension(:,:),   intent(inout) :: e_sal !< The bottom geopotential anomaly due to self-attraction and loading
                                                 !! [Z ~> m].
  ! Local variables
  type(PF_Mont_host_ctx),      pointer :: ctx ! The objects behind host
  type(ocean_grid_type),       pointer :: G   ! The ocean's grid structure.
  type(unit_scale_type),       pointer :: US  ! A dimensional unit scaling type
  type(PressureForce_Mont_CS), pointer :: CS  ! The Montgomery potential PGF control structure

  call c_f_pointer(host, ctx)
  G => ctx%G ; US => ctx%US ; CS => ctx%CS

  call calc_SAL(eta, e_sal, G, CS%SAL_CSp, tmp_scale=US%Z_to_m)

end subroutine Mont_calc_SAL_host

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
  call get_param(param_file, mdl, "RHO_0", CS%PF_Mont_opts_CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINESQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%R_to_kg_m3)
  call get_param(param_file, mdl, "TIDES", CS%PF_Mont_opts_CS%tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  call get_param(param_file, mdl, "CALCULATE_SAL", CS%PF_Mont_opts_CS%calculate_SAL, &
                 "If true, calculate self-attraction and loading.", default=CS%PF_Mont_opts_CS%tides)
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

  if (CS%PF_Mont_opts_CS%calculate_SAL) then
    CS%id_e_sal = register_diag_field('ocean_model', 'e_sal', diag%axesT1, Time, &
        'Self-attraction and loading height anomaly', 'meter', conversion=US%Z_to_m)
  endif
  if (CS%PF_Mont_opts_CS%tides) then
    CS%id_e_tide = register_diag_field('ocean_model', 'e_tidal', diag%axesT1, Time, &
        'Tidal Forcing Astronomical and SAL Height Anomaly', 'meter', conversion=US%Z_to_m)
    CS%id_e_tide_eq  = register_diag_field('ocean_model', 'e_tide_eq', diag%axesT1, Time, &
        'Equilibrium tides height anomaly', 'meter', conversion=US%Z_to_m)
    CS%id_e_tide_sal = register_diag_field('ocean_model', 'e_tide_sal', diag%axesT1, Time, &
        'Read-in tidal self-attraction and loading height anomaly', 'meter', conversion=US%Z_to_m)
  endif

  CS%PF_Mont_opts_CS%GFS_scale = 1.0
  if (GV%g_prime(1) /= GV%g_Earth) CS%PF_Mont_opts_CS%GFS_scale = GV%g_prime(1) / GV%g_Earth

  call PF_vgrid_consts_build(CS%PF_vgrid_consts_CS, GV)

  call log_param(param_file, mdl, "GFS / G_EARTH", CS%PF_Mont_opts_CS%GFS_scale, units="nondim")

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
