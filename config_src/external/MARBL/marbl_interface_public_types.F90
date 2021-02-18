module marbl_interface_public_types

    implicit none
    private ! Only want a few types to be public

    type :: marbl_single_diagnostic_type
        character(len=0)                  :: long_name
        character(len=0)                  :: short_name
        character(len=0)                  :: units
        character(len=0)                  :: vertical_grid
        logical                           :: compute_now
        logical                           :: ltruncated_vertical_extent
        integer                           :: ref_depth
        real, allocatable, dimension(:)   :: field_2d
        real, allocatable, dimension(:,:) :: field_3d
    end type marbl_single_diagnostic_type

    type, public :: marbl_diagnostics_type
        type(marbl_single_diagnostic_type), dimension(:), pointer :: diags
    end type marbl_diagnostics_type

    type :: marbl_single_saved_state_type
        integer :: rank
        character(len=0) :: short_name
        character(len=0) :: units
        character(len=0) :: vertical_grid
        real, allocatable :: field_2d(:)
        real, allocatable :: field_3d(:,:)
    end type marbl_single_saved_state_type

    type, public :: marbl_saved_state_type
        integer :: saved_state_cnt
        type(marbl_single_saved_state_type), dimension(:), pointer :: state => NULL()
    end type marbl_saved_state_type

    type :: marbl_forcing_fields_metadata_type
        character(len=0) :: varname
        ! character(len=0) :: field_units
        ! integer          :: rank
        ! integer,  allocatable :: extent(:)
    end type marbl_forcing_fields_metadata_type

    type, public :: marbl_forcing_fields_type
        type(marbl_forcing_fields_metadata_type) :: metadata
        real, pointer :: field_0d(:)   => NULL()
        real, pointer :: field_1d(:,:) => NULL()
    end type marbl_forcing_fields_type

    type, public :: marbl_tracer_metadata_type
        character(len=0) :: short_name
        character(len=0) :: long_name
        character(len=0) :: units
        ! character(len=0) :: tend_units
        ! character(len=0) :: flux_units
        ! logical          :: lfull_depth_tavg
        ! character(len=0) :: tracer_module_name
    end type marbl_tracer_metadata_type

type, public :: marbl_domain_type
    integer :: kmt
    real, allocatable :: zt(:)
    real, allocatable :: zw(:)
    real, allocatable :: delta_z(:)
end type marbl_domain_type

end module marbl_interface_public_types