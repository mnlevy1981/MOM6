module marbl_interface

    use MOM_error_handler,            only : MOM_error, FATAL
    use marbl_logging,                only : marbl_log_type
    use marbl_interface_public_types, only : marbl_forcing_fields_type
    use marbl_interface_public_types, only : marbl_tracer_metadata_type
    use marbl_interface_public_types, only : marbl_saved_state_type
    use marbl_interface_public_types, only : marbl_diagnostics_type
    use marbl_interface_public_types, only : marbl_domain_type
    implicit none
    private ! Only want marbl_interface_class to be public, not supporting functions

    type, public :: marbl_interface_class
        type(marbl_log_type) :: StatusLog
        type(marbl_forcing_fields_type), allocatable :: surface_flux_forcings(:)
        type(marbl_forcing_fields_type), allocatable :: interior_tendency_forcings(:)
        type(marbl_tracer_metadata_type), allocatable :: tracer_metadata(:)
        type(marbl_domain_type) :: domain
        type(marbl_saved_state_type) :: surface_flux_saved_state
        type(marbl_saved_state_type) :: interior_tendency_saved_state
        type(marbl_diagnostics_type) :: interior_tendency_diags
        type(marbl_diagnostics_type) :: surface_flux_diags
        real, allocatable :: tracers(:,:)
        real, allocatable :: tracers_at_surface(:,:)
        real, allocatable :: surface_fluxes(:,:)
        real, allocatable :: interior_tendencies(:,:)
       contains
        procedure, public  :: init
        procedure, public  :: surface_flux_compute
        procedure, public  :: interior_tendency_compute
        procedure, public  :: shutdown
        procedure, public :: put_setting
    end type marbl_interface_class

    character(len=*), parameter :: error_msg = "MOM6 built the MARBL stubs rather than the full library"

contains

    subroutine put_setting(self, str_in)
        class(marbl_interface_class), intent(in) :: self
        character(len=*),             intent(in) :: str_in

        call MOM_error(FATAL, error_msg)
    end subroutine put_setting

    subroutine init(self,                  &
        gcm_num_levels,                    &
        gcm_num_PAR_subcols,               &
        gcm_num_elements_surface_flux,     &
        gcm_delta_z,                       &
        gcm_zw,                            &
        gcm_zt,                            &
        lgcm_has_global_ops)

        class(marbl_interface_class), intent(inout) :: self
        integer,                      intent(in)    :: gcm_num_levels
        integer,                      intent(in)    :: gcm_num_PAR_subcols
        integer,                      intent(in)    :: gcm_num_elements_surface_flux
        real,                         intent(in)    :: gcm_delta_z(gcm_num_levels)
        real,                         intent(in)    :: gcm_zw(gcm_num_levels)
        real,                         intent(in)    :: gcm_zt(gcm_num_levels)
        logical,                      intent(in)    :: lgcm_has_global_ops

        call MOM_error(FATAL, error_msg)
    end subroutine init

    subroutine surface_flux_compute(self)

        class(marbl_interface_class), intent(inout) :: self

        call MOM_error(FATAL, error_msg)

    end subroutine surface_flux_compute

    subroutine interior_tendency_compute(self)

        class(marbl_interface_class), intent(inout) :: self

        call MOM_error(FATAL, error_msg)

    end subroutine interior_tendency_compute

    subroutine shutdown(self)

        class(marbl_interface_class), intent(inout) :: self

        call MOM_error(FATAL, error_msg)

    end subroutine shutdown

end module marbl_interface