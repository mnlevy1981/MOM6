module marbl_logging

    implicit none
    private

    type, public :: marbl_status_log_entry_type
        integer :: ElementInd
        logical :: lonly_master_writes
        character(len=0) :: LogMessage
        type(marbl_status_log_entry_type), pointer :: next
    end type marbl_status_log_entry_type

    type, public :: marbl_log_type
        logical, public  :: labort_marbl
        type(marbl_status_log_entry_type), pointer :: FullLog
    contains
        procedure, public :: log_error_trace
        procedure, public :: erase
    end type marbl_log_type

contains

    subroutine log_error_trace(self, RoutineName, CodeLoc, ElemInd)
        class(marbl_log_type), intent(inout) :: self
        character(len=*),      intent(in)    :: RoutineName, CodeLoc
        integer, optional,     intent(in)    :: ElemInd
    end subroutine log_error_trace

    subroutine erase(self)
        class(marbl_log_type), intent(inout) :: self
    end subroutine erase

end module marbl_logging