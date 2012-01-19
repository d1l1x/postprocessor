MODULE TIMER_CLASS
USE NRTYPE
IMPLICIT NONE
! Type definition
TYPE,PUBLIC :: TIMER ! This will be the name we instantiate
    PRIVATE
    REAL(DP) :: saved_time
    CONTAINS
    PROCEDURE,PUBLIC :: START_TIMER => start_timer_sub
    PROCEDURE,PUBLIC :: ELAPSED_TIME => elapsed_time_fn
END TYPE TIMER
!
PRIVATE :: start_timer_sub, elapsed_time_fn
!
CONTAINS
    SUBROUTINE start_timer_sub(this)
    IMPLICIT NONE
    CLASS(timer) :: this
    INTEGER,DIMENSION(8) :: value
    CALL date_and_time(VALUES=value)
    this%saved_time = 86400.D0 * value(3) + 3600.D0 * value(5) + 60.D0 * value(6) + value(7) + 0.001D0 * value(8)
    END SUBROUTINE start_timer_sub
    !
    REAL(DP) FUNCTION elapsed_time_fn(this)
    IMPLICIT NONE
    CLASS(timer) :: this
    INTEGER,DIMENSION(8) :: value
    REAL(DP) :: current_time
    CALL date_and_time(VALUES=value)
    current_time = 86400.D0 * value(3) + 3600.D0 * value(5) + 60.D0 * value(6) + value(7) + 0.001D0 * value(8)
    ! Get elapsed time
    elapsed_time_fn = current_time - this%saved_time
    END FUNCTION elapsed_time_fn
END MODULE TIMER_CLASS
