!=============================================================================
! MODULE: Module Name
!
!> @author
!> Felix Dietzsch
!
! DESCRIPTION
!> Contains definition of error messages
!>
! REVISION HISTORY
! 25 07 2011 - Initial Version
!=============================================================================
MODULE ONERROR
IMPLICIT NONE
SAVE
CONTAINS
    !=============================================================================
    !> @author Felix Dietzsch
    !
    ! DESCRIPTION:
    !> Prints allocation error message
    !
    ! REVISION HISTORY:
    ! 25 07 2011 Initial Version
    !
    !> @param[in] IERR Error descriptor
    !=============================================================================
    SUBROUTINE ALLOCATION_ERROR(IERR)
        INTEGER :: IERR
        PRINT *,'An error occured during allocation: ',IERR
        CALL EXIT(1)
    END SUBROUTINE
    !
    !=============================================================================
    !> @author Felix Dietzsch
    !
    ! DESCRIPTION:
    !> Prints deallocation error message
    !
    ! REVISION HISTORY:
    ! 25 07 2011 Initial Version
    !
    !> @param[in] IERR Error descriptor
    !> @param[in] UNIT_VALUE Error descriptor\n
    !> 1 -- COORD\n
    !> 2 -- INDX\n
    !> 3 -- INDY\n
    !> 4 -- INDZ\n
    !> 5 -- UUX\n
    !> 6 -- UUY\n
    !> 7 -- UUZ\n
    !=============================================================================
    SUBROUTINE DEALLOCATION_ERROR(IERR,UNIT_VALUE)
        INTEGER :: IERR,UNIT_VALUE
        PRINT *,'An error occured during deallocation of: ',UNIT_VALUE,&
           'with error core',IERR
        CALL EXIT(1)
    END SUBROUTINE
END MODULE ONERROR
