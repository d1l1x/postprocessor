MODULE PREMIXED_CLASS
USE NRTYPE
    INTERFACE COMP_PROGRESS 
      !MODULE PROCEDURE COMP_PROGRESS
        SUBROUTINE COMP_PROGRESS(TEMPER,PROG_VAR,SAVEVAR,SPEC,SPECUB,SPECB,limit,c_limit_index)
        USE NRTYPE
        USE IO_CLASS
        REAL(DP) :: TEMP_MIN,TEMP_MAX
        REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: PROG_VAR
        REAL(DP),DIMENSION(:,:,:),INTENT(IN) :: TEMPER
        REAL(DP),DIMENSION(:,:,:),INTENT(IN),OPTIONAL :: SPEC
        LOGICAL,INTENT(IN),OPTIONAL :: SAVEVAR
        REAL(DP),INTENT(IN),OPTIONAL :: SPECUB
        REAL(DP),INTENT(IN),OPTIONAL :: SPECB
        REAL(SP),INTENT(INOUT),OPTIONAL :: limit 
        LOGICAL,DIMENSION(:,:,:),ALLOCATABLE :: temp
        LOGICAL,DIMENSION(:,:,:),OPTIONAL,INTENT(INOUT):: c_limit_index
        END SUBROUTINE COMP_PROGRESS
    END INTERFACE COMP_PROGRESS 
    INTERFACE GET_VAL_UNBURNT
        SUBROUTINE GET_VAL_UNBURNT(var,c_limit_index,output)
            USE NRTYPE
            IMPLICIT NONE
            INTEGER(SP) :: I,J,K
            REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: var
            LOGICAL,DIMENSION(:,:,:),INTENT(INOUT) :: c_limit_index
            REAL(DP),DIMENSION(:),INTENT(INOUT) :: output
        END SUBROUTINE GET_VAL_UNBURNT
    END INTERFACE GET_VAL_UNBURNT
END MODULE PREMIXED_CLASS
! subroutine for computing the progress variable
SUBROUTINE COMP_PROGRESS(TEMPER,PROG_VAR,SAVEVAR,SPEC,SPECUB,SPECB,&
                         limit,c_limit_index)
    USE NRTYPE
    USE IO_CLASS
    REAL(DP) :: TEMP_MIN,TEMP_MAX
    REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: PROG_VAR
    REAL(DP),DIMENSION(:,:,:),INTENT(IN) :: TEMPER
    REAL(DP),DIMENSION(:,:,:),INTENT(IN),OPTIONAL :: SPEC
    LOGICAL,INTENT(IN),OPTIONAL :: SAVEVAR
    REAL(DP),INTENT(IN),OPTIONAL :: SPECUB
    REAL(DP),INTENT(IN),OPTIONAL :: SPECB
    REAL(SP),INTENT(INOUT),OPTIONAL :: limit 
    LOGICAL,DIMENSION(:,:,:),OPTIONAL,INTENT(INOUT):: c_limit_index

    IF (PRESENT(SPEC))THEN
      ! for the moment this is only valid for CO2
        PROG_VAR(:,:,:) = (SPEC(:,:,:) - SPECUB)/(SPECB - SPECUB)
    ELSE
        TEMP_MIN = MINVAL(TEMPER)
        TEMP_MAX = MAXVAL(TEMPER)
        PROG_VAR(:,:,:) = (TEMPER(:,:,:)-TEMP_MIN)/(TEMP_MAX-TEMP_MIN)
    ENDIF
    IF (SAVEVAR.AND. PRESENT(SPEC)) THEN
        CALL WRITE_PROG('./OUTPUT/PROG_VAR_SPEC',PROG_VAR)
        CALL WRITE_PROG('./OUTPUT/SPEC',SPEC)
    ELSE
        CALL WRITE_PROG('./OUTPUT/PROG_VAR',PROG_VAR)
        CALL WRITE_PROG('./OUTPUT/TEMPER',TEMPER)
    ENDIF
    IF (PRESENT(limit).AND.PRESENT(c_limit_index)) THEN
        c_limit_index(:,:,:) = PROG_VAR(:,:,:) < limit  
    END IF
END SUBROUTINE COMP_PROGRESS
!=============================================================================
!> @author Felix Dietzsch
!
! DESCRIPTION:
!> Extracts the velocity of the unburnt mixture.
!
! REVISION HISTORY:
! 19 01 2012 Initial Version
!
!=============================================================================
SUBROUTINE GET_VAL_UNBURNT(var,c_limit_index,output)
    USE NRTYPE
    IMPLICIT NONE
    INTEGER(SP) :: i
    REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: var
    LOGICAL,DIMENSION(:,:,:),INTENT(INOUT) :: c_limit_index
    REAL(DP),DIMENSION(:),INTENT(INOUT) :: output
    LOGICAL,DIMENSION(:,:,:),ALLOCATABLE :: mask
    
    ALLOCATE(mask(size(var,1),size(var,2),size(var,3)))
    mask = (c_limit_index.EQV..TRUE.)
    i = COUNT(mask)
    output = PACK(var,mask)
    PRINT*,"count:",i
END SUBROUTINE GET_VAL_UNBURNT
