MODULE FIELD_CLASS
USE NRTYPE
IMPLICIT NONE
TYPE,PUBLIC :: field
    PRIVATE
    INTEGER(SP) :: dim1
    INTEGER(SP) :: dim2
    INTEGER(SP) :: dim3
    REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: array
    INTEGER(SP) :: ierr
    CONTAINS
    PROCEDURE,PUBLIC :: create_2d => alloc_array2d_f
    PROCEDURE,PUBLIC :: create_3d => alloc_array3d_f
END TYPE
PRIVATE :: alloc_array2d_f,alloc_array3d_f 
CONTAINS
    REAL FUNCTION alloc_array2d_f(this,dimen)
    USE NRTYPE
    IMPLICIT NONE
    CLASS(field) :: this
    REAL(DP),DIMENSION(3) :: dimen
    !REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: array
    !ALLOCATE(this%array(dimen(1),dimen(2),1),STAT=this%ierr)
    !IF (this%ierr.NE.0) CALL ALLOCATION_ERROR(this%ierr)
    alloc_array2d_f = 1
    END FUNCTION alloc_array2d_f
    !
    FUNCTION alloc_array3d_f(this,dimen)
    USE NRTYPE
    IMPLICIT NONE
    CLASS(field) :: this
    INTEGER(SP),DIMENSION(3) :: dimen
    REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: array
    REAL(DP),DIMENSION(dimen(1),dimen(2),dimen(3)) :: alloc_array3d_f
    this%dim1=dimen(1)
    this%dim2=dimen(2)
    this%dim3=dimen(3)
    ALLOCATE(this%array(this%dim1,this%dim2,this%dim3),STAT=this%ierr)
    IF (this%ierr.NE.0) THEN
        PRINT *,'An error occured during allocation: ',this%ierr
        CALL EXIT(1)
    END IF
    alloc_array3d_f = this%array
    END FUNCTION alloc_array3d_f
END MODULE FIELD_CLASS
