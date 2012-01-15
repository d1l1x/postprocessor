!=============================================================================
! MODULE: Module Name
!
!> @author
!> Felix Dietzsch
!
! DESCRIPTION
!> Controls the input and output for the underlying postprocessor.
!>
! REVISION HISTORY
! 25 07 2011 - Initial Version
! 05 12 2011 - Added reading and writing routines
!=============================================================================
MODULE IO_CLASS
USE NRTYPE
    INTERFACE INDEXING
        SUBROUTINE INDEXING(COORD,NUMBER_OF_NODES,INDX,NX,INDY,NY,INDZ,NZ)
        USE INIT
        USE NRTYPE
        REAL(DP) :: XMIN,XMAX,YMIN,YMAX,LX,LY
        INTEGER(SP) :: NUMBER_OF_NODES
        INTEGER(SP),INTENT(INOUT) :: NX,NY
        INTEGER(SP),INTENT(INOUT),OPTIONAL :: NZ
        INTEGER(SP) :: IERR
        INTEGER(SP) :: N
        INTEGER(SP) :: NXM,NYM
        INTEGER(SP),DIMENSION(:),INTENT(INOUT) :: INDX,INDY
        INTEGER(SP),DIMENSION(:),INTENT(INOUT),OPTIONAL :: INDZ
        REAL(DP),DIMENSION(:,:) :: COORD
        END SUBROUTINE INDEXING
    END INTERFACE INDEXING
    INTERFACE READCOORD
        SUBROUTINE READCOORD2D(FILENAME,COORD)
        USE INIT
        USE NRTYPE
        USE ONERROR
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        INTEGER(SP) :: NUMBER_OF_NODES
        INTEGER(SP) :: IERR,STATUS_READ
        REAL(DP),DIMENSION(:,:),INTENT(INOUT) :: COORD
        END SUBROUTINE READCOORD2D
    END INTERFACE READCOORD
    INTERFACE READ_DIM
      SUBROUTINE READ_DIM(FILEX,FILEY,FILEZ,DIMEN)
        USE NRTYPE
        INTEGER(SP),DIMENSION(:),INTENT(INOUT) :: DIMEN
        CHARACTER(LEN=*),INTENT(IN) :: FILEX
        CHARACTER(LEN=*),INTENT(IN) :: FILEY
        CHARACTER(LEN=*),INTENT(IN) :: FILEZ
        LOGICAL :: FILE_EXISTS
        INTEGER(SP) :: NPX,NPY,NPZ
        INTEGER(SP) :: I
        END SUBROUTINE READ_DIM
    END INTERFACE READ_DIM
    INTERFACE READ_GRID
      SUBROUTINE READ_GRID(FILEX,FILEY,FILEZ,GRIDX,GRIDY,GRIDZ)
        USE NRTYPE
        REAL(DP),DIMENSION(:),INTENT(INOUT) :: GRIDX
        REAL(DP),DIMENSION(:),INTENT(INOUT) :: GRIDY
        REAL(DP),DIMENSION(:),INTENT(INOUT) :: GRIDZ
        CHARACTER(LEN=*),INTENT(IN) :: FILEX
        CHARACTER(LEN=*),INTENT(IN) :: FILEY
        CHARACTER(LEN=*),INTENT(IN) :: FILEZ
        LOGICAL :: FILE_EXISTS
        INTEGER(SP) :: NPX,NPY,NPZ
        INTEGER(SP) :: I
      END SUBROUTINE READ_GRID
    END INTERFACE READ_GRID
    INTERFACE READVALUE
      SUBROUTINE READVALUE(FILENAME,VAR)
         USE INIT
         USE NRTYPE
         USE ONERROR
         INTEGER :: I,J,K
         INTEGER :: NPROC_Z,NPOINTS_Z
         REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: VAR
         CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE READVALUE 
    END INTERFACE READVALUE
    INTERFACE WRITE_VALUE
      SUBROUTINE WRITE_VALUE(FILENAME,value)
        USE NRTYPE
        INTEGER :: I,J,K
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        REAL(DP),DIMENSION(:,:,:),INTENT(IN) :: value 
        END SUBROUTINE WRITE_VALUE
    END INTERFACE WRITE_VALUE
    INTERFACE WRITE_PROG
      SUBROUTINE WRITE_PROG(FILENAME,PROG_VAR)
        USE NRTYPE
        INTEGER :: I,J,K
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        REAL(DP),DIMENSION(:,:,:),INTENT(IN) :: PROG_VAR
      END SUBROUTINE WRITE_PROG
    END INTERFACE WRITE_PROG
    INTERFACE READVEL
        SUBROUTINE READVEL(FILENAME,UUX,UUY,UUZ)
        USE INIT
        USE NRTYPE
        USE ONERROR
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        INTEGER(SP) :: NP
        INTEGER(SP) :: IERR
        INTEGER(SP) :: II,JJ,KK,N
        REAL(DP) :: UX,UY,UZ
        !INTEGER(SP),DIMENSION(:) :: INDX,INDY
        !INTEGER(SP),DIMENSION(:),OPTIONAL :: INDZ
        REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: UUX,UUY
        REAL(DP),DIMENSION(:),ALLOCATABLE :: VX,VY,VZ
        REAL(DP),DIMENSION(:,:,:),INTENT(INOUT),OPTIONAL :: UUZ
        END SUBROUTINE READVEL
        !SUBROUTINE READVEL(FILENAME,NUMBER_OF_NODES,INDX,UUX,INDY,UUY,INDZ,UUZ)
        !USE INIT
        !USE NRTYPE
        !USE ONERROR
        !REAL(SP),DIMENSION(:),ALLOCATABLE :: UUXX
        !CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        !INTEGER(SP) :: NUMBER_OF_NODES
        !INTEGER(SP) :: IERR
        !INTEGER(SP) :: II,JJ,KK,N
        !REAL(DP) :: UX,UY,UZ
        !INTEGER(SP),DIMENSION(:) :: INDX,INDY
        !INTEGER(SP),DIMENSION(:),OPTIONAL :: INDZ
        !REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: UUX,UUY
        !REAL(DP),DIMENSION(:,:,:),INTENT(INOUT),OPTIONAL :: UUZ
        !!REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: UU
        !!COMPLEX(SPC),DIMENSION(:),ALLOCATABLE :: UUXX
        !END SUBROUTINE READVEL
    END INTERFACE READVEL
    INTERFACE READVEL_BIN
        SUBROUTINE READVEL_BIN(FILENAME,UX,VY,WZ)
        USE INIT
        USE NRTYPE
        USE ONERROR
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        INTEGER(SP) :: NX,NY,NZ,I,J,K
        REAL(SP),DIMENSION(:,:,:) :: UX,VY,WZ
        END SUBROUTINE READVEL_BIN
    END INTERFACE READVEL_BIN
    !###################################################
    TYPE,PUBLIC :: IO
        PRIVATE
        INTEGER(SP) :: NOLINES
        INTEGER(SP) :: NONODES
        INTEGER(SP) :: NONX
        INTEGER(SP) :: NONY
        INTEGER(SP) :: NONZ
        INTEGER(SP) :: DIM
        REAL(DP),DIMENSION(:,:),POINTER :: dummy
        CONTAINS
        PROCEDURE,PUBLIC :: file_stats => read_file_sub
        PROCEDURE,PUBLIC :: NOL => number_of_lines_f
        PROCEDURE,PUBLIC :: NON => number_of_nodes_f
        PROCEDURE,PUBLIC :: NX => number_of_nodes_x_f
        PROCEDURE,PUBLIC :: NY => number_of_nodes_y_f
        PROCEDURE,PUBLIC :: NZ => number_of_nodes_z_f
        PROCEDURE,PUBLIC :: read_file => read_dummy_sub
        PROCEDURE,PUBLIC :: field => field_f
    END TYPE IO
    PRIVATE :: read_file_sub,number_of_nodes_f,number_of_lines_f
    !
    CONTAINS
        SUBROUTINE read_dummy_sub(this,filename,dimen)
        USE INIT
        USE NRTYPE
        USE ONERROR
        IMPLICIT NONE
        CLASS(IO) :: this
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER(SP) :: STATUS_READ
        INTEGER(SP) :: N
        INTEGER(SP) :: dimen
        REAL(DP),DIMENSION(3),TARGET :: DUM1
        REAL(DP),DIMENSION(:,:),ALLOCATABLE,TARGET :: DUM
        TYPE :: ptr
            REAL(DP),DIMENSION(:),POINTER :: p
        END TYPE ptr
        TYPE(ptr) :: p1
        STATUS_READ = 0
        ! First loop to get number of lines
        OPEN(1,FILE=filename,STATUS='UNKNOWN',FORM='FORMATTED')
        DO
            READ(1,*,IOSTAT=STATUS_READ)
            IF(STATUS_READ.NE.0) EXIT
            this%NOLINES = this%NOLINES + 1
        END DO
        this%NONODES = this%NOLINES
        !WRITE(*,*) this%NONODES
        this%NONX = REAL(this%NONODES+1)**(REAL(1.D0/dimen))
        this%NONY = REAL(this%NONODES+1)**(REAL(1.D0/dimen))
        this%NONZ = REAL(this%NONODES+1)**(REAL(1.D0/dimen))
        this%DIM = dimen
        IF (dimen.EQ.3)ALLOCATE(this%dummy(this%NOLINES,3))
        IF (dimen.EQ.2)ALLOCATE(this%dummy(this%NOLINES,2))
        !IF (dimen.EQ.3)ALLOCATE(p1(this%NOLINES))
        !IF (dimen.EQ.2)ALLOCATE(p1(this%NOLINES))
        !IF (ASSOCIATED(this%dummy)) THEN
            IF (dimen.EQ.3) THEN
                REWIND(1)
                DO N=1,this%NONODES-1
                    !ALLOCATE(p1(dimen,N))
                    READ(1,*,IOSTAT=STATUS_READ) DUM1(1),DUM1(2),DUM1(3)
                    p1%p => DUM1
                    this%dummy(N,:) = p1%p(:)
                END DO
            ELSE 
                DO N=1,this%NONODES-1
                    READ(1,*,IOSTAT=STATUS_READ) this%dummy(1,N),this%dummy(2,N)
                END DO
            END IF
        !END IF
        IF (dimen.EQ.3)ALLOCATE(this%dummy(3,this%NOLINES))
        IF (dimen.EQ.2)ALLOCATE(this%dummy(2,this%NOLINES))
        this%dummy => DUM
        CLOSE(1)
        END SUBROUTINE read_dummy_sub
        !
        SUBROUTINE read_file_sub(this,filename,dimen)
        USE NRTYPE
        IMPLICIT NONE
        CLASS(IO) :: this
        CHARACTER(LEN=*),INTENT(IN) :: filename
        INTEGER(SP) :: IERR,STATUS_READ
        INTEGER(SP) :: dimen
        !WRITE(*,*) 'Reading file', filename
        !OPEN(1,FILE=filename,STATUS='UNKNOWN',FORM='FORMATTED')
        !DO
            !READ(1,*,IOSTAT=STATUS_READ)
            !IF(STATUS_READ.NE.0) EXIT
            !this%NOLINES = this%NOLINES + 1
        !END DO
        !CLOSE(1)
        !this%NONODES = this%NOLINES
        !this%NONX = REAL(this%NONODES+1)**(REAL(1.D0/dimen))
        !this%NONY = REAL(this%NONODES+1)**(REAL(1.D0/dimen))
        !this%NONZ = REAL(this%NONODES+1)**(REAL(1.D0/dimen))
        !this%DIM = dimen
        END SUBROUTINE read_file_sub
        !
        FUNCTION field_f(this) RESULT (array)
        USE NRTYPE
        IMPLICIT NONE
        CLASS(IO) :: this
        REAL(DP),DIMENSION(:,:),POINTER :: array
        array => this%dummy
        END FUNCTION field_f
        !
        REAL FUNCTION number_of_lines_f(this)
        IMPLICIT NONE
        CLASS(IO) :: this
        number_of_lines_f = this%NOLINES
        END FUNCTION number_of_lines_f
        !
        REAL FUNCTION number_of_nodes_f(this)
        IMPLICIT NONE
        CLASS(IO) :: this
        number_of_nodes_f = this%NONODES
        END FUNCTION number_of_nodes_f
        !
        REAL FUNCTION number_of_nodes_x_f(this)
        IMPLICIT NONE
        CLASS(IO) :: this
        number_of_nodes_x_f =  this%NONX
        END FUNCTION number_of_nodes_x_f
        !
        REAL FUNCTION number_of_nodes_y_f(this)
        IMPLICIT NONE
        CLASS(IO) :: this
        number_of_nodes_y_f =  this%NONY
        END FUNCTION number_of_nodes_y_f
        !
        REAL FUNCTION number_of_nodes_z_f(this)
        IMPLICIT NONE
        CLASS(IO) :: this
        number_of_nodes_z_f =  this%NONZ
        END FUNCTION number_of_nodes_z_f
        !!
END MODULE IO_CLASS
!=============================================================================
!> @author Felix Dietzsch
!
! DESCRIPTION:
!> Reads coordinate data provided by file 'FILENAME'.
!> Returns array 'COORD'
!
! REVISION HISTORY:
! 25 07 2011 Initial Version
!
!> @param[in] FILENAME Specifies the name of the file containing the node data.
!> @param[out] COORD Array containing node coordinates
!> @param[out] NUMBER_OF_NODES Number of nodes
!=============================================================================
SUBROUTINE READCOORD2D(FILENAME,COORD)
   USE INIT
   USE NRTYPE
   USE ONERROR
   CHARACTER(LEN=*), INTENT(IN) :: FILENAME
   INTEGER(SP) :: NUMBER_OF_NODES
   INTEGER(SP) :: IERR,STATUS_READ
   REAL(DP),DIMENSION(:,:),INTENT(INOUT) :: COORD
   OPEN(1,FILE=FILENAME,FORM='FORMATTED')
   NUMBER_OF_NODES = SIZE(COORD,2)!NUMBER_OF_LINES
   REWIND(1)
   IF (SIZE(COORD,1).EQ.3) THEN
       DO N=1,NUMBER_OF_NODES-1
           READ(1,*,IOSTAT=STATUS_READ) COORD(1,N),COORD(2,N),COORD(3,N)
       END DO
   ELSE 
       DO N=1,NUMBER_OF_NODES-1
           READ(1,*,IOSTAT=STATUS_READ) COORD(1,N),COORD(2,N)
       END DO
   END IF
   CLOSE(1)
END SUBROUTINE READCOORD2D
!=============================================================================
!> @author Felix Dietzsch
!
! DESCRIPTION:
!> Reads dimension of the problem.
!
! REVISION HISTORY:
! 06 12 2011 Initial Version
!
!> @param[in] FILENAME Specifies the name of the file containing the data to be read.
!> @param[out] Array of dimensions of the problem
!=============================================================================
SUBROUTINE READ_DIM(FILEX,FILEY,FILEZ,DIMEN)
  USE NRTYPE
  INTEGER(SP),DIMENSION(:),INTENT(INOUT) :: DIMEN
  CHARACTER(LEN=*),INTENT(IN) :: FILEX
  CHARACTER(LEN=*),INTENT(IN) :: FILEY
  CHARACTER(LEN=*),INTENT(IN) :: FILEZ
  LOGICAL :: FILE_EXISTS
  INTEGER(SP) :: NPX,NPY,NPZ
  INTEGER(SP) :: I

  NPY=1
  NPZ=1
  NPX=1
  INQUIRE(FILE=FILEX,EXIST=FILE_EXISTS)
  IF (FILE_EXISTS) THEN
    OPEN(1,FILE=FILEX,FORM="UNFORMATTED",ACTION="READ")
    REWIND(1)
    READ(1) NPX
    CLOSE(1)
  ELSE
    PRINT*,"ERROR: The file",FILEX,"is missing"
    STOP
  END IF
  INQUIRE(FILE=FILEY,EXIST=FILE_EXISTS)
  IF (FILE_EXISTS) THEN
    OPEN(1,FILE=FILEY,FORM="UNFORMATTED",ACTION="READ")
    REWIND(1)
    READ(1) NPY
    CLOSE(1)
  END IF
  INQUIRE(FILE=FILEZ,EXIST=FILE_EXISTS)
  IF (FILE_EXISTS) THEN
    OPEN(1,FILE=FILEZ,FORM="UNFORMATTED",ACTION="READ")
    REWIND(1)
    READ(1) NPZ
    CLOSE(1)
  END IF
  DIMEN(1)=NPX
  DIMEN(2)=NPY
  DIMEN(3)=NPZ
END SUBROUTINE READ_DIM
!=============================================================================
!> @author Felix Dietzsch
!
! DESCRIPTION:
!> Reads coordinate vectors returns arrays X, Y ,Z.
!
! REVISION HISTORY:
! 06 12 2011 Initial Version
!
!> @param[in] FILENAME Specifies the name of the file containing the data to be read.
!> @param[out] X Array containing x coordinates
!> @param[out] Y Array containing y coordinates 
!> @param[out] Z Array containing z coordinates 
!=============================================================================
SUBROUTINE READ_GRID(FILEX,FILEY,FILEZ,GRIDX,GRIDY,GRIDZ)
  USE NRTYPE
  REAL(DP),DIMENSION(:),INTENT(INOUT) :: GRIDX
  REAL(DP),DIMENSION(:),INTENT(INOUT) :: GRIDY
  REAL(DP),DIMENSION(:),INTENT(INOUT) :: GRIDZ
  CHARACTER(LEN=*),INTENT(IN) :: FILEX
  CHARACTER(LEN=*),INTENT(IN) :: FILEY
  CHARACTER(LEN=*),INTENT(IN) :: FILEZ
  LOGICAL :: FILE_EXISTS
  INTEGER(SP) :: NPX,NPY,NPZ
  INTEGER(SP) :: I

  GRIDX(:)=1.0D0
  GRIDY(:)=1.0D0
  GRIDZ(:)=1.0D0
  INQUIRE(FILE=FILEX,EXIST=FILE_EXISTS)
  IF (FILE_EXISTS) THEN
    OPEN(1,FILE=FILEX,FORM="UNFORMATTED",ACTION="READ")
    REWIND(1)
    READ(1) NPX
    READ(1) (GRIDX(I),I=1,NPX)
    CLOSE(1)
  ELSE
    PRINT*,"ERROR: The file"//FILEX//" is missing"
    STOP
  END IF
  INQUIRE(FILE=FILEY,EXIST=FILE_EXISTS)
  IF (FILE_EXISTS) THEN
    OPEN(1,FILE=FILEY,FORM="UNFORMATTED",ACTION="READ")
    REWIND(1)
    READ(1) NPY
    READ(1) (GRIDY(I),I=1,NPY)
    CLOSE(1)
  END IF
  INQUIRE(FILE=FILEZ,EXIST=FILE_EXISTS)
  IF (FILE_EXISTS) THEN
    OPEN(1,FILE=FILEZ,FORM="UNFORMATTED",ACTION="READ")
    REWIND(1)
    READ(1) NPZ
    READ(1) (GRIDZ(I),I=1,NPZ)
    CLOSE(1)
  END IF
END SUBROUTINE READ_GRID
!=============================================================================
!> @author Felix Dietzsch
!
! DESCRIPTION:
!> Reads temperature data provided by file 'FILENAME'.
!> Returns array 'VAR'
!
! REVISION HISTORY:
! 05 12 2011 Initial Version
!
!> @param[in] FILENAME Specifies the name of the file containing the data to be read.
!> @param[out] VAR Array containing node temperatures
!=============================================================================
SUBROUTINE READVALUE(FILENAME,VAR)
   !USE INIT
   USE NRTYPE
   USE ONERROR
   INTEGER :: I,J,K
   INTEGER :: NPROC_Z,NPOINTS_Z
   REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: VAR
   CHARACTER(LEN=*), INTENT(IN) :: FILENAME
   OPEN(1,FILE=FILENAME,FORM="UNFORMATTED",ACTION="read")
   REWIND(1)
   READ(1) NPROC_Z
   READ(1) NPOINTS_Z
   READ(1) (((VAR(I,J,K),I=1,SIZE(VAR,1)),&
                         J=1,SIZE(VAR,2)),&
                         K=1,SIZE(VAR,3))
   CLOSE(1)
END SUBROUTINE READVALUE
!=============================================================================
!> @author Felix Dietzsch
!
! DESCRIPTION:
!> write out progress variable.
!
! REVISION HISTORY:
! 05 12 2011 Initial Version
!
!> @param[in] FILENAME Specifies the name of the file containing the node data.
!> @param[in] PROG_VAR Specifies the name of the file containing the progress variable data.
!=============================================================================
SUBROUTINE WRITE_PROG(FILENAME,PROG_VAR)
  USE NRTYPE
  INTEGER :: I,J,K
  CHARACTER(LEN=*), INTENT(IN) :: FILENAME
  REAL(DP),DIMENSION(:,:,:),INTENT(IN) :: PROG_VAR
  OPEN(1,FILE=FILENAME,FORM="FORMATTED",ACTION="WRITE")
  DO K=1,SIZE(PROG_VAR,3)
    DO J=1,SIZE(PROG_VAR,2)
      DO I=1,SIZE(PROG_VAR,1)
        WRITE(1,*) PROG_VAR(I,J,K)
      END DO
    END DO
  END DO
  CLOSE(1)
END SUBROUTINE WRITE_PROG
!=============================================================================
!> @author Felix Dietzsch
!
! DESCRIPTION:
!> write out temperature.
!
! REVISION HISTORY:
! 05 12 2011 Initial Version
!
!> @param[in] FILENAME Specifies the name of the file containing the temperature data.
!> @param[in] TEMPER Specifies the name of the file containing the temperature data.
!=============================================================================
SUBROUTINE WRITE_VALUE(FILENAME,value)
  USE NRTYPE
  INTEGER :: I,J,K
  CHARACTER(LEN=*), INTENT(IN) :: FILENAME
  REAL(DP),DIMENSION(:,:,:),INTENT(IN) :: value
  OPEN(1,FILE=FILENAME,FORM="FORMATTED",ACTION="WRITE")
  DO K=1,SIZE(value,3)
    DO J=1,SIZE(value,2)
      DO I=1,SIZE(value,1)
        WRITE(1,*) value(I,J,K)
      END DO
    END DO
  END DO
  CLOSE(1)
END SUBROUTINE WRITE_VALUE
!!=============================================================================
!> @author Felix Dietzsch
!
! DESCRIPTION:
!> Provides vectors that control the indexing into the three
!> velocity component arrays
!
! REVISION HISTORY:
! 25 07 2011 Initial Version
!
!> @param[in] COORD Array containing node coordinates
!> @param[in] NUMBER_OF_NODES Number of nodes
!> @param[out] INDX Vector containing x indices
!> @param[out] INDY Vector containing y indices
!> @param[out] INDZ Vector containing z indices
!> @param[out] NX Number of nodes in x direction 
!> @param[out] NY Number of nodes in y direction 
!> @param[out] NZ Number of nodes in z direction 
!=============================================================================
SUBROUTINE INDEXING(COORD,NUMBER_OF_NODES,INDX,NX,INDY,NY,INDZ,NZ)
    USE NRTYPE
    USE INIT
    REAL(DP) :: XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,LX,LY,LZ
    INTEGER(SP) :: NUMBER_OF_NODES
    INTEGER(SP),INTENT(INOUT) :: NX,NY
    INTEGER(SP),INTENT(INOUT),OPTIONAL :: NZ
    INTEGER(SP) :: IERR
    INTEGER(SP) :: N
    INTEGER(SP) :: NXM,NYM,NZM
    INTEGER(SP),DIMENSION(:),INTENT(INOUT) :: INDX,INDY
    INTEGER(SP),DIMENSION(:),INTENT(INOUT),OPTIONAL :: INDZ
    REAL(DP),DIMENSION(:,:) :: COORD
    IF (PRESENT(INDZ).AND.PRESENT(NZ)) THEN
        NX = REAL(NUMBER_OF_NODES+1)**(1/3.D0)
    ELSE
        NX = REAL(NUMBER_OF_NODES+1)**(1/2.D0)
    END IF
    NY = NX
    NXM = NX -1 ! number of cell in one direction
    NYM = NY -1
    XMIN = MINVAL(COORD(1,:),1)
    XMAX = MAXVAL(COORD(1,:),1)
    YMIN = MINVAL(COORD(2,:),1)
    YMAX = MAXVAL(COORD(2,:),1)
    LX = XMAX-XMIN
    LY = YMAX-YMIN
    IF (PRESENT(INDZ).AND.PRESENT(NZ)) THEN
        NZ = NX
        NZM = NZ - 1
        ZMIN = MINVAL(COORD(3,:),1)
        ZMAX = MAXVAL(COORD(3,:),1)
        LZ = ZMAX-ZMIN
    END IF
    WRITE(*,*) 'Number of cells in on direction', NXM
    DO N=1,NUMBER_OF_NODES-1
       X = COORD(1,N)
       INDX(N) = IDINT(NXM*(X-XMIN)/LX+0.5D0)+1
       Y = COORD(2,N)
       INDY(N) = IDINT(NYM*(Y-YMIN)/LY+0.5D0)+1
       IF (PRESENT(INDZ).AND.PRESENT(NZ)) THEN
           Z = COORD(3,N)
           INDZ(N) = IDINT(NZM*(Z-ZMIN)/LZ+0.5D0)+1
       END IF
    END DO
END SUBROUTINE INDEXING
!=============================================================================
!> @author Felix Dietzsch
!
! DESCRIPTION:
!> Provides arrays for the velocity components
!> 
!
! REVISION HISTORY:
! 25 07 2011 Initial Version
!
!> @param[in] FILENAME Specifies the name of the file containing the velocity data
!> @param[in] NUMBER_OF_NODES Number of nodes
!> @param[in] INDX Vector containing x indices
!> @param[in] INDY Vector containing y indices
!> @param[in] INDZ Vector containing z indices
!> @param[out] UUX Array of x-velocity components
!> @param[out] UUY Array of y-velocity components
!> @param[out] UUZ Array of z-velocity components
!=============================================================================
!SUBROUTINE READVEL(FILENAME,NUMBER_OF_NODES,INDX,UUX,INDY,UUY,INDZ,UUZ)
SUBROUTINE READVEL(FILENAME,UUX,UUY,UUZ)
    USE INIT
    USE NRTYPE
    USE ONERROR
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    INTEGER(SP) :: NP
    INTEGER(SP) :: IERR
    INTEGER(SP) :: II,JJ,KK,N
    REAL(DP) :: UX,UY,UZ
    !INTEGER(SP),DIMENSION(:) :: INDX,INDY
    !INTEGER(SP),DIMENSION(:),OPTIONAL :: INDZ
    REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: UUX,UUY
    REAL(DP),DIMENSION(:),ALLOCATABLE :: VX,VY,VZ
    REAL(DP),DIMENSION(:,:,:),INTENT(INOUT),OPTIONAL :: UUZ
    WRITE(*,*) 'Reading velocity field...'
    NP=SIZE(UUX,1)
    OPEN(2,FILE=FILENAME,STATUS='UNKNOWN')
    IF (PRESENT(UUZ)) THEN
        ALLOCATE(VX(NP**3))
        ALLOCATE(VY(NP**3))
        ALLOCATE(VZ(NP**3))
        DO N=1,NP**3
            READ(2,*) UX,UY,UZ
            VX(N)=UX
            VY(N)=UY
            VZ(N)=UZ
        END DO
        UUX=RESHAPE(VX,(/NP,NP,NP/))
        UUY=RESHAPE(VY,(/NP,NP,NP/))
        UUZ=RESHAPE(VZ,(/NP,NP,NP/))
    ELSE
        ALLOCATE(VX(NP**2))
        ALLOCATE(VY(NP**2))
        DO N=1,NP**2
            READ(2,*) UX,UY
            VX(N)=UX
            VY(N)=UY
        END DO
        UUX=RESHAPE(VX,(/NP,NP,1/))
        UUY=RESHAPE(VY,(/NP,NP,1/))
    END IF
    CLOSE(1)
    IF (PRESENT(UUZ)) THEN
        DEALLOCATE(VX)
        DEALLOCATE(VY)
        DEALLOCATE(VZ)
    ELSE
        DEALLOCATE(VX)
        DEALLOCATE(VY)
    END IF
    WRITE(*,*) 'SUCCESSFUL'
END SUBROUTINE READVEL
!=============================================================================
!> @author Felix Dietzsch
!
! DESCRIPTION:
!> Provides arrays for the velocity components. The data is read from a binary
!> input
!> 
!
! REVISION HISTORY:
! 25 07 2011 Initial Version
!
!> @param[in] FILENAME Specifies the name of the file containing the velocity data
!> @param[out] U Array of x-velocity components
!> @param[out] V Array of y-velocity components
!> @param[out] W Array of z-velocity components
!=============================================================================
SUBROUTINE READVEL_BIN(FILENAME,UX,VY,WZ)
    USE INIT
    USE NRTYPE
    USE ONERROR
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    INTEGER(SP) :: NX,NY,NZ,I,J,K
    REAL(SP),DIMENSION(:,:,:) :: UX,VY,WZ
    NX=128
    NY=128
    NZ=128
    open(1, file=FILENAME,ACCESS='direct',RECL=3*nx,FORM='unformatted')
      do k=1,nz
         do j=1,ny
            read(1,rec= 2+(j-1)+(k-1)*ny)(UX(i,j,k),VY(i,j,k),WZ(i,j,k),i=1,nx)
         end do
      end do
    CLOSE(1)
    WRITE(*,*) 'SUCCESSFUL'
END SUBROUTINE READVEL_BIN
