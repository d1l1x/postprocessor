PROGRAM SPECTRAL_ANALYSIS

USE FGSL
USE INIT
USE ONERROR
USE NRTYPE
USE STATISTICS
USE IO_CLASS
USE PREMIXED_CLASS
USE,INTRINSIC :: iso_c_binding

IMPLICIT NONE
!include "mpif.h"
include "fftw3.f03"

INTEGER :: I 

INTEGER(SP) :: SPEC_NUM
REAL(SP) :: int_temp
INTEGER(SP),DIMENSION(3) :: DIMEN
LOGICAL :: saving
CHARACTER(2) :: ITERATOR
REAL(DP),DIMENSION(:),ALLOCATABLE :: GRIDX,GRIDY,GRIDZ 
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: TEMPER
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: PROG_VAR
REAL(DP),DIMENSION(:,:,:,:),ALLOCATABLE :: SPEC
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: UVEL
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: VVEL
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: WVEL
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: DIFCOF
REAL(DP),DIMENSION(:,:,:,:),ALLOCATABLE :: velo
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: temp
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE ::autocor 
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE ::temp_3d
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: spectrum
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: kappa
REAL(DP),DIMENSION(:),ALLOCATABLE ::temp_1d
REAL(DP),DIMENSION(:),ALLOCATABLE ::c_limit_array
LOGICAL,DIMENSION(:,:,:),ALLOCATABLE ::c_limit_index
CHARACTER(12) :: OUTDIR
CHARACTER(11) :: INDIR
CHARACTER(8) :: ITERATION

saving=.FALSE.

!integer(fgsl_size_t)::test
!
!!real(fgsl_double) :: data(5) = (/17.2D0, 18.1D0, 16.5D0, 18.3D0, 12.6D0 /)
!real(fgsl_double) :: data(1) = (/1.0D0/)
!real(fgsl_double) :: meannn, variancen, largest, smallest
!test=1
  !meannn     = fgsl_stats_mean(data, 1_fgsl_size_t,test) 
  !variancen = fgsl_stats_variance(data, 1_fgsl_size_t, 5_fgsl_size_t)
  !largest  = fgsl_stats_max(data, 1_fgsl_size_t, 5_fgsl_size_t)
  !smallest = fgsl_stats_min(data, 1_fgsl_size_t, 5_fgsl_size_t)
 ! 
  !write(6, '(''The dataset is '',5(F9.5))') data
  !write(6, '(''The sample mean is '',F9.5)') meannn
  !write(6, '(''The estimated variance is '',F9.5)') variancen
  !write(6, '(''The largest value is '',F9.5)') largest
  !write(6, '(''The smallest value is '',F9.5)') smallest
!

!CALL MPI_INIT(IERR)
!CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
!CALL MPI_COMM_RANK(MPI_COMM_WORLD,PROCNUM,IERR)

OUTDIR = './OUTPUT/3D/'
INDIR = "./INPUT/3D/"
ITERATION = TRIM("00000000")
SPEC_NUM = 17
WRITE(*,*) '===================================='
WRITE(*,*) '======== SPECTRAL ANALYSIS ========='
WRITE(*,*) '===================================='

PRINT*,'READING DIMENSION...'
print*,indir
CALL READ_DIM(indir//"SID_GRIDX_UFORM"//iteration,&
              indir//"SID_GRIDY_UFORM"//iteration,&
              indir//"SID_GRIDZ_UFORM"//iteration,&
               DIMEN)
IF (dimen(3).EQ.1.AND.dimen(2).EQ.1) THEN
    PRINT*,"   We have a 1D problem"
    PRINT*,"   Number of grid points in x: ",dimen(1)
    DIM=1
ELSEIF (dimen(3).EQ.1) THEN
    PRINT*,"   We have a 2D problem"
    PRINT*,"   Number of grid points in x: ",dimen(1)
    PRINT*,"   Number of grid points in y: ",dimen(2)
    DIM=2
ELSE
    PRINT*,"   We have a 3D problem"
    PRINT*,"   Number of grid points in x: ",dimen(1)
    PRINT*,"   Number of grid points in y: ",dimen(1)
    PRINT*,"   Number of grid points in z: ",dimen(1)
    DIM=3
END IF
PRINT*,'READING GRID...'
ALLOCATE(GRIDX(DIMEN(1)))
ALLOCATE(GRIDY(DIMEN(2)))
ALLOCATE(GRIDZ(DIMEN(3)))
CALL READ_GRID(indir//"SID_GRIDX_UFORM"//iteration,&
               indir//"SID_GRIDY_UFORM"//iteration,&
               indir//"SID_GRIDZ_UFORM"//iteration,&
               GRIDX,GRIDY,GRIDZ)

ALLOCATE(TEMPER(DIMEN(1),DIMEN(2),DIMEN(3)))
ALLOCATE(UVEL(DIMEN(1),DIMEN(2),DIMEN(3)))
ALLOCATE(DIFCOF(DIMEN(1),DIMEN(2),DIMEN(3)))
ALLOCATE(PROG_VAR(DIMEN(1),DIMEN(2),DIMEN(3)))
ALLOCATE(SPEC(SPEC_NUM,DIMEN(1),DIMEN(2),DIMEN(3)))
ALLOCATE(c_limit_index(DIMEN(1),DIMEN(2),DIMEN(3)))

PRINT*,'READING VELOCITIES...'
IF (DIM.EQ.1) THEN
    ALLOCATE(velo(DIMEN(1),DIMEN(2),DIMEN(3),1))
    PRINT*,'   READINX X VEL'
    CALL READ_VAL(INDIR//'SID_UVEL_'//ITERATION,UVEL)
    velo(:,:,:,1)=uvel(:,:,:)
ELSEIF (DIM.EQ.2) THEN
    ALLOCATE(velo(DIMEN(1),DIMEN(2),DIMEN(3),2))
    PRINT*,'   READINX X VEL'
    CALL READ_VAL(INDIR//'SID_UVEL_'//ITERATION,UVEL)
    velo(:,:,:,1)=uvel(:,:,:)
    PRINT*,'   READING V VEL'
    ALLOCATE(VVEL(DIMEN(1),DIMEN(2),DIMEN(3)))
    CALL READ_VAL(INDIR//'SID_VVEL_'//ITERATION,VVEL)
    velo(:,:,:,2)=vvel(:,:,:)
ELSEIF (DIM.EQ.3) THEN
    ALLOCATE(velo(DIMEN(1),DIMEN(2),DIMEN(3),3))
    PRINT*,'   READINX X VEL'
    CALL READ_VAL(INDIR//'SID_UVEL_'//ITERATION,UVEL,'3D')
    velo(:,:,:,1)=uvel(:,:,:)
    PRINT*,'   READING V VEL'
    ALLOCATE(VVEL(DIMEN(1),DIMEN(2),DIMEN(3)))
    CALL READ_VAL(INDIR//'SID_VVEL_'//ITERATION,VVEL,'3D')
    velo(:,:,:,2)=vvel(:,:,:)
    PRINT*,'   READING W VEL'
    ALLOCATE(WVEL(DIMEN(1),DIMEN(2),DIMEN(3)))
    CALL READ_VAL(INDIR//'SID_WVEL_'//ITERATION,WVEL,'3D')
    velo(:,:,:,3)=wvel(:,:,:)
END IF
!
IF (DIM.EQ.1.OR.DIM.EQ.2) THEN
! read temperature data
    PRINT*,'READING TEMPERATUR...'
    CALL READ_VAL(INDIR//'SID_TEMPER_'//ITERATION,TEMPER)
! read diffusion coefficient data
    PRINT*,'READING DIFFUSION COEF. DATA...'
    CALL READ_VAL(INDIR//'SID_D11_'//ITERATION,DIFCOF)
! read species data
    PRINT*,"READING SPECIES DATA..."
    DO I=1,SPEC_NUM
        WRITE(ITERATOR,"(I0)") I 
        PRINT*,'   SPECIES',ITERATOR
        CALL READ_VAL(INDIR//'SID_Y'//TRIM(ITERATOR)//'_'//ITERATION,SPEC(I,:,:,:))
    END DO
ELSE IF (DIM.EQ.3) THEN
! read temperature data
    PRINT*,'READING TEMPERATUR...'
    CALL READ_VAL(INDIR//'SID_TEMPER_'//ITERATION,TEMPER,'3D')
! read diffusion coefficient data
    PRINT*,'READING DIFFUSION COEF. DATA...'
    CALL READ_VAL(INDIR//'SID_D11_'//ITERATION,DIFCOF,'3D')
! read species data
    PRINT*,"READING SPECIES DATA..."
    DO I=1,SPEC_NUM
        WRITE(ITERATOR,"(I0)") I 
        PRINT*,'   SPECIES',ITERATOR
        CALL READ_VAL(INDIR//'SID_Y'//TRIM(ITERATOR)//'_'//ITERATION,SPEC(I,:,:,:),'3D')
    END DO
END IF

! compute c with Y
int_temp = 0.1
!CALL COMP_PROGRESS(TEMPER,PROG_VAR,.TRUE.,SPEC(11,:,:,:),MINVAL(SPEC(11,:,:,:)),MAXVAL(SPEC(11,:,:,:)),int_temp,c_limit_index)

!ALLOCATE(c_limit_array(COUNT(c_limit_index.EQV..TRUE.)))
!CALL GET_VAL_UNBURNT(velo,c_limit_index,c_limit_array)

! compute c with T
IF (saving) THEN
    PRINT*,"Write out data..."
    CALL COMP_PROGRESS(TEMPER,PROG_VAR,.TRUE.)
    PRINT*,"   Temperature"
    CALL WRITE_VALUE(OUTDIR//'TEMPER',TEMPER)
    PRINT*,"   Diffusion coefficient"
    CALL WRITE_VALUE(OUTDIR//'DIFCOF11',DIFCOF)
END IF
    PRINT*,"   X velocity"
    CALL WRITE_VALUE(OUTDIR//'uvel',uvel)
    PRINT*,"   Y velocity"
    CALL WRITE_VALUE(OUTDIR//'vvel',vvel)
    PRINT*,"   Z velocity"
    CALL WRITE_VALUE(OUTDIR//'wvel',wvel)
    PRINT*,"successfull"

!!########## AUTO CORRELATION ############
PRINT*,"compute correlation"
ALLOCATE(autocor(dimen(1),dimen(2),dimen(3)))
!ALLOCATE(TEMP(size(c_limit_array,1),1,1))
!temp(:,1,1) = c_limit_array(:)
CALL CORREL(uvel,autocor)
IF (saving) THEN
    CALL WRITE_VALUE(OUTDIR//'AUTOCORU',autocor)
END IF
CALL CORREL(vvel,autocor)
IF (saving) THEN
    CALL WRITE_VALUE(OUTDIR//'AUTOCORV',autocor)
END IF
CALL CORREL(wvel,autocor)
IF (saving) THEN
    CALL WRITE_VALUE(OUTDIR//'AUTOCORW',autocor)
END IF
!!########## 3D Spectrum ############
ALLOCATE(spectrum(dimen(1)*dimen(2)*dimen(3)/8,1,1))
ALLOCATE(kappa(dimen(1)*dimen(2)*dimen(3)/8,1,1))
CALL SPECTR(velo,spectrum,kappa)
CALL WRITE_VALUE(outdir//'spectrum',spectrum)
CALL WRITE_VALUE(outdir//'kappa',kappa)

END PROGRAM

