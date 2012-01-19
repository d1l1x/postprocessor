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
INTEGER(SP),DIMENSION(3) :: DIMEN
CHARACTER(2) :: ITERATOR
REAL(DP),DIMENSION(:),ALLOCATABLE :: GRIDX,GRIDY,GRIDZ 
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: TEMPER
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: PROG_VAR
REAL(DP),DIMENSION(:,:,:,:),ALLOCATABLE :: SPEC
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: UVEL
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: VVEL
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: WVEL
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: DIFCOF
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: VELO
REAL(DP),DIMENSION(:,:,:),ALLOCATABLE ::autocor 
CHARACTER(9) :: OUTDIR
CHARACTER(8) :: INDIR
CHARACTER(8) :: ITERATION

integer(fgsl_size_t)::test

!real(fgsl_double) :: data(5) = (/17.2D0, 18.1D0, 16.5D0, 18.3D0, 12.6D0 /)
real(fgsl_double) :: data(1) = (/1.0D0/)
real(fgsl_double) :: meannn, variancen, largest, smallest
test=1
  meannn     = fgsl_stats_mean(data, 1_fgsl_size_t,test) 
  variancen = fgsl_stats_variance(data, 1_fgsl_size_t, 5_fgsl_size_t)
  largest  = fgsl_stats_max(data, 1_fgsl_size_t, 5_fgsl_size_t)
  smallest = fgsl_stats_min(data, 1_fgsl_size_t, 5_fgsl_size_t)
  
  write(6, '(''The dataset is '',5(F9.5))') data
  write(6, '(''The sample mean is '',F9.5)') meannn
  write(6, '(''The estimated variance is '',F9.5)') variancen
  write(6, '(''The largest value is '',F9.5)') largest
  write(6, '(''The smallest value is '',F9.5)') smallest


!CALL MPI_INIT(IERR)
!CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
!CALL MPI_COMM_RANK(MPI_COMM_WORLD,PROCNUM,IERR)

!IF (PROCNUM.EQ.ROOT) THEN
    OUTDIR = './OUTPUT/'
    INDIR = './INPUT/'
    ITERATION = TRIM("99")
    SPEC_NUM = 17
    WRITE(*,*) '===================================='
    WRITE(*,*) '======== SPECTRAL ANALYSIS ========='
    WRITE(*,*) '===================================='
    
    PRINT*,'READING DIMENSION...'
    CALL READ_DIM("./INPUT/GRIDX_UFORM",&
                  "./INPUT/GRIDY_UFORM",&
                  "./INPUT/GRIDZ_UFORM",&
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
    CALL READ_GRID("./INPUT/GRIDX_UFORM",&
                   "./INPUT/GRIDY_UFORM",&
                   "./INPUT/GRIDZ_UFORM",&
                   GRIDX,GRIDY,GRIDZ)

    ALLOCATE(TEMPER(DIMEN(1),DIMEN(2),DIMEN(3)))
    ALLOCATE(UVEL(DIMEN(1),DIMEN(2),DIMEN(3)))
    ALLOCATE(velo(DIMEN(1),DIMEN(2),DIMEN(3)))
    ALLOCATE(DIFCOF(DIMEN(1),DIMEN(2),DIMEN(3)))
    ALLOCATE(PROG_VAR(DIMEN(1),DIMEN(2),DIMEN(3)))
    ALLOCATE(SPEC(SPEC_NUM,DIMEN(1),DIMEN(2),DIMEN(3)))
    PRINT*,'READING TEMPERATUR...'
    CALL READVALUE(INDIR//'SID_TEMPER_'//ITERATION,TEMPER)
    PRINT*,'READING VELOCITIES...'
    PRINT*,'   READINX X VEL'
    IF (DIM.EQ.1) THEN
        CALL READVALUE(INDIR//'SID_UVEL_'//ITERATION,UVEL)
        VELO(:,:,:)=UVEL(:,:,:)
    ELSEIF (DIM.EQ.2) THEN
        PRINT*,'   READING V VEL'
        ALLOCATE(VVEL(DIMEN(1),DIMEN(2),DIMEN(3)))
        CALL READVALUE(INDIR//'SID_VVEL_'//ITERATION,VVEL)
        VELO(:,:,:) = DSQRT(UVEL(:,:,:)*UVEL(:,:,:)+VVEL(:,:,:)*VVEL(:,:,:))
    ELSEIF (DIM.EQ.3) THEN
        PRINT*,'   READING V VEL'
        ALLOCATE(VVEL(DIMEN(1),DIMEN(2),DIMEN(3)))
        CALL READVALUE(INDIR//'SID_VVEL_'//ITERATION,VVEL)
        PRINT*,'   READING W VEL'
        ALLOCATE(WVEL(DIMEN(1),DIMEN(2),DIMEN(3)))
        CALL READVALUE(INDIR//'SID_WVEL_'//ITERATION,WVEL)
        velo(:,:,:) = DSQRT(UVEL(:,:,:)*UVEL(:,:,:)+VVEL(:,:,:)*VVEL(:,:,:)+WVEL(:,:,:)*WVEL(:,:,:))
    END IF
! read diffusion coefficient data
    PRINT*,'READING DIFFUSION COEF. DATA...'
    CALL READVALUE(INDIR//'SID_D11_'//ITERATION,DIFCOF)
! read species data
    PRINT*,"READING SPECIES DATA..."
    DO I=1,SPEC_NUM
        WRITE(ITERATOR,"(I0)") I 
        PRINT*,'   SPECIES',ITERATOR
        CALL READVALUE(INDIR//'SID_Y'//TRIM(ITERATOR)//'_'//ITERATION,SPEC(I,:,:,:))
    END DO
    ! compute c with Y
    CALL COMP_PROGRESS(TEMPER,PROG_VAR,.TRUE.,SPEC(11,:,:,:),MINVAL(SPEC(11,:,:,:)),MAXVAL(SPEC(11,:,:,:)))
    ! compute c with T
    CALL COMP_PROGRESS(TEMPER,PROG_VAR,.TRUE.)
    !CALL WRITE_COORD('./OUTPUT/GRIDX',GRIDX)
    CALL WRITE_VALUE('./OUTPUT/TEMPER',TEMPER)
    CALL WRITE_VALUE('./OUTPUT/DIFCOF11',DIFCOF)
    CALL WRITE_VALUE('./OUTPUT/VELO',velo)
    PRINT*,"SUCCESSFUL"
    PRINT*,"READING DIFFUSION COEFFICIENTS..."

    
    !################################
    ! process the velocity file
    !IF (DIM.GT.2) THEN
        !NX=NX-1
        !NY=NY-1
        !NZ=NZ-1
        !ALLOCATE(UUX(NX,NY,NZ),STAT=IERR) 
        !IF (IERR.NE.0) CALL ALLOCATION_ERROR(IERR)
        !ALLOCATE(UUY(NX,NY,NZ),STAT=IERR) 
        !IF (IERR.NE.0) CALL ALLOCATION_ERROR(IERR)
        !ALLOCATE(UUZ(NX,NY,NZ),STAT=IERR) 
        !IF (IERR.NE.0) CALL ALLOCATION_ERROR(IERR)
    !ELSE
        !ALLOCATE(UUX(NX,NY,NZ),STAT=IERR) 
        !IF (IERR.NE.0) CALL ALLOCATION_ERROR(IERR)
        !ALLOCATE(UUY(NX,NY,NZ),STAT=IERR) 
        !IF (IERR.NE.0) CALL ALLOCATION_ERROR(IERR)
    !END IF
    !IF (DIM.GT.2) THEN
        !!CALL READVEL(INDIR//'SID_VEL_2D0',NUMBER_OF_NODES,INDX,UUX,INDY,UUY,INDZ,UUZ)
        !CALL READVEL(INDIR//'SID_VEL_2D0',UUX,UUY,UUZ)
    !ELSE
        !!CALL READVEL(INDIR//'SID_VEL_2D0',NUMBER_OF_NODES,INDX,UUX,INDY,UUY)
        !CALL READVEL(INDIR//'SID_VEL_2D0',UUX,UUY)
    !END IF
    !IF (DIM.GT.2) THEN
        !ALLOCATE(UU(NX,NY,NZ,3))
        !UU(:,:,:,1) = UUX(:,:,:)
        !UU(:,:,:,2) = UUY(:,:,:)
        !UU(:,:,:,3) = UUZ(:,:,:)
    !ELSE
        !ALLOCATE(UU(NX,NY,NZ,2))
        !UU(:,:,:,1) = UUX(:,:,:)
        !UU(:,:,:,2) = UUY(:,:,:)
    !END IF
    !################################
    ! computing statistical properties
!    CALL MEAN(UU,NX,NY,NZ,DIM,MEANVAL)
    !WRITE(*,*) 'MEAN: ',MEANVAL
    !!
!!    CALL DEVIAT(UUX,NX,NY,NZ,SIGMA)
    !WRITE(*,*) 'SIGMAX: ',SIGMA
!!    CALL DEVIAT(UUY,NX,NY,NZ,SIGMA)
    !WRITE(*,*) 'SIGMAY: ',SIGMA
    !IF (DIM.GT.2) THEN
!!        CALL DEVIAT(UUZ,NX,NY,NZ,SIGMA)
        !WRITE(*,*) 'SIGMAZ: ',SIGMA
        !WRITE(*,*) 'SUCCESSFUL'
    !END IF
    !!
    !ALLOCATE(DUMMY1(NX,NY,NZ))
    !ALLOCATE(DUMMY2(NX,NY,NZ))
    !ALLOCATE(TEST_OUT(NX,NY,NZ))
    !DUMMY1(:,:,:)=DCMPLX(velo(:,:,:))
    !DUMMY2(:,:,:)=DCMPLX(UUY(:,:,:),0.0)
    !#########################################
    !!########## AUTO CORRELATION ############
    PRINT*,"compute correlation"
    ALLOCATE(autocor(dimen(1),dimen(2),dimen(3)))
    CALL CORREL(velo,autocor)
    !!ALLOCATE(TEST_IN(NX,NY))
    !!TEST_IN(:,:)=TEST_OUT(:,:,1)
    !!CALL FFTSHIFT(TEST_IN,NX,NY) ! 2D shift of fourier transform
    CALL WRITE_VALUE('./OUTPUT/AUTOCOR',autocor)
    !DEALLOCATE(TEST_IN)
    !!#########################################
    !########## CROSS CORRELATION ############
    !DUMMY1(:,:,:)=DCMPLX(UUX(:,:,:),0.0)
    !DUMMY2(:,:,:)=DCMPLX(UUX(:,:,:),0.0)
    !PRINT*, "Calling correl function"
!!    CALL CORREL(DUMMY1,DUMMY2,NX,NY,NZ,TEST_OUT)
    !PRINT*, "SUCCESSFUL"
    !OPEN(1,FILE=OUTDIR//'XCORRELATION.OUT')
    !DO N=1,NY
        !WRITE(1,*) REAL(TEST_OUT(1,N,1))
    !END DO
    !CLOSE(1)
!END IF
!if(procnum.eq.root) then
    !write(*,*) NX,NY,NZ
!end if
!if(procnum.eq.root)then
    !allocate(sendbuf(nproc,3))
    !sendbuf(:,1) = NX
    !sendbuf(:,2) = NY
    !sendbuf(:,3) = NZ
!end if
!allocate(recvbuf(3))
!CALL MPI_SCATTER(sendbuf(:,1),1,MPI_INTEGER,recvbuf(1),1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
!CALL MPI_SCATTER(sendbuf(:,2),1,MPI_INTEGER,recvbuf(2),1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
!CALL MPI_SCATTER(sendbuf(:,3),1,MPI_INTEGER,recvbuf(3),1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
!NX = recvbuf(1)
!NY = recvbuf(2)
!NZ = recvbuf(3)
!deallocate(sendbuf)
!deallocate(recvbuf)
!if(procnum.ne.root) then
    !write(*,*) NX,NY,NZ
!end if
!IF (DIM.GT.2) THEN
    !NX=NX-1
    !NY=NY-1
    !NZ=NZ-1
    !ALLOCATE(UUX_loc(NX,NY,NZ),STAT=IERR) 
    !IF (IERR.NE.0) CALL ALLOCATION_ERROR(IERR)
    !ALLOCATE(UUY_loc(NX,NY,NZ),STAT=IERR) 
    !IF (IERR.NE.0) CALL ALLOCATION_ERROR(IERR)
    !ALLOCATE(UUZ_loc(NX,NY,NZ),STAT=IERR) 
    !IF (IERR.NE.0) CALL ALLOCATION_ERROR(IERR)
!ELSE
    !ALLOCATE(UUX_loc(NX,NY,NZ),STAT=IERR) 
    !IF (IERR.NE.0) CALL ALLOCATION_ERROR(IERR)
    !ALLOCATE(UUY_loc(NX,NY,NZ),STAT=IERR) 
    !IF (IERR.NE.0) CALL ALLOCATION_ERROR(IERR)
!END IF
!!!!IF (PROCNUM.EQ.ROOT) THEN
!!!!    CALL t%start_timer()
!!!!    if(procnum.eq.root) then
!!!!        if (allocated(uux)) then
!!!!            CALL MPI_SCATTER(UUX,NX*NY*NZ,MPI_INTEGER,UUX_loc,NX*NY*NZ,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
!!!!        end if
!!!!        if (allocated(uuy)) then
!!!!            CALL MPI_SCATTER(UUY,NX*NY*NZ,MPI_INTEGER,UUY_loc,NX*NY*NZ,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
!!!!        end if
!!!!        if (allocated(uuz)) then
!!!!            CALL MPI_SCATTER(UUZ,NX*NY*NZ,MPI_INTEGER,UUZ_loc,NX*NY*NZ,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
!!!!        end if
!!!!    end if
!!!!end if
!!!!ALLOCATE(OUT(NX/2))
!!!!!IF (DIM.GT.2) THEN
!!!!    CALL SPECTR(PROCNUM,UUX,UUY,UUZ,NX,NY,NZ,OUT)
    !OPEN(1,FILE=OUTDIR//'SPECTRUM.OUT')
        !DO I=1,NX/2
            !WRITE(1,*) REAL(OUT(I))
        !END DO
    !CLOSE(1)
!ELSE
    !CALL SPECTR(UUX,UUY,NX,NY,OUT)
    !OPEN(1,FILE=OUTDIR//'SPECTRUM.OUT')
        !DO I=1,NX/2
            !WRITE(1,*) REAL(OUT(I))
        !END DO
    !CLOSE(1)
!END IF
!DEALLOCATE(OUT)
!ALLOCATE(OUT(2))
!WRITE (*,'(A,F8.3,A)') 'Time =', t%ELAPSED_TIME(), 's'



!ALLOCATE(SPECTRUM(2,NXM+1))
!IF (IERR.NE.0) CALL ALLOCATION_ERROR(IERR)
!! Deallocation stage
!WRITE(*,*) 'Entering deallocation phase...'
!!WRITE(*,*) '    Deallocate COORD'
!!UNIT_VALUE=1
!!DEALLOCATE(COORD,STAT=IERR)
!!IF (IERR.NE.0) CALL DEALLOCATION_ERROR(IERR,UNIT_VALUE)
!WRITE(*,*) '    Deallocate INDX'
!UNIT_VALUE=2
!DEALLOCATE(INDX,STAT=IERR)
!IF (IERR.NE.0) CALL DEALLOCATION_ERROR(IERR,UNIT_VALUE)
!WRITE(*,*) '    Deallocate INDY'
!UNIT_VALUE=3
!DEALLOCATE(INDY,STAT=IERR)
!IF (IERR.NE.0) CALL DEALLOCATION_ERROR(IERR,UNIT_VALUE)
!WRITE(*,*) '    Deallocate UU'
!UNIT_VALUE=4
!DEALLOCATE(UU,STAT=IERR)
!IF (IERR.NE.0) CALL DEALLOCATION_ERROR(IERR,UNIT_VALUE)
!WRITE(*,*) '    Deallocate UUX'
!UNIT_VALUE=5
!DEALLOCATE(UUX,STAT=IERR)
!IF (IERR.NE.0) CALL DEALLOCATION_ERROR(IERR,UNIT_VALUE)
!WRITE(*,*) '    Deallocate UUY'
!UNIT_VALUE=6
!DEALLOCATE(UUY,STAT=IERR)
!IF (IERR.NE.0) CALL DEALLOCATION_ERROR(IERR,UNIT_VALUE)
!WRITE(*,*) '    Deallocate R11'
!UNIT_VALUE=7
!DEALLOCATE(R11,STAT=IERR)
!IF (IERR.NE.0) CALL DEALLOCATION_ERROR(IERR,UNIT_VALUE)
!WRITE(*,*) '    Deallocate SPECTRUM'
!UNIT_VALUE=8
!DEALLOCATE(SPECTRUM,STAT=IERR)
!IF (IERR.NE.0) CALL DEALLOCATION_ERROR(IERR,UNIT_VALUE)
!WRITE(*,*) 'Finished deallocation phase'
!CALL MPI_FINALIZE(IERR)
END PROGRAM

