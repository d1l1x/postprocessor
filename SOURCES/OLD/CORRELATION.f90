!=======================================
! numerical recipes page 1361
!MODULE NRTYPE
    !INTEGER,PARAMETER :: SPC =KIND((1.0,1.0))
    !INTEGER,PARAMETER :: DPC =KIND((1.0D0,1.0D0))
    !INTEGER,PARAMETER :: SP = KIND(1.0)
    !INTEGER,PARAMETER :: DP = KIND(1.0D0)
    !INTEGER,PARAMETER :: I4B = SELECTED_INT_KIND(9)
    !INTEGER, PARAMETER :: LGT = KIND(.true.)
    !REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
    !REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
!END MODULE NRTYPE
!!=======================================
! numerical recipes page 1364
MODULE NRUTIL
    USE NRTYPE
    INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
    INTERFACE swap
        MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
        swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
        masked_swap_rs,masked_swap_rv,masked_swap_rm
    END INTERFACE
    INTERFACE ASSERT
        MODULE PROCEDURE ASSERT1
    END INTERFACE
    INTERFACE ASSERT_EQ
        MODULE PROCEDURE ASSERT_EQ2
    END INTERFACE
    INTERFACE arth
        MODULE PROCEDURE arth_r, arth_d, arth_i
    END INTERFACE
    CONTAINS
        SUBROUTINE ASSERT1(N1,STRING)
            CHARACTER(LEN=*), INTENT(IN) :: STRING
            LOGICAL, INTENT(IN) :: N1
            IF (.NOT.N1) THEN
                WRITE(*,*) 'NERROR: An assertion failed with this tag: ',&
                    STRING
                STOP 'Program terminated by ASSERT2'
            END IF
        END SUBROUTINE ASSERT1
        !
        SUBROUTINE swap_i(a,b)
        INTEGER(I4B), INTENT(INOUT) :: a,b
        INTEGER(I4B) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_i
        SUBROUTINE swap_r(a,b)
        REAL(SP), INTENT(INOUT) :: a,b
        REAL(SP) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_r
        SUBROUTINE swap_rv(a,b)
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
        REAL(SP), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_rv
        SUBROUTINE swap_c(a,b)
        COMPLEX(SPC), INTENT(INOUT) :: a,b
        COMPLEX(SPC) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_c
        SUBROUTINE swap_cv(a,b)
        COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
        COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_cv
        SUBROUTINE swap_cm(a,b)
        COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
        COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_cm
        SUBROUTINE swap_z(a,b)
        COMPLEX(DPC), INTENT(INOUT) :: a,b
        COMPLEX(DPC) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_z
        SUBROUTINE swap_zv(a,b)
        COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
        COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_zv
        SUBROUTINE swap_zm(a,b)
        COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
        COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_zm
        SUBROUTINE masked_swap_rs(a,b,mask)
        REAL(SP), INTENT(INOUT) :: a,b
        LOGICAL(LGT), INTENT(IN) :: mask
        REAL(SP) :: swp
        if (mask) then
            swp=a
            a=b
            b=swp
        end if
        END SUBROUTINE masked_swap_rs
        SUBROUTINE masked_swap_rv(a,b,mask)
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(size(a)) :: swp
        where (mask)
            swp=a
            a=b
            b=swp
        end where
        END SUBROUTINE masked_swap_rv
        SUBROUTINE masked_swap_rm(a,b,mask)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
        where (mask)
            swp=a
            a=b
            b=swp
        end where
        END SUBROUTINE masked_swap_rm
        !
        FUNCTION ASSERT_EQ2(N1,N2,STRING)
            CHARACTER(LEN=*), INTENT(IN) :: string
            INTEGER, INTENT(IN) :: N1,N2
            INTEGER :: ASSERT_EQ2
            IF (N1.EQ.N2) THEN
                ASSERT_EQ2=N1
            ELSE
                WRITE(*,*) 'NRERROR: An ASSERT_EQ failed with this tag: ',&
                    STRING
                STOP 'Program terminated by ASSERT_EQ2'
            END IF
        END FUNCTION ASSERT_EQ2
        FUNCTION arth_r(first,increment,n)
        !Array function returning an arithmetic progression.
        REAL(SP), INTENT(IN) :: first,increment
        INTEGER(I4B), INTENT(IN) :: n
        REAL(SP), DIMENSION(n) :: arth_r
        INTEGER(I4B) :: k,k2
        REAL(SP) :: temp
        if (n > 0) arth_r(1)=first
        if (n <= NPAR_ARTH) then
            do k=2,n
            arth_r(k)=arth_r(k-1)+increment
            end do
        else
            do k=2,NPAR2_ARTH
                arth_r(k)=arth_r(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
                if (k >= n) exit
                k2=k+k
                arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
                temp=temp+temp
                k=k2
            end do
        end if
        END FUNCTION ARTH_R
        FUNCTION arth_d(first,increment,n)
        REAL(DP), INTENT(IN) :: first,increment
        INTEGER(I4B), INTENT(IN) :: n
        REAL(DP), DIMENSION(n) :: arth_d
        INTEGER(I4B) :: k,k2
        REAL(DP) :: temp
        if (n > 0) arth_d(1)=first
        if (n <= NPAR_ARTH) then
            do k=2,n
            arth_d(k)=arth_d(k-1)+increment
            end do
        else
            do k=2,NPAR2_ARTH
            arth_d(k)=arth_d(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
            if (k >= n) exit
            k2=k+k
            arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
            temp=temp+temp
            k=k2
            end do
        end if
        END FUNCTION arth_d
        FUNCTION arth_i(first,increment,n)
        INTEGER(I4B), INTENT(IN) :: first,increment,n
        INTEGER(I4B), DIMENSION(n) :: arth_i
        INTEGER(I4B) :: k,k2,temp
        if (n > 0) arth_i(1)=first
        if (n <= NPAR_ARTH) then
            do k=2,n
            arth_i(k)=arth_i(k-1)+increment
            end do
        else
            do k=2,NPAR2_ARTH
            arth_i(k)=arth_i(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
            if (k >= n) exit
            k2=k+k
            arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
            temp=temp+temp
            k=k2
            end do
        end if
        END FUNCTION arth_i
        FUNCTION zroots_unity(n,nn)
        !Complex function returning nn powers of the nth root of unity.
        INTEGER(I4B), INTENT(IN) :: n,nn
        COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
        INTEGER(I4B) :: k
        REAL(SP) :: theta
        zroots_unity(1)=1.0
        theta=TWOPI/n
        k=1
        do
        if (k >= nn) exit
        zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
        zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
        zroots_unity(2:min(k,nn-k))
        k=2*k
        end do
        END FUNCTION zroots_unity
END MODULE NRUTIL
!=================================
MODULE NR
    INTERFACE fourrow
        SUBROUTINE fourrow_dp(data,isign)
        USE nrtype
        COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
        INTEGER(I4B), INTENT(IN) :: isign
        END SUBROUTINE fourrow_dp
        SUBROUTINE fourrow_sp(data,isign)
        USE nrtype
        COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
        INTEGER(I4B), INTENT(IN) :: isign
        END SUBROUTINE fourrow_sp
    END INTERFACE
    INTERFACE REALFT
        SUBROUTINE REALFT_DP(DATA,ISIGN,ZDATA)
        USE NRTYPE
        REAL(DP),DIMENSION(:),INTENT(INOUT) :: DATA
        INTEGER(I4B),INTENT(IN) :: ISIGN
        COMPLEX(DPC),DIMENSION(:),OPTIONAL,TARGET :: ZDATA
        END SUBROUTINE REALFT_DP
        !
        SUBROUTINE REALFT_SP(DATA,ISIGN,ZDATA)
        USE NRTYPE
        REAL(SP),DIMENSION(:),INTENT(INOUT) :: DATA
        INTEGER(I4B),INTENT(IN) :: ISIGN
        COMPLEX(SPC),DIMENSION(:),OPTIONAL,TARGET :: ZDATA
        END SUBROUTINE REALFT_SP
    END INTERFACE REALFT
    INTERFACE four1
        SUBROUTINE four1_dp(data,isign)
        USE nrtype
        COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
        INTEGER(I4B), INTENT(IN) :: isign
        END SUBROUTINE four1_dp
        SUBROUTINE four1_sp(data,isign)
        USE nrtype
        COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
        INTEGER(I4B), INTENT(IN) :: isign
        END SUBROUTINE four1_sp
    END INTERFACE
    INTERFACE CORRELATION
        SUBROUTINE CORRELATION(DATA1,DATA2)
        USE NRTYPE
        INTEGER(I4B) :: NO2,N
        REAL,DIMENSION(:),INTENT(INOUT) :: DATA1,DATA2
        REAL(SP),DIMENSION(SIZE(DATA1)) :: CORR
        COMPLEX(SPC),DIMENSION(SIZE(DATA1)/2) :: CDAT1,CDAT2
        END SUBROUTINE
    END INTERFACE CORRELATION
END MODULE NR
!===================================
SUBROUTINE fourrow_sp(data,isign)
        USE nrtype; USE nrutil, ONLY: assert,swap
        IMPLICIT NONE
        COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
        INTEGER(I4B), INTENT(IN) :: isign
        INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
        REAL(DP) :: theta
        COMPLEX(SPC), DIMENSION(size(data,1)) :: temp
        COMPLEX(DPC) :: w,wp !Double precision for the trigonometric recurrences.
        COMPLEX(SPC) :: ws
        n=size(data,2)
        call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_sp')
        n2=n/2
        j=n2
        !This is the bit-reversal section of the routine.
        do i=1,n-2
        if (j > i) call swap(data(:,j+1),data(:,i+1))
        m=n2
        do
        if (m < 2 .or. j < m) exit
        j=j-m
        m=m/2
        end do
        j=j+m
        end do
        mmax=1
        !Here begins the Danielson-Lanczos section of the routine.
        do !Outer loop executed log2 N times.
        if (n <= mmax) exit
        istep=2*mmax
        theta=PI_D/(isign*mmax) !Initialize for the trigonometric recurrence.
        wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
        w=cmplx(1.0_dp,0.0_dp,kind=dpc)
        do m=1,mmax !Here are the two nested inner loops.
        ws=w
        do i=m,n,istep
        j=i+mmax
        temp=ws*data(:,j) !This is the Danielson-Lanczos formula.
        data(:,j)=data(:,i)-temp
        data(:,i)=data(:,i)+temp
        end do
        w=w*wp+w !Trigonometric recurrence.
        end do
        mmax=istep
        END DO
END SUBROUTINE fourrow_sp
SUBROUTINE four1_sp(data,isign)
    USE nrtype; USE nrutil, ONLY: arth,assert
    USE nr, ONLY: fourrow
    IMPLICIT NONE
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
    REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
    INTEGER(I4B) :: n,m1,m2,j
    n=size(data)
    call assert(iand(n,n-1)==0,'n must be a power of 2 in four1_sp')
    !Find dimensions as close to square as possible, allocate space, and reshape the
    !input array.
    m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
    m2=n/m1
    allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat=reshape(data,shape(dat))
    call fourrow(dat,isign)
    theta=arth(0,isign,m1)*TWOPI_D/n !Set up recurrence.
    wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
    w=cmplx(1.0_dp,0.0_dp,kind=dpc)
    do j=2,m2 !Multiply by the extra phase factor.
    w=w*wp+w
    dat(:,j)=dat(:,j)*w
    end do
    temp=transpose(dat) !Transpose, and transform on (original) first in
    call fourrow(temp,isign)
    data=reshape(temp,shape(data)) !Reshape the result back to one dimension.
    deallocate(dat,w,wp,theta,temp)
END SUBROUTINE four1_sp

SUBROUTINE realft_sp(data,isign,zdata)
USE nrtype
USE nrutil, ONLY: assert,assert_eq,zroots_unity
USE nr, ONLY: four1
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
INTEGER(I4B), INTENT(IN) :: isign
COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
!When isign = 1, calculates the Fourier transform of a set of N real-valued data
!points,
!input in the array data. If the optional argument zdata is not present, the data
!are replaced
!by the positive frequency half of its complex Fourier transform. The real-valued
!first and
!last components of the complex transform are returned as elements data(1) and
!data(2),
!respectively. If the complex array zdata of length N/2 is present, data is
!unchanged and
!the transform is returned in zdata. N must be a power of 2. If isign = −1, this
!routine
!calculates the inverse transform of a complex data array if it is the transform
!of real data.
!(Result in this case must be multiplied by 2/N.) The data can be supplied either
!in data,
!with zdata absent, or in zdata.
INTEGER(I4B) :: n,ndum,nh,nq
COMPLEX(SPC), DIMENSION(size(data)/4) :: w
COMPLEX(SPC), DIMENSION(size(data)/4-1) :: h1,h2
COMPLEX(SPC), DIMENSION(:), POINTER :: cdata !Used for internal complex computations
COMPLEX(SPC) :: z
REAL(SP) :: c1=0.5_sp,c2
n=size(data)
call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_sp')
nh=n/2
nq=n/4
if (present(zdata)) then
    ndum=assert_eq(n/2,size(zdata),'realft_sp')
    cdata=>zdata !Use zdata as cdata.
    if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
else
    allocate(cdata(n/2)) !Have to allocate storage ourselves.
    cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
end if
if (isign == 1) then
    c2=-0.5_sp
    call four1(cdata,+1) !The forward transform is here.
else !Otherwise set up for an inverse trans
    c2=0.5_sp !form.
end if
w=zroots_unity(sign(n,isign),n/4)
w=cmplx(-aimag(w),real(w),kind=spc)
h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1))) !The two separate transforms are
h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1))) !arated out of cdata.
!Next they are recombined to form the true transform of the original real data:
cdata(2:nq)=h1+w(2:nq)*h2
cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
z=cdata(1) !Squeeze the first and last data together
!to get them all within the
!original array.
if (isign == 1) then
    cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=spc)
else
    cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=spc)
    call four1(cdata,-1) !This is the inverse transform for the
end if                      !case isign=-1.
if (present(zdata)) then !Ship out answer in data if required.
    if (isign /= 1) then
        data(1:n-1:2)=real(cdata)
        data(2:n:2)=aimag(cdata)
    end if
else
    data(1:n-1:2)=real(cdata)
    data(2:n:2)=aimag(cdata)
    deallocate(cdata)
end if
END SUBROUTINE REALFT_SP

SUBROUTINE CORRELATION(DATA1,DATA2)
    USE NRTYPE
    USE NRUTIL,ONLY:assert,assert_eq
    USE NR,ONLY:REALFT
    IMPLICIT NONE
    INTEGER(I4B) :: NO2,N
    REAL,DIMENSION(:),INTENT(INOUT) :: DATA1,DATA2
    REAL(SP),DIMENSION(SIZE(DATA1)) :: CORR
    COMPLEX(SPC),DIMENSION(SIZE(DATA1)/2) :: CDAT1,CDAT2
    N=ASSERT_EQ(SIZE(DATA1),SIZE(DATA2),'CORRELATION')
    CALL ASSERT(IAND(N,N-1)==0,'N must be a power of 2 in CORRELATION')
    NO2=N/2
    CALL REALFT(DATA1,1,CDAT1)
    CALL REALFT(DATA2,1,CDAT2)
    CDAT1(1)=CMPLX(REAL(CDAT1(1))*REAL(CDAT2(1))/NO2,AIMAG(CDAT1(1))*AIMAG(CDAT2(1))/NO2,KIND=SPC)
    CDAT1(2:)=CDAT1(2:)*CONJG(CDAT2(2:))/NO2
    CALL REALFT(CORR,-1,CDAT1)
END SUBROUTINE CORRELATION
