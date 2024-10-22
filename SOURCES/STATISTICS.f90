!=============================================================================
! MODULE: Module Name
!
!> @author
!> Felix Dietzsch
!
! DESCRIPTION
!> Provides functions for statistical computations
!>
!=============================================================================
MODULE STATISTICS
    !INTERFACE FFT
       !MODULE PROCEDURE FFT1D
    !END INTERFACE
    INTERFACE MEAN
       MODULE PROCEDURE MEAN1D,MEAN3D,MEAN4D
    END INTERFACE
    INTERFACE CORREL
       MODULE PROCEDURE CORREL3D,XCORREL
    END INTERFACE CORREL
    INTERFACE SHIFT
        MODULE PROCEDURE FFTSHIFT
    END INTERFACE SHIFT
    INTERFACE SPECTR
        MODULE PROCEDURE SPEC3D,SPEC2D
    END INTERFACE
    INTERFACE DEVIAT
       MODULE PROCEDURE DEVIATION
    END INTERFACE DEVIAT
    CONTAINS
        !=============================================================================
        !> @author Felix Dietzsch
        !
        ! DESCRIPTION:
        !> @detail
        !> Computes the three dimensional energy spectrum as the fourier transformation\n
        !> of the auto-correlation function. In the first step all three velocity\n
        !> components are fourier transformed. The second step consists of multiplying\n
        !> each transform with its conjugate complex in order to get the elements of the\n
        !> main diagonal of the velocity spectrum tensor
        !> \f[\hat{R}_{ii}\left(\kappa,t\right)=
        !> \left<\hat{u}_i^{*}\left(\kappa,t\right)\,\hat{u}_i\left(\kappa,t\right)\right>=
        !> \Phi_{ii}\left(\kappa,t\right).\f]
        !> After the components of \f$\Phi\f$ have been determined an integration over\n
        !> spherical shells has to be computed in order to get the three dimensional\n
        !> velocity spectrum.
        !> \f[
        !> E(k)=\oint\frac{1}{2}\Phi_{ii}(\kappa)\,d\mathrm{S}(\kappa)
        !> \f]
        !> @brief
        !> Computes the two dimensional spectrum
        !
        ! REVISION HISTORY:
        ! 09 08 2011 Initial Version
        !
        !> @param[in] IN1 3D array of first velocity component
        !> @param[in] IN2 3D array of second velocity component
        !> @param[in] NX Number of nodes in x direction
        !> @param[in] NY Number of nodes in y direction
        !> @param[out] OUT Vector containing the spectrum for the input
        !> velocities
        !=============================================================================

        SUBROUTINE SPEC2D(IN1,IN2,NX,NY,OUT)
            USE NRTYPE
            USE INIT
            USE, INTRINSIC :: iso_c_binding
            IMPLICIT NONE
            include 'fftw3.f03'
            INTEGER(SP) :: NX,NY
            INTEGER(SP) :: I,J,KAPPA_POS
            REAL(DP) :: KAPPA,II,JJ,SCALING
            !INTEGER(SP) :: FFTW_DIRECTION !< Direction of the Fourier transform. Can be found in the  fftw3.f file, located in the include directory
            !INTEGER(SP) :: FFTW_ESTIMATE !< Perform FFT without any precomputation for an optimal algorithm. Can be found in the  fftw3.f file, located in the include directory

            TYPE(C_PTR) :: PLAN
            COMPLEX(DPC),DIMENSION(:,:),ALLOCATABLE :: WORKDATA1,WORKDATA2
            COMPLEX(DPC),DIMENSION(:,:),ALLOCATABLE :: PHI_X,PHI_Y
            COMPLEX(DPC),DIMENSION(:),INTENT(INOUT) :: OUT
            REAL(DP),DIMENSION(:,:,:) :: IN1,IN2
            COMPLEX(DPC),DIMENSION(NX,NY) :: TEMP1,TEMP2
            ALLOCATE(WORKDATA1(NX,NY))
            ALLOCATE(WORKDATA2(NX,NY))
            ALLOCATE(PHI_X(NX/2,NY/2))
            ALLOCATE(PHI_Y(NX/2,NY/2))
            TEMP1(:,:)=DCMPLX(IN1(:,:,1))
            TEMP2(:,:)=DCMPLX(IN2(:,:,1))
            OUT(:) = DCMPLX(0.D0,0.D0)
            !FFTW_DIRECTION = -1
            !FFTW_ESTIMATE = 64 
            !CALL FFTW_PLAN_DFT_2D(PLAN,NX,NY,TEMP1,WORKDATA1,FFTW_FORWARD,FFTW_ESTIMATE)
            plan = FFTW_PLAN_DFT_2D(NX,NY,TEMP1,WORKDATA1,FFTW_FORWARD,FFTW_ESTIMATE)
            CALL FFTW_EXECUTE_DFT(PLAN,TEMP1,WORKDATA1)
            CALL FFTW_DESTROY_PLAN(PLAN)
            plan = FFTW_PLAN_DFT_2D(NX,NY,TEMP2,WORKDATA2,FFTW_FORWARD,FFTW_ESTIMATE)
            CALL FFTW_EXECUTE_DFT(PLAN,TEMP1,WORKDATA2)
            CALL FFTW_DESTROY_PLAN(PLAN)
            SCALING = NX*NY
            PHI_X(:,:)= WORKDATA1(1:NX/2,1:NY/2)*CONJG(WORKDATA1(1:NX/2,1:NY/2))*1.D0/SCALING**2
            PHI_Y(:,:)= WORKDATA2(1:NX/2,1:NY/2)*CONJG(WORKDATA2(1:NX/2,1:NY/2))*1.D0/SCALING**2
            DEALLOCATE(WORKDATA1)
            DEALLOCATE(WORKDATA2)
            DO J=1,NY/2
                DO I=1,NX/2
                    II = I-1
                    JJ = J-1
                    KAPPA = DSQRT(II**2+JJ**2)
                    KAPPA_POS = INT(KAPPA+.5)
                    IF (KAPPA_POS.LE.SIZE(OUT)-1) THEN
                        !KAPPA = .5 *((REAL(PHI_X(I,J))**2 + AIMAG(PHI_X(I,J))**2) &
                                    !+(REAL(PHI_Y(I,J))**2 + AIMAG(PHI_Y(I,J))**2))
                        KAPPA = .5 *REAL(((PHI_X(I,J)+PHI_Y(I,J))))
                        OUT(KAPPA_POS) = OUT(KAPPA_POS) + 2.D0*KAPPA
                    END IF
                END DO
            END DO
            DEALLOCATE(PHI_X)
            DEALLOCATE(PHI_Y)
        END SUBROUTINE SPEC2D
        !=============================================================================
        !> @author Felix Dietzsch
        !
        ! DESCRIPTION:
        !> @detail
        !> Computes the three dimensional energy spectrum as the fourier transformation\n
        !> of the auto-correlation function. In the first step all three velocity\n
        !> components are fourier transformed. The second step consists of multiplying\n
        !> each transform with its conjugate complex in order to get the elements of the\n
        !> main diagonal of the velocity spectrum tensor
        !> \f[\hat{R}_{ii}\left(\kappa,t\right)=
        !> \left<\hat{u}_i^{*}\left(\kappa,t\right)\,\hat{u}_i\left(\kappa,t\right)\right>=
        !> \Phi_{ii}\left(\kappa,t\right).\f]
        !> After the components of \f$\Phi\f$ have been determined an integration over\n
        !> spherical shells has to be computed in order to get the three dimensional\n
        !> velocity spectrum.
        !> \f[
        !> E(k)=\oint\frac{1}{2}\Phi_{ii}(\kappa)\,d\mathrm{S}(\kappa)
        !> \f]
        !> @brief
        !> Computes the three dimensional spectrum
        !
        ! REVISION HISTORY:
        ! 25 07 2011 Initial Version
        !
        !> @param[in] IN1 3D array of first velocity component
        !> @param[in] IN2 3D array of second velocity component
        !> @param[in] IN3 3D array of third velocity component
        !> @param[in] NX Number of nodes in x direction
        !> @param[in] NY Number of nodes in y direction
        !> @param[in] NZ Number of nodes in z direction
        !> @param[out] Vector containing the spectrum for the input
        !> velocities
        !
        !> @todo
        !>  Implementation of the 3D spectrum computation\n
        !>  eMail from Michael Gauding from July 12th 2011\n
        !>  concerns the integration of spherical shells
        !=============================================================================
        SUBROUTINE SPEC3D(velo,out,kappa)
            USE NRTYPE
            USE INIT
            USE,INTRINSIC :: iso_c_binding
            IMPLICIT NONE
            include "fftw3.f03"
            TYPE(C_PTR) :: plan
            INTEGER(SP) :: nx,ny,nz
            INTEGER(DP) :: i,j,k
            INTEGER(DP) :: kappa_pos
            REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: kappa
            REAL(DP) :: scaling,kappa_abs
            REAL(DP),DIMENSION(:,:,:,:),INTENT(IN) :: velo
            REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: temp
            COMPLEX(DPC),DIMENSION(:,:,:),ALLOCATABLE :: workdata1,workdata2,workdata3
            COMPLEX(DPC),DIMENSION(:,:,:),ALLOCATABLE :: phi_x,phi_y,phi_z,phi
            ! the spectrm ('out') is 3D for reasons of simplicity during
            ! execution, since almost all routines rely on 3D data
            REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: out

            nx = SIZE(velo,1)
            ny = SIZE(velo,2)
            nz = SIZE(velo,3)

            ALLOCATE(workdata1(nx,ny,nz))
            ALLOCATE(workdata2(nx,ny,nz))
            ALLOCATE(workdata3(nx,ny,nz))
            ALLOCATE(temp(nx,ny,nz))
            ALLOCATE(phi_x(nx/2,ny/2,nz/2))
            ALLOCATE(phi_y(nx/2,ny/2,nz/2))
            ALLOCATE(phi_z(nx/2,ny/2,nz/2))
            ALLOCATE(phi(nx,ny,nz))
            !
            temp(:,:,:) = velo(:,:,:,1)
            plan = FFTW_PLAN_DFT_R2C_3D(nx,ny,nz,temp,workdata1,FFTW_ESTIMATE)
            CALL FFTW_EXECUTE_DFT_R2C(plan,temp,workdata1)
            CALL FFTW_DESTROY_PLAN(plan)
            temp(:,:,:) = velo(:,:,:,2)
            plan = FFTW_PLAN_DFT_R2C_3D(nx,ny,nz,temp,workdata2,FFTW_ESTIMATE)
            CALL FFTW_EXECUTE_DFT_R2C(plan,temp,workdata2)
            CALL FFTW_DESTROY_PLAN(plan)
            temp(:,:,:) = velo(:,:,:,3)
            plan = FFTW_PLAN_DFT_R2C_3D(nx,ny,nz,temp,workdata3,FFTW_ESTIMATE)
            CALL FFTW_EXECUTE_DFT_R2C(plan,temp,workdata3)
            CALL FFTW_DESTROY_PLAN(plan)
            scaling = nx*ny*nz ! comes from FFT
            ! the product of a*conj(a) is casted to real in order to avoid
            ! warnings during compilation
            phi_x(:,:,:)= (workdata1(nx/2,ny/2,nz/2)*CONJG(workdata1(nx/2,ny/2,nz/2)))/scaling/scaling
            phi_y(:,:,:)= (workdata2(nx/2,ny/2,nz/2)*CONJG(workdata2(nx/2,ny/2,nz/2)))/scaling/scaling
            phi_z(:,:,:)= (workdata3(nx/2,ny/2,nz/2)*CONJG(workdata3(nx/2,ny/2,nz/2)))/scaling/scaling
            !phi(:,:,:) = phi_x(:,:,:) + phi_y(:,:,:) + phi_z(:,:,:)
            DEALLOCATE(workdata1)
            DEALLOCATE(workdata2)
            DEALLOCATE(workdata3)
            DEALLOCATE(temp)
            PRINT*,"Compute spectrum"
            DO k=1,nz/2-1
                DO j=1,ny/2-1
                    DO i=1,nx/2-1
                        kappa_abs = SQRT(REAL(k*k)+REAL(j*j)+REAL(i*i))
                        kappa_pos = INT(kappa_abs)
                        !out(kappa_pos,1,1) = out(kappa_pos,1,1) + phi(i,j,k)*kappa_abs*kappa_abs
                        out(kappa_pos,1,1) = out(kappa_pos,1,1) + &
                                            (REAL(phi_x(i,j,k))*REAL(phi_x(i,j,k))+AIMAG(phi_x(i,j,k))*AIMAG(phi_x(i,j,k)) &
                                            +REAL(phi_y(i,j,k))*REAL(phi_y(i,j,k))+AIMAG(phi_y(i,j,k))*AIMAG(phi_y(i,j,k)) &
                                            +REAL(phi_z(i,j,k))*REAL(phi_z(i,j,k))+AIMAG(phi_z(i,j,k))*AIMAG(phi_z(i,j,k)))
                        kappa(kappa_pos,1,1) = kappa_pos
                    ENDDO
                ENDDO
            ENDDO
            out(:,1,1) = 2*PI_D*out(:,1,1)
            DEALLOCATE(phi_x)
            DEALLOCATE(phi_y)
            DEALLOCATE(phi_z)

            !DO J=1,NY/2
                !DO I=1,NX/2
                    !DO K=1,NZ/2
            !            II = I
            !            JJ = J
            !            KK = K
            !            KAPPA = DSQRT(II**2+JJ**2+KK**2)
            !            KAPPA_POS = INT(KAPPA+.5)
            !            IF (KAPPA_POS.LE.SIZE(OUT)-1) THEN
            !                !KAPPA = 0.5d0 *((REAL(PHI_X(I,J,K))**2 + AIMAG(PHI_X(I,J,K))**2) &
            !                              !+(REAL(PHI_Y(I,J,K))**2 + AIMAG(PHI_Y(I,J,K))**2) &
            !                              !+(REAL(PHI_Z(I,J,K))**2 + AIMAG(PHI_Z(I,J,K))**2))
            !            KAPPA = .5 *(PHI_X(I,J,K)+PHI_Y(I,J,K)+PHI_Z(I,J,K))
            !            OUT(KAPPA_POS) = OUT(KAPPA_POS) + 2.D0*KAPPA
            !            END IF
            !        END DO
            !    END DO
            !END DO
            !OUT(:)=OUT(:)*1.D0
        END SUBROUTINE SPEC3D
        !=============================================================================
        !> @author Felix Dietzsch
        !
        ! DESCRIPTION:
        !> @detail
        !> Computes the standard deviation of the 3D input array. It is computed
        !>according to
        !> \f[\sigma = \sqrt{\frac{1}{N} \sum_{i=1}^N (x_i - \mu)^2}, {\rm \ \ where\
        !> \ } \mu = \frac{1}{N} \sum_{i=1}^N x_i.
        !> \f]
        !> @brief
        !> Computes the standard deviation of the 3D input array.
        ! REVISION HISTORY:
        ! 25 07 2011 Initial Version
        !
        !> @param[in] VELO 3D array of which the deviation is to be computed
        !> @param[in] NX Number of nodes in x direction
        !> @param[in] NY Number of nodes in y direction
        !> @param[in] NZ Number of nodes in z direction
        !> @param[out] SIGMA Standard deviation of the input array
        !> velocities
        !=============================================================================
        SUBROUTINE DEVIATION(VELO,NX,NY,NZ,SIGMA)
            USE INIT
            USE NRTYPE
            INTEGER(SP) :: NX,NY,I,J,K
            REAL(DP) :: SIGMA,MEAN
            REAL(DP),DIMENSION(:,:,:) :: VELO
            !OPEN(1,FILE='INIT',STATUS='UNKNOWN',FORM='FORMATTED')
            !READ(1,NML=NML_INIT)
            !CLOSE(1)
            SIGMA=0
            CALL MEAN3D(VELO,NX,NY,NZ,MEAN)
            DO I=1,NX
                DO J=1,NY
                    IF (NZ.GT.1) THEN
                        DO K=1,NZ-1
                            SIGMA=SIGMA+(VELO(I,J,K)-MEAN)**2
                        END DO
                    ELSE
                        DO K=1,NZ
                            SIGMA=SIGMA+(VELO(I,J,K)-MEAN)**2
                        END DO
                    END IF
                END DO
            END DO
            SIGMA=SQRT(1.0/(NX*NY*NZ)*SIGMA)
        END SUBROUTINE DEVIATION
        !=============================================================================
        !> @author Felix Dietzsch
        !
        ! DESCRIPTION:
        !> @detail
        !> Shifts the input array towards the zero frequencies.
        !> \image html fftshift.jpg "Illustration of the fftshift routine"
        !> @brief
        !> Shifts the input array towards the zero frequencies.
        ! REVISION HISTORY:
        ! 25 07 2011 Initial Version
        !
        !> @param[in] IN 2D array of shifted frequencies 
        !> @param[in] NX Number of nodes in y direction
        !> @param[in] NY Number of nodes in z direction
        !> @param[out] IN 2D array of shifted frequencies 
        !=============================================================================
        SUBROUTINE FFTSHIFT(IN,NX,NY)
            USE NRTYPE
            USE INIT
            IMPLICIT NONE
            INTEGER(SP) :: NX,NY
            COMPLEX(DPC),DIMENSION(NX,NY),INTENT(INOUT) :: IN
            !COMPLEX(DPC),DIMENSION(2*NX,NX/2) :: TEMP
            COMPLEX(DPC),DIMENSION(NX/2,NY/2) :: SUB1
            COMPLEX(DPC),DIMENSION(NX/2,NY/2) :: SUB2
            COMPLEX(DPC),DIMENSION(NX/2,NY/2) :: SUB3
            COMPLEX(DPC),DIMENSION(NX/2,NY/2) :: SUB4
            SUB1(1:NX/2,1:NY/2) = IN(1:NX/2,1:NY/2)
            SUB2(1:NX/2,1:NY/2) = IN(NX/2+1:NX,1:NY/2)
            SUB3(1:NX/2,1:NY/2) = IN(NX/2+1:NX,NY/2+1:NY)
            SUB4(1:NX/2,1:NY/2) = IN(1:NX/2,NY/2+1:NY)
            IN(1:NX/2,1:NY/2) = SUB3(:,:)
            IN(NX/2+1:NX,NY/2+1:NY)=SUB1(:,:)
            IN(NX/2+1:NX,1:NY/2) = SUB4(:,:)
            IN(1:NX/2,NY/2+1:NY) = SUB2(:,:)
        END SUBROUTINE FFTSHIFT
        !=============================================================================
        !> @author
        !> Felix Dietzsch
        !
        ! DESCRIPTION:
        !> @detail
        !> Computes the auto-correlation coefficients of the 3D input array. The
        !> general definition of the auto-correlation (often exressed as the
        !> covariance) function is
        !> \f[\mathrm{cov}(f,f)=\left(f\star f\right)(r)=\int_0^{\infty}f(x)\,f(x+r)\,dx
        !> \f]
        !> One important feature of the correlation function is that is satisfies
        !> \f[\mathcal{F}\{f\star f\}=(\mathcal{F}\{f\})^*\cdot\mathcal{F}\{f\}.
        !> \f]
        !> In the code this allows for fast computations of the correlation
        !> function using the FFT approach.
        !> In order to get the correlation coefficients, the previously mentioned
        !> function has to be devided by standard deviations of the input
        !> arrays.
        !> \f[r=\frac{\mathrm{cov}(f,f)}{\sigma_f\,\sigma_f}
        !> \f]
        !> @brief
        !> Computes the auto-correlation coefficients of the 3D input array 
        !
        ! REVISION HISTORY:
        ! 25 07 2011 Initial Version
        !
        !> @param[in] IN 3D array for the computation of the 3D auto-correlation
        !> coefficients
        !> @param[in] NX Number of nodes in x direction
        !> @param[in] NY Number of nodes in y direction
        !> @param[in] NZ Number of nodes in z direction
        !> @param[out] OUT 3D array of auto-correlation coeffcients 
        !=============================================================================
        SUBROUTINE CORREL3D(IN,OUT)
            USE NRTYPE
            USE INIT
            USE FGSL
            USE,INTRINSIC :: iso_c_binding
            IMPLICIT NONE
            include 'fftw3.f03'
            INTEGER(SP) :: nx,ny,nz
            TYPE(C_PTR) :: plan
            REAL(fgsl_double) :: sigma
            INTEGER(fgsl_size_t) :: dimen
            REAL(fgsl_double),DIMENSION(:),ALLOCATABLE :: velo
            REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: temp
            REAL(DP),DIMENSION(:,:,:),INTENT(IN) :: in
            REAL(DP),DIMENSION(:,:,:),INTENT(INOUT) :: out
            COMPLEX(DPC),DIMENSION(:,:,:),ALLOCATABLE :: workdata
            COMPLEX(DPC),DIMENSION(:,:,:),ALLOCATABLE :: workdata1
            
            nx = SIZE(in,1)
            ny = SIZE(in,2)
            nz = SIZE(in,3)
            dimen = nx*ny*nz
            ALLOCATE(velo(dimen))
            ALLOCATE(workdata(nx,ny,nz))
            ALLOCATE(workdata1(nx,ny,nz))
            ALLOCATE(temp(nx,ny,nz))

            velo = RESHAPE(in,(/nx*ny*nz/))
            velo(:) = velo(:) - fgsl_stats_mean(velo,1_fgsl_size_t,dimen)

            sigma = fgsl_stats_variance(velo,1_fgsl_size_t,dimen)
            sigma = DSQRT(sigma)
            !PRINT*,"SIGMA=",sigma
            !PRINT*,"MEAN=",fgsl_stats_mean(velo,1_fgsl_size_t,dimen)

            temp(:,:,:) = in(:,:,:) - fgsl_stats_mean(velo,1_fgsl_size_t,dimen)

            !!!An r2c transform produces the same output as a FFTW_FORWARD complex
            !!!DFT of the same input

            plan = FFTW_PLAN_DFT_R2C_3D(nx,ny,nz,temp,workdata,FFTW_ESTIMATE)
            CALL FFTW_EXECUTE_DFT_R2C(plan,temp,workdata)
            !compute autocorrelation in Fourier space
            workdata1(:,:,:) = workdata(:,:,:)*CONJG(workdata(:,:,:))/dimen/dimen/sigma/sigma
            !PRINT*,MAXVAL(REAL(workdata1))
            !PRINT*,MAXVAL(AIMAG(workdata1))
            CALL FFTW_DESTROY_PLAN(plan)
            ! Perform backward transformation to get real valued correlation
            plan = FFTW_PLAN_DFT_C2R_3D(nx,ny,nz,workdata1,out,FFTW_ESTIMATE)
            CALL FFTW_EXECUTE_DFT_C2R(plan,workdata1,out)
            DEALLOCATE(velo)
            DEALLOCATE(workdata)
            DEALLOCATE(workdata1)
            DEALLOCATE(temp)
            PRINT*,'DONE executing first transform'
        END SUBROUTINE CORREL3D
        !=============================================================================
        !> @author Felix Dietzsch
        !
        ! DESCRIPTION:
        !> @detail
        !> Computes the cross-correlation coefficients of the 3D input arrays. The
        !> general definition of the cross-correlation (often exressed as the
        !> covariance) function is
        !> \f[\mathrm{cov}(f,g)=\left(f\star g\right)(r)=\int_0^{\infty}f(x)\,g(x+r)\,dx
        !> \f]
        !> One important feature of the correlation function is that is satisfies
        !> \f[\mathcal{F}\{f\star g\}=(\mathcal{F}\{f\})^*\cdot\mathcal{F}\{g\}.
        !> \f]
        !> In the code this allows for fast computations of the correlation
        !> function using the FFT approach.
        !> In order to get the correlation coefficients, the previously mentioned
        !> function has to be devided by standard deviations of the input
        !> arrays.
        !> \f[r=\frac{\mathrm{cov}(f,g)}{\sigma_f\,\sigma_g}
        !> \f]
        !> Computes the cross-correlation of the 3D input array
        !> @brief
        !>
        !
        ! REVISION HISTORY:
        ! 25 07 2011 Initial Version
        ! 05 12 2011 Made VELO allocatable --> fixed error during compilation
        !
        !> @param[in] IN1 3D array for the computation of the 3D cross-correlation
        !> coefficients
        !> @param[in] IN2 3D array for the computation of the 3D cross-correlation
        !> coefficients
        !> @param[in] NX Number of nodes in x direction
        !> @param[in] NY Number of nodes in y direction
        !> @param[in] NZ Number of nodes in z direction
        !> @param[out] OUT 3D array of cross-correlation coeffcients 
        !=============================================================================
        SUBROUTINE XCORREL(IN1,IN2,NX,NY,NZ,OUT)
            USE NRTYPE
            USE INIT
            USE, INTRINSIC :: iso_c_binding
            IMPLICIT NONE
            include 'fftw3.f03'
            INTEGER(SP) :: NX,NY,NZ
            !INTEGER(SP) :: FFTW_DIRECTION,FFTW_ESTIMATE
            TYPE(C_PTR) :: PLAN
            REAL(DP) :: SIGMAX,SIGMAY
            REAL(SP) :: N
            REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: VELO
            COMPLEX(DPC),DIMENSION(:,:,:),ALLOCATABLE :: WORKDATA1
            COMPLEX(DPC),DIMENSION(:,:,:),ALLOCATABLE :: WORKDATA2
            COMPLEX(DPC),DIMENSION(:,:,:),INTENT(INOUT) :: IN1
            COMPLEX(DPC),DIMENSION(:,:,:),INTENT(INOUT) :: IN2
            COMPLEX(DPC),DIMENSION(:,:,:),INTENT(INOUT) :: OUT
            ALLOCATE(VELO(NX,NY,NZ))
            VELO(:,:,:) = REAL(IN1(:,:,:))
            CALL DEVIATION(VELO,NX,NY,NZ,SIGMAX)
            VELO(:,:,:) = REAL(IN2(:,:,:))
            CALL DEVIATION(VELO,NX,NY,NZ,SIGMAY)
            ALLOCATE(WORKDATA1(NX,NY,NZ))
            ALLOCATE(WORKDATA2(NX,NY,NZ))
            !FFTW_DIRECTION = -1 ! forward transform
            ! internal parameter of the fftw routines. Can be found in the
            ! fftw3.f file, located in the include directory
            !FFTW_ESTIMATE = 64 
            plan = FFTW_PLAN_DFT_3D(NX,NY,NZ,IN1,WORKDATA1,FFTW_FORWARD,FFTW_ESTIMATE)
            CALL FFTW_EXECUTE_DFT(PLAN,IN1,WORKDATA1)
            plan = FFTW_PLAN_DFT_3D(NX,NY,NZ,IN2,WORKDATA2,FFTW_FORWARD,FFTW_ESTIMATE)
            CALL FFTW_EXECUTE_DFT(PLAN,IN1,WORKDATA2)
            N=NX*NY*NZ
            OUT(:,:,:) = WORKDATA1(:,:,:)*CONJG(WORKDATA2(:,:,:))/SIGMAX/SIGMAY/N**2
            CALL FFTW_DESTROY_PLAN(PLAN)
            !FFTW_DIRECTION=1 ! backward transform
            plan = FFTW_PLAN_DFT_3D(NX,NY,NZ,OUT,OUT,FFTW_BACKWARD,FFTW_ESTIMATE)
            CALL FFTW_EXECUTE_DFT(PLAN,OUT,OUT)
            DEALLOCATE(WORKDATA1)
            DEALLOCATE(WORKDATA2)
        END SUBROUTINE XCORREL
        SUBROUTINE MEAN1D()
            USE NRTYPE
            CONTINUE
        END SUBROUTINE MEAN1D
        !=============================================================================
        !> @author Felix Dietzsch
        !
        ! DESCRIPTION:
        !> Shifts the input array towards the zero frequencies 
        !
        ! REVISION HISTORY:
        ! 25 07 2011 Initial Version
        !
        !> @param[in] M 3D array for the computation of the mean value
        !> @param[in] NX Number of nodes in x direction
        !> @param[in] NY Number of nodes in y direction
        !> @param[in] NZ Number of nodes in z direction
        !> @param[out] MEAN Computed mean value 
        !=============================================================================
        SUBROUTINE MEAN3D(M,NX,NY,NZ,MEAN)
            USE NRTYPE
            INTEGER(SP) :: NX,NY
            REAL(DP) :: MEAN
            REAL(DP),DIMENSION(:,:,:) :: M
            MEAN = 0.0
            MEAN = SUM(M)
            MEAN = MEAN/(NX*NY*NZ)
        END SUBROUTINE MEAN3D
        SUBROUTINE MEAN4D(M,NX,NY,NZ,MEAN)
            USE NRTYPE
            USE INIT
            INTEGER(SP) :: NX,NY
            REAL(DP) :: MEAN
            REAL(DP),DIMENSION(:,:,:,:) :: M
            MEAN = 0.0
            MEAN=SUM(M)
            MEAN = MEAN/(DIM*NX*NY*NZ)
        END SUBROUTINE MEAN4D
END MODULE STATISTICS
