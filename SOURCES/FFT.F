c  It runs and computes correctly although it is not in any 
c  sense sophisticated. However, it is free, and useful if 
c  you are in a hurry and don't want to think things out for yourself.
c   
c  Arthur Wouk (wouk@brl-vgr)
c  
c  *******************
c  Fast Fourier Transform
c  *******************
c  
	SUBROUTINE FFT (X, N, K)
C	FFT COMPUTES THE (FAST) FOURIER TRANSFORM OF THE VECTOR X
C	(A COMPLEX ARRAY OF DIMENSION N). SOURCE: Ferziger; Numerical
C       methods for engineering applications.
C
C	X = DATA TO BE TRANSFORMED; ON RETURN IT CONTAINS THE TRANSFORM.
C	N = SIZE OF VECTOR. MUST BE A POWER OF 2 (<32769).
C	K = 1 FOR FORWARD TRANSFORM.
C	K = -1 FOR INVERSE TRANSFORM.
C
	IMPLICIT INTEGER (A-Z)
	INTEGER SBY2,S
	REAL GAIN, PI2, ANG, RE, IM
	COMPLEX X(N), XTEMP, T, U(16), V, W
	LOGICAL NEW
	DATA	PI2,GAIN,NO,KO /6.283185307, 1., 0, 0/
C
C	TEST FIRST CALL?
C
	NEW = ( NO .NE. N)
	IF ( .NOT. NEW) GO TO 2
C
C	IF FIRST CALL COMPUTE LOG2 (N).
C
	L2N = 0
	NO = 1
    1	L2N = L2N + 1
	NO = NO + NO
	IF (NO .LT. N) GO TO 1
	GAIN = 1./N
	ANG = PI2*GAIN
	RE = COS (ANG)
	IM = SIN (ANG)
C
C	COMPUTE COMPLEX EXPONENTIALS IF NOT FIRST CALL
C
    2	IF (.NOT. NEW .AND. K*KO .GE. 1) GO TO 4
	U(1) = CMPLX (RE, -SIGN(IM, FLOAT(K)))
	DO 3 I = 2,L2N
  	   U(I) = U(I-1)*U(I-1)
    3	CONTINUE
	K0 = K
C
C	MAIN LOOP
C
    4	SBY2 = N
	DO 7 STAGE = 1,L2N
  	   V = U(STAGE)
  	   W = (1., 0.)
  	   S = SBY2
  	   SBY2 = S/2
  	   DO 6 L = 1,SBY2
		DO 5 I = 1,N,S
		   P = I + L- 1
		   Q = P + SBY2
		   T = X(P) + X(Q)
		   X(Q) = ( X(P) - X(Q))*W
		   X(P) =T
    5 		CONTINUE
		W = W*V
    6	   CONTINUE
    7	CONTINUE
C
C	REORDER THE ELEMENTS BY BIT REVERSAL
C
	DO 9 I = 1,N
	   INDEX = I-1
	   JNDEX = 0
	   DO 8 J = 1,L2N
		JNDEX = JNDEX+JNDEX
		ITEMP = INDEX/2
		IF (ITEMP+ITEMP .NE. INDEX) JNDEX = JNDEX + 1
		INDEX = ITEMP
    8	   CONTINUE
	   J = JNDEX + 1
	   IF (J .LT. I) GO TO 9
	   XTEMP = X(J)
	   X(J) = X(I)
	   X(I) = XTEMP
    9	CONTINUE
C
C	FORWARD TRANSFORM DONE
C
	IF (K .GT. 0) RETURN
C
C	INVERSE TRANSFORM
C
	DO 10 I = 1,N
	   X(I) = X(I)*GAIN
   10  	CONTINUE
	RETURN
	END
