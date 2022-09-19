* COPYRIGHT (c) 1988 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date 14 June 2001
C  June 2001: threadsafe version of MC41
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC71AD(N,KASE,X,EST,W,IW,KEEP)
C
C      MC71A/AD ESTIMATES THE 1-NORM OF A SQUARE MATRIX A.
C      REVERSE COMMUNICATION IS USED FOR EVALUATING
C      MATRIX-VECTOR PRODUCTS.
C
C
C         N       INTEGER
C                 THE ORDER OF THE MATRIX.  N .GE. 1.
C
C         KASE    INTEGER
C                 SET INITIALLY TO ZERO . IF N .LE. 0 SET TO -1
C                 ON INTERMEDIATE RETURN
C                 = 1 OR 2.
C                ON FINAL RETURN
C                 =  0  ,IF SUCCESS
C                 = -1  ,IF N .LE.0
C
C         X       DOUBLE PRECISION ARRAY OF DIMENSION (N)
C                 IF 1-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      A*X,             IF KASE=1,
C                      TRANSPOSE(A)*X,  IF KASE=2,
C
C                 AND MC71 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C                 IF INFINITY-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      TRANSPOSE(A)*X,  IF KASE=1,
C                      A*X,             IF KASE=2,
C
C                 AND MC71 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C
C         EST     DOUBLE PRECISION
C                 CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C         W       DOUBLE PRECISION ARRAY OF DIMENSION (N)
C                 = A*V,   WHERE  EST = NORM(W)/NORM(V)
C                          (V  IS NOT RETURNED).
C         IW      INTEGER(N) USED AS WORKSPACE.
C
C         KEEP    INTEGER ARRAY LENGTH 5 USED TO PRESERVE PRIVATE
C                 DATA, JUMP, ITER, J AND JLAST BETWEEN CALLS,
C                 KEEP(5) IS SPARE.
C
C      REFERENCE
C      N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C      THE ONE-NORM OF A
C      REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C      TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
C      UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
C
C      SUBROUTINES AND FUNCTIONS
C
C
C      INTERNAL VARIABLES
C
C     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EST
      INTEGER KASE,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION W(*),X(*)
      INTEGER IW(*),KEEP(5)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALTSGN,TEMP
      INTEGER I,ITER,J,JLAST,JUMP
C     ..
C     .. External Functions ..
      INTEGER IDAMAX
      EXTERNAL IDAMAX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN,NINT,DBLE
C     ..
C     .. Executable Statements ..
C
      IF (N.LE.0) THEN
        KASE = -1
        RETURN

      END IF

      IF (KASE.EQ.0) THEN
        DO 10 I = 1,N
          X(I) = ONE/DBLE(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        KEEP(1) = JUMP
        KEEP(2) = 0
        KEEP(3) = 0
        KEEP(4) = 0
        RETURN

      END IF
C
      JUMP  = KEEP(1)
      ITER  = KEEP(2)
      J     = KEEP(3)
      JLAST = KEEP(4)
C
      GO TO (100,200,300,400,500) JUMP
C
C      ................ ENTRY   (JUMP = 1)
C
  100 CONTINUE
      IF (N.EQ.1) THEN
        W(1) = X(1)
        EST = ABS(W(1))
C         ... QUIT
        GO TO 510

      END IF
C
      DO 110 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  110 CONTINUE
      KASE = 2
      JUMP = 2
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 2)
C
  200 CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
C
C      MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
  220 CONTINUE
      DO 230 I = 1,N
        X(I) = ZERO
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 3)
C
  300 CONTINUE
C
C      COPY X INTO W
C
      DO 310 I = 1,N
        W(I) = X(I)
  310 CONTINUE
      DO 320 I = 1,N
        IF (NINT(SIGN(ONE,X(I))).NE.IW(I)) GO TO 330
  320 CONTINUE
C
C      REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 410
C
  330 CONTINUE
      DO 340 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 4)
C
  400 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF ((ABS(X(JLAST)).NE.ABS(X(J))) .AND. (ITER.LT.ITMAX)) THEN
        ITER = ITER + 1
        GO TO 220

      END IF
C
C      ITERATION COMPLETE.  FINAL STAGE.
C
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
        EST = EST + ABS(W(I))
  420 CONTINUE
C
      ALTSGN = ONE
      DO 430 I = 1,N
        X(I) = ALTSGN* (ONE+DBLE(I-1)/DBLE(N-1))
        ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 5)
C
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
        TEMP = TEMP + ABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/DBLE(3*N)
      IF (TEMP.GT.EST) THEN
C
C      COPY X INTO W
C
        DO 530 I = 1,N
          W(I) = X(I)
  530   CONTINUE
        EST = TEMP
      END IF
C
  510 KASE = 0
C
 1010 CONTINUE
      KEEP(1) = JUMP
      KEEP(2) = ITER
      KEEP(3) = J
      KEEP(4) = JLAST
      RETURN
C
      END
