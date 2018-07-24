MODULE separable_leastsq
IMPLICIT NONE

! This is a conversion of the VARPRO package to F90 style.

! Code converted using TO_F90 by Alan Miller
! Date: 1999-03-15  Time: 11:15:27
! Latest version - 11 May 2001
! Changes from original:
! 1. Many variables needed to be saved in routine VPDPA.
! 2. A couple of arguments previously declared with INTENT(OUT) should
!    have been INTENT(IN OUT).

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(30,60)

CONTAINS


SUBROUTINE varpro (l, nl, n, nmax, lpp2, iv, t, y, w, ada, a,  &
                   iprint, alf, beta, ierr)

!        GIVEN A SET OF N OBSERVATIONS, CONSISTING OF VALUES Y(1), Y(2), ...,
!     Y(N) OF A DEPENDENT VARIABLE Y, WHERE Y(I) CORRESPONDS TO THE IV
!     INDEPENDENT VARIABLE(S) T(I,1), T(I,2), ..., T(I,IV), VARPRO ATTEMPTS TO
!     COMPUTE A WEIGHTED LEAST SQUARES FIT TO A FUNCTION ETA (THE 'MODEL')
!     WHICH IS A LINEAR COMBINATION

!                           L
!    ETA(ALF, BETA; T)  =  SUM  BETA  * PHI (ALF; T) + PHI   (ALF; T)
!                          J=1      J      J              L+1

!     OF NONLINEAR FUNCTIONS PHI(J) (E.G., A SUM OF EXPONENTIALS AND/
!     OR GAUSSIANS).  THAT IS, DETERMINE THE LINEAR PARAMETERS
!     BETA(J) AND THE VECTOR OF NONLINEAR PARAMETERS ALF BY MINIMIZING

!                      2     N                                 2
!        NORM(RESIDUAL)  =  SUM  W  * (Y  - ETA(ALF, BETA; T )) .
!                           I=1   I     I                   I

!     THE (L+1)-ST TERM IS OPTIONAL, AND IS USED WHEN IT IS DESIRED
!     TO FIX ONE OR MORE OF THE BETA'S (RATHER THAN LET THEM BE
!     DETERMINED).  VARPRO REQUIRES FIRST DERIVATIVES OF THE PHI'S.

!                             NOTES:

!     A)  THE ABOVE PROBLEM IS ALSO REFERRED TO AS 'MULTIPLE NONLINEAR
!     REGRESSION'.  FOR USE IN STATISTICAL ESTIMATION, VARPRO RETURNS THE
!     RESIDUALS, THE COVARIANCE MATRIX OF THE LINEAR AND NONLINEAR PARAMETERS,
!     AND THE ESTIMATED VARIANCE OF THE OBSERVATIONS.

!     B) AN ETA OF THE ABOVE FORM IS CALLED 'SEPARABLE'.  THE CASE OF A
!     NONSEPARABLE ETA CAN BE HANDLED BY SETTING L = 0 AND USING PHI(L+1).

!     C) VARPRO MAY ALSO BE USED TO SOLVE LINEAR LEAST SQUARES PROBLEMS
!     (IN THAT CASE NO ITERATIONS ARE PERFORMED).  SET NL = 0.

!     D)  THE MAIN ADVANTAGE OF VARPRO OVER OTHER LEAST SQUARES PROGRAMS IS
!     THAT NO INITIAL GUESSES ARE NEEDED FOR THE LINEAR PARAMETERS.  NOT ONLY
!     DOES THIS MAKE IT EASIER TO USE, BUT IT OFTEN LEADS TO FASTER CONVERGENCE.


!  DESCRIPTION OF PARAMETERS

!     L       NUMBER OF LINEAR PARAMETERS BETA (MUST BE >= 0).
!     NL      NUMBER OF NONLINEAR PARAMETERS ALF (MUST BE >= 0).
!     N       NUMBER OF OBSERVATIONS.  N MUST BE GREATER THAN L + NL
!             (I.E., THE NUMBER OF OBSERVATIONS MUST EXCEED THE
!             NUMBER OF PARAMETERS).
!     IV      NUMBER OF INDEPENDENT VARIABLES T.
!     T       REAL N BY IV MATRIX OF INDEPENDENT VARIABLES.  T(I, J) CONTAINS
!             THE VALUE OF THE I-TH OBSERVATION OF THE J-TH INDEPENDENT
!             VARIABLE.
!     Y       N-VECTOR OF OBSERVATIONS, ONE FOR EACH ROW OF T.
!     W       N-VECTOR OF NONNEGATIVE WEIGHTS.  SHOULD BE SET TO 1'S IF WEIGHTS
!             ARE NOT DESIRED.  IF VARIANCES OF THE INDIVIDUAL OBSERVATIONS
!             ARE KNOWN, W(I) SHOULD BE SET TO 1./VARIANCE(I).
!     INC     NL X (L+1) INTEGER INCIDENCE MATRIX.  INC(K, J) = 1 IF
!             NON-LINEAR PARAMETER ALF(K) APPEARS IN THE J-TH
!             FUNCTION PHI(J).  (THE PROGRAM SETS ALL OTHER INC(K, J)
!             TO ZERO.)  IF PHI(L+1) IS INCLUDED IN THE MODEL,
!             THE APPROPRIATE ELEMENTS OF THE (L+1)-ST COLUMN SHOULD
!             BE SET TO 1'S.  INC IS NOT NEEDED WHEN L = 0 OR NL = 0.
!             CAUTION:  THE DECLARED ROW DIMENSION OF INC (IN ADA)
!             MUST CURRENTLY BE SET TO 12.  SEE 'RESTRICTIONS' BELOW.
!     NMAX    THE DECLARED ROW DIMENSION OF THE MATRICES A AND T.
!             IT MUST BE AT LEAST MAX(N, 2*NL+3).
!     LPP2    L+P+2, WHERE P IS THE NUMBER OF ONES IN THE MATRIX INC.
!             THE DECLARED COLUMN DIMENSION OF A MUST BE AT LEAST
!             LPP2.  (IF L = 0, SET LPP2 = NL+2. IF NL = 0, SET LPP2 L+2.)
!     A       REAL MATRIX OF SIZE MAX(N, 2*NL+3) BY L+P+2.  ON INPUT
!             IT CONTAINS THE PHI(J)'S AND THEIR DERIVATIVES (SEE
!             BELOW).  ON OUTPUT, THE FIRST L+NL ROWS AND COLUMNS OF
!             A WILL CONTAIN AN APPROXIMATION TO THE (WEIGHTED)
!             COVARIANCE MATRIX AT THE SOLUTION (THE FIRST L ROWS
!             CORRESPOND TO THE LINEAR PARAMETERS, THE LAST NL TO THE
!             NONLINEAR ONES), COLUMN L+NL+1 WILL CONTAIN THE
!             WEIGHTED RESIDUALS (Y - ETA), A(1, L+NL+2) WILL CONTAIN
!             THE (EUCLIDEAN) NORM OF THE WEIGHTED RESIDUAL, AND
!             A(2, L+NL+2) WILL CONTAIN AN ESTIMATE OF THE (WEIGHTED)
!             VARIANCE OF THE OBSERVATIONS, NORM(RESIDUAL)**2 / (N - L - NL).
!     IPRINT  INPUT INTEGER CONTROLLING PRINTED OUTPUT.  IF IPRINT IS
!             POSITIVE, THE NONLINEAR PARAMETERS, THE NORM OF THE
!             RESIDUAL, AND THE MARQUARDT PARAMETER WILL BE OUTPUT
!             EVERY IPRINT-TH ITERATION (AND INITIALLY, AND AT THE
!             FINAL ITERATION).  THE LINEAR PARAMETERS WILL BE
!             PRINTED AT THE FINAL ITERATION.  ANY ERROR MESSAGES
!             WILL ALSO BE PRINTED.  (IPRINT = 1 IS RECOMMENDED AT
!             FIRST.) IF IPRINT = 0, ONLY THE FINAL QUANTITIES WILL
!             BE PRINTED, AS WELL AS ANY ERROR MESSAGES.  IF IPRINT =
!             -1, NO PRINTING WILL BE DONE.  THE USER IS THEN
!             RESPONSIBLE FOR CHECKING THE PARAMETER IERR FOR ERRORS.
!     ALF     NL-VECTOR OF ESTIMATES OF NONLINEAR PARAMETERS
!             (INPUT).  ON OUTPUT IT WILL CONTAIN OPTIMAL VALUES OF
!             THE NONLINEAR PARAMETERS.
!     BETA    L-VECTOR OF LINEAR PARAMETERS (OUTPUT ONLY).
!     IERR    INTEGER ERROR FLAG (OUTPUT):
!             > 0 - SUCCESSFUL CONVERGENCE, IERR IS THE NUMBER OF ITERATIONS
!                 TAKEN.
!             -1  TERMINATED FOR TOO MANY ITERATIONS.
!             -2  TERMINATED FOR ILL-CONDITIONING (MARQUARDT
!                 PARAMETER TOO LARGE.)  ALSO SEE IERR = -8 BELOW.
!             -4  INPUT ERROR IN PARAMETER N, L, NL, LPP2, OR NMAX.
!             -5  INC MATRIX IMPROPERLY SPECIFIED, OR P DISAGREES WITH LPP2.
!             -6  A WEIGHT WAS NEGATIVE.
!             -7  'CONSTANT' COLUMN WAS COMPUTED MORE THAN ONCE.
!             -8  CATASTROPHIC FAILURE - A COLUMN OF THE A MATRIX HAS
!                 BECOME ZERO.  SEE 'CONVERGENCE FAILURES' BELOW.

!             (IF IERR .LE. -4, THE LINEAR PARAMETERS, COVARIANCE
!             MATRIX, ETC. ARE NOT RETURNED.)

!  SUBROUTINES REQUIRED

!        NINE SUBROUTINES, VPDPA, VPFAC1, VPFAC2, VPBSOL, VPPOST, VPCOV,
!     VPNORM, VPINIT, AND VPERR ARE PROVIDED.  IN ADDITION, THE USER MUST
!     PROVIDE A SUBROUTINE (CORRESPONDING TO THE ARGUMENT ADA) WHICH,
!     GIVEN ALF, WILL EVALUATE THE FUNCTIONS PHI(J) AND THEIR PARTIAL
!     DERIVATIVES D PHI(J)/D ALF(K), AT THE SAMPLE POINTS T(I).
!     ITS CALLING SEQUENCE IS

!     SUBROUTINE ADA (L+1, NL, N, NMAX, LPP2, IV, A, INC, T, ALF, ISEL)

!        THE USER SHOULD MODIFY THE EXAMPLE SUBROUTINE 'ADA' (GIVEN
!     ELSEWHERE) FOR HIS OWN FUNCTIONS.

!        THE VECTOR SAMPLED FUNCTIONS PHI(J) SHOULD BE STORED IN THE
!     FIRST N ROWS AND FIRST L+1 COLUMNS OF THE MATRIX A, I.E.,
!     A(I, J) SHOULD CONTAIN PHI(J, ALF; T(I,1), T(I,2), ...,
!     T(I,IV)), I = 1, ..., N; J = 1, ..., L (OR L+1).  THE (L+1)-ST
!     COLUMN OF A CONTAINS PHI(L+1) IF PHI(L+1) IS IN THE MODEL,
!     OTHERWISE IT IS RESERVED FOR WORKSPACE.  THE 'CONSTANT' FUNC-
!     TIONS (THESE ARE FUNCTIONS PHI(J) WHICH DO NOT DEPEND UPON ANY
!     NONLINEAR PARAMETERS ALF, E.G., T(I)**J) (IF ANY) MUST APPEAR
!     FIRST, STARTING IN COLUMN 1.  THE COLUMN N-VECTORS OF NONZERO
!     PARTIAL DERIVATIVES D PHI(J) / D ALF(K) SHOULD BE STORED
!     SEQUENTIALLY IN THE MATRIX A IN COLUMNS L+2 THROUGH L+P+1.
!     THE ORDER IS

!       D PHI(1)  D PHI(2)        D PHI(L)  D PHI(L+1)  D PHI(1)
!       --------, --------, ...,  --------, ----------, --------,
!       D ALF(1)  D ALF(1)        D ALF(1)   D ALF(1)   D ALF(2)

!       D PHI(2)       D PHI(L+1)       D PHI(1)        D PHI(L+1)
!       --------, ..., ----------, ..., ---------, ..., ----------,
!       D ALF(2)        D ALF(2)        D ALF(NL)       D ALF(NL)

!     OMITTING COLUMNS OF DERIVATIVES WHICH ARE ZERO, AND OMITTING
!     PHI(L+1) COLUMNS IF PHI(L+1) IS NOT IN THE MODEL.  NOTE THAT
!     THE LINEAR PARAMETERS BETA ARE NOT USED IN THE MATRIX A.
!     COLUMN L+P+2 IS RESERVED FOR WORKSPACE.

!     THE CODING OF ADA SHOULD BE ARRANGED SO THAT:

!     ISEL = 1  (WHICH OCCURS THE FIRST TIME ADA IS CALLED) MEANS:
!               A.  FILL IN THE INCIDENCE MATRIX INC
!               B.  STORE ANY CONSTANT PHI'S IN A.
!               C.  COMPUTE NONCONSTANT PHI'S AND PARTIAL DERIVATIVES.
!          = 2  MEANS COMPUTE ONLY THE NONCONSTANT FUNCTIONS PHI
!          = 3  MEANS COMPUTE ONLY THE DERIVATIVES

!     (WHEN THE PROBLEM IS LINEAR (NL = 0) ONLY ISEL = 1 IS USED, AND
!     DERIVATIVES ARE NOT NEEDED.)

!  RESTRICTIONS

!        THE SUBROUTINES VPDPA, VPINIT (AND ADA) CONTAIN THE LOCALLY
!     DIMENSIONED MATRIX INC, WHOSE DIMENSIONS ARE CURRENTLY SET FOR
!     MAXIMA OF L+1 = 8, NL = 12.  THEY MUST BE CHANGED FOR LARGER
!     PROBLEMS.  DATA PLACED IN ARRAY A IS OVERWRITTEN ('DESTROYED').
!     DATA PLACED IN ARRAYS T, Y AND INC IS LEFT INTACT.  THE PROGRAM
!     RUNS IN WATFIV, EXCEPT WHEN L = 0 OR NL = 0.

!        IT IS ASSUMED THAT THE MATRIX PHI(J, ALF; T(I)) HAS FULL
!     COLUMN RANK.  THIS MEANS THAT THE FIRST L COLUMNS OF THE MATRIX
!     A MUST BE LINEARLY INDEPENDENT.

!        OPTIONAL NOTE:  AS WILL BE NOTED FROM THE SAMPLE SUBPROGRAM ADA,
!     THE DERIVATIVES D PHI(J)/D ALF(K) (ISEL = 3) MUST BE COMPUTED
!     INDEPENDENTLY OF THE FUNCTIONS PHI(J) (ISEL = 2), SINCE THE FUNCTION
!     VALUES ARE OVERWRITTEN AFTER ADA IS CALLED WITH ISEL = 2.
!     THIS IS DONE TO MINIMIZE STORAGE, AT THE POSSIBLE EXPENSE OF SOME
!     RECOMPUTATION (SINCE THE FUNCTIONS AND DERIVATIVES FREQUENTLY HAVE SOME
!     COMMON SUBEXPRESSIONS).  TO REDUCE THE AMOUNT OF COMPUTATION AT THE
!     EXPENSE OF SOME STORAGE, CREATE A MATRIX B OF DIMENSION NMAX BY L+1
!     IN ADA, AND AFTER THE COMPUTATION OF THE PHI'S (ISEL = 2), COPY THE
!     VALUES INTO B.  THESE VALUES CAN THEN BE USED TO CALCULATE THE
!     DERIVATIVES (ISEL = 3).  (THIS MAKES USE OF THE FACT THAT WHEN A CALL TO
!     ADA WITH ISEL = 3 FOLLOWS A CALL WITH ISEL = 2, THE ALFS ARE THE SAME.)

!        TO CONVERT TO OTHER MACHINES, CHANGE THE OUTPUT UNIT IN THE
!     DATA STATEMENTS IN VARPRO, VPDPA, VPPOST, AND VPERR.  THE
!     PROGRAM HAS BEEN CHECKED FOR PORTABILITY BY THE BELL LABS PFORT
!     VERIFIER.  FOR MACHINES WITHOUT DOUBLE PRECISION HARDWARE, IT
!     MAY BE DESIRABLE TO CONVERT TO SINGLE PRECISION.  THIS CAN BE
!     DONE BY CHANGING (A) THE DECLARATIONS 'DOUBLE PRECISION' TO
!     'REAL', (B) THE PATTERN '.D' TO '.E' IN THE 'DATA' STATEMENT IN
!     VARPRO, (C) DSIGN, SQRT AND ABS TO SIGN, SQRT AND ABS,
!     RESPECTIVELY, AND (D) DEXP TO EXP IN THE SAMPLE PROGRAMS ONLY.

!  NOTE ON INTERPRETATION OF COVARIANCE MATRIX

!        FOR USE IN STATISTICAL ESTIMATION (MULTIPLE NONLINEAR REGRESSION)
!     VARPRO RETURNS THE COVARIANCE MATRIX OF THE LINEAR AND NONLINEAR
!     PARAMETERS.  THIS MATRIX WILL BE USEFUL ONLY IF THE USUAL STATISTICAL
!     ASSUMPTIONS HOLD:  AFTER WEIGHTING, THE ERRORS IN THE OBSERVATIONS ARE
!     INDEPENDENT AND NORMALLY DISTRIBUTED, WITH MEAN ZERO AND THE SAME
!     VARIANCE.  IF THE ERRORS DO NOT HAVE MEAN ZERO (OR ARE UNKNOWN), THE
!     PROGRAM WILL ISSUE A WARNING MESSAGE (UNLESS IPRINT .LT. 0) AND THE
!     COVARIANCE MATRIX WILL NOT BE VALID.  IN THAT CASE, THE MODEL SHOULD BE
!     ALTERED TO INCLUDE A CONSTANT TERM (SET PHI(1) = 1.).

!        NOTE ALSO THAT, IN ORDER FOR THE USUAL ASSUMPTIONS TO HOLD,
!     THE OBSERVATIONS MUST ALL BE OF APPROXIMATELY THE SAME MAGNITUDE
!     (IN THE ABSENCE OF INFORMATION ABOUT THE ERROR OF EACH OBSERVATION),
!     OTHERWISE THE VARIANCES WILL NOT BE THE SAME.  IF THE OBSERVATIONS ARE
!     NOT THE SAME SIZE, THIS CAN BE CURED BY WEIGHTING.

!        IF THE USUAL ASSUMPTIONS HOLD, THE SQUARE ROOTS OF THE DIAGONALS OF
!     THE COVARIANCE MATRIX A GIVE THE STANDARD ERROR S(I) OF EACH PARAMETER.
!     DIVIDING A(I,J) BY S(I)*S(J) YIELDS THE CORRELATION MATRIX OF THE
!     PARAMETERS.  PRINCIPAL AXES AND CONFIDENCE ELLIPSOIDS CAN BE OBTAINED BY
!     PERFORMING AN EIGENVALUE/EIGENVECTOR ANALYSIS ON A.  ONE SHOULD CALL THE
!     EISPACK PROGRAM TRED2, FOLLOWED BY TQL2 (OR USE THE EISPAC CONTROL
!     PROGRAM).

!  CONVERGENCE FAILURES

!        IF CONVERGENCE FAILURES OCCUR, FIRST CHECK FOR INCORRECT
!     CODING OF THE SUBROUTINE ADA.  CHECK ESPECIALLY THE ACTION OF
!     ISEL, AND THE COMPUTATION OF THE PARTIAL DERIVATIVES.  IF THESE
!     ARE CORRECT, TRY SEVERAL STARTING GUESSES FOR ALF.  IF ADA
!     IS CODED CORRECTLY, AND IF ERROR RETURNS IERR = -2 OR -8
!     PERSISTENTLY OCCUR, THIS IS A SIGN OF ILL-CONDITIONING, WHICH
!     MAY BE CAUSED BY SEVERAL THINGS.  ONE IS POOR SCALING OF THE
!     PARAMETERS; ANOTHER IS AN UNFORTUNATE INITIAL GUESS FOR THE
!     PARAMETERS, STILL ANOTHER IS A POOR CHOICE OF THE MODEL.

!  ALGORITHM

!        THE RESIDUAL R IS MODIFIED TO INCORPORATE, FOR ANY FIXED ALF, THE
!     OPTIMAL LINEAR PARAMETERS FOR THAT ALF.  IT IS THEN POSSIBLE TO MINIMIZE
!     ONLY ON THE NONLINEAR PARAMETERS.  AFTER THE OPTIMAL VALUES OF THE
!     NONLINEAR PARAMETERS HAVE BEEN DETERMINED, THE LINEAR PARAMETERS CAN BE
!     RECOVERED BY LINEAR LEAST SQUARES TECHNIQUES (SEE REF. 1).

!        THE MINIMIZATION IS BY A MODIFICATION OF OSBORNE'S (REF. 3)
!     MODIFICATION OF THE LEVENBERG-MARQUARDT ALGORITHM.  INSTEAD OF
!     SOLVING THE NORMAL EQUATIONS WITH MATRIX

!              T      2
!            (J J + NU  * D),      WHERE  J = D(ETA)/D(ALF),

!     STABLE ORTHOGONAL (HOUSEHOLDER) REFLECTIONS ARE USED ON A
!     MODIFICATION OF THE MATRIX
!                                (   J  )
!                                (------) ,
!                                ( NU*D )

!     WHERE D IS A DIAGONAL MATRIX CONSISTING OF THE LENGTHS OF THE
!     COLUMNS OF J.  THIS MARQUARDT STABILIZATION ALLOWS THE ROUTINE
!     TO RECOVER FROM SOME RANK DEFICIENCIES IN THE JACOBIAN.
!     OSBORNE'S EMPIRICAL STRATEGY FOR CHOOSING THE MARQUARDT PARAM-
!     ETER HAS PROVEN REASONABLY SUCCESSFUL IN PRACTICE.  (GAUSS-
!     NEWTON WITH STEP CONTROL CAN BE OBTAINED BY MAKING THE CHANGE
!     INDICATED BEFORE THE INSTRUCTION LABELED 5).  A DESCRIPTION CAN
!     BE FOUND IN REF. (3), AND A FLOW CHART IN (2), P. 22.

!     FOR REFERENCE, SEE

!     1.  GENE H. GOLUB AND V. PEREYRA, 'THE DIFFERENTIATION OF
!         PSEUDO-INVERSES AND NONLINEAR LEAST SQUARES PROBLEMS WHOSE
!         VARIABLES SEPARATE,' SIAM J. NUMER. ANAL. 10, 413-432 (1973).
!     2.  ------, SAME TITLE, STANFORD C.S. REPORT 72-261, FEB. 1972.
!     3.  OSBORNE, MICHAEL R., 'SOME ASPECTS OF NON-LINEAR LEAST
!         SQUARES CALCULATIONS,' IN LOOTSMA, ED., 'NUMERICAL METHODS
!         FOR NON-LINEAR OPTIMIZATION,' ACADEMIC PRESS, LONDON, 1972.
!     4.  KROGH, FRED, 'EFFICIENT IMPLEMENTATION OF A VARIABLE PROJECTION
!         ALGORITHM FOR NONLINEAR LEAST SQUARES PROBLEMS,'
!         COMM. ACM 17, PP. 167-169 (MARCH, 1974).
!     5.  KAUFMAN, LINDA, 'A VARIABLE PROJECTION METHOD FOR SOLVING SEPARABLE
!         NONLINEAR LEAST SQUARES PROBLEMS', B.I.T. 15, 49-57 (1975).
!     6.  DRAPER, N., AND SMITH, H., APPLIED REGRESSION ANALYSIS,
!         WILEY, N.Y., 1966 (FOR STATISTICAL INFORMATION ONLY).
!     7.  C. LAWSON AND R. HANSON, SOLVING LEAST SQUARES PROBLEMS,
!         PRENTICE-HALL, ENGLEWOOD CLIFFS, N. J., 1974.

!                   JOHN BOLSTAD
!                   COMPUTER SCIENCE DEPT., SERRA HOUSE
!                   STANFORD UNIVERSITY
!                   JANUARY, 1977

!     ..................................................................

INTEGER, INTENT(IN)        :: l
INTEGER, INTENT(IN)        :: nl
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nmax
INTEGER, INTENT(IN)        :: lpp2
INTEGER, INTENT(IN)        :: iv
REAL (dp), INTENT(IN)      :: t(:,:)       ! t(nmax,iv)
REAL (dp), INTENT(IN)      :: y(:)
REAL (dp), INTENT(IN OUT)  :: w(:)
REAL (dp), INTENT(OUT)     :: a(:,:)       ! a(nmax,lpp2)
INTEGER, INTENT(IN)        :: iprint
REAL (dp), INTENT(IN OUT)  :: alf(:)       ! alf(nl)
REAL (dp), INTENT(IN OUT)  :: beta(:)
INTEGER, INTENT(OUT)       :: ierr

! Local variables

REAL (dp) :: acum, gnstep, nu, prjres, r, rnew
INTEGER   :: b1, isel, isub, iter, iterin, j, jsub, k, ksub, lnl2, lp1,  &
             modit, nlp1
LOGICAL   :: skip

INTERFACE
  SUBROUTINE ada (lp1, nl, n, nmax, lpp2, iv, a, inc, t, alf, isel)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(30, 60)
    INTEGER, INTENT(IN)     :: lp1
    INTEGER, INTENT(IN)     :: nl
    INTEGER, INTENT(IN)     :: n
    INTEGER, INTENT(IN)     :: nmax
    INTEGER, INTENT(IN)     :: lpp2
    INTEGER, INTENT(IN)     :: iv
    REAL (dp), INTENT(OUT)  :: a(:,:)
    INTEGER, INTENT(OUT)    :: inc(:,:)
    REAL (dp), INTENT(IN)   :: t(:,:)
    REAL (dp), INTENT(IN)   :: alf(:)
    INTEGER, INTENT(IN)     :: isel
  END SUBROUTINE ada
END INTERFACE

!     THE FOLLOWING TWO PARAMETERS ARE USED IN THE CONVERGENCE TEST:
!     EPS1 IS AN ABSOLUTE AND RELATIVE TOLERANCE FOR THE NORM OF THE PROJECTION
!     OF THE RESIDUAL ONTO THE RANGE OF THE JACOBIAN OF THE VARIABLE PROJECTION
!     FUNCTIONAL.
!     ITMAX IS THE MAXIMUM NUMBER OF FUNCTION AND DERIVATIVE EVALUATIONS
!     ALLOWED.
!     CAUTION:  EPS1 MUST NOT BE < 10 TIMES THE UNIT ROUND-OFF OF THE MACHINE.

REAL (dp), PARAMETER :: eps1 = 100._dp * EPSILON(1.0_dp)
INTEGER, PARAMETER   :: itmax = 50, output = 6

!-----------------------------------------------------------------
! CALL LIB MONITOR FROM VARPRO, MAINTENANCE NUMBER 509, DATE 77178
!        CALL LIBMON(509)
!***PLEASE DON'T REMOVE OR CHANGE THE ABOVE CALL. IT IS YOUR ONLY
!***PROTECTION AGAINST YOUR USING AN OUT-OF-DATE OR INCORRECT
!***VERSION OF THE ROUTINE. THE LIBRARY MONITOR REMOVES THIS CALL,
!***SO IT ONLY OCCURS ONCE, ON THE FIRST ENTRY TO THIS ROUTINE.
!-----------------------------------------------------------------
ierr = 1
iter = 0
lp1 = l + 1
b1 = l + 2
lnl2 = l + nl + 2
nlp1 = nl + 1
skip = .false.
modit = iprint
IF (iprint <= 0) modit = itmax + 2
nu = 0.
!              IF GAUSS-NEWTON IS DESIRED REMOVE THE NEXT STATEMENT.
nu = 1.

!              BEGIN OUTER ITERATION LOOP TO UPDATE ALF.
!              CALCULATE THE NORM OF THE RESIDUAL AND THE DERIVATIVE OF
!              THE MODIFIED RESIDUAL THE FIRST TIME, BUT ONLY THE
!              DERIVATIVE IN SUBSEQUENT ITERATIONS.

10 CALL vpdpa (l, nl, n, nmax, lpp2, iv, t, y, w, alf, ada, ierr,  &
               iprint, a, beta, a(:,lp1), r)
gnstep = 1.0
iterin = 0
IF (iter > 0) GO TO 20
IF (nl == 0) GO TO 160
IF (ierr /= 1) GO TO 180

IF (iprint <= 0) GO TO 20
WRITE (output,220) iterin,r
WRITE (output,190) nu
!                              BEGIN TWO-STAGE ORTHOGONAL FACTORIZATION
20 CALL vpfac1 (nlp1, n, l, iprint, a(:,b1:), prjres, ierr)
IF (ierr < 0) GO TO 180
ierr = 2
IF (nu == 0.) GO TO 40

!              BEGIN INNER ITERATION LOOP FOR GENERATING NEW ALF AND
!              TESTING IT FOR ACCEPTANCE.

30 CALL vpfac2 (nlp1, nu, a(:,b1:))

!              SOLVE A NL X NL UPPER TRIANGULAR SYSTEM FOR DELTA-ALF.
!              THE TRANSFORMED RESIDUAL (IN COL. LNL2 OF A) IS OVER-
!              WRITTEN BY THE RESULT DELTA-ALF.

40 CALL vpbsol (nl, a(:,b1:), a(:,lnl2))
a(1:nl,b1) = alf(1:nl) + a(1:nl,lnl2)
!           NEW ALF(K) = ALF(K) + DELTA ALF(K)

!              STEP TO THE NEW POINT NEW ALF, AND COMPUTE THE NEW
!              NORM OF RESIDUAL.  NEW ALF IS STORED IN COLUMN B1 OF A.

60 CALL vpdpa (l, nl, n, nmax, lpp2, iv, t, y, w, a(:,b1), ada, ierr,  &
               iprint, a, beta, a(:,lp1), rnew)
IF (ierr /= 2) GO TO 180
iter = iter + 1
iterin = iterin + 1
skip = MOD(iter,modit) /= 0
IF (skip) GO TO 70
WRITE (output,200) iter
WRITE (output,240) a(1:nl,b1)
WRITE (output,220) iterin, rnew

70 IF (iter < itmax) GO TO 80
ierr = -1
CALL vperr (iprint, ierr, 1)

! 1000       FORMAT(15X,'IERR AT 45 =',I7)
GO TO 170
80 IF (rnew-r < eps1*(r+1.0_dp)) GO TO 130

!              RETRACT THE STEP JUST TAKEN

IF (nu /= 0.) GO TO 100
!                                             GAUSS-NEWTON OPTION ONLY
gnstep = 0.5*gnstep
IF (gnstep < eps1) GO TO 170
DO  k=1,nl
  a(k,b1) = alf(k) + gnstep*a(k,lnl2)
END DO
GO TO 60
!                                        ENLARGE THE MARQUARDT PARAMETER
100 nu = 1.5*nu
IF (.NOT.skip) WRITE (output,210) nu
IF (nu <= 100.) GO TO 110
ierr = -2
CALL vperr (iprint, ierr, 1)
!       WRITE(6,1000)IERR
! 1000  FORMAT(15X,'IERR AT 60 =',I7)

GO TO 170
!                                        RETRIEVE UPPER TRIANGULAR FORM
!                                        AND RESIDUAL OF FIRST STAGE.
110 DO  k=1,nl
  ksub = lp1 + k
  DO  j=k,nlp1
    jsub = lp1 + j
    isub = nlp1 + j
    a(k,jsub) = a(isub,ksub)
  END DO
END DO
GO TO 30
!                                        END OF INNER ITERATION LOOP
!           ACCEPT THE STEP JUST TAKEN

130 r = rnew
DO  k=1,nl
  alf(k) = a(k,b1)
END DO
!                                        CALC. NORM(DELTA ALF)/NORM(ALF)
acum = gnstep*vpnorm(nl,a(:,lnl2)) / vpnorm(nl,alf)

!           IF ITERIN IS > 1, A STEP WAS RETRACTED DURING THIS OUTER ITERATION.

IF (iterin == 1) nu = 0.5*nu
IF (skip) GO TO 150
WRITE (output,190) nu
WRITE (output,230) acum
150 ierr = 3
IF (prjres > eps1*(r+1.0_dp)) GO TO 10
!           END OF OUTER ITERATION LOOP

!           CALCULATE FINAL QUANTITIES -- LINEAR PARAMETERS, RESIDUALS,
!           COVARIANCE MATRIX, ETC.

160 ierr = iter
170 IF (nl > 0) THEN
  isel = 4
  CALL vpdpa (l, nl, n, nmax, lpp2, iv, t, y, w, alf, ada, isel, iprint, a,  &
              beta, a(:,lp1), r)
END IF
CALL vppost (l, nl, n, lnl2, eps1, r, iprint, alf, w, a, a(:,lp1), beta)
180 RETURN

190 FORMAT ('     NU =',e15.7)
200 FORMAT ('0  ITERATION',i4,'    NONLINEAR PARAMETERS')
210 FORMAT ('     STEP RETRACTED, NU =',e15.7)
220 FORMAT ('0',i5,'  NORM OF RESIDUAL =',e15.7)
230 FORMAT ('     NORM(DELTA-ALF) / NORM(ALF) =',e12.3)
240 FORMAT ('0',7E15.7)
END SUBROUTINE varpro



SUBROUTINE vpfac1 (nlp1, n, l, iprint, b, prjres, ierr)

!            STAGE 1:  HOUSEHOLDER REDUCTION OF

!                   (    .    )      ( DR'. R3 )    NL
!                   ( DR . R2 )  TO  (----. -- ),
!                   (    .    )      (  0 . R4 )  N-L-NL

!                     NL    1          NL   1

!         WHERE DR = -D(Q2)*Y IS THE DERIVATIVE OF THE MODIFIED RESIDUAL
!         PRODUCED BY VPDPA, R2 IS THE TRANSFORMED RESIDUAL FROM DPA,
!         AND DR' IS IN UPPER TRIANGULAR FORM (AS IN REF. (2), P. 18).
!         DR IS STORED IN ROWS L+1 TO N AND COLUMNS L+2 TO L + NL + 1 OF
!         THE MATRIX A (I.E., COLUMNS 1 TO NL OF THE MATRIX B).  R2 IS
!         STORED IN COLUMN L + NL + 2 OF THE MATRIX A (COLUMN NL + 1 OF
!         B).  FOR K = 1, 2, ..., NL, FIND REFLECTION I - U * U' / BETA
!         WHICH ZEROES B(I, K), I = L+K+1, ..., N.

!     ..................................................................

INTEGER, INTENT(IN)        :: nlp1
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: l
INTEGER, INTENT(IN)        :: iprint
REAL (dp), INTENT(IN OUT)  :: b(:,:)     ! b(nmax,nlp1)
REAL (dp), INTENT(OUT)     :: prjres
INTEGER, INTENT(IN OUT)    :: ierr

! Local variables

REAL (dp) :: acum, alpha, beta, u
INTEGER   :: i, j, jsub, k, kp1, lp1, lpk, nl, nl23

nl = nlp1 - 1
nl23 = 2*nl + 3
lp1 = l + 1

DO  k=1,nl
  lpk = l + k
  alpha = SIGN(vpnorm(n+1-lpk, b(lpk:,k)), b(lpk,k))
  u = b(lpk,k) + alpha
  b(lpk,k) = u
  beta = alpha*u
  IF (alpha /= 0.0_dp) GO TO 10
!                                                   COLUMN WAS ZERO
  ierr = -8
  CALL vperr (iprint, ierr, lp1+k)
  GO TO 70
!                                APPLY REFLECTIONS TO REMAINING COLUMNS
!                                OF B AND TO RESIDUAL VECTOR.
  10 kp1 = k+1
  DO  j = kp1,nlp1
    acum = DOT_PRODUCT( b(lpk:n,k), b(lpk:n,j) )
    acum = acum/beta
    DO  i=lpk,n
      b(i,j) = b(i,j) - b(i,k)*acum
    END DO
  END DO
  b(lpk,k) = -alpha
END DO

prjres = vpnorm(nl, b(lp1:,nlp1))

!           SAVE UPPER TRIANGULAR FORM AND TRANSFORMED RESIDUAL, FOR USE
!           IN CASE A STEP IS RETRACTED.  ALSO COMPUTE COLUMN LENGTHS.

IF (ierr == 4) GO TO 70
DO  k=1,nl
  lpk = l + k
  DO  j=k,nlp1
    jsub = nlp1 + j
    b(k,j) = b(lpk,j)
    b(jsub,k) = b(lpk,j)
  END DO
  b(nl23,k) = vpnorm(k, b(lp1:,k))
END DO

70 RETURN
END SUBROUTINE vpfac1



SUBROUTINE vpfac2 (nlp1, nu, b)

!        STAGE 2:  SPECIAL HOUSEHOLDER REDUCTION OF

!                      NL     ( DR' . R3 )      (DR'' . R5 )
!                             (-----. -- )      (-----. -- )
!                  N-L-NL     (  0  . R4 )  TO  (  0  . R4 )
!                             (-----. -- )      (-----. -- )
!                      NL     (NU*D . 0  )      (  0  . R6 )

!                                NL    1          NL    1

!        WHERE DR', R3, AND R4 ARE AS IN VPFAC1, NU IS THE MARQUARDT
!        PARAMETER, D IS A DIAGONAL MATRIX CONSISTING OF THE LENGTHS OF
!        THE COLUMNS OF DR', AND DR'' IS IN UPPER TRIANGULAR FORM.
!        DETAILS IN (1), PP. 423-424.  NOTE THAT THE (N-L-NL) BAND OF
!        ZEROES, AND R4, ARE OMITTED IN STORAGE.

!     ..................................................................

INTEGER, INTENT(IN)     :: nlp1
REAL (dp), INTENT(IN)   :: nu
REAL (dp), INTENT(OUT)  :: b(:,:)     ! b(nmax,nlp1)

! Local variables

REAL (dp) :: acum, alpha, beta, u
INTEGER   :: i, j, k, kp1, nl, nl2, nl23, nlpk, nlpkm1

nl = nlp1 - 1
nl2 = 2*nl
nl23 = nl2 + 3
DO  k=1,nl
  kp1 = k + 1
  nlpk = nl + k
  nlpkm1 = nlpk - 1
  b(nlpk,k) = nu*b(nl23,k)
  b(nl,k) = b(k,k)
  alpha = SIGN(vpnorm(k+1, b(nl:,k)), b(k,k))
  u = b(k,k) + alpha
  beta = alpha*u
  b(k,k) = -alpha
!                        THE K-TH REFLECTION MODIFIES ONLY ROWS K,
!                        NL+1, NL+2, ..., NL+K, AND COLUMNS K TO NL+1.
  DO  j=kp1,nlp1
    b(nlpk,j) = 0.0_dp
    acum = u*b(k,j)
    DO  i=nlp1,nlpkm1
      acum = acum + b(i,k)*b(i,j)
    END DO
    acum = acum/beta
    b(k,j) = b(k,j) - u*acum
    DO  i=nlp1,nlpk
      b(i,j) = b(i,j) - b(i,k)*acum
    END DO
  END DO
END DO

RETURN
END SUBROUTINE vpfac2



SUBROUTINE vpdpa (l, nl, n, nmax, lpp2, iv, t, y, w, alf, ada,  &
                  isel, iprint, a, u, r, rnorm)

!        COMPUTE THE NORM OF THE RESIDUAL (IF ISEL = 1 OR 2), OR THE
!        (N-L) X NL DERIVATIVE OF THE MODIFIED RESIDUAL (N-L) VECTOR
!        Q2*Y (IF ISEL = 1 OR 3).  HERE Q * PHI = S, I.E.,

!         L     ( Q1 ) (     .   .        )   ( S  . R1 .  F1  )
!               (----) ( PHI . Y . D(PHI) ) = (--- . -- . ---- )
!         N-L   ( Q2 ) (     .   .        )   ( 0  . R2 .  F2  )

!                 N       L    1      P         L     1     P

!        WHERE Q IS N X N ORTHOGONAL, AND S IS L X L UPPER TRIANGULAR.
!        THE NORM OF THE RESIDUAL = NORM(R2), AND THE DESIRED DERIVATIVE
!        ACCORDING TO REF. (5), IS
!                                               -1
!                    D(Q2 * Y) = -Q2 * D(PHI)* S  * Q1* Y.

!     ..................................................................

INTEGER, INTENT(IN)        :: l
INTEGER, INTENT(IN)        :: nl
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nmax
INTEGER, INTENT(IN)        :: lpp2
INTEGER, INTENT(IN)        :: iv
REAL (dp), INTENT(IN)      :: t(:,:)    ! t(nmax,iv)
REAL (dp), INTENT(IN)      :: y(:)
REAL (dp), INTENT(IN OUT)  :: w(:)
REAL (dp), INTENT(IN OUT)  :: alf(:)
INTEGER, INTENT(IN OUT)    :: isel
INTEGER, INTENT(IN)        :: iprint
REAL (dp), INTENT(OUT)     :: a(:,:)    ! a(nmax,lpp2)
REAL (dp), INTENT(IN OUT)  :: u(:)
REAL (dp), INTENT(OUT)     :: r(:)
REAL (dp), INTENT(OUT)     :: rnorm

INTERFACE
  SUBROUTINE ada (lp1, nl, n, nmax, lpp2, iv, a, inc, t, alf, isel)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(30, 60)
    INTEGER, INTENT(IN)     :: lp1
    INTEGER, INTENT(IN)     :: nl
    INTEGER, INTENT(IN)     :: n
    INTEGER, INTENT(IN)     :: nmax
    INTEGER, INTENT(IN)     :: lpp2
    INTEGER, INTENT(IN)     :: iv
    REAL (dp), INTENT(OUT)  :: a(:,:)
    INTEGER, INTENT(OUT)    :: inc(:,:)
    REAL (dp), INTENT(IN)   :: t(:,:)
    REAL (dp), INTENT(IN)   :: alf(:)
    INTEGER, INTENT(IN)     :: isel
  END SUBROUTINE ada
END INTERFACE

! Local variables

REAL (dp) :: acum, alpha, beta
INTEGER   :: i, inc(12,8), j, k, kp1, ksub, m
LOGICAL, SAVE   :: nowate, philp1
REAL (dp), SAVE :: saved
INTEGER, SAVE   :: firstc, firstr, lastc, lnl2, lp1, lp2, lpp1, ncon, nconp1

IF (isel /= 1) GO TO 10
lp1 = l + 1
lnl2 = l + 2 + nl
lp2 = l + 2
lpp1 = lpp2 - 1
firstc = 1
lastc = lpp1
firstr = lp1
CALL vpinit (l, nl, n, nmax, lpp2, iv, t, w, alf, ada, isel,  &
             iprint, a, inc, ncon, nconp1, philp1, nowate)
IF (isel /= 1) GO TO 180
GO TO 30

10 CALL ada (lp1, nl, n, nmax, lpp2, iv, a, inc, t, alf, MIN(isel,3))
IF (isel == 2) GO TO 20
!                                                 ISEL = 3 OR 4
firstc = lp2
lastc = lpp1
firstr = (4-isel)*l + 1
GO TO 70
!                                                 ISEL = 2
20 firstc = nconp1
lastc = lp1
IF (ncon == 0) GO TO 30
IF (a(1,ncon) == saved) GO TO 30
isel = -7
CALL vperr (iprint, isel, ncon)
GO TO 180
!                                                  ISEL = 1 OR 2
30 IF (philp1) GO TO 50
r(1:n) = y(1:n)
GO TO 70
50 r(1:n) = y(1:n) - r(1:n)
!                                             WEIGHT APPROPRIATE COLUMNS
70 IF (nowate) GO TO 90
DO  i=1,n
  acum = w(i)
  DO  j=firstc,lastc
    a(i,j) = a(i,j)*acum
  END DO
END DO

!     COMPUTE ORTHOGONAL FACTORIZATIONS BY HOUSEHOLDER REFLECTIONS.
!     IF ISEL = 1 OR 2, REDUCE PHI (STORED IN THE FIRST L COLUMNS OF THE
!     MATRIX A) TO UPPER TRIANGULAR FORM, (Q*PHI = S), AND TRANSFORM Y (STORED
!     IN COLUMN L+1), GETTING Q*Y = R.
!     IF ISEL = 1, ALSO TRANSFORM J = D PHI (STORED IN COLUMNS L+2 THROUGH
!     L+P+1 OF THE MATRIX A), GETTING Q*J = F.
!     IF ISEL = 3 OR 4, PHI HAS ALREADY BEEN REDUCED, TRANSFORM ONLY J.
!     S, R, AND F OVERWRITE PHI, Y, AND J, RESPECTIVELY, AND A FACTORED FORM
!     OF Q IS SAVED IN U AND THE LOWER TRIANGLE OF PHI.

90 IF (l == 0) GO TO 130
DO  k=1,l
  kp1 = k + 1
  IF (isel >= 3 .OR. (isel == 2 .AND. k < nconp1)) GO TO 100
  alpha = SIGN(vpnorm(n+1-k, a(k:,k)), a(k,k))
  u(k) = a(k,k) + alpha
  a(k,k) = -alpha
  firstc = kp1
  IF (alpha /= 0.0_dp) GO TO 100
  isel = -8
  CALL vperr (iprint, isel, k)
  GO TO 180
!                                        APPLY REFLECTIONS TO COLUMNS
!                                        FIRSTC TO LASTC.
  100 beta = -a(k,k)*u(k)
  DO  j=firstc,lastc
    acum = u(k)*a(k,j) + DOT_PRODUCT( a(kp1:n,k), a(kp1:n,j) )
    acum = acum/beta
    a(k,j) = a(k,j) - u(k)*acum
    DO  i=kp1,n
      a(i,j) = a(i,j) - a(i,k)*acum
    END DO
  END DO
END DO

130 IF (isel >= 3) GO TO 140
rnorm = vpnorm(n-l,r(lp1:))
IF (isel == 2) GO TO 180
IF (ncon > 0) saved = a(1,ncon)

!     F2 IS NOW CONTAINED IN ROWS L+1 TO N AND COLUMNS L+2 TO L+P+1 OF THE
!     MATRIX A.   NOW SOLVE THE L X L UPPER TRIANGULAR SYSTEM S*BETA = R1 FOR
!     THE LINEAR PARAMETERS BETA.  BETA OVERWRITES R1.

140 IF (l > 0) CALL vpbsol (l, a, r)

!           MAJOR PART OF KAUFMAN'S SIMPLIFICATION OCCURS HERE.  COMPUTE
!           THE DERIVATIVE OF ETA WITH RESPECT TO THE NONLINEAR PARAMETERS

!   T   D ETA        T    L          D PHI(J)    D PHI(L+1)
!  Q * --------  =  Q * (SUM BETA(J) --------  + ----------)  =  F2*BETA
!      D ALF(K)          J=1         D ALF(K)     D ALF(K)

!           AND STORE THE RESULT IN COLUMNS L+2 TO L+NL+1.  IF ISEL NOT
!           = 4, THE FIRST L ROWS ARE OMITTED.  THIS IS -D(Q2)*Y.  IF
!           ISEL NOT = 4 THE RESIDUAL R2 = Q2*Y (IN COL. L+1) IS COPIED
!           TO COLUMN L+NL+2.  OTHERWISE ALL OF COLUMN L+1 IS COPIED.

DO  i=firstr,n
  IF (l == ncon) GO TO 170
  m = lp1
  DO  k=1,nl
    acum = 0.0_dp
    DO  j=nconp1,l
      IF (inc(k,j) == 0) CYCLE
      m = m + 1
      acum = acum + a(i,m)*r(j)
    END DO
    ksub = lp1+k
    IF (inc(k,lp1) == 0) GO TO 160
    m = m+1
    acum = acum + a(i,m)
    160 a(i,ksub) = acum
  END DO
  170 a(i,lnl2) = r(i)
END DO

180 RETURN
END SUBROUTINE vpdpa



SUBROUTINE vpinit (l, nl, n, nmax, lpp2, iv, t, w, alf, ada, isel,  &
                   iprint, a, inc, ncon, nconp1, philp1, nowate)

!        CHECK VALIDITY OF INPUT PARAMETERS, AND DETERMINE NUMBER OF
!        CONSTANT FUNCTIONS.

!     ..................................................................

INTEGER, INTENT(IN)        :: l
INTEGER, INTENT(IN)        :: nl
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nmax
INTEGER, INTENT(IN)        :: lpp2
INTEGER, INTENT(IN)        :: iv
REAL (dp), INTENT(IN)      :: t(:,:)   ! t(nmax,iv)
REAL (dp), INTENT(IN OUT)  :: w(:)
REAL (dp), INTENT(IN OUT)  :: alf(:)
INTEGER, INTENT(IN OUT)    :: isel
INTEGER, INTENT(IN)        :: iprint
REAL (dp), INTENT(IN OUT)  :: a(:,:)   ! a(nmax,lpp2)
INTEGER, INTENT(OUT)       :: inc(12,8)
INTEGER, INTENT(OUT)       :: ncon
INTEGER, INTENT(OUT)       :: nconp1
LOGICAL, INTENT(OUT)       :: philp1
LOGICAL, INTENT(OUT)       :: nowate

INTERFACE
  SUBROUTINE ada (lp1, nl, n, nmax, lpp2, iv, a, inc, t, alf, isel)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(30, 60)
    INTEGER, INTENT(IN)     :: lp1
    INTEGER, INTENT(IN)     :: nl
    INTEGER, INTENT(IN)     :: n
    INTEGER, INTENT(IN)     :: nmax
    INTEGER, INTENT(IN)     :: lpp2
    INTEGER, INTENT(IN)     :: iv
    REAL (dp), INTENT(OUT)  :: a(:,:)
    INTEGER, INTENT(OUT)    :: inc(:,:)
    REAL (dp), INTENT(IN)   :: t(:,:)
    REAL (dp), INTENT(IN)   :: alf(:)
    INTEGER, INTENT(IN)     :: isel
  END SUBROUTINE ada
END INTERFACE

! Local variables

INTEGER  :: i, inckj, j, k, lnl2, lp1, p
INTEGER, PARAMETER  :: output = 6

lp1 = l+1
nconp1 = lp1
lnl2 = l+2+nl
!                                          CHECK FOR VALID INPUT
IF (l >= 0 .AND. nl >= 0 .AND. l+nl < n .AND. lnl2 <= lpp2 .AND.   &
    2*nl + 3 <= nmax .AND. n <= nmax .AND. iv > 0 .AND. .NOT.(nl == 0  &
    .AND. l == 0)) GO TO 10
isel = -4
CALL vperr (iprint, isel, 1)
GO TO 90

10 IF (l == 0 .OR. nl == 0) GO TO 30
DO  j=1,lp1
  DO  k=1,nl
    inc(k,j) = 0
  END DO
END DO

30 CALL ada (lp1, nl, n, nmax, lpp2, iv, a, inc, t, alf, isel)

nowate = .true.
DO  i=1,n
  nowate = nowate .AND. (w(i) == 1.0_dp)
  IF (w(i) >= 0.0_dp) GO TO 40
!                                                ERROR IN WEIGHTS
  isel = -6
  CALL vperr (iprint, isel, i)
  GO TO 90
  40 w(i) = SQRT(w(i))
END DO

ncon = l
philp1 = (l == 0)
IF (philp1 .OR. nl == 0) GO TO 90
!                                   CHECK INC MATRIX FOR VALID INPUT AND
!                                   DETERMINE NUMBER OF CONSTANT FCNS.
p = 0
DO  j=1,lp1
  IF (p  == 0) nconp1 = j
  DO  k=1,nl
    inckj = inc(k,j)
    IF (inckj /= 0 .AND. inckj /= 1) GO TO 60
    IF (inckj == 1) p = p + 1
  END DO
END DO

ncon = nconp1-1
IF (iprint >= 0) WRITE (output,100) ncon
IF (l+p+2 == lpp2) GO TO 70
!                                              INPUT ERROR IN INC MATRIX
60 isel = -5
CALL vperr (iprint, isel, 1)
GO TO 90
!                                 DETERMINE IF PHI(L+1) IS IN THE MODEL.
70 DO  k=1,nl
  IF (inc(k,lp1) == 1) philp1 = .true.
END DO

90 RETURN
100 FORMAT ('0  NUMBER OF CONSTANT FUNCTIONS =', i4/)
END SUBROUTINE vpinit



SUBROUTINE vpbsol (n, a, x)

!        BACKSOLVE THE N X N UPPER TRIANGULAR SYSTEM A*X = B.
!        THE SOLUTION X OVERWRITES THE RIGHT SIDE B.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: a(:,:)
REAL (dp), INTENT(IN OUT)  :: x(:)

! Local variables

REAL (dp) :: acum
INTEGER   :: i, iback, ip1, j, np1

x(n) = x(n)/a(n,n)
IF (n == 1) GO TO 30
np1 = n + 1
DO  iback=2,n
  i = np1 - iback
!           I = N-1, N-2, ..., 2, 1
  ip1 = i + 1
  acum = x(i)
  DO  j=ip1,n
    acum = acum - a(i,j)*x(j)
  END DO
  x(i) = acum/a(i,i)
END DO

30 RETURN
END SUBROUTINE vpbsol



SUBROUTINE vppost (l, nl, n, lnl2, eps, rnorm, iprint, alf,  &
                   w, a, r, u)

!        CALCULATE RESIDUALS, SAMPLE VARIANCE, AND COVARIANCE MATRIX.
!        ON INPUT, U CONTAINS INFORMATION ABOUT HOUSEHOLDER REFLECTIONS
!        FROM VPDPA.  ON OUTPUT, IT CONTAINS THE LINEAR PARAMETERS.

INTEGER, INTENT(IN)        :: l
INTEGER, INTENT(IN)        :: nl
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: lnl2
REAL (dp), INTENT(IN)      :: eps
REAL (dp), INTENT(IN)      :: rnorm
INTEGER, INTENT(IN)        :: iprint
REAL (dp), INTENT(IN OUT)  :: alf(:)
REAL (dp), INTENT(IN OUT)  :: w(:)
REAL (dp), INTENT(IN OUT)  :: a(:,:)
REAL (dp), INTENT(IN OUT)  :: r(:)
REAL (dp), INTENT(IN OUT)  :: u(:)

! Local variables

REAL (dp) :: acum, prjres, SAVE
INTEGER   :: i, isel, k, kback, kp1, lp1, lpnl, lnl1
INTEGER, PARAMETER :: output = 6

lp1 = l + 1
lpnl = lnl2 - 2
lnl1 = lpnl + 1
DO  i=1,n
  w(i) = w(i)**2
END DO

!              UNWIND HOUSEHOLDER TRANSFORMATIONS TO GET RESIDUALS,
!              AND MOVE THE LINEAR PARAMETERS FROM R TO U.

IF (l == 0) GO TO 40
DO  kback=1,l
  k = lp1 - kback
  kp1 = k + 1
  acum = DOT_PRODUCT( a(kp1:n,k), r(kp1:n) )
  SAVE = r(k)
  r(k) = acum/a(k,k)
  acum = -acum/(u(k)*a(k,k))
  u(k) = SAVE
  DO  i=kp1,n
    r(i) = r(i) - a(i,k)*acum
  END DO
END DO
!                                              COMPUTE MEAN ERROR
40 acum = SUM( r(1:n) )
SAVE = acum/n

!              THE FIRST L COLUMNS OF THE MATRIX HAVE BEEN REDUCED TO
!              UPPER TRIANGULAR FORM IN VPDPA.  FINISH BY REDUCING ROWS
!              L+1 TO N AND COLUMNS L+2 THROUGH L+NL+1 TO TRIANGULAR
!              FORM.  THEN SHIFT COLUMNS OF DERIVATIVE MATRIX OVER ONE
!              TO THE LEFT TO BE ADJACENT TO THE FIRST L COLUMNS.

IF (nl == 0) GO TO 70
isel = 4
CALL vpfac1 (nl+1, n, l, iprint, a(:,l+2:), prjres, isel)
DO  i=1,n
  a(i,lnl2) = r(i)
  DO  k=lp1,lnl1
    a(i,k) = a(i,k+1)
  END DO
END DO
!                                              COMPUTE COVARIANCE MATRIX
70 a(1,lnl2) = rnorm
acum = rnorm*rnorm/(n-l-nl)
a(2,lnl2) = acum
CALL vpcov (lpnl, acum, a)

IF (iprint < 0) GO TO 80
WRITE (output,90)
IF (l > 0) WRITE (output,100) u(1:l)
IF (nl > 0) WRITE (output,110) alf(1:nl)
WRITE (output,120) rnorm, SAVE, acum
IF (ABS(SAVE) > eps) WRITE (output,130)
WRITE (output,90)
80 RETURN

90 FORMAT ('0',50(''''))
100 FORMAT ('0  LINEAR PARAMETERS'//(7E15.7))
110 FORMAT ('0  NONLINEAR PARAMETERS'//(7E15.7))
120 FORMAT ('0  NORM OF RESIDUAL =',e15.7,  &
    ' EXPECTED ERROR OF OBSERVATIONS =',e15.7,/  &
    '   ESTIMATED VARIANCE OF OBSERVATIONS =',e15.7)
130 FORMAT('  WARNING -- EXPECTED ERROR OF OBSERVATIONS IS NOT ZERO.',  &
    ' COVARIANCE MATRIX MAY BE MEANINGLESS.'/)
END SUBROUTINE vppost



SUBROUTINE vpcov (n, sigma2, a)

!           COMPUTE THE SCALED COVARIANCE MATRIX OF THE L + NL
!        PARAMETERS.  THIS INVOLVES COMPUTING

!                               2     -1    -T
!                          SIGMA  *  T   * T

!        WHERE THE (L+NL) X (L+NL) UPPER TRIANGULAR MATRIX T IS
!        DESCRIBED IN SUBROUTINE VPPOST.  THE RESULT OVERWRITES THE
!        FIRST L+NL ROWS AND COLUMNS OF THE MATRIX A.  THE RESULTING
!        MATRIX IS SYMMETRIC.  SEE REF. 7, PP. 67-70, 281.

!     ..................................................................

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: sigma2
REAL (dp), INTENT(OUT)  :: a(:,:)

! Local variables

REAL (dp) :: sum
INTEGER   :: i, ip1, j, jm1, nm1

DO  j=1,n
  a(j,j) = 1._dp / a(j,j)
END DO

!                 INVERT T UPON ITSELF

IF (n == 1) GO TO 40
nm1 = n - 1
DO  i=1,nm1
  ip1 = i + 1
  DO  j=ip1,n
    jm1 = j - 1
    sum = DOT_PRODUCT( a(i,i:jm1), a(i:jm1,j) )
    a(i,j) = -sum*a(j,j)
  END DO
END DO

!                 NOW FORM THE MATRIX PRODUCT

40 DO  i=1,n
  DO  j=i,n
    sum = DOT_PRODUCT( a(i,j:n), a(j,j:n) )
    sum = sum*sigma2
    a(i,j) = sum
    a(j,i) = sum
  END DO
END DO

RETURN
END SUBROUTINE vpcov



SUBROUTINE vperr (iprint, ierr, k)

!        PRINT ERROR MESSAGES

INTEGER, INTENT(IN)  :: iprint
INTEGER, INTENT(IN)  :: ierr
INTEGER, INTENT(IN)  :: k

INTEGER :: errno
INTEGER, PARAMETER :: output = 6

IF (iprint < 0) GO TO 80
errno = ABS(ierr)
SELECT CASE ( errno )
  CASE (    1)
    WRITE (output,90)
  CASE (    2)
    WRITE (output,100)
  CASE (    3)
    RETURN
  CASE (    4)
    WRITE (output,110)
  CASE (    5)
    WRITE (output,120)
  CASE (    6)
    WRITE (output,130) k
  CASE (    7)
    WRITE (output,140) k
  CASE (    8)
    WRITE (output,150) k
END SELECT

80 RETURN

90 FORMAT ('0  PROBLEM TERMINATED FOR EXCESSIVE ITERATIONS'//)
100 FORMAT ('0  PROBLEM TERMINATED BECAUSE OF ILL-CONDITIONING'//)
110 FORMAT (/' INPUT ERROR IN PARAMETER L, NL, N, LPP2, OR NMAX.'/)
120 FORMAT ('0  ERROR -- INC MATRIX IMPROPERLY SPECIFIED, OR ',  &
    'DISAGREES WITH LPP2.'/)
130 FORMAT ('0  ERROR -- WEIGHT(',i4,') IS NEGATIVE.'/)
140 FORMAT ('0  ERROR -- CONSTANT COLUMN ', i3,  &
    ' MUST BE COMPUTED ONLY WHEN ISEL = 1.'/)
150 FORMAT ('0  CATASTROPHIC FAILURE -- COLUMN', i4,  &
    ' IS ZERO, SEE DOCUMENTATION.'/)
END SUBROUTINE vperr



FUNCTION vpnorm (n, x) RESULT(fn_val)

!        COMPUTE THE L2 (EUCLIDEAN) NORM OF A VECTOR, MAKING SURE TO AVOID
!        UNNECESSARY UNDERFLOWS.  NO ATTEMPT IS MADE TO SUPPRESS OVERFLOWS.

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x(:)
REAL (dp)              :: fn_val

! Local variables
REAL (dp) :: rmax, sum, term
INTEGER   :: i

!           FIND LARGEST (IN ABSOLUTE VALUE) ELEMENT
rmax = 0._dp
DO  i=1,n
  IF (ABS(x(i)) > rmax) rmax = ABS(x(i))
END DO

sum = 0._dp
IF (rmax == 0._dp) GO TO 30
DO  i=1,n
  term = 0._dp
  IF (rmax+ABS(x(i)) /= rmax) term = x(i)/rmax
  sum = sum + term*term
END DO

30 fn_val = rmax*SQRT(sum)
RETURN
END FUNCTION vpnorm

END MODULE separable_leastsq
