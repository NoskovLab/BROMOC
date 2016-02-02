!    PB-PNP - Poisson-Boltzmann and Poisson-Nernst-Planck Equations Solver
!    Copyright (C) 2014 Pablo M. De Biase (pablodebiase@gmail.com)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

      SUBROUTINE EIGRS(A,N,JOBN,D,Z,IZ,WK,IER)
C---------------------------------------------------------------

      implicit none

C                                  SPECIFICATIONS FOR ARGUMENTS
      integer*4            N,JOBN,IZ,IER
      REAL*8   A(*),D(*),WK(*),Z(IZ,*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      integer*4            IJOB,IR,JR,IJ,JI,NP1
      integer*4            JER,NA,ND,IIZ,IBEG,IL,KK,LK,I,J,K,L
      REAL*8   ANORM,ASUM,PI,SUMZ,SUMR,AN,S,RDELP
      real*8   one,zero,thosnd,ten
CCC      DATA               RDELP/2.775557562D-17/
C                                  INITIALIZE ERROR PARAMETERS
C                                  FIRST EXECUTABLE STATEMENT
      RDELP=2.22045D-16
      IER = 0
      JER = 0
      zero   =    0.0D0
      one    =    1.0D0
      ten    =   10.0D0
      thosnd = 1000.0D0
      IF (JOBN.LT.10) GO TO 15
C                                  CONVERT TO SYMMETRIC STORAGE MODE
      K = 1
      JI = N-1
      IJ = 1
      DO 10 J=1,N
         DO 5 I=1,J
            A(K) = A(IJ)
            IJ = IJ+1
            K = K+1
    5    CONTINUE
         IJ = IJ + JI
         JI = JI - 1
   10 CONTINUE
   15 IJOB = MOD(JOBN,10)
      IF (IJOB.GE.0.AND.IJOB.LE.3) GO TO 20
C                                  WARNING ERROR - IJOB IS NOT IN THE
C                                    RANGE
      IER = 66
      IJOB = 1
      GO TO 25
   20 IF (IJOB.EQ.0) GO TO 35
   25 IF (IZ.GE.N) GO TO 30
C                                  WARNING ERROR - IZ IS LESS THAN N
C                                    EIGENVECTORS CAN NOT BE COMPUTED,
C                                    IJOB SET TO ZERO
      IER = 67
      IJOB = 0
   30 IF (IJOB.EQ.3) GO TO 75
   35 NA = (N*(N+1))/2
      IF (IJOB.NE.2) GO TO 45
      DO 40 I=1,NA
         WK(I) = A(I)
   40 CONTINUE
C                                  SAVE INPUT A IF IJOB = 2
   45 ND = 1
      IF (IJOB.EQ.2) ND = NA+1
C                                  REDUCE A TO SYMMETRIC TRIDIAGONAL
C                                    FORM
Cmu...25-Jul-93 Remove the second WK(ND) in the call to EHOUSS
Cmu   CALL EHOUSS(A,N,D,WK(ND),WK(ND))
      CALL EHOUSS(A,N,D,WK(ND))
      IIZ = 1
      IF (IJOB.EQ.0) GO TO 60
      IIZ = IZ
C                                  SET Z TO THE IDENTITY MATRIX
      DO 55 I=1,N
         DO 50 J=1,N
            Z(I,J) = ZERO
   50    CONTINUE
         Z(I,I) = ONE
   55 CONTINUE
C                                  COMPUTE EIGENVALUES AND EIGENVECTORS
   60 CALL EQRT2S(D,WK(ND),N,Z,IIZ,JER)
      IF (IJOB.EQ.0) GO TO 9000
      IF (JER.GT.128) GO TO 65
C                                  BACK TRANSFORM EIGENVECTORS
      CALL EHOBKS(A,N,1,N,Z,IZ)
   65 IF (IJOB.LE.1) GO TO 9000
C                                  MOVE INPUT MATRIX BACK TO A
      DO 70 I=1,NA
         A(I) = WK(I)
   70 CONTINUE
      WK(1) = THOSND
      IF (JER.NE.0) GO TO 9000
C                                  COMPUTE 1 - NORM OF A
   75 ANORM = ZERO
      IBEG = 1
      DO 85 I=1,N
         ASUM = ZERO
         IL = IBEG
         KK = 1
         DO 80 L=1,N
            ASUM = ASUM+ABS(A(IL))
            IF (L.GE.I) KK = L
            IL = IL+KK
   80    CONTINUE
         ANORM = MAX(ANORM,ASUM)
         IBEG = IBEG+I
   85 CONTINUE
      IF (ANORM.EQ.ZERO) ANORM = ONE
C                                  COMPUTE PERFORMANCE INDEX
      PI = ZERO
      DO 100 I=1,N
         IBEG = 1
         S = ZERO
         SUMZ = ZERO
         DO 95 L=1,N
            LK = IBEG
            KK = 1
            SUMZ = SUMZ+ABS(Z(L,I))
            SUMR = -D(I)*Z(L,I)
            DO 90 K=1,N
               SUMR = SUMR+A(LK)*Z(K,I)
               IF (K.GE.L) KK = K
               LK = LK+KK
   90       CONTINUE
            S = S+ABS(SUMR)
            IBEG = IBEG+L
   95    CONTINUE
         IF (SUMZ.EQ.ZERO) GO TO 100
         PI = MAX(PI,S/SUMZ)
  100 CONTINUE
      AN = N
      PI = PI/(ANORM*TEN*AN*RDELP)
      WK(1) = PI
      IF (JOBN.LT.10) GO TO 9000
C                                  CONVERT BACK TO FULL STORAGE MODE
      NP1 = N+1
      IJ = (N-1)*NP1 + 2
      K = (N*(NP1))/2
      DO 110 JR=1,N
         J = NP1-JR
         DO 105 IR=1,J
            IJ = IJ-1
            A(IJ) = A(K)
            K = K-1
  105    CONTINUE
         IJ = IJ-JR
  110 CONTINUE
      JI = 0
      K = N-1
      DO 120 I=1,N
         IJ = I-N
         DO 115 J=1,I
            IJ = IJ+N
            JI = JI+1
            A(IJ) = A(JI)
  115    CONTINUE
         JI = JI + K
         K = K-1
  120 CONTINUE
 9000 CONTINUE
CMJF  IF (IER.NE.0) CALL UERTST(IER,6HEIGRS )
      IF (JER.EQ.0) GO TO 9005
      IER = JER
CMJF  CALL UERTST(IER,6HEIGRS )
 9005 RETURN
      END


      SUBROUTINE EHOBKS(A,N,M1,M2,Z,IZ)
C
      implicit none
      integer*4 N,M1,M2,IZ
      REAL*8  A(*),Z(IZ,*)
C
      integer*4 I,J,K,L,IA
      REAL*8   H,S,one,zero

      zero   =  0.0D0
      one    =  1.0D0
C                                  FIRST EXECUTABLE STATEMENT
      IF (N .EQ. 1) GO TO 30
      DO 25 I=2,N
         L = I-1
         IA = (I*L)/2
         H = A(IA+I)
         IF (H.EQ.0.D0) GO TO 25
C                                  DERIVES EIGENVECTORS M1 TO M2 OF
C                                  THE ORIGINAL MATRIX FROM EIGENVECTORS
C                                  M1 TO M2 OF THE SYMMETRIC
C                                  TRIDIAGONAL MATRIX
         DO 20 J = M1,M2
            S = 0.0D0
            DO 10 K = 1,L
               S = S+A(IA+K)*Z(K,J)
   10       CONTINUE
            S = S/H
            DO 15 K=1,L
               Z(K,J) = Z(K,J)-S*A(IA+K)
   15       CONTINUE
   20    CONTINUE
   25 CONTINUE
   30 RETURN
      END


      SUBROUTINE EHOUSS(A,N,D,E)
C-----------------------------------------------------------------------
C     called by EIGRS
Cmu...25-Jul-93 M.E. Karpen.  The original E2 vector is the square
C     of the E vector.  The EIGRS routine passed the same vector (WK)
C     for both E & E2, which caused conflicting, machine-dependent results. 
C     Since EIGRS uses only the E vector, and is currently the only routine
C     that calls EHOUSS, E2 has been eliminated to prevent this ambiguity.
C

      implicit none

      integer*4  N
      REAL*8   A(*),D(N),E(N),one
C
      integer*4  JK1,JP1,NP1,I,J,K,L,NBEG,II,IK,JK,NK,NN
      REAL*8   ZERO,H,SCALE,F,G,HH
      DATA     ZERO/0.0D0/
c      zero   =  0.0
      one    =  1.0
C                                  FIRST EXECUTABLE STATEMENT
      NP1 = N+1
      NN = (N*NP1)/2-1
      NBEG = NN+1-N
      DO 70 II = 1,N
         I = NP1-II
         L = I-1
         H = ZERO
         SCALE = ZERO
         IF (L .LT. 1) GO TO 10
C                                  SCALE ROW (ALGOL TOL THEN NOT NEEDED)
         NK = NN
         DO 5 K = 1,L
            SCALE = SCALE+ABS(A(NK))
            NK = NK-1
    5    CONTINUE
         IF (SCALE .NE. ZERO) GO TO 15
   10    E(I) = ZERO
Cmu...25-Jul-93 remove E2 usage (NOT used)
Cmu         E2(I) = ZERO
         GO TO 65
   15    NK = NN
         DO 20 K = 1,L
            A(NK) = A(NK)/SCALE
            H = H+A(NK)*A(NK)
            NK = NK-1
   20    CONTINUE
Cmu         E2(I) = SCALE*SCALE*H
         F = A(NN)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE*G
         H = H-F*G
         A(NN) = F-G
         IF (L .EQ. 1) GO TO 55
         F = ZERO
         JK1 = 1
         DO 40 J = 1,L
            G = ZERO
            IK = NBEG+1
            JK = JK1
C                                  FORM ELEMENT OF A*U
            DO 25 K = 1,J
               G = G+A(JK)*A(IK)
               JK = JK+1
               IK = IK+1
   25       CONTINUE
            JP1 = J+1
            IF (L .LT. JP1) GO TO 35
            JK = JK+J-1
            DO 30 K = JP1,L
               G = G+A(JK)*A(IK)
               JK = JK+K
               IK = IK+1
   30       CONTINUE
C                                  FORM ELEMENT OF P
   35       E(J) = G/H
            F = F+E(J)*A(NBEG+J)
            JK1 = JK1+J
   40    CONTINUE
         HH = F/(H+H)
C                                  FORM REDUCED A
         JK = 1
         DO 50 J = 1,L
            F = A(NBEG+J)
            G = E(J)-HH*F
            E(J) = G
            DO 45 K = 1,J
               A(JK) = A(JK)-F*E(K)-G*A(NBEG+K)
               JK = JK+1
   45       CONTINUE
   50    CONTINUE
   55    DO 60 K = 1,L
            A(NBEG+K) = SCALE*A(NBEG+K)
   60    CONTINUE
   65    D(I) = A(NBEG+I)
         A(NBEG+I) = H*SCALE*SCALE
         NBEG = NBEG-I+1
         NN = NN-I
   70 CONTINUE
      RETURN
      END


      SUBROUTINE EQRT2S(D,E,N,Z,IZ,IER)
C

      implicit none

      integer*4  N,IZ,IER
      REAL*8   D(*),E(*),Z(IZ,*)
C
      integer*4  IP1,MM1,I,J,K,L,M,L1,II,MM1PL
      REAL*8   B,C,F,G,H,P,R,S,RDELP,one,zero
CCC      DATA               RDELP/2.775557562D-17/
C                                  MOVE THE LAST N-1 ELEMENTS
C                                  OF E INTO THE FIRST N-1 LOCATIONS
C                                  FIRST EXECUTABLE STATEMENT
      zero   =  0.0D0
      one    =  1.0D0
      RDELP=2.22045D-16
      IER  = 0
      IF (N .EQ. 1) GO TO 9005
      DO 5  I=2,N
         E(I-1) = E(I)
    5 CONTINUE
      E(N) = ZERO
      B = ZERO
      F = ZERO
      DO  60  L=1,N
         J = 0
         H = RDELP*(ABS(D(L))+ABS(E(L)))
         IF (B.LT.H) B = H
C                                  LOOK FOR SMALL SUB-DIAGONAL ELEMENT
         DO 10  M=L,N
            K=M
            IF (ABS(E(K)) .LE. B) GO TO 15
   10    CONTINUE
   15    M = K
         IF (M.EQ.L) GO TO 55
   20    IF (J .EQ. 30) GO TO 85
         J = J+1
         L1 = L+1
         G = D(L)
         P = (D(L1)-G)/(E(L)+E(L))
         R = ABS(P)
         IF (RDELP*ABS(P) .LT. 1.0D0) R = SQRT(P*P+ONE)
         D(L) = E(L)/(P+SIGN(R,P))
         H = G-D(L)
         DO 25 I = L1,N
            D(I) = D(I)-H
   25    CONTINUE
         F = F+H
C                                  QL TRANSFORMATION
         P = D(M)
         C = ONE
         S = ZERO
         MM1 = M-1
         MM1PL = MM1+L
         IF (L.GT.MM1) GO TO 50
         DO 45 II=L,MM1
            I = MM1PL-II
            G = C*E(I)
            H = C*P
            IF (ABS(P).LT.ABS(E(I))) GO TO 30
            C = E(I)/P
            R = SQRT(C*C+ONE)
            E(I+1) = S*P*R
            S = C/R
            C = ONE/R
            GO TO 35
   30       C = P/E(I)
            R = SQRT(C*C+ONE)
            E(I+1) = S*E(I)*R
            S = ONE/R
            C = C*S
   35       P = C*D(I)-S*G
            D(I+1) = H+S*(C*G+S*D(I))
            IF (IZ .LT. N) GO TO 45
C                                  FORM VECTOR
            DO 40 K=1,N
               H = Z(K,I+1)
               Z(K,I+1) = S*Z(K,I)+C*H
               Z(K,I) = C*Z(K,I)-S*H
   40       CONTINUE
   45    CONTINUE
   50    E(L) = S*P
         D(L) = C*P
         IF (ABS(E(L)) .GT.B) GO TO 20
   55    D(L) = D(L) + F
   60 CONTINUE
C                                  ORDER EIGENVALUES AND EIGENVECTORS
      DO  80  I=1,N
         K = I
         P = D(I)
         IP1 = I+1
         IF (IP1.GT.N) GO TO 70
         DO 65  J=IP1,N
            IF (D(J) .GE. P) GO TO 65
            K = J
            P = D(J)
   65    CONTINUE
   70    IF (K.EQ.I) GO TO 80
         D(K) = D(I)
         D(I) = P
         IF (IZ .LT. N) GO TO 80
         DO 75 J = 1,N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
   75    CONTINUE
   80 CONTINUE
      GO TO 9005
   85 IER = 128+L
CMJF  CALL UERTST(IER,6HEQRT2S)
 9005 RETURN
      END
