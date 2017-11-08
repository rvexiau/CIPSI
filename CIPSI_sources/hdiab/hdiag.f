c----------------------------------------------------------
      SUBROUTINE HDIAG(C,T,N,NDIM)
C     METHODE DE JACOBI
C     CALCUL DES VALEURS PROPRES ET DES VECTEURS PROPRES
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      DIMENSION C(NDIM,NDIM),T(NDIM,NDIM)
      R2=DSQRT(.5D0)
      DELTA=0.1D0
      DO 10 I=1,N
      DO 5 J=I,N
      T(I,J)=0.
5     T(J,I)=0.
10    T(I,I)=1.D0
35    N1=N-1
      I4=0
      DO 85 I=1,N1
      I2=I+1
      DO 85 J=I2,N
      A7=C(I,J)
      A9=A7*A7
      IF(A9-DELTA)85,85,40
40    I4=1
      A1=C(I,I)
      A2=C(J,J)
      A3=A2-A1
      IF(A3)50,45,50
45    A3=R2
      A4=A3
      A9=.5D0
      A10=A9
      A6=A7
      GO TO 57
50    A7=2.D0*A7
      ET=A7/A3
      A=ET/(DSQRT(1.D0+ET*ET)+1.D0)
      A9=1.D0/(1.D0+A*A)
      A10=1.D0-A9
      A3=DSQRT(A9)
      A4=A*A3
      A6=A7*A*A9
57    A5=A1*A9+A2*A10
      A8=A1*A10+A2*A9
      DO 84 K=1,N
      IF(I-K)60,70,60
60    IF(J-K)65,75,65
65    A9=C(K,I)
      A14=C(K,J)
      C(K,I)=A9*A3-A14*A4
      C(K,J)=A9*A4+A14*A3
      C(I,K)=C(K,I)
      C(J,K)=C(K,J)
      GO TO 80
70    C(I,K)=A5-A6
      C(I,J)=0.D0
      C(J,I)=0.D0
      GO TO 80
75    C(J,K)=A8+A6
80    A15=T(K,I)
      A16=T(K,J)
      T(K,I)=A15*A3-A16*A4
      T(K,J)=A15*A4+A16*A3
84    CONTINUE
85    CONTINUE
      IF(I4)35,90,35
90    IF(DELTA-1.D-18)100,100,95
95    DELTA=DELTA/10.
      GO TO 35
100   N1=N-1
      DO 115 I=1,N1
      IF(C(I,I)-C(I+1,I+1))115,115,105
105   A=C(I,I)
      C(I,I)=C(I+1,I+1)
      C(I+1,I+1)=A
      DO 110 K=1,N
      A=T(K,I)
      T(K,I)=T(K,I+1)
110   T(K,I+1)=A
      GO TO 100
115   CONTINUE
      DO 777 I=1,N
      IF(T(1,I)) 770,777,777
770   DO 778 J=1,N
778   T(J,I)=-T(J,I)
777   CONTINUE
      RETURN
      END
