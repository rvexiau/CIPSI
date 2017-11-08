       FUNCTION DETERM(N,NA,A)
      IMPLICIT REAL*8(A-H,O,P,R-Z)
      DIMENSION A(NA,NA)
      TEST=1.D-10
       NP1=N+1
       DETERM=1.D0
       DO 1 L=1,N
       LP1=L+1
       IF(DABS(A(L,L)).GT.TEST) GO TO 2
C
C CAS OU LE PIVOT EST TROP PETIT
C
      IF(L.EQ.N) GO TO 18
      DO 5 I=LP1,N
      IF(DABS(A(I,L)).GT.TEST) GO TO 3
  5   CONTINUE
 18   DETERM=0.D0
  11  RETURN
   3  DO 4 J=1,N
       ALJ=A(L,J)
      A(L,J)=A(I,J)
  4   A(I,J)=ALJ
      DETERM=-DETERM
C
C  CALCUL NORMAL
C
 2     LM1=L-1
      DETERM=DETERM*A(L,L)
      IF(L.EQ.N) GO TO 11
       DO 6 J=LP1,N
   6   A(L,J)=A(L,J)/A(L,L)
       DO 1 I=LP1,N
       DO 70 J=LP1,N
 70   A(I,J)=A(I,J)-A(I,L)*A(L,J)
  1    CONTINUE
      RETURN
       END
