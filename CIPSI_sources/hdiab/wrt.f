c-----------------------------------------------------------
      SUBROUTINE WRT(A,N,K1,K2,MZ)
C
C  IMPRESSION D UNE MATRICE
C  K1,K2 :COLONNES
C  N DIMENSION DES LIGNES
C
      REAL*8 A(MZ,MZ)
          MTER=10
      DO 200 K=K1,K2,MTER
      KP=MIN0(K2,K+MTER-1)
      WRITE(6,210)
      DO 200 I=1,N
200   WRITE(6,210) (A(I,J),J=K,KP)
210   FORMAT(1X,12F10.6)
      RETURN
      END
