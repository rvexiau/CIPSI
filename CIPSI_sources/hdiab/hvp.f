c--------------------------------------------------------
      SUBROUTINE HVP(H,NV,NDIM)
      IMPLICIT REAL*8(A-H,O,P,R-Z),LOGICAL*1(Q)
      include 'pshf.prm'
      DIMENSION H(metz,metz),HD(metz,metz),V(metz,metz)
      write(6,*)
      write(6,*) ' VERIF'
      DO 10 I=1,NV
      DO 15 J=1,I
      HD(I,J)=H(I,J)
15    HD(J,I)=H(I,J)
10    WRITE(6,100) (HD(I,J),J=1,I)
100   FORMAT(11F12.6)
      CALL HDIAG(HD,V,NV,NDIM)
      WRITE(6,105)
105   FORMAT(/,' VALEURS PROPRES ET VECTEURS PROPRES',/)
      WRITE(6,100) (HD(I,I),I=1,NV)
       WRITE(6,110)
110   FORMAT(/)
      DO 20 I=1,NV
20    WRITE(6,100) (V(I,J),J=1,NV)
      RETURN
      END
