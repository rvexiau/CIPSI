      IMPLICIT REAL*8 (A-H,O-Z)
      character*80 apsd
      DIMENSION A(40),B(1600),C(40)
      OPEN(UNIT=20,FORM='UNFORMATTED',STATUS='OLD',
     &FILE='MWBPOT')
      OPEN(UNIT=10,FORM='FORMATTED',STATUS='unknown',FILE='MWBPOT_FORM')
    1 READ(20,END=9000) APSD,NP
      write(6,*) apsd,np
      WRITE(10,1000) APSD,NP
      DO 10 Isym=1,NP+1
      READ(20)lsym,N,M,(c(i),i=1,m),(a(i),i=1,n),(b(i),i=1,n*m)
      write(6,*) isym,lsym,n,m   
      WRITE(10,1001)apsd,lsym,N,M
      WRITE(10,1002) (c(i),i=1,m)
      WRITE(10,1003) (a(i),i=1,n)
      WRITE(10,1004) (b(i),i=1,m*n)
   10 CONTINUE
      GO TO 1
 1000 FORMAT(A80,I4)
 1001 FORMAT('label =',a80, ' sym.max =', i3,/,' nombre d''exposants =',
     * i3,' nombre de combinaisons lineairement independantes =',i3)
 1002 FORMAT(' coefficients des operateurs',/, (4D20.12))
 1003 FORMAT(' exposants des operateurs',/, (4D20.12))
 1004 FORMAT(' coefficients des combinaisons lineaires ',/, (4D20.12))
 9000 CONTINUE
      STOP
      END
