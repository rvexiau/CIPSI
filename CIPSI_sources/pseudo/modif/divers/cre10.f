
      IMPLICIT REAL*8 (A-H,O-Z)
      character*80 apsd
      DIMENSION A(40),B(1600),C(40)
      DIMENSION A1(20),B1(400),C1(20)
      OPEN(UNIT=20,FORM='UNFORMATTED',STATUS='OLD',
     &FILE='PSNL')
      OPEN(UNIT=10,FORM='FORMATTED',STATUS='unknown',FILE='PSNL_FORM')
    1 READ(20,END=9000) APSD,NP
      write(6,*) apsd,np
      WRITE(10,1000) APSD,NP
      DO 10 Isym=1,NP
      READ(20)N,M
      read(20)a1
      read(20)b1
      read(20)c1
      WRITE(10,1001)amolcas,lsym,N,M
      WRITE(10,1002) c
      WRITE(10,1003) a
      WRITE(10,1004) b
   10 CONTINUE
      GO TO 1
 1000 FORMAT(A80,I4)
 1001 FORMAT('label =',a80, ' sym.max =', i3,/,' nombre d''exposants =',
     * i3,' nombre de combinaisons lineairement independantes =',i3)
 1002 FORMAT(' coefficients des operateurs',/, 4D20.12)
 1003 FORMAT(' exposants des operateurs',/, 4D20.12)
 1004 FORMAT(' coefficients des combinaisons lineaires ',/, 4D20.12)
 9000 CONTINUE
      STOP
      END
