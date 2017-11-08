      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(20),B(400),C(20)
      OPEN(UNIT=20,SHARED,FORM='UNFORMATTED',STATUS='OLD',
     &FILE='PSNL.TEI')
      OPEN(UNIT=10,FORM='FORMATTED',STATUS='NEW',FILE='PSNL.TEI_FORM')
    1 READ(20,END=9000) APSD,NP
      WRITE(6,2000) APSD,NP,ZN
      WRITE(10,1000) APSD,NP,ZN
      DO 10 I=1,NP
      READ(20)N,M,A,B,C
      WRITE(10,1001)N,M
      WRITE(10,1002) A
      WRITE(10,1002) B
      WRITE(10,1002) C
   10 CONTINUE
      GO TO 1
 1000 FORMAT(A8,I4,F10.5)
 1001 FORMAT(2I3)
 1002 FORMAT(4D20.12)
 9000 CONTINUE
 2000 FORMAT(' PSEUDO NAME ',A8,
     * ' NUMBER OF SYMMETRIES',I4,' NUM OF VALENCE ELECTRONS ',F5.2)  
      STOP
      END