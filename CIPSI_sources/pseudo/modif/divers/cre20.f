      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*8 APSD
      DIMENSION A(20),B(400),C(20)
C
C     THIS PROGRAM READS PSEUDOPOTENTIAL DATA ON FILE 10 (FORMATED)
C     AND CREATES A FILE 20 (UNFORMATED) TO BE USED IN PSGRAD OR PSHF
C      OR MONSTER
      OPEN(UNIT=10,FORM='FORMATTED',STATUS='OLD',FILE='PSNL_FORM')
      OPEN(UNIT=20,FORM='UNFORMATTED',STATUS='new',FILE='PSNL')
   1  READ(10,1000,END=9000) APSD,NP
      WRITE(20) APSD,NP,ZN
      WRITE(6,999) APSD,NP,ZN
 999  FORMAT(2X,'ATOMIC NAME=',1X,A8,' MAX.SYMM=',I3,
     *' NUMBER OF VALENCE ELECTRONS',F5.2)
      DO I=1,NP
        READ  (10,1001) N,M
        READ(10,1002) A
        READ(10,1002) B
        READ(10,1002) C
        WRITE(20) N,M,A,B,C
      END DO
      GO TO 1
1000  FORMAT(A8,I4,F10.5)
1001  FORMAT(2I3)
1002  FORMAT(4D20.12)
9000  CONTINUE
      STOP
      END
