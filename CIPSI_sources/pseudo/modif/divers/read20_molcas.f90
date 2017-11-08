      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*80 APSD
      DIMENSION A(20),B(400),C(20)
!
!     THIS PROGRAM READS PSEUDOPOTENTIAL DATA ON FILE 20 (UNFORMATTED)
!     USED IN PSHF OR MONSTER
!     AND CREATES A FILE 10 (FORMATTED) FOR READABILITY OR
!     PORTABILITY
! 
999  FORMAT(2X,'ATOMIC NAME=',1X,A8,' MAX.SYMM=',I3,' NUMBER OF VALENCE ELECTRONS',F5.2)
1000  FORMAT(A8,I4,F10.5)
1001  FORMAT(2I3)
1002  FORMAT(4D20.12) 
      OPEN(UNIT=10,FORM='FORMATTED',STATUS='NEW',FILE='PSNL_FORM')
      OPEN(UNIT=20,FORM='UNFORMATTED',STATUS='OLD',FILE='PSNL_MOLCAS_AR1E')
    1  READ(20,END=9000) APSD,NP,ZN
      NP=NP+1
      WRITE(10,1000) APSD,NP
      WRITE(6,999) APSD,NP,ZN

      DO I=1,NP
        READ(20) II,N,M,(A(j),j=1,M),(B(j),j=1,N),(C(j),j=1,N*M)
        WRITE  (10,1001) N,M
        WRITE(10,1002) A
        WRITE(10,1002) B
        WRITE(10,1002) C
      END DO
      GO TO 1

9000  CONTINUE
      STOP
      END