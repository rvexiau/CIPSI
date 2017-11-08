      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*80 APSO
      DIMENSION Apsl(60),Cpslo(3600),Cpsl(600)
C
C     cre23
C
C     THIS PROGRAM READS PSEUDOPOTENTIAL DATA ON FILE 12 (FORMATTED)
C     "MOLCAS" FORMAT
C     AND CREATES A FILE 23 (UNFORMATTED) TO BE USED IN PSHF
C      
      OPEN(UNIT=12,FORM='FORMATTED',STATUS='OLD',
     *                   FILE='psnl_gaussian_txt')
      OPEN(UNIT=23,FORM='UNFORMATTED',STATUS='NEW',
     *                     FILE='PSNL_MOLCAS_AR1E')
C
   1  READ(12,1000,END=9000) APSO,MNO,ZATO
      WRITE(23) APSO,MNO-1,ZATO
      WRITE(6,999) APSO,MNO,ZATO
      DO lsy=1,MNO
        READ(12,*)    Lsym,Nops,Nvois
        READ(12,1113) (Cpsl(I),I=1,Nvois)
        read(12,1113) (Apsl(I),I=1,Nops)
        do Itl=0,Nops-1
        read(12,1113) (Cpslo(I),I=1+Itl,Nops*Nvois,Nops) 
        enddo 
        WRITE(23) Lsym,Nops,Nvois,(Cpsl(J),J=1,Nvois),
     *                               (Apsl(J),J=1,Nops),
     *                        (Cpslo(J),J=1,Nops*Nvois)
      END DO
      GO TO 1
C
  999 FORMAT(2X,'ATOMIC NAME=',1X,A12,' MAX.SYMM=',I4,
     *' NUMBER OF VALENCE ELECTRONS',F10.5)
 1000 FORMAT(A12,I4,F10.5)
 1002 FORMAT(4D20.12)
 1113 FORMAT(3D21.13)
 9000 CONTINUE
      STOP
      END
