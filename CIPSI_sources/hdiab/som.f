C***********************************************************
C     RECOUVREMENT ENTR OM DE REFERENCES ET/OU OM UTILISEES
C     FILE  8 OM DE REFERENCES
C     FILE  2 DONNEES BASE PSHF
C     FILE 12 MATRICES DE RECOUVREMENT FINALES
C     FILE 11 OM ISSUES DE IJKL
C***********************************************************
C     
C     MODE D EMPLOI
C
C     NAOR    NOMBRE D OA DE REFERENCE
C     NOMR    NOMBRE D OM DE REFERENCE
C     VR      COEFFICIENTS DES OM DE REFERENCE
C     IREF    NUMERO DES OA DE REF UTILISEES SI NAOR<NAO
C     JREF    NUMERO DES OM DE REF UTILISEES
C     QFR     LECTURE DES OM DE REF SUR FILE 8
C     QPRT    IMPRESSIONS INTERMEDIAIRES
C     QRUHF   OM REFERENCE ALPHA NON EGAL A BETA     
C     NID     NOMBRE D OM DE REF IDENTIQUES AUX OM DU CALCUL
C************************************************************

      IMPLICIT REAL*8(A-H,O,P,R-Z),LOGICAL*1(Q)
      PARAMETER(MD=1000)
      DIMENSION SAO(MD,MD),OM(MD,MD),OMR(MD,MD),SMO(MD,MD),G(MD)
      DIMENSION IREF(MD),VR(MD*MD),JREF(MD)
      NAMELIST/DATOM/NOMR,IREF,JREF,VR,NAOR,QFR,QPRT,NID
      ndim=md
      QPRT=.FALSE.
      QFR=.FALSE.
      QRUHF=.FALSE.
      call openf
C  IDENTITE DE CERTAINES OM DE REF AVEC LES OM SCF
      NID=0
C
      DO 50 I=1,MD
      IREF(I)=I
      JREF(I)=I
      DO 50 J=1,MD
      OMR(I,J)=0.D0
50    OM(I,J)=0.D0
      READ(2) (IODA,I=1,19),NT,NSHELL,NUM
      NAO=NUM
      NOM=NUM
      READ(5,DATOM)
      KK=0
      IF(NID.LT.NOMR) THEN
      DO 250 I=NID+1,NOMR
      DO 250 J=1,NAOR
      KK=KK+1
250   OM(J,I)=VR(KK)
      ENDIF
C      REWIND 11
C      READ(11)
C      READ(11)
      READ(2) ((SAO(I,J),J=1,I),I=1,NAO)
      DO 20 I=1,NAO
      DO 20 J=1,I
   20 SAO(J,I)=SAO(I,J)
      IF(QPRT) THEN
      WRITE(6,*) ' RECOUVREMENT ENTRE OA'
      CALL WRT(SAO,NAO,1,NAO,NDIM)
      ENDIF
C   LECTURE DES OM DE REFERENCE
      IF(NID.NE.0) THEN
      REWIND 11
      READ(11) ((OM(I,J),I=1,NAO),J=1,NID)
      ENDIF
      QTR=QFR.AND.(NID.LT.NOMR)
      IF(QTR) THEN
      OPEN(UNIT=8,FORM='UNFORMATTED',STATUS='OLD')
      READ(8) ((OM(I,J),I=1,NAOR),J=NID+1,NOMR)
      IF(QRUHF) THEN
      READ(8) ((OM(I,J),I=1,NAOR),J=NID+1+NOMR,NOMR)
      NOMR=NOMR+NOMR
      ENDIF
      ENDIF
      DO 60 J=1,NOMR
      JJ=JREF(J)
      DO 60 I=1,NAOR
      II=IREF(I)
60    OMR(II,J)=OM(I,JJ)
      IF(QPRT) THEN
      WRITE(6,*) ' OM DE REFERENCE'
      CALL WRT(OMR,NAO,1,NOMR,NDIM)
      ENDIF
C
C RECOUVREMENT ENTRE OM DE REFERENCE
C
      DO 500 L=1,NOMR
      DO 550 J=1,NAO
      A=0.D0
      DO 560 I=1,NAO
560   A=A+OMR(I,L)*SAO(I,J)
550   G(J)=A
      DO 500 K=1,L
      A=0.D0
      DO 580 I=1,NAO
580   A=A+OMR(I,K)*G(I)
      SMO(L,K)=A
500   SMO(K,L)=A
      WRITE(12) NOMR,NOMR,((SMO(K,L),K=1,NOMR),L=1,NOMR)
      IF(QPRT) THEN
      WRITE(6,*) ' RECOUVREMENT ENTRE REFERENCES'
      CALL WRT(SMO,NOMR,1,NOMR,NDIM)
      ENDIF
C  LECTURE DES OM
      REWIND 11
      READ(11)((OM(I,J),I=1,NAO),J=1,NOM)
      IF(QPRT) THEN
      WRITE(6,*) ' OM DU CALCUL'
      CALL WRT(OM,NAO,1,NOM,NDIM)
      ENDIF
C
C  RECOUVREMENT ENTRE OM ET OM DE REFERENCE
C
      DO 5500 L=1,NOM
      DO 5550 J=1,NAO
      A=0.D0
      DO 5560 I=1,NAO
5560   A=A+OM(I,L)*SAO(I,J)
5550   G(J)=A
      DO 5500 K=1,NOMR
      A=0.D0
      DO 5580 I=1,NAO
5580   A=A+OMR(I,K)*G(I)
5500   SMO(K,L)=A
      WRITE(12) NOMR,NOM,((SMO(K,L),K=1,NOMR),L=1,NOM)
      IF(QPRT) THEN
      WRITE(6,*) ' RECOUVREMENT REFERENCES/OM'
      CALL WRT(SMO,NOMR,1,NOM,NDIM)
      ENDIF
      STOP
      END
      SUBROUTINE WRT(A,N,K1,K2,MD)
C
C  IMPRESSION D UNE MATRICE
C  K1,K2 :COLONNES
C  N DIMENSION DES LIGNES
C
      REAL*8 A(MD,MD)
          MTER=10
      DO 200 K=K1,K2,MTER
      KP=MIN0(K2,K+MTER-1)
      WRITE(6,210)
      DO 200 I=1,N
200   WRITE(6,210) (A(I,J),J=K,KP)
210   FORMAT(1X,12F10.6)
      RETURN
      END

C     ******************************************************************

      subroutine openf
      character*10 prefix,info,som,ijom,omref
      namelist/SOMFIL/PREFIX,info,ijom,som,omref
      info='info       '
      ijom='ijom       '
      som='som     '
      omref='omref       '
      prefix= '          '
      read(5,SOMFIL)
      write(6,somfil)
      call nomfil(8,omref,prefix,'UNKNOWN   ')
      call nomfil(2,info,prefix,'OLD       ')
      call nomfil(11,ijom,prefix,'OLD       ')
      call nomfil(12,som,prefix,'NEW       ')
      return
      end
      subroutine nomfil(iunit,racine,prefix,status)
      character*1 prefix(10),racine(10)
      character*10 status
      character*20 file
      character*1  fil(20)
      equivalence (fil(1),file)
      file='                                           '
C      OPEN FILES
      call nomf(fil,prefix,racine)
      write(6,*)'ouverture de la file',iunit
      write(6,*)fil,iunit,racine,status
      OPEN (UNIT=iunit,FORM='UNFORMATTED',STATUS=status,file=file)
      write(6,*)'fin ouverture des files'
      return
      end
      subroutine nomf(fil,prefix,racine)
      character*1 prefix(10),racine(10),blanc,fil(20)
      blanc= ' '
      do 1 i=1,10
	 if(prefix(i).eq.blanc)then
	    do 2 j=1,i-1
	       fil(j)=prefix(j)
 2          continue
	    k=i-1
	    goto 3
         endif
 1    continue
 3    continue
      do 4 i=1,10
	 if(racine(i).eq.blanc) then
	    do 5 j=1,i-1
	       fil(k+j)=racine(j)
 5          continue
	    k=k+i-1
	    goto 6
         endif
 4    continue
 6    continue
      write(6,*)fil
      return
      end
                                                                                
