c---------------------------------------------------------------
C   Programme d extraction de pseudoorbitales atomiques
c
c
c   tache
c
c   Extraction d une pseudoorbitale normalisee sans noeud 
c   a partir d 'une orbitale avec noeuds
c    r > rc       fit de l'orbitale vraie 
c    r > rc       prolongement  par une forme reguliere 
c
c   
c   makes a normalized  nodeless pseudoorbital
c   from an actual spinorbital including nodes
c   r > rc        fit of the actual orbital
c   r < rc        continuous and regular extension with smooth shape

c   La fonction de slater est r**np * exp(-br)
c   tabulation de la densite de probabilite radiale p(r)=r*f(r)
c
c   Donnees
c
c   namelist inporb
c
c   n             nombre de fonctions de Slater de base 
c   nn
c   (np(i),i=1,n) exposant de r
c   (b(i),i=1,n)  exposants de l orbitale
c   (c(i),i=1,n)  coefficients de l'orbitale 
c   (d(i),i=1,n)  coefficients de la derniere  orbitale de coeur
c                 de meme symetrie
c   im            nombre de points de la grille numerique
c   ircn = 0      rc est determine par intersection avec 
c                 la derniere orbitale de coeur
c        =1       rc est impose
c   rc            rayon de coupure
c   zat           charge atomique reelle
c   rmax          position du maximum de l orbitale
c   nnn           exposant du coefficient de r
c                 donne la forme a l origine
c   rx
c   rxeta
c   rmu           equilibre dans le fit poids partie interne/externe
c
c---------------------------------------------------------------

      program psorba
      IMPLICIT REAL*8 (A-H,O-Z)                                         00000030
      DIMENSION RK(120),XK(120),XU(30),EPREC(30),W(120) ,DP(30)         00000040
      DIMENSION RXE(30)                                                 00000050
      DIMENSION XPAR(30)                                                00000060
      DIMENSION B(30),P(120),Q(120),DIF(120),VAR(4),FONC(4),DERIV(4),   00000070
     1TRAB(4),XVAR(50),XFONC(50),U(50)                                  00000080
      COMMON/PSINF/D(30),C(30),CP(30),X(30),R,TEST,A1,NP(30),NNP(30),N, 00000090
     1NN                                                                00000100
      COMMON/DEVI/ RPRIM,RMU,RXETA,RX  ,CONSC,NNN                       00000110
      COMMON /OPTI/ DP                                                  00000120
      COMMON /FAC/ FACTOR(16)                                           00000130
      NAMELIST/INPORB/TITRE,N,NN,NP,D,B,C,IM,ZAT,NNN,IRCRN,RC,RMAX,
     1RMU,RX,RXETA
      character*80  TITRE
c
      READ(5,INPORB)
c
c  longueur du nom du fichier de sortie
c  le nom du fichier de sortie orbitale.plot
c  orbitale: debut du titre jusqu au premier caractere blanc    
c
      lo=0
      do i=1,80
       if(titre(i:i).eq.' ') go to 1644
       lo=lo+1
      enddo
1644  continue
     
  102 FORMAT(4I2)                                                       00000200
  103 FORMAT(I3,3D15.8)                                                 00000220
      DO99I=1,NN                                                        00000230
      INDEX=N-NN+I                                                      00000240
      NNP(I)=NP(INDEX)                                                  00000250
      DP(I)=D(INDEX)                                                    00000260
  99   CP(I)=C(INDEX)                                                   00000270
  105 FORMAT(2I3)                                                       00000290
  104 FORMAT(4D15.8)                                                    00000330
      NPREC=10                                                          00000340
      TEST=10.D0**(-NPREC)                                              00000350
C                                                                       00000360
C     CALCUL DES XK ET DES RK                                           00000370
C                                                                       00000380
      Z=ZAT                                                             00000390
      PI=3.14159265358979D0                                             00000400
      RK(1)=1.D-6                                                       00000410
      XK(1)=1.D-6                                                       00000420
      PAX=0.005D0                                                       00000430
      MIN=-8                                                            00000440
      MAX=1                                                             00000450
      DO600I=1,11                                                       00000460
      PAX=PAX*2.D0                                                      00000470
      MIN=MIN+10                                                        00000480
      MAX=MAX+10                                                        00000490
      DO600J=MIN,MAX                                                    00000500
      XK(J)=XK(J-1)+PAX                                                 00000510
      RK(J)=XK(J)*(9.D0*PI*PI/(128.D0*Z))**(1.D0/3.D0)                  00000520
  600 CONTINUE                                                          00000530
C                                                                       00000540
      S=0.D0                                                            00000550
      SS=0.D0                                                           00000560
      DO 40 I=1,N                                                       00000570
      S=S+DK(NP(I),NP(I),D(I),D(I))*DJ(2*NP(I),2.D0*D(I))*C(I)**2       00000580
      M=I-1                                                             00000590
      IF(I.EQ.1) GO TO 40                                               00000600
      DO 41 J=1,M                                                       00000610
      SS=SS+C(I)*C(J)*DK(NP(I),NP(J),D(I),D(J))*DJ(NP(I)+NP(J),D(I)+D(J)00000620
     1)                                                                 00000630
   41 CONTINUE                                                          00000640
   40 CONTINUE                                                          00000650
      REC=S+2.D0*SS                                                     00000660
      DO 42 I=1,N                                                       00000670
   42 C(I)=C(I)/DSQRT(REC)                                              00000680
C                                                                       00000690
C      CALCUL DE RC ET DE RMAX                                          00000700
C                                                                       00000710
      IF(IRCRN.EQ.0)GOTO601                                             00000720
C                                                                       00000730
C     CALCUL  DE  P  ET  Q                                              00000740
C                                                                       00000750
      DO810I=1,IM                                                       00000760
       S=0.D0                                                           00000770
      T=0.D0                                                            00000780
      DO820J=1,N                                                        00000790
      RNP=NP(J)                                                         00000800
      S=S+B(J)*DEXP(-D(J)*RK(I))*RK(I)**NP(J)*(2.D0*D(J))**(RNP+0.5D0)/ 00000810
     1DSQRT(FACTOR(2*NP(J)))                                            00000820
      T=T+C(J)*DEXP(-D(J)*RK(I))*RK(I)**NP(J)*(2.D0*D(J))**(RNP+0.5D0)/ 00000830
     1DSQRT(FACTOR(2*NP(J)))                                            00000840
  820 CONTINUE                                                          00000850
      P(I)=S                                                            00000860
      Q(I)=T                                                            00000870
      DIF(I)=DABS(P(I))-DABS(Q(I))                                      00000880
  810 CONTINUE                                                          00000890
C                                                                       00000900
C       MAXIMUM                                                         00000910
C                                                                       00000920
      DO850I=1,IM                                                       00000930
      P(I)=DABS(P(I))                                                   00000940
  850 Q(I)=DABS(Q(I))                                                   00000950
      DO500I=1,IM                                                       00000960
      INPAS=IM-I                                                        00000970
      IF((Q(INPAS)-Q(INPAS-1)).GT.0.D0.AND.Q(INPAS).GT.0.1D0)GOTO520    00000980
  500 CONTINUE                                                          00000990
  520 CONTINUE                                                          00001000
      T1=RK(INPAS-1)                                                    00001010
      T2=RK(INPAS)                                                      00001020
      T3=RK(INPAS+1)                                                    00001030
      R1=Q(INPAS-1)                                                     00001040
      R2=Q(INPAS)                                                       00001050
      R3=Q(INPAS+1)                                                     00001060
      TA=((R1-R2)/(T1-T2)-(R3-R2)/(T3-T2))/(T1-T3)                      00001070
      TB=-TA*(T1+T2)+(R1-R2)/(T1-T2)                                    00001080
      RMAX=-TB/(2.D0*TA)                                                00001090
      IRC=IM-INPAS                                                      00001100
C                                                                       00001110
C      RAYON DE COUPURE                                                 00001120
C                                                                       00001130
      IF(DIF(INPAS))250,250,260                                         00001140
  250 DO200I=IRC,IM                                                     00001150
      INPAS=IM-I                                                        00001160
      IF(DIF(INPAS))200,200,210                                         00001170
  200 CONTINUE                                                          00001180
  260 DO220I=IRC,IM                                                     00001190
      INPAS=IM-I                                                        00001200
      IF(DIF(INPAS))210,210,220                                         00001210
  220 CONTINUE                                                          00001220
  210 CONTINUE                                                          00001230
      S1=P(INPAS+1)-P(INPAS)-Q(INPAS+1)+Q(INPAS)                        00001240
      S2=(Q(INPAS)-P(INPAS))*(RK(INPAS+1)-RK(INPAS))                    00001250
      RC=(RK(INPAS)*S1+S2)/S1                                           00001260
  601 CONTINUE                                                          00001270
      R=RMAX                                                            00001280
      RPRIM=RC/2                                                        00001290
C                                                                       00001300
C     CALCUL DE A                                                       00001310
C                                                                       00001320
      S=0.D0                                                            00001330
      SS=0.D0                                                           00001340
      DO10I=1,N                                                         00001350
      S=S+C(I)**2*DK(NP(I),NP(I),D(I),D(I))*DI(2*NP(I),2.D0*D(I),R)     00001360
      MAX=I-1                                                           00001370
      IF(I.EQ.1)GOTO10                                                  00001380
      DO15J=1,MAX                                                       00001390
      SS=SS+C(I)*C(J)*DK(NP(I),NP(J),D(I),D(J))*DI(NP(I)+NP(J),D(I)+D(J)00001400
     1,R)                                                               00001410
 15   CONTINUE                                                          00001420
  10  CONTINUE                                                          00001430
      A1=S+2.D0*SS                                                      00001440
C                                                                       00001450
C   MODIF DU 21.11.78  CALCUL DE P(RC)                                  00001460
C                                                                       00001470
      RNNN=NNN                                                          00001480
      S=0.D0                                                            00001490
      DO700I=1,N                                                        00001500
      RNP=NP(I)                                                         00001510
      S=S+C(I)*(2.D0*D(I))**(RNP+0.5D0)*RC**RNP*DEXP(-D(I)*RC)/DSQRT(   00001520
     1FACTOR(2*NP(I)))                                                  00001530
  700 CONTINUE                                                          00001540
      PRC=S                                                             00001550
      IF(S.LT.0.D0)GOTO710                                              00001560
      GOTO730                                                           00001570
  710 PRC=-PRC                                                          00001580
      DO720I=1,N                                                        00001590
  720 C(I)=-C(I)                                                        00001600
      DO725I=1,NN                                                       00001610
  725 CP(I)=-CP(I)                                                      00001620
  730 CONTINUE                                                          00001630
      NM1=NNN-1                                                         00001640
      CONSC=0.6D0*RX*PRC/(RC*RPRIM**NM1)                                00001650
C                                                                       00001660
C      IMPRESSIONS DONNEES                                              00001670
C                                                                       00001680
      WRITE(6,1)                                                        00001690
  1   FORMAT(1H1)                                                       00001700
      WRITE(6,17)                                                       00001710
  17  FORMAT(5X,' PSEUDO ORBITALES ATOMIQUES',///)                      00001720
      WRITE(6,2)TITRE
  2   FORMAT(A80)        
      WRITE(6,3)                                                        00001750
  3   FORMAT(1X, //)                                                    00001760
      WRITE(6,4)N,NN,IM,ZAT                                             00001770
  4   FORMAT(5X,'N =',I2,3X,'NP =',I2,3X,'M =',I3,3X,'Z =', F6.2,//)    00001780
      WRITE(6,5)                                                        00001790
  5   FORMAT(' COEFFICIENTS DES ORBITALES ATOMIQUES',/)                 00001800
      WRITE(6,6)(NP(I),D(I),C(I),I=1,N)                                 00001810
    6 FORMAT(1X,I3,5X,D15.9,5X,D15.9)                                   00001820
      WRITE(6,3)                                                        00001830
      WRITE(6,7)                                                        00001840
  7   FORMAT(' COEFFICIENTS DE DEPART DES PSEUDO ORBITALES ATOMIQUES ') 00001850
      WRITE(6,6)(NNP(I),DP(I),CP(I),I=1,NN)                             00001870
      WRITE(6,3)                                                        00001880
      WRITE(6,8)RC,RMAX,R,RPRIM                                         00001890
  8   FORMAT( ' RC = ',D15.8,3X,' RMAX = ',D15.8,3X, ' R = ',D15.8,3X, '00001900
     1 RPRIM = ',D15.8//)                                               00001910
      WRITE(6,14)NNN,RMU,RX,RXETA                                       00001920
  14  FORMAT( ' NNN = ',I3,3X, ' MU = ',D15.8,3X,' RX = ',D15.8,3X, ' DZ00001930
     1ETA = ',D15.8)                                                    00001940
      WRITE(6,3)                                                        00001950
      WRITE(6,9)A1,CONSC                                                00001960
  9   FORMAT(5X, ' A = ',D15.9,3X, ' CONS.C = ',D15.9//)                00001970
      WRITE(6,3)                                                        00001980
C                                                                       00001990
      CALL CALCFX(NN,XPAR,FUNC)                                         00002000
      DO 20 I=1,NN                                                      00002010
      XU(I)=DP(I)                                                       00002020
  20  CONTINUE                                                          00002030
C                                                                       00002040
      WRITE(6,3)                                                        00002050
      WRITE(6,635)                                                      00002060
  635 FORMAT(3X,'K',9X,'RK',10X,'P(ORBITALES)',5X,'P(PSEUDO ORB.)')     00002070
      open(1,form='FORMATTED',status='UNKNOWN',
     *     file=titre(1:lo)//'.plot')
      write(1,1642)  titre
1642  format('#',a80)

C                                                                       00002080
C      CALCUL DE P1 ET DE PP                                            00002090
C                                                                       00002100
      DO610I=1,IM                                                       00002110
        S=0.D0
      DO620J=1,N                                                        00002130
      RNP=NP(J)                                                         00002140
      S=S+C(J)*DEXP(-D(J)*RK(I))*RK(I)**NP(J)*(2.D0*D(J))**(RNP+0.5D0)/ 00002150
     1DSQRT(FACTOR(2*NP(J)))                                            00002160
  620 CONTINUE                                                          00002170
      P1=S                                                              00002180
      S=0.D0                                                            00002190
      TP=0.D0                                                           00002200
      TS=0.D0                                                           00002210
      DO630J=1,NN                                                       00002220
      RNP=NNP(J)                                                        00002230
      PCO=X(J)*(2.D0*XU(J))**(RNP+0.5D0)/DSQRT(FACTOR(2*NNP(J)))        00002240
      S=S+PCO*DEXP(-XU(J)*RK(I))*RK(I)**NNP(J)                          00002250
      TP=TP+PCO*DEXP(-XU(J)*RK(I))*(NNP(J)*RK(I)**(NNP(J)-1)-XU(J)*     00002260
     1RK(I)**NNP(J))                                                    00002270
      TS=TS+PCO*DEXP(-XU(J)*RK(I))*(NNP(J)*(NNP(J)-1)*RK(I)**(NNP(J)-2) 00002280
     1-2.D0*NNP(J)*XU(J)*RK(I)**(NNP(J)-1)+XU(J)**2*RK(I)**NNP(J))      00002290
  630 CONTINUE                                                          00002300
      PP=S                                                              00002310
      WRITE(6,640)I,RK(I),P1,PP,TP,TS                                   00002320
      WRITE(1,1640)RK(I),P1,PP
  640 FORMAT(2X,I3,5(2X,D15.9))                                         00002330
 1640 FORMAT(5(2X,D15.9))        
      XK(I)=P1                                                          00002340
       W(I)=PP                                                          00002350
  610 CONTINUE                                                          00002360
C         CALCUL DES VALEURS MOYENNES DE R**P                           00002370
      WRITE(6,3)                                                        00002380
      WRITE(6,11)                                                       00002390
   11 FORMAT(5X,'VALEURS MOYENNES DE R**P',/)                           00002400
      WRITE(6,12)                                                       00002410
   12 FORMAT(5X,'P',10X,'ORBITALE',9X,'PSEUDO-ORBITALE')                00002420
      DO 50 IP=1,7                                                      00002430
      IPP=IP-3                                                          00002440
      S=0.D0                                                            00002450
      SS=0.D0                                                           00002460
      DO 51 I=1,N                                                       00002470
      S=S+DK(NP(I),NP(I),D(I),D(I))*DJ(2*NP(I)+IPP,2.D0*D(I))*C(I)**2   00002480
      M=I-1                                                             00002490
      IF(I.EQ.1) GO TO 51                                               00002500
      DO 52 J=1,M                                                       00002510
      SS=SS+C(I)*C(J)*DK(NP(I),NP(J),D(I),D(J))*DJ(NP(I)+NP(J)+IPP,D(I)+00002520
     1D(J))                                                             00002530
   52 CONTINUE                                                          00002540
   51 CONTINUE                                                          00002550
      VMR=S+2.D0*SS                                                     00002560
      S=0.D0                                                            00002570
      SS=0.D0                                                           00002580
      DO 53 I=1,NN                                                      00002590
      S=S+DK(NNP(I),NNP(I),XU(I),XU(I))*DJ(2*NNP(I)+IPP,2.D0*XU(I))*X(I)00002600
     1**2                                                               00002610
      M=I-1                                                             00002620
      IF(I.EQ.1) GO TO 53                                               00002630
      DO 54 J=1,M                                                       00002640
      SS=SS+X(I)*X(J)*DK(NNP(I),NNP(J),XU(I),XU(J))*DJ(NNP(I)+NNP(J)+IPP00002650
     1,XU(I)+XU(J))                                                     00002660
   54 CONTINUE                                                          00002670
   53 CONTINUE                                                          00002680
      VMRP=S+2.D0*SS                                                    00002690
      WRITE(6,13) IPP,VMR,VMRP                                          00002700
   13 FORMAT(4X,I2,2(6X,D15.8))                                         00002710
   50 CONTINUE                                                          00002720
 2000 STOP                                                              00002740
      END                                                               00002750
      SUBROUTINE CALCFX(NPAR,XPAR,RNORU)                                00002760
      IMPLICIT REAL*8 (A-H,O-Z)                                         00002770
      DIMENSION A(30,30),CC(30,30),Y(30),B(30,30),ARES(30,30),          00002780
     1F(100),RLAMB(100),DP(30)                                          00002790
      DIMENSION  P(30,30),AU(30,30),A2(30,30),Y1(30),Y2(30)             00002800
      DIMENSION XPAR(30)                                                00002810
      COMMON/PSINF/D(30),C(30),CP(30),X(30),R,TEST,A1,NP(30),NNP(30),N, 00002820
     1NN                                                                00002830
      COMMON/DEVI/ RPRIM,RMU,RXETA,RX  ,CONSC,NNN                       00002840
      COMMON/OPTI/DP                                                    00002850
      COMMON /FAC/ FACTOR(16)                                           00002860
C                                                                       00002870
C      CALCUL DES COEFFICIENTS A(I,J)                                   00002880
C                                                                       00002890
      DO20I=1,NN                                                        00002900
      DO20J=1,I                                                         00002910
      P(I,J)=DK(NNP(I),NNP(J),DP(I),DP(J))                              00002920
      AU(I,J)=P(I,J)*DI(NNP(I)+NNP(J),DP(I)+DP(J),R)                    00002930
      A2(I,J)=RMU*P(I,J)*(DJ(NNP(I)+NNP(J),DP(I)+DP(J))-DI(NNP(I)+NNP(J)00002940
     1,DP(I)+DP(J),RPRIM))                                              00002950
      A(I,J)=AU(I,J)+A2(I,J)                                            00002960
      A(J,I)=A(I,J)                                                     00002970
  20  CONTINUE                                                          00002980
C                                                                       00002990
C     CALCUL DES COEFFICIENTS C(I,J)                                    00003000
C                                                                       00003010
      DO25I=1,NN                                                        00003020
      DO25J=1,N                                                         00003030
      CC(I,J)=DK(NNP(I),NP(J),DP(I),D(J))*DI(NNP(I)+NP(J),DP(I)+D(J),R) 00003040
  25  CONTINUE                                                          00003050
C                                                                       00003060
C     CALCUL DES COEFFICIENTS Y(I)                                      00003070
C                                                                       00003080
      DO30I=1,NN                                                        00003090
      S=0.D0                                                            00003100
      DO35J=1,N                                                         00003110
      S=S+CC(I,J)*C(J)                                                  00003120
  35  CONTINUE                                                          00003130
      Y1(I)=S                                                           00003140
      RNNP=NNP(I)                                                       00003150
      DKMOD=(2.D0*DP(I))**(RNNP+0.5D0)/DSQRT(FACTOR(2*NNP(I)))          00003160
      Y2(I)=RMU*CONSC*   DKMOD                  *(DJ(NNP(I)+NNN,DP(I)+  00003170
     1RXETA)-DI(NNP(I)+NNN,DP(I)+RXETA,RPRIM))                          00003180
      Y(I)=Y1(I)+Y2(I)                                                  00003190
  30  CONTINUE                                                          00003200
C                                                                       00003210
C     CALCUL DES COEFFICIENTS B(I,J)                                    00003220
C                                                                       00003230
      DO40I=1,NN                                                        00003240
      DO40J=1,I                                                         00003250
      B(I,J)=DK(NNP(I),NNP(J),DP(I),DP(J))*DJ(NNP(I)+NNP(J),DP(I)+DP(J))00003260
      B(J,I)=B(I,J)                                                     00003270
  40  CONTINUE                                                          00003280
C                                                                       00003290
C      CALCUL DE LAMBDA(0) ET DE F(0)                                   00003300
C                                                                       00003310
      DO45I=1,NN                                                        00003320
      X(I)=Y(I)                                                         00003330
      DO45J=1,NN                                                        00003340
      ARES(I,J)=A(I,J)                                                  00003350
  45  CONTINUE                                                          00003360
      CALL DRESYL(ARES,30,30,X,30,1,NN,1,1.D-15,KK)                     00003370
      S=0.D0                                                            00003380
      SS=0.D0                                                           00003390
      DO50I=1,NN                                                        00003400
      S=S+X(I)**2*B(I,I)                                                00003410
      MAX=I-1                                                           00003420
      IF(I.EQ.1)GOTO50                                                  00003430
      DO55J=1,MAX                                                       00003440
      SS=SS+X(I)*X(J)*B(I,J)                                            00003450
  55  CONTINUE                                                          00003460
  50  CONTINUE                                                          00003470
      F(1)=-1.D0+S+2.D0*SS                                              00003480
      RLAMB(1)=0.D0                                                     00003490
      K=1                                                               00003500
C                                                                       00003510
C     CALCUL DE LAMBDA(1)                                               00003520
C                                                                       00003530
      RLAMB(2)=1.D-5                                                    00003540
  300 K=K+1                                                             00003550
C                                                                       00003560
C      CALCUL DU MEILLEUR LAMBDA                                        00003570
C                                                                       00003580
      DO75I=1,NN                                                        00003590
      X(I)=Y(I)                                                         00003600
      DO75J=1,NN                                                        00003610
      ARES(I,J)=A(I,J)-RLAMB(K)*B(I,J)                                  00003620
  75  CONTINUE                                                          00003630
      CALL DRESYL(ARES,30,30,X,30,1,NN,1,1.D-15,KK)                     00003640
C                                                                       00003650
C      CALCUL DE F(K)                                                   00003660
C                                                                       00003670
      S=0.D0                                                            00003680
      SS=0.D0                                                           00003690
      DO80I=1,NN                                                        00003700
      S=S+X(I)**2*B(I,I)                                                00003710
      MAX=I-1                                                           00003720
      IF(I.EQ.1)GOTO80                                                  00003730
      DO85J=1,MAX                                                       00003740
      SS=SS+X(I)*X(J)*B(I,J)                                            00003750
  85  CONTINUE                                                          00003760
  80  CONTINUE                                                          00003770
      F(K)=-1.D0+S+2.D0*SS                                              00003780
      IF(DABS(F(K)).LE.TEST)GOTO200                                     00003790
      IF(K.GE.50)GOTO200                                                00003800
      RLAMB(K+1)=RLAMB(K)-F(K)*(RLAMB(K)-RLAMB(K-1))/(F(K)-F(K-1))      00003810
      GOTO300                                                           00003820
  200 CONTINUE                                                          00003830
C                                                                       00003840
C      CALCUL DE FONCTION                                               00003850
C                                                                       00003860
      S=0.D0                                                            00003870
      SS=0.D0                                                           00003880
      SSS=0.D0                                                          00003890
      SD=0.D0                                                           00003900
      SSD=0.D0                                                          00003910
      DO90I=1,NN                                                        00003920
      SD=SD+X(I)*(X(I)*AU(I,I)-2.D0*Y1(I))                              00003930
      S=S+X(I)*(X(I)*A(I,I)-2.D0*Y(I))                                  00003940
      SSS=SSS+X(I)*Y(I)                                                 00003950
      MAX=I-1                                                           00003960
      IF(I.EQ.1)GOTO90                                                  00003970
      DO91J=1,MAX                                                       00003980
      SS=SS+X(I)*X(J)*A(I,J)                                            00003990
      SSD=SSD+X(I)*X(J)*AU(I,J)                                         00004000
  91  CONTINUE                                                          00004010
  90  CONTINUE                                                          00004020
      IF(RXETA.EQ.0.D0)GOTO92                                           00004030
      STZ=CONSC**2*    1.D0               *(DJ(2*NNN,2.D0*RXETA)-DI(2*  00004040
     1NNN,2.D0*RXETA,RPRIM))                                            00004050
      GOTO93                                                            00004060
  92  STZ=CONSC**2*1.D0*RPRIM**(2*NNN+1)/(2*NNN+1)                      00004070
  93  CONTINUE                                                          00004080
      FUNC1=A1+RMU*STZ+S+2.D0*SS-RLAMB(K)                               00004090
      FUNC2=A1+RMU*STZ-SSS                                              00004100
      RNORT=DSQRT(FUNC2+RLAMB(K))                                       00004110
      RNORU=DSQRT(A1+SD+2.D0*SSD)                                       00004120
      WRITE(6,13)(DP(I),I=1,NN)                                         00004130
   13 FORMAT(1x,/,' EXP. ',6D15.8)                                      00004140
      WRITE(6,16)(X(I),I=1,NN)                                          00004150
   16 FORMAT(' COEF ',6D15.8)                                           00004160
      WRITE(6,501)RLAMB(K)                                              00004170
  501 FORMAT(' LAMBDA ',D15.8)                                          00004180
      WRITE(6,502)FUNC1,FUNC2                                               
  502 FORMAT( ' FONC 1 ',D15.8, 'FONC 2 ',D15.8)
      WRITE(6,503)RNORT,RNORU
  503 FORMAT( ' NORME TOTALE ',D15.8, 'NORME UTILE ',D15.8)
      RETURN                                                            00004220
      END                                                               00004230
      DOUBLE PRECISION FUNCTION DI(M,A,B)                               00004240
      DOUBLE PRECISION A,B,S,FACTOR,DEXP                                00004250
      COMMON /FAC/ FACTOR(16)                                           00004260
      S=0.D0                                                            00004270
      DO5J=1,M                                                          00004280
      I=J-1                                                             00004290
      MR=M-I                                                            00004300
  5   S=S+FACTOR(M)*B**MR/(FACTOR(MR)*A**J)                             00004310
      DI=(S+FACTOR(M)/A**(M+1))*DEXP(-A*B)                              00004320
      RETURN                                                            00004330
      END                                                               00004340
      DOUBLE PRECISION FUNCTION DJ(M,A)                                 00004350
      DOUBLE PRECISION A,FACTOR                                         00004360
      COMMON /FAC/ FACTOR(16)                                           00004370
      IF(M.EQ.0)GOTO10                                                  00004380
      DJ=FACTOR(M)/A**(M+1)                                             00004390
      GOTO20                                                            00004400
  10  DJ=1.D0/A                                                         00004410
  20  CONTINUE                                                          00004420
      RETURN                                                            00004430
      END                                                               00004440
      DOUBLE PRECISION FUNCTION DK(M,N,A,B)                             00004450
      DOUBLE PRECISION A,B,FACTOR,S,D,RM,RN,DSQRT                       00004460
      COMMON /FAC/ FACTOR(16)                                           00004470
      RM=M                                                              00004480
      RN=N                                                              00004490
      S=(2.D0*A)**(RM+0.5D0)*(2.D0*B)**(RN+0.5D0)                       00004500
      D=DSQRT(FACTOR(2*M)*FACTOR(2*N))                                  00004510
      DK=S/D                                                            00004520
      RETURN                                                            00004530
      END                                                               00004540
      BLOCK DATA                                                        00004550
      DOUBLE PRECISION FACTOR(16)                                       00004560
      COMMON /FAC/ FACTOR                                               00004570
      DATA FACTOR/                  .10000000000D+01, .20000000000D+01, 00004580
     1            .60000000000D+01, .24000000000D+02, .12000000000D+03, 00004590
     2            .72000000000D+03, .50400000000D+04, .40320000000D+05, 00004600
     3            .36288000000D+06, .36288000000D+07, .39916800000D+08, 00004610
     4            .47900160000D+09, .62270208000D+10, .87178291200D+11, 00004620
     5            .13076743680D+13, .20922789888D+14/                   00004630
      END                                                               00004640
      SUBROUTINE DRESYL (A,NA,MA,B,NB,MB,N,M,TEST,K1)                   00004650
      DOUBLE PRECISION A(NA,MA),B(NB,MB),TEST,DABS,ALJ                  00004660
      DOUBLE PRECISION BLJ                                              00004670
       NP1=N+1                                                          00004680
       K1=0                                                             00004690
       DO 1 L=1,N                                                       00004700
       LP1=L+1                                                          00004710
      IF(DABS(A(L,L)).GT.TEST) GO TO 2                                  00004720
C                                                                       00004730
C CAS OU LE PIVOT EST TROP PETIT                                        00004740
C                                                                       00004750
       I=L                                                              00004760
 5     I=I+1                                                            00004770
       IF(I.GT.N) GO TO 101                                             00004780
      IF(DABS(A(I,L)).LT.TEST) GO TO 5                                  00004790
        GO TO 3                                                         00004800
C                                                                       00004810
C CAS OU LA RESOLUTION EST IMPOSSIBLE                                   00004820
C                                                                       00004830
 101   K1=1                                                             00004840
       WRITE(6,*) ' ERREUR DANS LE PROGRAMME DRESYL'
       WRITE(6,7)TEST   
7      FORMAT(2X,'IL N EXISTE PAS DE PIVOT SUPERIEUR A',D24.17)
       write(6,*) 'LE DETERMINANT EST NUL A LA PRECISION DEMANDEE'
       go to 11
   3  DO 4 J=1,N                                                        00004900
       ALJ=A(L,J)                                                       00004910
      A(L,J)=A(I,J)                                                     00004920
  4   A(I,J)=ALJ                                                        00004930
       DO 50 J=1,M                                                      00004940
       BLJ=B(L,J)                                                       00004950
       B(L,J)=B(I,J)                                                    00004960
  50  B(I,J)=BLJ                                                        00004970
 2     IF(L.EQ.N) GO TO 100                                             00004980
       DO 6 J=LP1,N                                                     00004990
   6   A(L,J)=A(L,J)/A(L,L)                                             00005000
 100    DO 8 J=1,M                                                      00005010
   8   B(L,J)=B(L,J)/A(L,L)                                             00005020
       IF(L.EQ.N) GO TO 1                                               00005030
       DO 19 I=LP1,N                                                    00005040
       DO 70 J=LP1,N                                                    00005050
  70  A(I,J)=A(I,J)-A(I,L)*A(L,J)                                       00005060
       DO 80 J=1,M                                                      00005070
  80  B(I,J)=B(I,J)-A(I,L)*B(L,J)                                       00005080
 19     CONTINUE                                                        00005090
  1    CONTINUE                                                         00005100
       I=N                                                              00005110
 10    I=I-1                                                            00005120
c  
c   fin
c
       IF(I.EQ.0) GO TO 11                                              00005130
c
       IP1=I+1                                                          00005140
       DO 12 K=IP1,N                                                    00005150
        DO 12 J=1,M                                                     00005160
  12  B(I,J)=B(I,J)-A(I,K)*B(K,J)                                       00005170
       GO TO 10                                                         00005180
  11  RETURN                                                            00004890
      END                                                               00005190
