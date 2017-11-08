      PROGRAM PSDF87                                                    PSD0003
      IMPLICIT REAL*8 (A-H,O-Z)                                         PSD00040
      DIMENSION SFI(22,22),AL(22),D(22,22),A(60),H(22)                  PSD00050
      COMMON/SUP/FSUP(22,10),ESUP(10),VHSUPV(10),R000,NSUP              PSD00060
      DIMENSION XU(6),EPREC(6),W(100)                                   PSD00070
      COMMON W                                                          PSD00080
      NAMELIST/PSEUDO/PRECIS,ESCALE,MAXIT,NETAT,YRSTRT                  PSD00090
     &,MVAR,PASVAR,ORIVAR,POIDS                                         PSD00100
     &,MAXV,ISTOP,DMNXP,YNLOC                                           PSD00110
     &,DMNXP                                                            PSD00120
     &,CSPO,NGAM,ESO,PHISO                                              PSD00130
C      EXTERNAL INT,ININT
      LOGICAL YRSTRT,YNLOC                                              PSD00140
      COMMON/SPOC/ CSPO(10),ESO(2),PHISO(22,3),NGAM                     PSD00150
      COMMON/RSVC/FATT(30),FACTO(16),PQN(60),ZETA(60),AJMN(15),         PSD00160
     1OCCUCS(4),OCCUP(4),C(22,8),FC(500),FO(500),S(500),U(500),T(500),  PSD00170
     2P(1627),Q(1627),DT(500),DOS(500),ZN,NBAS(4),NCSH(4),NOSH(4),NSYM, PSD00180
     3NFLAG1,NOSHT,NBAST,NFI                                            PSD00190
      COMMON/PSEUD/ALFCT(30),CACT(30),ALFA(30,4),CA(30,4),PSE(4),       PSD00200
     1ICT(30),ISG(30,4),NT(4),NTCT,NCT(30),NX(30,4)                     PSD00210
      COMMON/UN/TAF(22,22,3),APQ(22,22,3),ZO(22,3),AO(22,3),PSEO(3),    PSD00220
     1CPSO(22,3),NO(22,3),HV(22,3),VHV(3),L,NBAS1,NETAT,NEND            PSD00230
      DIMENSION NFIL(4)                                                 PSD00240
      DATA NFIL/1,2,3,4/                                                PSD00250
       NEND=1                                                           PSD00260
      call openf
      YRSTRT=.FALSE.                                                    PSD00300
      YNLOC=.FALSE.                                                     PSD00310
      DMNXP=1.D-8                                                       PSD00320
      PRECIS=1.D-5                                                      PSD00330
      DMNXP=1.D-8                                                       PSD00340
      ESCALE=1.D05                                                      PSD00350
      MAXIT=30                                                          PSD00360
      NETAT=1                                                           PSD00370
      CSPO(1)=0.                                                        PSD00380
      CSPO(2)=0.                                                        PSD00390
      CSPO(3)=0.                                                        PSD00400
      MVAR=1                                                            PSD00410
      PASVAR=0.1                                                        PSD00420
      ORIVAR=0.1                                                        PSD00430
      MAXV=1                                                            PSD00440
      ISTOP=0                                                           PSD00450
      NSUP=0                                                            PSD00460
      R000=0.D0                                                         PSD00470
      READ(5,PSEUDO)                                                    PSD00480
      WRITE(6,PSEUDO)                                                   PSD00490
      KET=1                                                             PSD00500
      FACTO(1)=1.                                                       PSD00510
      KVAR=MVAR                                                         PSD00520
      FACTO(2)=1.                                                       PSD00530
      FATT(1)=1.                                                        PSD00540
      FATT(2)=1.                                                        PSD00550
      DO 1301 I=3,16                                                    PSD00560
      IM1=I-1                                                           PSD00570
      FACTO(I)=IM1*FACTO(I-2)                                           PSD00580
 1301  FATT(I)=IM1*FATT(IM1)                                            PSD00590
      DO 1201 I=17,30                                                   PSD00600
      IM1=I-1                                                           PSD00610
 1201 FATT(I)=IM1*FATT(IM1)                                             PSD00620
  255 CALL INPUT                                                        PSD00630
      REWIND 1                                                          PSD00640
      DO 2 I=1,NBAST                                                    PSD00650
      N1=PQN(I)                                                         PSD00660
      DEP=PURX(N1,N1,ZETA(I),ZETA(I),0,0.D0)                            PSD00670
      A(I)=DSQRT(1.D0/DEP)                                              PSD00680
    2 CONTINUE                                                          PSD00690
      IST=0                                                             PSD00700
      DO 50 L=1,NSYM                                                    PSD00710
      NBAS1=NBAS(L)                                                     PSD00720
      CNOR=0.D0                                                         PSD00730
      DO 51 I=1,NBAS1                                                   PSD00740
      IPZA=IST+I                                                        PSD00750
      N1=PQN(IPZA)                                                      PSD00760
      Z1=ZETA(IPZA)                                                     PSD00770
      DO 51 J=1,NBAS1                                                   PSD00780
      JPZA=IST+J                                                        PSD00790
      N2=PQN(JPZA)                                                      PSD00800
      Z2=ZETA(JPZA)                                                     PSD00810
      COEF=A(IPZA)*A(JPZA)*C(I,L)*C(J,L)                                PSD00820
   51 CNOR=CNOR+COEF*PURX(N1,N2,Z1,Z2,0,0.D0)                           PSD00830
      DO 52 II=1,NBAS1                                                  PSD00840
      C(II,L)=C(II,L)/DSQRT(CNOR)                                       PSD00850
   52 CONTINUE                                                          PSD00860
      WRITE(6,61) (C(II,L),II=1,NBAS1)                                  PSD00870
61    FORMAT('0ORBITALE RENORMALISEE:', 9F11.6,/(23X, 9F11.6))          PSD00880
   50 IST=IST+NBAS1                                                     PSD00890
C                                                                       PSD00900
       CALL OEISG                                                       PSD00910
      NFI=NFIL(KET)                                                     PSD00920
      IF(.NOT.YRSTRT) CALL TEIPSG                                       PSD00930
      CALL DENSI                                                        PSD00940
      CALL HAMIL                                                        PSD00950
      MVAR=KVAR                                                         PSD00960
      NSTEP=0                                                           PSD00970
      NB1=0                                                             PSD00980
      DO 10 L=1,NSYM                                                    PSD00990
      NTL=NT(L)                                                         PSD01000
      NBAS1=NBAS(L)                                                     PSD01010
      NCSHL=NCSH(L)                                                     PSD01020
      IF(NTL.NE.0) GO TO 100                                            PSD01030
      NSTEP=NSTEP+NBAS1*(NBAS1+1)/2                                     PSD01040
      NB1=NB1+NBAS1                                                     PSD01050
   10 CONTINUE                                                          PSD01060
  100 WRITE(6,20) L                                                     PSD01070
      IF(CSPO(KET).EQ.0) GO TO 57                                       PSD01080
      CNOR=0.                                                           PSD01090
      DO 58 I=1,NBAS1                                                   PSD01100
      II=I+NB1                                                          PSD01110
      N1=PQN(II)                                                        PSD01120
      DO 58 J=1,NBAS1                                                   PSD01130
      JJ=J+NB1                                                          PSD01140
      N2=PQN(JJ)                                                        PSD01150
58    CNOR=CNOR+PHISO(I,KET)*PHISO(J,KET)*A(II)*A(JJ)                   PSD01160
     &*PURX(N1,N2,ZETA(II),ZETA(JJ),0,0.D0)                             PSD01170
      CNOR=1.D0/DSQRT(CNOR)                                             PSD01180
      DO 53 I=1,NBAS1                                                   PSD01190
53    PHISO(I,KET)=PHISO(I,KET)*CNOR                                    PSD01200
      WRITE(6,56) (PHISO(I,KET),I=1,NBAS1)                              PSD01210
56    FORMAT('0,PHISO NORMALISEE',9F11.6,/,(17X,9F11.6))                PSD01220
      GO TO 60                                                          PSD01230
57    DO 59 I=1,NBAS1                                                   PSD01240
59    PHISO(I,KET)=C(I,L)                                               PSD01250
60    CONTINUE                                                          PSD01260
   20 FORMAT(1H1,50X,'SYMMETRY',I2,///)                                 PSD01270
      DO 55 I=1,NBAS1                                                   PSD01280
   55 CPSO(I,KET)=PHISO(I,KET)                                          PSD01290
      PSEO(KET)=PSE(L)                                                  PSD01300
      DO 11 I=1,NBAS1                                                   PSD01310
      DO 11 J=1,I                                                       PSD01320
      K=NSTEP+(I*(I-1)/2)+J                                             PSD01330
      SFI(I,J)=S(K)                                                     PSD01340
      SFI(J,I)=SFI(I,J)                                                 PSD01350
      IF(NCSHL.EQ.0) GO TO 12                                           PSD01360
      TAF(I,J,KET)=FC(K)                                                PSD01370
      GO TO 11                                                          PSD01380
   12 TAF(I,J,KET)=FO(K)                                                PSD01390
11    FC(K-NSTEP)=T(K)-ZN*U(K)                                          PSD01400
      CALL DVLVC1(SFI,22,22,AL,22,D,22,22,NBAS1,II)                     PSD01410
      DO 805 IP=1,NBAS1                                                 PSD01420
      IF(AL(IP).LT.DMNXP) GO TO 805                                     PSD01430
      ALIP=AL(IP)**(-.25)                                               PSD01440
      DO 804 IQ=1,NBAS1                                                 PSD01450
804   D(IP,IQ)=D(IP,IQ)*ALIP                                            PSD01460
805   CONTINUE                                                          PSD01470
      DO 17 IP=1,NBAS1                                                  PSD01480
      DO 17 IQ=1,IP                                                     PSD01490
      SPQ=0.D0                                                          PSD01500
      DO 18 IR=1,NBAS1                                                  PSD01510
18    SPQ=SPQ+D(IR,IP)*D(IR,IQ)                                         PSD01520
      APQ(IP,IQ,KET)=SPQ                                                PSD01530
      APQ(IQ,IP,KET)=APQ(IP,IQ,KET)                                     PSD01540
   17 CONTINUE                                                          PSD01550
C                                                                       PSD01560
      NB1=1                                                             PSD01570
      DO 1 LL=1,NSYM                                                    PSD01580
      IF(LL.GE.L) GO TO 1                                               PSD01590
      NB1=NB1+NBAS(LL)                                                  PSD01600
    1 CONTINUE                                                          PSD01610
      NB2=NB1+NBAS1-1                                                   PSD01620
      NSB=NB1-1                                                         PSD01630
      DO 4 I=NB1,NB2                                                    PSD01640
      NEW=I-NSB                                                         PSD01650
      ZO(NEW,KET)=ZETA(I)                                               PSD01660
      NO(NEW,KET)=PQN(I)                                                PSD01670
    4 AO(NEW,KET)=A(I)                                                  PSD01680
C                                                                       PSD01690
      IF(NTCT.EQ.0) GO TO 6                                             PSD01700
      DO 7 I=1,NBAS1                                                    PSD01710
      DO 7 J=1,I                                                        PSD01720
      SOM=0.D0                                                          PSD01730
      DO 8 JCT=1,NTCT                                                   PSD01740
      IX=ICT(JCT)                                                       PSD01750
      KX=NCT(JCT)                                                       PSD01760
      AX=ALFCT(JCT)                                                     PSD01770
    8 SOM=SOM+CACT(JCT)*WPS(ZO(I,KET),AO(I,KET),NO(I,KET),ZO(J,KET),    PSD01780
     1AO(J,KET),NO(J,KET),IX,KX,AX)                                     PSD01790
      TAF(I,J,KET)=TAF(I,J,KET)+SOM                                     PSD01800
    7 CONTINUE                                                          PSD01810
6     IF(MVAR.EQ.1) GO TO 115                                           PSD01820
      ALFA(1,L)=ORIVAR                                                  PSD01830
      DO 114 I=2,NTL                                                    PSD01840
114   ALFA(I,L)=ALFA(1,L)                                               PSD01850
115   DNOR=0.1D+30
      DO 210 I=1,NBAS1                                                  PSD01880
      DO 210 J=1,I                                                      PSD01890
  210 TAF(J,I,KET)=TAF(I,J,KET)                                         PSD01900
      SV=0.D0                                                           PSD01910
      DO 200 IP=1,NBAS1                                                 PSD01920
      SOM=0.D0                                                          PSD01930
      DO 201 IQ=1,NBAS1                                                 PSD01940
  201 SOM=SOM+CPSO(IQ,KET)*TAF(IP,IQ,KET)                               PSD01950
      H(IP)=SOM                                                         PSD01960
200   SV=SV+CPSO(IP,KET)*SOM                                            PSD01970
      VHV(KET)=SV                                                       PSD01980
      WRITE(6,216) VHV(KET)                                             PSD01990
216   FORMAT('0VALEUR MOYENNE SUR H MOYEN',F20.6)                       PSD02000
      DO 202 IP=1,NBAS1                                                 PSD02010
      SOM=0.D0                                                          PSD02020
      DO 203 IQ=1,NBAS1                                                 PSD02030
  203 SOM=SOM+APQ(IP,IQ,KET)*H(IQ)                                      PSD02040
  202 HV(IP,KET)=SOM                                                    PSD02050
      IF(NSUP.EQ.0) GO TO 285                                           PSD02060
      DO 280 I=1,NSUP                                                   PSD02070
      TOT=0.                                                            PSD02080
      DO 270 IP=1,NBAS1                                                 PSD02090
      SOM=0.                                                            PSD02100
      DO 260 IQ=1,NBAS1                                                 PSD02110
  260 SOM=SOM+FSUP(IQ,I)*TAF(IQ,IP,KET)                                 PSD02120
  270 TOT=TOT+SOM*FSUP(IP,I)                                            PSD02130
  280 VHSUPV(I)=TOT                                                     PSD02140
  285 CONTINUE                                                          PSD02150
      KET=KET+1                                                         PSD02160
      IF(KET.LE.NETAT) GO TO 255                                        PSD02170
      IF(YNLOC) CALL EXTRAC                                             PSD02180
116   CALL ICALCF                                                       PSD02190
      IF(ISTOP.GE.3) GO TO 610                                          PSD02200
      DO 120 NVAR=1,MVAR                                                PSD02210
      NU=1                                                              PSD02220
      XU(1)=ALFA(1,L)                                                   PSD02230
      EPREC(1)=PRECIS*XU(1)                                             PSD02240
      IF(NTL.LT.2) GO TO 1114                                           PSD02250
      DO 14 I=2,NTL                                                     PSD02260
      IF(ALFA(I,L).EQ.ALFA(I-1, L)) GO TO 14                            PSD02270
      NU=NU+1                                                           PSD02280
      XU(NU)=ALFA(I,L)                                                  PSD02290
      EPREC(NU)=PRECIS*XU(NU)                                           PSD02300
   14 CONTINUE                                                          PSD02310
1114  CALL CALCFX(NU,XU,GAM)                                            PSD02320
      NEND=0                                                            PSD02330
      IF(DNOR.LT.GAM) GO TO 117                                         PSD02340
      DNOR=GAM                                                          PSD02350
      ALNOR=ALFA(1,L)                                                   PSD02360
117   DO 120 I=1,NTL                                                    PSD02370
120   ALFA(I,L)=ALFA(I,L)+PASVAR                                        PSD02380
      IF(MVAR.EQ.1) GO TO 125                                           PSD02390
      DO 121 I=1,NTL                                                    PSD02400
121   ALFA(I,L)=ALNOR                                                   PSD02410
      IF(DNOR.GT.2.) GO TO 110                                          PSD02420
      IF(ISTOP.EQ.2) GO TO 110                                          PSD02430
      MVAR=1                                                            PSD02440
      GO TO 116                                                         PSD02450
125   CONTINUE                                                          PSD02460
      IF(ISTOP.EQ.2) GO TO 110                                          PSD02470
       NEND=0                                                           PSD02480
      CALL DVA0XX(XU,EPREC,NU,GAM,ESCALE,1,1,MAXIT)                     PSD02490
       NEND=1                                                           PSD02500
      CALL CALCFX(NU,XU,GAM)                                            PSD02510
      IF(ISTOP.EQ.1) GO TO 110                                          PSD02520
      NEND=0                                                            PSD02530
      DO 510 I=1,NTL                                                    PSD02540
      XU(I)=XU(I+1)                                                     PSD02550
510   EPREC(I)=PRECIS*XU(I)                                             PSD02560
      NU=NTL                                                            PSD02570
      CALL DVA0XX(XU,EPREC,NU,GAM,ESCALE,1,1,MAXV)                      PSD02580
      NEND=1                                                            PSD02590
      CALL CALCFX(NU,XU,GAM)                                            PSD02600
      IF(ISTOP.EQ.0) GO TO 110                                          PSD02610
610   J=NTL+1                                                           PSD02620
      XU(J)=ALFA(1,L)                                                   PSD02630
      EPREC(J)=XU(J)*PRECIS                                             PSD02640
      NU=0                                                              PSD02650
      DO 620 I=1,NTL                                                    PSD02660
      NU=NU+1                                                           PSD02670
      XU(NU)=CA(I,L)                                                    PSD02680
      IF(ALFA(I,L).EQ.XU(J)) GO TO 620                                  PSD02690
      J=J+1                                                             PSD02700
      XU(J)=ALFA(I,L)                                                   PSD02710
      EPREC(J)=XU(J)*PRECIS                                             PSD02720
620   EPREC(NU)=XU(NU)*PRECIS                                           PSD02730
      NU=J                                                              PSD02740
      CALL CALCFX(NU,XU,GAM)                                            PSD02750
      IF(GAM.GT.10.) GO TO 110                                          PSD02760
      NEND=0                                                            PSD02770
      CALL DVA0XX(XU,EPREC,NU,GAM,ESCALE,1,1,MAXIT)                     PSD02780
      NEND=1                                                            PSD02790
      CALL CALCFX(NU,XU,GAM)                                            PSD02800
110   STOP                                                              PSD02810
      END                                                               PSD02820
      SUBROUTINE INPUT                                                  PSD02830
      IMPLICIT REAL*8 (A-H,O-Z)                                         PSD02840
      CHARACTER*8 SLATER,IBLA,NAME,PSDTYP,PSDCST                        PSD02850
      COMMON/SUP/FSUP(22,10),ESUP(10),VHSUPV(10),R000,NSUP            
      DIMENSION NCCUP(4),NAME(18)                                       PSD02870
      COMMON/RSVC/FATT(30),FACTO(16),PQN(60),ZETA(60),AJMN(15),         PSD02880
     1OCCUCS(4),OCCUP(4),C(22,8),FC(500),FO(500),S(500),U(500),T(500),  PSD02890
     2P(1627),Q(1627),DT(500),DOS(500),ZN,NBAS(4),NCSH(4),NOSH(4),NSYM, PSD02900
     3NFLAG1,NOSHT,NBAST                                                PSD02910
      LOGICAL YW,YPSD                                                   PSD02920
      DIMENSIONNBPOT(4),APOT(10),CPOT(10),NPOT(10)                      PSD02930
      COMMON/PSEUD/ALFCT(30),CACT(30),ALFA(30,4),CA(30,4),PSE(4),       PSD02940
     1ICT(30),ISG(30,4),NT(4),NTCT,NCT(30),NX(30,4)                     PSD02950
      DIMENSION TITLE(9),VCC(15),NOOC(4),NFLAG(4)                       PSD02960
      EQUIVALENCE (TITLE(1),NAME(1)),(VCC(1),AJMN(1)),(NOOC(1),NCCUP(1))PSD02970
     &,(NFLAG1,NFLAG(1))                                                PSD02980
      NAMELIST/ATOM/ TITLE,CHARGE,NBAS,NCSH,NOOC,PQN,ZETA,VCC,NFLAG     PSD02990
     &,PSE,YPSD,C,NSUP,ESUP,FSUP,PSDTYP,NTCT,CACT,ALFCT,PSDCST,         PSD03000
     &NBPOT,NPOT,APOT,CPOT,NCT                                          PSD03010
      IBLA='    '                                                       PSD03020
      SLATER='SLATER'                                                   PSD03030
      READ(5,10,END=1000,ERR=1000) (NAME(I),I=1,18)                     PSD03040
 10   FORMAT(18A4)                                                      PSD03050
      NSUP=0                                                            PSD03060
      NTCT=0                                                            PSD03070
      YPSD=.TRUE.                                                       PSD03080
      WRITE(6,80) (NAME(I),I=1,18)                                      PSD03090
      IF(NAME(1).EQ.IBLA) GO TO 400                                     PSD03100
 80   FORMAT   (1H1,5X,'PSEUDOPOTENTIAL OPTIMIZATION'//5X,18A4)         PSD03110
      READ (5,20) NFLAG1                                                PSD03120
 20   FORMAT(36I2)                                                      PSD03130
      IF(NFLAG1.EQ.1) GOTO 11                                           PSD03140
      WRITE(6,12)                                                       PSD03150
 12   FORMAT(1H0,5X,20HTYPE SLATER ORBITALS)                            PSD03160
      GOTO 14                                                           PSD03170
 11   WRITE(6,13)                                                       PSD03180
 13   FORMAT (1H0,5X,19HTYPE GAUSS ORBITALS)                            PSD03190
   14 READ(5,15) ZN                                                     PSD03200
 15   FORMAT(F5.1)                                                      PSD03210
      WRITE(6,85) ZN                                                    PSD03220
 85    FORMAT(1H0,5X,' CHARGE=',F10.6)                                  PSD03230
      READ(5,45)NSYM                                                    PSD03240
 45   FORMAT(7I1)                                                       PSD03250
      READ(5,20) (NBAS(I),I=1,NSYM)                                     PSD03260
      READ(5,45) (NCSH(I),I=1,NSYM)                                     PSD03270
      READ(5,45) (NOSH(I),I=1,NSYM)                                     PSD03280
      READ(5,45)(NCCUP(I),I=1,NSYM)                                     PSD03290
      DO 600 I=1,4                                                      PSD03300
 600  OCCUP(I)=NCCUP(I)                                                 PSD03310
 25   FORMAT(18F4.1)                                                    PSD03320
      WRITE(6,115) (NBAS(I),I=1,NSYM)                                   PSD03330
 115  FORMAT(1H0,5X,'SYMMETRY SPECIES',11X,'S',5X,'P',5X,'D',5X,'F'/6X, PSD03340
     1'NUMBER OF BASIS FUNCTIONS=',4(I2,4X))                            PSD03350
      WRITE(6,120) (NCSH(I),I=1,NSYM)                                   PSD03360
 120  FORMAT(6X,' NUMBER OF CLOSED SHELLS = ',I1,5X,I1,5X,I1,5X,I1)           PS
      WRITE(6,125) (NCCUP(I),I=1,NSYM)                                  PSD03380
 125  FORMAT(6X,' OPEN SHELL OCCUPATION   = ',I1,5X,I1,5X,I1,5X,I1)           PS
      READ(5,30)  AJMN                                                  PSD03400
 30   FORMAT(5D15.8)                                                    PSD03410
      WRITE(6,135)                                                      PSD03420
 135  FORMAT(1H0,5X,'VECTOR COUPLING COEFFICIENTS K'//)                 PSD03430
      WRITE(6,140)  AJMN                                                PSD03440
 140  FORMAT(5F16.8)                                                    PSD03450
      NBAST=0                                                           PSD03460
      NOSHT=0                                                           PSD03470
      NCSHT=0                                                           PSD03480
      DO 22 I=1,NSYM                                                    PSD03490
      NBAS1=NBAS(I)                                                     PSD03500
      NOSHT=NOSHT+NOSH(I)                                               PSD03510
      NCSHT=NCSHT+NCSH(I)                                               PSD03520
 22   NBAST=NBAST+NBAS1                                                 PSD03530
      NSHT=NCSHT+NOSHT                                                  PSD03540
      READ(5,25) (PQN(J),J=1,NBAST)                                     PSD03550
      READ(5,39) (ZETA(J),J=1,NBAST)                                    PSD03560
      GO TO 450                                                         PSD03570
400   CONTINUE                                                          PSD03580
      DO 405 I=1,4                                                      PSD03590
      NCSH(I)=0                                                         PSD03600
      NBAS(I)=0                                                         PSD03610
      NOSH(I)=0                                                         PSD03620
      NOOC(I)=0                                                         PSD03630
      VCC(I)=0.                                                         PSD03640
      VCC(I+4)=0.                                                       PSD03650
      VCC(I+8)=0.                                                       PSD03660
      NBPOT(I)=0                                                        PSD03670
405   NFLAG(I)=0                                                        PSD03680
      READ(5,ATOM,END=1000)                                             PSD03690
      NSYM=0                                                            PSD03700
      NBAST=0                                                           PSD03710
      NOSHT=0                                                           PSD03720
      NCSHT=0                                                           PSD03730
      DO 410 I=1,4                                                      PSD03740
      IF(NBAS(I).EQ.0) GO TO 420                                        PSD03750
      NSYM=NSYM+1                                                       PSD03760
      NBAST=NBAST+NBAS(I)                                               PSD03770
      NCSHT=NCSHT+NCSH(I)                                               PSD03780
      OCCUP(I)=NOOC(I)                                                  PSD03790
      IF(NCSH(I).EQ.0.AND.NOOC(I).EQ.0) NOSH(I)=1                       PSD03800
      IF(NOOC(I).NE.0) NOSH(I)=1                                        PSD03810
410   NOSHT=NOSHT+NOSH(I)                                               PSD03820
420   NSHT=NCSHT+NOSHT                                                  PSD03830
      ZN=CHARGE                                                         PSD03840
 39   FORMAT(5E14.6)                                                    PSD03850
      IF(NFLAG1.EQ.1) WRITE(6,13)                                       PSD03860
      IF(NFLAG1.EQ.0) WRITE(6,12)                                       PSD03870
      WRITE(6,85) ZN                                                    PSD03880
      WRITE(6,115) (NBAS(I),I=1,NSYM)                                   PSD03890
      WRITE(6,120) (NCSH(I),I=1,NSYM)                                   PSD03900
      WRITE(6,125)(NCCUP(I),I=1,NSYM)                                   PSD03910
      WRITE(6,135)                                                      PSD03920
      WRITE(6,140) (AJMN(I),I=1,15)                                     PSD03930
450   CONTINUE                                                          PSD03940
       WRITE(6,145)                                                     PSD03950
145   FORMAT('0PRINCIPAL QUANTUM NUMBER ORBITALS EXPONENTS AND COEFFICIEPSD03960
     &NTS DEFINING PSEUDO ORBITALS OF S P D F SYMMETRY',/)              PSD03970
      IF(YPSD) GO TO 149                                                PSD03980
      READ(5,30) (PSE(I),I=1,NSHT)                                      PSD03990
C                                                                       PSD04000
      DO 51 J=1,NSHT                                                    PSD04010
   51 READ(5,39) (C(I,J),I=1,NBASM)                                     PSD04020
149   I1=NBAS(1)                                                        PSD04030
      I2=I1+NBAS(2)                                                     PSD04040
      I3=I2+NBAS(3)                                                     PSD04050
      DO 230 I=1,30                                                     PSD04060
      WRITE(6,*)                                                        PSD04070
      YW=.TRUE.                                                         PSD04080
      IF(I1.LT.I) GO TO 160                                             PSD04090
      WRITE(6,150) PQN(I),ZETA(I),C(I,1)                                PSD04100
      YW=.FALSE.                                                        PSD04110
150   FORMAT('+',F3.0,F14.6,F11.6)                                      PSD04120
160   IF(NBAS(2).LT.I) GO TO 180                                        PSD04130
      WRITE(6,170) PQN(I1+I),ZETA(I1+I),C(I,2)                          PSD04140
      YW=.FALSE.                                                        PSD04150
170   FORMAT('+',31X,F3.0,F14.6,F11.6)                                  PSD04160
180   IF(NBAS(3).LT.I) GO TO 200                                        PSD04170
      WRITE(6,190) PQN(I2+I),ZETA(I2+I),C(I,3)                          PSD04180
      YW=.FALSE.                                                        PSD04190
190   FORMAT('+',62X,F3.0,F14.6,F11.6)                                  PSD04200
200   IF(NBAS(4).LT.I) GO TO 220                                        PSD04210
      WRITE(6,210) PQN(I3+I),ZETA(I3+I),C(I,4)                          PSD04220
      GO TO 230                                                         PSD04230
210   FORMAT('+',93X,F3.0,F14.6,F11.6)                                  PSD04240
220   IF(YW) GO TO 240                                                  PSD04250
230   CONTINUE                                                          PSD04260
240   CONTINUE                                                          PSD04270
      WRITE(6,250) (PSE(I),I=1,NSYM)                                    PSD04280
250   FORMAT('0',30X,' PSEUDOORBITAL ENERGIES',/,'0',4(10X,F15.6,5X))   PSD04290
      IF(NSUP.EQ.0) GO TO 500                                           PSD04300
      WRITE(6,320)                                                      PSD04310
  320 FORMAT('0 FONCTIONS DE CONTRAINTE')                               PSD04320
      DO 330 I=1,NBASM                                                  PSD04330
  330 WRITE(6,*)   (FSUP(I,J),J=1,NSUP)                                 PSD04340
C                                                                       PSD04350
500   CONTINUE                                                          PSD04360
C                                                                       PSD04370
      IF(NTCT.EQ.0) GO TO 99                                            PSD04380
      WRITE(6,251)                                                      PSD04390
251   FORMAT('0,    PSEUDOPOTENTIEL CONSTANT POUR TOUTES LES SYMETRIES')PSD04400
      ITY=1                                                             PSD04410
      IF(PSDCST.EQ.SLATER)ITY=0                                         PSD04420
      IF(NFLAG1.EQ.1) ITY=0                                             PSD04430
      DO 630 I=1,NTCT                                                   PSD04440
630   ICT(I)=ITY                                                        PSD04450
      WRITE(6,96)                                                       PSD04460
   96 FORMAT(1H0,5X,3HICT,2X,3HNCT,6X,5HALFCT,12X,4HCACT,/)             PSD04470
      WRITE(6,97)(ICT(I),NCT(I),ALFCT(I),CACT(I),I=1,NTCT)              PSD04480
   97 FORMAT(I7,I6,2D17.8)                                              PSD04490
99    WRITE(6,252)                                                      PSD04500
      WRITE(6,95)                                                       PSD04510
252   FORMAT('0 STARTING PARAMETERS FOR THE OPTIMIZATION IN THE L SYMMETPSD04520
     &RY')                                                              PSD04530
   95 FORMAT(1H0,3H  L,2X,3HISG,2X,3H NX,6X,5H ALFA,12X,4H  CA,/)       PSD04540
      ITY=1                                                             PSD04550
      IF(PSDTYP.EQ.SLATER) ITY=0                                        PSD04560
      IF(NFLAG1.EQ.1) ITY=0                                             PSD04570
      DO 98 L=1,NSYM                                                    PSD04580
      NT(L)=NBPOT(L)                                                    PSD04590
      NTL=NT(L)                                                         PSD04600
      IF(NTL.EQ.0) GO TO 98                                             PSD04610
      DO 610 I=1,NTL                                                    PSD04620
      ISG(I,L)=ITY                                                      PSD04630
      NX(I,L)=NPOT(I)                                                   PSD04640
      ALFA(I,L)=APOT(I)                                                 PSD04650
610   CA(I,L)=CPOT(I)                                                   PSD04660
      WRITE(6,94)(L,ISG(I,L),NX(I,L),ALFA(I,L),CA(I,L),I=1,NTL)         PSD04670
   94 FORMAT(I3,I7,I6,2D17.8)                                           PSD04680
   98 CONTINUE                                                          PSD04690
   21 FORMAT(2I2,2D15.8)                                                PSD04700
       RETURN                                                           PSD04710
 1000 STOP                                                              PSD04720
       END                                                              PSD04730
      SUBROUTINE OEISG                                                  PSD04740
      IMPLICIT REAL*8 (A-H,O-Z)                                         PSD04750
      COMMON/RSVC/FATT(30),FACTO(16),PQN(60),ZETA(60),AJMN(15),         PSD04760
     1OCCUCS(4),OCCUP(4),C(22,8),FC(500),FO(500),S(500),U(500),T(500),  PSD04770
     2P(1627),Q(1627),DT(500),DOS(500),ZN,NBAS(4),NCSH(4),NOSH(4),NSYM, PSD04780
     3NFLAG1,NOSHT                                                      PSD04790
      NSTEP1=0                                                          PSD04800
      NSTEP=0                                                           PSD04810
      DO 15 L=1,NSYM                                                    PSD04820
       NBAS1=NBAS(L)                                                    PSD04830
       DO 10 I=1,NBAS1                                                  PSD04840
       I1=NSTEP1+I                                                      PSD04850
       ZP=ZETA(I1)                                                      PSD04860
       NP=PQN(I1)                                                       PSD04870
       WP=2.*(NP-L)/ZP                                                  PSD04880
       DO 10 J=1,I                                                      PSD04890
       J1=NSTEP1+J                                                      PSD04900
       K=NSTEP+(I*(I-1)/2)+J                                            PSD04910
       ZQ=ZETA(J1)                                                      PSD04920
       ZPQ=0.5*(ZP+ZQ)                                                  PSD04930
       NQ=PQN(J1)                                                       PSD04940
       NPQ=NP+NQ+1                                                      PSD04950
       WQ=2.*(NQ-L)/ZQ                                                  PSD04960
       NPQ1=NPQ-1                                                       PSD04970
       NPQ2=NPQ-2                                                       PSD04980
       IF (NFLAG1.NE.0) GO TO 1                                         PSD04990
       NN=2*NP+1                                                        PSD05000
       NM=2*NQ+1                                                        PSD05010
       VP=FATT(NN)/ZP**NN                                               PSD05020
       VQ=FATT(NM)/ZQ**NM                                               PSD05030
       VPQ=FATT(NPQ)/ZPQ**NPQ                                           PSD05040
       VPQ1=2.*FATT(NPQ1)/ZPQ**NPQ1                                     PSD05050
      VPQ2=FATT(NPQ2)/ZPQ**NPQ2                                         PSD05060
      TERM2=VPQ-(WP+WQ)*VPQ1/2.+WP*WQ*VPQ2                              PSD05070
       GO TO 5                                                          PSD05080
    1  VPQ=FACTO(NPQ1)/ZPQ**(0.5*NPQ)                                   PSD05090
       VP=FACTO(2*NP)/ZP**(NP+0.5)                                      PSD05100
       VQ=FACTO(2*NQ)/ZQ**(NQ+0.5)                                      PSD05110
       IF (NPQ1.GT.2) GO TO 3                                           PSD05120
       VPQM2=1.                                                         PSD05130
        GO TO 4                                                         PSD05140
    3  VPQM2=FACTO(NPQ2-1)/ZPQ**(0.5*NPQ2)                              PSD05150
    4  VPQ1= FACTO(NPQ2)*1.59576912/ZPQ**(0.5*NPQ1)                     PSD05160
       VPQP2=FACTO(NPQ+1)/ZPQ**(0.5*(NPQ+2))                            PSD05170
       TERM2=VPQP2-VPQ*(WP+WQ)+WP*WQ*VPQM2                              PSD05180
    5  TERM1=1./DSQRT(VP*VQ)                                            PSD05190
       S(K)=TERM1*VPQ                                                   PSD05200
       U(K)=TERM1*VPQ1                                                  PSD05210
       T(K)=0.5*ZP*ZQ*TERM1*TERM2                                       PSD05220
   10  CONTINUE                                                         PSD05230
       NSTEP1=NSTEP1+NBAS1                                              PSD05240
   15  NSTEP=NSTEP+NBAS1*(NBAS1+1)/2                                    PSD05250
       RETURN                                                           PSD05260
       END                                                              PSD05270
      SUBROUTINE TEIPSG                                                 PSD05280
      IMPLICIT REAL*8 (A-H,O-Z)                                         PSD05290
      COMMON/RSVC/FATT(30),FACTO(16),PQN(60),ZETA(60),AJMN(15),         PSD05300
     1OCCUCS(4),OCCUP(4),C(22,8),FC(500),FO(500),S(500),U(500),T(500),  PSD05310
     2P(1627),Q(1627),DT(500),DOS(500),ZN,NBAS(4),NCSH(4),NOSH(4),NSYM, PSD05320
     3NFLAG1,NOSHT,NBAST,NFI                                            PSD05330
      J=0                                                               PSD05340
      KMX=0                                                             PSD05350
      DO 100 I=1,NSYM                                                   PSD05360
      KIN=KMX+1                                                         PSD05370
      KMX=KIN+NBAS(I)-1                                                 PSD05380
      DO 100 K=KIN,KMX                                                  PSD05390
      NP=PQN(K)                                                         PSD05400
      ZP=ZETA(K)                                                        PSD05410
      NP2P1=2*NP+1                                                      PSD05420
      DO 100 L=KIN,K                                                    PSD05430
      NQ=PQN(L)                                                         PSD05440
      ZQ=ZETA(L)                                                        PSD05450
      NPQ=NP+NQ                                                         PSD05460
      ZPQ=0.5*(ZP+ZQ)                                                   PSD05470
      NQ2P1=2*NQ+1                                                      PSD05480
      MMX=0                                                             PSD05490
      DO 100 IM=1,I                                                     PSD05500
      MIN=MMX+1                                                         PSD05510
      MMX=MIN-1+NBAS(IM)                                                PSD05520
      IF(IM-I) 20,30,30                                                 PSD05530
 20   MMXP=MMX                                                          PSD05540
      GOTO 31                                                           PSD05550
  30  MMXP=K                                                            PSD05560
 31   DO 100 M=MIN,MMXP                                                 PSD05570
      NR=PQN(M)                                                         PSD05580
      NR2P1=2*NR+1                                                      PSD05590
      ZR=ZETA(M)                                                        PSD05600
      NPR=NP+NR                                                         PSD05610
      NQR=NQ+NR                                                         PSD05620
      ZPR=0.5*(ZP+ZR)                                                   PSD05630
      ZQR=0.5*(ZQ+ZR)                                                   PSD05640
      IF(M-K) 40,41,41                                                  PSD05650
 40   NMX=M                                                             PSD05660
      GOTO 42                                                           PSD05670
   41 NMX=L                                                             PSD05680
 42   DO 100 N=MIN,NMX                                                  PSD05690
      IF(J.LT.1627) GO TO 43                                            PSD05700
      WRITE(NFI) P,Q                                                    PSD05710
      J=0                                                               PSD05730
   43 J=J+1                                                             PSD05740
      NS=PQN(N)                                                         PSD05750
      ZS=ZETA(N)                                                        PSD05760
       NRS=NR+NS                                                        PSD05770
      ZRS=0.5*(ZR+ZS)                                                   PSD05780
      NS2P1=2*NS+1                                                      PSD05790
      NQS=NQ+NS                                                         PSD05800
      NPS=NP+NS                                                         PSD05810
      ZQS=0.5*(ZQ+ZS)                                                   PSD05820
      ZPS=0.5*(ZP+ZS)                                                   PSD05830
      IF(NFLAG1.EQ.1) GOTO 6                                            PSD05840

C      VP=FATT(NP2P1)/ZP**NP2P1
C      VP=VP*FATT(NQ2P1)/ZQ**NQ2P1
C      VP=VP*FATT(NR2P1)/ZR**NR2P1
C      VP=VP*FATT(NS2P1)/ZS**NS2P1
      vp=dsqrt(fatt(np2p1))
      vp=vp*dsqrt(fatt(nq2p1))
      vp=vp*dsqrt(fatt(nr2p1))
      vp=vp*dsqrt(fatt(ns2p1))
      deno=dsqrt(zp**np2p1)
      deno=deno*dsqrt(zq**nq2p1)
      deno=deno*dsqrt(zr**nr2p1)
      deno=deno*dsqrt(zs**ns2p1)
      SSS=1.                                                            PSD05890
      XJ=2.*(SINY(0,NPQ,NRS,ZPQ,ZRS)+SINY(0,NRS,NPQ,ZRS,ZPQ))           PSD05900
      GOTO 7                                                            PSD05910
 6    VP=FACTO(NP2P1-1)/ZP**(NP+0.5)                                    PSD05920
      VP=VP*FACTO(NQ2P1-1)/ZQ**(NQ+0.5)                                 PSD05930
      VP=VP*FACTO(NR2P1-1)/ZR**(NR+0.5)                                 PSD05940
      VP=VP*FACTO(NS2P1-1)/ZS**(NS+0.5)                                 PSD05950
      SSS=1.59576912                                                    PSD05960
      XJ=GINY(0,NPQ,NRS,ZPQ,ZRS)+GINY(0,NRS,NPQ,ZRS,ZPQ)                PSD05970
 7    IF(NFLAG1.EQ.1) THEN
      TERM2=SSS/DSQRT(VP) 
      ELSE
      TERM2=SSS*DENO/VP
      ENDIF  
      XJ=XJ*TERM2                                                       PSD05990
      SSS=0.                                                            PSD06000
      SSU=0.                                                            PSD06010
      NULO=I-IM+1                                                       PSD06020
      NUHI=I+IM-1                                                       PSD06030
      DO 50 KIK=NULO,NUHI,2                                             PSD06040
      NU=KIK-NULO                                                       PSD06050
       BUM1=FATT(NU+1)/FATT((NU+2)/2)**2                                PSD06060
      NU=KIK+NULO                                                       PSD06070
      BUM1=BUM1*FATT(NU-1)/FATT(NU/2)**2                                PSD06080
      NU=NUHI-KIK                                                       PSD06090
      BUM1=BUM1*FATT(NU+1)/FATT((NU+2)/2)**2                            PSD06100
      NU=NUHI+KIK                                                       PSD06110
      BUM1=BUM1/(FATT(NU-1)/FATT(NU/2)**2)                              PSD06120
      ALMN=BUM1/(NU-1)                                                  PSD06130
      NU=KIK-1                                                          PSD06140
      IF(NFLAG1.EQ.1) GOTO 21                                           PSD06150
      YJNU=SINY(NU,NPR,NQS,ZPR,ZQS)+SINY(NU,NQS,NPR,ZQS,ZPR)+SINY(NU,NPSPSD06160
     1,NQR,ZPS,ZQR)+SINY(NU,NQR,NPS,ZQR,ZPS)                            PSD06170
      YJNU=YJNU*TERM2                                                   PSD06180
      GOTO 22                                                           PSD06190
 21   YJNU=0.5*(GINY(NU,NPR,NQS,ZPR,ZQS)+GINY(NU,NQS,NPR,ZQS,ZPR)+GINY( PSD06200
     1NU,NPS,NQR,ZPS,ZQR)+GINY(NU,NQR,NPS,ZQR,ZPS))*TERM2               PSD06210
   22 SSU=SSU+YJNU*ALMN                                                 PSD06220
      IF(NOSHT) 50,50,52                                                PSD06230
 52   IF(I-IM) 50,59,53                                                 PSD06240
 53   IF(NU-2) 56,55,54                                                 PSD06250
  54  AKOO=AJMN(7)                                                      PSD06260
      GOTO 68                                                           PSD06270
 55   AKOO=AJMN(5)                                                      PSD06280
      GOTO 68                                                           PSD06290
 56   IF(I-2) 50,57,58                                                  PSD06300
 57   AKOO=AJMN(2)                                                      PSD06310
      GOTO 68                                                           PSD06320
 58   AKOO=AJMN(6)                                                      PSD06330
      GOTO 68                                                           PSD06340
 59   IF(I.GE.4)GO TO 80                                                PSD06350
      IF(NU-2) 64,61,60                                                 PSD06360
 60   AKOO=AJMN(10)                                                     PSD06370
      GOTO 68                                                           PSD06380
 61   IF(I-2) 50,63,62                                                  PSD06390
 62   AKOO=AJMN(9)                                                      PSD06400
      GOTO 68                                                           PSD06410
 63   AKOO=AJMN(4)                                                      PSD06420
      GOTO 68                                                           PSD06430
 64   IF(I-2) 65,66,67                                                  PSD06440
 65   AKOO=AJMN(1)                                                      PSD06450
      GOTO 68                                                           PSD06460
 66   AKOO=AJMN(3)                                                      PSD06470
      GOTO 68                                                           PSD06480
 67   AKOO=AJMN(8)                                                      PSD06490
      GO TO 68                                                          PSD06500
 80   IFNU= NU/2 + 1                                                    PSD06510
      GOTO(81,82,83,84),IFNU                                            PSD06520
 81   AKOO=AJMN(12)                                                     PSD06530
      GO TO 68                                                          PSD06540
 82   AKOO=AJMN(13)                                                     PSD06550
      GO TO 68                                                          PSD06560
 83   AKOO=AJMN(14)                                                     PSD06570
      GO TO 68                                                          PSD06580
 84   AKOO=AJMN(15)                                                     PSD06590
 68   SSS=SSS+AKOO*YJNU                                                 PSD06600
 50   CONTINUE                                                          PSD06610
       P(J)=XJ-0.5*SSU                                                  PSD06620
      Q(J)=-0.5*SSS                                                     PSD06630
      IF(IM.NE.2) GOTO 100                                              PSD06640
      IF(I.NE.3) GOTO 100                                               PSD06650
      IF(AJMN(11).EQ.0.) GOTO 100                                       PSD06660
      IF(NFLAG1.EQ.1) GOTO 73                                           PSD06670
      XJNU=2.*(SINY(2,NPQ,NRS,ZPQ,ZRS)+SINY(2,NRS,NPQ,ZRS,ZPQ))         PSD06680
      GOTO 74                                                           PSD06690
 73   XJNU=GINY(2,NPQ,NRS,ZPQ,ZRS)+GINY(2,NRS,NPQ,ZRS,ZPQ)              PSD06700
 74   Q(J)=Q(J)+AJMN(11)*TERM2*XJNU                                     PSD06710
 100  CONTINUE                                                          PSD06720
C                                                                       PSD06730
      WRITE(NFI) P,Q                                                    PSD06740
      RETURN                                                            PSD06750
      END                                                               PSD06760
      FUNCTION SINY(NY,NAB,NCD,ZAB,ZCD)                                 PSD06770
      IMPLICIT REAL*8 (A-H,O-Z)                                         PSD06780
      COMMON/RSVC/FATT(30),FACTO(16)                                    PSD06790
      TABCD=ZAB/ZCD                                                     PSD06800
      NABNY=NAB-NY-1                                                    PSD06810
      NCDNY=NCD+NY                                                      PSD06820
      VNAB=FATT(NABNY+1)/ZAB**(NABNY+1)                                 PSD06830
      VNCD=FATT(NCDNY+1)/ZCD**(NCDNY+1)                                 PSD06840
      TERM1=0.                                                          PSD06850
      NALF=NABNY+1                                                      PSD06860
      DO1I=1,NALF                                                       PSD06870
      TERM1=TERM1+FATT(NABNY+NCDNY+2)*TABCD**(I-1)/(FATT(I)*FATT(NABNY+NPSD06880
     1CDNY+3-I))                                                        PSD06890
  1   CONTINUE                                                          PSD06900
      CABCD=TERM1*(TABCD+1.)**(-NABNY-NCDNY-1)                          PSD06910
      SINY=VNAB*VNCD*CABCD                                              PSD06920
      RETURN                                                            PSD06930
      END                                                               PSD06940
      FUNCTION GINY(NY,NAB,NCD,ZAB,ZCD)                                 PSD06950
      IMPLICIT REAL*8 (A-H,O-Z)                                         PSD06960
      COMMON/RSVC/FATT(30),FACTO(16)                                    PSD06970
      TABCD=ZAB/ZCD                                                     PSD06980
      NABNY=NAB-NY-1                                                    PSD06990
      NCDNY=NCD+NY                                                      PSD07000
      VNAB=FACTO(NABNY)/ZAB**(0.5*(NABNY+1))                            PSD07010
      VNCD=FACTO(NCDNY)/ZCD**(0.5*(NCDNY+1))                            PSD07020
      TERM1=0.                                                          PSD07030
      DO1I=1,NABNY,2                                                    PSD07040
      TERM1=TERM1+FACTO(I+NCDNY-1)/(FACTO(I)*FACTO(NCDNY))*(TABCD/(1.+TAPSD07050
     1BCD))**(0.5*(I-1))                                                PSD07060
  1   CONTINUE                                                          PSD07070
      CABCD=TERM1*(1.+TABCD)**(-0.5*(NCDNY+1))                          PSD07080
      GINY=VNAB*VNCD*CABCD                                              PSD07090
      RETURN                                                            PSD07100
      END                                                               PSD07110
      SUBROUTINE DENSI                                                  PSD07120
      IMPLICIT REAL*8 (A-H,O-Z)                                         PSD07130
      COMMON/RSVC/FATT(30),FACTO(16),PQN(60),ZETA(60),AJMN(15),         PSD07140
     1OCCUCS(4),OCCUP(4),C(22,8),FC(500),FO(500),S(500),U(500),T(500),  PSD07150
     2P(1627),Q(1627),DT(500),DOS(500),ZN,NBAS(4),NCSH(4),NOSH(4),NSYM, PSD07160
     3NFLAG1,NOSHT                                                      PSD07170
      OCCUCS(1)=2.                                                      PSD07180
      OCCUCS(2)=6.                                                      PSD07190
      OCCUCS(3)=10.                                                     PSD07200
      OCCUCS(4)=14.                                                     PSD07210
      K=1                                                               PSD07220
      NSTEP=0                                                           PSD07230
      J1=0                                                              PSD07240
      DO45I=1,NSYM                                                      PSD07250
      IS=0                                                              PSD07260
      IK=0                                                              PSD07270
      NBAS1=NBAS(I)                                                     PSD07280
      NOSHIC=NOSH(I)                                                    PSD07290
      NCSHIC=NCSH(I)                                                    PSD07300
      NOCSH=NOSHIC+NCSHIC                                               PSD07310
      J1=J1+NOCSH                                                       PSD07320
      IF(NOSHIC.GT.0)IS=1                                               PSD07330
      IF(NCSHIC.GT.0)IK=1                                               PSD07340
      DO75M=1,NBAS1                                                     PSD07350
      DO75N=1,M                                                         PSD07360
      TERM1=OCCUP(I)*C(M,J1)*C(N,J1)                                    PSD07370
      TERM2=0.                                                          PSD07380
      DO5J=1,NCSHIC                                                     PSD07390
      J2=J+NSTEP                                                        PSD07400
  5   TERM2=TERM2+C(M,J2)*C(N,J2)                                       PSD07410
      TERM2=TERM2*OCCUCS(I)                                             PSD07420
      IF(M.EQ.N)GOTO35                                                  PSD07430
      TERM1=2.*TERM1                                                    PSD07440
      TERM2=2.*TERM2                                                    PSD07450
  35  DOSA=TERM1*IS                                                     PSD07460
      DT(K)=TERM2*IK+DOSA                                               PSD07470
      DOS(K)=DOSA                                                       PSD07480
  75  K=K+1                                                             PSD07490
  45  NSTEP=J1                                                          PSD07500
      RETURN                                                            PSD07510
      END                                                               PSD07520
      SUBROUTINE HAMIL                                                  PSD07530
      IMPLICIT REAL*8 (A-H,O-Z)                                         PSD07540
      DIMENSION N1(6),PCAP(500),QCAP(500),SMIN(22,8),GMIN(22,8)         PSD07550
      COMMON/RSVC/FATT(30),FACTO(16),PQN(60),ZETA(60),AJMN(15),         PSD07560
     1OCCUCS(4),OCCUP(4),C(22,8),FC(500),FO(500),S(500),U(500),T(500),  PSD07570
     2P(1627),Q(1627),DT(500),DOS(500),ZN,NBAS(4),NCSH(4),NOSH(4),NSYM, PSD07580
     3NFLAG1,NOSHT,NBAST,NFI                                            PSD07590
      N1T=0                                                             PSD07600
      ITEST=0                                                           PSD07610
      DO5I=1,NSYM                                                       PSD07620
      NBAS1=NBAS(I)                                                     PSD07630
      ITEST=ITEST+NCSH(I)*2+OCCUP(I)+.1
      WRITE(6,*) ' ITEST ',ITEST
      WRITE(6,*) NCSH
      WRITE(6,*) OCCUP
      NSTEP1=NBAS1*(NBAS1+1)/2                                          PSD07650
      N1T=N1T+NSTEP1                                                    PSD07660
  5   N1(I)=NSTEP1                                                      PSD07670
C                                                                       PSD07680
      IPQ=1627                                                          PSD07690
      DO 3 I=1,N1T                                                      PSD07700
      PCAP(I)=0.D0                                                      PSD07710
    3 QCAP(I)=0.D0                                                      PSD07720
      IF(ITEST.LE.1) GO TO 40                                           PSD07730
      REWIND NFI                                                        PSD07740
      DO 35 I=1,N1T                                                     PSD07750
      DO 30 J=1,I                                                       PSD07760
      IPQ=IPQ+1                                                         PSD07770
      IF(IPQ.LE.1627) GO TO 36                                          PSD07780
      READ (NFI) P,Q                                                    PSD07790
C                                                                       PSD07800
      IPQ=1                                                             PSD07810
   36 PCAP(I)=PCAP(I)+P(IPQ)*DT(J)                                      PSD07820
      QCAP(I)=QCAP(I)+Q(IPQ)*DOS(J)                                     PSD07830
      IF(I.EQ.J) GO TO 30                                               PSD07840
      PCAP(J)=PCAP(J)+P(IPQ)*DT(I)                                      PSD07850
      QCAP(J)=QCAP(J)+Q(IPQ)*DOS(I)                                     PSD07860
   30 CONTINUE                                                          PSD07870
   35 CONTINUE                                                          PSD07880
40    CONTINUE                                                          PSD07890
      IF(ITEST.LE.1) WRITE(6,601)                                       PSD07900
601   FORMAT(' ATTENTION PROBLEME A 0 OU 1 ELECTRON PAS DE PQRS')       PSD07910
      NSTEP1=0                                                          PSD07920
      NSTEP2=0                                                          PSD07930
      DO70I=1,NSYM                                                      PSD07940
      NSH=NCSH(I)+NOSH(I)                                               PSD07950
      NBAS1=NBAS(I)                                                     PSD07960
      DO65J=1,NSH                                                       PSD07970
      J1=NSTEP2+J                                                       PSD07980
      DO65M=1,NBAS1                                                     PSD07990
      FACT1=0.                                                          PSD08000
      FACT2=0.                                                          PSD08010
      DO2N=1,NBAS1                                                      PSD08020
      IF(N-M)50,50,55                                                   PSD08030
  50  K=(M*(M-1)/2)+N+NSTEP1                                            PSD08040
      GOTO60                                                            PSD08050
  55  K=(N*(N-1)/2)+M+NSTEP1                                            PSD08060
  60  FACT1=FACT1+S(K)*C(N,J1)                                          PSD08070
  2   FACT2=FACT2+QCAP(K)*C(N,J1)                                       PSD08080
      SMIN(M,J1)=FACT1                                                  PSD08090
  65  GMIN(M,J1)=FACT2                                                  PSD08100
      NSTEP2=NSTEP2+NSH                                                 PSD08110
  70  NSTEP1=NSTEP1+N1(I)                                               PSD08120
      K=1                                                               PSD08130
      NSTEP=0                                                           PSD08140
      DO105I=1,NSYM                                                     PSD08150
      NOSH1=NOSH(I)                                                     PSD08160
      NSH=NCSH(I)+NOSH1                                                 PSD08170
      FACT2=OCCUP(I)                                                    PSD08180
      FACT1=OCCUCS(I)-FACT2                                             PSD08190
      NBAS1=NBAS(I)                                                     PSD08200
      DO100M=1,NBAS1                                                    PSD08210
      DO100N=1,M                                                        PSD08220
      TERM1=0.                                                          PSD08230
      TERM2=0.                                                          PSD08240
      DO95J=1,NSH                                                       PSD08250
      J1=J+NSTEP                                                        PSD08260
      IF(J-NSH)90,80,90                                                 PSD08270
  80  IF(NOSH1)85,90,85                                                 PSD08280
  85  TERM2=SMIN(M,J1)*GMIN(N,J1)+GMIN(M,J1)*SMIN(N,J1)                 PSD08290
      GOTO95                                                            PSD08300
  90  TERM1=TERM1+SMIN(M,J1)*GMIN(N,J1)+GMIN(M,J1)*SMIN(N,J1)           PSD08310
  95  CONTINUE                                                          PSD08320
      HPPC=PCAP(K)+T(K)-ZN*U(K)                                         PSD08330
      FC(K)=HPPC+TERM2*FACT2/FACT1                                      PSD08340
      FO(K)=HPPC-QCAP(K)+TERM1*OCCUCS(I)/FACT1                          PSD08350
  100 K=K+1                                                             PSD08360
  105 NSTEP=NSTEP+NSH                                                   PSD08370
      RETURN                                                            PSD08380
      END                                                               PSD08390
      SUBROUTINE HPSEUD(TAF,APQ,CIK,VIK,WIJ,AIH,NTL,NBAS1,NEND)         PSD08400
       IMPLICIT REAL*8 (A-H,O-Z)                                        PSD08410
      DIMENSION HPS(22,22),HPQ(22,22),VP(22),HD(22,22),VECT(22,22)      PSD08420
      COMMON/RSVC/FATT(30),FACTO(16),PQN(60),ZETA(60),AJMN(15),         PSD08430
     1OCCUCS(4),OCCUP(4),C(22,8),FC(500),FO(500),S(500),U(500),T(500),  PSD08440
     2P(1627),Q(1627),DT(500),DOS(500),ZN,NBAS(4),NCSH(4),NOSH(4),NSYM, PSD08450
     3NFLAG1,NOSHT                                                      PSD08460
      COMMON/SUP/FSUP(22,10),ESUP(10),VHSUPV(10),R000,NSUP                   PSD
C                                                                       PSD08480
      DIMENSION TAF(22,22),APQ(22,22),CIK(22),AIH(10)                   PSD08490
      DIMENSION WIJ(NBAS1,NBAS1,NTL)                                    PSD08500
      DIMENSION TPM(22)                                                 PSD08510
      DO 75 I=1,NBAS1                                                   PSD08520
      DO 75 J=1,I                                                       PSD08530
      SOM=0.D0                                                          PSD08540
      DO 76 IT=1,NTL                                                    PSD08550
   76 SOM=SOM+AIH(IT)*WIJ(I,J,IT)                                       PSD08560
      HPS(I,J)=TAF(I,J)+SOM                                             PSD08570
   75 HPS(J,I)=HPS(I,J)                                                 PSD08580
      DO 23 I=1,NBAS1                                                   PSD08590
      DO 60 IP=1,NBAS1                                                  PSD08600
      SOM=0.                                                            PSD08610
      DO 50 IQ=1,NBAS1                                                  PSD08620
50    SOM=SOM+APQ(I,IQ)*HPS(IQ,IP)                                      PSD08630
60    TPM(IP)=SOM                                                       PSD08640
      DO 23 J=1,I                                                       PSD08650
      SOM=0.                                                            PSD08660
      DO 24 IP=1,NBAS1                                                  PSD08670
C                                                                       PSD08680
24    SOM=SOM+TPM(IP)*APQ(IP,J)                                         PSD08690
      HPQ(I,J)=SOM                                                      PSD08700
23    HPQ(J,I)=SOM                                                      PSD08710
      CALL DVLVC1(HPQ,22,22,VP,22,HD,22,22,NBAS1,IIT)                   PSD08720
       IF(NEND.EQ.0) GOTO 40                                            PSD08730
      WRITE(6,631) (VP(I),I=1,NBAS1)                                    PSD08740
      REWIND 8                                                          PSD08750
      WRITE(8) HPS                                                      PSD08760
      KKM=0                                                             PSD08770
      DO 85 I=1,NBAS1                                                   PSD08780
      DO 85 J=1,I                                                       PSD08790
      KKM=KKM+1                                                         PSD08800
      HPS(I,J)=HPS(I,J)-TAF(I,J)+FC(KKM)                                PSD08810
85    HPS(J,I)=HPS(I,J)                                                 PSD08820
      WRITE(8) HPS                                                      PSD08830
  631 FORMAT(/' VAL.PROP. ',4D15.8)                                     PSD08840
  40  DO 26 I=1,NBAS1                                                   PSD08850
      DO 26 J=1,NBAS1                                                   PSD08860
      SC=0.D0                                                           PSD08870
      DO 27 IR=1,NBAS1                                                  PSD08880
   27 SC=SC+HD(I,IR)*APQ(IR,J)                                          PSD08890
      VECT(I,J)=SC                                                      PSD08900
   26 CONTINUE                                                          PSD08910
      VIK=VP(1)                                                         PSD08920
      DO 30 I=1,NBAS1                                                   PSD08930
   30 CIK(I)=VECT(1,I)                                                  PSD08940
      IF(NEND.EQ.0) RETURN                                              PSD08950
      SOM=0.                                                            PSD08960
      DO 90 I=1,NBAS1                                                   PSD08970
      DO 90 J=1,NBAS1                                                   PSD08980
90    SOM=SOM+CIK(I)*HPS(I,J)*CIK(J)                                    PSD08990
      WRITE(6,91) SOM                                                   PSD09000
91    FORMAT('0 VALEUR MOYENNE MONOELECTRONIQUE',D15.8)                 PSD09010
      X=0.                                                              PSD09020
      DO 510 I=1,NBAS1                                                  PSD09030
      DO 510 J=1,NBAS1                                                  PSD09040
510   X=X+CIK(I)*TAF(I,J)*CIK(J)                                        PSD09050
      WRITE(6,520) X                                                    PSD09060
520   FORMAT(' VALEUR MOYENNE FONCTION SUR HPS MOYEN',F20.8)            PSD09070
      WRITE(6,530) (CIK(I),I=1,NBAS1)                                   PSD09080
530   FORMAT(' VECTEUR PROPRE',10F10.6,/,(1X,10F10.6))                  PSD09090
      RETURN                                                            PSD09100
      END                                                               PSD09110
      FUNCTION WPS(Z1,A1,N1,Z2,A2,N2,IX,KX,AX)                          PSD09120
       IMPLICIT REAL*8 (A-H,O-Z)                                        PSD09130
C                                                                       PSD09140
C                                                                       PSD09150
C                                                                       PSD09160
C                                                                       PSD09170
C                                                                       PSD09180
      COEF=A1*A2                                                        PSD09190
      IF(IX.EQ.1) GO TO 1                                               PSD09200
      WPS=COEF*PURX(N1,N2,Z1,Z2,KX,AX)                                  PSD09210
      RETURN                                                            PSD09220
    1 NIJ=N1+N2+KX                                                      PSD09230
      RAX=DSQRT(AX)                                                     PSD09240
      NIJ1=NIJ+1                                                        PSD09250
      ZIJ=(Z1+Z2)/RAX                                                   PSD09260
      PALF=1.D0/(RAX**NIJ1)                                             PSD09270
      WPS=COEF*PALF*DSLAG(NIJ,ZIJ)                                      PSD09280
      ENTRY UDRFLW                                                      PSD09290
      RETURN                                                            PSD09300
      END                                                               PSD09310
      DOUBLE PRECISION FUNCTION PURX(N1,N2,Z1,Z2,IX,ALF)                PSD09320
      IMPLICIT REAL*8 (A-H,O-Z)                                         PSD09330
      COMMON/RSVC/FATT(30),FACTO(16),PQN(60),ZETA(60),AJMN(15),         PSD09340
     &OCCUCS(4),OCCUP(4),CC(22,8),FC(2500),P(1627),Q(1627),DT(1000),    PSD09350
     &ZN,NBAS(4),NCSH(4),NOSH(4),NSYM,NFLAG1                            PSD09360
      DOUBLE PRECISION C,Z1,Z2,DEP,ALF                                  PSD09370
      DATA PI/3.14159265358979/                                         PSD09380
      PURX=0.D0                                                         PSD09390
      M=N1+N2+IX                                                        PSD09400
      C=1.D0/(Z1+Z2+ALF)                                                PSD09410
      IF(NFLAG1.NE.0) GO TO 20                                          PSD09420
      IF(M) 99,11,12                                                    PSD09430
   11 PURX=C                                                            PSD09440
      GO TO 99                                                          PSD09450
12    PURX=C**(M+1)*FATT(M+1)                                           PSD09460
      GO TO 99                                                          PSD09470
20    M2=M/2                                                            PSD09480
      IF((M2*2).NE.M) GO TO 60                                          PSD09490
      IF(M.EQ.0) GO TO 50                                               PSD09500
      PURX=DSQRT(PI*C)*FACTO(M)*.5*(.5*C)**M2                           PSD09510
      GO TO 99                                                          PSD09520
50    PURX=.5D0*DSQRT(PI*C)                                             PSD09530
      GO TO 99                                                          PSD09540
60    M2=M2+1                                                           PSD09550
      PURX=FATT(M2)*.5D0*C**M2                                          PSD09560
   99 RETURN                                                            PSD09570
      END                                                               PSD09580
      SUBROUTINE CALCFX(IN,XU,GAM)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/RSVC/FATT(30),FACTO(16),PQN(60),ZETA(60),AJMN(15),
     &OCCUCS(4),OCCUP(4),CC(22,8),FC(2500),P(1627),Q(1627),DT(1000),
     &ZN,NBAS(4),NCSH(4),NOSH(4),NSYM,NFLAG1
      DIMENSION CIK(22),XU(6)
      COMMON/SUP/FSUP(22,10),ESUP(10),VHSUPV(10),R000,NSUP
      DIMENSION WIJ(9680),AIAJ(22,22),AIH(22)
      COMMON/PSEUD/ALFCT(30),CACT(30),ALFA(30,4),CA(30,4),PSE(4),
     1ICT(30),ISG(30,4),NT(4),NTCT,NCT(30),NX(30,4)
      COMMON/UN/TAF(22,22,3),APQ(22,22,3),ZO(22,3),AO(22,3),PSEO(3),
     1CPSO(22,3),NO(22,3),HV(22,3),VHV(3),L,NBAS1,NETAT,NEND
      DIMENSION SSOP(22,22),SSVP(22),SSV(22,22),SSD(22)
      DIMENSION VAV(22,3),AV(22,22,3)
      LOGICAL*1 YCAL
      NTL=NT(L)
      YCAL=IN.GT.NTL
      IF(YCAL) GO TO 110
      DO 2 I=1,NTL
      AIH(I)=0.D0
      DO 2 J=1,NTL
    2 AIAJ(I,J)=0.D0
      DO 3 IK=1,NETAT
      IWIJ=NBAS1*NBAS1*NTL*(IK-1)+1
      CALL CAIH(APQ(1,1,IK),ZO(1,IK),AO(1,IK),NO(1,IK),XU(1),
     1CPSO(1,IK),PSEO(IK),HV(1,IK),WIJ(IWIJ),AIAJ(1,1),AIH(1),L,NTL,
     2NBAS1,VAV(1,IK),AV(1,1,IK),VHV(IK))
    3 CONTINUE
      IF(IN.LE.1) GOTO2101
      CALL COUTR(AIAJ,AIH,ALFA(1,L),NX(1,L),NTL,VAV,VHV,AV(1,1,1))
      GO TO 130
2101  CONTINUE
      DO 19 I=1,NTL
19    SSD(I)=DSQRT(AIAJ(I,I))
      DO 20 I=1,NTL
      AIH(I)=AIH(I)/SSD(I)
      DO 20 J=1,I
      SSOP(I,J)=AIAJ(I,J)/(SSD(I)*SSD(J))
20    SSOP(J,I)=SSOP(I,J)
      IF(NTL.EQ.1) GO TO 210
      CALL DVLVC1(SSOP,22,22,SSVP,22,SSV,22,22,NTL,IIT)
      GO TO 220
210    CONTINUE
      SSVP(1)=SSOP(1,1)
      SSV(1,1)=1.
220   CONTINUE
30    FORMAT(' VALEUR PROPRE DE LA MATRICE DES OPERATEURS:',(1X,6E12.2))
      DO 40 I=1,NTL
      SSVP(I)=1.D0/SSVP(I)
40    IF(SSVP(I).GT.1.D10) SSVP(I)=0.D0
      DO 60 I=1,NTL
      DO 60 J=1,I
      TOT=0.
      DO 50 IK=1,NTL
50    TOT=TOT+SSV(IK,I)*SSV(IK,J)*SSVP(IK)
      AIAJ(I,J)=TOT
60    AIAJ(J,I)=TOT
      DO 80 I=1,NTL
      TOT=0.
      DO 70 J=1,NTL
70    TOT=TOT+AIAJ(I,J)*AIH(J)
80    SSD(I)=TOT/SSD(I)
      DO 90 I=1,NTL
90    AIH(I)=SSD(I)
      GO TO 130
110     DO 120 I=1,NTL
120    AIH(I)=XU(I)
130    CONTINUE
      PGAM=1.D0
      GAM=0.D0
      DO 1 IK=1,NETAT
      IWIJ=NBAS1*NBAS1*NTL*(IK-1)+1
      IF(YCAL) CALL CAWIJ(WIJ(IWIJ),ZO(1,IK),NO(1,IK),XU(NTL+1),L,NTL,NB
     &AS1,AO(1,IK))
      CALL HPSEUD(TAF(1,1,IK),APQ(1,1,IK),CIK(1),VIK,WIJ(IWIJ),
     1AIH(1),NTL,NBAS1,NEND)
      REC=0.D0
      REC1=0.D0
      DO 10 J1=1,NBAS1
      N1=NO(J1,IK)
      Z1=ZO(J1,IK)
      DO 10 J2=1,NBAS1
      N2=NO(J2,IK)
      Z2=ZO(J2,IK)
      COEF=AO(J1,IK)*CIK(J1)*AO(J2,IK)*CPSO(J2,IK)
      TOT=AO(J1,IK)*AO(J2,IK)*CPSO(J2,IK)*PURX(N1,N2,Z1,Z2,0,0.D0)
      REC=REC+CIK(J1)*TOT
      REC1=TOT*CPSO(J1,IK)+REC1
   10 CONTINUE
      GNOR=DABS(VIK*REC*REC-PSEO(IK)*REC1*REC1)+DABS(VIK*REC1*REC1-REC*R
     &EC*PSEO(IK))
      IF(NEND.EQ.0) GO TO 41
      WRITE(8) NBAS1,ZO,NO
      WRITE(6,652) REC,REC1
652    FORMAT(' RECOUVREMENT PHI PHIP',D15.8,' NORME PHI',D15.8)
      REC=REC-REC1
      WRITE(6,634) REC,GNOR,(XU(I),I=1,IN)
  634 FORMAT(' REC=',D18.12,'  NORM=',D18.12,/,3X,6D18.8)
      WRITE(6,30) (SSVP(I),I=1,NTL)
      WRITE(6,651)(CIK(I),I=1,NBAS1)
651    FORMAT('0 FONCTION PROPRE',(1X,10F11.7))
      IF(GNOR.GT.1.D-2) GO TO 770
      WRITE(6,730)
730   FORMAT(10X,'R',10X,'PHI',7X,'PHI PRIME',7X,
     &'DERIVEES DE PHI  PRIME',
     &12X,'PSEUDOPOTENTIEL',3X,'L(L+1)/2R2',7X,'1/R')
      RA=1.D-6
       RAIN=.1
      DO 750 I=1,500
      T=0.
      TOTO=0.
      TD=0.
      TOT=0.
      TOTG=0.
      TDG=0.
      RU=RA
      IF(NFLAG1.NE.0) RU=RA*RA
      DO 740 J=1,NBAS1
      EXPON=DMIN1(173.D0,RU*ZO(J,IK))
      TA=AO(J,IK)*RA**NO(J,IK)*DEXP(-EXPON)
      T=T+TA*CPSO(J,IK)
      TA=TA*CIK(J)
      TOT=TOT+TA
      TOTO=TOTO+TA*(NO(J,IK)/ RA-ZO(J,IK))
      TOTS= NO(J,IK)/RA-(RA+RA)*ZO(J,IK)
      TOTG=TOTG+TA*TOTS
      TDG=TDG+TA*(TOTS*TOTS-NO(J,IK)/RU-2.D0*ZO(J,IK))
      TD=TD+TA*(NO(J,IK)*(NO(J,IK)-1)/(RA*RA)-2.D0*ZO(J,IK)*NO(J,IK)
     &/RA+ZO(J,IK)**2)
740    CONTINUE
      TUT=0.D0
      DO 741 J=1,NTL
      INDAL=MIN0(J,IN)
      EXPON=DMIN1(173.D0,RA*RA*XU(INDAL))
741   TUT=TUT+AIH(J)*RA** NX(J,L)*DEXP(-EXPON)
      IF(NFLAG1.EQ.0) GO TO 747
      TOTO=TOTG
      TD=TDG
747   CONTINUE
      TAT=1.D0/RA
       TATA=L*(L-1)/(RA*RA)*0.5D0
      WRITE(6,760) RA,T,TOT,TOTO,TD,TUT,TATA,TAT
      IF(DABS(T).LT.1.D-2.AND.RA.GT.5.D0) GO TO 770
      IF(RA.GT.2.D0) RAIN=.1D0*RA
750   RA=RA+RAIN
770    CONTINUE
760   FORMAT(1X,9F14.6)
      IF(YCAL) GO TO 41
      DO 150 I=1,NTL
150    XU(I+1)=AIH(I)
41      GAM=GAM+GNOR
    1 CONTINUE
      PGAM=PGAM*GNOR
      WRITE(6,640) (AIH(I),I=1,NTL),GNOR,XU(1)
  640 FORMAT(5X,8D15.8)
      IF(NETAT.EQ.1) RETURN
      GAM=GAM/DFLOAT(NETAT)
      GAM=GAM*GAM*GAM/PGAM
      RETURN
      ENTRY ICALCF
      TOT=FACTOR(I)
      R00=R000
      IF(R000.GE.0.D0) RETURN
      R00=50.
      DEL=1.
      SIG=1.
      DO 810 I=1,500
      TOTO=TOT
      TOT=0.
      IF(R00.LT.0.D0) GO TO 840
      DO 800 J=1,NBAS1
      R0A=R00*ZO(J,1)
      IF(R0A.GT.173.D0)GO TO 800
      TOT=TOT+AO(J,1)*CPSO(J,1)*R00**NO(J,1)*DEXP(-R0A)
800    CONTINUE
      TOT=DABS(TOT)
      IF(R00.GT.51.D0) GO TO 850
      IF(I.EQ.1) GO TO 810
         IF(TOT.GT.TOTO) GO TO 810
      IF(TOT.LT.0.05) GO TO 810
      IF(DEL.LT.0.001) GO TO 820
      GO TO 860
850   R00=50.
      GO TO 860
840    R00=0.
860   CONTINUE
      SIG=-SIG
      DEL=DEL/10
810   R00=R00-DEL*SIG
820   WRITE(6,830) R00
830   FORMAT(1X,'MAXI ORBITALE',F15.8)
      RETURN
      END
      FUNCTION OVLAP(NI,NJ,ZI,ZJ)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 F(30),ZI,ZJ,R
      OVLAP=F(NI+NJ+1)*(ZI+ZJ)**(-NI-NJ-1)
      RETURN
      ENTRY RECOUV(NI,NJ,ZI,ZJ,R)
      K=NI+NJ
      T=0.
      RR=R*(ZI+ZJ)
      IF(RR.GT.174.) GO TO 5
      T=F(K+1)*(ZI+ZJ)**(-K-1)*DEXP(-RR)
5     CONTINUE
      S=T
      DO 10 L=1,K
      T=T*RR/DFLOAT(L)
10    S=S+T
      RECOUV=S
      RETURN
       ENTRY FACTOR (IL )
      F(1)=1.D0
      DO 20 I=2,30
20     F(I)=DFLOAT(I-1)*F(I-1)
      L=IL
      FACTOR=L
      RETURN
      END
      SUBROUTINE CAIH(APQ,ZO,AO,NO,XU,CPSO,PSEO,HV,WIJ,AIAJ,AIH,L,
     1NTL,NBAS1,VAV,AV,VHV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VAV(22)
      COMMON/SUP/FSUP(22,22),ESUP(22),VHSUPV(22),NSUP
      DIMENSION WIJ(NBAS1,NBAS1,NTL)
      DIMENSION XU(6),AVP(22),CPSO(22),APQ(22,22),ZO(22),AO(22),NO(22)
      DIMENSION AV(22,22),HV(22),AIH(22),AIAJ(22,22)
      COMMON/PSEUD/ALFCT(30),CACT(30),ALFA(30,4),CA(30,4),PSE(4),
     1ICT(30),ISG(30,4),NT(4),NTCT,NCT(30),NX(30,4)
      INB=1
      DO 50 IT=1,NTL
      ISGIL=ISG(IT,L)
      KX=NX(IT,L)
      IF(IT.EQ.1) GO TO 51
      IF(ALFA(IT,L).EQ.ALFA(IT-1,L)) GO TO 52
      INB=INB+1
   51 AX=XU(INB)
52    DO 60 I=1,NBAS1
      DO 60 J=1,I
      WIJ(I,J,IT)=WPS(ZO(I),AO(I),NO(I),ZO(J),AO(J),NO(J),ISGIL,KX,AX)
   60 WIJ(J,I,IT)=WIJ(I,J,IT)
      SS=0.D0
      DO 65 IP=1,NBAS1
      SOM=0.D0
      DO 66 IQ=1,NBAS1
   66 SOM=SOM+CPSO(IQ)*WIJ(IP,IQ,IT)
      AVP(IP)=SOM
   65 SS=SS+CPSO(IP)*AVP(IP)
      VAV(IT)=SS
      DO 67 IP=1,NBAS1
      SOM=0.D0
      DO 68 IQ=1,NBAS1
   68 SOM=SOM+AVP(IQ)*APQ(IQ,IP)
   67 AV(IP,IT)=SOM
      SOM=0.D0
      DO 69 IP=1,NBAS1
      SOM=SOM+AV(IP,IT)*HV(IP)
69    CONTINUE
      AIH(IT)=AIH(IT)+PSEO*VAV(IT)-SOM
   50 CONTINUE
      DO 70 IT=1,NTL
      DO 70 IS=1,IT
      SOM=0.D0
      DO 72 IP=1,NBAS1
      SOM=SOM+AV(IP,IT)*AV(IP,IS)
72    CONTINUE
      AIAJ(IT,IS)=AIAJ(IT,IS)+SOM
   70 AIAJ(IS,IT)=AIAJ(IT,IS)
      IF(NSUP.EQ.0) RETURN
      NTL1=NTL+1
      NTL2=NTL+NSUP
      DO 130 I=1,NSUP
      DO 130 J=1,NTL
      TOT=0.
      DO 120 IP=1,NBAS1
      SOM=0.
      DO 110 IQ=1,NBAS1
  110 SOM=SOM+FSUP(IQ,I)*WIJ(IQ,IP,J)
  120 TOT=TOT+SOM*FSUP(IP,I)
      TOT=TOT-VAV(J)
      AIAJ(NTL+I,J)=TOT
  130 AIAJ(J,NTL+I)=TOT
      DO 140 I=NTL1,NTL2
      II=I-NTL
      DO 135 J=NTL1,NTL2
  135 AIAJ(J,I)=0.
  140 AIH(I)=(ESUP(II)-PSEO)-(VHSUPV(II)-VHV)
      WRITE(6,*) (VAV(I),I=1,NTL)
       WRITE(6,*) PSEO,VHV,VHSUPV(1),VHSUPV(2)
      DO 150 I=1,NTL2
  150 WRITE(6,*) (AIAJ(I,J),J=1,NTL2),AIH(I)
      RETURN
      END
      SUBROUTINE CAWIJ(WIJ,ZO,NO,XU,L,NTL,NBAS1,AO)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WIJ(NBAS1,NBAS1,NTL),ZO(1),NO(1),AO(1),XU(1)
      COMMON/PSEUD/ALFCT(30),CACT(30),ALFA(30,4),CA(30,4),PSE(4),
     &ICT(30),ISG(30,4),NT(4),NTCT,NCT(30),NX(30,4)
      INB=1
      DO 50 IT=1,NTL
      ISGIL=ISG(IT,L)
      KX=NX(IT,L)
      IF(IT.EQ.1) GO TO 51
      IF(ALFA(IT,L).EQ.ALFA(IT-1,L))GO TO 52
      INB=INB+1
51    AX=XU(INB)
52    DO 60 I=1,NBAS1
      DO 60 J=1,I
      WIJ(I,J,IT)=WPS(ZO(I),AO(I),NO(I),ZO(J),AO(J),NO(J),ISGIL,KX,AX)
60    WIJ(J,I,IT)=WIJ(I,J,IT)
50     CONTINUE
      RETURN
      END
      SUBROUTINE EXTRAC
      IMPLICIT REAL*8(A-H,O-Z)
      WRITE(6,1)
1     FORMAT(' ATTENTION PSEUDOS NON LOCAUX')
      STOP
      END
      SUBROUTINE DVA0XX (X,E,N,F,ESCALE,IPRINT,ICON,MAXIT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(1),X(1),E(1)
      COMMON W
      DDMAG=0.1*ESCALE
      SCER=0.05/ESCALE
      JJ=N*N+N
      JJJ=JJ+N
      K=N+1
      NFCC=1
      DO 1 I=1,N
      DO 2 J=1,N
      W(K)=0.
      IF(I-J)4,3,4
    3 W(K)=DABS(E(I))
      W(I)=ESCALE
    4 K=K+1
    2 CONTINUE
    1 CONTINUE
      ITERC=1
      ISGRAD=2
      CALL CALCFX(N,X,F)
      FKEEP=DABS(F)+DABS(F)
    5 ITONE=1
      FP=F
      SUM=0.
      IXP=JJ
      DO 6 I=1,N
      IXP=IXP+1
      W(IXP)=X(I)
    6 CONTINUE
      IDIRN=N+1
      ILINE=1
    7 DMAX=W(ILINE)
      DACC=DMAX*SCER
      DMAG= DMIN1(DDMAG,0.1*DMAX)
       DMAG = DMAX1(DMAG,20.*DACC)
      DDMAX=10.*DMAG
      GO TO (70,70,71),ITONE
   70 DL=0.
      D = DMAG
      FPREV=F
      IS=5
      FA=F
      DA=DL
    8 DD=D-DL
      DL=D
   58 K=IDIRN
      DO 9 I=1,N
      X(I)=X(I)+DD*W(K)
      K=K+1
    9 CONTINUE
      CALL CALCFX(N,X,F)
      NFCC=NFCC+1
      GO TO (10,11,12,13,14,112),IS
   14 IF(F-FA)15,16,24
   16 IF(DABS(D)-DMAX)  17,17,18
   17 D=D+D
      GO TO 8
   18 WRITE(6,19)
   19 FORMAT(5X,44HVA04A MAXIMUM CHANGE DOES NOT ALTER FUNCTION)
      GO TO 20
   15 FB=F
      DB=D
      GO TO 21
   24 FB=FA
      DB=DA
      FA=F
      DA=D
   21 GO TO (83,23),ISGRAD
   23 D=DB+DB-DA
      IS=1
      GO TO 8
   83 D=0.5*(DA+DB-(FA-FB)/(DA-DB))
      IS=4
      IF((DA-D)*(D-DB))25,8,8
   25 IS=1
      IF(DABS(D-DB)-DDMAX) 8,8,26
   26 D=DB+DSIGN(DDMAX,DB-DA)
      IS=1
      DDMAX=DDMAX+DDMAX
      DDMAG=DDMAG+DDMAG
      IF(DDMAX-DMAX)8,8,27
   27 DDMAX=DMAX
      GO TO 8
   13 IF(F-FA)28,23,23
   28 FC=FB
      DC=DB
   29 FB=F
      DB=D
      GO TO 30
   12 IF(F-FB)28,28,31
   31 FA=F
      DA=D
      GO TO 30
   11 IF(F-FB)32,10,10
   32 FA=FB
      DA=DB
      GO TO 29
   71 DL=1.
      DDMAX=5.
      FA=FP
      DA=-1.
      FB=FHOLD
      DB=0.
      D=1.
   10 FC=F
      DC=D
   30 A=(DB-DC)*(FA-FC)
      B=(DC-DA)*(FB-FC)
      IF((A+B)*(DA-DC))33,33,34
   33 FA=FB
      DA=DB
      FB=FC
      DB=DC
      GO TO 26
   34 D=0.5*(A*(DB+DC)+B*(DA+DC))/(A+B)
      DI=DB
      FI=FB
      IF(FB-FC)44,44,43
   43 DI=DC
      FI=FC
   44 GO TO (86,86,85),ITONE
   85 ITONE=2
      GO TO 45
  86  IF(DABS(D-DI)-DACC)41,41,93
  93  IF(DABS(D-DI)-0.03*DABS(D))41,41,45
   45 IF ((DA-DC)*(DC-D)) 47,46,46
   46 FA=FB
      DA=DB
      FB=FC
      DB=DC
      GO TO 25
   47 IS=2
      IF ((DB-D)*(D-DC)) 48,8,8
   48 IS=3
      GO TO 8
   41 F=FI
      D=DI-DL
      DD=DSQRT((DC-DB)*(DC-DA)*(DA-DB)/(A+B))
      DO 49 I=1,N
      X(I)=X(I)+D*W(IDIRN)
      W(IDIRN)=DD*W(IDIRN)
      IDIRN=IDIRN+1
   49 CONTINUE
      W(ILINE)=W(ILINE)/DD
      ILINE=ILINE+1
      IF(IPRINT-1)51,50,51
   50 WRITE(6,52)ITERC,NFCC,F,(X(I),I=1,N)
   52 FORMAT (/1X,9HITERATION,I5,I15,16H FUNCTION VALUES,
     110X,3HF =,E21.14/(5E24.14))
      GO TO (51,109),IPRINT
   51 GO TO (55,38),ITONE
   55 IF (FPREV-F-SUM) 94,95,95
   95 SUM=FPREV-F
      JIL=ILINE
   94 IF (IDIRN-JJ) 7,7,92
   92 FHOLD=F
      IS=6
      IXP=JJ
      DO 59 I=1,N
      IXP=IXP+1
      W(IXP)=X(I)-W(IXP)
   59 CONTINUE
      DD=1.
      GO TO 58
  112 IF (FP-F) 37,91,91
   91 D=2.*(FP+F-2.*FHOLD)/(FP-F)**2
      IF (D*(FP-FHOLD-SUM)**2-SUM) 87,37,37
   87 J=JIL*N+1
      IF (J-JJ) 60,60,61
   60 DO 62 I=J,JJ
      K=I-N
      W(K)=W(I)
   62 CONTINUE
      DO 97 I=JIL,N
      W(I-1)=W(I)
   97 CONTINUE
   61 IDIRN=IDIRN-N
      ITONE=3
      K=IDIRN
      IXP=JJ
      AAA=0.
      DO 65 I=1,N
      IXP=IXP+1
      W(K)=W(IXP)
      IF (AAA-DABS(W(K)/E(I))) 66,67,67
   66 AAA=DABS(W(K)/E(I))
   67 K=K+1
   65 CONTINUE
      DDMAG=1.
      W(N)=ESCALE/AAA
      ILINE=N
      GO TO 7
   37 IXP=JJ
      AAA=0.
      F=FHOLD
      DO 99 I=1,N
      IXP=IXP+1
      X(I)=X(I)-W(IXP)
      IF (AAA*DABS(E(I))-DABS(W(IXP))) 98,99,99
 98   AAA=DABS(W(IXP)/E(I))
   99 CONTINUE
      GO TO 72
   38 AAA=AAA*(1.+DI)
  72  IF(IPRINT-2)109,50,50
  109 IF(AAA-0.1)20,20,76
76    IF(F-FP)35,78,78
   78 WRITE(6,80)
   80 FORMAT (5X,37HVA04A ACCURACY LIMITED BY ERRORS IN F)
      GO TO 20
   35 DDMAG=0.4*DSQRT(FP-F)
      ISGRAD=1
  108 ITERC=ITERC+1
      IF (ITERC-MAXIT) 5,5,81
   81 WRITE(6,82)MAXIT
   82 FORMAT(I5,30H ITERATIONS COMPLETED BY VA04A)
      IF (F-FKEEP) 20,20,110
  110 F=FKEEP
      DO 111 I=1,N
      JJJ=JJJ+1
      X(I)=W(JJJ)
  111 CONTINUE
   20 RETURN
      END
      SUBROUTINE DVLVC1(A,NA,MA,VP,NVP,C,NC,MC,N,IA)
      DOUBLE PRECISION A(NA,MA),VP(NVP),C(NC,MC)
      DOUBLE PRECISION E,SIN45,COS45,S45SQ,C45SQ,AMAX,TEMPA,TEST,DSQRT
      DOUBLE PRECISION AIJMAX,TAIJ, APP,AQQ,APQ,TMT,Z,SINT,COST,Z2
      DOUBLE PRECISION SINSQ ,COSSQ,APK,AQK,CPK,CQK,VPPROV
      INTEGER P,Q
C
C                    METHODE DE  JACOBI.
C
C    DIAGONALISATION D'UNE MATRICE SYMETRIQUE A ,D'ORDRE N
C      RECHERCHE DU PLUS GRAND ELEMENT  EXTRADIAGONAL
C A LA   FIN   DU CALCUL
C      LA  MATRICE  A  EST  DETRUITE
C IA  CONTIENT  LE NOMBRE D'ITERATIONS EFFECTUEES
C LES VALEURS PROPRES SONT DANS  VP(I).
C        ELLES SONT RANGEES DANS L'ORDRE ALGEBRIQUE CROISSANT
C LES VECTEURS PROPRES SONT DANS  C(I,J)
C        LA 1 ERE LIGNE DE C  REPRESENTE  LE VECTEUR
C    PROPRE ASSOCIE  A LA PLUS PETITE  VALEUR PROPRE
C
      E=1.D-20
      ITER=0
      SIN45= ( DSQRT(2.D0) )/2.D0
      COS45=SIN45
      S45SQ=0.5D0
      C45SQ=S45SQ
C
C       INITIALISATION
C
      DO 70 I = 1,N
       DO 7 J=1,N
 7    C(I,J)=0.D0
 70   C(I,I)=1.D0
C
C      DETERMINATION DU TEST D'ARRET
C
      AMAX=0.D0
      DO 1 I = 1,N
      DO 1 J = 1,I
      TEMPA=DABS(A(I,J))
      IF (AMAX-TEMPA) 2,1,1
    2 AMAX = TEMPA
    1 CONTINUE
      TEST = AMAX*E
C
C       RECHERCHE DU PLUS GRAND ELEMENT EXTRADIAGONAL
C
 72   AIJMAX=0.D0
      ITER=ITER+1
      DO 3 I = 2,N
      LIM = I-1
      DO 3 J = 1,LIM
       TAIJ=DABS(A(I,J))
       IF (AIJMAX-TAIJ) 4,3,3
    4 AIJMAX = TAIJ
      P=I
      Q=J
    3 CONTINUE
      IF(AIJMAX-TEST)300,300,5
C
C       CALCUL  DES ELEMENTS DE LA MATRICE T
C
 5    APP=A(P,P)
      AQQ=A(Q,Q)
      APQ=A(P,Q)
      TMT=APP-AQQ
      IF(DABS(TMT/APQ) - 1.D-9) 60,60,6
 60   IF(APQ) 10,10,11
 6    Z=APQ/(2.D0*TMT)
 90   IF(DABS(Z) - 0.38268D0) 8,8,9
 9    IF(Z) 10,10,11
   10 SINT = -SIN45
      GO TO 12
   11 SINT = SIN45
   12 COST = COS45
      SINSQ = S45SQ
      COSSQ = C45SQ
      GO TO 120
 8    Z2=Z*Z
      SINT=2.D0*Z/( 1.D0+Z2)
      COST=(1.D0 - Z2) / (1.D0 +Z2 )
      SINSQ=SINT*SINT
      COSSQ=COST*COST
C
C     CALCUL DES NOUVEAUX ELEMENTS DE LA MATRICE
C    CALCUL DES VECTEURS PROPRES
C
  120 DO 13 K = 1,N
      APK=A(P,K)
      AQK=A(Q,K)
      A(P,K)=APK*COST + AQK*SINT
      A(Q,K)=AQK*COST - APK*SINT
      CPK=C(P,K)
      CQK=C(Q,K)
      C(P,K)=CPK*COST + CQK*SINT
 13   C(Q,K)=CQK*COST - CPK*SINT
      A(P,P)=APP*COSSQ + AQQ*SINSQ +2.D0*APQ*SINT*COST
      A(Q,Q)=APP*SINSQ + AQQ*COSSQ - 2.D0*APQ*SINT*COST
      A(P,Q)=APQ*(COSSQ - SINSQ) - SINT*COST*TMT
      A(Q,P)=A(P,Q)
      DO 30 K = 1,N
      A(K,P)=A(P,K)
 30   A(K,Q)=A(Q,K)
      GO TO 72
C
C     CLASSEMENT DES VALEURS PROPRES PAR ODRE CROISSANT
C
  300 DO 14 I=1,N
 14   VP(I)=A(I,I)
      L=N-1
      DO 154 K=1,L
      KK=K+1
      DO 152 KJ=KK,N
      IF( VP(KJ). GE .VP(K) )   GOTO 152
      DO 153 J=1,N
      A(1,J)=C(K,J)
      C(K,J)=C(KJ,J)
      C(KJ,J)=A(1,J)
 153  CONTINUE
      VPPROV=VP(K)
      VP(K)=VP(KJ)
      VP(KJ)=VPPROV
 152  CONTINUE
 154  CONTINUE
      IA=ITER
      RETURN
      END
      DOUBLE PRECISION FUNCTION DGLEG4(F,A,B)
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION H(20),X(20)
      DATA X( 1)/  -0.9982377097105589D 00/
      DATA H( 1)/   0.4521277098533191D-02/
      DATA X( 2)/  -0.9907262386994567D 00/
      DATA H( 2)/   0.1049828453115281D-01/
      DATA X( 3)/  -0.9772599499837742D 00/
      DATA H( 3)/   0.1642105838190788D-01/
      DATA X( 4)/  -0.9579168192137913D 00/
      DATA H( 4)/   0.2224584919416695D-01/
      DATA X( 5)/  -0.9328128082786765D 00/
      DATA H( 5)/   0.2793700698002340D-01/
      DATA X( 6)/  -0.9020988069688741D 00/
      DATA H( 6)/   0.3346019528254784D-01/
      DATA X( 7)/  -0.8659595032122592D 00/
      DATA H( 7)/   0.3878216797447202D-01/
      DATA X( 8)/  -0.8246122308333116D 00/
      DATA H( 8)/   0.4387090818567327D-01/
      DATA X( 9)/  -0.7783056514265193D 00/
      DATA H( 9)/   0.4869580763507223D-01/
      DATA X(10)/  -0.7273182551899268D 00/
      DATA H(10)/   0.5322784698393678D-01/
      DATA X(11)/  -0.6719566846141792D 00/
      DATA H(11)/   0.5743976909939152D-01/
      DATA X(12)/  -0.6125538896679799D 00/
      DATA H(12)/   0.6130624249292891D-01/
      DATA X(13)/  -0.5494671250951279D 00/
      DATA H(13)/   0.6480401345660100D-01/
      DATA X(14)/  -0.4830758016861787D 00/
      DATA H(14)/   0.6791204581523383D-01/
      DATA X(15)/  -0.4137792043716050D 00/
      DATA H(15)/   0.7061164739128673D-01/
      DATA X(16)/  -0.3419940908257584D 00/
      DATA H(16)/   0.7288658239580400D-01/
      DATA X(17)/  -0.2681521850072536D 00/
      DATA H(17)/   0.7472316905796821D-01/
      DATA X(18)/  -0.1926975807013711D 00/
      DATA H(18)/   0.7611036190062619D-01/
      DATA X(19)/  -0.1160840706752552D 00/
      DATA H(19)/   0.7703981816424793D-01/
      DATA X(20)/  -0.3877241750605079D-01/
      DATA H(20)/   0.7750594797842478D-01/
      FXI=0.D0
       BMA=(B-A)*0.5
       BPA=(B+A)*0.5
       DO 1 I=1,20
       XI=BMA*X(I)+BPA
 1    FXI=F(XI)*H(I)+FXI
       DO 2 I=21,40
       J=41-I
       XI=BPA-BMA*X(J)
 2     FXI=F(XI)*H(J)+FXI
       DGLEG4=BMA*FXI
      RETURN
      END
       DOUBLE PRECISION FUNCTION DSLAG(N,A)
C
C
C     LA FONCTION DSLAG CALCULE L'INTEGRALE :
C     X**N*EXP(-X**2-A*X)
C     ENTRE  ZERO ET +INFINI.
C     A DOIT ETRE SUP.OU EGAL A -26.
C     N DOIT ETRE INF.OU EGAL A 20
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 J0,J(20)
      DIMENSION BORN1(21),BORN2(21),C(30),FACT(16)
       DATA QL2    /    .2772588722239782D+01/
      DATA FACT                   / .10000000000D+01, .20000000000D+01,
     1            .60000000000D+01, .24000000000D+02, .12000000000D+03,
     2            .72000000000D+03, .50400000000D+04, .40320000000D+05,
     3            .36288000000D+06, .36288000000D+07, .39916800000D+08,
     4            .47900160000D+09, .62270208000D+10, .87178291200D+11,
     5            .13076743680D+13, .20922789888D+14/
      DATA BORN1/10.0,10.0, 7.0, 5.0, 5.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
     1            4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0/
      DATA BORN2/35.0,28.0,16.0,12.0, 9.0, 8.0, 8.0, 8.0, 8.0, 8.0, 5.0,
     1            5.0, 5.0, 5.0, 5.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0/
      DATA C/    .5000000000000000D+00,    .4431134627263790D+00,
     1           .5000000000000000D+00,    .6646701940895684D+00,
     2           .1000000000000000D+01,    .1661675485223921D+01,
     3           .3000000000000000D+01,    .5815864198283724D+01,
     4           .1200000000000000D+02,    .2617138889227676D+02,
     5           .6000000000000000D+02,    .1439426389075221D+03,
     6           .3600000000000000D+03,    .9356271528988939D+03,
     7           .2520000000000000D+04,    .7017203646741704D+04,
     8           .2016000000000000D+05,    .5964623099730448D+05,
     9           .1814400000000000D+06,    .5666391944743925D+06,
     *           .1814400000000000D+07,    .5949711541981121D+07,
     *           .1995840000000000D+08,    .6842168273278289D+08,
     *           .2395008000000000D+09,    .8552710341597860D+09,
     *           .3113510400000000D+10,    .1154615896115711D+11,
     *           .4358914560000000D+11,    .1674193049367781D+12/
      DATA SQPIS2/0.886226925452758D0/
        COMMON/ETDSL/AA,NN
       NN=N
       AA=A
       AN=N
      NP1=N+1
      B1=BORN1(NP1)
      B2=BORN2(NP1)
       DSLAG=0.D0
      IF(A.LE.0.) GOTO 2030
C
C------------------ CAS OU A EST POSITIF
C
       IF(A.LE.0.1D0) GOTO 2021
      DL=A*A+QL2
       XMAX=(-A+DSQRT(A*A+8.D0*AN))*0.25
       IF(N.EQ.0) XMAX=( DSQRT(DL)-A )*0.5
      IF(A.GT.1.) GOTO 3
      BSUP=B1*XMAX
       GOTO 1
 3    IF(A.GE.20.) GOTO 5
      BSUP=( B1+(B2-B1)*DLOG10(A)*0.77 )*XMAX
       GOTO 1
 5    BSUP=B2*XMAX
 1     DSLAG=DGLEGX(BSUP)
       RETURN
C
C------------------  A.LE.0.1
C
 2021  CONTINUE
       SIG=-1.D0
       DO 200 K=1,10
       DSLAG=DSLAG+SIG*A**K*C(K+N)/        FACT(K)
       SIG=-SIG
 200   CONTINUE
       IF(N-1)11,12,13
C
C      CAS N=0
C
 11    DSLAG=DSLAG+SQPIS2
       GOTO 300
C
C      CAS N=1
C
 12    DSLAG=DSLAG+0.5D0
       GOTO 300
C
C      CAS N.GT.1
C
 13    DSLAG=DSLAG+C(N)
 300   CONTINUE
       RETURN
C
C---------------------  CAS OU A EST NEGATIF OU NUL
C
 2030  CONTINUE
      EX=DEXP(A*A*0.25)
      AS2=DABS(A*0.5)
       DERFCA=1.D0+DERFA(AS2)
       J0=SQPIS2*EX*DERFCA
       IF(N-1)211,212,213
C
C      CAS N=0
C
 211   DSLAG=J0
       GOTO 301
C
C      CAS N=1
C
 212  J(1)=(1.D0-A*J0)*0.5
       DSLAG=J(1)
       GOTO 301
C
C      CAS N.GT.1
C
 213  J(1)=(1.D0-A*J0)*0.5
      J(2)=(-A*J(1)+J0)*0.5
       DO 2130 I=3,N
      DIM1=I-1
      J(I)=( -A*J(I-1)+DIM1*J(I-2) )*0.5
 2130  CONTINUE
       DSLAG=J(N)
 301   CONTINUE
 30    CONTINUE
       RETURN
       END
      DOUBLE PRECISION FUNCTION DGLEGX(BSUP)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(20),X(20)
       COMMON/ETDSL/AA,NN
      DATA X( 1)/  -0.9982377097105589D 00/
      DATA H( 1)/   0.4521277098533191D-02/
      DATA X( 2)/  -0.9907262386994567D 00/
      DATA H( 2)/   0.1049828453115281D-01/
      DATA X( 3)/  -0.9772599499837742D 00/
      DATA H( 3)/   0.1642105838190788D-01/
      DATA X( 4)/  -0.9579168192137913D 00/
      DATA H( 4)/   0.2224584919416695D-01/
      DATA X( 5)/  -0.9328128082786765D 00/
      DATA H( 5)/   0.2793700698002340D-01/
      DATA X( 6)/  -0.9020988069688741D 00/
      DATA H( 6)/   0.3346019528254784D-01/
      DATA X( 7)/  -0.8659595032122592D 00/
      DATA H( 7)/   0.3878216797447202D-01/
      DATA X( 8)/  -0.8246122308333116D 00/
      DATA H( 8)/   0.4387090818567327D-01/
      DATA X( 9)/  -0.7783056514265193D 00/
      DATA H( 9)/   0.4869580763507223D-01/
      DATA X(10)/  -0.7273182551899268D 00/
      DATA H(10)/   0.5322784698393678D-01/
      DATA X(11)/  -0.6719566846141792D 00/
      DATA H(11)/   0.5743976909939152D-01/
      DATA X(12)/  -0.6125538896679799D 00/
      DATA H(12)/   0.6130624249292891D-01/
      DATA X(13)/  -0.5494671250951279D 00/
      DATA H(13)/   0.6480401345660100D-01/
      DATA X(14)/  -0.4830758016861787D 00/
      DATA H(14)/   0.6791204581523383D-01/
      DATA X(15)/  -0.4137792043716050D 00/
      DATA H(15)/   0.7061164739128673D-01/
      DATA X(16)/  -0.3419940908257584D 00/
      DATA H(16)/   0.7288658239580400D-01/
      DATA X(17)/  -0.2681521850072536D 00/
      DATA H(17)/   0.7472316905796821D-01/
      DATA X(18)/  -0.1926975807013711D 00/
      DATA H(18)/   0.7611036190062619D-01/
      DATA X(19)/  -0.1160840706752552D 00/
      DATA H(19)/   0.7703981816424793D-01/
      DATA X(20)/  -0.3877241750605079D-01/
      DATA H(20)/   0.7750594797842478D-01/
      BMA=BSUP*0.5
      FXI=0.D0
      DO 1 I=1,20
      XI=BMA*(X(I)+1.)
 1    FXI=FXI+XI**NN*DEXP( -XI*(XI+AA) )*H(I)
      DO 2 I=21,40
       J=41-I
      XI=BMA*(-X(J)+1.)
 2    FXI=FXI+XI**NN*DEXP( -XI*(XI+AA) )*H(J)
      DGLEGX=BMA*FXI
      RETURN
      END
      SUBROUTINE COUTR(AIAJ,AIH,APOT,NPOT,NTL,VAV,VXV,AV)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/UN/TAF(22,22,3),APQ(22,22,3),ZO(22,3),AO(22,3),PSEO(3),
     1CPSO(22,3),NO(22,3),HV(22,3),VHV(3),L,NBAS1,NETAT,NEND
      DIMENSION SJEU(22,22),SSAV(22,22),HSAV(22),AIAJ(22,22),AIH(22),
     1APOT(22),NPOT(22),SSD(22),SSVP(22),SSV(22,22),CIK(22),TAB(22),
     2VAV(22,3),AV(22,22,3),SS(22),SVQ(22,3)
      DIMENSION TJEU(22,22),ZJEU(22),XJEU(22)
      LOGICAL YCAS
      LOGICAL YPAS
      YPAS=.FALSE.
      DO 705 I=1,NBPOT
      DO 705 J=1,I
705   SJEU(I,J)=0.D0
      DO 16 IK=1,NETAT
      DO 11 J1=1,NBAS1
      N1=NO(J1,IK)
      Z1=ZO(J1,IK)
      REC1=0.D0
      DO 10 J2=1,NBAS1
      N2=NO(J2,IK)
      Z2=ZO(J2,IK)
      TOT=AO(J2,IK)*AO(J1,IK)*CPSO(J2,IK)*PURX(N1,N2,Z1,Z2,0,0.D0)
      REC1=TOT+REC1
10    CONTINUE
      SS(J1)=REC1
11    CONTINUE
      DO 15 J1=1,NBAS1
      T=0.
      DO 14 J2=1,NBAS1
14    T=T+SS(J2)*APQ(J2,J1,IK)
15    SVQ(J1,IK)=T
16       CONTINUE
      CALL LININT(R0,RDIV,NDIV,PASLAM,FMUL,NITER,ORILAM)
      NBPOT=NTL
      DO 770 IK=1,NETAT
      DO 750 J=1,NTL
      DO 730 K=1,NBAS1
      DO 720 I=1,NBAS1
      NPOT1=NPOT(J)+NO(I,IK)+2+NO(K,IK)
      APOT1=APOT(J)+ZO(I,IK)+ZO(K,IK)
      CALL LINT(NPOT1,APOT1,S,R0,RDIV,NDIV)
      S=S*AO(I,IK)*AO(K,IK)
720   ZJEU(I)=S
      S=0.
      DO 725 I=1,NBAS1
725    S=S+CPSO(I,IK)*ZJEU(I)
730     XJEU(K)=S
      DO 740 I=1,NBAS1
      S=0.D0
      DO 735 K=1,NBAS1
735   S=S+APQ(K,I,IK)*XJEU(K)
740   TJEU(I,J)=S
750   CONTINUE
      SMAX=0.D0
      IMAX=0
      DO 191 I=1,NBPOT
      DO 190 J=1,I
       S=0.
      DO 765 K=1,NBAS1
765   S=S+TJEU(K,I)*TJEU(K,J)
      SJEU(I,J)=SJEU(I,J)+S
190   CONTINUE
      WRITE(6,*) (SJEU(I,J),J=1,I)
      WRITE(6,*) (AIAJ(I,J),J=1,I)
      IF(SMAX.GT.S) GO TO 191
      SMAX=S
      IMAX=I
191   CONTINUE
770   CONTINUE
      YCAS=.TRUE.
      FMUL=1.
      IF(IMAX.NE.0) FMUL=AIAJ(IMAX,IMAX)/SMAX
      WRITE(6,*) IMAX,SMAX,FMUL
      GNORMI=100.
305      CONTINUE
      DLAMB=ORILAM
      DO 300 ITER=1,NITER
      QLAMB=DLAMB*FMUL
      DO 211 I=1,NTL
      HSAV(I)=AIH(I)
      DO 210 J=1,I
210   SSAV(I,J)=AIAJ(I,J)+QLAMB*SJEU(I,J)
211     CONTINUE
      DO 19 I=1,NTL
19    SSD(I)=DSQRT(SSAV(I,I))
      DO 21 I=1,NTL
      HSAV(I)=HSAV(I)/SSD(I)
      DO 20 J=1,I
      SSAV(I,J)=SSAV(I,J)/(SSD(I)*SSD(J))
20    SSAV(J,I)=SSAV(I,J)
21           CONTINUE
      CALL DVLVC1(SSAV,22,22,SSVP,22,SSV,22,22,NTL,IIT)
      DO 40 I=1,NTL
      SSVP(I)=1.D0/SSVP(I)
      IF(SSVP(I).LT.0.D0) SSVP(I)=0.D0
40    IF(SSVP(I).GT.1.D10) SSVP(I)=0.D0
      DO 60 I=1,NTL
      DO 60 J=1,I
      TOT=0.D0
      DO 50 IK=1,NTL
50    TOT=TOT+SSV(IK,I)*SSV(IK,J)*SSVP(IK)
      SSAV(I,J)=TOT
60    SSAV(J,I)=TOT
      DO 80 I=1,NTL
      TOT=0.D0
      DO 70 J=1,NTL
70    TOT=TOT+SSAV(I,J)*HSAV(J)
80    SSD(I)=TOT/SSD(I)
      DO 90 I=1,NTL
90    HSAV(I)=SSD(I)
      TT=0.D0
      DO 220 I=1,NTL
      T=0.D0
      DO 201 J=1,I
201   T=T+HSAV(J)*SJEU(I,J)
220   TT=TT+T*HSAV(I)
      GNOR=0.
      DO 130 IK=1,NETAT
      GNOR=PSEO(IK)**2+GNOR
      TOT=0.
      SUM=0.
      TUT=0.
      DO 401 IP=1,NBAS1
      T=0.
      DO 402 IT=1,NTL
402   T=T+AV(IP,IT,IK)*HSAV(IT)
      SUM=SUM+T*SVQ(IP,IK)
      TOT=TOT+HV(IP,IK)*HV(IP,IK)
401   TUT=TUT+HV(IP,IK)*T
      GNOR=GNOR+TOT-2.D0*PSEO(IK)*VHV(IK)
      GNOR=GNOR-SUM*PSEO(IK)+TUT
130   CONTINUE
      TNOR=1.
      IF(TNOR.GT.1.D-3) TNOR=100.
      IF(GNOR.LT.1.D-5) TNOR=GNOR*1.D5
      ANOR=TT+GNOR*TNOR
      WRITE(6,*) DLAMB,TT,GNOR,ANOR
140   ANORO=ANOR
      OLAMB=DLAMB
      DLAMB=DLAMB*PASLAM
      YCAS=.NOT.YCAS
      IF(ANOR.GT.GNORMI) GO TO 300
      ORILAM=DLAMB/PASLAM
       GNORMI=ANOR
300   CONTINUE
      IF(YPAS) GO TO 320
      YPAS=.TRUE.
      NITER=1
      GO TO 305
320    CONTINUE
330   DO 350 I=1,NTL
350   AIH(I)=HSAV(I)
      WRITE(6,*) (HSAV(I),I=1,NTL)
      WRITE(6,30) (SSVP(I),I=1,NTL)
30    FORMAT(' VALEUR PROPRE DE LA MATRICE DES OPERATEUR:',(1X,6E12.2))
      WRITE(6,*) (SSVP(I),I=1,NTL)
      RA=1.D-6
      DO 250 I=1,500,10
      DO 245 K=1,10
      T=0.
      DO 240 J=1,NTL
      EXPON=DMIN1(173.D0,RA*RA*APOT(J))
240   T=T+HSAV(J)*RA**NPOT(J)*DEXP(-EXPON)
      IF(DABS(T).LT.1.D-6.AND.RA.GT.5) GO TO 250
      RA=RA+.2
245   TAB(K)=T
      WRITE(6,246) TAB
246   FORMAT(1X,10F12.6)
226   FORMAT(' TABULATION R=R0 A R0+2.')
250   CONTINUE
      RETURN
      END
      SUBROUTINE LINT(NPOT,APOT,S,R0,RDIV,NDIV)
      IMPLICIT REAL*8(A-H,O-Z)
      S=0.D0
      Z=R0
      ZZ=Z*Z
      X=R0+RDIV
      XX=X*X
      TPOT=DMIN1(173.D0,APOT*ZZ)
      Y1=Z**NPOT*DEXP(-TPOT)
      XDIV=RDIV
      DO 11 K=1,NDIV
      SPOT=DMIN1(173.D0,APOT*XX)
      Y2=X**NPOT*DEXP(-SPOT)
      S=S+(Y1+Y2)*XDIV
      XDIV=XDIV*1.1
      X=X+XDIV
      XX=X*X
      Y1=Y2
11    CONTINUE
      S=S*.5D0
      RETURN
      ENTRY LININT(R0,RDIV,NDIV,PASLAM,FMUL,NITER,ORILAM)
      READ(5,*) R0,RDIV,NDIV,FMUL,PASLAM,NITER,ORILAM
      WRITE(6,*) R0,RDIV,NDIV,FMUL,PASLAM,NITER,ORILAM
      RETURN
      END
      DOUBLE PRECISION FUNCTION DERFA(X)
      DOUBLE PRECISION X,C(15),A(21),T(43),XS3,F,S,DSQRT,DEXP
      DATA C( 1)/    0.975083423708555D+00/
      DATA C( 2)/   -0.240493938504140D-01/
      DATA C( 3)/    0.820452240880000D-03/
      DATA C( 4)/   -0.434293081300000D-04/
      DATA C( 5)/    0.301844703400000D-05/
      DATA C( 6)/   -0.254473319000000D-06/
      DATA C( 7)/    0.248583530000000D-07/
      DATA C( 8)/   -0.273172000000000D-08/
      DATA C( 9)/    0.330847000000000D-09/
      DATA C(10)/   -0.435050000000000D-10/
      DATA C(11)/    0.614100000000000D-11/
      DATA C(12)/   -0.922000000000000D-12/
      DATA C(13)/    0.146000000000000D-12/
      DATA C(14)/   -0.240000000000000D-13/
      DATA C(15)/    0.400000000000000D-14/
      DATA A( 1)/    0.109547129977762D+01/
      DATA A( 2)/   -0.289175401126989D+00/
      DATA A( 3)/    0.110456398633795D+00/
      DATA A( 4)/   -0.412531882278560D-01/
      DATA A( 5)/    0.140828380706510D-01/
      DATA A( 6)/   -0.432929544743100D-02/
      DATA A( 7)/    0.119827190159200D-02/
      DATA A( 8)/   -0.299972962353000D-03/
      DATA A( 9)/    0.683258603780000D-04/
      DATA A(10)/   -0.142469884540000D-04/
      DATA A(11)/    0.273540877200000D-05/
      DATA A(12)/   -0.486191287000000D-06/
      DATA A(13)/    0.803872760000000D-07/
      DATA A(14)/   -0.124184180000000D-07/
      DATA A(15)/    0.179953200000000D-08/
      DATA A(16)/   -0.245479000000000D-09/
      DATA A(17)/    0.316250000000000D-10/
      DATA A(18)/   -0.385900000000000D-11/
      DATA A(19)/    0.447000000000000D-12/
      DATA A(20)/   -0.490000000000000D-13/
      DATA A(21)/    0.500000000000000D-14/
      IF(X.GT.3.D0) GO TO 1
      XS3=X/3.D0
      T(1)=1.D0
      T(2)=XS3
      F=XS3+XS3
      DO 2 I=2,42
  2   T(I+1)=F*T(I)-T(I-1)
      S=0.D0
      DO 3 N=1,21
         IN=N-1
  3   S=S+A(N)*T(2*IN+2)
      DERFA=S*1.128379167095512D0
      GO TO 4
  1   XS3=3.00D0/X
      T(1)=1.D0
      T(2)=XS3
      F=XS3+XS3
      DO 5 I=2,29
  5   T(I+1)=F*T(I)-T(I-1)
      S=0.D0
      DO 6 N=1,15
         IN=N-1
  6   S=S+C(N)*T(2*IN+1)
      DERFA=1.D0-DEXP(-X**2)*S/(X*1.772453850905516D0)
  4   RETURN
      END
C     ******************************************************************

      subroutine openf
      character*10 prefix,f01,f08
      namelist/PSFIL/PREFIX,F01,F08
      f01='restart    '
      f08='work       '
      prefix= '          '
      read(5,PSFIL)
      write(6,psfil)
      call nomfil(1,f01,prefix,'UNKNOWN       ')
      call nomfil(8,f08,prefix,'UNKNOWN       ')
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
                                                                                
