C 
C PROGRAMME NLOCVS
C
C
      OPEN (UNIT=21,form='UNFORMATTED',status='SCRATCH')
      OPEN (UNIT=20,form='UNFORMATTED',file='PSNL.TEI',status='OLD')
      CALL NOLOC                                                        00000100
      STOP                                                              00000200
      END                                                               00000300
      SUBROUTINE NOLOC                                                  00000400
      IMPLICIT REAL*8(A-H,O-T,V-Z)                                      00000500
      DIMENSION INUM(40),EPSDS(40),Z(40),ZOP(40),RNOR(40),OPNOR(40)     00000600
      DIMENSION TOPE(1600),TORT(1600),REC(1600),TEMP(1600),VKL(840)     00000700
      DIMENSION SW(400),SV(400), SN(400),SC(400)                        00000800
      DIMENSION SL(1600)                                                00000900
      DIMENSION ZOP1(20),TEMP1(400),EPSDS1(20)                          00001000
      DIMENSION ZOP2(20),TEMP2(400),EPSDS2(20)                          00001100
      EQUIVALENCE ( ZOP (1),ZOP2(1) ), ( TEMP (1),TEMP2(1) )            00001200
      EQUIVALENCE ( EPSDS (1),EPSDS2(1) )                               00001300
      DIMENSION OLDPS(840),VERI(840)                                    00001400
      COMMON VL(44500)                                                  00001500
      DIMENSION NPOT(20),APOT(20),CPOT(20)                              00001600
      COMMON/GIV/NBLOCK(400),NSTART(400),NSIZE,IMNM(400)                00001700
      CHARACTER*8 FMT,apsd,anam                                              00001800
      DIMENSION RZ(4),OVZ(4),EX(4),OZ(2)                                00001900
      NAMELIST/INPUT/ ZN,L,NBAS,NOP,ARAF,BRAF,NVOIS,AOP,BOP,NBPOT,NPOT, 00002000
     &APOT,CPOT,Z,ZOP,NEWBAS,DMNXP2,DMNXP,OPTES,BASTES,APSD             00002100
     &,FMT,IDEN,UPUNCH,USO,NENR,UREPL,ULISTE                            00002200
     &,RZ,OVZ,UCUT,OZ                                                   00002300
     U,UPROJ                                                            00002400
     U,UOLDPS,UVERI                                                     00002500
     U,USO,NENR,NOM1,NOM2                                               00002600
      LOGICAL   NEWBAS ,UCUT,USO, UREPL,ULISTE                          00002700
      LOGICAL IVERI,UPUNCH                                              00002800
      LOGICAL   UPROJ,UOLDPS,UVERI                                      00002900
      I=INIPUR(1)                                                       00003000
      UPROJ=.FALSE.                                                     00003100
      UOLDPS=.FALSE.                                                    00003200
      UPUNCH=.FALSE.                                                    00003300
      UVERI=.FALSE.                                                     00003400
      USO=.FALSE.                                                       00003500
      UREPL=.FALSE.                                                     00003600
      DATA FMT/'(5Z16)'/                                                00003700
      MAXBAS=19                                                         00003800
      MAXBAS=20                                                         00003900
      UCUT=.FALSE.                                                      00004000
       ULISTE=.FALSE.                                                   00004100
      MDIM=400                                                          00004200
      RZ(1)=0.                                                          00004300
      OZ(1)=0.                                                          00004400
      OPTES=1.D-6                                                       00004500
      BASTES=1.D-6                                                      00004600
      DMNXP=1.D-10                                                      00004700
      ARAF=0.                                                           00004800
      IVERI=.TRUE.                                                      00004900
      DO 2300 I=1,400                                                   00005000
2300     IMNM(I)=I*(I-1)/2                                              00005100
      AOP=0.                                                            00005200
      OZ(1)=0.                                                          00005300
1     CONTINUE                                                          00005400
      DMNXP2=1.D-9                                                      00005500
      NEWBAS=.FALSE.                                                    00005600
      READ(5,input,END=999)
      WRITE(6,INPUT)                                                    00005700
C     DEFINE FILE 20 (94,882,U,NAW)                                     00005800
C     DEFINE FILE 30 (30,885,U,NAW)                                     00005900
      IF(NEWBAS) GO TO 8884                                             00006000
      WRITE(6,8888) IDEN,APSD,L,NOP,RZ,OVZ,OZ                           00006100
 8888 FORMAT(//////////,130(1H*),/,' EXTRACTION PSEUDO NON-LOCAUX POUR  00006200
     1 NUMERO ATOMIQUE',I4,10X,'ATOME',3X,A8,'SYMETRIE',I3,/,           00006300
     2 '  NOMBRE D''OPERATEURS ',I3,                                    00006400
     3 /,' EXTENSION DE LA BASE D''EXTRACTION',4F10.4,/,' RECOUVREMENTS'00006500
     4,4F10.4,/,' EXTENSION DE LA BASE D''OPERATEURS',2F10.4)           00006600
      WRITE(6,8887) NBPOT                                               00006700
      WRITE(6,8886) (NPOT(I),APOT(I),CPOT(I),I=1,NBPOT)                 00006800
 8887 FORMAT(/,'  OPERATEUR SEMI-LOCAL  NOMBRE DE TERMES=',I4)          00006900
 8886 FORMAT(I4,2F12.6)                                                 00007000
 8884 CONTINUE                                                          00007100
      LL=L+L                                                            00007200
      IF(RZ(1).EQ.0.D0) GO TO 2200                                      00007300
      DO 2050 K=1,4                                                     00007400
      EX(K)=L/(2.D0*RZ(K)*RZ(K))                                        00007500
2050   CONTINUE                                                         00007600
      J=1                                                               00007700
      Z(1)=EX(1)                                                        00007800
      DO 2150 K=1,3                                                     00007900
      T=OVZ(K)**(2.D0/(DFLOAT(L)+0.5D0))                                00008000
      T=(2.D0-T+DSQRT((2.D0-T)**2-T))/T                                 00008100
      DO 2100 I=1,39                                                    00008200
      J=J+1                                                             00008300
      Z(J)=Z(J-1)*T                                                     00008400
      IF(Z(J).GT.EX(K+1)) GO TO 2110                                    00008500
2100    CONTINUE                                                        00008600
2110   J=J-1                                                            00008700
2150    CONTINUE                                                        00008800
      NBAS=J                                                            00008900
      IJ=0                                                              00009000
      DO 40 I=1,NBAS                                                    00009100
      RNOR(I)=1.D0/DSQRT(PURX(LL,Z(I)+Z(I)))                            00009200
      DO 40 J=1,I                                                       00009300
      ZZ=Z(I)+Z(J)                                                      00009400
      T=0.D0                                                            00009500
      DO 35 K=1,NBPOT                                                   00009600
   35 T=T+PURX(LL+NPOT(K),ZZ+APOT(K))*CPOT(K)                           00009700
      T=T*RNOR(I)*RNOR(J)                                               00009800
      IJ=IJ+1                                                           00009900
      VKL(IJ)=T                                                         00010000
   40 CONTINUE                                                          00010100
      IJ=0                                                              00010200
      WRITE(6,8883)                                                     00010300
 8883 FORMAT(/,'  ELEMTS DE MATRICE DS LA BASE D''EXTRACTION')          00010400
      DO 45 I=1,NBAS                                                    00010500
      WRITE(6,42) (VKL(IJ+J),J=1,I)                                     00010600
  45  IJ=IJ+I                                                           00010700
C                                                                       00010800
C     VERIFICATION A PARTIR DE DONNEES SUR FILE 20 (NEWBAS=T ET UPUNCH=T00010900
C                                                                       00011000
      IF(.NOT.NEWBAS) GOTO 3900                                         00011100
      REWIND 20                                                         00011200
 102  READ(20,END=1750) ANAM,NPMAX                                      00011300
      IF(ANAM.EQ.APSD) GOTO 110                                         00011400
      DO 105 J=1,NPMAX                                                  00011500
 105  READ(20)                                                          00011600
      GOTO 102                                                          00011700
 1750 WRITE(6,9995) APSD                                                00011800
 9995 FORMAT(//,'   PAS DE PSEUDO-POTENTIEL SUR FT20 POUR =',A6,//)     00011900
      GOTO 1                                                            00012000
 110  IF(L.EQ.1) GOTO 120                                               00012100
      L1=L-1                                                            00012200
      DO 115 I=1,L1                                                     00012300
 115  READ(20)                                                          00012400
 120  READ(20) NOP,NVOIS,ZOP2,TEMP2,EPSDS2                              00012500
      IF (NOP.NE.0) GOTO 125                                            00012600
      WRITE(6,9996) APSD,L                                              00012700
      GOTO 1                                                            00012800
 9996 FORMAT(//,'   PAS DE PSEUDO-POTENTIEL SUR FT20 POUR=',A6,         00012900
     1'  DE SYMETRIE L=',I3)                                            00013000
 125  CONTINUE                                                          00013100
      WRITE(6,3678) APSD,NOP,NVOIS,(ZOP(I),I=1,NOP)                     00013200
      WRITE(6,3681)                                                     00013300
      WRITE(6,3679) (EPSDS(I),I=1,NVOIS)                                00013400
      WRITE(6,3682)                                                     00013500
      IJ=0                                                              00013600
      DO 3654 I=1,NVOIS                                                 00013700
      WRITE(6,3679) (TEMP(IJ+J),J=1,NOP)                                00013800
 3654 IJ=IJ+NOP                                                         00013900
 3678 FORMAT(////,' OP. NON LOCAL POUR ATOME',5X,A8,//,                 00014000
     *' NOMBRE DE PRIMITIVES',I3,/,                                     00014100
     *' NOMBRE DE CONTRACTEES',I3,/,                                    00014200
     *' EXPOSANTS',/,(10F12.6))                                         00014300
 3679 FORMAT(1X,8D15.8)                                                 00014400
 3681 FORMAT(' COEFFICIENTS DE L''OP NON LOCAL')                        00014500
 3682 FORMAT(' COEFFICIENTS DE CONTRACTION')                            00014600
      DO 3100 I=1,NOP                                                   00014700
 3100 OPNOR(I)=1.D0/DSQRT(PURX(LL,ZOP(I)+ZOP(I)))                       00014800
      IJ=0                                                              00014900
      DO 3200 I=1,NBAS                                                  00015000
      DO 3200 J=1,NOP                                                   00015100
      IJ=IJ+1                                                           00015200
 3200 REC(IJ)=PURX(LL,Z(I)+ZOP(J))*RNOR(I)*OPNOR(J)                     00015300
      CALL ORTHOE(REC,SV,TEMP,NOP,NVOIS,NBAS)                           00015400
      IJ=0                                                              00015500
      DO 3300 I=1,NBAS                                                  00015600
      NI=(I-1)*NVOIS                                                    00015700
      DO 3300 J=1,I                                                     00015800
      NJ=(J-1)*NVOIS                                                    00015900
      IJ=IJ+1                                                           00016000
      T=0.D0                                                            00016100
      DO 3350 K=1,NVOIS                                                 00016200
 3350 T=T+REC(NI+K)*REC(NJ+K)*EPSDS(K)                                  00016300
 3300 VERI(IJ)=T-VKL(IJ)                                                00016400
      WRITE(6,6)                                                        00016500
      WRITE(6,3400)                                                     00016600
 3400 FORMAT(' :::::::: VERIFICATION :::::::')                          00016700
      IJ=0                                                              00016800
      DO 3410 I=1,NBAS                                                  00016900
      WRITE(6,585) (VERI(IJ+J),J=1,I)                                   00017000
 3410 IJ=IJ+I                                                           00017100
      GOTO 1                                                            00017200
 3900 CONTINUE                                                          00017300
      IF(NEWBAS)GO TO 2160                                              00017400
      IF(OZ(1).EQ.0.D0) GO TO  2160                                     00017500
      IF(NOP.EQ.0) NOP=NBAS                                             00017600
      AOP=L/(2.D0*OZ(1)*OZ(1))                                          00017700
      BOP=L/(2.D0*OZ(2)*OZ(2))                                          00017800
      BOP=(BOP/AOP)**(1.D0/DFLOAT(NOP))                                 00017900
2160     CONTINUE                                                       00018000
      IF(NBAS.LT.40) GO TO 2250                                         00018100
      WRITE(6,2210)                                                     00018200
2210   FORMAT(' NOMBRE D EXP GENERES >40')                              00018300
      STOP                                                              00018400
2200   CONTINUE                                                         00018500
      IF(ARAF.EQ.0.) GO TO 11                                           00018600
      Z(1)=ARAF                                                         00018700
      DO 10 I=2,NBAS                                                    00018800
10    Z(I)=Z(I-1)*BRAF                                                  00018900
2250     CONTINUE                                                       00019000
11     IF(AOP.EQ.0.) GO TO 12                                           00019100
      IF(NEWBAS) GO TO 12                                               00019200
      ZOP(1)=AOP                                                        00019300
      DO 20 I=2,NOP                                                     00019400
20    ZOP(I)=ZOP(I-1)*BOP                                               00019500
      IF(NVOIS.EQ.0) NVOIS=NOP                                          00019600
12      CONTINUE                                                        00019700
      WRITE(6,6003)(Z(I),I=1,NBAS)                                      00019800
6003  FORMAT(/,' EXPOSANTS DE LA BASE D''EXTRACTION',/,(10F12.6))       00019900
      WRITE(6,6004)(ZOP(I),I=1,NOP)                                     00020000
6004  FORMAT(/,' EXPOSANTS DE LA BASE D''OPERATEURS',/,(10F12.6))       00020100
      I=INIPUR(1)                                                       00020200
      IF(NEWBAS) GO TO 721                                              00020300
      DO 21  I=1,NOP                                                    00020400
      OPNOR(I)=1.D0/DSQRT(PURX(LL,ZOP(I)+ZOP(I)))                       00020500
21     CONTINUE                                                         00020600
       IJ=0                                                             00020700
      DO 710 I=1,NOP                                                    00020800
      DO 710 J=1,I                                                      00020900
      IJ=IJ+1                                                           00021000
710    SL(IJ)=OPNOR(I)*OPNOR(J)*PURX(LL,  ZOP(I)+ZOP(J))                00021100
      WRITE(6,759)                                                      00021200
      WRITE(6,31) SL(2),SL(4),SL(7)                                     00021300
759     FORMAT('0 ORTHONORM DES  OPERATEURS')                           00021400
      CALL DIAGOS(SL,TOPE,TEMP,EPSDS,NOP,NOPR,DMNXP2)                   00021500
       DO 720 J=1,NOPR                                                  00021600
      IF(EPSDS(J).GT.OPTES) NVOISI=NVOISI+1                             00021700
      INUM(J)=J                                                         00021800
720   CONTINUE                                                          00021900
761   FORMAT(1X,10F12.7)                                                00022000
721       CONTINUE                                                      00022100
      IJ=0                                                              00022200
      JI=0                                                              00022300
      DO 30 I=1,NBAS                                                    00022400
      ZZ=Z(I)                                                           00022500
      RNOR(I)=1.D0/DSQRT(PURX(LL,ZZ+ZZ))                                00022600
      DO 25 J=1,I                                                       00022700
      T=RNOR(I)*RNOR(J)*PURX(LL,ZZ+Z(J))                                00022800
      JI=JI+1                                                           00022900
25    SL(JI)=T                                                          00023000
      DO 30 J=1,NOP                                                     00023100
      IJ=IJ+1                                                           00023200
30    REC(IJ)=PURX(LL,ZZ+ZOP(J))*RNOR(I)*OPNOR(J)                       00023300
      CALL DIAGOS(SL,TORT,TEMP,EPSDS,NBAS,NBASLD,DMNXP2)                00023400
      WRITE(6,31) SL(2),SL(4),SL(7)                                     00023500
31      FORMAT(' RECOUVREMENT ENTRE VOISINS:',3F15.8)                   00023600
      FNORMO=1.D10                                                      00023700
      NBAI=0                                                            00023800
      DO 730 I=1,NBASLD                                                 00023900
      IF(EPSDS(I).GT.BASTES) NBAI=NBAI+1                                00024000
730    CONTINUE                                                         00024100
        NBASLD=NBAI                                                     00024200
      IF(UPROJ       ) GO TO 620                                        00024300
      IF(UCUT.AND.NOPR.GT.NBASLD) NOPR=NBASLD                           00024400
      IF(NOPR.GT.MAXBAS) NOPR=MAXBAS                                    00024500
      GO TO 640                                                         00024600
620    IJ=0                                                             00024700
      IF(NEWBAS) GO TO 640                                              00024800
      DO 630 I=1,NOP                                                    00024900
      DO 630 J=1,I                                                      00025000
      ZZ=ZOP(I)+ZOP(J)                                                  00025100
      T=0.                                                              00025200
      DO 625 K=1,NBPOT                                                  00025300
625    T=T+PURX(LL+NPOT(K),ZZ+APOT(K))*CPOT(K)                          00025400
      T=T*OPNOR(I)*OPNOR(J)                                             00025500
      IJ=IJ+1                                                           00025600
      SC(IJ)=T                                                          00025700
630       CONTINUE                                                      00025800
      CALL ORTHO(SC,TEMP,SV,TOPE,NOP,NOPR)                              00025900
       IJ=0                                                             00026000
       DO 635 I=1,NOP                                                   00026100
       IJ=IJ+I                                                          00026200
635     SC(IJ)=0.5D0*SC(IJ)                                             00026300
      WRITE(6,6)                                                        00026400
      WRITE(6,*) (SC(I),I=1,IJ)                                         00026500
640    CONTINUE                                                         00026600
      IF(NVOISI.GT.NBAI)NVOISI=NBAI                                     00026700
      CALL ORTHO(VKL,TEMP,SV,TORT,NBAS,NBASLD)                          00026800
        IJ=NBAS*(NBAS+1)/2                                              00026900
      WRITE(6,6)                                                        00027000
6     FORMAT('0',20(' * '),/,'0')                                       00027100
      CALL ORTHOD    (REC,SV,TORT,NBAS,NBASLD,NOP)                      00027200
      CALL ORTHOE(REC,SV,TOPE,NOP,NOPR,NBASLD)                          00027300
42    FORMAT(' VKL=',12F10.4)                                           00027400
      NBASS=NBAS                                                        00027500
      NBAS=NBASLD                                                       00027600
      NOPS=NOP                                                          00027700
      NOP=NOPR                                                          00027800
        NOPEF=0                                                         00027900
      JI=0                                                              00028000
      NVOIS=NOP                                                         00028100
      IF(.NOT.UOLDPS) GO TO 680                                         00028200
      IJ=0                                                              00028300
      DO 670 I=1,NBAS                                                   00028400
      DO 670 J=1,I                                                      00028500
      IJ=IJ+1                                                           00028600
670   VKL(IJ)=VKL(IJ)-OLDPS(IJ)                                         00028700
680     CONTINUE                                                        00028800
      IF(.NOT. UVERI) GO TO 690                                         00028900
      IF(IVERI) GO TO 690                                               00029000
      IJ=0                                                              00029100
      DO 685 I=1,NBAS                                                   00029200
      DO 685 J=1,I                                                      00029300
      IJ=IJ+1                                                           00029400
685   VKL(IJ)=VKL(IJ)-VERI(IJ)                                          00029500
690   CONTINUE                                                          00029600
      WRITE(6,53) NVOIS                                                 00029700
53    FORMAT(' NB BANDES COMPATIBLES AVEC MDIM',I4)                     00029800
      IF(NEWBAS) GO TO 349                                              00029900
      IF(UPROJ) GO TO 349                                               00030000
      I2=NOP                                                            00030100
      IJJI=0                                                            00030200
      IJ=0                                                              00030300
      DO 200 II=1,NVOIS                                                 00030400
      I=INUM(II)                                                        00030500
      DO 190 IIV=1,II                                                   00030600
      IV=INUM(IIV)                                                      00030700
      IP=IV-I                                                           00030800
      T=REC(I)*REC(IV)                                                  00030900
      TEMP(1)=T+T                                                       00031000
      MN=1                                                              00031100
      TT=T*VKL(1)                                                       00031200
      MI=NOP+I                                                          00031300
      DO 60 M=2,NBAS                                                    00031400
      MN1=M-1                                                           00031500
       REI=REC(MI)                                                      00031600
      REV=REC(MI+IP)                                                    00031700
      NI=I                                                              00031800
      DO 55 N=1,MN1                                                     00031900
      MN=MN+1                                                           00032000
      T=REI*REC(NI+IP)+REV*REC(NI)                                      00032100
      TT=TT+T*VKL(MN)                                                   00032200
      NI=NI+NOP                                                         00032300
55    TEMP(MN)=T                                                        00032400
      MI=MI+NOP                                                         00032500
      MN=MN+1                                                           00032600
      T=REI*REV                                                         00032700
      TEMP(MN)=T+T                                                      00032800
60    TT=TT+T*VKL(MN)                                                   00032900
      IJ=IJ+1                                                           00033000
      SV(IJ)=TT                                                         00033100
      JI=0                                                              00033200
      DO 180 JJ=1,NVOIS                                                 00033300
      J=INUM(JJ)                                                        00033400
      DO 100 JJV=1,JJ                                                   00033500
      JV=INUM(JJV)                                                      00033600
      JP=JV-J                                                           00033700
      JI=JI+1                                                           00033800
      TT=TEMP(1)*REC(J)*REC(JV)                                         00033900
      MN=1                                                              00034000
      MJ=NOP+J                                                          00034100
      DO 80 M=2,NBAS                                                    00034200
      MN1=M-1                                                           00034300
      REJ=REC(MJ)                                                       00034400
      REV=REC(MJ+JP)                                                    00034500
      NJ=J                                                              00034600
      DO 70 N=1,MN1                                                     00034700
      MN=MN+1                                                           00034800
      TT=TT+TEMP(MN)*(REJ*REC(NJ+JP)+REV*REC(NJ))                       00034900
70    NJ=NJ+NOP                                                         00035000
      MJ=MJ+NOP                                                         00035100
      MN=MN+1                                                           00035200
      TT=TT+TEMP(MN)*REJ*REV                                            00035300
80    CONTINUE                                                          00035400
      IJJI=IJJI+1                                                       00035500
      VL(IJJI)=TT                                                       00035600
      IF(JI.NE.IJ) GO TO 100                                            00035700
      TT=1.D0/DSQRT(TT)                                                 00035800
      SV(IJ)=SV(IJ)*TT                                                  00035900
      SN(IJ)=TT                                                         00036000
      N=IJJI-IJ                                                         00036100
      DO 90 M=1,IJ                                                      00036200
90    VL(N+M)=TT*SN(M)*VL(N+M)                                          00036300
      GO TO 185                                                         00036400
100   CONTINUE                                                          00036500
180   CONTINUE                                                          00036600
185   CONTINUE                                                          00036700
190   CONTINUE                                                          00036800
200   I2=I2-1                                                           00036900
      CALL DIAGOS(VL,VL,VL,TEMP,IJ,IJP,DMNXP)                           00037000
      WRITE(6,206) IJ,IJP                                               00037100
206   FORMAT(' NB TERME',I4,' NB REDUIT',I4)                            00037200
      DO 300 I=1,IJ                                                     00037300
      TT=0.                                                             00037400
      DO 260 J=1,IJ                                                     00037500
      T=0.                                                              00037600
      JM=J                                                              00037700
      IM=I                                                              00037800
      DO 250 M=1,IJP                                                    00037900
      T=T+VL(IM)*VL(JM)                                                 00038000
      JM=JM+IJ                                                          00038100
250   IM=IM+IJ                                                          00038200
260   TT=TT+T*SV(J)                                                     00038300
300   SC(I)=TT                                                          00038400
        DO 850 I=1,IJ                                                   00038500
 850    SC(I)=SC(I)*SN(I)                                               00038600
      WRITE(6,*) (SC(I),I=1,IJ)                                         00038700
349     CONTINUE                                                        00038800
      IJ=NVOIS*(NVOIS+1)/2                                              00038900
348     CONTINUE                                                        00039000
      MN=0                                                              00039100
      DO 350 I=1,NBAS                                                   00039200
      DO 350 J=1,I                                                      00039300
      MN=MN+1                                                           00039400
350   TEMP(MN)=0.                                                       00039500
      IJ=0                                                              00039600
      IP=0                                                              00039700
      I2=NOP                                                            00039800
      I1=1                                                              00039900
      DO 500 II=1,NVOIS                                                 00040000
      I=INUM(II)                                                        00040100
      DO 490 IIV=1,II                                                   00040200
      IV=INUM(IIV)                                                      00040300
      IP=IV-I                                                           00040400
        IJ=IJ+1                                                         00040500
      T=REC(I)*REC(IV)                                                  00040600
      TEMP(1)=(T+T)*SC(IJ)+TEMP(1)                                      00040700
      MN=1                                                              00040800
      MI=NOP+I                                                          00040900
      DO460 M=2,NBAS                                                    00041000
      MN1=M-1                                                           00041100
       REI=REC(MI)                                                      00041200
      REV=REC(MI+IP)                                                    00041300
      NI=I                                                              00041400
      DO455 N=1,MN1                                                     00041500
      MN=MN+1                                                           00041600
      T=REI*REC(NI+IP)+REV*REC(NI)                                      00041700
      NI=NI+NOP                                                         00041800
      TEMP(MN)=T   *SC(IJ)+TEMP(MN)                                     00041900
455    CONTINUE                                                         00042000
      MN=MN+1                                                           00042100
      T=REI*REV                                                         00042200
      TEMP(MN)=(T+T)*SC(IJ)+TEMP(MN)                                    00042300
      MI=MI+NOP                                                         00042400
460   CONTINUE                                                          00042500
490       CONTINUE                                                      00042600
      IP=IP+1                                                           00042700
      I2=I2-1                                                           00042800
500        CONTINUE                                                     00042900
      WRITE(6,6)                                                        00043000
      WRITE(6,6)                                                        00043100
      IJ=0                                                              00043200
      DO 610 I=1,NBAS                                                   00043300
      DO 610 J=1,I                                                      00043400
      IJ=IJ+1                                                           00043500
      OLDPS(IJ)=TEMP(IJ)                                                00043600
610   CONTINUE                                                          00043700
      IF(.NOT.UVERI) GO TO 615                                          00043800
      IVERI=.NOT.IVERI                                                  00043900
      IJ=0                                                              00044000
      DO 614 I=1,NBAS                                                   00044100
      DO 614 J=1,I                                                      00044200
      IJ=IJ+1                                                           00044300
614    VERI(IJ)=TEMP(IJ)                                                00044400
615   CONTINUE                                                          00044500
        IJ=NBAS*(NBAS+1)/2                                              00044600
      T=TEMP(1)-VKL(1)                                                  00044700
      TEMP(1)=T                                                         00044800
      MN=1                                                              00044900
      FNORM=T*T                                                         00045000
      DO 520 M=2,NBAS                                                   00045100
      MN1=M-1                                                           00045200
      DO 515 N=1,MN1                                                    00045300
      MN=MN+1                                                           00045400
      T=TEMP(MN)-VKL(MN)                                                00045500
      TEMP(MN)=T                                                        00045600
515     FNORM=FNORM+2.D0*T*T                                            00045700
      MN=MN+1                                                           00045800
      T=TEMP(MN)-VKL(MN)                                                00045900
      TEMP(MN)=T                                                        00046000
520    FNORM=FNORM+T*T                                                  00046100
510     FORMAT(' NORME',F20.8)                                          00046200
         WRITE(6,510) FNORM                                             00046300
      IJ=0                                                              00046400
      DO 550 J=1,NBAS                                                   00046500
      JI=J                                                              00046600
      T=EPSDS(J)                                                        00046700
      DO 550 I=1,NBASS                                                  00046800
      IJ=IJ+1                                                           00046900
      VL(JI)=TORT(IJ)*T                                                 00047000
550   JI=JI+NBAS                                                        00047100
      CALL ORTHO(TEMP,SL,SV,VL,NBAS,NBASS)                              00047200
      IJ=0                                                              00047300
      WRITE(6,584)                                                      00047400
584   FORMAT('0  DIFFERENCE DES ELMTS DE MAT DU PSEUDO DS LA BASE INI') 00047500
      DO 580 I=1,NBASS                                                  00047600
      WRITE(6,585) (TEMP(IJ+J),J=1,I)                                   00047700
580   IJ=IJ+I                                                           00047800
      NOP=NOPS                                                          00047900
       IF(NEWBAS)GO TO 1                                                00048000
 997    NTER=NVOIS*(NVOIS+1)/2                                          00048100
585   FORMAT(1X,12D10.2)                                                00048200
       FNORM=FNORM/NBAS                                                 00048300
      IF(FNORM.GT.1.D0) GO TO 1                                         00048400
         IJV=0                                                          00048500
         DO 1950I=1,NVOIS                                               00048600
         DO1950J=1,I                                                    00048700
          IJV=IJV+1                                                     00048800
          T=SC(IJV)                                                     00048900
          IF(I.EQ.J)T=T+T                                               00049000
         TEMP(NVOIS*(I-1)+J)=T                                          00049100
1950     TEMP(NVOIS*(J-1)+I)=T                                          00049200
         CALL JACVAO(TEMP,NVOIS,EPSDS,TORT)                             00049300
        DO2000I=1,NVOIS                                                 00049400
        DO2000J=1,NOP                                                   00049500
         JK=J                                                           00049600
         KI=NVOIS*(I-1)                                                 00049700
         T=0.D0                                                         00049800
         DO1900K=1,NVOIS                                                00049900
         KI=KI+1                                                        00050000
         T=T+TORT(KI)*TOPE(JK)                                          00050100
1900      JK=JK+NOP                                                     00050200
2000       TEMP(NOP*(I-1)+J)=T                                          00050300
      IF(UPUNCH) GO TO 5647                                             00050400
C......STOCKAGE SUR FILE20                                              00050500
C......RELECTURE DE FILE20                                              00050600
C......TEST POUR SAVOIR SI CE CAS EXISTE DEJA                           00050700
 699  CONTINUE                                                          00050800
      REWIND 20                                                         00050900
      REWIND 21                                                         00051000
C 700  READ(20,END=750,ERR=787) ANAM,NPMAX                               00051100
 700  READ(20,END=750,ERR=750) ANAM,NPMAX                               00051100
      IF (ANAM.EQ.APSD) GOTO 705                                        00051200
      WRITE(21) ANAM,NPMAX                                              00051300
      DO 702 I=1,NPMAX                                                  00051400
      READ(20,END=780,ERR=780)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1           00051500
 702  WRITE(21)NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                            00051600
      GOTO 700                                                          00051700
 705  IF( L.LE.NPMAX)      GOTO 1720                                    00051800
      NP1=NPMAX+1                                                       00051900
C                                                                       00052000
C.....PSEUDO EXISTE DEJA. ON AJOUTE LA SYMETRIE SUIVANTE                00052100
C                                                                       00052200
      WRITE(21) ANAM,NP1                                                00052300
      DO 1710 I=1,NPMAX                                                 00052400
      READ(20,END=790,ERR=790)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1           00052500
 1710  WRITE(21) NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                          00052600
C                                                                       00052700
C.....SYMETRIE ACTUELLEMENT CALCULEE                                    00052800
C                                                                       00052900
      WRITE(21) NOP,NVOIS,ZOP2,TEMP2,EPSDS2                             00053000
C                                                                       00053100
C.....ON COMPLETE LA FILE21                                             00053200
C                                                                       00053300
 715  READ(20,END=760,ERR=789) ANAM,NPMAX                               00053400
         WRITE(6,783) ANAM,NPMAX                                        00053500
 783     FORMAT(' *****715****',A6,I6)                                  00053600
      WRITE(21) ANAM ,NPMAX                                             00053700
      DO 718 I=1,NPMAX                                                  00053800
      READ(20,END=795,ERR=795)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1           00053900
 718  WRITE(21) NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                           00054000
      GOTO 715                                                          00054100
 1720  IF(.NOT.UREPL) GOTO 800                                          00054200
C                                                                       00054300
C.....PSEUDO EXISTE DEJA. ON REMPLACE UNE SYMETRIE DEJA EXISTANTE       00054400
C                                                                       00054500
      WRITE(21) ANAM,NPMAX                                              00054600
      IF (L.EQ.1) GOTO  726                                             00054700
      L1=L-1                                                            00054800
      DO 725 I=1,L1                                                     00054900
      READ(20,END=796,ERR=796)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1           00055000
 725  WRITE(21) NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                           00055100
 726  CONTINUE                                                          00055200
      READ(20,END=797,ERR=797)                                          00055300
      WRITE(21) NOP ,NVOIS ,ZOP2,TEMP2,EPSDS2                           00055400
      IF(L.EQ.NPMAX)      GOTO 715                                      00055500
      L1=L+1                                                            00055600
      DO 1730 I=L1,NPMAX                                                00055700
      READ(20,END=798,ERR=798)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1           00055800
 1730  WRITE(21) NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                          00055900
C                                                                       00056000
C.....ON COMPLETE LA FILE21                                             00056100
C                                                                       00056200
      GOTO 715                                                          00056300
 750  CONTINUE                                                          00056400
C                                                                       00056500
C.....LE PSEUDO N'EXISTE PAS SUR LA FILE                                00056600
C.....ON LE RAJOUTE A LA FIN                                            00056700
C                                                                       00056800
      NPMAX=1                                                           00056900
      WRITE(21) APSD,NPMAX                                              00057000
      WRITE(21) NOP,NVOIS,ZOP2,TEMP2,EPSDS2                             00057100
C                                                                       00057200
C.....REECRITURE DE LA FILE20 MISE A JOUR                               00057300
C                                                                       00057400
 760  CONTINUE                                                          00057500
      REWIND 20                                                         00057600
      REWIND 21                                                         00057700
 765  READ(21,END=1) ANAM,NPMAX                                         00057800
      WRITE(20) ANAM,NPMAX                                              00057900
      WRITE(6,9843) ANAM,NPMAX                                          00058000
      DO 770 I=1,NPMAX                                                  00058100
      READ(21)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                           00058200
      WRITE(20) NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                           00058300
      IF(ULISTE) WRITE(6,9844) I,NOP1,NVOIS1,ZOP1                       00058400
 770  CONTINUE                                                          00058500
 9843 FORMAT(/'    ATOME=',A8,'   NB DE SYMETRIE ',I3)                  00058600
 9844 FORMAT('   SYMETRIE=',I3,'  NB D OPERATEURS=',I3,                 00058700
     1'  NB DE COMB.ORTHOG.=',I3 ,/,' EXPOSANTS=',12F10.6,(/,13F10.6))  00058800
      GOTO 765                                                          00058900
C                                                                       00059000
C ERREUR LECTURE                                                        00059100
C                                                                       00059200
 787    ICAS=700                                                        00059300
        WRITE(6,781) ANAM,NPMAX,ICAS                                    00059400
        STOP                                                            00059500
 789     ICAS=715                                                       00059600
        WRITE(6,781) ANAM,NPMAX,ICAS                                    00059700
        STOP                                                            00059800
 780    ICAS=702                                                        00059900
        WRITE(6,781) ANAM,NPMAX,ICAS                                    00060000
        STOP                                                            00060100
 790    ICAS=1710                                                       00060200
        WRITE(6,781) ANAM,NPMAX,ICAS                                    00060300
        STOP                                                            00060400
 795    ICAS=718                                                        00060500
        WRITE(6,781) ANAM,NPMAX,ICAS                                    00060600
        STOP                                                            00060700
 796     ICAS=725                                                       00060800
        WRITE(6,781) ANAM,NPMAX,ICAS                                    00060900
        STOP                                                            00061000
 797      ICAS=0                                                        00061100
        WRITE(6,781) ANAM,NPMAX,ICAS                                    00061200
        STOP                                                            00061300
 798      ICAS=1730                                                     00061400
        WRITE(6,781) ANAM,NPMAX,ICAS                                    00061500
        STOP                                                            00061600
 781   FORMAT(//' ERREUR SUR FILE20 APRES ',A6,I3,'  ICAS=',I5)         00061700
C                                                                       00061800
C..... ERREUR....ON ESSAIE DE REMPLACER UN PSEUDO EXISTANT              00061900
C      SANS SPECIFIER  UREPL=.TRUE.                                     00062000
C                                                                       00062100
 800  CONTINUE                                                          00062200
       WRITE(6,9845)                                                    00062300
 9845  FORMAT(//,'*****ESSAI DE REMPLACEMENT D UN PSEUDO SANS SPECIFIER 00062400
     1UREPL=.TRUE.********')                                            00062500
      STOP                                                              00062600
 5647 CONTINUE                                                          00062700
      WRITE(8,FMT) IDEN,L,NOP,NVOIS,NVOIS                               00062800
      WRITE(8,FMT) (ZOP(I),I=1,NOP)                                     00062900
      WRITE(8,FMT)((TEMP( I+(J-1)*NOP),I=1,NOP),J=1,NVOIS)              00063000
                        WRITE(8,FMT)(EPSDS(I),I=1,NVOIS    )            00063100
 5650 CONTINUE                                                          00063200
999    CONTINUE                                                         00063300
      STOP                                                              00063400
      END                                                               00063500
$NOSTANDARD SYSTEM
      SUBROUTINE GIVENS (h,E,V,N,NEV,NVEC)
      IMPLICIT REAL*8 (A-H,O-X,Z)
      PARAMETER(NDETZ=5000,METZ=37)
C
C                       NSTART(I) : POSITION DU PREMIER ELEMENT DE
C                       LIGNE I
C
C                  H    matrice a diagonaliser rangee sous la forme d'un
c                       demi matrice superieure
c                  E    VALEURS PROPRES EN ORDRE CROISSANT
C                       RANGEES DANS UN VECTEUR DE DIMENSION  N
C
C                  V    MATRICES DES VECTEURS PROPRES
C                        DE DIMENSION NVEC
C
C                  N    DIMENSION DE LA MATRICE
C                  NEV   NOMBRE DE VALAURS PROPRES DEMANDEES
C                  NVEC  NOMBRE DE VECTEURS  PROPRES DEMANDES
C
      INTEGER AG,RN,ALF
      REAL*8 NORM,LAMBD,L,MULT
      DIMENSION h(1),E(1),V(1),B(NDETZ),C(NDETZ),P(NDETZ),Q(NDETZ),
     1R(NDETZ),W(NDETZ),VY(NDETZ+2),IN(NDETZ),Z(NDETZ)
      DIMENSION HS1(NDETZ),HS2(NDETZ)
c***********************************************************************
c     nstart(i)=indice du premier element de la ligne i d'une demi-matri
c     superieure
      nstart(i)=n*(n+1)/2-(n-i+1)*(n-i+2)/2+1
c***********************************************************************
C  SOIT I  LE NOMBRE DE DECIMALES POUR LA MACHINE UTILISEE. ON PREND
C  EPSI=10**(-(I+3))   ET   EPSDG=EPSV= 10**(-(I-2))
      I=15
      EPSI=10.0**(-(I+3))
      EPSV=10.0**(-(I-2))
      EPSDG=EPSV
C     EPSI=1.D-10
C     EPSDG=1.0D-5
C     EPSV=1.D-5
      ITEVM=3
      JOYM=15
C       MODM1 = 2**15 - 1
      MODM1 = 32767
c dlgmax= puissance de 10 maximale
      dlgmax=36.d0
      RN=MODM1 - 2
      ALF=259
      AGNES=2.0**(-14)
      NP1=N+1
      NM1=N-1
      NM2=N-2
      NG=1
      IF(NM2)99,9,200
  200 CONTINUE
c transfert de la premiere ligne de h dans hs1
      DO 2 I=1,N
    2 HS1(I)=h(I)
C
C     DEBUT DE LA TRIDIAGONALISATION
C
C     HS1 CONTIENT LA LIGNE IR  A TRIDIAGONALISER
C
      DO 801 IR=1,NM2
      C(IR)=HS1(IR)
      IP1=IR+1
      IP2=IR+2
C
C     DETERMINATION DES ROTATIONS ANNULANT LES ELEMENTS IR+2,..,N
C     LE SINUS DE L'ANGLE DE ROTATION EST STOCKE A LA PLACE DE L'ELEMENT
C     ANNULE
      DO 802 I=IP2,N
      T=HS1(I)
      IF(dabs(t).lt.1.d-8) then
      HS1(I)=1.
      HS2(I-IP1)=0.
      else
      S=HS1(IP1)
      SS=DSQRT(S*S+T*T)
      HS1(IP1)=SS
      S=S/SS
      T=T/SS
      HS1(I)=S
      HS2(I-IP1)=T
      endif
  802 CONTINUE
C
C
C     ISH3 : INDICE-IP1 DU PREMIER ELEMENT DE LA LIGNE IP1
C
      ISH3=NSTART(IP1)-IP1
C
      DO 805 K=IP2,N
      KP1=K+1
C
C     ISH4 :INDICE-K DU PREMIER ELEMENT DE LA LIGNE K
      ISH4=NSTART(K)-K
C
C
C     ON FAIT LES ROTATIONS ENTRE LES LIGNES IR+1 ET K K=IR+2,...,N
C
      S=HS1(K)
      T=HS2(K-IP1)
      IF(T) 820,821,820
  820 CONTINUE
      IH4K=ISH4+K
      IH3K=ISH3+K
      IH3P=ISH3+IP1
      U=H(IH3P)
      L=H(IH3K)
      TT=h(IH4K)
C
      h(IH3P)= (U*S+L*T)*S+(L*S+TT*T)*T
      h(IH3K) =-(U*S+L*T)*T+(L*S+TT*T)*S
      h(IH4K) = (U*T-L*S)*T+ (-L*T+TT*S)*S
C
C
      IF(KP1-N) 840,840,841
  840 CONTINUE
C
      DO 807 LL=KP1,N
      IH3L=ISH3+LL
      IH4L=ISH4+LL
      U=h(IH3L)
      L=h(IH4L)
      h(IH3L)= U*S+L*T
      h(IH4L)=-U*T+L*S
  807 CONTINUE
C
  821 CONTINUE
C
      IF(KP1-N) 823,823,841
  823 CONTINUE
      IH3K=ISH3+K
      DO 808 LL=KP1,N
      S=HS1(LL)
      T=HS2(LL-IP1)
      IF(T) 822,808,822
  822 CONTINUE
      IH4L=ISH4+LL
      U=h(IH3K)
      L=h(IH4L)
      h(IH3K)= U*S+L*T
      h(IH4L)= -U*T+L*S
  808 CONTINUE
  841 CONTINUE
C
C     FIN DE ROTATION DE LA LIGNE K
C
  805 CONTINUE
      B(IR)=HS1(IP1)
C
C     ON STOCKE DANS HS2 LE COSINUS DE LA ROTATION OU LE COSINUS PLUS
C     3. (SI LE SINUS EST NEGATIF)
      DO 3410 K=IP2,N
      S=HS1(K)
      IF(S.GE.0.D0) GO TO 3410
      HS2(K-IP1)=HS2(K-IP1)+3.0D00
 3410 CONTINUE
C
C     ON TRANSFERE LA LIGNE SUIVANTE DANS HS1
      DO 3050 K=IP1,N
 3050 HS1(K)=h(ISH3+K)
C   LES COSINUS DES ROTATIONS STOCKES DANS HS2 SONT TRANSFERES DANS h
      ISH1=NSTART(IR)-IR
      DO 3020 K=IP2,N
 3020 h(ISH1+K)=HS2(K-IP1)
  801 CONTINUE
C
      C(NM1)=h(NSTART(NM1))
      B(NM1)=h(NSTART(NM1)+1)
      C(N)=h(NSTART(N))
      B(N)=0.D0
      GO TO 850
    9 continue
      iiii=0
      C(NM1)=h(iiii+1)
      B(NM1)=h(iiii+2)
      C(N)=h(iiii+3)
      B(N)=0.D0
C
  850 CONTINUE
      NORM=DABS(C(1))+DABS(B(1))
      DO 10 I=2,N
      T=DABS(C(I))+DABS(B(I))+DABS(B(I-1))
      IF(NORM-T) 104,10,10
  104 NORM=T
   10 CONTINUE
      DO 11 I=1,N
   11 W(I)=B(I)*B(I)
      K=1
      L=-NORM
      MULT=NORM*EPSI
      DO 12 I=1,NEV
   12 E(I)=NORM
   13 U=E(K)
   14 LAMBD=0.5*(L+U)
      IF(LAMBD-L-MULT) 30,30,106
  106 IF(U-LAMBD-MULT) 30,30,107
  107 AG=0
      I=1
   16 S=C(I)-LAMBD
   18 IF(DABS(S)-1.D-21) 20,108,108
  108 IF(S) 110,110,109
  109 AG=AG+1
  110 I=I+1
      IF(I-N) 111,111,22
  111 S=C(I)-LAMBD-W(I-1)/S
      GO TO 18
   20 AG=AG+1
      I=I+2
      IF(I-N) 16,16,22
   22 IF(AG-N+K) 24,24,112
  112 L=LAMBD
      GO TO 14
   24 U=LAMBD
      IF(N-AG-NEV) 113,113,114
  113 M=N-AG
      GO TO 115
  114 M=NEV
  115 DO 26 I=K,M
   26 E(I)=LAMBD
      GO TO 14
   30 E(K)=LAMBD
      K=K+1
      IF(K-NEV) 13,13,116
  116 IF(NVEC) 40,999,40
   40 II=0
      INN=-N
      EPSIN=NORM*EPSDG
      DO 90 I=1,NVEC
      INN=INN+N
      IEPS=0
      VY(N+1)=0.
      VY(N+2)=0.
  402 T=E(I)
      IEPS=IEPS+1
      GO TO (430,404),IEPS
  404 J=I-II
      IF(J) 406,406,407
  406 S=NORM
      GO TO 408
  407 S=T-E(J)
  408 L=T
      DO 414 K=I,NEV
      U=E(K)
      IF(U-L - EPSIN) 410,416,416
  410 L=U
  414 CONTINUE
  416 TT=U-T
      IF(TT-S) 418,420,420
  418 S=TT
  420 U=DABS(T)
      IF(U-S) 422,424,424
  422 U=S
  424 T=T + U*EPSDG*0.1
      IR=0
  430 CONTINUE
      DO 44 J=1,N
      P(J)=0.
      Q(J)=B(J)
      R(J)=C(J)-T
      GO TO (42,44),IEPS
   42 VY(J)=1.0
   44 CONTINUE
      DO 50 J=1,NM1
      IF(DABS(R(J))+DABS(B(J))) 152,152,154
  152 IN(J)=0
      W(J)=0.
      R(J)=1.D-30
      GO TO 50
  154 CONTINUE
      IF(DABS(R(J))-DABS(B(J)))49,117,117
  117 MULT=B(J)/R(J)
      IN(J)=0
      GO TO 48
   49 MULT=R(J)/B(J)
      IN(J)=1
      R(J)=B(J)
      T=R(J+1)
      R(J+1)=Q(J)
      Q(J)=T
      P(J)=Q(J+1)
      Q(J+1)=0.
   48 W(J)=MULT
      Q(J+1)=Q(J+1)-MULT*P(J)
      R(J+1)=R(J+1)-MULT*Q(J)
      IF(DABS(R(J))-1.D-30) 118,118,50
  118 R(J)=1.D-30
   50 CONTINUE
      IF(DABS(R(N))-1.D-30)  119,119,120
  119 R(N)=1.D-30
  120 GO TO (1202,145),IEPS
 1202 IF(I-NVEC) 155,121,121
  155 IF(DABS(E(I+1)-E(I)) - EPSIN) 53,121,121
  121 IF(II) 122,55,122
  122 KK=1
      GO TO 54
   53 KK=2
   54 II=II+1
      IF(II-1) 123,55,123
  123 JOY=0
      IR=0
   51 JOY=JOY+1
      IF(JOY-JOYM) 124,95,95
  124 DO 52 J=1,N
      VY(J)=Ran(alf)
   52 CONTINUE
   55 ITEV=0
   56 ITEV=ITEV+1
      DO 66 JI=1,N
      K=N-JI+1
   62 T=VY(K)
      tnum= T-VY(K+1)*Q(K)-VY(K+2)*P(K)
      tden=r(k)
      if(dabs(tnum).gt.1.d-30) then
      dltnum=dlog(dabs(tnum))
      else
      dltnum=-30.d0
      end if
      if((dltnum-dlog(dabs(tden))).gt.dlgmax) then
      DO 64 J=1,N
   64 VY(J)=VY(J)*1.D-5
      GO TO 62
      end if
      vy(k)=tnum/tden
   66 CONTINUE
      IF(ITEV-1) 145,145,131
  131 S=E(I)
      TT=DABS((C(1)-S)*VY(1) + B(1)*VY(2))
      T=DABS(VY(1))
      DO 135 J=2,N
      U =DABS(B(J-1)*VY(J-1) + (C(J)-S)*VY(J) + B(J)*VY(J+1))
      IF(TT-U ) 132,133,133
  132 TT=U
  133 U =DABS(VY(J))
      IF(T-U ) 134,135,135
  134 T=U
  135 CONTINUE
      S=TT/T
      IF(S-EPSV) 136,136,139
  136 IF(ITEV-2) 69,69,137
  137 WRITE(6,100) I,ITEV
      GO TO 69
  139 WRITE(6,101) I,EPSV,ITEV ,TT,T,S
      IF(ITEV-ITEVM) 145,69,69
  145 DO 68 J=1,NM1
      IF(IN(J)) 144,144,67
  144 VY(J+1)=VY(J+1)-W(J)*VY(J)
      GO TO 68
   67 T=VY(J)
      VY(J)=VY(J+1)
      VY(J+1)=T-W(J)*VY(J+1)
   68 CONTINUE
      GO TO 56
   69 CONTINUE
      T=1./T
      DO 98 J=1,N
   98 VY(J)=VY(J)*T
      IF(II-1) 77,147,72
   72 JI=I-II+1
      M=I-1
      T=0.
      DO 70 J=1,N
      U=VY(J)
   70 T=T + U*U
      T=1./T
      AG=1
      IKK=N*(JI-2)
      DO 75 K=JI,M
      IKK=IKK+N
      S=0.
      DO 73 J=1,N
      JVK=IKK+J
   73 S=S+VY(J)*V(JVK)
      TT=S*Z(K)
      DO 74 J=1,N
      JVK=IKK+J
   74 VY(J) = VY(J)-TT*V(JVK)
      U=S*TT*T
      T=T/(1.-U)
      IF(U-0.75) 75,75,76
   76 AG=2
   75 CONTINUE
      GO TO (160,143),AG
  143 ITEV=1
      IR=IR+1
      GO TO (145,1432),IR
 1432 GO TO (402,51),IEPS
  160 IF(JOY-2) 147,147,146
  146 WRITE(6,100) I,II,JOY,IR
  147 S=0.
      DO 162 J=1,N
      JVI=INN+J
      U=VY(J)
      S=S + U*U
  162 V(JVI)=U
      Z(I)=1./S
C
C
      GO TO (78,90),KK
   77 DO 777 J=1,N
      JVI=INN+J
  777 V(JVI)=VY(J)
   78 CONTINUE
      II=0
   90 CONTINUE
      NBL1=0
      IF(NM2)85,85,81
   81 DO 84 J=1,NM2
      K=N-J-1
      KS=J+2
      KP=K+1
      M=K+2
      IR=N+1
      IS=J+1
      ISH1=NSTART(K)-K
C
      DO 83 KK=M,N
      IS=IS-1
      IR=IR-1
      T=h(ISH1+IR)
      IF(T.GT.1.5D0) GO TO 3460
      IF(T.EQ.0.D0) GO TO 83
      S=DSQRT(1.D0-T*T)
      GO TO 3470
 3460 T=T-3.0D0
      S=-DSQRT(1.0D0-T*T)
 3470 CONTINUE
      INN=-N
      DO 89 JOY=1,NVEC
      INN=INN+N
      KPJOY=INN+KP
      IRJOY=INN+IR
      U=V(KP JOY)
      L=V(IR JOY)
      V(KP JOY)=U*S-L*T
      V(IR JOY)=U*T+L*S
   89 CONTINUE
   83 CONTINUE
   84 CONTINUE
   85 CONTINUE
      INN=-N
      DO 87 JOY=1,NVEC
      INN=INN+N
      S=0.
      DO 86 J=1,N
      JJOY=INN+J
   86 S=S+V(J JOY)*V(J JOY)
      S=1./DSQRT(S)
      DO 88 J=1,N
      JJOY=INN+J
   88 V(J JOY)=V(J JOY)*S
   87 CONTINUE
  999 RETURN
   95 WRITE(6,102)I,II,JOYM
      CALL EXIT
  100 FORMAT(1H0,10I5)
  101 FORMAT(1H0,I4,E16.8,I3,3E18.8)
  102 FORMAT(1H0,I4,I3)
      RETURN
 99   continue
      E(1)=h(1)
      V(1)=1.0
      RETURN
      END
      FUNCTION PURX(M,ZZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FACTO(20),FATT(20)
      DATA PI/3.14159265358979/
      C=1.D0/ZZ
 20   M2=M/2
      IF((M2*2).NE.M) GOTO 60
      IF(M.EQ.0) GOTO 50
      PURX=DSQRT(PI*C)*FACTO(M)*.5*(.5*C)**M2
      GOTO 99
 50   PURX=.5D0*DSQRT(PI*C)
      GOTO 99
 60   M2=M2+1
      PURX=FATT(M2)*.5D0*C**M2
 99   RETURN
      ENTRY INIPUR(J)
      FACTO(1)=1.
      FACTO(2)=1.
      FATT(1)=1.
      FATT(2)=1.
      DO 200 I=2,19
      FACTO(I+1)=I*FACTO(I-1)
 200  FATT(I+1)=I*FATT(I)
      INIPUR=I
      RETURN
      END
C--------------------------------------------------------------------------     
      SUBROUTINE DIAGOS(S,DT,V,EPS,NBAS1,NBASLD,DMNXP)
      IMPLICIT REAL*8 (A-H,O-T,V-Z),LOGICAL  (U)
      COMMON/GIV/NBLOCK(400),NSTART(400),NSIZE,IMNM(400)
      parameter (nbmax=400)
      DIMENSION EPSI(1000)
      DIMENSION TAB(500)
      EQUIVALENCE ( TAB(1),EPSI(501) )
      DIMENSION DT(NBAS1,NBAS1),S(1),V(NBAS1,NBAS1),EPS(1)
      dimension h(nbmax*(nbmax+1)/2)
      ij=0
      do i=1,nbas1
      do j=i,nbas1
      ij=ij+1
      h(ij)=s(j*(j-1)/2+i)
      end do
      end do
      CALL GIVENS(h,EPSI,V,NBAS1,NBAS1,NBAS1)
      WRITE(6,63) ( EPSI(I),I=1,NBAS1 )
 63   FORMAT('  VP DE S ',12F10.6)
      K=0
      DO 50 II=1,NBAS1
      I=II
      IF(EPSI(I).LT.DMNXP) GOTO 50
      TERM=1.D0/DSQRT(EPSI(I))
      K=K+1
      EPS(NBAS1-K+1)=EPSI(I)
      DO 40 J=1,NBAS1
 40   DT(J,K)=V(J,I)*TERM
 50   CONTINUE
      NBASLD=K
      NBI=K/2
      DO 100 I=1,NBI
      II=K-I+1
      DO 80 J=1,NBAS1
 80   TAB(J)=DT(J,I)
      DO 90 J=1,NBAS1
      DT(J,I)=DT(J,II)
 90   DT(J,II)=TAB(J)
 100  CONTINUE
      WRITE(6,60) NBASLD
 60   FORMAT('  BASE REDUITE DIMENSION:',I4)
      RETURN
      END
      SUBROUTINE ORTHOE(V,EPSI,W,NOP,NR,NB)      
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(1),EPSI(1),W(NOP,NR)
      JK=0
      JL=0
      DO 80 I=1,NB
      DO 70 J=1,NOP
      JK=JK+1
 70   EPSI(J)=V(JK)
      DO 80 J=1,NR
      T=0.
      DO 60 K=1,NOP
 60   T=T+W(K,J)*EPSI(K)
      JL=JL+1
 80   V(JL)=T
      RETURN
      END
      SUBROUTINE ORTHOD(V,EPSI,W,NB,NR,NOP)
      IMPLICIT REAL*8 (A-H,O-T,V-Z),LOGICAL  (U)
      DIMENSION V(NOP,NB),W(NB,NB),EPSI(NB)
      DO 50 I=1,NOP
      DO 40 J=1,NB
 40   EPSI(J)=V(I,J)
      DO 50 J=1,NR
      T=0.
      DO 30 K=1,NB
 30   T=T+W(K,J)*EPSI(K)
 50   V(I,J)=T
      RETURN
      END
      SUBROUTINE ORTHO(VPHI,V,EPSI,W,NB,NR)
      IMPLICIT REAL*8 (A-H,O-T,V-Z),LOGICAL  (U)
      DIMENSION W(NB,NB),V(NB,NB),VPHI(1),EPSI(1)
      IJ=0
      DO 20 I=1,NB
      DO 20 J=1,I
      IJ=IJ+1
      V(I,J)=VPHI(IJ)
 20   V(J,I)=VPHI(IJ)
      KL=0
      DO 50 K=1,NR
      DO 35 J=1,NB
      T=0.
      DO 30 I=1,NB
 30   T=T+W(I,K)*V(I,J)
 35   EPSI(J)=T
      DO 50 L=1,K
      T=0.
      DO 40 J=1,NB
 40   T=T+EPSI(J)*W(J,L)
      KL=KL+1
 50   VPHI(KL)=T
      RETURN
      END
      SUBROUTINE JACVAO(A,N,EIVU,EIVR)                                  00000100
      IMPLICIT REAL*8 (A-H,O-Z)                                         00000200
      DIMENSION A(N,N),EIVU(N),EIVR(N,N)                                00000300
      IF(N-1) 2,2,1                                                     00000400
    2 EIVR(1,1)=1.0D0                                                   00000500
      EIVU(1)=A(1,1)                                                    00000600
      RETURN                                                            00000700
    1 DO 101 J=1,N                                                      00000800
      DO 100 I=1,N                                                      00000900
  100 EIVR(I,J)=0.0D0                                                   00001000
  101 EIVR(J,J)=1.0D0                                                   00001100
C        FIND THE ABSOLUTELY LARGEST ELEMENT OF A                       00001200
  102 ATOP=0.D0                                                         00001300
      DO 111 I=1,N                                                      00001400
      DO 111 J=I,N                                                      00001500
      IF (ATOP-DABS(A(I,J))) 104,111,111                                00001600
  104 ATOP=DABS(A(I,J))                                                 00001700
  111 CONTINUE                                                          00001800
      IF(ATOP)109,109,113                                               00001900
  109 RETURN                                                            00002000
C        CALCULATE THE STOPPING CRITERION -- DSTOP                      00002100
  113 AVGF=DFLOAT(N*(N-1))*.55D0                                        00002200
      D=0.0D0                                                           00002300
      DO 114 JJ=2,N                                                     00002400
      DO 114 II=2,JJ                                                    00002500
      S=A(II-1,JJ)/ATOP                                                 00002600
  114 D=S*S+D                                                           00002700
      DSTOP=(1.D-09)*D                                                  00002800
C        CALCULATE THE THRESHOLD, THRSH                                 00002900
C                                                                       00003000
      THRSH=DSQRT(D/AVGF)*ATOP                                          00003100
C                                                                       00003200
C        START A SWEEP                                                  00003300
C                                                                       00003400
  115 IFLAG=0                                                           00003500
      DO 130 JCOL=2,N                                                   00003600
      JCOL1=JCOL-1                                                      00003700
      DO 130 IROW=1,JCOL1                                               00003800
      AIJ=A(IROW,JCOL)                                                  00003900
C                                                                       00004000
C        COMPARE THE OFF-DIAGONAL ELEMENT WITH THRSH                    00004100
C                                                                       00004200
      IF (DABS(AIJ)-THRSH) 130,130,117                                  00004300
  117 AII=A(IROW,IROW)                                                  00004400
      AJJ=A(JCOL,JCOL)                                                  00004500
      S=AJJ-AII                                                         00004600
C                                                                       00004700
C        CHECK TO SEE IF THE CHOSEN ROTATION IS LESS THAN THE ROUNDING E00004800
C        IF SO , THEN DO NOT ROTATE.                                    00004900
C                                                                       00005000
      IF (DABS(AIJ)-1.D-09*DABS(S)) 130,130,118                         00005100
  118 IFLAG=1                                                           00005200
C                                                                       00005300
C        IF THE ROTATION IS VERY CLOSE TO 45 DEGREES, SET SIN AND COS   00005400
C     TO 1/(ROOT 2).                                                    00005500
C                                                                       00005600
      IF (1.D-10*DABS(AIJ)-DABS(S)) 116,119,119                         00005700
  119 S=.707106781186548D0                                              00005800
      C=S                                                               00005900
      GO TO 120                                                         00006000
C                                                                       00006100
C        CALCULATION OF SIN AND COS FOR ROTATION THAT IS NOT VERY CLOSE 00006200
C        TO 45 DEGREES                                                  00006300
C                                                                       00006400
  116 T=AIJ/S                                                           00006500
      S=0.25D0/DSQRT(0.25D0+T*T)                                        00006600
C                                                                       00006700
C     COS=C,  SIN=S                                                     00006800
C                                                                       00006900
      C=DSQRT(0.5D0+S)                                                  00007000
      S=2.D0*T*S/C                                                      00007100
C                                                                       00007200
C        CALCULATION OF THE NEW ELEMENTS OF MATRIX A                    00007300
C                                                                       00007400
  120 DO 121 I=1,IROW                                                   00007500
      T=A(I,IROW)                                                       00007600
      U=A(I,JCOL)                                                       00007700
      A(I,IROW)=C*T-S*U                                                 00007800
  121 A(I,JCOL)=S*T+C*U                                                 00007900
      I2=IROW+2                                                         00008000
      IF(I2.GT.JCOL) GO TO 123                                          00008100
      DO 122 I=I2,JCOL                                                  00008200
      T=A(I-1,JCOL)                                                     00008300
      U=A(IROW,I-1)                                                     00008400
      A(I-1,JCOL)=S*U+C*T                                               00008500
  122 A(IROW,I-1)=C*U-S*T                                               00008600
  123 A(JCOL,JCOL)=S*AIJ+C*AJJ                                          00008700
      A(IROW,IROW)=C*A(IROW,IROW)-S*(C*AIJ-S*AJJ)                       00008800
      DO 124 J=JCOL,N                                                   00008900
      T=A(IROW,J)                                                       00009000
      U=A(JCOL,J)                                                       00009100
      A(IROW,J)=C*T-S*U                                                 00009200
  124 A(JCOL,J)=S*T+C*U                                                 00009300
C                                                                       00009400
C        ROTATION COMPLETED.                                            00009500
C                                                                       00009600
  131 DO 125 I=1,N                                                      00009700
      T=EIVR(I,IROW)                                                    00009800
      EIVR(I,IROW)=C*T-EIVR(I,JCOL)*S                                   00009900
  125 EIVR(I,JCOL)=S*T+EIVR(I,JCOL)*C                                   00010000
C                                                                       00010100
C        CALCULATE THE NEW NORM D AND COMPARE WITH DSTOP                00010200
C                                                                       00010300
      S=AIJ/ATOP                                                        00010400
      D=D-S*S                                                           00010500
      IF(D.GE.DSTOP) GO TO 129                                          00010600
C                                                                       00010700
C        RECALCULATE DSTOP AND THRSH TO DISCARD ROUNDING ERRORS         00010800
C                                                                       00010900
      D=0.D0                                                            00011000
      DO 128 JJ=2,N                                                     00011100
      DO 128 II=2,JJ                                                    00011200
      S=A(II-1,JJ)/ATOP                                                 00011300
  128 D=S*S+D                                                           00011400
      DSTOP=(1.D-09)*D                                                  00011500
  129 THRSH=DSQRT(D/AVGF)*ATOP                                          00011600
  130 CONTINUE                                                          00011700
      IF(IFLAG.GT.0) GO TO 115                                          00011800
C                                                                       00011900
C     PLACE EIGENVALUES IN EIVU                                         00012000
      DO 132 J=1,N                                                      00012100
      EIVU(J)=A(J,J)                                                    00012200
  132 CONTINUE                                                          00012300
      RETURN                                                            00012400
      END                                                               00012500
      SUBROUTINE OVERFL(K)
      K=2
      RETURN
      END
