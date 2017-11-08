C 
C PROGRAMME NLOCVS
C
C
      OPEN (UNIT=21,form='UNFORMATTED',status='SCRATCH')
      OPEN (UNIT=20,form='UNFORMATTED',file='PSNL',status='OLD')
      CALL NOLOC                                                        
      STOP                                                              
      END                                                               
      SUBROUTINE NOLOC                                                  
      IMPLICIT REAL*8(A-H,O-T,V-Z)                                      
      DIMENSION INUM(40),EPSDS(40),Z(40),ZOP(40),RNOR(40),OPNOR(40)     
      DIMENSION TOPE(1600),TORT(1600),REC(1600),TEMP(1600),VKL(840)     
      DIMENSION SW(400),SV(400), SN(400),SC(400)                        
      DIMENSION SL(1600)                                                
      DIMENSION ZOP1(20),TEMP1(400),EPSDS1(20)                          
      DIMENSION ZOP2(20),TEMP2(400),EPSDS2(20)                          
      EQUIVALENCE ( ZOP (1),ZOP2(1) ), ( TEMP (1),TEMP2(1) )            
      EQUIVALENCE ( EPSDS (1),EPSDS2(1) )                               
      DIMENSION OLDPS(840),VERI(840)                                    
      COMMON VL(44500)                                                  
      DIMENSION NPOT(20),APOT(20),CPOT(20)                              
      COMMON/GIV/NBLOCK(400),NSTART(400),NSIZE,IMNM(400)                
      CHARACTER*8 FMT,apsd,anam                                         
      DIMENSION RZ(4),OVZ(4),EX(4),OZ(2)                                
      NAMELIST/INPUT/ ZN,L,NBAS,NOP,ARAF,BRAF,NVOIS,AOP,BOP,NBPOT,NPOT, 
     &APOT,CPOT,Z,ZOP,NEWBAS,DMNXP2,DMNXP,OPTES,BASTES,APSD             
     &,FMT,IDEN,UPUNCH,USO,NENR,UREPL,ULISTE                            
     &,RZ,OVZ,UCUT,OZ                                                   
     U,UPROJ                                                            
     U,UOLDPS,UVERI                                                     
     U,USO,NENR,NOM1,NOM2                                               
      LOGICAL   NEWBAS ,UCUT,USO, UREPL,ULISTE                          
      LOGICAL IVERI,UPUNCH                                              
      LOGICAL   UPROJ,UOLDPS,UVERI                                      
      I=INIPUR(1)                                                       
      UPROJ=.FALSE.                                                     
      UOLDPS=.FALSE.                                                    
      UPUNCH=.FALSE.                                                    
      UVERI=.FALSE.                                                     
      USO=.FALSE.                                                       
      UREPL=.FALSE.                                                     
      DATA FMT/'(5Z16)'/                                                
      MAXBAS=19                                                         
      MAXBAS=20                                                         
      UCUT=.FALSE.                                                      
       ULISTE=.FALSE.                                                   
      MDIM=400                                                          
      RZ(1)=0.                                                          
      OZ(1)=0.                                                          
      OPTES=1.D-6                                                       
      BASTES=1.D-6                                                      
      DMNXP=1.D-10                                                      
      ARAF=0.                                                           
      IVERI=.TRUE.                                                      
      DO 2300 I=1,400                                                   
2300     IMNM(I)=I*(I-1)/2                                              
      AOP=0.                                                            
      OZ(1)=0.                                                          
1     CONTINUE                                                          
      DMNXP2=1.D-9                                                      
      NEWBAS=.FALSE.                                                    
      READ(5,input,END=999)
      WRITE(6,INPUT)                                                    
C     DEFINE FILE 20 (94,882,U,NAW)                                     
C     DEFINE FILE 30 (30,885,U,NAW)                                     
      IF(NEWBAS) GO TO 8884                                             
      WRITE(6,8888) IDEN,APSD,L,NOP,RZ,OVZ,OZ                           
 8888 FORMAT(//////////,130(1H*),/,' EXTRACTION PSEUDO NON-LOCAUX POUR  
     1 NUMERO ATOMIQUE',I4,10X,'ATOME',3X,A8,'SYMETRIE',I3,/,           
     2 '  NOMBRE D''OPERATEURS ',I3,                                    
     3 /,' EXTENSION DE LA BASE D''EXTRACTION',4F10.4,/,' RECOUVREMENTS'
     4,4F10.4,/,' EXTENSION DE LA BASE D''OPERATEURS',2F10.4)           
      WRITE(6,8887) NBPOT                                               
      WRITE(6,8886) (NPOT(I),APOT(I),CPOT(I),I=1,NBPOT)                 
 8887 FORMAT(/,'  OPERATEUR SEMI-LOCAL  NOMBRE DE TERMES=',I4)          
 8886 FORMAT(I4,2F12.6)                                                 
 8884 CONTINUE                                                          
      LL=L+L                                                            
      IF(RZ(1).EQ.0.D0) GO TO 2200                                      
      DO 2050 K=1,4                                                     
      EX(K)=L/(2.D0*RZ(K)*RZ(K))                                        
2050   CONTINUE                                                         
      J=1                                                               
      Z(1)=EX(1)                                                        
      DO 2150 K=1,3                                                     
      T=OVZ(K)**(2.D0/(DFLOAT(L)+0.5D0))                                
      T=(2.D0-T+DSQRT((2.D0-T)**2-T))/T                                 
      DO 2100 I=1,39                                                    
      J=J+1                                                             
      Z(J)=Z(J-1)*T                                                     
      IF(Z(J).GT.EX(K+1)) GO TO 2110                                    
2100    CONTINUE                                                        
2110   J=J-1                                                            
2150    CONTINUE                                                        
      NBAS=J                                                            
      IJ=0                                                              
      DO 40 I=1,NBAS                                                    
      RNOR(I)=1.D0/DSQRT(PURX(LL,Z(I)+Z(I)))                            
      DO 40 J=1,I                                                       
      ZZ=Z(I)+Z(J)                                                      
      T=0.D0                                                            
      DO 35 K=1,NBPOT                                                   
   35 T=T+PURX(LL+NPOT(K),ZZ+APOT(K))*CPOT(K)                           
      T=T*RNOR(I)*RNOR(J)                                               
      IJ=IJ+1                                                           
      VKL(IJ)=T                                                         
   40 CONTINUE                                                          
      IJ=0                                                              
      WRITE(6,8883)                                                     
 8883 FORMAT(/,'  ELEMTS DE MATRICE DS LA BASE D''EXTRACTION')          
      DO 45 I=1,NBAS                                                    
      WRITE(6,42) (VKL(IJ+J),J=1,I)                                     
  45  IJ=IJ+I                                                           
C                                                                       
C     VERIFICATION A PARTIR DE DONNEES SUR FILE 20 (NEWBAS=T ET UPUNCH=T
C                                                                       
      IF(.NOT.NEWBAS) GOTO 3900                                         
      REWIND 20                                                         
 102  READ(20,END=1750) ANAM,NPMAX                                      
      IF(ANAM.EQ.APSD) GOTO 110                                         
      DO 105 J=1,NPMAX                                                  
 105  READ(20)                                                          
      GOTO 102                                                          
 1750 WRITE(6,9995) APSD                                                
 9995 FORMAT(//,'   PAS DE PSEUDO-POTENTIEL SUR FT20 POUR =',A6,//)     
      GOTO 1                                                            
 110  IF(L.EQ.1) GOTO 120                                               
      L1=L-1                                                            
      DO 115 I=1,L1                                                     
 115  READ(20)                                                          
 120  READ(20) NOP,NVOIS,ZOP2,TEMP2,EPSDS2                              
      IF (NOP.NE.0) GOTO 125                                            
      WRITE(6,9996) APSD,L                                              
      GOTO 1                                                            
 9996 FORMAT(//,'   PAS DE PSEUDO-POTENTIEL SUR FT20 POUR=',A6,         
     1'  DE SYMETRIE L=',I3)                                            
 125  CONTINUE                                                          
      WRITE(6,3678) APSD,NOP,NVOIS,(ZOP(I),I=1,NOP)                     
      WRITE(6,3681)                                                     
      WRITE(6,3679) (EPSDS(I),I=1,NVOIS)                                
      WRITE(6,3682)                                                     
      IJ=0                                                              
      DO 3654 I=1,NVOIS                                                 
      WRITE(6,3679) (TEMP(IJ+J),J=1,NOP)                                
 3654 IJ=IJ+NOP                                                         
 3678 FORMAT(////,' OP. NON LOCAL POUR ATOME',5X,A8,//,                 
     *' NOMBRE DE PRIMITIVES',I3,/,                                     
     *' NOMBRE DE CONTRACTEES',I3,/,                                    
     *' EXPOSANTS',/,(10F12.6))                                         
 3679 FORMAT(1X,8D15.8)                                                 
 3681 FORMAT(' COEFFICIENTS DE L''OP NON LOCAL')                        
 3682 FORMAT(' COEFFICIENTS DE CONTRACTION')                            
      DO 3100 I=1,NOP                                                   
 3100 OPNOR(I)=1.D0/DSQRT(PURX(LL,ZOP(I)+ZOP(I)))                       
      IJ=0                                                              
      DO 3200 I=1,NBAS                                                  
      DO 3200 J=1,NOP                                                   
      IJ=IJ+1                                                           
 3200 REC(IJ)=PURX(LL,Z(I)+ZOP(J))*RNOR(I)*OPNOR(J)                     
      CALL ORTHOE(REC,SV,TEMP,NOP,NVOIS,NBAS)                           
      IJ=0                                                              
      DO 3300 I=1,NBAS                                                  
      NI=(I-1)*NVOIS                                                    
      DO 3300 J=1,I                                                     
      NJ=(J-1)*NVOIS                                                    
      IJ=IJ+1                                                           
      T=0.D0                                                            
      DO 3350 K=1,NVOIS                                                 
 3350 T=T+REC(NI+K)*REC(NJ+K)*EPSDS(K)                                  
 3300 VERI(IJ)=T-VKL(IJ)                                                
      WRITE(6,6)                                                        
      WRITE(6,3400)                                                     
 3400 FORMAT(' :::::::: VERIFICATION :::::::')                          
      IJ=0                                                              
      DO 3410 I=1,NBAS                                                  
      WRITE(6,585) (VERI(IJ+J),J=1,I)                                   
 3410 IJ=IJ+I                                                           
      GOTO 1                                                            
 3900 CONTINUE                                                          
      IF(NEWBAS)GO TO 2160                                              
      IF(OZ(1).EQ.0.D0) GO TO  2160                                     
      IF(NOP.EQ.0) NOP=NBAS                                             
      AOP=L/(2.D0*OZ(1)*OZ(1))                                          
      BOP=L/(2.D0*OZ(2)*OZ(2))                                          
      BOP=(BOP/AOP)**(1.D0/DFLOAT(NOP))                                 
2160     CONTINUE                                                       
      IF(NBAS.LT.40) GO TO 2250                                         
      WRITE(6,2210)                                                     
2210   FORMAT(' NOMBRE D EXP GENERES >40')                              
      STOP                                                              
2200   CONTINUE                                                         
      IF(ARAF.EQ.0.) GO TO 11                                           
      Z(1)=ARAF                                                         
      DO 10 I=2,NBAS                                                    
10    Z(I)=Z(I-1)*BRAF                                                  
2250     CONTINUE                                                       
11     IF(AOP.EQ.0.) GO TO 12                                           
      IF(NEWBAS) GO TO 12                                               
      ZOP(1)=AOP                                                        
      DO 20 I=2,NOP                                                     
20    ZOP(I)=ZOP(I-1)*BOP                                               
      IF(NVOIS.EQ.0) NVOIS=NOP                                          
12      CONTINUE                                                        
      WRITE(6,6003)(Z(I),I=1,NBAS)                                      
6003  FORMAT(/,' EXPOSANTS DE LA BASE D''EXTRACTION',/,(10F12.6))       
      WRITE(6,6004)(ZOP(I),I=1,NOP)                                     
6004  FORMAT(/,' EXPOSANTS DE LA BASE D''OPERATEURS',/,(10F12.6))       
      I=INIPUR(1)                                                       
      IF(NEWBAS) GO TO 721                                              
      DO 21  I=1,NOP                                                    
      OPNOR(I)=1.D0/DSQRT(PURX(LL,ZOP(I)+ZOP(I)))                       
21     CONTINUE                                                         
       IJ=0                                                             
      DO 710 I=1,NOP                                                    
      DO 710 J=1,I                                                      
      IJ=IJ+1                                                           
710    SL(IJ)=OPNOR(I)*OPNOR(J)*PURX(LL,  ZOP(I)+ZOP(J))                
      WRITE(6,759)                                                      
      WRITE(6,31) SL(2),SL(4),SL(7)                                     
759     FORMAT('0 ORTHONORM DES  OPERATEURS')                           
      CALL DIAGOS(SL,TOPE,TEMP,EPSDS,NOP,NOPR,DMNXP2)                   
       DO 720 J=1,NOPR                                                  
      IF(EPSDS(J).GT.OPTES) NVOISI=NVOISI+1                             
      INUM(J)=J                                                         
720   CONTINUE                                                          
761   FORMAT(1X,10F12.7)                                                
721       CONTINUE                                                      
      IJ=0                                                              
      JI=0                                                              
      DO 30 I=1,NBAS                                                    
      ZZ=Z(I)                                                           
      RNOR(I)=1.D0/DSQRT(PURX(LL,ZZ+ZZ))                                
      DO 25 J=1,I                                                       
      T=RNOR(I)*RNOR(J)*PURX(LL,ZZ+Z(J))                                
      JI=JI+1                                                           
25    SL(JI)=T                                                          
      DO 30 J=1,NOP                                                     
      IJ=IJ+1                                                           
30    REC(IJ)=PURX(LL,ZZ+ZOP(J))*RNOR(I)*OPNOR(J)                       
      CALL DIAGOS(SL,TORT,TEMP,EPSDS,NBAS,NBASLD,DMNXP2)                
      WRITE(6,31) SL(2),SL(4),SL(7)                                     
31      FORMAT(' RECOUVREMENT ENTRE VOISINS:',3F15.8)                   
      FNORMO=1.D10                                                      
      NBAI=0                                                            
      DO 730 I=1,NBASLD                                                 
      IF(EPSDS(I).GT.BASTES) NBAI=NBAI+1                                
730    CONTINUE                                                         
        NBASLD=NBAI                                                     
      IF(UPROJ       ) GO TO 620                                        
      IF(UCUT.AND.NOPR.GT.NBASLD) NOPR=NBASLD                           
      IF(NOPR.GT.MAXBAS) NOPR=MAXBAS                                    
      GO TO 640                                                         
620    IJ=0                                                             
      IF(NEWBAS) GO TO 640                                              
      DO 630 I=1,NOP                                                    
      DO 630 J=1,I                                                      
      ZZ=ZOP(I)+ZOP(J)                                                  
      T=0.                                                              
      DO 625 K=1,NBPOT                                                  
625    T=T+PURX(LL+NPOT(K),ZZ+APOT(K))*CPOT(K)                          
      T=T*OPNOR(I)*OPNOR(J)                                             
      IJ=IJ+1                                                           
      SC(IJ)=T                                                          
630       CONTINUE                                                      
      CALL ORTHO(SC,TEMP,SV,TOPE,NOP,NOPR)                              
       IJ=0                                                             
       DO 635 I=1,NOP                                                   
       IJ=IJ+I                                                          
635     SC(IJ)=0.5D0*SC(IJ)                                             
      WRITE(6,6)                                                        
      WRITE(6,*) (SC(I),I=1,IJ)                                         
640    CONTINUE                                                         
      IF(NVOISI.GT.NBAI)NVOISI=NBAI                                     
      CALL ORTHO(VKL,TEMP,SV,TORT,NBAS,NBASLD)                          
        IJ=NBAS*(NBAS+1)/2                                              
      WRITE(6,6)                                                        
6     FORMAT('0',20(' * '),/,'0')                                       
      CALL ORTHOD    (REC,SV,TORT,NBAS,NBASLD,NOP)                      
      CALL ORTHOE(REC,SV,TOPE,NOP,NOPR,NBASLD)                          
42    FORMAT(' VKL=',12F10.4)                                           
      NBASS=NBAS                                                        
      NBAS=NBASLD                                                       
      NOPS=NOP                                                          
      NOP=NOPR                                                          
        NOPEF=0                                                         
      JI=0                                                              
      NVOIS=NOP                                                         
      IF(.NOT.UOLDPS) GO TO 680                                         
      IJ=0                                                              
      DO 670 I=1,NBAS                                                   
      DO 670 J=1,I                                                      
      IJ=IJ+1                                                           
670   VKL(IJ)=VKL(IJ)-OLDPS(IJ)                                         
680     CONTINUE                                                        
      IF(.NOT. UVERI) GO TO 690                                         
      IF(IVERI) GO TO 690                                               
      IJ=0                                                              
      DO 685 I=1,NBAS                                                   
      DO 685 J=1,I                                                      
      IJ=IJ+1                                                           
685   VKL(IJ)=VKL(IJ)-VERI(IJ)                                          
690   CONTINUE                                                          
      WRITE(6,53) NVOIS                                                 
53    FORMAT(' NB BANDES COMPATIBLES AVEC MDIM',I4)                     
      IF(NEWBAS) GO TO 349                                              
      IF(UPROJ) GO TO 349                                               
      I2=NOP                                                            
      IJJI=0                                                            
      IJ=0                                                              
      DO 200 II=1,NVOIS                                                 
      I=INUM(II)                                                        
      DO 190 IIV=1,II                                                   
      IV=INUM(IIV)                                                      
      IP=IV-I                                                           
      T=REC(I)*REC(IV)                                                  
      TEMP(1)=T+T                                                       
      MN=1                                                              
      TT=T*VKL(1)                                                       
      MI=NOP+I                                                          
      DO 60 M=2,NBAS                                                    
      MN1=M-1                                                           
       REI=REC(MI)                                                      
      REV=REC(MI+IP)                                                    
      NI=I                                                              
      DO 55 N=1,MN1                                                     
      MN=MN+1                                                           
      T=REI*REC(NI+IP)+REV*REC(NI)                                      
      TT=TT+T*VKL(MN)                                                   
      NI=NI+NOP                                                         
55    TEMP(MN)=T                                                        
      MI=MI+NOP                                                         
      MN=MN+1                                                           
      T=REI*REV                                                         
      TEMP(MN)=T+T                                                      
60    TT=TT+T*VKL(MN)                                                   
      IJ=IJ+1                                                           
      SV(IJ)=TT                                                         
      JI=0                                                              
      DO 180 JJ=1,NVOIS                                                 
      J=INUM(JJ)                                                        
      DO 100 JJV=1,JJ                                                   
      JV=INUM(JJV)                                                      
      JP=JV-J                                                           
      JI=JI+1                                                           
      TT=TEMP(1)*REC(J)*REC(JV)                                         
      MN=1                                                              
      MJ=NOP+J                                                          
      DO 80 M=2,NBAS                                                    
      MN1=M-1                                                           
      REJ=REC(MJ)                                                       
      REV=REC(MJ+JP)                                                    
      NJ=J                                                              
      DO 70 N=1,MN1                                                     
      MN=MN+1                                                           
      TT=TT+TEMP(MN)*(REJ*REC(NJ+JP)+REV*REC(NJ))                       
70    NJ=NJ+NOP                                                         
      MJ=MJ+NOP                                                         
      MN=MN+1                                                           
      TT=TT+TEMP(MN)*REJ*REV                                            
80    CONTINUE                                                          
      IJJI=IJJI+1                                                       
      VL(IJJI)=TT                                                       
      IF(JI.NE.IJ) GO TO 100                                            
      TT=1.D0/DSQRT(TT)                                                 
      SV(IJ)=SV(IJ)*TT                                                  
      SN(IJ)=TT                                                         
      N=IJJI-IJ                                                         
      DO 90 M=1,IJ                                                      
90    VL(N+M)=TT*SN(M)*VL(N+M)                                          
      GO TO 185                                                         
100   CONTINUE                                                          
180   CONTINUE                                                          
185   CONTINUE                                                          
190   CONTINUE                                                          
200   I2=I2-1                                                           
      CALL DIAGOS(VL,VL,VL,TEMP,IJ,IJP,DMNXP)                           
      WRITE(6,206) IJ,IJP                                               
206   FORMAT(' NB TERME',I4,' NB REDUIT',I4)                            
      DO 300 I=1,IJ                                                     
      TT=0.                                                             
      DO 260 J=1,IJ                                                     
      T=0.                                                              
      JM=J                                                              
      IM=I                                                              
      DO 250 M=1,IJP                                                    
      T=T+VL(IM)*VL(JM)                                                 
      JM=JM+IJ                                                          
250   IM=IM+IJ                                                          
260   TT=TT+T*SV(J)                                                     
300   SC(I)=TT                                                          
        DO 850 I=1,IJ                                                   
 850    SC(I)=SC(I)*SN(I)                                               
      WRITE(6,*) (SC(I),I=1,IJ)                                         
349     CONTINUE                                                        
      IJ=NVOIS*(NVOIS+1)/2                                              
348     CONTINUE                                                        
      MN=0                                                              
      DO 350 I=1,NBAS                                                   
      DO 350 J=1,I                                                      
      MN=MN+1                                                           
350   TEMP(MN)=0.                                                       
      IJ=0                                                              
      IP=0                                                              
      I2=NOP                                                            
      I1=1                                                              
      DO 500 II=1,NVOIS                                                 
      I=INUM(II)                                                        
      DO 490 IIV=1,II                                                   
      IV=INUM(IIV)                                                      
      IP=IV-I                                                           
        IJ=IJ+1                                                         
      T=REC(I)*REC(IV)                                                  
      TEMP(1)=(T+T)*SC(IJ)+TEMP(1)                                      
      MN=1                                                              
      MI=NOP+I                                                          
      DO460 M=2,NBAS                                                    
      MN1=M-1                                                           
       REI=REC(MI)                                                      
      REV=REC(MI+IP)                                                    
      NI=I                                                              
      DO455 N=1,MN1                                                     
      MN=MN+1                                                           
      T=REI*REC(NI+IP)+REV*REC(NI)                                      
      NI=NI+NOP                                                         
      TEMP(MN)=T   *SC(IJ)+TEMP(MN)                                     
455    CONTINUE                                                         
      MN=MN+1                                                           
      T=REI*REV                                                         
      TEMP(MN)=(T+T)*SC(IJ)+TEMP(MN)                                    
      MI=MI+NOP                                                         
460   CONTINUE                                                          
490       CONTINUE                                                      
      IP=IP+1                                                           
      I2=I2-1                                                           
500        CONTINUE                                                     
      WRITE(6,6)                                                        
      WRITE(6,6)                                                        
      IJ=0                                                              
      DO 610 I=1,NBAS                                                   
      DO 610 J=1,I                                                      
      IJ=IJ+1                                                           
      OLDPS(IJ)=TEMP(IJ)                                                
610   CONTINUE                                                          
      IF(.NOT.UVERI) GO TO 615                                          
      IVERI=.NOT.IVERI                                                  
      IJ=0                                                              
      DO 614 I=1,NBAS                                                   
      DO 614 J=1,I                                                      
      IJ=IJ+1                                                           
614    VERI(IJ)=TEMP(IJ)                                                
615   CONTINUE                                                          
        IJ=NBAS*(NBAS+1)/2                                              
      T=TEMP(1)-VKL(1)                                                  
      TEMP(1)=T                                                         
      MN=1                                                              
      FNORM=T*T                                                         
      DO 520 M=2,NBAS                                                   
      MN1=M-1                                                           
      DO 515 N=1,MN1                                                    
      MN=MN+1                                                           
      T=TEMP(MN)-VKL(MN)                                                
      TEMP(MN)=T                                                        
515     FNORM=FNORM+2.D0*T*T                                            
      MN=MN+1                                                           
      T=TEMP(MN)-VKL(MN)                                                
      TEMP(MN)=T                                                        
520    FNORM=FNORM+T*T                                                  
510     FORMAT(' NORME',F20.8)                                          
         WRITE(6,510) FNORM                                             
      IJ=0                                                              
      DO 550 J=1,NBAS                                                   
      JI=J                                                              
      T=EPSDS(J)                                                        
      DO 550 I=1,NBASS                                                  
      IJ=IJ+1                                                           
      VL(JI)=TORT(IJ)*T                                                 
550   JI=JI+NBAS                                                        
      CALL ORTHO(TEMP,SL,SV,VL,NBAS,NBASS)                              
      IJ=0                                                              
      WRITE(6,584)                                                      
584   FORMAT('0  DIFFERENCE DES ELMTS DE MAT DU PSEUDO DS LA BASE INI') 
      DO 580 I=1,NBASS                                                  
      WRITE(6,585) (TEMP(IJ+J),J=1,I)                                   
580   IJ=IJ+I                                                           
      NOP=NOPS                                                          
       IF(NEWBAS)GO TO 1                                                
 997    NTER=NVOIS*(NVOIS+1)/2                                          
585   FORMAT(1X,12D10.2)                                                
       FNORM=FNORM/NBAS                                                 
      IF(FNORM.GT.1.D0) GO TO 1                                         
         IJV=0                                                          
         DO 1950I=1,NVOIS                                               
         DO1950J=1,I                                                    
          IJV=IJV+1                                                     
          T=SC(IJV)                                                     
          IF(I.EQ.J)T=T+T                                               
         TEMP(NVOIS*(I-1)+J)=T                                          
1950     TEMP(NVOIS*(J-1)+I)=T                                          
         CALL JACVAO(TEMP,NVOIS,EPSDS,TORT)                             
        DO2000I=1,NVOIS                                                 
        DO2000J=1,NOP                                                   
         JK=J                                                           
         KI=NVOIS*(I-1)                                                 
         T=0.D0                                                         
         DO1900K=1,NVOIS                                                
         KI=KI+1                                                        
         T=T+TORT(KI)*TOPE(JK)                                          
1900      JK=JK+NOP                                                     
2000       TEMP(NOP*(I-1)+J)=T                                          
      IF(UPUNCH) GO TO 5647                                             
C......STOCKAGE SUR FILE20                                              
C......RELECTURE DE FILE20                                              
C......TEST POUR SAVOIR SI CE CAS EXISTE DEJA                           
 699  CONTINUE                                                          
      REWIND 20                                                         
      REWIND 21                                                         
C 700  READ(20,END=750,ERR=787) ANAM,NPMAX                              
 700  READ(20,END=750,ERR=750) ANAM,NPMAX                               
      IF (ANAM.EQ.APSD) GOTO 705                                        
      WRITE(21) ANAM,NPMAX                                              
      DO 702 I=1,NPMAX                                                  
      READ(20,END=780,ERR=780)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1           
 702  WRITE(21)NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                            
      GOTO 700                                                          
 705  IF( L.LE.NPMAX)      GOTO 1720                                    
      NP1=NPMAX+1                                                       
C                                                                       
C.....PSEUDO EXISTE DEJA. ON AJOUTE LA SYMETRIE SUIVANTE                
C                                                                       
      WRITE(21) ANAM,NP1                                                
      DO 1710 I=1,NPMAX                                                 
      READ(20,END=790,ERR=790)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1           
 1710  WRITE(21) NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                          
C                                                                       
C.....SYMETRIE ACTUELLEMENT CALCULEE                                    
C                                                                       
      WRITE(21) NOP,NVOIS,ZOP2,TEMP2,EPSDS2                             
C                                                                       
C.....ON COMPLETE LA FILE21                                             
C                                                                       
 715  READ(20,END=760,ERR=789) ANAM,NPMAX                               
         WRITE(6,783) ANAM,NPMAX                                        
 783     FORMAT(' *****715****',A6,I6)                                  
      WRITE(21) ANAM ,NPMAX                                             
      DO 718 I=1,NPMAX                                                  
      READ(20,END=795,ERR=795)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1           
 718  WRITE(21) NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                           
      GOTO 715                                                          
 1720  IF(.NOT.UREPL) GOTO 800                                          
C                                                                       
C.....PSEUDO EXISTE DEJA. ON REMPLACE UNE SYMETRIE DEJA EXISTANTE       
C                                                                       
      WRITE(21) ANAM,NPMAX                                              
      IF (L.EQ.1) GOTO  726                                             
      L1=L-1                                                            
      DO 725 I=1,L1                                                     
      READ(20,END=796,ERR=796)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1           
 725  WRITE(21) NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                           
 726  CONTINUE                                                          
      READ(20,END=797,ERR=797)                                          
      WRITE(21) NOP ,NVOIS ,ZOP2,TEMP2,EPSDS2                           
      IF(L.EQ.NPMAX)      GOTO 715                                      
      L1=L+1                                                            
      DO 1730 I=L1,NPMAX                                                
      READ(20,END=798,ERR=798)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1           
 1730  WRITE(21) NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                          
C                                                                       
C.....ON COMPLETE LA FILE21                                             
C                                                                       
      GOTO 715                                                          
 750  CONTINUE                                                          
C                                                                       
C.....LE PSEUDO N'EXISTE PAS SUR LA FILE                                
C.....ON LE RAJOUTE A LA FIN                                            
C                                                                       
      NPMAX=1                                                           
      WRITE(21) APSD,NPMAX                                              
      WRITE(21) NOP,NVOIS,ZOP2,TEMP2,EPSDS2                             
C                                                                       
C.....REECRITURE DE LA FILE20 MISE A JOUR                               
C                                                                       
 760  CONTINUE                                                          
      REWIND 20                                                         
      REWIND 21                                                         
 765  READ(21,END=1) ANAM,NPMAX                                         
      WRITE(20) ANAM,NPMAX                                              
      WRITE(6,9843) ANAM,NPMAX                                          
      DO 770 I=1,NPMAX                                                  
      READ(21)  NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                           
      WRITE(20) NOP1,NVOIS1,ZOP1,TEMP1,EPSDS1                           
      IF(ULISTE) WRITE(6,9844) I,NOP1,NVOIS1,ZOP1                       
 770  CONTINUE                                                          
 9843 FORMAT(/'    ATOME=',A8,'   NB DE SYMETRIE ',I3)                  
 9844 FORMAT('   SYMETRIE=',I3,'  NB D OPERATEURS=',I3,                 
     1'  NB DE COMB.ORTHOG.=',I3 ,/,' EXPOSANTS=',12F10.6,(/,13F10.6))  
      GOTO 765                                                          
C                                                                       
C ERREUR LECTURE                                                        
C                                                                       
 787    ICAS=700                                                        
        WRITE(6,781) ANAM,NPMAX,ICAS                                    
        STOP                                                            
 789     ICAS=715                                                       
        WRITE(6,781) ANAM,NPMAX,ICAS                                    
        STOP                                                            
 780    ICAS=702                                                        
        WRITE(6,781) ANAM,NPMAX,ICAS                                    
        STOP                                                            
 790    ICAS=1710                                                       
        WRITE(6,781) ANAM,NPMAX,ICAS                                    
        STOP                                                            
 795    ICAS=718                                                        
        WRITE(6,781) ANAM,NPMAX,ICAS                                    
        STOP                                                            
 796     ICAS=725                                                       
        WRITE(6,781) ANAM,NPMAX,ICAS                                    
        STOP                                                            
 797      ICAS=0                                                        
        WRITE(6,781) ANAM,NPMAX,ICAS                                    
        STOP                                                            
 798      ICAS=1730                                                     
        WRITE(6,781) ANAM,NPMAX,ICAS                                    
        STOP                                                            
 781   FORMAT(//' ERREUR SUR FILE20 APRES ',A6,I3,'  ICAS=',I5)         
C                                                                       
C..... ERREUR....ON ESSAIE DE REMPLACER UN PSEUDO EXISTANT              
C      SANS SPECIFIER  UREPL=.TRUE.                                     
C                                                                       
 800  CONTINUE                                                          
       WRITE(6,9845)                                                    
 9845  FORMAT(//,'*****ESSAI DE REMPLACEMENT D UN PSEUDO SANS SPECIFIER 
     1UREPL=.TRUE.********')                                            
      STOP                                                              
 5647 CONTINUE                                                          
      WRITE(8,FMT) IDEN,L,NOP,NVOIS,NVOIS                               
      WRITE(8,FMT) (ZOP(I),I=1,NOP)                                     
      WRITE(8,FMT)((TEMP( I+(J-1)*NOP),I=1,NOP),J=1,NVOIS)              
                        WRITE(8,FMT)(EPSDS(I),I=1,NVOIS    )            
 5650 CONTINUE                                                          
999    CONTINUE                                                         
      STOP                                                              
      END                                                               
      subroutine given (h,e,v,n,nev,nvec)
      implicit real*8 (a-h,o-x,z)
c
c                  h    matrice a diagonaliser rangee sous la forme d'un
c                       demi matrice superieure
c                  e    valeurs propres en ordre croissant
c                       rangees dans un vecteur de dimension  n
c
c                  v    matrices des vecteurs propres
c                        de dimension nvec
c
c                  n    dimension de la matrice
c                  nev   nombre de valaurs propres demandees
c                  nvec  nombre de vecteurs  propres demandes
c
      parameter ndetz=500000
      integer ag,rn,alf
      real*8 norm,lambd,l,mult
      dimension h(1),e(1),v(1),b(ndetz),c(ndetz),p(ndetz),q(ndetz),
     1r(ndetz),w(ndetz),vy(ndetz+2),in(ndetz),z(ndetz)
      dimension hs1(ndetz),hs2(ndetz)
c***********************************************************************
c     nstart(i)=indice du premier element de la ligne i d'une demi-matri
c     superieure
      nstart(i)=n*(n+1)/2-(n-i+1)*(n-i+2)/2+1
c***********************************************************************
c  soit i  le nombre de decimales pour la machine utilisee. on prend
c  epsi=10**(-(i+3))   et   epsdg=epsv= 10**(-(i-2))
      i=15
      epsi=10.0**(-(i+3))
      epsv=10.0**(-(i-2))
      epsdg=epsv
c     epsi=1.d-10
c     epsdg=1.0d-5
c     epsv=1.d-5
      itevm=3
      joym=15
c       modm1 = 2**15 - 1
      modm1 = 32767
c dlgmax= puissance de 10 maximale
      dlgmax=36.d0
      rn=modm1 - 2
      alf=259
      agnes=2.0**(-14)
      np1=n+1
      nm1=n-1
      nm2=n-2
      ng=1
      if(nm2)99,9,200
  200 continue
c transfert de la premiere ligne de h dans hs1
      do 2 i=1,n
    2 hs1(i)=h(i)
c
c     debut de la tridiagonalisation
c
c     hs1 contient la ligne ir  a tridiagonaliser
c
      do 801 ir=1,nm2
      c(ir)=hs1(ir)
      ip1=ir+1
      ip2=ir+2
c
c     determination des rotations annulant les elements ir+2,..,n
c     le sinus de l'angle de rotation est stocke a la place de l'element
c     annule
      do 802 i=ip2,n
      t=hs1(i)
      if(dabs(t).lt.1.d-8) then
      hs1(i)=1.
      hs2(i-ip1)=0.
      else
      s=hs1(ip1)
      ss=dsqrt(s*s+t*t)
      hs1(ip1)=ss
      s=s/ss
      t=t/ss
      hs1(i)=s
      hs2(i-ip1)=t
      endif
  802 continue
c
c
c     ish3 : indice-ip1 du premier element de la ligne ip1
c
      ish3=nstart(ip1)-ip1
c
      do 805 k=ip2,n
      kp1=k+1
c
c     ish4 :indice-k du premier element de la ligne k
      ish4=nstart(k)-k
c
c
c     on fait les rotations entre les lignes ir+1 et k k=ir+2,...,n
c
      s=hs1(k)
      t=hs2(k-ip1)
      if(t) 820,821,820
  820 continue
      ih4k=ish4+k
      ih3k=ish3+k
      ih3p=ish3+ip1
      u=h(ih3p)
      l=h(ih3k)
      tt=h(ih4k)
c
      h(ih3p)= (u*s+l*t)*s+(l*s+tt*t)*t
      h(ih3k) =-(u*s+l*t)*t+(l*s+tt*t)*s
      h(ih4k) = (u*t-l*s)*t+ (-l*t+tt*s)*s
c
c
      if(kp1-n) 840,840,841
  840 continue
c
      do 807 ll=kp1,n
      ih3l=ish3+ll
      ih4l=ish4+ll
      u=h(ih3l)
      l=h(ih4l)
      h(ih3l)= u*s+l*t
      h(ih4l)=-u*t+l*s
  807 continue
c
  821 continue
c
      if(kp1-n) 823,823,841
  823 continue
      ih3k=ish3+k
      do 808 ll=kp1,n
      s=hs1(ll)
      t=hs2(ll-ip1)
      if(t) 822,808,822
  822 continue
      ih4l=ish4+ll
      u=h(ih3k)
      l=h(ih4l)
      h(ih3k)= u*s+l*t
      h(ih4l)= -u*t+l*s
  808 continue
  841 continue
c
c     fin de rotation de la ligne k
c
  805 continue
      b(ir)=hs1(ip1)
c
c     on stocke dans hs2 le cosinus de la rotation ou le cosinus plus
c     3. (si le sinus est negatif)
      do 3410 k=ip2,n
      s=hs1(k)
      if(s.ge.0.d0) go to 3410
      hs2(k-ip1)=hs2(k-ip1)+3.0d00
 3410 continue
c
c     on transfere la ligne suivante dans hs1
      do 3050 k=ip1,n
 3050 hs1(k)=h(ish3+k)
c   les cosinus des rotations stockes dans hs2 sont transferes dans h
      ish1=nstart(ir)-ir
      do 3020 k=ip2,n
 3020 h(ish1+k)=hs2(k-ip1)
  801 continue
c
      c(nm1)=h(nstart(nm1))
      b(nm1)=h(nstart(nm1)+1)
      c(n)=h(nstart(n))
      b(n)=0.d0
      go to 850
    9 continue
      iiii=0
      c(nm1)=h(iiii+1)
      b(nm1)=h(iiii+2)
      c(n)=h(iiii+3)
      b(n)=0.d0
c
  850 continue
      norm=dabs(c(1))+dabs(b(1))
      do 10 i=2,n
      t=dabs(c(i))+dabs(b(i))+dabs(b(i-1))
      if(norm-t) 104,10,10
  104 norm=t
   10 continue
      do 11 i=1,n
   11 w(i)=b(i)*b(i)
      k=1
      l=-norm
      mult=norm*epsi
      do 12 i=1,nev
   12 e(i)=norm
   13 u=e(k)
   14 lambd=0.5*(l+u)
      if(lambd-l-mult) 30,30,106
  106 if(u-lambd-mult) 30,30,107
  107 ag=0
      i=1
   16 s=c(i)-lambd
   18 if(dabs(s)-1.d-21) 20,108,108
  108 if(s) 110,110,109
  109 ag=ag+1
  110 i=i+1
      if(i-n) 111,111,22
  111 s=c(i)-lambd-w(i-1)/s
      go to 18
   20 ag=ag+1
      i=i+2
      if(i-n) 16,16,22
   22 if(ag-n+k) 24,24,112
  112 l=lambd
      go to 14
   24 u=lambd
      if(n-ag-nev) 113,113,114
  113 m=n-ag
      go to 115
  114 m=nev
  115 do 26 i=k,m
   26 e(i)=lambd
      go to 14
   30 e(k)=lambd
      k=k+1
      if(k-nev) 13,13,116
  116 if(nvec) 40,999,40
   40 ii=0
      inn=-n
      epsin=norm*epsdg
      do 90 i=1,nvec
      inn=inn+n
      ieps=0
      vy(n+1)=0.
      vy(n+2)=0.
  402 t=e(i)
      ieps=ieps+1
      go to (430,404),ieps
  404 j=i-ii
      if(j) 406,406,407
  406 s=norm
      go to 408
  407 s=t-e(j)
  408 l=t
      do 414 k=i,nev
      u=e(k)
      if(u-l - epsin) 410,416,416
  410 l=u
  414 continue
  416 tt=u-t
      if(tt-s) 418,420,420
  418 s=tt
  420 u=dabs(t)
      if(u-s) 422,424,424
  422 u=s
  424 t=t + u*epsdg*0.1
      ir=0
  430 continue
      do 44 j=1,n
      p(j)=0.
      q(j)=b(j)
      r(j)=c(j)-t
      go to (42,44),ieps
   42 vy(j)=1.0
   44 continue
      do 50 j=1,nm1
      if(dabs(r(j))+dabs(b(j))) 152,152,154
  152 in(j)=0
      w(j)=0.
      r(j)=1.d-30
      go to 50
  154 continue
      if(dabs(r(j))-dabs(b(j)))49,117,117
  117 mult=b(j)/r(j)
      in(j)=0
      go to 48
   49 mult=r(j)/b(j)
      in(j)=1
      r(j)=b(j)
      t=r(j+1)
      r(j+1)=q(j)
      q(j)=t
      p(j)=q(j+1)
      q(j+1)=0.
   48 w(j)=mult
      q(j+1)=q(j+1)-mult*p(j)
      r(j+1)=r(j+1)-mult*q(j)
      if(dabs(r(j))-1.d-30) 118,118,50
  118 r(j)=1.d-30
   50 continue
      if(dabs(r(n))-1.d-30)  119,119,120
  119 r(n)=1.d-30
  120 go to (1202,145),ieps
 1202 if(i-nvec) 155,121,121
  155 if(dabs(e(i+1)-e(i)) - epsin) 53,121,121
  121 if(ii) 122,55,122
  122 kk=1
      go to 54
   53 kk=2
   54 ii=ii+1
      if(ii-1) 123,55,123
  123 joy=0
      ir=0
   51 joy=joy+1
      if(joy-joym) 124,95,95
  124 do 52 j=1,n
      vy(j)=ran(alf)
   52 continue
   55 itev=0
   56 itev=itev+1
      do 66 ji=1,n
      k=n-ji+1
   62 t=vy(k)
      tnum= t-vy(k+1)*q(k)-vy(k+2)*p(k)
      tden=r(k)
      if(dabs(tnum).gt.1.d-30) then
      dltnum=dlog(dabs(tnum))
      else
      dltnum=-30.d0
      end if
      if((dltnum-dlog(dabs(tden))).gt.dlgmax) then
      do 64 j=1,n
   64 vy(j)=vy(j)*1.d-5
      go to 62
      end if
      vy(k)=tnum/tden
   66 continue
      if(itev-1) 145,145,131
  131 s=e(i)
      tt=dabs((c(1)-s)*vy(1) + b(1)*vy(2))
      t=dabs(vy(1))
      do 135 j=2,n
      u =dabs(b(j-1)*vy(j-1) + (c(j)-s)*vy(j) + b(j)*vy(j+1))
      if(tt-u ) 132,133,133
  132 tt=u
  133 u =dabs(vy(j))
      if(t-u ) 134,135,135
  134 t=u
  135 continue
      s=tt/t
      if(s-epsv) 136,136,139
  136 if(itev-2) 69,69,137
  137 write(6,100) i,itev
      go to 69
  139 write(6,101) i,epsv,itev ,tt,t,s
      if(itev-itevm) 145,69,69
  145 do 68 j=1,nm1
      if(in(j)) 144,144,67
  144 vy(j+1)=vy(j+1)-w(j)*vy(j)
      go to 68
   67 t=vy(j)
      vy(j)=vy(j+1)
      vy(j+1)=t-w(j)*vy(j+1)
   68 continue
      go to 56
   69 continue
      t=1./t
      do 98 j=1,n
   98 vy(j)=vy(j)*t
      if(ii-1) 77,147,72
   72 ji=i-ii+1
      m=i-1
      t=0.
      do 70 j=1,n
      u=vy(j)
   70 t=t + u*u
      t=1./t
      ag=1
      ikk=n*(ji-2)
      do 75 k=ji,m
      ikk=ikk+n
      s=0.
      do 73 j=1,n
      jvk=ikk+j
   73 s=s+vy(j)*v(jvk)
      tt=s*z(k)
      do 74 j=1,n
      jvk=ikk+j
   74 vy(j) = vy(j)-tt*v(jvk)
      u=s*tt*t
      t=t/(1.-u)
      if(u-0.75) 75,75,76
   76 ag=2
   75 continue
      go to (160,143),ag
  143 itev=1
      ir=ir+1
      go to (145,1432),ir
 1432 go to (402,51),ieps
  160 if(joy-2) 147,147,146
  146 write(6,100) i,ii,joy,ir
  147 s=0.
      do 162 j=1,n
      jvi=inn+j
      u=vy(j)
      s=s + u*u
  162 v(jvi)=u
      z(i)=1./s
c
c
      go to (78,90),kk
   77 do 777 j=1,n
      jvi=inn+j
  777 v(jvi)=vy(j)
   78 continue
      ii=0
   90 continue
      nbl1=0
      if(nm2)85,85,81
   81 do 84 j=1,nm2
      k=n-j-1
      ks=j+2
      kp=k+1
      m=k+2
      ir=n+1
      is=j+1
      ish1=nstart(k)-k
c
      do 83 kk=m,n
      is=is-1
      ir=ir-1
      t=h(ish1+ir)
      if(t.gt.1.5d0) go to 3460
      if(t.eq.0.d0) go to 83
      s=dsqrt(1.d0-t*t)
      go to 3470
 3460 t=t-3.0d0
      s=-dsqrt(1.0d0-t*t)
 3470 continue
      inn=-n
      do 89 joy=1,nvec
      inn=inn+n
      kpjoy=inn+kp
      irjoy=inn+ir
      u=v(kp joy)
      l=v(ir joy)
      v(kp joy)=u*s-l*t
      v(ir joy)=u*t+l*s
   89 continue
   83 continue
   84 continue
   85 continue
      inn=-n
      do 87 joy=1,nvec
      inn=inn+n
      s=0.
      do 86 j=1,n
      jjoy=inn+j
   86 s=s+v(j joy)*v(j joy)
      s=1./dsqrt(s)
      do 88 j=1,n
      jjoy=inn+j
   88 v(j joy)=v(j joy)*s
   87 continue
  999 return
   95 write(6,102)i,ii,joym
      call exit
  100 format(1h0,10i5)
  101 format(1h0,i4,e16.8,i3,3e18.8)
  102 format(1h0,i4,i3)
      return
 99   continue
      e(1)=h(1)
      v(1)=1.0
      return
      end
c$NOSTANDARD SYSTEM
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
C-----------------------------------------------------------------------
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
      CALL GIVEN(h,EPSI,V,NBAS1,NBAS1,NBAS1)
      do i=1,nbas1
      do j=1,nbas1
      dum=0.d0
      do k=1,nbas1 
      kj=k*(k-1)/2+j
      if(k.lt.j)kj=j*(j-1)/2+k
      dum=dum+s(kj)*v(k,i)
      end do 
      dx=dum-epsi(i)*v(j,i)
      if(dabs(dx).gt.1.d-10)then
      write(6,*)'erreur diagos, i,epsi,dum,v(j,i)dx',i,epsi(i),dum,
     *v(j,i),dx
      end if
      end do
      end do

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
      SUBROUTINE JACVAO(A,N,EIVU,EIVR)                                  
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(N,N),EIVU(N),EIVR(N,N)                                
      IF(N-1) 2,2,1                                                     
    2 EIVR(1,1)=1.0D0                                                   
      EIVU(1)=A(1,1)                                                    
      RETURN                                                            
    1 DO 101 J=1,N                                                      
      DO 100 I=1,N                                                      
  100 EIVR(I,J)=0.0D0                                                   
  101 EIVR(J,J)=1.0D0                                                   
C        FIND THE ABSOLUTELY LARGEST ELEMENT OF A                       
  102 ATOP=0.D0                                                         
      DO 111 I=1,N                                                      
      DO 111 J=I,N                                                      
      IF (ATOP-DABS(A(I,J))) 104,111,111                                
  104 ATOP=DABS(A(I,J))                                                 
  111 CONTINUE                                                          
      IF(ATOP)109,109,113                                               
  109 RETURN                                                            
C        CALCULATE THE STOPPING CRITERION -- DSTOP                      
  113 AVGF=DFLOAT(N*(N-1))*.55D0                                        
      D=0.0D0                                                           
      DO 114 JJ=2,N                                                     
      DO 114 II=2,JJ                                                    
      S=A(II-1,JJ)/ATOP                                                 
  114 D=S*S+D                                                           
      DSTOP=(1.D-09)*D                                                  
C        CALCULATE THE THRESHOLD, THRSH                                 
C                                                                       
      THRSH=DSQRT(D/AVGF)*ATOP                                          
C                                                                       
C        START A SWEEP                                                  
C                                                                       
  115 IFLAG=0                                                           
      DO 130 JCOL=2,N                                                   
      JCOL1=JCOL-1                                                      
      DO 130 IROW=1,JCOL1                                               
      AIJ=A(IROW,JCOL)                                                  
C                                                                       
C        COMPARE THE OFF-DIAGONAL ELEMENT WITH THRSH                    
C                                                                       
      IF (DABS(AIJ)-THRSH) 130,130,117                                  
  117 AII=A(IROW,IROW)                                                  
      AJJ=A(JCOL,JCOL)                                                  
      S=AJJ-AII                                                         
C                                                                       
C        CHECK TO SEE IF THE CHOSEN ROTATION IS LESS THAN THE ROUNDING E
C        IF SO , THEN DO NOT ROTATE.                                    
C                                                                       
      IF (DABS(AIJ)-1.D-09*DABS(S)) 130,130,118                         
  118 IFLAG=1                                                           
C                                                                       
C        IF THE ROTATION IS VERY CLOSE TO 45 DEGREES, SET SIN AND COS   
C     TO 1/(ROOT 2).                                                    
C                                                                       
      IF (1.D-10*DABS(AIJ)-DABS(S)) 116,119,119                         
  119 S=.707106781186548D0                                              
      C=S                                                               
      GO TO 120                                                         
C                                                                       
C        CALCULATION OF SIN AND COS FOR ROTATION THAT IS NOT VERY CLOSE 
C        TO 45 DEGREES                                                  
C                                                                       
  116 T=AIJ/S                                                           
      S=0.25D0/DSQRT(0.25D0+T*T)                                        
C                                                                       
C     COS=C,  SIN=S                                                     
C                                                                       
      C=DSQRT(0.5D0+S)                                                  
      S=2.D0*T*S/C                                                      
C                                                                       
C        CALCULATION OF THE NEW ELEMENTS OF MATRIX A                    
C                                                                       
  120 DO 121 I=1,IROW                                                   
      T=A(I,IROW)                                                       
      U=A(I,JCOL)                                                       
      A(I,IROW)=C*T-S*U                                                 
  121 A(I,JCOL)=S*T+C*U                                                 
      I2=IROW+2                                                         
      IF(I2.GT.JCOL) GO TO 123                                          
      DO 122 I=I2,JCOL                                                  
      T=A(I-1,JCOL)                                                     
      U=A(IROW,I-1)                                                     
      A(I-1,JCOL)=S*U+C*T                                               
  122 A(IROW,I-1)=C*U-S*T                                               
  123 A(JCOL,JCOL)=S*AIJ+C*AJJ                                          
      A(IROW,IROW)=C*A(IROW,IROW)-S*(C*AIJ-S*AJJ)                       
      DO 124 J=JCOL,N                                                   
      T=A(IROW,J)                                                       
      U=A(JCOL,J)                                                       
      A(IROW,J)=C*T-S*U                                                 
  124 A(JCOL,J)=S*T+C*U                                                 
C                                                                       
C        ROTATION COMPLETED.                                            
C                                                                       
  131 DO 125 I=1,N                                                      
      T=EIVR(I,IROW)                                                    
      EIVR(I,IROW)=C*T-EIVR(I,JCOL)*S                                   
  125 EIVR(I,JCOL)=S*T+EIVR(I,JCOL)*C                                   
C                                                                       
C        CALCULATE THE NEW NORM D AND COMPARE WITH DSTOP                
C                                                                       
      S=AIJ/ATOP                                                        
      D=D-S*S                                                           
      IF(D.GE.DSTOP) GO TO 129                                          
C                                                                       
C        RECALCULATE DSTOP AND THRSH TO DISCARD ROUNDING ERRORS         
C                                                                       
      D=0.D0                                                            
      DO 128 JJ=2,N                                                     
      DO 128 II=2,JJ                                                    
      S=A(II-1,JJ)/ATOP                                                 
  128 D=S*S+D                                                           
      DSTOP=(1.D-09)*D                                                  
  129 THRSH=DSQRT(D/AVGF)*ATOP                                          
  130 CONTINUE                                                          
      IF(IFLAG.GT.0) GO TO 115                                          
C                                                                       
C     PLACE EIGENVALUES IN EIVU                                         
      DO 132 J=1,N                                                      
      EIVU(J)=A(J,J)                                                    
  132 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE OVERFL(K)
      K=2
      RETURN
      END
