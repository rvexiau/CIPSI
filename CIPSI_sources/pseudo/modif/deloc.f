       proGRAM POCELO
C ONE CENTER EXPANSION OF LOCAL OPERATORS
      IMPLICIT REAL*16(A-H,O-Z)
      NAMELIST /OCELO/LABEL,ZATO,NVAL,FPSTAR,FBSTAR,FPEND,FBEND,SORTIE
     * ,nce
      COMMON/CASE/ZATO,NVAL,NBID(3),RA,FPSTAR,FBSTAR,FPEND,FBEND
      COMMON/CHAR/LABEL,SORTIE
      DIMENSION A(40),B(1600),C(40) ,XXNOR(10),ISYM(10)
C FPSTAR AMPLITUDE OF THE PSEUDOPOT FOR STARTING THE REPRESENTATION AT R
C FBSTAR * RA :START OF COLLOCATION BASIS
C FPEND  * LSYM/NVAL :END   OF REPRESENTATION
C FBEND  * LSYM/NVAL :END   OF COLLOCATION BASIS
       CHARACTER*8LAB,SORTIE
       CHARACTER*80 LABEL
       SORTIE='MOLCAS'
*1       READ(5,686,END=99)LAB
*686       FORMAT(A8)
*       IF(LAB.EQ.'SYMETRIC')THEN
1        LABEL='        '
        do i=1,80
           label(i:i)=' '
        enddo
        DO10I=1,5
        ZATO=NVAL
        NVAL=0
        FPSTAR=1.D-4
        FBSTAR=4.
        FPEND=.4
        FBEND=.8
        open(unit=31,form='unformatted',status='unknown',
     &       file='psnl_gaussian_bin')
        open(unit=32,form='formatted',status='unknown',
     &       file='psnl_gaussian_txt')
        open(unit=33,form='formatted',status='unknown',
     &       file='MWBPOT_new')     
        READ(5,OCELO,end=99)
        IF(NVAL.EQ.0)GOTO99
          write(6,'(//,a)') Label
          write(7,'(//,a)') Label

10       CALL NOLOC (XXNOR(I),ISYM(I))
50       CONTINUE
        NP=I-1
60       CONTINUE
 1000  FORMAT(A8,I4,F10.6)
 1001  FORMAT(3I3,F12.8)
 1002  FORMAT(4D20.12)
*       ENDIF
      GOTO1
*       ELSEIF(LAB.EQ.'POLARISA')THEN
C        CALL NLOCG
*       WRITE(6,*)' NLOCG INHIBE'
*       ELSE
*      STOP
*        ENDIF
*        GOTO1
99       CONTINUE
           IF(index(sortie,'MOLCAS').ne.0)THEN
             LABEL=' fin de fichier'
             lcst=0
             nspot=0
           WRITE(6,*)' MOLCAS    SAVED ', LABEL,lcst-1

*             rewind 22
             write(22)LABEL,lcst-1,
     *   nspot,0

c nm1=nspot nm2=0
           endif

       close(22)
       close(31)
       close(32)
       close(33)
      STOP
      END
C MAIN ROUTINE FOR EXPANDING A LOCAL POTENTIAL IN A NON LOCAL OPERATOR
C THE TECHNIQUE USED IS SIMILAR TO THE  COLLOCATION TECHNIQUE BUT THE
C NUMERICAL GRID IS REPLACED BY A SET OF GAUSSIAN FUNCTIONS
      SUBROUTINE NOLOC(XNOR,ISYM)
      IMPLICIT REAL*16(A-H,O-Z)
      COMMON/CASE/ZATO,NVAL,NBID(3),RA,FPSTAR,FBSTAR,FPEND,FBEND
      DIMENSION EPSDS(60),Z(60),ZOP(60),RNOR(60),OPNOR(60)
      COMMON/POOL2/VL(4000000)
      COMMON/POOL/SL(3600),TEMP(3600),VKL(3600),INUM(60)
     *, TOPE(3600),TORT(3600),REC(3600)
     *, SW(600),SV(600), SN(600),SC(3600),POIDS(3600),okl(3600)
         real*8 dzop(60),dsc(3600),dsn(600),wapot(50),wcpot(50)
      COMMON/INFO/NPOT(80),APOT(80),CPOT(80),RZ(4),OVZ(3),EX(4),OZ(2),
     * ARAF,AOP,BRAF,BOP,DMNXP,DMNXP2,OPTES,BASTES,NSPOT,NEWBAS,NSTARP
     * ,ISUP,QVER
      LOGICAL*4 NEWBAS,QPAR,QVER,qcst,tabu,extrac
      COMMON/CHAR/LABEL,SORTIE
       CHARACTER*8LAB,SORTIE
       CHARACTER*80 LABEL,labelr
       CHARACTER*8FMT,APSD
      DIMENSION FMT(1),NBTYP(0:5),grid(3),xnort(100) 
      NAMELIST/INPUT/ NBTYP,NPOT,LSYM,
     &APOT,CPOT,DMNXP2,DMNXP,OPTES,BASTES,
     &OVP,OVB,QPAR,APSD,POIX,QVER,QCST,
     * icut,mcut,tabu,extrac,grid,facso

      DATA FMT/'(5D16.8)'/
      DATA APSD/'        '/
       icut=0
        mcut=0
       facso=1.d0
      MAXBAS=40
      MDIM=400
      OVB=.96
      OVP=.96
      OPTES=1.D-6
      BASTES=1.D-6
      DMNXP=1.D-12
      DMNXP2=1.D-12
      NEWBAS=.FALSE.
      QVER=.FALSE.
      tabu=.false.
       extrac=.true.
      grid(1)=7.
       grid(2)=.1
      grid(3)=0.1
      QPAR=.FALSE.
      Qcst=.FALSE.
      POIX=1.
                do i=0,5
       nbtyp(i)=0
          enddo

      READ(5,INPUT,END=999)

         write(6,'(//,a,i2,//)')' moment angulaire :',lsym-1
         write(7,'(a,i2)')' moment angulaire :',lsym-1

      if(nbtyp(lsym).eq.0)then
       write(6,*)' pas de terme semi local dans cette symetrie'
      return
      endif


      write(6,*)' Dans le cas du SO les operateurs Christiansen'
      write(6,*)' doivent etre multiplies par -1. variable facso'
       write(6,*)' ici facso vaut : ',facso
       IF(QVER)THEN
         if(lsym.eq.0)return
	rewind 25
c        do j=1,5
c           read(25) LABEL(1:8),Lsr,NOP,Nopr,(dZOP(I),I=1,NOP),
c     *    (dSC(I),I=1,NOP*nopr),(dSN(I),I=1,Nopr)
c         write(6,*)label,lsr,nop,nopr
c         if(lsr.eq.lsym)goto9
c        enddo
9001      continue
       IF(index(sortie,'MOLCAS').ne.0)THEN
             read(25)LABELr,lcstr,
     *   nspot,lsr,(apot(i),i=1,nspot+lsr),(cpot(i),i=1,nspot+lsr)
         lcstr=lcstr+1
           WRITE(6,*)' MOLCAS  was  SAVED ', LABELr

           if(labelr.eq.label)then
            do j=1,lcstr

           read(25) Lsr,NOP,Nopr,(DSN(I),I=1,Nopr),
     *    (DZOP(I),I=1,Nop),(DSC(I),I=1,NOP*Nopr)
           
         if(lsr.eq.lsym)goto9
        enddo
          else
           do j=1,lcstr
            read(25)
           enddo
          endif
      else
       write(6,*)' non molcas a revoir'
      endif
       go to 9001
 9         continue
         doi=1,nop
           zop(i)=dzop(i)
           sn(i)=dsn(i)
         enddo
         do i=1,nop*nopr
          sc(i)=dsc(i)
           enddo
       write(6,*) ' passage de verification sur lsym=',lsr,'nop',nop
       NEWBAS=.TRUE.
       ENDIF
       WRITE(6,*)'  APSD=',APSD
       IF(QPAR) WRITE (6,*) 'TEST ENGENDREMENT BASE '
      L=LSYM
      ISYM=LSYM
       lcst=0
       do i=1,5
        if(nbtyp(i).ne.0)then
         lcst=i
        endif
       enddo
      
      NSTARP=0
      DO11I=0,LSYM-1
11    NSTARP=NSTARP+NBTYP(I)
* le tableau des param psd est contracte pour avoir seulement la
* symetrie consideree
         j=1
       do i=nstarp+1,nstarp+nbtyp(lsym)
        apot(j)=apot(i)
        cpot(j)=cpot(i)*facso
        npot(j)=npot(i)
        j=j+1
       enddo
       if(qcst.and.lcst.ne.0)then
        DO I=LSYM,lcst-1
        NSTARP=NSTARP+NBTYP(I)
         enddo
        do i=nstarp+1,nstarp+nbtyp(lcst)
        apot(j)=apot(i)
        cpot(j)=cpot(i)*facso
        npot(j)=npot(i)
        j=j+1
        enddo
 
       endif

      NSPOT=j-1
      NSTARP=1
      OVZ(1)=OVB
      OZ(2)=OVP
686      FORMAT(A4)
c        WRITE(6,INPUT)
        write(6,*)' semiloc pot'
        write(6,'(1x,3(i3,2f15.6))')(npot(i),apot(i),cpot(i)
     * ,i=nstarp,nspot)
        if(lsym.ne.0)then
       do i=nstarp,nspot
       if(abs(apot(i)).lt.1.d-6.or.abs(cpot(i)).lt.1.d-6)then
       write(6,*)' terme semi-local nul verifier les donnees'
       stop
       endif
       enddo
        endif
       if(qcst.and.lcst.eq.0)then
       write(6,*)' qcst true        lcst nul verifier les donnees'
       stop
       endif

1     CONTINUE

*terme M1 des psds modeles a la suedoise sert de terme general
*pour les psds semi locaux en -Zeff/r* sum(Ai*exp(-ai r**2))
*terme M2 des psds modeles a la suedoise sert de terme general
*pour les psds semi locaux en -Zeff/r * r * sum (Ai * exp(-ai r**2))
* attention la serie de projecteurs s'arretant a f les projecteurs g
* sont traite comme un terme d'ecrantage pour les psd de stuttgart.
*ceci ne respecte pas la formulation originale de Stuttgart.
*Les differences induites sont negligeables sauf dans le cas de
* mauvais parametrage de la symetrie g. Il est souvent mal aise
* d'extraire des parametres pour les symetries ne comportant pas
* d'electrons de coeur, en particulier avec les programmes atomiques
* numeriques. Les pseudos de HAY et WADT sont peu surs de ce point de
* vue.
* attention il n'y a pas de possibilite pour des termes en 1/r*r
* dans MOLCAS3

        if(LSYM.eq.0.and..not.newbas)THEN
           IF(index(sortie,'MOLCAS').ne.0)THEN
          if(label(1:1).ne.' ')then
          do i=80,2,-1
           label(i:i)=label(i-1:i-1)
          enddo
          label(1:1)=' '
        endif

           WRITE(6,*)' MOLCAS    SAVED ', LABEL,lcst-1
          nspot1=0
          do i=1,nspot
          if(npot(i).eq.-2)nspot1=nspot1+1
          enddo
         if(nspot1.ne.0)then
           write(6,*)' il y a des termes en  r-2  '
     &   ,'dans le pot ecran non supportes dans molcas3'
          stop 20
          endif

          nspot2=0
         do i=1,nspot
          if(npot(i).eq.-1)then
           nspot2=nspot2+1
           wapot(nspot2)=apot(i)
           wcpot(nspot2)=-cpot(i)/float(nval)
           endif
          enddo

          nspot3=nspot2
         do i=1,nspot
          if(npot(i).eq.0)then
          nspot3=nspot3+1
          wapot(nspot3)=apot(i)
          wcpot(nspot3)=-cpot(i)/float(nval)
           endif
         enddo



*             rewind 22
             write(22)LABEL,lcst-1,
     *   nspot2,nspot3-nspot2,(wapot(i),i=1,nspot3),
     &                          (wcpot(i),i=1,nspot3)
              rewind 31
              rewind 32
              rewind 33
              write(31) label,lcst-1
              write(32,1112) label
  
1112          format(a80)
              write(32,*) lcst-1
c nm1=nspot nm2=0
           endif
            return
         else
C DEFINITION OF THE EXPANSION FUNCTIONS FOR THE POTENTIAL
C DEFINITION OF THE COLLOCATION BASIS SET
      CALL SSLOC(NBAS,NBASLD,NOP,NOPR,Z,ZOP,L,L,TORT,TOPE,REC,
     * EPSDS,RNOR,OPNOR)

         endif


       IF(QPAR)RETURN
      LL=L+L


      IF(.NOT.NEWBAS)THEN
        if(tabu)then
        iflag=10 
       RA=grid(2)
       PAS=grid(3)
        ntf=(grid(1)-grid(2))/grid(3)
       DO NT=1,ntf+1
       RV=0.
       Do J=NSTARP,NSPOT
        RV=RV+CPOT(J)*RA**(NPOT(J))*EXP(-RA*RA*APOT(J))
       enddo
        if(abs(rv).lt.1.d-6)iflag=iflag-1
       if(iflag.gt.0)WRITE(6,'(5(a,d15.6))')
     * ' r= ',RA,' pot ',RV-float(nval)/ra+l*(l+1)/(2.*ra*ra),
     *' pspot=',RV,' -Z/r',-float(nval)/ra,' KIN',
     * l*(l+1)/(2.*ra*ra)
        RA=ra+pas
        enddo
        endif

       IJ=0
      DO  I=1,nop
      DO  J=1,I
       ZZ=Zop(I)+Zop(J)
       T=0.
       DO  K=NSTARP,NSPOT
       T=T+PURX(LL+NPOT(K),ZZ+APOT(K))*CPOT(K)
       enddo
        T=T*OPNOR(I)*OPNOR(J)
        tt=OPNOR(I)*OPNOR(J)*PURX(LL-1,ZZ)*float(nval)
        IJ=IJ+1
        okl(IJ)=T
         dsc(ij)=opnor(i)*opnor(j)*purx(ll,zz)
         temp(ij)=t-tt
       enddo
       enddo

c        write(6,*)' test densite / pot'
c       ij=1
c       do j=2,nop
c       ij=ij+j
c       write(6,'(1x,5e15.6)') okl(ij-j),okl(ij),okl(ij-1),dsc(ij-1)
c     * ,okl(ij-1)-(okl(ij)+okl(ij-j))*.5*dsc(ij-1)
c       enddo


c      WRITE(6,*)' ope sl -Z/R ds base de projection'
c      IJ=0
c       DO  I=1,Nop
c       WRITE(6,42) zop(i),(temp(ij+j),J=1,I)
c       IJ=IJ+I
c        enddo

      CALL ORTHO(okl,TEMP,SV,TOpe,Nop,Nopr)
        nvois=nopr
     
       CALL RSP(NVOIS,NVOIS,(NVOIS+1)*NVOIS/2,okl,Sv,1,poids,VL(1),
     * VL(NVOIS+1),IERR)

        DOI=1,NVOIS
        DOJ=1,NOP
         JK=J
         KI=NVOIS*(I-1)
         T=0.D0
         DO K=1,NVOIS
         KI=KI+1
         T=T+poids(KI)*TOPE(JK)
         JK=JK+NOP
          enddo
         temp(NOP*(I-1)+J)=T
          enddo
           enddo
           ij=0
            ijp=0
            nopr=0
         write(6,'(a,6(1x,e11.4))')' vp du pot',
     *     (sv(i),i=1,nvois)

             emax=max(abs(sv(1)),abs(sv(nvois)))+1.
             ncut=0
            do i=1,nvois
          if(abs(sv(i)).lt.dmnxp.and.icut.eq.0)then
           ncut=ncut+1
c          sv(i)=emax
          endif
          enddo
            write(6,*)ncut,' VP eliminees < ',dmnxp
          icut=ncut
             
         do i=1,nvois
c          if(abs(sv(i)).lt.emax)then
            emin=emax
            ii=0
            do j=1,nvois
            if(abs(sv(j)).lt.emin)then
             emin=abs(sv(j))
             ii=j
             endif
             enddo
                  if(ii.eq.0)then
                   write(6,*)'bug dans classement vp du pot'
                    stop
                   endif
            nopr=nopr+1
          do j=1,nop
          ijp=ijp+1
          tope(ijp)=temp(NOP*(Ii-1)+J)
          enddo
          do j=1,nopr
           ij=ij+1
          okl(ij)=0.d0
          enddo
            okl(ij)=sv(ii)
             sv(ii)=emax
c         endif
         enddo
         write(6,'(a,6(1x,e11.4))')' vp du pot',
     *     (okl(i*(i+1)/2),i=1,nvois)
         write(6,*)' base reduite du potentiel dim:',nopr

         write(6,*) ' coupure dans le spectre de l''operateur:'
         write(6,*) icut,' valeurs propres petites en ||'
          write(6,*) mcut,' valeurs propres grandes en ||'

c      CALL ORTHOD    (REC,SV,TORT,NBAS,NBASLD,NOP)

      CALL ORTHOE(REC,SV,TOPE,NOP,NOPR,NBAS)
       do i=1,nopr*(nopr+1)/2
        vkl(i)=okl(i)
       enddo
       call ortho(vkl,temp,sv,rec,nopr,nbas)
c     WRITE(6,*)' ope projete ds base extraction'
c      IJ=0
c      DO  I=1,NBAS
c      WRITE(6,42) z(i),(vKL(IJ+J),J=1,I)
c      IJ=IJ+I
c        enddo
       else
       do i=1,nbas*(nbas+1)/2
        vkl(i)=0.d0
       enddo

       endif

          fnorm=0.
      IJ=0
      DO 40 I=1,NBAS
      DO 40 J=1,I
      ZZ=Z(I)+Z(J)
      T=0.
      DO 35 K=NSTARP,NSPOT
35    T=T+PURX(LL+NPOT(K),ZZ+APOT(K))*CPOT(K)
      T=T*RNOR(I)*RNOR(J)
      IJ=IJ+1
      VKL(IJ)=T-vkl(ij)
       fnorm=fnorm+vkl(ij)**2
40    CONTINUE

       if(newbas)then
      WRITE(6,*)'  SL dans base verification'
        else
      WRITE(6,*)' difference SL NL approx 0 base extraction'
         WRITE(6,510) sqrt(FNORM)
        endif
      IJ=0
      DO 45 I=1,NBAS
      WRITE(6,42) z(i),(VKL(IJ+J),J=1,I)
45    IJ=IJ+I
      IF(.NOT.NEWBAS)THEN
      CALL ORTHO(VKL,TEMP,SV,TORT,NBAS,NBASLD)
        IJ=NBAS*(NBAS+1)/2
      WRITE(6,6)
6     FORMAT('0',20(' * '),/,'0')
      CALL ORTHOD    (REC,SV,TORT,NBAS,NBASLD,NOP)
c      CALL ORTHOE(REC,SV,TOPE,NOP,NOPR,NBASLD)
           ELSE
      CALL ORTHOE(REC,SV,SC,NOP,NOPR,NBASLD)
C       WRITE(6,*)' COEF'
C         WRITE(6,*)(SN(I),I=1,NOP)
C         WRITE(6,*)' CONTRACTION'
C          DO4023J=1,NOPR
C4023       WRITE(6,2024)(SC(J+NOP*(I-1)),I=1,NOP)
C          WRITE(6,*)' BASE OPER  '
C           DO4024J=1,NOPR
C4024       WRITE(6,2024)(REC(J+NOPR*(I-1)),I=1,NBASLD)
C      WRITE(6,584)
         write(6,*)' qualite de la reproduction'
        XNOR=0.
        IK=0
        IJ=0
        DO3000I=1,NBASLD
         JK=0
       DO3100J=1,I
         IJ=IJ+1
       IK0=IK
        TT=0.D0
       DO3200K=1,NOP
       JK=JK+1
       IK0=IK0+1
       TT=TT+REC(IK0)*REC(JK)*SN(K)
3200   CONTINUE
        TEMP(J)=TT-VKL(IJ)
        XNOR=XNOR+(TEMP(J)*TEMP(J))
3100    CONTINUE
         IK=IK+NOP
         xnort(i)=sqrt(xnor)
         if(i.le.nbasld-6)xnorg=xnor
3000       WRITE(6,585)(TEMP(J),J=1,I)
        XNORg=SQRT(XNORg)
       WRITE(6,*)' NORME VERIF (1,NBAS-6) ',XNORg,xnort(nbasld)
       write(6,*) ' norme cumulee en partant de  l''exposant '
       write(6,*)' le + externe : exposant norme '
       write(6,'(4(f11.4,d10.2),/)')(z(i),xnort(i),i=1,nbasld)
       WRITE(7,*)' NORME VERIF (1,NBAS-6) ',XNORg,sqrt(xnort(nbasld))
       write(7,*) ' norme cumulee en partant de  l''exposant '
       write(7,*)' le + externe : exposant norme '
       write(7,'(4(f11.4,d10.2),/)')(z(i),xnort(i),i=1,nbasld)

          RETURN
          ENDIF
42    FORMAT(' exp=',f11.4,' VKL=',8F11.4,(1x,10f11.4))

       if(.not.extrac)return
      NBASS=NBAS
      NBAS=NBASLD
      NOPS=NOP
      JI=0

**************************

      ivois=icut+1
      NVOIS=min(nbasld,NOPR-icut-mcut)
      WRITE(6,53) NVOIS
53    FORMAT(' NB DE TERMES DU PSEUDO NON SINGULIERS ',I4)
       nvois=nvois+icut
      IJJI=0
C       CALL PRTS('VKL',VKL,NBAS)
C       CALL PRT(' REC',REC,NOP)
C LEAST SQUARE FIT OF THE POTENTIAL
C SIGMA       CI,IV <XM|UI> <UIV|XN> = <XM|PSPOT|XN>
      WRITE(6,*)' ISUP',ISUP,QVER,POIX
      MN=0
      DO 54 M=1,NBAS
      DO 54 N=1,M
       MN=MN+1
c      IF(M.LE.ISUP)THEN
c      POIDS(MN)=POIX*POIX
c      ELSE IF(N.LE.ISUP)THEN
c      POIDS(MN)=POIX
c c     ELSE
      POIDS(MN)=1.
c      ENDIF
54     CONTINUE
c      WRITE(6,*)' NVOIS',NVOIS
      IJ=0
      DO 200 I =ivois,NVOIS
      DO 190  IV=ivois,I

C OPERATEUR |I> <IV|
C LIGNE = D/ DCI,IV
      TT=0.D0
      MI=I
      MIV=IV
      MN=0
      DO 60 M=1,NBAS
       REI=REC(MI)
       REV=REC(MIV)
       NIV=IV
       NI=I
      DO 55 N=1,M
       MN=MN+1
       T=REI*REC(NIV)+REV*REC(NI)
       IF(M.EQ.N)T=T*.5D0
       TEMP(MN)=T
       TT=TT+T*VKL(MN)*POIDS(MN)
       NIV=NIV+NOP
55    NI=NI+NOP
       MI=MI+NOP
       MIV=MIV+NOP
60     CONTINUE

      IJ=IJ+1
      SV(IJ)=TT

C     WRITE(6,*) I,IV,IJ
      JI=0
      DO 180  J=ivois,NVOIS
      DO 100  JV=ivois,J
C OPERATEUR |J> <JV|
      JI=JI+1

      TT=0.D0
      MN=0
      MJ=J
      MJV=JV
       DO 80 M=1,NBAS
       REJ=REC(MJ)
        REV=REC(MJV)
      NJ=J
       NJV=JV
      DO 70 N=1,M
       MN=MN+1
      TTT=(REJ*REC(NJV)+REV*REC(NJ))
C05         IF(N.EQ.M)TTT=TTT*.5D0
      TT=TT+TTT*TEMP(MN)*POIDS(MN)
      NJV=NJV+NOP
70    NJ=NJ+NOP
      MJ=MJ+NOP
      MJV=MJV+NOP
80    CONTINUE

      IJJI=IJJI+1
      VL(IJJI)=TT
      IF(JI.NE.IJ) GO TO 100
C           CALL PRTC(' VVVV',VL(IJJI-IJ+1),IJ)
      TT=1.D0/SQRT(TT)
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
200   CONTINUE

C CANONICAL SOLUTION OF THE LINEAR SYSTEM OF EQUATION
      WRITE(6,206) IJ,IJ
      CALL DIAGOS(VL,VL,VL(IJJI+1),TEMP,IJ,IJP,DMNXP)
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

c        ij=0
c        do i=1,nvois
c         ij=ij+i
c        sc(ij)=sc(ij)-.5d0*okl(ij)
c        enddo

c      WRITE(6,'(1x,6d20.12)') (SC(I),I=1,IJ)

      MN=0
      DO 350 I=1,NBAS
      DO 350 J=1,I
      MN=MN+1
350   TEMP(MN)=0.

      IJ=0
      IP=0
      I2=NOP
      I1=1
      DO 500 I=ivois,NVOIS
      DO 490 IV=ivois,I
        IJ=IJ+1
      T=REC(I)*REC(IV)
      MN=0
      MI=I
       MIV=IV
      DO460 M=1,NBAS
       REI=REC(MI)
      REV=REC(MIV)
      NI=I
      NIV=IV
      DO455 N=1,M
      MN=MN+1
      T=REI*REC(NIV)+REV*REC(NI)
      NI=NI+NOP
       NIV=NIV+NOP
      TEMP(MN)=T   *SC(IJ)+TEMP(MN)
455    CONTINUE
      MI=MI+NOP
       MIV=MIV+NOP
460   CONTINUE
490       CONTINUE
500        CONTINUE

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
       IF(NEWBAS)RETURN

997    NTER=NVOIS*(NVOIS+1)/2
585   FORMAT(1X,12D10.2)
       FNORM=FNORM/NBAS
        IJV=0
         ij=0 
       doI=1,ivois-1
       DOJ=1,I
        IJV=IJV+1
        TEMP(IJV)=okl(ijv)
        enddo
        enddo

       DO1950I=ivois,NVOIS
        doj=1,ivois-1
         IJV=IJV+1
        TEMP(IJV)=okl(ijv)
        enddo
       DO1945J=ivois,I
         ij=ij+1
        IJV=IJV+1
1945    TEMP(IJV)=SC(IJ)+okl(ijv)
         TEMP(IJV)=TEMP(IJV)+sc(ij)
1950         CONTINUE
       doI=NVOIS,nopr
       DOJ=1,I
        IJV=IJV+1
        TEMP(IJV)=okl(ijv)
        enddo
        enddo


        nvois=nopr

C        CALL GIVEIS(NVOIS,NVOIS,NVOIS,TEMP,VL,SN,TORT,IERR)
       CALL RSP(NVOIS,NVOIS,(NVOIS+1)*NVOIS/2,TEMP,SN,1,TORT,VL(1),
     * VL(NVOIS+1),IERR)
        DO2000I=1,NVOIS
        DO2000J=1,NOPS
         JK=J
         KI=NVOIS*(I-1)
         T=0.D0
         DO1900K=1,NVOIS
         KI=KI+1
         T=T+TORT(KI)*TOPE(JK)
1900     JK=JK+NOPS
2000      SC(NOPS*(I-1)+J)=T
 

         IF(index(sortie,'MELD').ne.0)then
           WRITE(6,*)' MELD    SAVED ', LABEL
         doi=1,nvois
            dsn(i)=sn(i)
        enddo
              do i=1,nops*nvois
          dsc(i)=sc(i)
           enddo

           doi=1,nops
             dzop(i)=zop(i)
              enddo

            
           WRITE(22) LABEL(1:8),L,NOPS,NVOIS,(DZOP(I),I=1,NOPS),
     *    (DSC(I),I=1,NOPS*NVOIS),(DSN(I),I=1,NVOIS)
           WRITE(23) LABEL(1:8),L,NOPS,NVOIS,(ZOP(I),I=1,NOPS),
     *    (SC(I),I=1,NOPS*NVOIS),(SN(I),I=1,NVOIS)

          ELSE  IF(index(sortie,'MOLCAS').ne.0)THEN
           WRITE(6,*)' MOLCAS    SAVED ', LABEL,nops,nvois
         doi=1,nvois
            dsn(i)=sn(i)
        enddo
              do i=1,nops*nvois
          dsc(i)=sc(i)
           enddo

           doi=1,nops
             dzop(i)=zop(i)
              enddo

           WRITE(22) L,NOPS,NVOIS,(DSN(I),I=1,NVOIS),
     *    (DZOP(I),I=1,NOPS),(DSC(I),I=1,NOPS*NVOIS)
           WRITE(31) L,NOPS,NVOIS,(DSN(I),I=1,NVOIS),
     *    (DZOP(I),I=1,NOPS),(DSC(I),I=1,NOPS*NVOIS)
           WRITE(32,*) L,NOPS,NVOIS
           WRITE(32,1113) (DSN(I),I=1,NVOIS)
           WRITE(32,1113) (DZOP(I),I=1,NOPS)
           DO Itl=0,NOPS-1
              WRITE(32,1113) 
     >        (DSC(I),I=1+Itl,NVOIS*NOPS,NOPS)
           ENDDO
1113       format(3d21.13)

           if(l.eq.1) write(33,4000) label,lcst-1
4000  FORMAT(A80,I4)  
           WRITE(33,4001) label,L,NOPS,NVOIS
           WRITE(33,4002) (DSN(I),I=1,NVOIS)
           WRITE(33,4003) (DZOP(I),I=1,NOPS)
           write(33,4004)
           write(33,4005) (DSC(I),I=1,NVOIS*NOPS)
c           DO Itl=0,NOPS-1
c              WRITE(33,4005) 
c     >        (DSC(I),I=1+Itl,NVOIS*NOPS,NOPS)
c           ENDDO
           
4001  FORMAT('label =',a80, ' sym.max =', i3,/,' nombre d''exposants =',
     * i3,' nombre de combinaisons lineairement independantes =',i3)
4002  FORMAT(' coefficients des operateurs',/, (4D20.12))
4003  FORMAT(' exposants des operateurs',/, (4D20.12))
4004  FORMAT(' coefficients des combinaisons lineaires ') 
4005  format(4d20.12)

          ELSE
          WRITE(6,*)SORTIE, ' FORMAT DE SORTIE NON DEFINI '
          ENDIF
cc        WRITE(6,*)' COEF'
c           WRITE(6,'(1x,6d20.8)')(SN(I),I=1,NOP)
c          WRITE(6,*)(SN(I),I=1,NOP)
c          WRITE(6,*)' CONTRACTION'
c           DO2023J=1,NOPS
c2023       WRITE(6,2024)(SC(J+NOP*(I-1)),I=1,NOP)
2024         FORMAT(1X,10F12.6)
       NEWBAS=.TRUE.
      GO TO 1
999    CONTINUE
      STOP
      END
      SUBROUTINE SSLOC
     &(NBASA,NBASRA,NOPA,NOPRA,ZA,ZOPA,LA,LBASA,TORTA,TOPA,RECA,EPSDSA,
     &RNORA,OPNORA)
      PARAMETER(MAXOP=60)
C DETERMINATION OF THE BASIS SET
      IMPLICIT REAL*16(A-H,O-Z)
      COMMON/CASE/ZATO,NVAL,NBID(3),RA,FPSTAR,FBSTAR,FPEND,FBEND
      DIMENSION      ZA(60),ZOPA(60),RNORA(60),OPNORA(60)
      DIMENSION TOPA(3600),TORTA(3600),RECA(3600)
      DIMENSION EPSDSA(60)
      COMMON/POOL2/VL(4000000)
      COMMON/POOL/SL(3600),TEMP(3600),VKL(3600),INUM(60)
      COMMON/INFO/NPOT(80),APOT(80),CPOT(80),RZ(4),OVZ(3),EX(4),OZ(2),
     * ARAF,AOP,BRAF,BOP,DMNXP,DMNXP2,OPTES,BASTES,NSPOT,NEWBAS,NSTARP
     * ,ISUP,QVER
      LOGICAL*4 NEWBAS
C DETERMINATION DE LA BASE DE REPRESENTATION

       WRITE(6,*)
       WRITE(6,*)' *********** POTENTIEL *********** '
       IF(.NOT.NEWBAS)WRITE(6,*)' TABULATION (GROSSIERE)  '
       RA=15.
       PAS=-.5
       DO530 NT=1,30
       RV=0.
       DO510 J=NSTARP,NSPOT
510       RV=RV+CPOT(J)*RA**(NPOT(J))*EXP(-RA*RA*APOT(J))
       IF(.NOT.NEWBAS)WRITE(6,*) RA,RV
        RA=RA+PAS
530     CONTINUE
       RA=15.
       PAS=-.1
       DO430 NT=1,150
       RV=0.
       DO410 J=NSTARP,NSPOT
410       RV=RV+CPOT(J)*RA**(NPOT(J))*EXP(-RA*RA*APOT(J))
        IF(ABS(RV).GT.FPSTAR)THEN
        IF(ABS(RV)-0.1*FPSTAR.LT..01*FPSTAR)GO TO435
        RA=RA-PAS
        PAS=PAS*.1
        ENDIF
        RA=RA+PAS
430     CONTINUE
435      CONTINUE
       RA=RA+RA
      I=ANIPUR(1)
C     IF(OZ(1).EQ.0.D0) GO TO  2160
      IF(NOPA.EQ.0) NOPA=NBASA
      BOP=OZ(2)
      RB=LA*FPEND/FLOAT(NVAL)
      RAP=RA
       IF(.NOT.NEWBAS)THEN
      WRITE(6,*) ' DEPART PS ',RA,' FIN ',RB,' * ',FPSTAR,FPEND,LA,NVAL
      CALL GENB(NOPA,ZOPA,OPNORA,LA,LBASA,RA,RB,BOP,TOPA,0.D0,ISUP)
      WRITE(6,31) SL(2),SL(4),SL(7)
759     FORMAT('0 ORTHONORM DES  OPERATEURS')
      CALL DIAGOS(SL,TOPA,TEMP,EPSDSA,NOPA,NOPRA,DMNXP2)
761   FORMAT(1X,10F12.7)
721        CONTINUE
      WRITE(6,6002)(NPOT(I),APOT(I),CPOT(I),I=NSTARP,NSPOT)
6002   FORMAT(' PSD S L',5(1X,I2,F10.5,F10.5))
      WRITE(6,6004)(ZOPA(I),I=1,NOPA)
6004  FORMAT('OPE N. L.',10F12.6)
      ELSE
       DO450J=1,NOPA
450   OPNORA(J)=1.D0/SQRT(PURX(LA+LA,ZOPA(J)+ZOPA(J)))
      ENDIF
C DETERMINATION DE LA BASE D EXTRACTION
       WRITE(6,*)
       WRITE(6,*)' *********** BASE D EXTRACTION *********** '
      RB=LBASA*FBEND/FLOAT(NVAL)
      EX(4)=RA*FBSTAR
      IF(NEWBAS)EX(4)=2.*EX(4)
      IF(NEWBAS)OVZ(1)=(1.-OVZ(1))*.2+OVZ(1)
      IF(NEWBAS)WRITE(6,*)' NEW OVERLAP ',OVZ(1)
      WRITE(6,*) ' DEPART BA ',EX(4),' FIN ',RB,' * ',FBSTAR,FBEND,LBASA
      CALL GENB(NBASA,ZA,RNORA,LBASA,LBASA,EX(4),RB,OVZ(1),
     *TORTA,RAP,ISUP)
      CALL DIAGOS(SL,TORTA,TEMP,EPSDSA,NBASA,NBASRA,DMNXP2)
c base operateur
      IJ=0
      DO 30 I=1,NBASA
      ZZ=ZA(I)
      DO 30 J=1,NOPA
      IJ=IJ+1
30    RECA(IJ)=PURX(LBASA+LBASA,ZZ+ZOPA(J))*RNORA(I)*OPNORA(J)
      IF(NEWBAS)WRITE(6,*)' VERIFICATION '
31      FORMAT(' RECAOUVREMENT ENTRE VOISINS:',3F15.8)
6001   FORMAT(' SYMETRIE NL',I3,'SYMETRIE BASE',I3,' NB GAUSSIENN
     &ES:',I3,'NB OPE N. L.',I3)
      WRITE(6,6003)(ZA(I),I=1,NBASA)
6003  FORMAT(' BASE:',10F12.6)
        RETURN
         END
      SUBROUTINE GENB
     &(NOPA,ZOPA,OPNORA,LA,LBASA,RA,RB,BOP,TOPA,RP,ISUP)
      PARAMETER(MAXOP=60)
C DETERMINATION OF THE BASIS SET
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION       ZOPA(60),OPNORA(60)
      DIMENSION TOPA(3600)
      COMMON/POOL2/VL(4000000)
      COMMON/POOL/SL(3600),TEMP(3600),VKL(3600),INUM(60)
C DETERMINATION DE LA BASE DE REPRESENTATION
      AOP=LA/(2.D0*RA*RA)
      COP=LA/(2.D0*RB*RB)
      DOP=0.
      IF(RP.GT.1.D-6)DOP=LA/(2.D0*RP*RP)
      T=BOP**(2.D0/(DFLOAT(LA)+0.5D0))
      TOP=(2.D0-T+SQRT((2.D0-T)**2-T))/T
      WRITE(6,*)' genb TOP AOP',TOP,AOP,LA,RA,BOP,T
2160     CONTINUE
       IF(AOP.EQ.0.) STOP     ' NO BASIS FOR PSPOT OPERATOR'
        DO110 I=1,MAXOP
        DO109 J=1,MAXOP
109     TOPA(J+MAXOP*(I-1))=0.
110     TOPA(I+MAXOP*(I-1))=1.
      ZOPA(1)=AOP
       I=1
      OPNORA(I)=1.D0/SQRT(PURX(LBASA+LBASA,ZOPA(I)+ZOPA(I)))
       ISUP=0
      SL(I)=1.

      DO 121 I=2,MAXOP
       IF(ZOPA(I-1).GT.COP)GOTO125
       IF(ZOPA(I-1).LT.DOP)ISUP=I

      DO 119 ITE=1,20
        DO111 J=1,I
111     TOPA(J+MAXOP*(I-1))=0.
        TOPA(I+MAXOP*(I-1))=1.
      ZOPA(I)=ZOPA(I-1)*TOP
c      WRITE(6,*) ' EXP ',ZOPA(I),I,ITE
      OPNORA(I)=1.D0/SQRT(PURX(LBASA+LBASA,ZOPA(I)+ZOPA(I)))
      IJ=(I-1)*I/2
      DO 710 J=1,I
      IJ=IJ+1
710    SL(IJ)=OPNORA(I)*OPNORA(J)*PURX(LBASA+LBASA,  ZOPA(I)+ZOPA(J))
      REC     = SPROJ(TOPA,SL,MAXOP,1,I-1,I)
c      WRITE(6,*)' REC',REC
       IF(ABS(REC-BOP).LT.1.D-3)GOTO120
C       IF(ITE.EQ.1)THEN
       IF(REC.LT.BOP)ROP=TOP*0.9
       IF(REC.GT.BOP)ROP=TOP*1.05
C      ELSE
C       ROP=TOP*BOP/REC(ITE)*TOP/TOPM*REC(ITE-1)/REC(ITE)
C      ENDIF
        TOPM=TOP
        TOP=ROP
119    CONTINUE

120    CONTINUE
        DO122 J=1,I
122     TOPA(J+MAXOP*(I-1))=0.
        TOPA(I+MAXOP*(I-1))=1.
      REC     = SPROJ(TOPA,SL,MAXOP,1,I-1,I)
      CALL OPROJ(TOPA,SL,MAXOP,I,1,I-1,I,1.D-21)
121    CONTINUE
125    NOPA=I-1
       WRITE(6,*)' NOMBRE D EXP GENERES ',NOPA
C      IF(ISUP.EQ.0)RETURN
C     DO 421 I=1,ISUP
C     ZOPA(NOPA+I)=SQRT(ZOPA(I)+ZOPA(I+1))
C     OPNORA(NOPA+I)=1.D0/SQRT(PURX(LBASA+LBASA,ZOPA(NOPA+I)*2))
C     IJ=(I+NOPA-1)*(I+NOPA)/2
C     DO 410 J=1,I+NOPA
C     IJ=IJ+1
C410    SL(IJ)=OPNORA(I+NOPA)*OPNORA(J)
C     *        *PURX(LBASA+LBASA,  ZOPA(I+NOPA)+ZOPA(J) )
C421    CONTINUE
C       NOPA=NOPA+ISUP
C       WRITE(6,*)' NOMBRE D EXP GENERES ',NOPA
        RETURN
         END
      SUBROUTINE TTP(S1,S2,N,M)
      IMPLICIT REAL*16 (A-H,O-Z)
      DIMENSION S1(N),S2(M,N)
       IJ=0
      DO 50 I=1,N
      DO 50 J=1,I
      IJ=IJ+1
      S2(I,J)=S1(IJ)
50    S2(J,I)=S1(IJ)
      RETURN
      ENTRY TPT(S1,S2,N,M)
       IJ=0
      DO 60 I=1,N
      S2(I,I)=SQRT(S2(I,I))
      T=1.D0/S2(I,I)
      DO 60 J=1,I
      IJ=IJ+1
60    S1(IJ)=S2(I,J)      *T
      RETURN
      END
      FUNCTION SPROJ(CV,S,N,N1,N2,N3)
      IMPLICIT REAL*16 (A-H,O-Z)
      DIMENSION CV(N,N),S(N*N)
       DIMENSION V(100),U(100)
       I=N3
      REC=0.D0
      DO 50 KP=1,N3
      T=0.D0
      DO 40 KQ=1,N3
      IF(KP.GE.KQ) KPKQ=KP*(KP-1)/2+KQ
      IF(KP.LT.KQ) KPKQ=KQ*(KQ-1)/2+KP
40    T=T+CV(KQ,I)*S(KPKQ)
      REC=REC+T*CV(KP,I)
50    V(KP)=T
C      WRITE(6,*)' OPROJ',(V(KP),KP=1,N3),REC
      DO 70 J=N1,N2
      T=0.D0
      DO 60 KP=1,N3
60    T=T+V(KP)*CV(KP,J)
      REC=REC-T*T
70    U(J)=T
C      WRITE(6,*)' OPROJ',(U(KP),KP=1,N3),REC
      REC=REC+T*T
       IF(REC.GT.1.D-12) THEN
      REC=1.D0/SQRT(REC)
       ELSE
       WRITE(6,*)' ELIMINE',N3,REC
       NC=NC-1
        RETURN
       ENDIF
      DO 95 KP=1,N3
      DO 90 J=N1,N2-1
        T=U(J)
90    CV(KP,I)=CV(KP,I)-T*CV(KP,J)
95    CV(KP,I)=CV(KP,I)*REC
100   CONTINUE
      DO150 KP=1,N3
C     WRITE(6,'(I3,5F8.4)')KP,(CV(KP,J),J=1,N3)
      T=0.D0
      DO140 KQ=1,N3
      IF(KP.GE.KQ) KPKQ=KP*(KP-1)/2+KQ
      IF(KP.LT.KQ) KPKQ=KQ*(KQ-1)/2+KP
140    T=T+CV(KQ,I)*S(KPKQ)
150    V(KP)=T
C      WRITE(6,*)' OPROJ',(V(KP),KP=1,N3)
      DO190 J=N1,N3
      T=0.D0
      DO 160 KP=1,N3
160    T=T+V(KP)*CV(KP,J)
190   U(J)=T
      SPROJ=U(N2)
C      WRITE(6,*)' OPROJ',(U(KP),KP=1,N3)
       REC=U(N3)*U(N3)-U(N2)*U(N2)
       IF(REC.GT.1.D-12) THEN
      REC=1.D0/SQRT(REC)
       ELSE
       WRITE(6,*)' ELIMINE',N3,REC
       NC=NC-1
        RETURN
       ENDIF
      DO295 KP=1,N3
            J=N2
        T=U(J)
      CV(KP,I)=CV(KP,I)-T*CV(KP,J)
295    CV(KP,I)=CV(KP,I)*REC
      RETURN
      END
      SUBROUTINE OPROJ(CV,S,ND,N,NA,NB,NC,PREC)
      IMPLICIT REAL*16 (A-H,O-Z)
      DIMENSION CV(ND,ND),S(N*N)
       DIMENSION V(100),U(100)
       N1=NA
      N2=NB
      N3=NC
       I=N3
      DO100ITER=1,3
      REC=0.D0
      DO 50 KP=1,N
      T=0.D0
      DO 40 KQ=1,N
      IF(KP.GE.KQ) KPKQ=KP*(KP-1)/2+KQ
      IF(KP.LT.KQ) KPKQ=KQ*(KQ-1)/2+KP
40    T=T+CV(KQ,I)*S(KPKQ)
      REC=REC+T*CV(KP,I)
50    V(KP)=T
      DO 70 J=N1,N2
      T=0.D0
      DO 60 KP=1,N
60    T=T+V(KP)*CV(KP,J)
      REC=REC-T*T
70    U(J)=T
       IF(REC.GT.PREC) THEN
      REC=1.D0/SQRT(REC)
       ELSE
       WRITE(6,*)' ELIMINE',N3,REC
       NC=NC-1
        RETURN
       ENDIF
      DO 95 KP=1,N
      DO 90 J=N1,N2
        T=U(J)
90    CV(KP,I)=CV(KP,I)-T*CV(KP,J)
95    CV(KP,I)=CV(KP,I)*REC
100   CONTINUE
      RETURN
      END
      FUNCTION PURX(M,ZZ)
C  OVERLAP FOR GAUSSIAN ORBITALS
      IMPLICIT REAL*16 (A-H,O-Z)
      common/ff/FACTO(20),FATT(20)
      DATA PI/3.14159265358979/
      C=1.D0/ZZ
20    M2=M/2
      IF((M2*2).NE.M) GO TO 60
      IF(M.EQ.0) GO TO 50
      PURX=SQRT(PI*C)*FACTO(M)*.5D0*(.5D0*C)**M2
      GO TO 99
50    PURX=.5D0*SQRT(PI*C)
      GO TO 99
60    M2=M2+1
      PURX=FATT(M2)*.5D0*C**M2
99    RETURN
      end
      function RPURX(M,ZZ,RC)
      IMPLICIT REAL*16 (A-H,O-Z)
      common/ff/FACTO(20),FATT(20)
      DATA PI/3.14159265358979/
            C=1.D0/ZZ
            CC=C*.5D0
             M2=(M+1)/2
             RI=0.D0
              RC2=RC*RC
              TU=RC2*ZZ
             IF(TU.LT.72.D0)THEN
            TT=CC
            DO 121I=1,M2-1
            RI=RI+2.D0
            TT=(TT*RI+RC2)*CC
121         CONTINUE
             TT=TT*EXP(-TU)
              ELSE
               TT=0.D0
              ENDIF
            RPURX=FATT(M2)*.5D0*C**M2-TT
             RETURN
                end
         function anipur(j)
       IMPLICIT REAL*16 (A-H,O-Z)
       common/ff/FACTO(20),FATT(20)
       FACTO(1)=1.D0
       FACTO(2)=1.D0
       FATT(1)=1.D0
       FATT(2)=1.D0
        RI=1.D0
       DO 200 I=2,19
         RI=RI+1.D0
       FACTO(I+1)=RI*FACTO(I-1)
 200    FATT(I+1)=RI*FATT(I)
       ANIPUR=RI
       RETURN
       END
      SUBROUTINE DIAGOS(S,DT,V,EPS,NBAS1,NBASLD,DMNXP)
C LOWDIN ORTHONORMALIZATION ROUTINE (S-1/2)
C USE OF PSEUDOINVERSE
C WITH ELIMINATION OF LINEAR DEPENDENCIES (THRESHOLD DMNXP)
      IMPLICITREAL*16(A-H,O-P,R-Z),LOGICAL*4(Q)
      PARAMETER(NCPM=1000)
      DIMENSION DT(NBAS1,NBAS1),S(NCPM),V(NBAS1,NBAS1),EPS(NCPM)
        DIMENSION C(10*NCPM)
C      CALL GIVEIS(NBAS1,NBAS1,NBAS1,S,C,EPS,V,IERR)
          MD=NBAS1*(NBAS1+1)/2
           IKP=1
       CALL RSP(NBAS1,NBAS1,MD,S,EPS,IKP,V,C,
     * C(NBAS1+1),IERR)
        WRITE(6,*)' IERR',IERR
      WRITE(6,63) (EPS(I),I=1,NBAS1)
63    FORMAT(' VP DE S',10E12.5)
64     FORMAT(1X,12F10.5)
      K=0
        DO 50 II=1,NBAS1
      I=II
      IF(EPS(I).LT.DMNXP) GO TO 50
      TERM=1.D0/SQRT(EPS(I))
      K=K+1
      EPS(K)=EPS(I)
C      WRITE(6,*)
C      WRITE(6,64)(V(J,I),J=1,NBAS1)
      DO 40 J=1,NBAS1
40    DT(J,K)=V(J,I)*TERM
50    CONTINUE
      NBASLD=K
      WRITE(6,60) NBASLD
60    FORMAT(' BASE REDUITE DIMENSION:',I4)
      RETURN
      END
      SUBROUTINE ORTHOE(V,EPSI,W,NOP,NR,NB)
      IMPLICIT REAL*16 (A-H,O-Z)
      DIMENSION V(1000)     ,EPSI(1000),W(NOP,NR)
      JK=0
       JL=0
       DO 80 I=1,NB
        DO 70 J=1,NOP
      JK=JK+1
70      EPSI(J)=V(JK)
      DO 80 J=1,NR
      T=0.D0
      DO 60 K=1,NOP
60      T=T+W(K,J)*EPSI(K)
      JL=JL+1
80       V(JL)=T
      RETURN
      END
      SUBROUTINE ORTHOD(V,EPSI,W,NB,NR,NOP)
C V(I,J)= V(I,K)*W(K,J)
      IMPLICIT REAL*16 (A-H,O-P,R-Z),LOGICAL*4(Q)
      DIMENSION V(NOP,NB),W(NB,NB),EPSI(NB)
      DO 50 I=1,NOP
      DO 40 J=1,NB
40    EPSI(J)=V(I,J)
      DO 50 J=1,NR
      T=0.D0
      DO 30 K=1,NB
30    T=T+W(K,J)*EPSI(K)
50    V(I,J)=T
      RETURN
      END
      SUBROUTINE ORTHO(VPHI,V,EPSI,W,NB,NR)
      IMPLICIT REAL*16 (A-H,O-P,R-Z),LOGICAL*4(Q)
      DIMENSION W(NB,NB),V(NB,NB),VPHI(1000),EPSI(1000)
      IJ=0
      DO 20 I=1,NB
      DO 20 J=1,I
      IJ=IJ+1
      V(I,J)=VPHI(IJ)
20    V(J,I)=VPHI(IJ)
      KL=0
      DO 50 K=1,NR
      DO 35 J=1,NB
      T=0.D0
      DO 30 I=1,NB
30    T=T+W(I,K)*V(I,J)
35    EPSI(J)=T
      DO 50 L=1,K
      T=0.D0
      DO 40 J=1,NB
40    T=T+EPSI(J)*W(J,L)
      KL=KL+1
50    VPHI(KL)=T
      RETURN
       END
       SUBROUTINE RSP(N1,ND,NZ,WORK,EW,IK,CW,C1,C2,IER)
       IMPLICIT REAL*16 (A-H,O-Z)
       DIMENSION WORK(1),EW(1),CW(1),C1(1),C2(1)
        CALL GIVEIS(ND,ND,ND,WORK,C1,EW,CW,IERR)
        RETURN
        END
*      DEBUG SUBCHK
*     ENDDEBUG
CEDIT  GIVEIS
C*****************************************************
      SUBROUTINE GIVEIS(N,NVECT,NV,A,B,ROOT,VECT,IERR)
C*****************************************************
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION A(1),B(N,9),ROOT(N),VECT(NV,1)
C
C  FIND ALL EIGENVALUES AND SOME EIGENVECTORS OF A REAL SYMMETRIC MATRIX
C               AUTHORS: C. MOLER AND D. SPANGLER, N.R.C.C., 4/1/79.
C
C     INPUT..
C     N     = ORDER OF MATRIX .
C     NVECT = NUMBER OF VECTORS DESIRED.  0 .LE. NVECT .LE. N .
C     NV    = ROW DIMENSION OF VECT .
C     A     = INPUT MATRIX, COLUMNS OF THE UPPER TRIANGLE PACKED INTO
C             LINEAR ARRAY OF DIMENSION N*(N+1)/2 .
C     B     = SCRATCH ARRAY, 9*N ELEMENTS (NOTE THIS IS MORE THAN
C             PREVIOUS VERSIONS OF GIVENS.)
C
C     OUTPUT..
C     A       DESTORYED .
C     ROOT  = ALL EIGENVALUES, ROOT(1) .LE. ... .LE. ROOT(N) .
C             (FOR OTHER ORDERINGS, SEE BELOW.)
C     VECT  = EIGENVECTORS FOR ROOT(1),..., ROOT(NVECT) .
C     IERR  = 0 IF NO ERROR DETECTED,
C           =  K IF ITERATION FOR K-TH EIGENVALUE FAILED,
C           = -K IF ITERATION FOR K-TH EIGENVECTOR FAILED.
C
C         IF TINVIT FAILS TO CONVERGE, TQL2 IS CALLED
C
C     SEE EISPACK USERS GUIDE, B. T. SMITH ET AL, SPRINGER-VERLAG
C     LECTURE NOTES IN COMPUTER SCIENCE, VOL. 6, 2-ND EDITION, 1976 .
C     IMTQLV AND TINVTB HAVE INTERNAL MACHINE DEPENDENT CONSTANTS.
C
      PARAMETER( ZERO=0.0D0, ONE=1.0D0)
CMULTI      SAVE ZERO,ONE
C
      CALL TRED3B(N,(N*N+N)/2,A,B(1,1),B(1,2),B(1,3))
      CALL IMTQLV(N,B(1,1),B(1,2),B(1,3),ROOT,B(1,9),IERR,B(1,4))
      IF (IERR .NE. 0) RETURN
C
      IF (NVECT .LE. 0) RETURN
      CALL TINVTB(NV,N,B(1,1),B(1,2),B(1,3),NVECT,ROOT,B(1,9),VECT,IERR,
     +     B(1,4),B(1,5),B(1,6),B(1,7),B(1,8))
      IF (IERR .EQ. 0) GO TO 160
C
C     IF INVERSE ITERATION GIVES AN ERROR IN DETERMINING THE VECTORS
C     TRY THE QL ALGORITHM IF ALL THE EIGENVECTORS ARE DESIRED.
C
      IF (NVECT .NE. N) RETURN
      DO 120 I = 1, NVECT
      DO 100 J = 1, N
      VECT(I,J) = ZERO
  100 CONTINUE
      VECT(I,I) = ONE
  120 CONTINUE
      CALL TQL2 (NV,N,B(1,1),B(1,2),VECT,IERR)
      DO 140 I = 1, NVECT
      ROOT(I) = B(I,1)
  140 CONTINUE
      IF (IERR .NE. 0) RETURN
  160 CALL TRBK3B(NV,N,(N*N+N)/2,A,NVECT,VECT)
      RETURN
      END
*      DEBUG SUBCHK
*     ENDDEBUG
CEDIT  DAXPY
C*******************************************
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C*******************************************
      IMPLICIT REAL*16(A-H,O-Z)
C
C     OVERWRITE DY(*) WITH DA*DX + DY.4
C
      DIMENSION DX(1),DY(1)
      PARAMETER(ZERO=0.D0)
CMULTI       SAVE ZERO
      IF(N.LE.0.OR.DA.EQ.ZERO) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = N - (N/4)*4
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DA*DX(I) + DY(I)
   70     CONTINUE
      RETURN
      END
*      DEBUG SUBCHK
*     ENDDEBUG
CEDIT eDOT
C******************************************************
CTL    DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
       REAL*16          FUNCTION DDOT(N,DX,INCX,DY,INCY)
CSINGLE    FUNCTION DDOT(N,DX,INCX,DY,INCY)
CCDC  FUNCTION DDOT(N,DX,INCX,DY,INCY)
C******************************************************
      IMPLICIT REAL*16(A-H,O-Z)
C
C     RETURNS THE DOT PRODUCT OF DX AND DY
C
      DIMENSION DX(1),DY(1)
      PARAMETER(ZERO=0.D0)
CMULTI        SAVE ZERO
      DDOT = ZERO
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + DX(IX)*DY(IY)
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = N - (N/5)*5
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     $   DX(I + 2)*DY(I + 2)+ DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE. 1
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + DX(I)*DY(I)
70    CONTINUE
      RETURN
      END
*      DEBUG SUBCHK
*     ENDDEBUG
CEDIT  TRED3B
C*************************************
      SUBROUTINE TRED3B(N,NV,A,D,E,E2)
C*************************************
      IMPLICIT REAL*16(A-H,O-Z)
C
      DIMENSION A(NV),D(N),E(N),E2(N)
      PARAMETER(ZERO=0.D0)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
C
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
C
C     ON OUTPUT-
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
C          TRANSFORMATIONS USED IN THE REDUCTION,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DO 300 II = 1, N
      I = N + 1 - II
      L = I - 1
      IZ = (I * L) / 2
      H = ZERO
      SCALE = ZERO
      IF (L .LT. 1) GO TO 120
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
      DO 100 K = 1, L
      IZ = IZ + 1
      D(K) = A(IZ)
      SCALE = SCALE + ABS(D(K))
  100 CONTINUE
C
      IF (SCALE .NE. ZERO) GO TO 140
  120 E(I) = ZERO
      E2(I) = ZERO
      GO TO 280
C
  140 DO 160 K = 1, L
      D(K) = D(K) / SCALE
      H = H + D(K) * D(K)
  160 CONTINUE
C
      E2(I) = SCALE * SCALE * H
      F = D(L)
      G = -SIGN(SQRT(H),F)
      E(I) = SCALE * G
      H = H - F * G
      D(L) = F - G
      A(IZ) = SCALE * D(L)
      IF (L .EQ. 1) GO TO 280
      F = ZERO
C
      JK = 1
      DO 220 J = 1, L
      JM1 = J - 1
      DT = D(J)
      G = ZERO
C     ********** FORM ELEMENT OF A*U **********
      IF (JM1 .EQ. 0) GO TO 200
      DO 180 K = 1, JM1
      E(K) = E(K) + DT * A(JK)
      G = G + D(K) * A(JK)
      JK = JK + 1
  180 CONTINUE
  200 E(J) = G + A(JK) * DT
      JK = JK + 1
C     ********** FORM ELEMENT OF P **********
  220 CONTINUE
      F = ZERO
      DO 240 J = 1, L
      E(J) = E(J) / H
      F = F + E(J) * D(J)
  240 CONTINUE
C
      HH = F / (H + H)
      JK = 0
C     ********** FORM REDUCED A **********
      DO 260 J = 1, L
      F = D(J)
      G = E(J) - HH * F
      E(J) = G
C
      DO 260 K = 1, J
      JK = JK + 1
      A(JK) = A(JK) - F * E(K) - G * D(K)
  260 CONTINUE
C
  280 D(I) = A(IZ+1)
      A(IZ+1) = SCALE * SQRT(H)
  300 CONTINUE
C
      RETURN
      END
*      DEBUG SUBCHK
*     ENDDEBUG
CEDIT  IMTQLV
C***********************************************
      SUBROUTINE IMTQLV(N,D,E,E2,W,IND,IERR,RV1)
C***********************************************
      IMPLICIT REAL*16(A-H,O-Z)
      INTEGER TAG
      DIMENSION D(N),E(N),E2(N),W(N),RV1(N),IND(N)
      PARAMETER( ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
C     THIS SUBROUTINE IS A VARIANT OF  IMTQL1  WHICH IS A TRANSLATION OF
C     ALGOL PROCEDURE IMTQL1, NUM. MATH. 12, 377-383(1968) BY MARTIN AND
C     WILKINSON, AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC TRIDIAGONAL
C     MATRIX BY THE IMPLICIT QL METHOD AND ASSOCIATES WITH THEM
C     THEIR CORRESPONDING SUBMATRIX INDICES.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C     ON OUTPUT-
C
C        D AND E ARE UNALTERED,
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO,
C
C        W CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES,
C
C        IND CONTAINS THE SUBMATRIX INDICES ASSOCIATED WITH THE
C          CORRESPONDING EIGENVALUES IN W -- 1 FOR EIGENVALUES
C          BELONGING TO THE FIRST SUBMATRIX FROM THE TOP,
C          2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS,
C
C        RV1 IS A TEMPORARY STORAGE ARRAY.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     XMCHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING RELATIVE PREC.
C
            XMCHEP=TWO**(-216)
c       xmchep=.929q-323
CCDC  XMCHEP=TWO**(-42)
C
      IERR = 0
      K = 0
      TAG = 0
C
      DO 100 I = 1, N
      W(I) = D(I)
      IF (I .NE. 1) RV1(I-1) = E(I)
  100 CONTINUE
C
      E2(1) = ZERO
      RV1(N) = ZERO
C
      DO 360 L = 1, N
      J = 0
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
  120 DO 140 M = L, N
      IF (M .EQ. N) GO TO 160
      IF (ABS(RV1(M)) .LE. XMCHEP * (ABS(W(M)) + ABS(W(M+1)))) GO TO
     +     160
C     ********** GUARD AGAINST UNDERFLOWED ELEMENT OF E2 **********
      IF (E2(M+1) .EQ. ZERO) GO TO 180
  140 CONTINUE
C
  160 IF (M .LE. K) GO TO 200
      IF (M .NE. N) E2(M+1) = ZERO
  180 K = M
      TAG = TAG + 1
  200 P = W(L)
      IF (M .EQ. L) GO TO 280
      IF (J .EQ. 30) GO TO 380
      J = J + 1
C     ********** FORM SHIFT **********
      G = (W(L+1) - P) / (TWO * RV1(L))
      R = SQRT(G*G+ONE)
      G = W(M) - P + RV1(L) / (G + SIGN(R,G))
      S = ONE
      C = ONE
      P = ZERO
      MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
      DO 260 II = 1, MML
      I = M - II
      F = S * RV1(I)
      B = C * RV1(I)
      IF (ABS(F) .LT. ABS(G)) GO TO 220
      C = G / F
      R = SQRT(C*C+ONE)
      RV1(I+1) = F * R
      S = ONE / R
      C = C * S
      GO TO 240
  220 S = F / G
      R = SQRT(S*S+ONE)
      RV1(I+1) = G * R
      C = ONE / R
      S = S * C
  240 G = W(I+1) - P
      R = (W(I) - G) * S + TWO * C * B
      P = S * R
      W(I+1) = G + P
      G = C * R - B
  260 CONTINUE
C
      W(L) = W(L) - P
      RV1(L) = G
      RV1(M) = ZERO
      GO TO 120
C     ********** ORDER EIGENVALUES **********
  280 IF (L .EQ. 1) GO TO 320
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
      DO 300 II = 2, L
      I = L + 2 - II
      IF (P .GE. W(I-1)) GO TO 340
      W(I) = W(I-1)
      IND(I) = IND(I-1)
  300 CONTINUE
C
  320 I = 1
  340 W(I) = P
      IND(I) = TAG
  360 CONTINUE
C
      GO TO 400
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
  380 IERR = L
  400 RETURN
      END
*      DEBUG SUBCHK
*     ENDDEBUG
CEDIT  TINVTB
C**********************************************************************
      SUBROUTINE TINVTB(NM,N,D,E,E2,M,W,IND,Z,IERR,RV1,RV2,RV3,RV4,RV6)
C**********************************************************************
      IMPLICIT REAL*16(A-H,O-Z)
      INTEGER P,Q,R,S,GROUP
      DIMENSION D(N),E(N),E2(N),W(M),Z(NM,M),
     +       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N),IND(M)
      PARAMETER( ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, TENM5=1.0D-5)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C          0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE,
C
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES,
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER,
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C
C     ON OUTPUT-
C
C        ALL INPUT ARRAYS ARE UNALTERED,
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS,
C
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** XMCHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
            XMCHEP=TWO**(-216)
CCDC  XMCHEP=TWO**(-42)
C
      IERR = 0
      IF (M .EQ. 0) GO TO 680
      ITAG = 0
      ORDER = ONE - E2(1)
      Q = 0
C     ********** ESTABLISH AND PROCESS NEXT SUBMATRIX **********
  100 P = Q + 1
C
      DO 120 Q = P, N
      IF (Q .EQ. N) GO TO 140
      IF (E2(Q+1) .EQ. ZERO) GO TO 140
  120 CONTINUE
C     ********** FIND VECTORS BY INVERSE ITERATION **********
  140 ITAG = ITAG + 1
      S = 0
C
      DO 660 R = 1, M
      IF (IND(R) .NE. ITAG) GO TO 660
      ITS = 1
      X1 = W(R)
      IF (S .NE. 0) GO TO 220
C     ********** CHECK FOR ISOLATED ROOT **********
      XU = ONE
      IF (P .NE. Q) GO TO 160
      RV6(P) = ONE
      GO TO 600
  160 XNORM = ABS(D(P))
      IP = P + 1
C
      DO 180 I = IP, Q
  180 XNORM = XNORM + ABS(D(I)) + ABS(E(I))
C     ********** EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW **********
      EPS2 = TENM5 * XNORM/(Q-P+1)
      EPS3 = XMCHEP * XNORM
             UK = DBLE(Q-P+1)
CCDC  UK =  FLOAT(Q-P+1)
      EPS4 = UK * EPS3
      UK = EPS4 / SQRT(UK)
      S = P
  200 GROUP = 0
      GO TO 240
C     ********** LOOK FOR CLOSE OR COINCIDENT ROOTS **********
  220 IF (ABS(X1-X0) .GE. EPS2) GO TO 200
      GROUP = GROUP + 1
      IF (ORDER * (X1 - X0) .LE. ZERO) X1 = X0 + ORDER * EPS3
C     ********** ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR **********
  240 V = ZERO
C
      DO 300 I = P, Q
      RV6(I) = UK
      IF (I .EQ. P) GO TO 280
      IF (ABS(E(I)) .LT. ABS(U)) GO TO 260
C     ********** WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY **********
      XU = U / E(I)
      RV4(I) = XU
      RV1(I-1) = E(I)
      RV2(I-1) = D(I) - X1
      RV3(I-1) = ZERO
      IF (I .NE. Q) RV3(I-1) = E(I+1)
      U = V - XU * RV2(I-1)
      V = -XU * RV3(I-1)
      GO TO 300
  260 XU = E(I) / U
      RV4(I) = XU
      RV1(I-1) = U
      RV2(I-1) = V
      RV3(I-1) = ZERO
  280 U = D(I) - X1 - XU * V
      IF (I .NE. Q) V = E(I+1)
  300 CONTINUE
C
      IF (U .EQ. ZERO) U = EPS3
      RV1(Q) = U
      RV2(Q) = ZERO
      RV3(Q) = ZERO
C     ********** BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- **********
  320 DO 340 II = P, Q
      I = P + Q - II
      RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
      V = U
      U = RV6(I)
  340 CONTINUE
C     ********** ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP **********
      IF (GROUP .EQ. 0) GO TO 400
      J = R
C
      DO 380 JJ = 1, GROUP
  360 J = J - 1
      IF (IND(J) .NE. ITAG) GO TO 360
      XU = DDOT(Q-P+1,RV6(P),1,Z(P,J),1)
C
      CALL DAXPY(Q-P+1,-XU,Z(P,J),1,RV6(P),1)
C
  380 CONTINUE
C
  400 XNORM = ZERO
C
      DO 420 I = P, Q
  420 XNORM = XNORM + ABS(RV6(I))
C
      IF (XNORM .GE. ONE) GO TO 560
C     ********** FORWARD SUBSTITUTION **********
      IF (ITS .EQ. 5) GO TO 540
      IF (XNORM .NE. ZERO) GO TO 440
      RV6(S) = EPS4
      S = S + 1
      IF (S .GT. Q) S = P
      GO TO 480
  440 XU = EPS4 / XNORM
C
      DO 460 I = P, Q
  460 RV6(I) = RV6(I) * XU
C     ********** ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE **********
  480 DO 520 I = IP, Q
      U = RV6(I)
C     ********** IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS **********
      IF (RV1(I-1) .NE. E(I)) GO TO 500
      U = RV6(I-1)
      RV6(I-1) = RV6(I)
  500 RV6(I) = U - RV4(I) * RV6(I-1)
  520 CONTINUE
C
      ITS = ITS + 1
      GO TO 320
C     ********** SET ERROR -- NON-CONVERGED EIGENVECTOR **********
  540 IERR = -R
      XU = ZERO
      GO TO 600
CDEL
C     ********** NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER **********
  560 U = ZERO
C
      DO 580 I = P, Q
  580 U = U + RV6(I)**2
C
      XU = ONE / SQRT(U)
C
  600 DO 620 I = 1, N
  620 Z(I,R) = ZERO
C
      DO 640 I = P, Q
  640 Z(I,R) = RV6(I) * XU
C
      X0 = X1
  660 CONTINUE
C
      IF (Q .LT. N) GO TO 100
  680 RETURN
      END
*      DEBUG SUBCHK
*     ENDDEBUG
CEDIT  TQL2
C*************************************
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C*************************************
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION D(N),E(N),Z(NM,N)
      PARAMETER(ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1,
C
C        E HAS BEEN DESTROYED,
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** XMCHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
            XMCHEP=TWO**(-216)
c       xmchep=.989 q-323
CCDC  XMCHEP=TWO**(-42)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 400
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = ZERO
      B = ZERO
      E(N) = ZERO
C
      DO 300 L = 1, N
      J = 0
      H = XMCHEP * (ABS(D(L)) + ABS(E(L)))
      IF (B .LT. H) B = H
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
      DO 120 M = L, N
      IF (ABS(E(M)) .LE. B) GO TO 140
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  120 CONTINUE
C
  140 IF (M .EQ. L) GO TO 280
  160 IF (J .EQ. 30) GO TO 380
      J = J + 1
C     ********** FORM SHIFT **********
      L1 = L + 1
      G = D(L)
      P = (D(L1) - G) / (TWO * E(L))
      R = SQRT(P*P+ONE)
      D(L) = E(L) / (P + SIGN(R,P))
      H = G - D(L)
C
      DO 180 I = L1, N
  180 D(I) = D(I) - H
C
      F = F + H
C     ********** QL TRANSFORMATION **********
      P = D(M)
      C = ONE
      S = ZERO
      MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
      DO 260 II = 1, MML
      I = M - II
      G = C * E(I)
      H = C * P
      IF (ABS(P) .LT. ABS(E(I))) GO TO 200
      C = E(I) / P
      R = SQRT(C*C+ONE)
      E(I+1) = S * P * R
      S = C / R
      C = ONE / R
      GO TO 220
  200 C = P / E(I)
      R = SQRT(C*C+ONE)
      E(I+1) = S * E(I) * R
      S = ONE / R
      C = C * S
  220 P = C * D(I) - S * G
      D(I+1) = H + S * (C * G + S * D(I))
C     ********** FORM VECTOR **********
      DO 240 K = 1, N
      H = Z(K,I+1)
      Z(K,I+1) = S * Z(K,I) + C * H
      Z(K,I) = C * Z(K,I) - S * H
  240 CONTINUE
C
  260 CONTINUE
C
      E(L) = S * P
      D(L) = C * P
      IF (ABS(E(L)) .GT. B) GO TO 160
  280 D(L) = D(L) + F
  300 CONTINUE
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
      DO 360 II = 2, N
      I = II - 1
      K = I
      P = D(I)
C
      DO 320 J = II, N
      IF (D(J) .GE. P) GO TO 320
      K = J
      P = D(J)
  320 CONTINUE
C
      IF (K .EQ. I) GO TO 360
      D(K) = D(I)
      D(I) = P
C
      DO 340 J = 1, N
      P = Z(J,I)
      Z(J,I) = Z(J,K)
      Z(J,K) = P
  340 CONTINUE
C
  360 CONTINUE
C
      GO TO 400
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
  380 IERR = L
  400 RETURN
      END
*      DEBUG SUBCHK
*     ENDDEBUG
CEDIT  TRBK3B
C*************************************
      SUBROUTINE TRBK3B(NM,N,NV,A,M,Z)
C*************************************
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION A(NV),Z(NM,M)
      DATA ZERO/0.0D0/, ONE/1.0D0/, TWO/2.0D0/
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3B.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS
C          USED IN THE REDUCTION BY  TRED3B IN ITS FIRST
C          N*(N+1)/2 POSITIONS,
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT-
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 140
      IF (N .EQ. 1) GO TO 140
C
      DO 120 I = 2, N
      L = I - 1
      IZ = (I * L) / 2
      IK = IZ + I
      H = A(IK)
      IF (H .EQ. ZERO) GO TO 120
C
      DO 100 J = 1, M
      S = -DDOT(L,A(IZ+1),1,Z(1,J),1)
C
C     ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********
      S = (S / H) / H
C
      CALL DAXPY(L,S,A(IZ+1),1,Z(1,J),1)
C
  100 CONTINUE
C
  120 CONTINUE
C
  140 RETURN
      END
