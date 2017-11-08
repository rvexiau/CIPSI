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
         real*8 dzop(60),dsc(3600),dsn(600)
      COMMON/INFO/NPOT(80),APOT(80),CPOT(80),RZ(4),OVZ(3),EX(4),OZ(2),
     * ARAF,AOP,BRAF,BOP,DMNXP,DMNXP2,OPTES,BASTES,NSPOT,NEWBAS,NSTARP
     * ,ISUP,QVER
      LOGICAL*4 NEWBAS,QPAR,QVER,qcst,tabu,extrac
      COMMON/CHAR/LABEL,SORTIE
       CHARACTER*8LAB,SORTIE
       CHARACTER*32 LABEL
       CHARACTER*8FMT,APSD
      DIMENSION FMT(1),NBTYP(5),grid(3)
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
                do i=1,5
       nbtyp(i)=0
          enddo

      READ(5,INPUT,END=999)
      write(6,*)' Dans le cas du SO les operateurs Christiansen'
      write(6,*)' doivent etre multiplies par -1. variable facso'
       write(6,*)' ici facso vaut : ',facso
       IF(QVER)THEN
	rewind 22
         do j=1,5
           read(22) LABEL(1:8),Lsr,NOP,Nopr,(dZOP(I),I=1,NOP),
     *    (dSC(I),I=1,NOP*nopr),(dSN(I),I=1,Nopr)
         write(6,*)label,lsr,nop,nopr
         if(lsr.eq.lsym)goto9
        enddo
 9         continue
         doi=1,nop
           zop(i)=dzop(i)
           sn(i)=dsn(i)
         enddo
         do i=1,nop*nopr
          sc(i)=dsc(i)
           enddo

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
      DO11I=1,LSYM-1
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
       do i=nstarp,nspot
       if(abs(apot(i)).lt.1.d-6.or.abs(cpot(i)).lt.1.d-6)then
       write(6,*)' terme semi-local nul verifier les donnees'
       stop
       endif
       enddo
       if(qcst.and.lcst.eq.0)then
       write(6,*)' qcst true        lcst nul verifier les donnees'
       stop
       endif

1     CONTINUE
C DEFINITION OF THE EXPANSION FUNCTIONS FOR THE POTENTIAL
C DEFINITION OF THE COLLOCATION BASIS SET
      CALL SSLOC(NBAS,NBASLD,NOP,NOPR,Z,ZOP,L,L,TORT,TOPE,REC,
     * EPSDS,RNOR,OPNOR)
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
        RV=RV+CPOT(J)*RA**(NPOT(J)+2)*EXP(-RA*RA*APOT(J))
       enddo
        if(abs(rv).lt.1.d-6)iflag=iflag-1
       if(iflag.gt.0)WRITE(6,'(3(a,d20.6))')
     * ' r= ',RA,' pspot=',RV,' -Z/r',-float(nval)/ra
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

        write(6,*)' test densite / pot'
       ij=1
       do j=2,nop
       ij=ij+j
       write(6,'(1x,5e15.6)') okl(ij-j),okl(ij),okl(ij-1),dsc(ij-1)
     * ,okl(ij-1)-(okl(ij)+okl(ij-j))*.5*dsc(ij-1)
       enddo


      WRITE(6,*)' ope sl -Z/R ds base de projection'
      IJ=0
       DO  I=1,Nop
       WRITE(6,42) zop(i),(temp(ij+j),J=1,I)
       IJ=IJ+I
        enddo

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
      WRITE(6,*)' ope projete ds base extraction'
      IJ=0
      DO  I=1,NBAS
      WRITE(6,42) z(i),(vKL(IJ+J),J=1,I)
      IJ=IJ+I
        enddo
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
       IF(J.LT.NBASLD-6) XNOR=XNOR+(TEMP(J)*TEMP(J))
3100    CONTINUE
         IK=IK+NOP
3000       WRITE(6,585)(TEMP(J),J=1,I)
        XNOR=SQRT(XNOR)
       WRITE(6,*)' NORME VERIF (1,NBAS-6) ',XNOR
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
      WRITE(6,206) IJ,IJP
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
 
          IF(SORTIE(1:7).EQ.'MONSTER')THEN
           WRITE(6,*)' MONSTER SAVED '
      WRITE(8,1001)NOPS,NVOIS
      WRITE(8,1002) (ZOP(I),I=1,NOPS)
      WRITE(8,1002) (SC(I),I=1,NOPS*NVOIS)
      WRITE(8,1002) (SN(I),I=1,NVOIS)
 1001 FORMAT(2I3)
 1002 FORMAT(4D20.12)
          ELSE  IF(SORTIE(1:4).EQ.'MELD')THEN
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
