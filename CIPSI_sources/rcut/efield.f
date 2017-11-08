      subroutine efield
Ccc   RG : modified on 25/11/08 : Differential calculation modifications      
Ccc   RG : modified on 26/11/08 : Output of the operators disabled
Ccc   RG : modified on 06/05/09 : Now reads the CPP operator from rpol
      implicit double precision (a-h,o-z)
      double precision i1sr4
      integer ent,sor,pl,pl1,plmax,plder,casp
      parameter (itypmax=3,plder=3,npar=2000000,lbloc=16384)
C
      include 'pshf.prm'
      logical*1 iandj,iecrit
      logical ijji
      common/ecrit/iecrit
c     common/output/nprint,itol,icut,normf,normp
      common/rayona/rcut1(dc,0:plder)
      common/rayonc/rcut(0:plder)
      common/infoa/nat,nnp,c(3,dc)
      common/infpot/iatno(dc),ipseud(dc)
      common/recomb/ s(225,20),nvar,itpmax
      include 'NSHEL.cmm'
      common/nocas/casp
      common/parinv/xi,yi,zi,xj,yj,zj,cx,cy,cz,lit,ljt,i1,j1,i2,j2,
     1              mini,maxi,mazi,minj,maxj,mazj,loci,locj
      common/alf/ai,aj
      common/project/pl
      common/nexpres/nint,noij
      common/par/ijmax,npartab,plmax
      common/iselec/ival,ivalf
      common/pmat/ pspt((itypmax+1)*(itypmax+2)*(itypmax+1)*
     1                             (itypmax+2)/4,0:4)
      common/ndebfin/ndeb(2,0:3,(itypmax+1)*(itypmax+2)/2),
     1               nfin(2,0:3,(itypmax+1)*(itypmax+2)/2)
      common /viandj/iandj
      integer*2 tabnew
      common/tabdat/tabnew(npar)
      common/lecrit/ent,sor
      dimension dij(225)
      dimension i1sr4(doas),
     1          elx(doas),
     1          ely(doas),
     1          elz(doas)
      data sqrt3/1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/

      logical diffCalc,firstRun,exist
      common /diff/ diffCalc,firstRun
      dimension tmpvek(doas),tmpvek2(doas)
      dimension tmpop(4*nat*nnp)
      character*2 ac2
c
c
      ia(i)=i*(i-1)/2
c     norm=normf.ne.1.or.normp.ne.1
c
c
c
      ioffset=0
      do 10000 ic=1,nat
      write(ac2,'(i2)') ic
      iato=iatno(ic)
      if(ipseud(ic).eq.0) go to 10000
      do 2 pl1=0,plmax
      rcut(pl1)=rcut1(ic,pl1)
    2 continue
      rewind 10
      write(sor,*) 'ic=',ic,'pl,rcut(pl)=',(pl1,rcut(pl1),pl1=0,plmax)
      do 5 i=1,nnp
      i1sr4(i)=0.0d0
      elx(i)=0.d0
      ely(i)=0.d0
      elz(i)=0.d0
      tmpvek=0.d0
      tmpvek2=0.d0
    5 continue
      cx=c(1,ic)
      cy=c(2,ic)
      cz=c(3,ic)
      do 9500 pl=0,plmax
      write(sor,*) 'pl de efield =',pl
      read(10) (((ndeb(i,j,k),nfin(i,j,k),k=1,ijmax),j=0,3),i=1,2)
        nbloc=(npartab*nfin(2,3,ijmax))/lbloc
        if((npartab*nfin(2,3,ijmax)).ne.lbloc*nbloc) nbloc=nbloc+1
        do 3788 i=1,nbloc
        mdeb=(i-1)*lbloc+1
        mfin=mdeb+lbloc
        if(i.eq.nbloc) mfin=npartab*nfin(2,3,ijmax)
        read(10) (tabnew(l),l=mdeb,mfin)
 3788   continue
 
c
c     ----- ishell
c
      do 9000 ii=1,nshell
      i=katom(ii)
      xi=c(1,i)
      yi=c(2,i)
      zi=c(3,i)
      i1=kstart(ii)
      i2=i1+kng(ii)-1
      lit=ktype(ii)-1
      if((lit+1).gt.5) then
      write(sor,*) 'cas imprevu : lit+1 > 5 , incompatible avec hondo'
                       stop
      endif
      mini=kmin(ii)
      maxi=kmax(ii)
      mazi=kmaz(ii)
      loci=kloc(ii)-mini
c
c     ----- jshell
c
      do 8000 jj=1,ii
c     if(iecrit) write(sor,*) 'ii=',ii,'jj=',jj
      nvar=1
      itpmax=4
      ival=1
      ivalf=2
      ijji=.false.
      j=katom(jj)
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      j1=kstart(jj)
      j2=j1+kng(jj)-1
      ljt=ktype(jj)-1
      if((ljt+1).gt.5) then
      write(sor,*) 'cas imprevu : ljt+1 > 5 ,incompatible avec hondo'
                       stop
      endif
      minj=kmin(jj)
      maxj=kmax(jj)
      mazj=kmaz(jj)
      locj=kloc(jj)-minj
c     iandj=ii.eq.jj
      iandj=.false.
c
c
      if(lit.lt.ljt)then
                ijji=.true.
                call invers
      endif
      call preintn
c
c
      do 40 j=1,4
      do 20 i=1,nint
      s(i,j)=0.0d0
 20   continue
 40   continue
c
c
      if(casp.ne.0)then
c
             if(ndeb(1,casp,noij).gt.nfin(1,casp,noij))then
c
c
                        if(ndeb(2,casp,noij).gt.nfin(2,casp,noij))then
                                if(ijji) call invers
                                go to 8000
                                                                  else
                                nvar=2
                                ival=2
                                go to 50
                        endif
c
c
                                                       else
c
c
                        if(ndeb(2,casp,noij).gt.nfin(2,casp,noij))then
                                itpmax=1
                                ivalf=1
                        endif
c
c
             endif
c
c
      endif
c
c
c
c     ----- i primitive
c
c
 50   jgmax=j2
c     if(iecrit) write(sor,*)'ig=i1,i2 :',i1,i2
      do 7000 ig=i1,i2
      ai=ex(ig)
      csi=cs(ig)
      cpi=cp(ig)
      cdi=cd(ig)
      cfi=cf(ig)
      cgi=cg(ig)
c
c
c     ----- j primitive
c
c     if(iandj) jgmax=ig
      do 6000 jg=j1,jgmax
c     if(iecrit) write(sor,*) 'ig=',ig,'jg=',jg
      aj=ex(jg)
c     if(iecrit) write(sor,*) 'ai=',ai,'aj=',aj
      csj=cs(jg)
      cpj=cp(jg)
      cdj=cd(jg)
      cfj=cf(jg)
      cgj=cg(jg)
c
c     ----- density factor
c
c     double=iandj.and.ig.ne.jg
      max=maxj
      nn=0
      do 310 i=mini,maxi
      go to ( 70, 80,180,180, 90,180,180,100,180,180,
     1       110,180,180,120,180,180,180,180,180,130,
     2       140,180,180,150,180,180,180,180,180,160,
     3       180,180,170,180,180),i
   70 dum1=csi
      go to 180
   80 dum1=cpi
      go to 180
   90 dum1=cdi
      go to 180
c 100 if (norm) dum1=dum1*sqrt3
  100 dum1=dum1*sqrt3
      go to 180
  110 dum1=cfi
      go to 180
c 120 if(norm) dum1=dum1*sqrt5
  120 dum1=dum1*sqrt5
      go to 180
c 130 if(norm) dum1=dum1*sqrt3
  130 dum1=dum1*sqrt3
      go to 180
  140 dum1=cgi
      go to 180
c 150 if(norm) dum1=dum1*sqrt7
  150 dum1=dum1*sqrt7
      go to 180
c 160 if(norm) dum1=dum1*sqrt5/sqrt3
  160 dum1=dum1*sqrt5/sqrt3
      go to 180
c 170 if(norm) dum1=dum1*sqrt3
  170 dum1=dum1*sqrt3
c 180 if(iandj) max=i
  180 do 310 j=minj,max
      go to (190,200,300,300,210,300,300,220,300,300,
     1       230,300,300,240,300,300,300,300,300,250,
     2       260,300,300,270,300,300,300,300,300,280,
     3       300,300,290,300,300),j
  190 dum2=dum1*csj
c     if(.not.double) go to 300
c     if(i.gt.1) go to 195
c     dum2=dum2+dum2
c     go to 300
c 195 dum2=dum2+csi*cpj
      go to 300
  200 dum2=dum1*cpj
c     if(double) dum2=dum2+dum2
      go to 300
  210 dum2=dum1*cdj
c     if(double) dum2=dum2+dum2
      go to 300
c 220 if(norm) dum2=dum2*sqrt3
  220 dum2=dum2*sqrt3
      go to 300
  230 dum2=dum1*cfj
c     if(double) dum2=dum2+dum2
      go to 300
c 240 if(norm) dum2=dum2*sqrt5
  240 dum2=dum2*sqrt5
      go to 300
c 250 if(norm) dum2=dum2*sqrt3
  250 dum2=dum2*sqrt3
      go to 300
  260 dum2=dum1*cgj
c     if(double) dum2=dum2+dum2
      go to 300
c 270 if(norm) dum2=dum2*sqrt7
  270 dum2=dum2*sqrt7
      go to 300
c 280 if(norm) dum2=dum2*sqrt5/sqrt3
  280 dum2=dum2*sqrt5/sqrt3
      go to 300
c 290 if(norm) dum2=dum2*sqrt3
  290 dum2=dum2*sqrt3
  300 nn=nn+1
  310 dij(nn)=dum2
c
c     ..... electric field
c
c
      call intnumx
c
c
      ivalj=4
      if(ivalf.eq.1) ivalj=1
      do 450 j=ival,ivalj
      do 440 i=1,nint
      s(i,j)=s(i,j)+pspt(i,j-1)*dij(i)
  440 continue
  450 continue
c
c
 6000 continue
 7000 continue
c
c
      if((lit.ge.2.and.lit.le.4).or.(ljt.ge.2.and.ljt.le.4))
     1 call comps(lit+1,ljt+1,iandj)
      max=mazj
      ij=0
      do 7100 i=mini,mazi
      if(iandj) max=i
      do 7100 j=minj,max
      ij=ij+1
      if(ii.eq.jj .and. j.gt.i) go to 7100
      nij=ia(loci+i)+locj+j
      if(ijji)nij=ia(locj+j)+loci+i
      i1sr4(nij)=i1sr4(nij)+s(ij,1)
      elx(nij)=elx(nij)+s(ij,2)
      ely(nij)=ely(nij)+s(ij,3)
      elz(nij)=elz(nij)+s(ij,4)
 7100 continue
      if(ijji)then
                call invers
      endif
 8000 continue
 9000 continue
c
cCc       write(sor,7602)(i,elx(i),ely(i),elz(i),i=1,nnp)
c
      if(pl.eq.3)then
cCc         write(sor,*) 'on ecrit : i, i1sr4'
cCc         write(sor,7604)(i,i1sr4(i),i=1,nnp)
7602    format(/,(i5,3d15.8))
7604    format(/,4(i5,f8.5))
      endif
c
c
 9500 continue
cCc ============================
cCc = Differential calculation =
cCc ============================
      if(diffCalc) then
        if(firstRun) then
cCc =================================================
cCc = Here, save the calculated operators in a file =
cCc =================================================
          open(unit=20,file='TMPOPS'//trim(adjustl(ac2)),
     +         form='unformatted',access='sequential')
          write(20) elx(1:nnp)
          write(20) ely(1:nnp)
          write(20) elz(1:nnp)
          write(20) i1sr4(1:nnp)
          close(20)
          write(*,*) 'Diff calc. Saving the ops. from the first run'
          write(*,*) 'Diff calc. Reading CPP operator from rpol'
!!$  >>> debug <<<          
          write(*,*) 'flag1 ',4*nat*nnp 
!!$  >>> debug <<<          
          inquire(file='VCPP_RPOL',exist=exist)
          if(.not.exist) then
            stop 'file VCPP_RPOL does not exist for the differential'
     +' calculation'
          endif
          open(unit=20,file='VCPP_RPOL',form='unformatted')
          iii=0
          do i=1,4*nat
            read(20) tmpop(1+iii:nnp+iii)
            iii=iii+nnp
          enddo
          close(20)
          write(*,*) 'Have read VCPP_RPOL'
        else
cCc ============================================
cCc = Here comes the differential calculation. =
cCc ============================================ 
          open(unit=20,file='TMPOPS'//trim(adjustl(ac2)),
     +         form='unformatted',access='sequential')
          open(unit=21,file='COMB_VCPP',form='unformatted',
     +         access='append')
          read(20) tmpvek(1:nnp)   ! that's elx
          write(*,*) 'flag1 ',1+ioffset,nnp+ioffset
          tmpvek2(1:nnp)=tmpop(1+ioffset:nnp+ioffset)
          elx=elx-tmpvek
          write(21) (elx(i)+tmpvek2(i),i=1,nnp) 
          ioffset=ioffset+nnp
          read(20) tmpvek(1:nnp)   ! that's ely
          write(*,*) 'flag2 ',1+ioffset,nnp+ioffset 
          tmpvek2(1:nnp)=tmpop(1+ioffset:nnp+ioffset)
          ely=ely-tmpvek
          write(21) (ely(i)+tmpvek2(i),i=1,nnp) 
          ioffset=ioffset+nnp
          read(20) tmpvek(1:nnp)   ! that's elz
          write(*,*) 'flag3 ',1+ioffset,nnp+ioffset
          tmpvek2(1:nnp)=tmpop(1+ioffset:nnp+ioffset)
          elz=elz-tmpvek
          write(21) (elz(i)+tmpvek2(i),i=1,nnp)
          ioffset=ioffset+nnp
          read(20) tmpvek(1:nnp)   ! that's i1sr4
          write(*,*) 'flag4 ',1+ioffset,nnp+ioffset
          tmpvek2(1:nnp)=tmpop(1+ioffset:nnp+ioffset)
          i1sr4=i1sr4-tmpvek
          write(21) (i1sr4(i)+tmpvek2(i),i=1,nnp)
          ioffset=ioffset+nnp
          close(20)
          close(21)
          write(*,*) 'Diff calc'
          write(*,*) 
     +'Have substracted the contribution from the first run'
          write(*,*) '**********************************************'
        endif
      endif
c
      if(.not.firstRun) then
        write(17) (elx(i),i=1,nnp)
        write(17) (ely(i),i=1,nnp)
        write(17) (elz(i),i=1,nnp)
        write(17) (i1sr4(i),i=1,nnp)
      endif
10000 continue
c
c
      write(sor,9400)
 9400 format(' ......  end of 1/r**4, end of electric field ......')
      return
      end
