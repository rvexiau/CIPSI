      subroutine atoms(iunt,pseud)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical redux,pseud
      character*4 pname,dname,fname,gname
      character*2 label
      character*4 dfg(35),dfg1(35),dfg2(35),sname,iflab
      character*3 number(dc),atnam(88)
      character*1 blk
      character*4 dn1,dn2,fn1,fn2,gn1,gn2
      character*80 apsd,apseud,amolcas
      character*8 anam
      character*8 line,linev,atom,atid,a
      integer ddmax,ffmax,ggmax
      common/output/iop1,itol,icut,normf,normp
      common/iofile/ir,iw,ip,is
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc)
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell,kmaz(ds)
      common/runlab/iflab(3,doa)
      common/symtry/invt(48),iso(ds,12),ptr(3,144),dtr(6,288)
     1,ftr(10,480),nt,ict(dc,48)
      dimension gtr(15,720)
      common/transf/xsmal,ysmal,zsmal,xnew,ynew,znew,xp,yp,zp
      common/fil2/ a(dc),ns(dc),ks(dc),newsh(ds,48),ngauss
      common/frame/u1,u2,u3,v1,v2,v3,w1,w2,w3,x0,y0,z0
      logical geoto,gestop
      common/geot/geoto,gestop
      common/psnloc/apseud(dc),iatno(dc)
      common/isopac/indin(48),indout(12)
c
      namelist/cent/atom,atid,znuc,iatn,x,y,z,ntypg,ncong,exg,csg,cpg,
     1 cdg,cfg,cgg,apsd
      dimension ntypg(100),ncong(100),exg(150),csg(100),cpg(100),
     +cdg(100),cfg(100),cgg(100)
cCc       dimension ntypg(62),ncong(62),exg(150),csg(62),cpg(62),cdg(62),
cCc      1 cfg(62),cgg(62)
cCc       dimension ntypg(50),ncong(50),exg(150),csg(50),cpg(50),cdg(50),
cCc      1 cfg(50),cgg(50)
      dimension b(dc),nbfs(10)
c
      dimension ind(ds),label(10)
      dimension pname(3),dname(6),fname(10),gname(15)
      dimension dn1(6),dn2(5),fn1(10),fn2(7),gn1(15),gn2(9)
c
      dimension csinp(dgp),cpinp(dgp),cdinp(dgp),cfinp(dgp),cginp(dgp)
      data sname/' s'/,blk/' '/
      data atnam/'h ','he','li','be','b ','c ','n ','o ','f ','ne',
     1                     'na','mg','al','si','p ','s ','cl','ar',
     2 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu','zn','ga',
     3 'ge','as','se','br','kr',
     4 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd','in',
     5 'sn','sb','te','i ','xe','cs','ba','la','ce','pr','nd','pm','sm',
     6 'eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta','w ','re',
     7 'os','ir','pt','au','hg','tl','pb','bi','po','at','rn','fr','ra'/
      data pname,dn1/'x ','y ','z ','xx','yy','zz','xy','xz','yz'/
      data dn2/'dsig','dpix','dpiy','ddxx','ddxy'/
      data dname/' ',' ',' ',' ',' ',' '/
      data fname/' ',' ',' ',' ',' ',' ',' ',' ',' ',' '/
      data fn1   /'xxx ','yyy ','zzz ','xxy ','xxz ','yyx ','yyz ',
     1 'zzx ','zzy ','xyz '/
      data fn2   /'fsig','fpix','fpiy','fdxx','fdxy','ffix','ffiy'/
      data gname /'xxxx','yyyy','zzzz','xxxy','xxxz','yyyx','yyyz ',
     1 'zzzx','zzzy','xxyy','xxzz','yyzz','xxyz','yyxz','zzxy'/
      data nbfs/1,3,5,7,9,1,4,6,10,15/
      data label/'s ','p ','d ','f ','g ','k ','l ','m ','n ','o '/
      data sqrt3 /1.73205080756888d+00/
      data sqrt5 /2.23606797749979d+00/
      data sqrt7 /2.64575131106459d+00/
      data pi32 /5.56832799683170d+00/
      data tol /1.0d-10/
      data line,linev /' * * * *','       *'/
C changed for LINUX compatibility 26/11/97 Malte GROSS
c     equivalence (dfg(1),sname),(dfg(2),pname(1))
c     equivalence(dfg(5),dname(1)),(dfg(11),fname(1)),(dfg(21),gname(1))
c     equivalence(dfg1(5),dn1(1)),(dfg1(11),fn1(1)),(dfg1(21),gn1(1))
c     equivalence(dfg2(5),dn2(1)),(dfg2(11),fn2(1)),(dfg2(21),gn2(1))
      dfg2=' '
      do i=1,dc
        write(number(i),'(i3)') i
       end do
      dfg(1)=sname
      do i=1,3
         dfg(i+1)=pname(i)
      end do
      do i=1,6
         dfg(i+4)=dname(i)
      end do
      do i=1,10
         dfg(i+10)=fname(i)
      end do
      do i=1,15
         dfg(i+20)=gname(i)
      end do
      do i=1,6
         dfg1(i+4)=dn1(i)
      end do
      do i=1,10
         dfg1(i+10)=fn1(i)
      end do
      do i=1,15
         dfg1(i+20)=gn1(i)
      end do
      do i=1,5
         dfg2(i+4)=dn2(i)
      end do
      do i=1,7
         dfg2(i+10)=fn2(i)
      end do
      do i=1,9
         dfg2(i+20)=gn2(i)
      end do
C changed for LINUX compatibility 26/11/97 Malte GROSS

 8888 format(i5,4x,a1,i5)
 8887 format(a8,i6,f5.0,3f20.12)
 8886 format(i5,e15.9,2e20.12)
 8885 format(/,20x,20(1h*),//,20x,'molecular basis set',//,
     1 20x,20(1h*))
 8884 format(//,20x,43(2h* ),//,20x,' atom pseudo    atomic   atomic',
     116x,'coordinates',17x,'number of',/,36x,'number charge',7x,'x',
     213x,'y',14x,'z',7x,'shells',//,20x,43(2h* ))
 8883 format(/,18x,2h* ,2a6,4x,i3,5x,f5.1,3(f12.7,3x),i5,10x,1h*)
 8882 format(/,28x,78(1h*),//,8x,' contracted primitive functions',//,
     1 18x,'shell type prim      exponents   ',
     2 9x,'contraction coefficients')
 8881 format(18x,i3,4x,a1,3x,i3,1x,2f15.6,3h  (,f12.6,3h  ))
 8880 format(/' there is no identical center ',a8,' for center',a8)
 8877 format(/,' total number of shells ',11x,'=',i5,/
     1 ,' total number of basis functions   =',i5,/,
     1 ' number of electrons ',14x,'=',i5,/,
     1 ' charge of molecule ',15x,'=',i5,/,
     1 ' sz                 ',15x,'=',i5,'/2',/,
     1 ' number of alpha orbitals',6x,'=',i5,/,
     1 ' number of beta  orbitals',6x,'=',i5,/,
     1 ' total number of atoms',13x,'=',i5)
 8876 format(' excessive number of atoms')
 8875 format(' excessive number of shells')
 8874 format(' excessive number of basis functions')
 8873 format(' excessive number of contracted primitive functions')
 8872 format(' excessive contraction into shell')
 8871 format(92x,5a2)
 8870 format(/,1x,5h*****,' transformation table of atoms',
     1 5h*****,/,30x,' rows are atoms',/,30x,
     1 ' columns are symmetry operations')
 8869 format(//)
 8868 format(1x,16a8)
 8867 format(1x,a8,15(2x,i2,3x,1h*))
 8866 format(1x,16(2x,i2,3x,1h*))
 8865 format(/,1x,5h*****,' transformation table of shells',
     1 5h*****,/,30x,' rows are shells',/,30x,
     1 ' columns are symmetry operations')
 8864 format(' the contracted primitive functions have been ',
     1 'unnormalized')
 8863 format(/,' transformation of the basis functions',/)
 8862 format(8x,10(3x,a4,3x))
 8861 format(2x,a4,2x,10f10.6)
 8860 format(1h1)
 8859 format(3i5)
 8858 format(/,21x,'transformation number',i4,/)
 8855 format(/)
 8851 format(/,1x,5h*****,' inverse transformations ',5h*****,/)
 8850 format(' the contracted basis functions are now normalized to',
     1 ' unity')
 8849 format(18x,i3,4x,a1,3x,i3,1x,2f15.6,3h  (,f12.6,3h  ),
     2 f15.6,3h  (,f12.6,3h  ))
 8848 format(///,' namelist cent is missing or incorrect.. program '
     +'stops')
c
c     read in molecule information (charge,multiplicity,units)
c     read in unique centers and atomic basis sets grouped
c     in shells
c     generate all new centers
c     set table ( centers versus transformations )
c     set table ( shells versus transformations )
c     set matrix of transformation of the basis functions
c
      units=1.0d+00/0.529167d+00
      redux=.true.
      if(iop1.ne.4) write(iw,8885)
      do 1010 i=1,dgp
      ex(i)=0.0d+00
      cs(i)=0.0d+00
      cp(i)=0.0d+00
      cd(i)=0.0d+00
      cf(i)=0.0d+00
      cg(i)=0.0d+00
      csinp(i)=0.0d+00
      cpinp(i)=0.0d+00
      cdinp(i)=0.0d+00
      cfinp(i)=0.0d+00
      cginp(i)=0.0d+00
 1010 continue
      nat=0
      ne=0
      nshell=0
      loc=0
      ngauss=0
      if(geoto) call geoda
      if(geoto.and.gestop) stop
 1020 atid=' '
      znuc=-100000.d0
      atom='        '
      apsd=''
      x=0.d0
      y=0.d0
      z=0.d0
      do  i=1,100
      csg(i)=0.d0
      cpg(i)=0.d0
      cdg(i)=0.d0
      cfg(i)=0.d0
      cgg(i)=0.d0
      ncong(i)=0
      ntypg(i)=0
      end do
      do i=1,150
      exg(i)=0.d0
      end do
c
      read(ir,cent,end=6020)
      ishell=0
      igt=0
      igs=0
      igp=0
      igd=0
      igf=0
      igg=0
      if(atom.eq.'        ') go to 1700
      nat=nat+1
      a(nat)=atom
      ns(nat)=0
      ks(nat)=nshell+1
      iatno(nat)=iatn
      if(geoto) then
         x=c(1,nat)
         y=c(2,nat)
         z=c(3,nat)
      endif
      if(iunt.le.0) go to 1030
      if(geoto) then
         c(1,nat)=c(1,nat)*units
         c(2,nat)=c(2,nat)*units
         c(3,nat)=c(3,nat)*units
      else
         x=x*units
         y=y*units
         z=z*units
      endif
 1030 continue
      if(znuc.gt.0.d0) ne=ne+znuc
      zan(nat)=znuc
      if(.not.geoto) then
         c(1,nat)=x
         c(2,nat)=y
         c(3,nat)=z
      endif
C changed for LINUX compatibility 26/11/97 Malte GROSS
      if(iatn.le.0.or.iatn.eq.int(znuc+.5))go to 1130
c     if(iatn.le.0.or.iatn.eq.ifix(znuc))go to 1130
C changed for LINUX compatibility 26/11/97 Malte GROSS
c     check if there is pseudop. data for this atom on file 20
      rewind 20
      if(apsd.eq.'') apsd=atom
      if(apsd.eq.'pssl'.or.apsd.eq.'PSSL')then 
      pseud=.true.
      go to 1130
      end if
  102 read(20,end=100) anam,npmax
      if(anam.eq.apsd(1:8)) go to 110
      do  i=1,npmax
      read(20)
      end do
      go to 102
  100 continue
  103 read(23,end=108) amolcas,npmax
      if(amolcas(1:12).eq.apsd(1:12)) go to 110
      npmax=npmax+1
      do i=1,npmax
      read(23)
      end do
      go to 103
  108 write(iw,9995) apsd
 9995 format(//,'** no pseudopotential data on file for apsd=',a12,//)
      stop
  110 continue
       write(6,*) 'pseudo trouve'
 1130 apseud(nat)=apsd
      if(atid.eq.' ') go to 1040
      natm=nat-1
      do 2010 iat=1,natm
      if(atid.eq.a(iat)) go to 2020
 2010 continue
      write(iw,8880) atid,atom
      stop
 2020 nato=iat
      ns(nat)=ns(nato)
      ks(nat)=nshell+1
      nshell=nshell+ns(nat)
      ns1=ns(nat)
      if(ns1.eq.0) go to 1060
      do 2030 k=1,ns1
      j=ks(nato)+k-1
      jj=ks(nat)+k-1
      kstart(jj)=kstart(j)
      ktype(jj)=ktype(j)
      kmax(jj)=kmax(j)
      kmin(jj)=kmin(j)
      kmaz(jj)=kmaz(j)
       kng(jj)=kng(j)
      katom(jj)=nat
      kloc(jj)=loc+1
      loc=loc+nbfs(ktype(jj))
 2030 continue
      go to 1060
 1040 ishell=ishell+1
      igauss=ncong(ishell)
c
      if(igauss.le.0) go to 1060
      if(igauss.le.20) go to 1042
      write(iw,8872)
      write(*,*) 'igauss = ',igauss,ishell
      write(*,*) 'ncong = ',ncong
      stop
 1042 continue
      ityp=ntypg(ishell)
      nshell=nshell+1
      go to (1043,1044,1045,1046,1047,1043,1048,1031,1032,1033),ityp
 1043 min=1
      max=1
      maz=1
      go to 1049
 1044 min=2
      max=4
      maz=4
      go to 1049
 1045 min=5
      max=10
      maz=9
      go to 1049
 1046 min=11
      max=20
      maz=17
      go to 1049
 1047 min=21
      max=35
      maz=29
      go to 1049
 1048 min=1
      max=4
      maz=4
 1031 min=5
      max=10
      maz=10
      go to 1049
 1032 min=11
      max=20
      maz=20
      go to 1049
 1033 min=21
      max=35
      maz=35
 1049 kmin(nshell)=min
      kmax(nshell)=max
      kmaz(nshell)=maz
      kstart(nshell)=ngauss+1
      katom(nshell)=nat
      ktype(nshell)=ityp
      if(ityp.ge.8) redux=.false.
      kng(nshell)=igauss
      kloc(nshell)=loc+1
      ngauss=ngauss+igauss
      loc=loc+nbfs(ityp)
      k1=kstart(nshell)
      k2=k1+kng(nshell)-1
      ns(nat)=ns(nat)+1
      do 1050 k=k1,k2
c     read(ir,8886) kdum,ex(k),c1,c2
      igt=igt+1
      ex(k)=exg(igt)
      go to (2041,2042,2043,2044,2045,2041,2046,2043,2044,2045),ityp
 2041 igs=igs+1
      csinp(k)=csg(igs)
      go to 1034
 2042 igp=igp+1
      cpinp(k)=cpg(igp)
      go to 1034
 2043 igd=igd+1
      cdinp(k)=cdg(igd)
      go to 1034
 2044 igf=igf+1
      cfinp(k)=cfg(igf)
      go to 1034
 2045 igg=igg+1
      cginp(k)=cgg(igg)
      go to 1034
 2046 igs=igs+1
      csinp(k)=csg(igs)
      igp=igp+1
       cpinp(k)=cpg(igp)
 1034 cs(k)=csinp(k)
      cp(k)=cpinp(k)
      cd(k)=cdinp(k)
      cf(k)=cfinp(k)
      cg(k)=cginp(k)
 1050 continue
c
c     if(normp.ne.1) ... unnormalization of the primitive functions.
c     if contraction coefficients are given in terms of normalized
c     primitive functions, change them to go with unnormalized
c     primitives .
c     for d shells, the input coefficients cd must be the coefficients
c     corresponding to the normalized primitive x**2 *exp(-a*r**2).
c     for f shells, the input coefficients cf must be the coefficients
c     corresponding to the normalized primitive x**3 *exp(-a*r**2)
c     for g shells, the input coefficients cg must be the coefficients
c     corresponding to the normalized primitive x**4 *exp(-a*r**2)
c
      if(normp.eq.1) go to 1052
      do 1051 ig=k1,k2
      ee=ex(ig)+ex(ig)
      facs=pi32/(ee*sqrt(ee))
      facp=0.5d+00*facs/ee
      facd=0.75d+00*facs/(ee*ee)
      facf=1.875d+00*facs/(ee**3)
      facg=6.5625d+00*facs/(ee**4)
      cs(ig)=cs(ig)/sqrt(facs)
      cp(ig)=cp(ig)/sqrt(facp)
      cd(ig)=cd(ig)/sqrt(facd)
      cf(ig)=cf(ig)/sqrt(facf)
      cg(ig)=cg(ig)/sqrt(facg)
 1051 continue
 1052 continue
c
c     if(normf.ne.1) normalize the contracted basis functions.
c
      if(normf.eq.1) go to 1040
      facs=0.0d+00
      facp=0.0d+00
      facd=0.0d+00
      facf=0.0d+00
      facg=0.0d+00
      do 1054 ig=k1,k2
      do 1054 jg=k1,ig
      ee=ex(ig)+ex(jg)
      fac=ee*sqrt(ee)
      dums=cs(ig)*cs(jg)/fac
      dump=0.5d+00*cp(ig)*cp(jg)/(ee*fac)
      dumd=0.75d+00*cd(ig)*cd(jg)/(ee**2*fac)
      dumf=1.875d+00*cf(ig)*cf(jg)/(ee**3*fac)
      dumg=6.5625d+00*cg(ig)*cg(jg)/(ee**4*fac)
      if(ig.eq.jg) go to 1053
      dums=dums+dums
      dump=dump+dump
      dumd=dumd+dumd
      dumf=dumf+dumf
      dumg=dumg+dumg
 1053 facs=facs+dums
      facp=facp+dump
      facd=facd+dumd
      facf=facf+dumf
      facg=facg+dumg
 1054 continue
      do 1055 ig=k1,k2
      if(facs.gt.tol) cs(ig)=cs(ig)/sqrt(facs*pi32)
      if(facp.gt.tol) cp(ig)=cp(ig)/sqrt(facp*pi32)
      if(facd.gt.tol) cd(ig)=cd(ig)/sqrt(facd*pi32)
      if(facf.gt.tol) cf(ig)=cf(ig)/sqrt(facf*pi32)
      if(facg.gt.tol) cg(ig)=cg(ig)/sqrt(facg*pi32)
 1055 continue
      go to 1040
 1060 if(nt.eq.1) go to 1020
c
c     generate equivalent centers
c
      call local(x,y,z,xs,ys,zs)
      xsmal=xs
      ysmal=ys
      zsmal=zs
      nat0=nat
      do 1250 it=1,nt
      if(it.eq.1) go to 1250
      nn=9*(it-1)
      call trans(nn)
      call rot
      do 1150 iat=1,nat
      test=(xp-c(1,iat))**2+(yp-c(2,iat))**2+(zp-c(3,iat))**2
      if(test.le.tol) go to 1250
 1150 continue
      nat=nat+1
      c(1,nat)=xp
      c(2,nat)=yp
      c(3,nat)=zp
      ns(nat)=ns(nat0)
      ks(nat)=ks(nat-1)+ns(nat0)
      a(nat)=a(nat0)
      apseud(nat)=apseud(nat0)
      iatno(nat)=iatno(nat0)
      zan(nat)=zan(nat0)
      if(zan(nat).gt.0.d0) ne=ne+zan(nat)
      nshell=nshell+ns(nat)
      ns1=ns(nat)
      do 1200 k=1,ns1
      j=ks(nat0)+k-1
      jj=ks(nat)+k-1
      kmin(jj)=kmin(j)
      kmax(jj)=kmax(j)
      kmaz(jj)=kmaz(j)
      kstart(jj)=kstart(j)
      ktype(jj)=ktype(j)
      kng(jj)=kng(j)
      katom(jj)=nat
      kloc(jj)=loc+1
 1200 loc=loc+nbfs(ktype(jj))
 1250 continue
      go to 1020
 1700 continue
c
c     if(normp.ne.1) the contracted basis functions have been
c     expressed in terms of unnormalized primitive functions.
c
      if(normp.ne.1) write(iw,8864)
c
c     if(normf.ne.1) the contracted basis functions have been
c     normalized to unity
c
      if(normf.ne.1) write(iw,8850)
c
c     form transformation tables for atoms and shells.
c
      do 1750 iat=1,nat
      x=c(1,iat)
      y=c(2,iat)
      z=c(3,iat)
      ns1=ks(iat)-1
      ns2=ns(iat)
      call local(x,y,z,xs,ys,zs)
      xsmal=xs
      ysmal=ys
      zsmal=zs
      do 1749 it=1,nt
      nn=9*(it-1)
      call trans(nn)
      call rot
      do 1746 i=1,nat
      test=(xp-c(1,i))**2+(yp-c(2,i))**2+(zp-c(3,i))**2
      if(test.gt.tol) go to 1746
      ictr=i
      go to 1747
 1746 continue
 1747 ict(iat,it)=ictr
      ns3=ks(ictr)-1
      do 1748 ish=1,ns2
 1748 newsh(ns1+ish,it)=ns3+ish
 1749 continue
 1750 continue
      ntwd=(nt+3)/4
      do 1754 i=1,nshell
      do 1751 it=1,ntwd
 1751 indout(it)=0
      do 1752 it=1,nt
      indin(it)=0
 1752 indin(it)=newsh(i,it)
      call isoin(nt)
      do 1753 it=1,ntwd
 1753 iso(i,it)=indout(it)
 1754 continue
c
c     calculate transforms of p,d,f, and g functions
c     for all symetry operations.
c
      x=x0+1.0d+00
      y=y0
      z=z0
      call local(x,y,z,xs,ys,zs)
      xsmal=xs
      ysmal=ys
      zsmal=zs
      do 1760 it=1,nt
      nn=9*(it-1)
      call trans(nn)
      call rot
      n=3*(it-1)
      ptr(1,n+1)=xp-x0
      ptr(2,n+1)=yp-y0
      ptr(3,n+1)=zp-z0
 1760 continue
      x=x0
      y=y0+1.0d+00
      z=z0
      call local(x,y,z,xs,ys,zs)
      xsmal=xs
      ysmal=ys
      zsmal=zs
      do 1770 it=1,nt
      nn=9*(it-1)
      call trans(nn)
      call rot
      n=3*(it-1)
      ptr(1,n+2)=xp-x0
      ptr(2,n+2)=yp-y0
      ptr(3,n+2)=zp-z0
 1770 continue
      x=x0
      y=y0
      z=z0+1.0d+00
      call local(x,y,z,xs,ys,zs)
      xsmal=xs
      ysmal=ys
      zsmal=zs
      do 1780 it=1,nt
      nn=9*(it-1)
      call trans(nn)
      call rot
      n=3*(it-1)
      ptr(1,n+3)=xp-x0
      ptr(2,n+3)=yp-y0
      ptr(3,n+3)=zp-z0
 1780 continue
      do 1820 it=1,nt
      np=3*(it-1)
      nd=6*(it-1)
      nf=10*(it-1)
      ng=15*(it-1)
      do 1788 i=1,6
      go to (1781,1782,1783,1784,1785,1786) ,i
 1781 j=1
      k=1
      go to 1787
 1782 j=2
      k=2
      go to 1787
 1783 j=3
      k=3
      go to 1787
 1784 j=1
      k=2
      go to 1787
 1785 j=1
      k=3
      go to 1787
 1786 j=2
      k=3
 1787 dtr(1,nd+i)=ptr(1,np+j)*ptr(1,np+k)
      dtr(2,nd+i)=ptr(2,np+j)*ptr(2,np+k)
      dtr(3,nd+i)=ptr(3,np+j)*ptr(3,np+k)
      dtr(4,nd+i)=ptr(1,np+j)*ptr(2,np+k)
     1           +ptr(2,np+j)*ptr(1,np+k)
      dtr(5,nd+i)=ptr(1,np+j)*ptr(3,np+k)
     1           +ptr(3,np+j)*ptr(1,np+k)
      dtr(6,nd+i)=ptr(2,np+j)*ptr(3,np+k)
     1           +ptr(3,np+j)*ptr(2,np+k)
 1788 continue
      do 1801 i=1,10
      go to (1790,1791,1792,1793,1794,1795,1796,1797,1798,1799),i
 1790 j=1
      k=1
      go to 1800
 1791 j=2
      k=2
      go to 1800
 1792 j=3
      k=3
      go to 1800
 1793 j=1
      k=2
      go to 1800
 1794 j=1
      k=3
      go to 1800
 1795 j=2
      k=1
      go to 1800
 1796 j=2
      k=3
      go to 1800
 1797 j=3
      k=1
      go to 1800
 1798 j=3
      k=2
      go to 1800
 1799 j=4
      k=3
 1800 ftr(1,nf+i)=dtr(1,nd+j)*ptr(1,np+k)
      ftr(2,nf+i)=dtr(2,nd+j)*ptr(2,np+k)
      ftr(3,nf+i)=dtr(3,nd+j)*ptr(3,np+k)
      ftr(4,nf+i)=dtr(1,nd+j)*ptr(2,np+k)
     1           +dtr(4,nd+j)*ptr(1,np+k)
      ftr(5,nf+i)=dtr(1,nd+j)*ptr(3,np+k)
     1           +dtr(5,nd+j)*ptr(1,np+k)
      ftr(6,nf+i)=dtr(2,nd+j)*ptr(1,np+k)
     1           +dtr(4,nd+j)*ptr(2,np+k)
      ftr(7,nf+i)=dtr(2,nd+j)*ptr(3,np+k)
     1           +dtr(6,nd+j)*ptr(2,np+k)
      ftr(8,nf+i)=dtr(3,nd+j)*ptr(1,np+k)
     1           +dtr(5,nd+j)*ptr(3,np+k)
      ftr(9,nf+i)=dtr(3,nd+j)*ptr(2,np+k)
     1           +dtr(6,nd+j)*ptr(3,np+k)
      ftr(10,nf+i)=dtr(4,nd+j)*ptr(3,np+k)
     1           + dtr(5,nd+j)*ptr(2,np+k)
     2           + dtr(6,nd+j)*ptr(1,np+k)
 1801 continue
      do 1819 i=1,15
      go to (1803,1804,1805,1806,1807,1808,1809,1810,1811,1812,
     1       1813,1814,1815,1816,1817),i
 1803 j=1
      k=1
      go to 1818
 1804 j=2
      k=2
      go to 1818
 1805 j=3
      k=3
      go to 1818
 1806 j=1
      k=2
      go to 1818
 1807 j=1
      k=3
      go to 1818
 1808 j=2
      k=1
      go to 1818
 1809 j=2
      k=3
      go to 1818
 1810 j=3
      k=1
      go to 1818
 1811 j=3
      k=2
      go to 1818
 1812 j=4
      k=2
      go to 1818
 1813 j=5
      k=3
      go to 1818
 1814 j=7
      k=3
      go to 1818
 1815 j=4
      k=3
      go to 1818
 1816 j=6
      k=3
      go to 1818
 1817 j=8
      k=2
 1818 gtr(1,ng+i)=ftr(1,nf+j)*ptr(1,np+k)
      gtr(2,ng+i)=ftr(2,nf+j)*ptr(2,np+k)
      gtr(3,ng+i)=ftr(3,nf+j)*ptr(3,np+k)
      gtr(4,ng+i)=ftr(1,nf+j)*ptr(2,np+k)
     1           +ftr(4,nf+j)*ptr(1,np+k)
      gtr(5,ng+i)=ftr(1,nf+j)*ptr(3,np+k)
     1           +ftr(5,nf+j)*ptr(1,np+k)
      gtr(6,ng+i)=ftr(2,nf+j)*ptr(1,np+k)
     1           +ftr(6,nf+j)*ptr(2,np+k)
      gtr(7,ng+i)=ftr(2,nf+j)*ptr(3,np+k)
     1           +ftr(7,nf+j)*ptr(2,np+k)
      gtr(8,ng+i)=ftr(3,nf+j)*ptr(1,np+k)
     1           +ftr(8,nf+j)*ptr(3,np+k)
      gtr(9,ng+i)=ftr(3,nf+j)*ptr(2,np+k)
     1           +ftr(9,nf+j)*ptr(3,np+k)
      gtr(10,ng+i)=ftr(4,nf+j)*ptr(2,np+k)
     1           + ftr(6,nf+j)*ptr(1,np+k)
      gtr(11,ng+i)=ftr(5,nf+j)*ptr(3,np+k)
     1           + ftr(8,nf+j)*ptr(1,np+k)
      gtr(12,ng+i)=ftr(7,nf+j)*ptr(3,np+k)
     1           + ftr(9,nf+j)*ptr(2,np+k)
      gtr(13,ng+i)=ftr(4,nf+j)*ptr(3,np+k)
     1           + ftr(5,nf+j)*ptr(2,np+k)
     2           +ftr(10,nf+j)*ptr(1,np+k)
      gtr(14,ng+i)=ftr(6,nf+j)*ptr(3,np+k)
     1           + ftr(7,nf+j)*ptr(1,np+k)
     2           +ftr(10,nf+j)*ptr(2,np+k)
      gtr(15,ng+i)=ftr(8,nf+j)*ptr(2,np+k)
     1           + ftr(9,nf+j)*ptr(1,np+k)
     2           +ftr(10,nf+j)*ptr(3,np+k)
 1819 continue
 1820 continue
      if(normf.eq.1.and.normp.eq.1) go to 1831
      do 1830 it=1,nt
      nd=6*(it-1)
      nf=10*(it-1)
      ng=15*(it-1)
      do 1822 i=1,6
      if(i.gt.3) go to 1821
      dtr(4,nd+i)=dtr(4,nd+i)/sqrt3
      dtr(5,nd+i)=dtr(5,nd+i)/sqrt3
      dtr(6,nd+i)=dtr(6,nd+i)/sqrt3
      go to 1822
 1821 dtr(1,nd+i)=dtr(1,nd+i)*sqrt3
      dtr(2,nd+i)=dtr(2,nd+i)*sqrt3
      dtr(3,nd+i)=dtr(3,nd+i)*sqrt3
 1822 continue
      do 1825 i=1,10
      if(i.gt.3) go to 1823
      ftr(4,nf+i)=ftr(4,nf+i)/sqrt5
      ftr(5,nf+i)=ftr(5,nf+i)/sqrt5
      ftr(6,nf+i)=ftr(6,nf+i)/sqrt5
      ftr(7,nf+i)=ftr(7,nf+i)/sqrt5
      ftr(8,nf+i)=ftr(8,nf+i)/sqrt5
      ftr(9,nf+i)=ftr(9,nf+i)/sqrt5
      ftr(10,nf+i)=ftr(10,nf+i)/(sqrt5*sqrt3)
      go to 1825
 1823 if(i.gt.9) go to 1824
      ftr(1,nf+i)=ftr(1,nf+i)*sqrt5
      ftr(2,nf+i)=ftr(2,nf+i)*sqrt5
      ftr(3,nf+i)=ftr(3,nf+i)*sqrt5
      ftr(10,nf+i)=ftr(10,nf+i)/sqrt3
      go to 1825
 1824 ftr(1,nf+i)=ftr(1,nf+i)*sqrt5*sqrt3
      ftr(2,nf+i)=ftr(2,nf+i)*sqrt5*sqrt3
      ftr(3,nf+i)=ftr(3,nf+i)*sqrt5*sqrt3
      ftr(4,nf+i)=ftr(4,nf+i)*sqrt3
      ftr(5,nf+i)=ftr(5,nf+i)*sqrt3
      ftr(6,nf+i)=ftr(6,nf+i)*sqrt3
      ftr(7,nf+i)=ftr(7,nf+i)*sqrt3
      ftr(8,nf+i)=ftr(8,nf+i)*sqrt3
      ftr(9,nf+i)=ftr(9,nf+i)*sqrt3
 1825 continue
      do 1829 i=1,15
      if(i.gt.3) go to 1826
      gtr(4,ng+i)=gtr(4,ng+i)/sqrt7
      gtr(5,ng+i)=gtr(5,ng+i)/sqrt7
      gtr(6,ng+i)=gtr(6,ng+i)/sqrt7
      gtr(7,ng+i)=gtr(7,ng+i)/sqrt7
      gtr(8,ng+i)=gtr(8,ng+i)/sqrt7
      gtr(9,ng+i)=gtr(9,ng+i)/sqrt7
      gtr(10,ng+i)=gtr(10,ng+i)*sqrt3/(sqrt5*sqrt7)
      gtr(11,ng+i)=gtr(11,ng+i)*sqrt3/(sqrt5*sqrt7)
      gtr(12,ng+i)=gtr(12,ng+i)*sqrt3/(sqrt5*sqrt7)
      gtr(13,ng+i)=gtr(13,ng+i)/(sqrt5*sqrt7)
      gtr(14,ng+i)=gtr(14,ng+i)/(sqrt5*sqrt7)
      gtr(15,ng+i)=gtr(15,ng+i)/(sqrt5*sqrt7)
      go to 1829
 1826 if(i.gt.9) go to 1827
      gtr(1,ng+i)=gtr(1,ng+i)*sqrt7
      gtr(2,ng+i)=gtr(2,ng+i)*sqrt7
      gtr(3,ng+i)=gtr(3,ng+i)*sqrt7
      gtr(10,ng+i)=gtr(10,ng+i)*sqrt3/sqrt5
      gtr(11,ng+i)=gtr(11,ng+i)*sqrt3/sqrt5
      gtr(12,ng+i)=gtr(12,ng+i)*sqrt3/sqrt5
      gtr(13,ng+i)=gtr(13,ng+i)/sqrt5
      gtr(14,ng+i)=gtr(14,ng+i)/sqrt5
      gtr(15,ng+i)=gtr(15,ng+i)/sqrt5
      go to 1829
 1827 if(i.gt.12) go to 1828
      gtr(1,ng+i)=gtr(1,ng+i)*sqrt7*sqrt5/sqrt3
      gtr(2,ng+i)=gtr(2,ng+i)*sqrt7*sqrt5/sqrt3
      gtr(3,ng+i)=gtr(3,ng+i)*sqrt7*sqrt5/sqrt3
      gtr(4,ng+i)=gtr(4,ng+i)*sqrt5/sqrt3
      gtr(5,ng+i)=gtr(5,ng+i)*sqrt5/sqrt3
      gtr(6,ng+i)=gtr(6,ng+i)*sqrt5/sqrt3
      gtr(7,ng+i)=gtr(7,ng+i)*sqrt5/sqrt3
      gtr(8,ng+i)=gtr(8,ng+i)*sqrt5/sqrt3
      gtr(9,ng+i)=gtr(9,ng+i)*sqrt5/sqrt3
      gtr(13,ng+i)=gtr(13,ng+i)/sqrt3
      gtr(14,ng+i)=gtr(14,ng+i)/sqrt3
      gtr(15,ng+i)=gtr(15,ng+i)/sqrt3
      go to 1829
 1828 gtr(1,ng+i)=gtr(1,ng+i)*sqrt7*sqrt5
      gtr(2,ng+i)=gtr(2,ng+i)*sqrt7*sqrt5
      gtr(3,ng+i)=gtr(3,ng+i)*sqrt7*sqrt5
      gtr(4,ng+i)=gtr(4,ng+i)*sqrt5
      gtr(5,ng+i)=gtr(5,ng+i)*sqrt5
      gtr(6,ng+i)=gtr(6,ng+i)*sqrt5
      gtr(7,ng+i)=gtr(7,ng+i)*sqrt5
      gtr(8,ng+i)=gtr(8,ng+i)*sqrt5
      gtr(9,ng+i)=gtr(9,ng+i)*sqrt5
      gtr(10,ng+i)=gtr(10,ng+i)*sqrt3
      gtr(11,ng+i)=gtr(11,ng+i)*sqrt3
      gtr(12,ng+i)=gtr(12,ng+i)*sqrt3
 1829 continue
 1830 continue
 1831 continue
c
c     iop1=5   ... print convergence data for each cycle
c                  + extensive printing
c     iop1=4   ... telex output
c     iop1=3   ... extensive printing
c     iop1=2   ... extensive printing of basis set only
c     iop1=1   ... normal printing + 2e-and 1e-integrals
c     iop1=0   ... normal printing
c
      if(redux) call reduc
      do 1840 i=5,35
      dfg(i)=dfg1(i)
      if(redux) dfg(i)=dfg2(i)
 1840 continue
      if(iop1.eq.4) go to 2500
      write(iw,8884)
      do 1851 iat=1,nat
      write(iw,8883) a(iat),apseud(iat),iatno(iat),zan(iat),c(1,iat),
     1 c(2,iat),c(3,iat),ns(iat)
      if(ns(iat).eq.0) go to 1851
      ns1=ks(iat)
      ns2=ns1+ns(iat)-1
      do 1850 ish=ns1,ns2
 1850 ind(ish)=ktype(ish)
      write(iw,8871) (label(ind(ish)),ish=ns1,ns2)
 1851 continue
      write(iw,8882)
      do 1860 ish=1,nshell
      write(iw,8855)
      i1=kstart(ish)
      i2=i1+kng(ish)-1
      ityp=ktype(ish)
      do 1859 ig=i1,i2
      go to (1852,1853,1854,1855,1856,1852,1857,1854,1855,1856),ityp
 1852 c1=cs(ig)
      c2=csinp(ig)
      go to 1858
 1853 c1=cp(ig)
      c2=cpinp(ig)
      go to 1858
 1854 c1=cd(ig)
      c2=cdinp(ig)
      go to 1858
 1855 c1=cf(ig)
      c2=cfinp(ig)
      go to 1858
 1856 c1=cg(ig)
      c2=cginp(ig)
      go to 1858
 1857 c1=cs(ig)
      c2=csinp(ig)
      c3=cp(ig)
      c4=cpinp(ig)
      write(iw,8849) ish,label(ityp),ig,ex(ig),c1,c2,c3,c4
      go to 1859
 1858 write(iw,8881) ish,label(ityp),ig,ex(ig),c1,c2
 1859 continue
 1860 continue
      if(iop1.lt.2) go to 2500
      write(iw,8860)
      write(iw,8870)
      imax=0
 1900 imin=imax+1
      imax=imax+15
      if(imax.gt.nt) imax=nt
      imax1=imax+1
      write(iw,8869)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8867) linev,(i,i=imin,imax)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (line ,i=imin,imax1)
      do 1950 iat=1,nat
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8866) iat,(ict(iat,i),i=imin,imax)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (line ,i=imin,imax1)
 1950 continue
      if(imax.lt.nt) go to 1900
      write(iw,8860)
      write(iw,8865)
      imax=0
 2000 imin=imax+1
      imax=imax+15
      if(imax.gt.nt) imax=nt
      imax1=imax+1
      write(iw,8869)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8867) linev,(i,i=imin,imax)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (line ,i=imin,imax1)
      do 2050 ish=1,nshell
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8866) ish,(newsh(ish,i),i=imin,imax)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (line ,i=imin,imax1)
 2050 continue
      if(imax.lt.nt) go to 2000
      write(iw,8860)
      write(iw,8851)
      imax=0
 2060 imin=imax+1
      imax=imax+15
      if(imax.gt.nt) imax=nt
      imax1=imax+1
      write(iw,8869)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8867) linev,(i,i=imin,imax)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (line ,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8867) linev,(invt(i),i=imin,imax)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (linev,i=imin,imax1)
      write(iw,8868) (line ,i=imin,imax1)
      if(imax.lt.nt) go to 2060
      if(iop1.ne.3) go to 2500
      write(iw,8863)
      do 2155 it=1,nt
      write(iw,8860)
      write(iw,8858) it
      np=3*(it-1)
      write(iw,8862) (pname(j),j=1,3)
      write(iw,8855)
      do 2150 i=1,3
 2150 write(iw,8861) pname(i),(ptr(i,np+j),j=1,3)
      write(iw,8869)
      nd=6*(it-1)
      ddmax=6
      if(redux)ddmax=5
      write(iw,8862) (dname(j),j=1,ddmax)
      write(iw,8855)
      do 2151 i=1,ddmax
 2151 write(iw,8861) dname(i),(dtr(i,nd+j),j=1,ddmax)
      write(iw,8869)
      nf=10*(it-1)
      ffmax=10
      if(redux) ffmax=7
      write(iw,8862) (fname(j),j=1,ffmax)
      write(iw,8855)
      do 2152 i=1,ffmax
 2152 write(iw,8861) fname(i),(ftr(i,nf+j),j=1,ffmax)
      write(iw,8869)
      ng=15*(it-1)
      ggmax=15
      if(redux) ggmax=9
      jmax=0
 2153 jmin=jmax+1
      jmax=jmax+ggmax
      if(jmax.gt.ggmax) jmax=ggmax
      write(iw,8862) (gname(j),j=jmin,jmax)
      write(iw,8855)
      do 2154 i=1,ggmax
 2154 write(iw,8861) gname(i),(gtr(i,ng+j),j=jmin,jmax)
      write(iw,8855)
      if(jmax.lt.ggmax) go to 2153
      write(iw,8869)
 2155 continue
 2500 continue
      if(nat.le.dc) go to 2510
      write(iw,8876)
      stop
 2510 if(nshell.le.ds) go to 2520
      write(iw,8875)
      stop
c2520 if(loc.le.80) go to 2530
 2520 if(loc.le.doa) go to 2530
      write(iw,8874)
      stop
 2530 if(ngauss.le.dgp) go to 2540
      write(iw,8873)
      stop
 6020 write(iw,8848)
      stop
 2540 continue
      num=loc
      ne=ne-ich
      na=(ne+mul)/2
      nb=(ne-mul)/2
      nx=num*(num+1)/2
      write(iw,8877) nshell,num,ne,ich,mul,na,nb,nat
      do 5200 i=1,nshell
      ic=katom(i)
      min=kmin(i)
      max=kmaz(i)
      iatn=iatno(ic)
      loc=kloc(i)-min
      do 5200 k=min,max
      lck=loc+k
      iflab(1,lck)=blk//atnam(iatn)//blk
      iflab(2,lck)=number(ic)//blk//blk
 5200 iflab(3,lck)=dfg(k)
      return
      end
