      subroutine poltwo
      implicit real*8(a-h,o-z)
      include 'pshf.prm'                                                   
      logical iandj,kandl,same
      logical out,ymono
      parameter(EPSI=1.D-7)
      real*8,dimension(:),allocatable :: xx_pqrs
      integer,dimension(:),allocatable :: ix_pqrs        
      common/infnp/nprint
      common/iofile/ir,iw,ip,is,iq
      common/shlt/ tol,cutoff,icount,out
      common/polar/calfa(dc)
      common/fil2/ a(dc),ns(dc),ks(dc),newsh(ds,48),ngauss
      common/symtry/invt(48),iso(ds,12),ptr(3,144),dtr(6,288),
     1 ftr(10,480),nt
      common/nshel/exg(dgp),cs(dgp),cp(dgp),cd(dgp),
     1 cf(dgp),cg(dgp),nshell,kstart(ds),katom(ds),
     2 ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     3 kmaz(ds)
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(dc),c(3,dc),
     1 ymono,mrec

      dimension m0(48),m1(48),m2(48),m3(48)
      dimension mm0(48),mm1(48),mm2(48),mm3(48)                         
      dimension ex(doas,nblk),ey(doas,nblk),ez(doas,nblk)
      dimension ix(lsize),xx(lsize)
      dimension ixout(lsize),xxout(lsize)
      dimension ish(doa),ia(doa)
      dimension listat(dc+1)

      do 5 i=1,doa
      ia(i)=(i*(i-1))/2
    5 continue
      write(6,*) 'dimension des buffers dans poltwo', lsize
      write(6,*) 'nombre de blocs', mrec
      out=nprint.eq.1
      ntwd=(nt+3)/4
      tol=1.d-12
      cutoff=1.d-09
      iin=is
      iout=iq
      rewind 17
      do 1 i=1,nshell
      lit=ktype(i)
      mini=kmin(i)
      maxi=kmaz(i)
      loci=kloc(i)-mini
      do 2 j=mini,maxi
   2  ish(loci+j)=i
   1  continue
      write (iw,*) ' polarisabilites atomiques des coeurs'
      write (iw,*) (calfa(i),i=1,nat)
c     boucle sur les atomes
C M. GROSS 16/6/97  par blocks de nblk atomes
C compte les atomes polarisables
      npat=0
      do i=1,nat
        if (calfa(i).gt.EPSI) then
           npat= npat+1
           listat(npat) = i
        end if
      end do
C index pour fin de tableau
      listat(npat+1) = 0
      do 15000 iblk=1,(npat-1)/nblk+1
      do ic=1,nblk
       iatom = listat((iblk-1)*nblk+ic)
       if (iatom.eq.0) goto 20
       write (iw,*) 'Block, atome',iblk,iatom,calfa(iatom)
       read(17) (ex(i,ic),i=1,nnp)
       read(17) (ey(i,ic),i=1,nnp)
       read(17) (ez(i,ic),i=1,nnp)
       read(17)
      end do
20    continue

      write(6,*) 'iin: ',iin
      rewind iin
      rewind iout
      nrec=0
      mind=1
      mout=0
      
c     stockage of file iin in memory 
      imin=1
      imax=lsize  

      allocate(xx_pqrs(lsize*(mrec+1)),ix_pqrs(lsize*(mrec+1)))
      do i=1,mrec     
      read(iin)xx_pqrs(imin:imax),ix_pqrs(imin:imax)     
      imin=imin+lsize
      imax=imin+lsize-1
      enddo
      read(iin,end=105)xx_pqrs(imin),ix_pqrs(imin)
105   continue
      imin=1
      imax=lsize
c      rewind iin
      
      xx=xx_pqrs(imin:imax)
      ix=ix_pqrs(imin:imax)
      imin=imin+lsize
      imax=imin+lsize-1
c      read(iin,end=100) xx,ix
      do 13000 in=1,nshell
      do 1020 it=1,nt
      mm0(it)=newsh(in,it)
      if (mm0(it).gt.in) go to 13000
 1020 continue
c
      do 12000 jn=1,in
      do 1030 it=1,nt
      mm1(it)=newsh(jn,it)
      if (mm1(it).gt.jn) go to 12000
 1030 continue
c
      do 11000 kn=1,in
      do 1040 it=1,nt
      mm2(it)=newsh(kn,it)
      if (mm2(it).gt.kn) go to 11000
 1040 continue
c
      maxln=kn
      if (in.eq.kn) maxln=jn
c
      do 10000 ln=1,maxln
      do 1050 it=1,nt
      mm3(it)=newsh(ln,it)
      if (mm3(it).gt.ln) go to 10000
 1050 continue
c
c      ......boucle sur les quadruplets de transf. de sym.
c
c     ----- ishell -----                                                
c                                                                       
      do 9000 iop=1,nt                                                  
      iii=mm0(iop)
      ip=iop
 1100 ip=ip+1
      if (ip.gt.nt) go to 1110
      if (iii.eq.mm0(ip)) go to 9000
      go to 1100
c                                                                       
c     ----- jshell -----                                                
c                                                                       
 1110 do 8000 jop=1,nt                                                  
      jjj=mm1(jop)
      jp=jop
 1120 jp=jp+1
      if (jp.gt.nt) go to 1130
      if (jjj.eq.mm1(jp)) go to 8000
      go to 1120
c
 1130 iii=mm0(iop)
      if (iii.ge.jjj) go to 1140
      if (in.eq.jn) go to 8000
      nd=iii
      iii=jjj
      jjj=nd
c                                                                       
c     ----- kshell -----                                                
c                                                                       
 1140 do 7000 kop=1,nt                                                  
      kk=mm2(kop)
      kp=kop
 1150 kp=kp+1
      if (kp.gt.nt) go to 1160
      if (kk.eq.mm2(kp)) go to 7000
      go to 1150
c                                                                       
c     ----- lshell -----                                                
c                                                                       
 1160 do 6000 lop=1,nt                                                  
      ll=mm3(lop)
      lp=lop
 1170 lp=lp+1
      if (lp.gt.nt) go to 1200
      if (ll.eq.mm3(lp)) go to 6000
      go to 1170
 1200 ii=iii
c
      jj=jjj
      kk=mm2(kop)
c
      if (kk.ge.ll) go to 1210
      if (kn.eq.ln) go to 6000
      nd=ll
      ll=kk
      kk=nd
 1210 if (ii.gt.kk) go to 1240
      if (ii.eq.kk) go to 1220
      nd=ii
      ii=kk
      kk=nd
      go to 1230
 1220 if (jj.ge.ll) go to 1240
 1230 if (in.eq.kn.and.jn.eq.ln) go to 6000
      nd=jj
      jj=ll
      ll=nd
c
 1240 continue
      do  it=1,nt
      m0(it)=newsh(ii,it)
      m1(it)=newsh(jj,it)
      m2(it)=newsh(kk,it)
      m3(it)=newsh(ll,it)
      end do
c
      n4=0
      do 1570 it=1,nt
      ld=m3(it)
      kd=m2(it)
      if (kd.ge.ld) go to 1520
      nd=kd
      kd=ld
      ld=nd
 1520 id=m0(it)
      jd=m1(it)
      if (id.ge.jd) go to 1530
      nd=id
      id=jd
      jd=nd
 1530 if (id.gt.kd) go to 1560
      if (id.eq.kd) go to 1540
      nd=id
      id=kd
      kd=nd
      go to 1550
 1540 if (jd.ge.ld) go to 1560
 1550 nd=jd
      jd=ld
      ld=nd
 1560 if (id.lt.ii) go to 1570
      if (id.gt.ii) go to 6000
      if (jd.lt.jj) go to 1570
      if (jd.gt.jj) go to 6000
      if (kd.lt.kk) go to 1570
      if (kd.gt.kk) go to 6000
      if (ld.lt.ll) go to 1570
      if (ld.gt.ll) go to 6000
      n4=n4+1
 1570 continue
c     couches en ordre
      io=ii
      jo=jj
      if(ktype(io).lt.ktype(jo)) then
      n=io
      io=jo
      jo=n
      endif
      ko=kk
      lo=ll
      if(ktype(ko).lt.ktype(lo)) then
      n=ko
      ko=lo
      lo=n
      endif
c
c
c     ----- calculate q4 factor for this group of shells
c
      q4=dfloat(nt)/dfloat(n4)
c
      iandj=ii.eq.jj
      kandl=kk.eq.ll
      same=ii.eq.kk.and.jj.eq.ll
      ijoa=0
      do 5500 ioa=kloc(io),kloc(io)+kmaz(io)-kmin(io)
      joamax=kloc(jo)+kmaz(jo)-kmin(jo)
      if(iandj) joamax=ioa
      do 5400 joa=kloc(jo),joamax
      ijoa=ijoa+1
      kloa=0
      do 5300 koa=kloc(ko),kloc(ko)+kmaz(ko)-kmin(ko)
      loamax=kloc(lo)+kmaz(lo)-kmin(lo)
      if(kandl) loamax=koa
      do 5200 loa=kloc(lo),loamax
      kloa=kloa+1
      if(same.and.kloa.gt.ijoa) go to 5200
      qq4=q4
c     ordre canonique
      i1=ioa
      i2=joa
      i3=koa
      i4=loa
      if(i1.ge.i2) go to 700
      n=i1
      i1=i2
      i2=n
  700 if(i3.ge.i4) go to 800
      n=i3
      i3=i4
      i4=n
  800 if(i1-i3) 900,1000,1111
  900 n=i1
      i1=i3
      i3=n
      n=i2
      i2=i4
      i4=n
      go to 1111
 1000 if(i2.lt.i4) go to 900
 1111 continue
      if(i1.eq.i2)qq4=qq4/2.d0
      if(i3.eq.i4)qq4=qq4/2.d0
      if(i1.eq.i3.and.i2.eq.i4) qq4=qq4/2.d0
      i12=ia(i1)+i2
      i34=ia(i3)+i4
      label=iword(i1,i2,i3,i4)
      dum1 = 0.0
      do ic=1,nblk
      iatom = listat((iblk-1)*nblk+ic)
      if (iatom.eq.0) goto 1300
        dum1=dum1-(ex(i12,ic)*ex(i34,ic)+ey(i12,ic)*ey(i34,ic)+
     &             ez(i12,ic)*ez(i34,ic))*calfa(iatom)*qq4
      end do
1300  dum=dum1
      labl2=ix(mind)
      if(label.eq.labl2) then
        xxm=xx(mind)
        dum=xxm+dum1
        mind=mind+1
        if(mind.gt.lsize) then
        xx=xx_pqrs(imin:imax)
        ix=ix_pqrs(imin:imax)
        imin=imin+lsize
        imax=imin+lsize-1
c        read(iin) xx,ix        
          mind=1
        endif
      endif
      if(dabs(dum).gt.cutoff) then
c      write(6,*) ic,i1,i2,i3,i4,label,labl2,xxm,dum1,dum,mout,mind
      mout=mout+1
      xxout(mout)=dum
      ixout(mout)=label
             if(mout.eq.lsize) then
             write(iout) xxout,ixout
             nrec=nrec+1
             mout=0
             endif
      endif
 5200 continue
 5300 continue
 5400 continue
 5500 continue
 6000 continue
 7000 continue
 8000 continue
 9000 continue
10000 continue
11000 continue
12000 continue
13000 continue
c
      ixout(mout+1)=0
      write(iout) xxout,ixout
      nrec=nrec+1
      nd1=iin
      nd2=iout
      nd=iin
      iin=iout
      iout=nd
      write (6,*) 'ic is iq iin iout nd1 nd2',ic,is,iq,iin,iout,nd1,nd2
c      write(6,*) (ixout(j),xxout(j),j=1,500)
15000 continue
c fin de boucle sur les atomes
c
c
c     verifier sur quelle unite se trouvent les integrales
c     attention dans cette version, on garde la file definitive sur f09
c
      if (iout.eq.is) go to 1015
      rewind is
      rewind iq
      do 1010 i=1,nrec
      read(is) xx,ix
 1010 write(iq) xx,ix
 1015 rewind is
      rewind iq
100   continue
      if(nrec.eq.0) write(6,*) 'attention file pqrs vide'
      write(6,'(a27)')'nombre de blocs sur la file'
      write(6,*) iq,nrec
      
      deallocate(xx_pqrs,ix_pqrs)
      return
      end
c     RV 01/16 rewrite iword function in fortran 
      function iword(iarg1,iarg2,iarg3,iarg4)
        implicit none
        integer :: iword
        integer :: iarg1,iarg2,iarg3,iarg4
        integer :: s1,s2,s3,s4
        
        s1 = ishft(iarg1,24)         
        s2 = ishft(iarg2,16)         
        s3 = ishft(iarg3,8)         
        s4 = ishft(iarg4,0)                                      
        iword= ior(s1,ior(s2,ior(s3,s4)))                              
        return                                         
      end function