      subroutine pert(f,v,d,e,ia,in,ipass,ending)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      logical ending,scfp
      logical convg
      character*8 aname,bname
      dimension igel(doa),g(doa),h(doa)
      common/iofile/ ir,iw,ip,is,iq,ih,iv
      common/scfit/ aname,bname,maxit,nconv,npunch,npi,cextrp,amix
      common/output/ nprint
      common/conv/acurcy,ehf,ehf0,iter
      dimension nam(10),gg(doa,5)
      equivalence (gg(1,1),g(1)),(gg(1,2),h(1))
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc)
      dimension f(*),v(*),d(*),e(*),ia(*),in(*)
      common/pseudo/qa(doa),qb(doa),ngla(doa),nglb(doa),numgla,numglb,
     * scfp
c
      data thresh,delmax/1.d-8,0.2/
 9999 format(/5x,' the dimension of the basis equals the number of occup
     *ied orbitals .... program stops')
 9998 format(10x,22(1h-),/,10x,'fock matrix(m.o. basis)',/,10x,22(1h-))
 9997 format(10x,22(1h-),/,10x,'fock matrix(after mixing)',/,10x,22(1h-)
     *)
 9994 format(/'  doubly occupied orbitals =',i3,
     *      /,'  singly occupied orbitals =',i3,
     *      /,'          virtual orbitals =',i3)
 9993 format(' maximum perturbation coefficient =',2i4,2d15.8)
c
      nxsq=numscf*numscf
      alp=amix
c      delmax=delmax*alp
c      if(iter.ne.1.or.bname.ne.'rhf') go to 3
      nda=0
      nvt=0
      nos=0
      if(bname.ne.'uhf') then
      do 2 i=1,numscf
      if(qa(i).eq.2.d0) nda=nda+1
      if(qa(i).eq.0.d0) nvt=nvt+1
      if(qa(i).eq.1.d0) nos=nos+1
      if(qa(i).eq.1.d0) nam(nos)=i
    2 continue
      if(iter.eq.1)write(iw,9994) nda,nos,nvt
      noc=nda
      else
      noc=na
      if(ipass.eq.2)noc=nb
      end if
    3 continue
      if(na.lt.numscf) go to 5
      write(iw,9999)
      stop
    5 continue
c
      do 8 i=1,numscf
    8 igel(i)=1
      if(ipass.eq.1) then
      nav1=5
      do 6 i=1,numgla
      igel(ngla(i))=0
 6    continue
      else
      nav1=8
      noc=nb
      do 7 i=1,numglb
      igel(nglb(i))=0
 7    continue
      end if
c
c
c     calculation of molecular fock matrix
c
c     get ao fock matrix
      call reada(f,nx,nav1)
c     get vectors
      call reada(v,nxsq,nav1+2)
      do 550 i=1,numscf
      ni=in(i)
      do 532 j=1,numscf
  532 g(j)=0.d0
      do 540 k=1,numscf
      dum=v(ni+k)
      if(dum.eq.0.d0) go to 540
      nk=ia(k)
      do 534 j=1,k
  534 g(j)=g(j)+f(nk+j)*dum
      if(k.eq.numscf) go to 540
      kp=k+1
      do 538 j=kp,numscf
  538 g(j)=g(j)+f(ia(j)+k)*dum
  540 continue
      do 545 j=i,numscf
      h(j)=0.d0
      nj=in(j)
      do 545 k=1,numscf
      dum=v(nj+k)
      if(dum.eq.0.d0) go to 545
      h(j)=h(j)+dum*g(k)
  545 continue
      do 550 j=i,numscf
  550 v(ni+j)=h(j)
      do 560 i=1,numscf
      do 560 j=i,numscf
      if(abs(v(in(i)+j)).lt.1.d-8) v(in(i)+j)=0.d0
      if(i.eq.j) e(i)=v(in(i)+j)
  560 f(ia(j)+i)=v(in(i)+j)
      if(nprint.lt.5) go to 590
      write(iw,9998)
      call fout(f,numscf)
  590 continue
      if(nos.eq.0) go to 670
      rewind 3
      write(3) (f(i),i=1,nx)
      rewind ih
      read(ih)
      read(ih) (v(i),i=1,nxsq)
      do 601 i=1,nx
  601 f(i)=0.d0
      do 600 ios=1,nos
      i1=nam(ios)
      do 600 i=1,numscf
      do 600 j=1,i
      dum      =  v(in(i1)+i)*v(in(i1)+j)
  600 f(ia(i)+j)=dum+f(ia(i)+j)
      rewind is
      call hstar(.false.)
      call symh
      do 605 i=1,nx
  605 f(i)=v(i)
      rewind ih
      read(ih)
      read(ih) (v(i),i=1,nxsq)
      do 610 ios=1,nos
      i1=nam(ios)
      do 610 i=1,numscf
      gg(i,ios)=0.d0
      do 610 j=1,numscf
      do 610 k=1,numscf
      jk=ia(j)+k
      if(j.lt.k) jk=ia(k)+j
      gg(i,ios)=gg(i,ios)+v(in(i)+k)*v(in(i1)+j)*f(jk)
  610 continue
      rewind 3
      read(3) (f(i),i=1,nx)
      do 620 ios=1,nos
      i1=nam(ios)
      do 620 i=1,numscf
      ii1=ia(i)+i1
      if(i.lt.i1) ii1=ia(i1)+i
      if(i.le.nda) go to 615
      f(ii1)=f(ii1)+gg(i,ios)
      go to 620
  615 f(ii1)=f(ii1)-gg(i,ios)
  620 continue
  650 continue
c      correction a l'energie
      dum=0.d0
      do 660 ios=1,nos
      i1=nam(ios)
  660 dum=dum+gg(i1,ios)
      ehf=ehf+0.5d0*dum
  670 continue
c
c     first order perturbation
c
      call reada(v,nxsq,nav1+2)
      delmm=0.d0
c      write(6,*) 'noc',noc,(igel(i),i=1,numscf)
      do 10 i=1,noc
      if(igel(i).eq.0) go to 10
      ei=e(i)
      do 11 j=noc+1,numscf
      if(igel(j).eq.0) go to 11
c      if(i.gt.nda.and.j.le.na) go to 11
      delta=f(ia(j)+i)/(ei-e(j))
c      write(6,*) 'i,j,f(ia(j)+i),e(i),e(j),delta,delmm',
c     8i,j,f(ia(j)+i),e(i),e(j),delta,delmm
      if(abs(delta).le.delmm) go to 9
      delmm=abs(delta)
      imax=i
      jmax=j
    9 f(ia(j)+i)=delta
   11 continue
   10 continue
c      call fout(f,numscf)
      write(iw,9993) imax,jmax,delmm,delmax
      if(delmm.lt.acurcy) ending=.true.
      do 212 i=1,numscf
  212 e(i)=f(ia(i)+i)
      red=1.d0
      if(delmm.gt.delmax) then
      red=delmax/delmm
      do 210 i=1,nx
  210 f(i)=f(i)*red
c      write(6,*) 'red',red
      end if
c      if(nprint. eq.5) then
c      call fout(f,numscf)
c      end if
c
      do 16 k=1,numscf
      do 12 i=1,numscf
   12 g(i)=v(in(i)+k)
      do 15 i=1,noc
      if(igel(i).eq.0) go to 15
      do 14 j=noc+1,numscf
c      if(i.gt.nda.and.j.le.na) go to 14
      if(igel(j).eq.0) go to 14
      fij=f(ia(j)+i)
      g(i)=g(i)+fij*v(in(j)+k)
      g(j)=g(j)-fij*v(in(i)+k)
   14 continue
   15 continue
      do 16 i=1,numscf
   16 v(in(i)+k)=g(i)
c
      if(nprint.eq.5) call vout(v,in,numscf,numscf)
      if(ipass.eq.1)call ortit(f,v,g,h,ia,in,ngla,numgla,nav1+2)
      if(ipass.eq.2)call ortit(f,v,g,h,ia,in,nglb,numglb,nav1+2)
      return
      end
