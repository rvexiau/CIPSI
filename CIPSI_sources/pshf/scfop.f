      subroutine scfop
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical backt,ending,scfp,uhf
      character*8 aname,bname,type
c
c     -----              scf calculations -----
c               bname='rhf'   restricted scf calculation
c               bname='osrhf' restricted open shell scf calculation
c     on tape (is)       two-electron integrals
c     on file (ih)   -1- general informations(basis set,geometry,..)
c                    -2- overlap
c                    -3- one electron integral matrix
c                    -4- fock transformation matrix  s**-1/2
c                    -5- closed shell fock matrix, r operator
c                    -6- total density matrix
c                    -7- molecular orbitals
c                    -8- a closed shell matrix
c                    -9- b closed shell matrix
c                    -10- shell 2 fock matrix
c                    -11- shell 2 a matrix
c                    -12- shell 2 b matrix
c                    -13- shell 3 fock matrix
c                    -14- shell 3 a matrix
c                    -15- shell 3 b matrix
c                    -16- closed shell projector (c+c)
c                    -17- shell 2      "
c                    -18- shell 3      "
c                    -19- virtual mos  "
c
c     if aname .eq. hcore , a set of starting orbitals will
c     be generated from the bare nucleus hamiltonian matrix
c     if aname .eq. vector , a set of trial vector are readed
c     from cards
c
c     irest = 7 ..... start scf
c     irest = 8 ..... restart scf
c
c     nprint = 5   ... convergence data printed after each iteration
c              4   ... energy data only printed for telex
c            < 4   ... normal printing after convergenge or at exit
c
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
c
      common/iofile/ir,iw,ip,is,iq,ih,iv
      common/open/nc,noc(3),oc(3),alpha(3,3),beta(3,3),rland
      common/times/ti,tx,tim,to
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc),en
c
      common/scfit/ aname,bname,maxit,nconv,npunch,npi,cextrp,amix
      common/conv/acurcy,ehf,ehf0,iter,idiis,ijump
      common/output/nprint,itol,icut,normf,normp,minp
      common/pseudo/qa(doa),qb(doa),ngla(doa),nglb(doa),numgla,numglb,
     1 scfp
      dimension f1(doas),f2(doas),f3(doas),da1(doas),da2(doas),da3(doas)
      dimension db1(doas),db2(doas),db3(doas),f(doas),h(doas)
      dimension v(doas1),g(doa),t(doa),e(doa),ia(doa),in(doa),q(doa)
      dimension ea(doa),eb(doa),d(doas)
      equivalence (f1(1),f(1)),(da1(1),d(1))
      dimension type(3)
      data type/'closed s','open s 1','open s 2'/
      namelist/vecinp/v,ivec,nvec
c
 9998 format(/' cycle ',i3,' energy ',f16.9,' elec. energy',f16.9,
     1 ' time ',3(f6.3,3x))
 9997 format(10x,20(1h-),/,10x,'fock matrix ',a8,/,10x,20(1h-))
 9983 format(10x,18(1h-),/,10x,'r operator matrix',/,10x,18(1h-))
 9984 format(/,10x,26(1h-),/,10x,'fock transformation matrix',/,
     1 10x,26(1h-))
 9996 format(/,10x,12(1h-),/,10x,'eigenvalues',/,10x,12(1h-))
 9995 format(8(i4,f12.6))
 9989 format(i2,i3,5f15.8)
 9994 format(/,10x,13(1h-),/,10x,'eigenvectors',/,10x,13(1h-))
 9993 format(/,10x,15(1h-),/,10x,'density matrix',/,10x,15(1h-))
 9992 format(i4,2x,3f16.9,f10.3)
 9991 format(/10x,16(1h-),/,10x,'energy converged',/,10x,16(1h-))
 9985 format('... scf has converged but program stopped',/,
     1 ' for timlim before ultimate cycle ...')
 9987 format(' ..... scf has not converged ...',
     1 ' ... job to be restarted ....')
 9986 format(/' irest=',i3,' nrec=',i5)
 9990 format(/,10x,30(1h-),/,10x,'excessive number of iterations',
     1 /,10x,30(1h-))
 9988 format(i2,i3,5e15.8)
 9980 format(//,10x,17(1h*),//,10x,' rhf calculation ',//,
     1 10x,17(1h*))
      rewind is
      do 50 i=1,numscf
      in(i)=(i-1)*numscf
   50 ia(i)=i*(i-1)/2
      ehf=99999.0d+00
      ehf0=99999.0d+00
      acurcy=10.0d+00**(-nconv)
c
      uhf=.false.
      iter=0
      idiis=0
      ijump=0
      ending=.false.
      nxsq=numscf*numscf
      nvec=0
      call secnd(tim)
      tim0=tim
      tims=tim
c
c
c     if restart job, go to 300
c
      if(irest.eq.8) go to 300
c
c     if there is a density matrix, go to 300
c
      if(aname.ne.'vector  ') go to 210
c     if there are trial vectors ,get them
c     first pass
      nvec=na
      ivec=0
      ipass=1
      nav1=6
      nav2=7
      do 201 i=1,numscf
  201 q(i)=qa(i)
c
  202 continue
      do 205 i=1,nxsq
  205 v(i)=0.0d0
c     nvec=na(or nb) starting occupied orbitals for standard scf
c     nvec=numscf starting occupied+virtual orbitals for os scf
      read(ir,vecinp,err=3500,end=3500)
      if(ivec.eq.0) go to 208
      nasq=nvec*numscf
      read(ivec) (v(i),i=1,nasq)
      if(ivec.ne.0.and.nvec.eq.numscf) then
         call wrtda(v,nxsq,nav2)
         call dmop(v,f1,da1,db1,da2,db2,da3,db3,in,ia)
         go to 350
      endif
  208 continue
      if(maxit.ne.0) go to 245
c     calculate density matrix with input vectors
c     (eventually non-orthogonal)
      call dmtx(v,f,q,numscf)
      call wrtda(f,nx,nav1)
      call wrtda(v,nxsq,nav2)
      go to 300
 245  continue
c     orthonormalize trial vectors and calculate density matrix
      call orthon (f,v,e,g,ia,in,ngel,q,numgel,nav2,nvec,scfp)
c     store trial density matrix and trial vectors
c
      if(nprint.ne.5) go to 250
      write(iw,9993)
      call fout (f,numscf)
  250 continue
c
      call wrtda(f,nx,nav1)
      call wrtda(v,nxsq,nav2)
c
  210 if(aname.ne.'hcore   ') go to 300
c
c     set h-core matrix to generate initial guess
c
      call reada(f,nx,3)
      call wrtda(f,nx,5)
      go to 510
  300 continue
      call reada(da1,nx,6)
      do 310 i=1,nx
  310 db1(i)=da1(i)/2.d0
      call wrtda(da1,nx,8)
      call wrtda(db1,nx,9)
      call wrtda(db1,nx,16)
      go to 400
350   call reada(da1,nx,8)
      call reada(db1,nx,9)
      if(nc.gt.1) then
         call reada(da2,nx,11)
         call reada(db2,nx,12)
         if(nc.eq.2) go to 355
         call reada(da3,nx,14)
         call reada(db3,nx,15)
355      continue
      endif
400   continue
      iter=iter+1
      call secnd(tim)
      timf=tim
      if(iter.eq.1.and.nvec.ne.numscf.and.aname.ne.'hcore')then
         nes=1
         else
         nes=nc
      endif
      call hstar2(f1,da1,db1,f2,da2,db2,f3,da3,db3,ia,in,nes)
      if(nes.gt.1) call wrtda(f2,nx,10)
      if(nes.gt.2) call wrtda(f3,nx,13)
      call secnd(tim)
      timf=tim-timf
      timh=tim
      ehf0=ehf
      ehf=0.d0
c     symmetrize fock matrix and calculate total energy
c     f1 and f are on the same location
c
      call reada(da1,nx,6)
      call reada(h,nx,3)
      ehf=tracep(da1,h,numscf)*0.5d0
      nav1=5
      nav2=16
      do 410 nk=1,nes
      call symh(f1,f2,ia)
      do 402 i=1,nx
 402  f1(i)=f1(i)+oc(nk)*h(i)/2.0d0
      call wrtda(f1,nx,nav1)
      call reada(da1,nx,nav2)
      ehf=ehf+tracep(f1,da1,numscf)
      if(nprint.ne.5) go to 405
      write(iw,9997)type(nk)
      call fout(f1,numscf)
 405  if(nk.eq.nes) go to 410
      nav1=10+(nk-1)*3
      nav2=nav2+1
      call reada(f1,nx,nav1)
 410  continue
      if(nes.eq.1) go to 411
      call roper(f,da1,da2,da3,f2,f3,db1,db2,ia)
      call cale (f,db2,v,t,1,ending,uhf,ia,in)
      go to 412
 411  call cale (f,da1,v,t,1,ending,uhf,ia,in)
      if(bname.eq.'osrhf'.and.iter.eq.1) ending=.false.
 412  call secnd(tim)
      timh=tim-timh
c
c
      if(maxit.eq.0) go to 1400
  510 continue
c
      ipass=1
      nav1=5
  515 continue
      call reada(d,nx,4)
      if(aname.eq.'hcore') go to 525
      if(nprint.ne.5) go to 525
      if(iter.gt.1) go to 525
      write(iw,9984)
      call fout(d,numscf)
  525 continue
      if(iter.eq.1.or.idiis.eq.(ijump+1)) call reada(f,nx,nav1)
      if(nprint.ne.5) go to 550
      write(iw,9983)
      call fout(f,numscf)
  550 continue
c
c     fock transformation     .... q*f*q
c
      call focktr (h,d,f,t,ia,numscf)
c
c     diagonalize fock matrix
c     coefficient back transformation = q*c
c
c      call ligen(h,v,e,ia,in,numscf)
      call diagonaliser(numscf,h,e,v,numscf)       
      do 560 i=1,numscf
  560 ea(i)=e(i)
c
      if(nprint.ne.4)write(iw,9995) (i,e(i),i=1,numscf)
      call reada(f,nx,4)
      call backtr(v,f,t,ia,in,numscf)
c
 1200 continue
      if(nprint.ne.5) go to 1250
      if(ending) go to 1250
      write(iw,9994)
      call vecout(v,e,in,numscf,numscf)
 1250 continue
c
c
c     form density matrix
c
      call dmop(v,h,da1,db1,da2,db2,da3,db3,in,ia)
      if(nprint.ne.5) go to 1260
c
      write(iw,9993)
      call fout(da2,numscf)
 1260 continue
      ni=numscf*numscf
      do 1270 i=1,numscf
 1270 v(ni+i)=e(i)
      nxst=nxsq+numscf
      call wrtda(v,nxst,7)
 1400 continue
      if(iter.eq.0) ehf0=0.d0
      if(iter.eq.0.and.aname.eq.'hcore') go to 1402
      diff=abs(ehf-ehf0)
c      if(diff.lt.acurcy) ending=.true.
      etot=ehf+en
      call secnd(tim)
      tim1=tim
      delt=tim1-tim0
      tim0=tim1
      if(ending) go to 1500
      write(iw,9998) iter,etot,ehf,timf,timh,delt
      if(maxit.eq.0) go to 1450
c
      if(ending) write(iw,9991)
c
c
 1402 continue
 1420 continue
 1450 continue
c     exit in case of time limit
c     punch out restart data
c
      if((tim-tims).lt.timlim) go to 1475
      if(ending) go to 1460
      write(iw,9987)
      go to 1470
 1460 write(iw,9985)
 1470 continue
      irest=8
      write(iw,9986) irest,nrec
c
      return
 1475 if(iter.lt.maxit) go to 350
      write(iw,9990)
c
c
 1500 continue
      if(ending) write(iw,9968) etot
      if(ending) write(iw,9967) etot
      if(ending) write(iw,9967) etot
      write(iw,9994)
      call reada(v,nxsq,7)
      call vecout(v,ea,in,numscf,numscf)
      if(npunch.ge.0) go to 1502
      rewind 11
      write(11) (v(i),i=1,nxsq),(ea(i),i=1,numscf)
 1502 continue
 1401 continue
      call secnd(tim)
      delt=tim-tims
      write(iw,9972) delt,tim-to
      return
 3500 write(iw,9974)
      stop
 9975 format(/' file (is) contains the results of iteration',i3)
 9974 format(/' namelist vecinp is missing or incorrect')
 9972 format(//,' ... end of scf ...',//,' elapsed time =',f8.3,5x,
     1 ' total time =',f8.3)
 9968 format(/10x, 14hfinal energy =,f13.6)
 9967 format(1h+,23x,f13.6)
      end
