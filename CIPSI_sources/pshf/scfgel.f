      subroutine scfgel
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical backt,ending,scfp,uhf
      character*8 aname,bname,type
c
c     -----              scf calculations -----
c               bname='rhf'   restricted scf calculation
c               bname='uhf' unrestricted scf calculation
c               bname='uhf' unrestricted scf calculation
c     on tape (is)       two-electron integrals
c     on file (ih)   -1- general informations(basis set,geometry,..)
c                    -2- overlap
c                    -3- one electron integral matrix
c                    -4- fock transformation matrix
c                    -5- alpha fock matrix
c                    -6- alpha density matrix
c                    -7- alpha molecular orbitals
c                    -8  beta fock matrix
c                    -9- beta density matrix
c                    -10- beta molecular orbitals
c
c     if aname .eq. hcore , a set of starting orbitals will
c     be generated from the bare nucleus hamiltonian matrix
c     if aname .eq. vector , a set of trial vector are readed
c     from cards
c
c     if scfp=.true. perturbative scf calculation
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
      common/times/ti,tx,tim,to
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc),en
c
      common/scfit/ aname,bname,maxit,nconv,npunch,npi,cextrp,amix
      common/conv/acurcy,ehf,ehf0,iter,idiis,ijump
      common/output/nprint,itol,icut,normf,normp,minp
      common/pseudo/qa(doa),qb(doa),ngla(doa),nglb(doa),numgla,numglb,
     1 scfp
      dimension fa(doas),da(doas),fb(doas),db(doas),f(doas),h(doas)
      dimension v(doas1),g(doa),t(doa),e(doa),ia(doa),in(doa),q(doa)
      dimension cc(doas1),cnew(doas1),ign(doa),ig(doa)
      dimension iagel(doa),iga(doa),egel(doa)
      dimension ea(doa),eb(doa),d(doas)
      equivalence (fa(1),f(1)),(da(1),d(1)),(db(1),h(1)),(v(1),fb(1))
      dimension type(2)
      data type/' alpha  ','  beta  '/
      namelist/vecinp/v,ivec,nvec
c
 9999 format(/,' ----- nuclear energy ----- = ',f20.12)
 9998 format(/' cycle ',i3,' energy ',f16.9,' elec. energy',f16.9,
     1 ' time ',3(f6.3,3x))
 9997 format(10x,19(1h-),/,10x,'fock matrix',a8,/,10x,19(1h-))
 9984 format(/,10x,26(1h-),/,10x,'fock transformation matrix',/,
     1 10x,26(1h-))
 9996 format(/,10x,19(1h-),/,10x,'eigenvalues',a8,/,10x,19(1h-))
 9995 format(8(i4,f12.6))
 9989 format(i2,i3,5f15.8)
 9994 format(/,10x,20(1h-),/,10x,'eigenvectors',a8,/,10x,20(1h-))
 9993 format(/,10x,22(1h-),/,10x,'density matrix',a8,/,10x,22(1h-))
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
 9981 format(//,10x,17(1h*),//,'perturbative scf calculation ',//,
     1 10x,17(1h*))
 9982 format(//,10x,17(1h*),//,10x,' uhf calculation ',//,
     1 10x,17(1h*))
      uhf=.false.
      nbgel=numgla
      if(bname.eq.'uhf') uhf=.true.
      if(scfp) write(iw,9981)
      if(.not.scfp.and..not.uhf) write(iw,9980)
      if(.not.scfp.and.     uhf) write(iw,9982)
      rewind is
      do 50 i=1,numscf
      in(i)=(i-1)*numscf
      ia(i)=i*(i-1)/2
50    continue
      do 51 i=1,numscf-numgla
          iga(i)= i*(i-1)/2
          ign(i)=(i-1)*(numscf-numgla)
51    continue
      do 52 i=1,numgla
          iagel(i) =i*(i-1)
          ig(i)=(i-1)*numgla
52    continue
      ehf=99999.0d+00
      ehf0=99999.0d+00
      acurcy=10.0d+00**(-nconv)
c
      iter=0
      idiis=0
      ijump=0
      ending=.false.
      nxsq=numscf*numscf
      nvec=0
c     call second(tim)
      tim0=tim
      tims=tim
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
c     nvec=numscf starting occupied+virtual orbitals for perturb. scf
      read(ir,vecinp,err=3500,end=3500)
      if(ivec.eq.0) go to 208
      nasq=nvec*numscf
      read(ivec) (v(i),i=1,nxsq)
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
      if(nvec.ne.0.and.ipass.eq.1)
     1 call orthon (f,v,e,g,ia,in,ngla,q,numgla,nav2,nvec,scfp)
      if(nvec.ne.0.and.ipass.eq.2)
     1 call orthon (f,v,e,g,ia,in,nglb,q,numglb,nav2,nvec,scfp)
c     store trial density matrix and trial vectors
c
      if(nprint.ne.5) go to 250
      write(iw,9993) type(ipass)
      call fout (f,numscf)
  250 continue
c
      call wrtda(f,nx,nav1)
      call wrtda(v,nxsq,nav2)
      if(ipass.eq.2.or..not.uhf) go to 300
      ipass=2
      nvec=nb
      nav1=9
      nav2=10
      do 203 i=1,numscf
  203 q(i)=qb(i)
      go to 202
c
  210 if(aname.ne.'hcore   ') go to 300
c
c     set h-core matrix to generate initial guess
c
      call reada(f,nx,3)
      call wrtda(f,nx,5)
      if(uhf) call wrtda(f,nx,8)
      go to 511
  300 continue
      call reada(v,nxsq,7)
      call dmtx(v,da,qa,numscf)
      call wrtda(da,nx,6)
        if(uhf) go to 320
      do 310 i=1,nx
  310 db(i)=0.d0
      go to 400
  320 call reada(db,nx,9)
  400 continue
      iter=iter+1
c     call second(tim)
      timf=tim
      call hstar(da,fa,db,fb,ia,uhf)
      do 785 i=1,nx
      if(abs(fa(i)).lt.1.d-8) fa(i)=0.d0
      if(abs(fb(i)).lt.1.d-8) fb(i)=0.d0
  785 continue
      if(uhf) call wrtda(fb,nx,8)
c     call second(tim)
      timf=tim-timf
      timh=tim
      ehf0=ehf
      ehf=0.d0
c     symmetrize fock matrix and calculate total energy
c     da and d,fa and f are on the same location
c
      call reada(h,nx,3)
      ehf=tracep(d,h,numscf)
      call symh(f,fb,ia)
      call addmat(f,h,f,nx)
      call wrtda(f,nx,5)
      ehf=ehf+tracep(f,d,numscf)
c       if(.not.scfp)call cale (f,d,v,t,1,ending,uhf,ia,in)
c     for beta system d,f,and h have been lost
      if(.not.uhf) go to 420
      call reada(f,nx,8)
      call symh(f,fb,ia)
      call reada(d,nx,9)
      call reada(h,nx,3)
      ehf=ehf+tracep(d,h,numscf)
      call addmat(f,h,f,nx)
      call wrtda(f,nx,8)
      ehf=ehf+tracep(f,d,numscf)
      if(.not.scfp)call cale (f,d,v,t,2,ending,uhf,ia,in)
  420 continue
      ehf=ehf*0.5d0
      if(iter.eq.1) oldehf=ehf
      if(iter.gt.1) then
          diff=oldehf-ehf
          if(diff.lt.acurcy.and.diff.ge.0.) ending=.true.
          oldehf=ehf
      end if
c     call second(tim)
      timh=tim-timh
c
c
      if(maxit.eq.0) go to 1400
c
  511 continue
      ipass=1
      nav1=5
      do 512 i=1,numscf
  512 q(i)=qa(i)
  515 continue
c     go to perturbative scf calculations if required
      if(.not.scfp) go to 510
c     if there are occupied only starting orbitals a first standard
c     iteration is necessary to generate virtual orbital
      if(aname.eq.'vector'.and.iter.eq.1.and.nvec.eq.na) go to 510
      call pert(f,v,d,e,ia,in,ipass,ending)
      go to 1200
  510 continue
      call reada(d,nx,4)
      if(aname.eq.'hcore') go to 525
      if(nprint.ne.5) go to 525
      if(iter.gt.1) go to 525
      if(ipass.eq.2) go to 525
      write(iw,9984)
      call fout(d,numscf)
  525 continue
      if(uhf.or.iter.eq.1.or.idiis.eq.(ijump+1)) call reada(f,nx,nav1)
      if(nprint.ne.5) go to 550
      write(iw,9997) type(ipass)
      call fout(f,numscf)
  550 continue
c
c
c      call focktr (h,d,f,t,ia,numscf)
       if(iter.eq.1) call reada(cc,nxsq,7)
c       write(6,*) 'om en base d oa avant diag'
c       call vecout(cc,e,in,numscf,numscf)
c creation de la matrice v contenant les om gelees
       call getgel(cc,v,numscf,numgla,ngla,in)
c       write(6,*) 'om gelees:'
c       call vecout(v,e,in,numscf,numgla)
c transformation de la matrice de fock correspondant aux om gelees
c dans la base d'om
       call ftrgel(v,f,h,numscf,numgla,in,ia,iagel,ngla)
c diagonalisation de cette matrice pour obtenir els energies orbitales
c des om gelees egel
       call hdiag(h,v,egel,iagel,ig,numgla)
c        write(6,*) 'egel :'
c        write(6,*) (egel(i),i=1,numgla)
c       write(6,*) 'apres diag de fgel:'
c       call vecout(v,egel,ig,numgla,numgla)
c transformation de la matrice de fock en base d'om
       call ftrom(cc,f,h,numscf,in,ia)
c       write(6,*) 'matrice de fock en base d om'
c       call fout(h,numscf)
c creation de la matrtice de fock correspondant aux om non gelees
       call elimat(h,numgla,ngla,numscf,ia,iga)
c       write(6,*) 'matrice de fock en base d om apres elimination'
c       call fout(h,numscf-numgla)
c creation de la matrice contenant les om non gelees dans la base d'oa
       call elimv(cc,v,numgla,ngla,in,numscf)
c       write(6,*) 'om non gelees en base d oa'
c       call vecout(v,e,in,numscf,numscf-numgla)
c
c     diagonalize fock matrix
c
c diagonalisation de la matrica de fock correspondant aux om non gelees
      call hdiag(h,cnew,t,iga,ign,numscf-numgla)
c      write(6,*) 'diagonalisation de la matrice fnongel'
c      call vecout(cnew,t,ign,numscf-numgla,numscf-numgla)
      k=1
      ii=0
      do 559 i=1,numscf
          if(i.ne.ngla(k)) then
            ii=ii+1
            e(i)=t(ii)
          else
            e(i)=egel(k)
            k=k+1
          end if
559   continue
      if(uhf.and.ipass.eq.2) go to 570
      do 560 i=1,numscf
  560 ea(i)=e(i)
      go to 580
  570 do 565 i=1,numscf
  565 eb(i)=e(i)
  580 continue
c
      if(nprint.ne.4.and..not.uhf)write(iw,9995) (i,e(i),i=1,numscf)
      call reada(f,nx,4)
c      call backtr(v,f,t,ia,in,numscf)
c creation des nouvelles om dans la base d'oa
      call backoa(v,cnew,numscf,numscf-numgla,in,ign)
c      write(6,*) 'nouvelles om(non gelee) dans la base d oa'
c      call vecout(cnew,e,in,numscf,numscf-numgla)
c creation de la matrice d'om totale dans la base d'oa
      call insert(cc,cnew,numscf,numgla,in,ngla)
c      write(6,*) ' om ds la base d oa apres diag'
c      call vecout(cc,e,in,numscf,numscf)
      do 1 j=1,numscf
          do 2 i=1,numscf
              ij=in(j)+i
              v(ij)=cc(ij)
2         continue
1     continue
c
 1200 continue
      if(nprint.ne.5) go to 1250
      if(ending) go to 1250
      write(iw,9994) type(ipass)
      call vecout(v,e,in,numscf,numscf)
 1250 continue
c
c
c     form density matrix
c
      call dmtx(v,f,q,numscf)
      if(nprint.ne.5) go to 1260
c
      write(iw,9993)type(ipass)
      call fout(f,numscf)
 1260 continue
      call wrtda(f,nx,nav1+1)
      ni=numscf*numscf
      do 1270 i=1,numscf
 1270 v(ni+i)=e(i)
      nxst=nxsq+numscf
      call wrtda(v,nxst,nav1+2)
      if(ipass.eq.2.or..not.uhf) go to 1400
      ipass=2
      nav1=8
      do 1300 i=1,numscf
 1300 q(i)=qb(i)
      go to 515
 1400 continue
      if(iter.eq.0) ehf0=0.d0
      if(iter.eq.0.and.aname.eq.'hcore') go to 1402
      diff=abs(ehf-ehf0)
      etot=ehf+en
c     call second(tim)
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
 1475 if(iter.lt.maxit) go to 300
      write(iw,9990)
c
c
 1500 continue
      if(ending) write(iw,9968) etot
      if(ending) write(iw,9967) etot
      if(ending) write(iw,9967) etot
c
c
      write(iw,9994)type(1)
      call reada(v,nxst,7)
      call vecout(v,v(nxsq+1),in,numscf,numscf)
c      if(.not.scfp) call gelom(v,e,fa,fb,da,db,ia,in,uhf)
      if(npunch.ge.0) go to 1502
      rewind 11
      write(11) (v(i),i=1,nxsq),(ea(i),i=1,numscf)
 1502 continue
      if(.not.uhf) go to 1401
      write(iw,9994) type(2)
      call reada(v,nxst,10)
      call vecout (v,v(nxsq+1),in,numscf,numscf)
      if(npunch.ge.0) go to 1503
      write(11) (v(i),i=1,nxsq),(eb(i),i=1,numscf)
 1503 continue
 1401 continue
c     if(.not.uhf) call nesbet
c
c     if(npi.ne.0) call potion
c     call second(tim)
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
