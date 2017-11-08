      subroutine mole(istp)
cdudu
C   Introduction du Born-Mayer 24/07/02 (modifs encadrees par cdudu)
C   Les parametres du Born-Mayer (abm,rhobm) sont lus dans la namelist
C   &HONDO
C   L'energie Born-Mayer est ajputee a l'energie nucleaire et stokee
C   par stkener (typener='energie nucleaire')
cdudu
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      character*40 typegp,typener
      logical scfp,egp,geoto,gestop
      logical pseud,wijkl,ymono
      common /infpot/apot(400),cpot(400),npot(400),nbtyp(200),
     1 ipseud(103),nsom(50),pseud
      common/geot/geoto,gestop
      character*8 aname,bname,group,direct,title
      common/times/ti,tx,tim
      common/elchp/elec(3)
      common/frame/u1,u2,u3,v1,v2,v3,w1,w2,w3,x0,y0,z0,x1,y1,z1,
     1 x2,y2,z2,direct,group,title(10),naxis
      common/output/iop1,itol,icut,normf,normp
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc),en
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      common/psgrx/negp,egp,typegp(50)
      common/scfit/ aname,bname,maxit,nconv,npunch,npi,cextrp,amix,ymono
      common/pseudo/qa(doa),qb(doa),ngla(doa),nglb(doa),numgla,numglb,
     1 scfp
      common/open/nc,noc(3),oc(3),alpha(3,3),beta(3,3),rland
      dimension qn(doa),ngel(doa)
      equivalence (ngla(1),ngel(1))
      namelist/hondo/timlim,nprint,itol,icut,irest,ist,jst,lst,kst,nrec,
     * intloc,group,naxis,mcharg,iunt,x0,y0,z0,x1,y1,z1,x2,y2,z2,direct,
Cdudu1istop,sz,egp,geoto,gestop,elec
     1istop,sz,egp,geoto,gestop,elec,abm,rhobm,ixyz
      namelist/scfinp/maxit,nconv,npunch,cextrp,amix,qa,qb,bname,aname,
     1 ngla,nglb,scfp,qn,ngel,istat,nc,oc,noc,alpha,beta,rland,ymono
c
 7777 format( 20x,10a8)
 7776 format(f10.0,10i3,i10,i5)
 7775 format(/,' time limit = ',d12.5,/,' output option = ',i5,/,
     1 ' prefactor tolerance for integrals = 1.0e-',i2,/,
     1 ' integral cutoff = 1.0e-',i2,/,' option for normalization of',
     1 ' basis functions = ',i5,/,' option for unnormalization of',
     1 ' primitives = ',i5,/,' restart option = ',i5,/,
     1 ' starting shells for 2e-integrals = ',4i5,/,
     1 ' record number = ',i10,/,' location of first integral in this'
     1 , ' record = ',i5)
 7774 format(10a8)
 8859 format(' champ electrique uniforme sur x : ',f12.10)
 8860 format(' champ electrique uniforme sur y : ',f12.10)
 8861 format(' champ electrique uniforme sur z : ',f12.10)
 8858 format(/,' occupation numbers of molecular orbitals')
 8857 format(3i5,a5,a6)
 8856 format(8f10.5)
 8855 format(/)
 8854 format(/,' occupation numbers for alpha orbitals')
 8853 format(/,' occupation numbers for beta  orbitals')
 8852 format(/,(1x,10(i4,f8.5)))
 8851 format(//' ****scf data ****',
     1       /,' maximum number of iterations =',i4,
     2       /,' convergenge threshold        =1.0e-',i2,
     3       /,' starting option              =',a6,
     4       /,' punch out option             =',i4,
     5       /,' scf option                     =',a6,
     6       /,' extrapolation threshold      =',f6.3,
     7       /,' mixing coefficient            =',f6.3)
c
 8850 format(/,' frozen alpha orbitals',/,(30i4))
 8849 format(/,' frozen beta  orbitals',/,(30i4))
 8848 format(/,' number of occupied alpha orbitals   =',i4,
     *       /,' number of occupied beta  orbitals   =',i4)
 8847 format(/,' number of occupied  orbitals   =',i4)
 8999 format(/,' hermiticity parameter  rland =',f8.5)
 9000 format(/,' state parameters. nc =',i2)
 9001 format(/,' occupation numbers of each shell',3f11.7)
 9002 format(3f11.7)
 9003 format(/,' alpha')
 9004 format(/,' beta')
 9005 format(/,' number of occupied orbitals of each shell',3i5)
 9999 format(/,' ----- nuclear energy ----- = ',f20.12)
cdudu
 9998 format(/,' ----- Born Mayer energy ----- = ',f20.12)
 9997 format(/,' ----- Total core energy ----- = ',f20.12)
 9990 format(/,' Parametres Born-Mayer: Abm= ',f10.6,' RHObm=',f10.6)
cdudu
      ixyz=0
      ir=5
      iw=6
      ip=7
      is=8
      iq=9
      ih=10
      iv=11
      if2=2
      if3=3
      if4=4
      rewind is
      ti=0.0d+00
      tx=0.0d+00
      call timit(1)
c
c     read title and options
c
      read(ir,7774) (title(i),i=1,10)
      wijkl=.false.
      geoto=.false.
      gestop=.false.
      timlim=1800.
      nprint=4
      pseud=.false.
      group='c1'
      istop=0
      do 5 i=1,doa
      ngla(i)=0
    5 nglb(i)=0
      numgla=0
      numglb=0
      naxis=0
      x0=0.d0
      y0=0.d0
      z0=0.d0
      x1=0.d0
      y1=0.d0
      z1=0.d0
      x2=0.d0
      y2=0.d0
      z2=0.d0
      mcharg=0
      iunt=0
      sz=0.d0
      direct='parall'
      itol=9
      icut=9
      irest=0
      nrec=0
      egp=.false.
      elec(1)=0.d0
      elec(2)=0.d0
      elec(3)=0.d0
cdudu
      abm=0.d0
      rhobm=0.d0
cdudu
      read(ir,hondo)
      istp = istop
      iop1=nprint
      ich=mcharg
      mul=idint(2.d0*sz)
      if(iop1.eq.4) go to 1
      write(iw,7777) (title(i),i=1,10)
    1 continue
      if(itol.le.0) itol=9
      if(icut.le.0) icut=9
      if(irest.gt.1) go to 3
      ist=1
      jst=1
      kst=1
      lst=1
      nrec=1
      intloc=1
    3 continue
      write(iw,7775) timlim,iop1,itol,icut,normf,normp,irest,
     1 ist,jst,kst,lst,nrec,intloc
c
c     read symmetry point group of the molecule,
c     and generate the transformation matrices.
c
      call ptgrp
c
c     read the basis set for the unique centers,
c     and generate the molecular basis set.
c
      call atoms(iunt,pseud)
      if(ixyz.eq.1) then
      open(unit=15,form='FORMATTED',status='OLD',file='xyz')
      write(6,*)
      write(6,*) ' coordonnees lues en format xyz'
      write(6,*)
        do i=1,nat
         read(15,*) (c(k,i),k=1,3)
         write(6,*) (c(k,i),k=1,3)
         
        enddo
      close(15)
      endif
c
c     nuclear energy
c
      en=0.0d+00
cdudu
      ebm=0.0d+00
cdudu
      if(nat.eq.1) go to 200
      do 150 i=2,nat
      ni=i-1
      do 150 j=1,ni
      rr=0.0d+00
      do 100 k=1,3
  100 rr=rr+(c(k,i)-c(k,j))**2
cdudu
      en=en+zan(i)*zan(j)/sqrt(rr)
        if(abm.ne.0.d0) then
          if(zan(i)*zan(j).lt.0.0) then
          ttt=abm*dexp(-sqrt(rr)/rhobm)
          ebm=ebm+ttt
          endif
        endif
  150 continue
cdudu 150 en=en+zan(i)*zan(j)/sqrt(rr)

c   RV 05/2016 : contribution nucleaire dans le champ electrique        
      do i=1,nat
        rr=elec(1)*c(1,i)+elec(2)*c(2,i)+elec(3)*c(3,i)
        en=en-zan(i)*rr
      end do  
      
  200 write(iw,9999) en
cdudu
        if(abm.ne.0.d0) then
        write(iw,9990)abm,rhobm
        write(iw,9998) ebm
        en=en+ebm
        write(iw,9997) en
        endif
cdudu
c
      typener='energie nucleaire'
      call stkener(typener,en,1)
c     read some scf parameters.
c
      do 2000 i=1,num
      qa(i)=0.d0
      qb(i)=0.d0
 2000 qn(i)=0.d0
      do 2006 i=1,na
      qa(i)=1.d0
 2006 qn(i)=2.d0
      do 2007 i=1,nb
 2007 qb(i)=1.d0
      nc=0
      do 2003 k=1,3
      oc(k)=0.d0
      noc(k)=0
      do 2003 l=1,3
      alpha(k,l)=0.d0
 2003 beta(k,l)=0.d0
      rland=1.0d0
      maxit=30
      nconv=4
      npunch=0
      cextrp=0.1d0
      amix=0.25d0
      aname='hcore'
      bname='rhf'
      ymono=.FALSE.
      read(ir,scfinp)
      if(bname.eq.'RHF')bname='rhf'
      if(bname.eq.'UHF')bname='uhf'
      if(bname.eq.'OSRHF')bname='osrhf'
      if(aname.eq.'HCORE')aname='hcore'
      if(aname.eq.'VECTOR')aname='vector'
      if(bname.eq.'uhf') then
      na=0
      nb=0
      do 2008 i=1,num
      if(qa(i).ne.0.d0) na=na+1
      if(qb(i).ne.0.d0) nb=nb+1
 2008 continue
      write(iw,8848)na,nb
      else
      na=0
      do 2009 i=1,num
      if(qn(i).ne.0.d0) na=na+1
      qa(i)=qn(i)
 2009 continue
      write(iw,8847)na
      end if
c      if(maxit.le.0) maxit=30
      if(nconv.le.0) nconv=6
c
      write(iw,8851)  maxit,nconv,aname,npunch,bname,cextrp,amix
c
c                 print out occupation numbers of the orbitals
c
c
      if(bname.ne.'osrhf') go to 2700
      if(nc.eq.0) then
         call instat(istat)
         else
         write(iw,9000) nc
         write(iw,9001)(oc(k),k=1,nc)
         write(iw,9005)(noc(k),k=1,nc)
         write(iw,9003)
         do 2610 k=1,nc
 2610    write(iw,9002)(alpha(k,l),l=1,nc)
         write(iw,9004)
         do 2620 k=1,nc
 2620    write(iw,9002)(beta(k,l),l=1,nc)
      endif
         write(iw,8999) rland
 2700 continue
      if(bname.ne.'uhf')then
         write(iw,8858)
         write(iw,8852)(i,qa(i),i=1,num)
         else
         write(iw,8854)
         write(iw,8852)(i,qa(i),i=1,num)
         write(iw,8853)
         write(iw,8852) (i,qb(i),i=1,num)
      endif
      do 2760 i=1,doa
      if(ngla(i).ne.0) numgla=numgla+1
      if(nglb(i).ne.0) numglb=numglb+1
 2760 continue
      if(numgla.ne.0) write(iw,8850) (ngla(i),i=1,numgla)
      if(numglb.ne.0) write(iw,8851) (nglb(i),i=1,numglb)
      write(iw,8859) elec(1)
      write(iw,8860) elec(2)
      write(iw,8861) elec(3)
      call dafile
      call timit(1)
      if(istop.ne.1) return
      stop
      end
