      subroutine set
Ccc   RG : modified on 25/11/08 : Make sure the cutoff radii are positive for the calculation      
Ccc   RG : modified on 25/11/08 : Added a rewind 5
c
c
c     ------------------------------------------------------------------
c
c     lecture des donnees du present programme :
c
c
c     - plmax
c     - (no(i),i=1,2)  i=1 : n (pour 1/r**4), i=2 : n (pour x,y,z/r**3).
c     - (rcut(i),i=0,plmax) : rayon de coupure pour chaque pl.
c                 l'unite 2.
c
c
c     - lecture des donnees generees par [hondom].
c
c
c     n.b. : les expressions des derivees generees par le programme :
c     ****   [derivee.f] sont lues sur l'unite 10 (derpl22) et les re-
c            sultats numeriques des integrales 1/r**4 et x,y,z/r**3 sont
c            ecrits sur la file : 1sr4xsr3 (unite 17).
c
c
c     ------------------------------------------------------------------
c
c
c
 
      implicit  double precision  (a-h,o-z)
      integer ent,sor,plmax,plder
      include 'pshf.prm'
      parameter (itypmax=3,plder=3)
      character*8 bname,group,title,a
      dimension ioda(19),title(10)
      dimension mmin(10),mmax(10),mmaz(10)
      dimension a(dc),ns(dc),ks(dc)
      dimension zan(dc)

      logical ysoft                              !OD, 13/2/2004

      include 'NSHEL.cmm'
      common/infoa/nat,nx,c(3,dc)
      common/infpot/iatno(dc),ipseud(dc)
      common/par/ijmax,npartab,plmax
      common/rayona/rcut1(dc,0:plder)
      common/puisno/no(2)
      common/lecrit/ent,sor
      
      logical diffCalc,firstRun
      common /diff/ diffCalc,firstRun
      
      namelist/rcutval/plmax,rcut1,ysoft,no           !OD, 13/2/2004

      data mmin/1,2,5,11,21,1,1,5,11,21/
      data mmax/1,4,10,20,35,1,4,10,20,35/
      data mmaz/1,4,9,17,30,1,4,10,20,35/
c
c
 5000 format(1x,'donnees lues sur la file :',a20,'generee par hondom')
c
c
      ent=5
      sor=6
c
c
c      read(ent,*) plmax			   !OD, 13/2/2004
c
c     plder : pl maximal dans le programme [derivee.f]. ce dernier
c             programme genere actuellement les derivees des pl=0-3
c             (inclus).
c

      ysoft =.FALSE.
      read(ent,rcutval)
      
      if(plmax.gt.plder) then
                write(sor,*) 'plmax > plder, incompatibilite'
                stop
      endif
c
c
c
      npartab=15+itypmax*(2*itypmax+5)
      ijmax=(itypmax+1)*(itypmax+2)/2
c
c
c     normf=0
c     normp=0
c
c
c      read(ent,*) (no(i),i=1,2)           !OD, 13/2/2004
cvax
c      open(10,form='unformatted',
c     1             status='old',shared)
c      open(17,form='unformatted',status='unknown')
c      open(2,form='unformatted',status='old')
cvax
      rewind 10
      rewind 17
      rewind 2
      rewind 5
c
c
c
 9000 read(2) ioda,nt,nshell,num,nveca,nvecb,bname,ngauss,
     *ne,ich,mul,na,nb,nat,group,naxis,title,
     *(ktype(k),kstart(k),kloc(k),katom(k),kng(k),k=1,nshell),
     *(ex(k),cs(k),cp(k),cd(k),cf(k),cg(k),k=1,ngauss),
     *((c(i,k),i=1,3),a(k),ns(k),ks(k),iatno(k),zan(k),k=1,nat)
cCc ============================
cCc = Differential calculation =
cCc ============================
      if(diffCalc) then
        rcut1=abs(rcut1)
        if(firstRun) then
          do i=1,nat
            xtemp=rcut1(i,0)
            rcut1(i,0:plmax)=xtemp
          enddo
          write(*,*) 'Diff calc. First run. rcut1() array : '
          write(*,*) rcut1(1:nat,0:plmax)
        else
          write(*,*) 'Diff calc. second run. rcut1() array : '
          write(*,*) rcut1(1:nat,0:plmax)
        endif
      endif
c
      do 50 i=1,nat
      iatn=iatno(i)
      znuc=zan(i)
      write(6,*) 'i,iatn,znuc',i,iatn,znuc
      ipseud(i)=0
c      if((iatn.le.0) .or. (iatn.eq.idint(znuc))) go to 50
      ipseud(i)=1
c      read(ent,*) (rcut1(i,j),j=0,plmax)           !OD, 13/02/2004
      write(sor,*)'rcut',(rcut1(i,j),j=0,plmax)
      if(rcut1(i,0).eq.0.d0) then
       ipseud(i)=0
      else
      write (sor,*) ' i rcut1(i,j)=',i,(rcut1(i,j),j=0,plmax)
      endif
 50   continue
c
c
c
      write(sor,*) 'i=','ipseud(i)=',(i,ipseud(i),i=1,nat)
c
c
c
      do 9005 i=1,nshell
      kmin(i)=mmin(ktype(i))
      kmax(i)=mmax(ktype(i))
      kmaz(i)=mmaz(ktype(i))
 9005 continue
      nx=(num*(num+1))/2
      return
      end
