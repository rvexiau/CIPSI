      subroutine set
c
c
c     ------------------------------------------------------------------
c
c     lecture des donnees du present programme :
c
c
c     - (no(i),i=1,2)  i=1 : n (pour 1/r**4), i=2 : n (pour x,y,z/r**3).
c     - rcut(i) : rayon de coupure pour chaque atome.
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
Ccc   RG : modified on 02/11/07 : Added a rewind statement to the stdin.
 
      implicit  double precision  (a-h,o-z)
      integer ent,sor
      include 'pshf.prm'
      parameter (itypmax=3)
      real*8 rcut(dc)
      logical ysoft
      character*8 bname,group,title,a

      dimension ioda(19),title(10)
      dimension mmin(10),mmax(10),mmaz(10)
      dimension a(dc),ns(dc),ks(dc)
      dimension zan(dc)
      include 'NSHEL.cmm'
      common/infoa/nat,nx,c(3,dc)
      common/infpot/iatno(dc),ipseud(dc)
      common/par/ijmax,npartab
      common/puisno/no(2)
      common/lecrit/ent,sor
      common /rayon/rcut,ysoft

      namelist /rcval/rcut,ysoft


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
      rewind 10
      rewind 17
      rewind 2
      rewind ent  ! RG 2 novembre 2007
c
c
c
 9000 read(2) ioda,nt,nshell,num,nveca,nvecb,bname,ngauss,
     *ne,ich,mul,na,nb,nat,group,naxis,title,
     *(ktype(k),kstart(k),kloc(k),katom(k),kng(k),k=1,nshell),
     *(ex(k),cs(k),cp(k),cd(k),cf(k),cg(k),k=1,ngauss),
     *((c(i,k),i=1,3),a(k),ns(k),ks(k),iatno(k),zan(k),k=1,nat)
     
     
c
      ysoft =.FALSE.
      read(ent,rcval)

      do 50 i=1,nat
      iatn=iatno(i)
      znuc=zan(i)
      write(6,*) 'i,iatn,znuc',i,iatn,znuc
      ipseud(i)=0
c      if((iatn.le.0) .or. (iatn.eq.idint(znuc))) go to 50
      ipseud(i)=1
c      rcut(i) = 1.455
      if(rcut(i).eq.0.d0) then
       ipseud(i)=0
      else
      write (sor,*) ' i rcut(i)=',i,rcut(i)
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
