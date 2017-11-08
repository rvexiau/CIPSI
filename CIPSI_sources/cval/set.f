      subroutine set
      implicit  double precision  (a-h,o-z)
      include 'pshf.prm'
      character*8 bname,group,title,a
      character*6 ante
c      character*50 nom
      logical ymono
      common/fil2/ a(dc),ns(dc),ks(dc),newsh(ds,48),ngauss
      common/symtry/invt(48),iso(ds,12),ptr(3,144),dtr(6,288),
     1 ftr(10,480),nt
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
c      common /nomfil/nom
      common/runlab/title(10)
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc),
     1 ymono
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),
     1 cf(dgp),cg(dgp),nshell,kstart(ds),katom(ds),
     2 ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     3 kmaz(ds)
      common/elchp/elec(3)
      dimension axx(doa*doa)
      dimension iatno(dc),iflab(3,doa)
      dimension mmin(10),mmax(10),mmaz(10)
      data mmin/1,2,5,11,21,1,1,5,11,21/
      data mmax/1,4,10,20,35,1,4,10,20,35/
      data mmaz/1,4,9,17,30,1,4,10,20,35/
      data pi32/5.5683279968317d0/
      
      open (unit=ih,access='direct',recl=16*lsize,status='scratch')
c
c
c     relecture des donnees (geometrie et base atomique sur 2)
c

      read(if2) ioda,nt,nshell,num,nveca,nvecb,bname,ngauss,
     1ne,ich,mul,na,nb,nat,group,naxis,title,
     2(ktype(k),kstart(k),kloc(k),katom(k),kng(k),k=1,nshell),
     3(ex(k),cs(k),cp(k),cd(k),cf(k),cg(k),k=1,ngauss),
     4((c(i,k),i=1,3),a(k),ns(k),ks(k),iatno(k),zan(k),k=1,nat)
     5,((newsh(i,k),i=1,nshell),k=1,nt)
     5,((ptr(i,k),i=1,3),k=1,3*nt),
     6((dtr(i,k),i=1,6),k=1,6*nt),((ftr(i,k),i=1,10),k=1,10*nt)
c    7                ,((gtr(i,k),i=1,15),k=1,15*nt)
     8  ,((iflab(i,k),i=1,3),k=1,num),((iso(i,j),i=1,nshell),j=1,12)
     9  ,invt,en,ante,ymono,(elec(i),i=1,3)

      do 9005 i=1,nshell
      kmin(i)=mmin(ktype(i))
      kmax(i)=mmax(ktype(i))
      kmaz(i)=mmaz(ktype(i))
 9005 continue
 9010 continue
      write(ih,rec=1) ioda,nt,nshell,num,nveca,nvecb,bname,ngauss,
     1ne,ich,mul,na,nb,nat,group,naxis,title,
     2(ktype(k),kstart(k),kloc(k),katom(k),kng(k),k=1,nshell),
     3(ex(k),cs(k),cp(k),cd(k),cf(k),cg(k),k=1,ngauss),
     4((c(i,k),i=1,3),a(k),ns(k),ks(k),iatno(k),zan(k),k=1,nat)
     5,((newsh(i,k),i=1,nshell),k=1,nt)
     5,((ptr(i,k),i=1,3),k=1,3*nt),
     6((dtr(i,k),i=1,6),k=1,6*nt),((ftr(i,k),i=1,10),k=1,10*nt)
c    7                ,((gtr(i,k),i=1,15),k=1,15*nt)
     8  ,((iflab(i,k),i=1,3),k=1,num),((iso(i,j),i=1,nshell),j=1,12)
     9  ,invt,en,ante,ymono,(elec(i),i=1,3)
      nx=num*(num+1)/2
      nxsq=num*num
      do 10 k=2,19
      nn=nx
      if(k.eq.7) nn=nxsq
      if(bname.eq.'uhf'.and.k.eq.10) nn=nxsq
      read(if2) (axx(i),i=1,nn)
      call wrtda(axx,nn,k)
   10 continue
      write(iw,1000) ih,if2
 1000 format(//,' direct access file ',i3,' restored from sequential
     1 file',i3)
      write(6,*) ' champ electrique externe'
      do k=1,3
      write(6,*) elec(k)
      enddo
      return
      end
