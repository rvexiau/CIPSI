      subroutine daclos
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      character*8 bname,group,title,a
      character*6 ante
      logical ymono
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      common/fil2/ a(dc),ns(dc),ks(dc),newsh(ds,48),ngauss
      common/runlab/title(10)
      common/symtry/invt(48),iso(ds,12),ptr(3,144),dtr(6,288),
     1 ftr(10,480),nt
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc),
     1 ymono,mrec
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),
     1 cf(dgp),cg(dgp),nshell,kstart(ds),katom(ds),
     2 ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     3 kmaz(ds)
      common/polxx/polz
      dimension iatno(dc),iflab(3,doa)
      dimension axx(doa*doa)
      rewind if3
      read(ih,rec=1) ioda,nt,nshell,num,norba,norbb,bname,ngauss,ne,
     1 ich,mul,na,nb,nat,group,naxis,title,
     2(ktype(k),kstart(k),kloc(k),katom(k),kng(k),k=1,nshell),
     3 (ex(k),cs(k),cp(k),cd(k),cf(k),cg(k),k=1,ngauss),
     4((c(i,k),i=1,3),a(k),ns(k),ks(k),iatno(k),zan(k),k=1,nat),
     5((newsh(i,k),i=1,nshell),k=1,nt),((ptr(i,k),i=1,3),k=1,3*nt),
     6((dtr(i,k),i=1,6),k=1,6*nt),((ftr(i,k),i=1,10),k=1,10*nt)
c    7                ,((gtr(i,k),i=1,15),k=1,15*nt)
     8  ,((iflab(i,k),i=1,3),k=1,num),((iso(i,j),i=1,nshell),j=1,12)
     9  ,invt,en,ante
c
      write(if3) ioda,nt,nshell,num,norba,norbb,bname,ngauss,ne,
     1 ich,mul,na,nb,nat,group,naxis,title,
     2(ktype(k),kstart(k),kloc(k),katom(k),kng(k),k=1,nshell),
     3 (ex(k),cs(k),cp(k),cd(k),cf(k),cg(k),k=1,ngauss),
     4((c(i,k),i=1,3),a(k),ns(k),ks(k),iatno(k),zan(k),k=1,nat),
     5((newsh(i,k),i=1,nshell),k=1,nt),((ptr(i,k),i=1,3),k=1,3*nt),
     6((dtr(i,k),i=1,6),k=1,6*nt),((ftr(i,k),i=1,10),k=1,10*nt)
c    7                ,((gtr(i,k),i=1,15),k=1,15*nt)
     8  ,((iflab(i,k),i=1,3),k=1,num),((iso(i,j),i=1,nshell),j=1,12)
     9  ,invt,en,ante,ymono

      nx=num*(num+1)/2
      nxsq=num*num
      do 10 k=2,19
      nn=nx
      if(k.eq.7) nn=nxsq
      if(bname.eq.'uhf'.and.k.eq.10) nn=nxsq
      call reada(axx,nn,k)
      write(if3) (axx(i),i=1,nn)
   10 continue
      write(iw,1000)ih,if3
 1000 format(//,' direct access file ',i3,' saved on sequential file',
     1 i3)
      return
      end
