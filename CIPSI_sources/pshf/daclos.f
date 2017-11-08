      subroutine daclos
c     RV 05/2016 : don't overwrite energy 'en' when reading the older record        
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      character*8 group,aname,bname,title,direct
      character*80 apseud
      character*6 ante
      character*4 iflab
      logical ymono
      common/frame/u1,u2,u3,v1,v2,v3,w1,w2,w3,x0,y0,z0,x1,y1,z1,
     1 x2,y2,z2,direct,group,title(10),naxis
      common/fil2/ a(dc),ns(dc),ks(dc),newsh(ds,48),ngauss
      common/output/iop1,itol,icut,normf,normp
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
      common/symtry/invt(48),iso(ds,12),ptr(3,144),dtr(6,288),
     1 ftr(10,480),nt
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc),en
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      common/scfit/ aname,bname,maxit,nconv,npunch,npi,cextrp,amix,ymono
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell,kmaz(ds)
      common/psnloc/apseud(dc),iatno(dc)
      common/pseudo/qa(doa),qb(doa),ngla(doa),nglb(doa),numgla,numglb
      common/runlab/iflab(3,doa)
      common/elchp/elec(3)
      dimension axx(doasq)
c
c
      rewind if2
      read(ih,rec=1) ioda,nt,nshell,num,norba,norbb,bname,ngauss,ne,
     1 ich,mul,na,nb,nat,group,naxis,title,
     2(ktype(k),kstart(k),kloc(k),katom(k),kng(k),k=1,nshell),
     3 (ex(k),cs(k),cp(k),cd(k),cf(k),cg(k),k=1,ngauss),
     4((c(i,k),i=1,3),a(k),ns(k),ks(k),iatno(k),zan(k),k=1,nat),
     5((newsh(i,k),i=1,nshell),k=1,nt),((ptr(i,k),i=1,3),k=1,3*nt),
     6((dtr(i,k),i=1,6),k=1,6*nt),((ftr(i,k),i=1,10),k=1,10*nt)
c    7                ,((gtr(i,k),i=1,15),k=1,15*nt)
     8  ,((iflab(i,k),i=1,3),k=1,num),((iso(i,j),i=1,nshell),j=1,12)
     9  ,invt,entmp,ante,ymono
      write(if2) ioda,nt,nshell,num,norba,norbb,bname,ngauss,ne,
     1 ich,mul,na,nb,nat,group,naxis,title,
     2(ktype(k),kstart(k),kloc(k),katom(k),kng(k),k=1,nshell),
     3 (ex(k),cs(k),cp(k),cd(k),cf(k),cg(k),k=1,ngauss),
     4((c(i,k),i=1,3),a(k),ns(k),ks(k),iatno(k),zan(k),k=1,nat),
     5((newsh(i,k),i=1,nshell),k=1,nt),((ptr(i,k),i=1,3),k=1,3*nt),
     6((dtr(i,k),i=1,6),k=1,6*nt),((ftr(i,k),i=1,10),k=1,10*nt)
c    7                ,((gtr(i,k),i=1,15),k=1,15*nt)
     8  ,((iflab(i,k),i=1,3),k=1,num),((iso(i,j),i=1,nshell),j=1,12)
     9  ,invt,en,ante,ymono,(elec(k),k=1,3)
      nx=num*(num+1)/2
      nxsq=num*num
      do 10 k=2,19
      nn=nx
      if(k.eq.7) nn=nxsq
      if(bname.eq.'uhf'.and.k.eq.10) nn=nxsq
      call reada(axx,nn,k)
      write(if2) (axx(i),i=1,nn)
   10 continue
      write(iw,1000)ih,if2
 1000 format(//,' direct access file ',i3,' saved on sequential file',
     1 i3)
      return
      end
