      subroutine dafile
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      character*8 aname,direct,title
      character*8 group,bname
      character*80 apseud
      character*6 ante
      character *4 iflab
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
c     common/pseudo/qn(80)
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell,kmaz(ds)
      common/psnloc/apseud(dc),iatno(dc)
      common/pseudo/qa(doa),qb(doa),ngla(doa),nglb(doa),numgla,numglb
      common/runlab/iflab(3,doa)
      dimension axx(doasq)
c
c      determine length of records
c      open (unit=ih,access='direct',recl=32*lsize,status='scratch')
       open (unit=ih,access='direct',recl=32*lsize,status='unknown',
     1    file='F10_SCR')
c
c      determine length of first record
      n1=(23+5*nshell)*4+(12+6*ngauss)*8
      n1=n1+(3*4+5*8)*nat
      n1=n1+nt*nshell*4+(9+36+100)*nt*8
c     n1=n1+225*nt*8
      ioda(1)=1
      n1m=n1/(lsize*8)
      if(n1.gt.lsize*n1m*8) n1m=n1m+1
      ioda(2)=ioda(1)+n1m
c      determine length of records for a half matrix of dimension num
      n1m=nx/lsize
c     determine length of record for a complete matrix+ one line
      nmn= num*(num+1)
      n2m=nmn/lsize
      if(nmn.gt.lsize*n2m) n2m=n2m+1
      if(nx.gt.lsize*n1m) n1m=n1m+1
c     determine location for overlap,one-electron,s-1/2,density and
c      vector matrices
      ioda(3)=ioda(2)+n1m
      ioda(4)=ioda(3)+n1m
      ioda(5)=ioda(4)+n1m
      ioda(6)=ioda(5)+n1m
      ioda(7)=ioda(6)+n1m
      ioda(8)=ioda(7)+n2m
      ioda(9)=ioda(8)+n1m
      ioda(10)=ioda(9)+n1m
      nmm=n1m
      if(bname.eq.'uhf')nmm=n2m
      ioda(11)=ioda(10)+nmm
      ioda(12)=ioda(11)+n1m
      ioda(13)=ioda(12)+n1m
      ioda(14)=ioda(13)+n1m
      ioda(15)=ioda(14)+n1m
      ioda(16)=ioda(15)+n1m
      ioda(17)=ioda(16)+n1m
      ioda(18)=ioda(17)+n1m
      ioda(19)=ioda(18)+n1m
      kp=3*nt
      kd=6*nt
      kf=10*nt
      kg=15*nt
c     write general informations on first record
      nav=ioda(1)
c      norba=num-numgla
c      norbb=num-numglb
      norba=num
      norbb=num
      ante='  pshf'      
      write(ih,rec=nav) ioda,nt,nshell,num,norba,norbb,bname,ngauss,ne,
     1 ich,mul,na,nb,nat,group,naxis,title,
     2(ktype(k),kstart(k),kloc(k),katom(k),kng(k),k=1,nshell),
     3 (ex(k),cs(k),cp(k),cd(k),cf(k),cg(k),k=1,ngauss),
     4((c(i,k),i=1,3),a(k),ns(k),ks(k),iatno(k),zan(k),k=1,nat),
     5((newsh(i,k),i=1,nshell),k=1,nt),((ptr(i,k),i=1,3),k=1,kp),
     6((dtr(i,k),i=1,6),k=1,kd),((ftr(i,k),i=1,10),k=1,kf)
c    7                ,((gtr(i,k),i=1,15),k=1,kg)
     8  ,((iflab(i,k),i=1,3),k=1,num),((iso(i,j),i=1,nshell),j=1,12)
     9  ,invt,en,ante,ymono
      if(irest.ge.3) return
c     write on last record to check if file ih is long enough
      do 50 k=2,19
      nmm=nx
      if(k.eq.7)nmm=nmn
      if(bname.eq.'uhf'.and.k.eq.10)nmm=nmn
      call wrtda(axx,nmm,k)
 50   continue
      nred=ioda(19)+n1m-1
      write(iw,9999) nred ,lsize
 9999 format(//,' number of records needed on file ih',2i6)
      return
      end
