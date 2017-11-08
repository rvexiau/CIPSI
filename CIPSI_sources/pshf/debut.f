      subroutine debut
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 out
      common/output/nprint,itol,icut,normf,normp
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
      common/buf/tei(lsize),ila(lsize)
      common/shlt/tol,cutoff,icount,out
      common/iofile/ir,iw,ip,is
      write(iw,1000)
 1000 format(//,10x,20(1h*),//,10x,'2 electron integrals',
     1 //,10x,20(1h*))
      out=nprint.eq.1
      if(icut.le.0) icut=9
      cutoff=1.0d+00/(1.0d+01**icut)
c     modification jpd novembre 90
c     if(icut.eq.20) on stocke toutes les integrales
      if(icut.eq.20) cutoff=-1.d0
c
      if(itol.le.0) itol=20
      tol=2.30258d+00*itol
      icount=1
      if(irest.le.1) go to 40
      icount=intloc
      n=nrec-1
      do 20 i=1,n
   20 read(is)
      read(is) tei,ila
      rewind is
      do 30 i=1,n
   30 read(is)
   40 continue
      return
      end
