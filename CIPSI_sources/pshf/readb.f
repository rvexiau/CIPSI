      subroutine readb (a,n,iter,index)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      common/iofile/ir,iw,ip,is,iq,ih
      common/extrp/bb(20,20),d(20),lrec,nxr,iodb(5000)
      dimension a(n)
      dimension b(lsize)
      nav=iodb(iter)+lrec*(index-1)
      nrec=(n+lsize-1)/lsize
      js=0
      imax=lsize
      do 10 irec=1,nrec
      read (iq,rec=nav) b
      if(js+lsize.gt.n) imax=n-js
      do 20 i=1,imax
      a(js+i)=b(i)
   20 continue
      js=js+lsize
      nav=nav+1
   10 continue
      return
      end
