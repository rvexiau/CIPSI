      subroutine wrtda (a,n,idaf)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      dimension a(n)
      dimension b(lsize)
      nav=ioda(idaf)
      nrec=(n+lsize-1)/lsize
      js=0
      imax=lsize
      do 10 irec=1,nrec
      if(js+lsize.gt.n) imax=n-js
      do 20 i=1,imax
   20 b(i)=a(i+js)
      write(ih,rec=nav) b
      js=js+lsize
      nav=nav+1
   10 continue
      return
      end
