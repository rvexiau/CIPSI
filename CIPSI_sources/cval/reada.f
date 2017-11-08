      subroutine reada (a,n,idaf)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      dimension b(lsize)
      dimension a(n)
      nav=ioda(idaf)
      nrec=(n+lsize-1)/lsize
      js=0
      imax=lsize
      do 10 irec=1,nrec
      read(ih,err=1000,rec=nav) b
      if(js+lsize.gt.n) imax=n-js
      do 20 i=1,imax
   20 a(i+js)=b(i)
      js=js+lsize
      nav=nav+1
   10 continue
      return
 1000 write(iw,9999) ih,idaf,n
      stop
 9999 format(//' error in processing direct access on file',i3,' for rec
     1ord number',i4,' of length',i6)
      end
