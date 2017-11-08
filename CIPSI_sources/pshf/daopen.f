      subroutine daopen
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      character*8 aname,bname
      character*6 ante
      character*4 iflab
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      common/scfit/aname,bname
      common/restar/timlim,irest
      common/runlab/iflab(3,doa)
      dimension a(doasq)
      read(if2) ioda,nt,nshell,num
      if(irest.le.3) return
      nx=num*(num+1)/2
      nxsq=num*num
      do 10 k=2,19
      nn=nx
      if(k.eq.7) nn=nxsq
      if(bname.eq.'uhf'.and.k.eq.10) nn=nxsq
      read(if2) (a(i),i=1,nn)
      call wrtda(a,nn,k)
   10 continue
      write(iw,1000) ih,if2
 1000 format(//,' direct access file ',i3,' restored from sequential
     1 file',i3)
      return
      end
