      subroutine matout(f,numscf)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
c
c     print out f-matrix
c
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      common/output/nprint,itol,icut,normf,normp,nopk
      common/runlab/iflab(3,doa)
      dimension f(2),n(doa),ia(doa)
 1000 format(/)
 1001 format(18x,7(2x,3a4,1x))
 1002 format(i5,1x,3a4,7e15.8)
      do 50 i=1,numscf
   50 ia(i)=i*(i-1)/2
      max=7
      if(nprint.eq.6) max=3
      imax=0
  100 imin=imax+1
      imax=imax+max
      if(imax.gt.numscf) imax=numscf
      write(iw,1000)
      write(iw,1001) ((iflab(k,i),k=1,3),i=imin,imax)
      write(iw,1000)
      do 400 j=1,numscf
      do 300 i=imin,imax
      if(i.gt.j) go to 200
      n(i)=ia(j)+i
      go to 300
  200 n(i)=ia(i)+j
  300 continue
      write(iw,1002) j,(iflab(k,j),k=1,3),(f(n(i)),i=imin,imax)
  400 continue
      if(imax.lt.numscf) go to 100
      return
      end
