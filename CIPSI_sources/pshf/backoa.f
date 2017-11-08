      subroutine backoa(tr,v,noa,nvec,in,ign)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      dimension tr(*),v(*),ign(*),in(*)
      dimension temp(doas)
      do 10 j=1,nvec
         do 20 i=1,noa
             dum=0.d0
             do 30 k=1,nvec
                 ik=in(k)+i
                 kj=ign(j)+k
                 dum=dum+tr(ik)*v(kj)
30           continue
             ij=in(j)+i
             temp(ij)=dum
20        continue
10    continue
      do 40 j=1,nvec
          do 50 i=1,noa
              ij=in(j)+i
              v(ij)=temp(ij)
50        continue
40    continue
      return
      end
