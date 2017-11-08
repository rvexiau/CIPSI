      subroutine ftrgel(tr,f,h,noa,nvec,in,ia,iagel,ind)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
       dimension tr(*),f(*),h(*),in(*),ia(*),ind(nvec),iagel(*)
       dimension t(doa)
       do 10 j=1,nvec
          do 20 k=1,noa
             dum=0.d0
             do 30 l=1,noa
                lj=in(j)+l
                kl=ia(l)+k
                if(k.gt.l) kl=ia(k)+l
                dum=dum+f(kl)*tr(lj)
30           continue
             t(k)=dum
20        continue
          do 40 i=1,nvec
             dum=0.d0
             do 50 k=1,noa
                ki=in(i)+k
                dum=dum+tr(ki)*t(k)
50           continue
             ij=iagel(j)+i
             if(i.gt.j) ij=iagel(i)+j
             h(ij)=dum
40         continue
10    continue
      return
      end
