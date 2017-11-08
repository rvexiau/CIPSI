      subroutine ecrita(f,n,ia)
      implicit real*8(a-h,o-z)
      dimension f(*),ia(*)
      do 10 j=1,n
         write(6,*) ' '
         write(6,*) 'vecteur: ',j
         do 20 i=1,j
            ij = ia(j)+i
            write(6,*) f(ij)
20       continue
10    continue
      return
      end
