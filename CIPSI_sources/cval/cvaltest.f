      program cvaltest
      implicit double precision (a-h,o-z)
      include 'pshf.prm'
      dimension x1(lsize),x2(lsize)
      integer ix1(lsize),ix2(lsize)
      open(11,file='/disc3/TEMP/gross/noar_pqrs_cv',
     &      form='unformatted',status='OLD') 
      open(12,file='/disc3/TEMP/gross/noar_pqrs_cvo',
     &      form='unformatted',status='OLD') 
      err =0.0
      n=0
      nd=0
1     read (11,end=999) x1,ix1
      read (12) x2,ix2
      do i=1,lsize
        if (ix1(i).ne.ix2(i)) then
	  nd=nd+1
c  print *,"ix1 and ix2 different:",ix1(i),ix2(i)
        end if
        err = err + (x1(i)-x2(i))**2
        n=n+1
      end do
      goto 1
999   continue
      print *,"number of lines:", n
      print *,"total error:", err,nd
      stop
      end
