      subroutine rhr
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      common/symtry/invt(48),iso(ds,12),ptr(3,144),dtr(6,288)
     1,ftr(10,480),nt
      common/hsym/t(35,35),mink,maxk,lkt,minl,maxl,llt,ntr
      dimension v(35)
c
c     ----- right multiply  t  by  r,
c           result back in  t
c
    6 go to (500,400,300,200,100,500,400,300,200,100),llt
c
c     ----- g shell
c
c 100 ng=15*(ntr-1)-20
c     do 130 k=mink,maxk
c     do 120 l=minl,maxl
c     dum=0.0d+00
c     do 110 n=minl,maxl
c 110 dum=dum+t(k,n)*gtr(n-20,ng+l)
c 120 v(l)=dum
c     do 130 l=minl,maxl
c 130 t(k,l)=v(l)
c     go to 500
  100 go to 500
c
c     ----- f shell
c
  200 nf=10*(ntr-1)-10
      do 230 k=mink,maxk
      do 220 l=minl,maxl
      dum=0.0d+00
      do 210 n=minl,maxl
  210 dum=dum+t(k,n)*ftr(n-10,nf+l)
  220 v(l)=dum
      do 230 l=minl,maxl
  230 t(k,l)=v(l)
      go to 500
c
c     ----- d shell
c
  300 nd= 6*(ntr-1)- 4
      do 330 k=mink,maxk
      do 320 l=minl,maxl
      dum=0.0d+00
      do 310 n=minl,maxl
  310 dum=dum+t(k,n)*dtr(n-4,nd+l)
  320 v(l)=dum
      do 330 l=minl,maxl
  330 t(k,l)=v(l)
      go to 500
c
c     ----- p shell
c
  400 np= 3*(ntr-1)- 1
      do 430 k=mink,maxk
      do 420 l=2,4
      dum=0.0d+00
      do 410 n=2,4
  410 dum=dum+t(k,n)*ptr(n-1,np+l)
  420 v(l)=dum
      do 430 l=2,4
  430 t(k,l)=v(l)
  500 continue
c
c     ----- left multiply  t  by r
c           result back in  t
c
   38 go to (1000,900,800,700,600,1000,900,800,700,600),lkt
c
c     ----- g shell
c
c 600 ng=15*(ntr-1)-20
c     do 630 l=minl,maxl
c     do 620 k=mink,maxk
c     dum=0.0d+00
c     do 610 n=mink,maxk
c 610 dum=dum+gtr(n-20,ng+k)*t(n,l)
c 620 v(k)=dum
c     do 630 k=mink,maxk
c 630 t(k,l)=v(k)
c     go to 1000
  600 go to 1000
c
c     ----- f shell
c
  700 nf=10*(ntr-1)-10
      do 730 l=minl,maxl
      do 720 k=mink,maxk
      dum=0.0d+00
      do 710 n=mink,maxk
  710 dum=dum+ftr(n-10,nf+k)*t(n,l)
  720 v(k)=dum
      do 730 k=mink,maxk
  730 t(k,l)=v(k)
      go to 1000
c
c     ----- d shell
c
  800 nd= 6*(ntr-1)-4
      do 830 l=minl,maxl
      do 820 k=mink,maxk
      dum=0.0d+00
      do 810 n=mink,maxk
  810 dum=dum+dtr(n-4,nd+k)*t(n,l)
  820 v(k)=dum
      do 830 k=mink,maxk
  830 t(k,l)=v(k)
      go to 1000
c
c     ----- p shell
c
  900 np= 3*(ntr-1)- 1
      do 930 l=minl,maxl
      do 920 k=2,4
      dum=0.0d+00
      do 910 n=2,4
  910 dum=dum+ptr(n-1,np+k)*t(n,l)
  920 v(k)=dum
      do 930 k=2,4
  930 t(k,l)=v(k)
 1000 continue
      return
      end
