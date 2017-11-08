      subroutine comtwo
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 iandj,kandl,same,out
      common/shlnos/qq4,ii,jj,kk,ll,lit,ljt,lkt,llt,mini,minj,mink,minl,
     1 maxi,maxj,maxk,maxl,loci,locj,lock,locl,ij,kl,ijkl,nij
      common/misc/maxll,iandj,kandl,same
      common/gout/g(10000)
      dimension mmax(10),mmaz(10),tab(15,15)
      dimension g1(5050)
      data mmax/1,3,6,10,15,1,4,6,10,15/,mmaz/1,3,5,7,9,1,4,6,10,15/
      imax=mmax(lit)
      jmax=mmax(ljt)
      kmax=mmax(lkt)
      lmax=mmax(llt)
      imaz=mmaz(lit)
      jmaz=mmaz(ljt)
      kmaz=mmaz(lkt)
      lmaz=mmaz(llt)
c
      if(.not.same) go to 100
c
      do 10 inv=1,ijkl
   10 g1(inv)=g(inv)
      nn=0
      nnp=0
      jm=jmax
      ijn=0
      do 40 i=1,imax
      if(iandj) jm=i
      do 40 j=1,jm
      ijn=ijn+1
      lm=lmax
      kln=0
      do 30 k=1,kmax
      if(kandl) lm=k
      do 30 l=1,lm
      kln=kln+1
      nn=nn+1
      if(kln-ijn) 21,20,22
   20 m=nnp+1
   21 nnp=nnp+1
      g(nn)=g1(nnp)
      go to 30
   22 g(nn)=g1(m+ijn)
      m=m+kln
   30 continue
   40 continue
c
c
  100 continue
      if(lkt.le.2.and.ljt.le.2) go to 210
      nn=0
      nm=0
      jm=jmax
      do 200 i=1,imax
      if(iandj) jm=i
      do 200 j=1,jm
      do 150 k=1,kmax
      if(kandl) go to 130
      do 120 l=1,lmax
      nn=nn+1
  120 tab(k,l)=g(nn)
      go to 150
  130 do 140 l=1,k
      nn=nn+1
      tab(k,l)=g(nn)
  140 tab(l,k)=g(nn)
  150 continue
c
      call comtab(lkt,llt,kandl,tab)
c
      lm=lmaz
      do 160 k=1,kmaz
      if(kandl) lm=k
      do 160 l=1,lm
      nm=nm+1
  160 g1(nm)=tab(k,l)
  200 continue
      go to 230
  210 continue
      nijkl=ij*kl
      do 220 i=1,nijkl
  220 g1(i)=g(i)
  230 continue
c
      if((lit.le.2.or.lit.gt.5).and.(ljt.le.2.or.ljt.gt.5)) go to 500
      nkl=0
      lkl=kmaz*lmaz
      if(kandl) lkl=(lkl+kmaz)/2
      lm=lmaz
      do 400 k=1,kmaz
      if(kandl) lm=k
      do 400 l=1,lm
      nkl=nkl+1
      nn=nkl-lkl
      do 350 i=1,imax
      if(iandj) go to 330
      do 320 j=1,jmax
      nn=nn+lkl
  320 tab(i,j)=g1(nn)
      go to 350
  330 do 340 j=1,i
      nn=nn+lkl
      tab(i,j)=g1(nn)
  340 tab(j,i)=g1(nn)
  350 continue
c
      call comtab(lit,ljt,iandj,tab)
c
      nn=nkl-lkl
      jm=jmaz
      do 360 i=1,imaz
      if(iandj) jm=i
      do 360 j=1,jm
      nn=nn+lkl
  360 g(nn)=tab(i,j)
  400 continue
c
      if(.not.same) return
c
      nn=0
      nm=0
      ijn=0
      jm=jmaz
      do 710 i=1,imaz
      if(iandj) jm=i
      do 710 j=1,jm
      ijn=ijn+1
      kln=0
      lm=lmaz
      do 720 k=1,kmaz
      if(kandl) lm=k
      do 720 l=1,lm
      kln=kln+1
      nn=nn+1
      if(kln.gt.ijn) go to 720
      nm=nm+1
      g1(nm)=g(nn)
  720 continue
  710 continue
c
  500 continue
      do 510 i=1,nm
  510 g(i)=g1(i)
      return
      end
