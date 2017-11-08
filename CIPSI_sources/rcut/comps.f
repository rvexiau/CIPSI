      subroutine comps (lit,ljt,iandj)
      implicit real*8(a-h,o-z)
      logical*1 iandj
      common/recomb/v(225,20),nvar,itpmax
      dimension mmax(10),mmaz(10),tab(15,15)
      data mmax/1,3,6,10,15,1,4,6,10,15/,mmaz/1,3,5,7,9,1,4,6,10,15/
c
      imax=mmax(lit)
      jmax=mmax(ljt)
      imaz=mmaz(lit)
      jmaz=mmaz(ljt)
c
      do 50 k=nvar,itpmax
c
      nn=0
      jm=jmax
      do 30 i=1,imax
      if(iandj) go to 15
      do 10 j=1,jmax
      nn=nn+1
   10 tab(i,j)=v(nn,k)
      go to 30
   15 do 20 j=1,i
      nn=nn+1
      tab(i,j)=v(nn,k)
   20 tab(j,i)=v(nn,k)
   30 continue
c
      call comtab(lit,ljt,iandj,tab)
c
      jm=jmaz
      nn=0
      do 40 i=1,imaz
      if(iandj) jm=i
      do 40 j=1,jm
      nn=nn+1
   40 v(nn,k)=tab(i,j)
c
   50 continue
      return
      end
