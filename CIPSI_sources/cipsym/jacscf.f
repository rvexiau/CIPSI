      subroutine jacscf(a,b,c,naa,nqq,epslon)
c********************************************
      implicit real*8(a-h,o-p,r-z)
c
c          find all eigenvalues and eigenvectors of a(*) by the
c          jacobi method.
c
c          this program treats the orthogonal case (s=1)
c      a(*) - input matrix to be diagonalized. (this matrix
c             is stored in lower triangular form and is
c             is destroyed during the calculation)
c      b(*) - upon exiting this routine this array contains the
c             eigenvectors stored in a rectangular compact form
c      c(*) - upon exiting this routine this array contains the
c             eigenvalues.
c             naa=dimension of a
c             nq=-1 causes b to be cleared to a unit matrix
c             nq>0 b(*) contains an naa by nq array which will be
c                  multiplied by the eigenvector matrix.
c             nq should not be dimensioned<naa if nq>0.
c            routine.
c      epslon is the convergence criteria for off diagonal elements
      dimension a(1),b(1),c(1)
      data zero/0.0d0/, half/0.5d0/, one/1.0d0/, two/2.0d0/,
     *     three/3.0d0/,scale/0.909090909/,sc2/0.95238095/
      nq=nqq
      na=naa
      nam=na-1
c   set b to unit matrix if nq .le. 0
      if(nq.gt.0) go to 120
         k=1
         nq=na
         do 115 i=1,na
            do 110 j=1,na
               if(i.eq.j) go to 100
                  b(k)=zero
                  go to 105
100            continue
               b(k)=one
105            continue
               k=k+1
110         continue
115      continue
c
c      --- determine initial and final thresholds ---
c
120   continue
      if(nam)325,310,125
125   continue
      k=1
      sum1=zero
      do 155 i=2,na
         im=i-1
         do 150 j=1,im
            sum1=sum1+abs(a(k))
            k=k+1
150      continue
         k=k+1
155   continue
      t=(na*(na-1))/2
      sum1=sum1/t
      k=1
      amax=a(1)
      amin=a(1)
      do 160 j=2,na
         k=k+j
         if(a(k).gt.amax) amax=a(k)
         if(a(k).lt.amin) amin=a(k)
160   continue
      t=(amax-amin)/dble(na)+sum1*two
      thrshg=epslon*t
      thresh=sum1
      if(thresh.lt.thrshg) thresh=thrshg
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< begin program loop 180
180   continue
      k=1
      n=0
      jd=1
      kd=0
      do 270 j=2,na
         jd=jd+j
         jj=j-1
         kd=kd+nq
         do 265 i=1,jj
c
c      --- test for large off diagonal element ---
c
            if(abs(a(k+i)).le.thresh) go to 265
c
c      --- compute rotation matrix ---
c
               n=n+1
               id=i*(i+1)/2
               alpha=(a(jd)-a(id))/(two*a(k+i))
               t=sign(one/(abs(alpha)+sqrt(one+alpha*alpha)),alpha)
               ca=sqrt(one+t*t)
               cc=one/ca
               s=-cc*t
c
c      --- apply rotation to elements i,i  i,j  j,j ---
c
               a(jd)=a(jd)+a(k+i)*t
               a(id)=a(id)-a(k+i)*t
               a(k+i)=zero
               ka=jd-j
               kb=id-i
               kc=(i-1)*nq
c
c      --- apply rotation to vector matrix ---
c
               do 215 l=1,nq
                  b(kd+l)=-s*b(kc+l)+cc*b(kd+l)
                  b(kc+l)=ca*b(kc+l)-t*b(kd+l)
215            continue
c
c      --- apply rotation to rest of operator matrix ---
c
               l1=i-1
               l3=j-2
               do 220 l=1,l1
                  a(ka+l)=-s*a(kb+l)+cc*a(ka+l)
                  a(kb+l)=ca*a(kb+l)-t*a(ka+l)
220            continue
               ka=ka+1
               kb=kb+i
               do 225 l=i,l3
                  kb=kb+l
                  a(ka+l)=-s*a(kb)+cc*a(ka+l)
                  a(kb)= ca*a(kb)-t*a(ka+l)
225            continue
               ka=ka+j-1
               kb=kb+j-1
               do 230 l=j,nam
                  ka=ka+l
                  kb=kb+l
                  a(ka)=-s*a(kb)+cc*a(ka)
                  a(kb)=ca*a(kb)-t*a(ka)
230            continue
265      continue
         k=k+j
270   continue
c
c      --- test for convergence (thresh=thrshg and n=0) ---
c
      if(thresh.eq.thrshg) go to 285
         thresh=thresh*scale
         scale=scale*sc2
         if(thresh.gt.thrshg) go to 290
         thresh=thrshg
         go to 290
  285 if(n.eq.0) go to 310
  290  continue
      go to 180
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> end program loop 180
  310  continue
      ll=0
315   continue
      do 320 l=1,na
         ll=ll+l
         c(l)=a(ll)
320   continue
325   continue
      return
      end
