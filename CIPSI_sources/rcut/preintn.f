      subroutine preintn
      implicit double precision (a-h,o-z)
      integer casp,pmax2a,pmax2b,pmaxab,pmaxabi,pmaxa,pmaxb,
     1        pmaxa12,pmaxb12,pl
      logical*1 iandj
      parameter(itypmax=3)
      common/norb/nf(0:itypmax)
      common/parinv/xi,yi,zi,xj,yj,zj,cx,cy,cz,lit,ljt,i1,j1,i2,j2,
     1              mini,maxi,mazi,minj,maxj,mazj,loci,locj
      common/project/pl
      common/nocas/casp
      common/ityp/itypg,itypd
      common/vintnum/pmax2a,pmax2b,pmaxab,pmaxabi,pmaxa,pmaxb
      common/nexpres/nint,noij
      common/dist/a(3),b(3),ac,bc
      common/distpar/ab(3),pab
      common/varvpl/vpl(0:3)
      common /viandj/iandj
      dimension pmaxa12(0:3,10),pmaxb12(0:3,10)
      data pmaxa12/1,1,1,1,2,2,3,3,2,3,3,4,3,3,4,4,3,4,4,5,3,4,5,5,
     1             4,4,5,5,4,5,5,6,4,5,6,6,4,5,6,7/
      data pmaxb12/1,1,1,1,1,2,3,2,2,3,3,4,1,2,4,3,2,3,4,4,3,4,5,5,
     1             1,2,5,4,2,3,5,5,3,4,6,6,4,5,6,7/
 
c
c
c
      itypg=lit
      itypd=ljt
c
c
c
      nint=nf(itypg)*nf(itypd)
      if(iandj) nint=(nf(itypg)*(nf(itypg)+1))/2
      noij=(itypg+1)*itypg/2+itypd+1
c
c
c
      pmax2a=itypg+min(pl,itypd)+1
      pmax2b=itypd+min(pl,itypg)+1
      pmaxab=itypg+itypd
      pmaxabi=1
      pmaxa=pmaxa12(pl,noij)
      pmaxb=pmaxb12(pl,noij)
c
c
c
      casp=0
      a(1)=xi-cx
      a(2)=yi-cy
      a(3)=zi-cz
      b(1)=xj-cx
      b(2)=yj-cy
      b(3)=zj-cz
c
c
      ac=dsqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      bc=dsqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))
c
c
      if(ac.lt.1.0d-09 .and. bc.lt.1.0d-09)then
                casp=3
                elseif(ac.lt.1.0d-09)then
                          casp=1
                          ab(1)=b(1)
                          ab(2)=b(2)
                          ab(3)=b(3)
                          pab=dsqrt(ab(1)*ab(1)+ab(2)*ab(2)+ab(3)*ab(3))
                elseif(bc.lt.1.0d-09)then
                          casp=2
                          ab(1)=a(1)
                          ab(2)=a(2)
                          ab(3)=a(3)
                          pab=dsqrt(ab(1)*ab(1)+ab(2)*ab(2)+ab(3)*ab(3))
      endif
c
c
c
      if(casp.eq.0)then
      ctta=(a(1)*b(1)+a(2)*b(2)+a(3)*b(3))/(ac*bc)
      pl1=ctta
      ctaa2=pl1*pl1
      pl2=((3.d0*ctaa2)-1.d0)/2.d0
      pl3=(2.5d0*ctaa2*pl1)-(1.5d0*pl1)
      vpl(0)=1.d0
      vpl(1)=pl1
      vpl(2)=pl2
      vpl(3)=pl3
      endif
c
c
c
      return
      end
