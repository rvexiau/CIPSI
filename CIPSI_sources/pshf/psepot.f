      subroutine psepot(nprim)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      logical pseud
      character*80 apseud
      common/psnloc/apseud(dc),iatno(dc)
      common /infpot/apot(400),cpot(400),npot(400),nbtyp(200),
     1 ipseud(103),nsom(50),pseud
      common/uncp/ax,ay,az,bx,by,bz,xc,yc,zc,ga,gb,gc,am,bm,cm,
     1 pt,alfa,alfb,alfc,ctta,r,f1,f2,f3,f4,f5,f6,f7,f8,f9,xbs,bs(10),
     2 xg(64),hg(64),rxg(40),h1(40),h23(40),h456(40),h78(40),h9(40),
     3 fgam(13),rab,l,n,k,inum,nbs
      common/pmat/ psdint(225),pspt(225)
      common/pdpot/ xaa,yaa,zaa,icaa,itaa,aa,
     1              xbb,ybb,zbb,icbb,itbb,bb,
     2              xnc,ync,znc,ic
      data pi4/0.1256637061435917d+2/
c......
      do 5 i=1,nprim
    5 pspt(i)=0.d0
      xc=xnc
      yc=ync
      zc=znc
      if(itaa.ge.itbb) go to 10
      alfa=bb
      alfb=aa
      ita=itbb
      itb=itaa
      ica=icbb
      icb=icaa
      xa=xbb
      ya=ybb
      za=zbb
      xb=xaa
      yb=yaa
      zb=zaa
      go to 20
   10 alfa=aa
      alfb=bb
      ita=itaa
      itb=itbb
      ica=icaa
      icb=icbb
      xa=xaa
      ya=yaa
      za=zaa
      xb=xbb
      yb=ybb
      zb=zbb
   20 inum=(ita-1)*ita/2+itb
      rab=alfa+alfb
      if(ica.eq.ic) go to 30
      if(icb.eq.ic) go to 40
      k=1
      ax =xa-xc
      ay= ya-yc
      az= za-zc
      pa=dsqrt(ax*ax+ay*ay+az*az)
      bx= xb-xc
      by= yb-yc
      bz= zb-zc
      pb=dsqrt(bx*bx+by*by+bz*bz)
      xc=(alfa*ax+alfb*bx)/rab
      yc=(alfa*ay+alfb*by)/rab
      zc=(alfa*az+alfb*bz)/rab
      pc=dsqrt(xc*xc+yc*yc+zc*zc)
      gaa=alfa*pa
      gaa=gaa+gaa
      gbb=alfb*pb
      gbb=gbb+gbb
      gcc=rab*pc
      gcc=gcc+gcc
      pt=alfa*pa*pa+alfb*pb*pb
      ctta=(ax*bx+ay*by+az*bz)/(pa*pb)
      go to 60
   30 if(icb.eq.ic) go to 50
      k=2
      bx =xb-xc
      by= yb-yc
      bz= zb-zc
      pb=dsqrt(bx*bx+by*by+bz*bz)
      xc=alfb*bx/rab
      yc=alfb*by/rab
      zc=alfb*bz/rab
      pc=dsqrt(xc*xc+yc*yc+zc*zc)
      gbb=alfb*pb
      gbb=gbb+gbb
      gcc=rab*pc
      gcc=gcc+gcc
      pt=alfb*pb*pb
      ctta=1.0
      go to 60
   40 k=3
      ax= xa-xc
      ay= ya-yc
      az =za-zc
      pa=dsqrt(ax*ax+ay*ay+az*az)
      xc=alfa*ax/rab
      yc=alfa*ay/rab
      zc=alfa*az/rab
      pc=dsqrt(xc*xc+yc*yc+zc*zc)
      gaa=alfa*pa
      gaa=gaa+gaa
      gcc=rab*pc
      gcc=gcc+gcc
      pt=alfa*pa*pa
      ctta=1.0
      go to 60
   50 k=4
      ctta=1.0
   60 continue
      iato=iatno(ic)
      iset=ipseud(iato)
      idx=nsom(iset)
      do 190 np=1,4
      itp=(iset-1)*4+np
      itpmax=nbtyp(itp)
      if(itpmax.eq.0) go to 190
      do 180 ipo=1,itpmax
      idx=idx+1
      alfc=apot(idx)
      cof=cpot(idx)
      n=npot(idx)
      l=np-2
      r=rab+alfc
      rr=dsqrt(r)
      go to (70,80,90,100),k
   70 ga=gaa/rr
      am=pa*rr
      gb=gbb/rr
      bm=pb*rr
      gc=gcc/rr
      cm=pc*rr
      go to 100
   80 gb=gbb/rr
      bm=pb*rr
      gc=gcc/rr
      cm=pc*rr
      go to 100
   90 ga=gaa/rr
      am=pa*rr
      gc=gcc/rr
      cm=pc*rr
  100 cte=pi4/(rr**(3-n))
      if(l.gt.0) cte=cte*dble(2*l+1)
      go to (110,120,130,140,150,160),inum
  110 call pgspls
      go to 170
  120 call pgppls
      go to 170
  130 call pgpplp
      go to 170
  140 call pgdpls
      go to 170
  150 call pgdplp
      go to 170
  160  call pgdpld
  170 continue
      do 175 i=1,nprim
  175 pspt(i)=pspt(i)+cof*cte*psdint(i)
  180 continue
  190 continue
      return
      end
