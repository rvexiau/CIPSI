      subroutine s0000
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 iandj,kandl,same,out
      common/shlinf/ag(20),csa(20),cpa(20),cda(20),cfa(20),cga(20),
     1              bg(20),csb(20),cpb(20),cdb(20),cfb(20),cgb(20),
     1              cg(20),csc(20),cpc(20),cdc(20),cfc(20),cgc(20),
     1              dg(20),csd(20),cpd(20),cdd(20),cfd(20),cgd(20),
     1              xi,yi,zi,xj,yj,zj,rri,xk,yk,zk,xl,yl,zl,rrk,
     1              nga,ngb,ngc,ngd
      common/shlt/tol,cutoff,icount,out
      common/ijpair/a(400),r(400),x1(400),y1(400),z1(400),dij(6400),
     1 ijd(225)
      common/shlnos/qq4,ii,jj,kk,ll,lit,ljt,lkt,llt,
     1 mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     1 loci,locj,lock,locl,ij,kl,ijkl,nij
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
      common/iofile/ir,iw,ip,is
      common/buf/ tei(lsize),ix(lsize)
      common/gout/gout
      common/misc/maxll,iandj,kandl,same
      common/wrtc/ in,jn,kn,ln,inx(4,5),q4,val ,q4x(5),valx(5),nwrt
      data pi252 /34.986836655250d+00/
      data pie4 /7.85398163397448d-01/

      label=0
      lgmax=ngd
      do 205 kg=1,ngc
      bk=cg(kg)
      brrk=bk*rrk
      bxk=bk*xk
      byk=bk*yk
      bzk=bk*zk
      csk=csc(kg)
      if(kandl) lgmax=kg
      do 200 lg=1,lgmax
      bl=dg(lg)
      bb=bk+bl
      dum=bl*brrk/bb
      if(dum.gt.tol) go to 200
      bbrrk=dum
      d2=csd(lg)*csk/bb
      if(kandl.and.lg.ne.kg) d2=d2+d2
      bbx=(bxk+bl*xl)/bb
      bby=(byk+bl*yl)/bb
      bbz=(bzk+bl*zl)/bb
      sum=0.0d+00
      nn=1
      do 150 n=1,nij
      dum=bbrrk+r(n)
      if(dum.gt.tol) go to 150
      expe=exp(-dum)
      aa=a(n)
      ab=aa+bb
      dum=x1(n)-bbx
      xx=dum*dum
      dum=y1(n)-bby
      xx=dum*dum+xx
      dum=z1(n)-bbz
      xx=dum*dum+xx
      x=xx*aa*bb/ab
      if(x.gt.5.0d+00) go to 50
      if(x.gt.1.0d+00) go to 30
      if(x.gt.3.0d-7) go to 20
      ww1=1.0d+00-x/3.0d+00
      go to 100
   20 f1=          ((((((((-8.36313918003957d-08*x+1.21222603512827d-06
     1)*x-1.15662609053481d-05 )*x+9.25197374512647d-05
     2)*x-6.40994113129432d-04 )*x+3.78787044215009d-03
     3)*x-1.85185172458485d-02 )*x+7.14285713298222d-02
     4)*x-1.99999999997023d-01 )*x+3.33333333333318d-01
      ww1=(x+x)*f1+exp(-x)
      go to 100
   30 if(x.gt.3.0d+00) go to 40
      y=x-2.0d+00
      f1=        ((((((((((-1.61702782425558d-10*y+1.96215250865776d-09
     1)*y-2.14234468198419d-08 )*y+2.17216556336318d-07
     2)*y-1.98850171329371d-06 )*y+1.62429321438911d-05
     3)*y-1.16740298039895d-04 )*y+7.24888732052332d-04
     4)*y-3.79490003707156d-03 )*y+1.61723488664661d-02
     5)*y-5.29428148329736d-02 )*y+1.15702180856167d-01
      ww1=(x+x)*f1+exp(-x)
      go to 100
   40 y=x-4.0d+00
      f1=        ((((((((((-2.62453564772299d-11*y+3.24031041623823d-10
     1)*y-3.614965656163d-09)*y+3.760256799971d-08)*y-3.553558319675d-07
     2)*y+3.022556449731d-06)*y-2.290098979647d-05)*y+1.526537461148d-04
     3)*y-8.81947375894379d-04 )*y+4.33207949514611d-03
     4)*y-1.75257821619926d-02 )*y+5.28406320615584d-02
      ww1=(x+x)*f1+exp(-x)
      go to 100
   50 if(x.gt.15.0d+00) go to 70
      e=exp(-x)
      if(x.gt.10.0d+00) go to 60
      ww1=    (((((( 4.6897511375022d-01/x-6.9955602298985d-01)/x
     1+5.3689283271887d-01)/x-3.2883030418398d-01)/x
     2+2.4645596956002d-01)/x-4.9984072848436d-01)/x
     3-3.1501078774085d-06)*e + sqrt(pie4/x)
      go to 100
   60 ww1=       (((-1.8784686463512d-01/x+2.2991849164985d-01)/x
     1-4.9893752514047d-01)/x-2.1916512131607d-05)*e + sqrt(pie4/x)
      go to 100
   70 if(x.gt.33.0d+00) go to 90
      e=exp(-x)
      ww1=        (( 1.9623264149430d-01/x-4.9695241464490d-01)/x
     1-6.0156581186481d-05)*e + sqrt(pie4/x)
      go to 100
   90 ww1=sqrt(pie4/x)
  100 sum=sum+dij(nn)*ww1*expe/sqrt(ab)
  150 nn=nn+16
      gout=gout+d2*sum
  200 continue
  205 continue
      gout=gout*pi252
      if(abs(gout).lt.cutoff) return
      in=loci+1
      jn=locj+1
      kn=lock+1
      ln=locl+1
      if(in.lt.jn) then
      idum=in
      in=jn
      jn=idum
      endif
      if(kn.lt.ln) then
      idum=kn
      kn=ln
      ln=idum
      endif
      if(in.lt.kn) then
      idum=in
      in=kn
      kn=idum
      idum=jn
      jn=ln
      ln=idum
      endif
      if(in.eq.kn.and.jn.lt.ln) then
      idum=jn
      jn=ln
      ln=idum
      endif
      q4=qq4
      if(in.eq.jn) q4=q4/2.0d+00
      if(kn.eq.ln) q4=q4/2.0d+00
      if(in.eq.kn.and.jn.eq.ln) q4=q4/2.0d+00
      gout=gout*q4
      val=gout
      if(out) call wrt(0)
c
c
      tei(icount)=gout
      label=iword(in,jn,kn,ln)
      ix(icount)=label
      icount=icount+1
      if(icount.le.lsize) return
      write(is) tei,ix
      icount=1
      nrec=nrec+1
      return
      end
