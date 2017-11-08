      subroutine qout
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 iandj,kandl,same,out                                                  
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
      common/iofile/ir,iw,ip,is
      common/shlt/tol,cutoff,icount,out
      common/output/nprint,itol,icut,normf,normp
      common/shlnos/qq4,ii,jj,kk,ll,lit,ljt,lkt,llt,
     1 mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     1 loci,locj,lock,locl,ij,kl,ijkl,nij
      common/misc/maxll,iandj,kandl,same
      common/gout/g(10000)
      common/buf/xx(lsize),ix(lsize)
      common/wrtc/ i1,i2,i3,i4,inx(4,5),q4,val,q4x(5),valx(5),nwrt
      dimension mmaz(10)
      data mmaz /1,4,9,17,29,1,4,10,20,35/

      if((lit.ge.3.and.lit.le.5).or.(ljt.ge.3.and.ljt.le.5).or.
     1   (lkt.ge.3.and.lkt.le.5).or.(llt.ge.3.and.llt.le.5))call comtwo
c
c     ----- pack the 4 indices of integral into one word
c     ----- write label + integral on tape (is)
c
      mazi=mmaz(lit)
      mazj=mmaz(ljt)
      mazk=mmaz(lkt)
      mazl=mmaz(llt)
c
      label=0
      nn=0
      ijn=0
      jmax=mazj
      do 1700 i=mini,mazi
      if(iandj) jmax=i
      do 1600 j=minj,jmax
      ijn=ijn+1
      lmax=mazl
      kln=0
      do 1500 k=mink,mazk
      if(kandl) lmax=k
      do 1400 l=minl,lmax
      kln=kln+1
      if(same.and.kln.gt.ijn) go to 1600
      nn=nn+1
      val=g(nn)
      if(abs(val).lt.cutoff) go to 1400
      i1=loci+i
      i2=locj+j
      i3=lock+k
      i4=locl+l
      if(i1.ge.i2) go to 700
      n=i1
      i1=i2
      i2=n
  700 if(i3.ge.i4) go to 800
      n=i3
      i3=i4
      i4=n
  800 if(i1-i3) 900,1000,1100
  900 n=i1
      i1=i3
      i3=n
      n=i2
      i2=i4
      i4=n
      go to 1100
 1000 if(i2.lt.i4) go to 900
 1100 continue
c
c     ----- compute q factor for each individual integral
c
      q4=qq4
      if(i1.eq.i2) q4=q4/2.0d+00
      if(i3.eq.i4) q4=q4/2.0d+00
      if(i1.eq.i3.and.i2.eq.i4) q4=q4/2.0d+00
      val=val*q4
      if(out) call wrt(0)
c
      xx(icount)=val
      label=iword(i1,i2,i3,i4)
c      write(6,*) i1,i2,i3,i4,label 
      ix(icount)=label
      icount=icount+1
c      write(6,'(2i12,4i4,e14.6)') ,icount,label,i1,i2,i3,i4,val
      if(icount.le.lsize) go to 1400
      write(is) xx,ix
      icount=1
      nrec=nrec+1
 1400 continue
 1500 continue
 1600 continue
 1700 continue
      return
      end

