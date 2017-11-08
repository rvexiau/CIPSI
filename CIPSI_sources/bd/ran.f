C---------------------------------------------------------------------
      function ran(idum)
      implicit real*8(a-h,m-z)
      parameter(mbig=4000000.,mseed=161803398.,mz=0.,fac=1./mbig)
      common/cran/ ma(55),inext,inextp
      data iff/0/
      if (idum.lt.0.or.iff.eq.0) then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
           ii=mod(21*i,55)
           ma(ii)=mk
           mk=mj-mk
           if (mk.lt.mz) mk=mk+mbig
           mj=ma(ii)
11       continue
        do 13 k=1,4
        do 12 i=1,55
           ma(i)=ma(i)-ma(1+mod(i+30,55))
           if (ma(i).lt.mz) ma(i)=ma(i)+mbig
12       continue
13       continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
cdudu      if (inext.eq.56) inext=1
      if (inext.eq.56) inext=2
      inextp=inextp+1
      if (inextp.eq.56) inextp=1
cdudu     if (inextp.eq.56) inextp=2
      mj=ma(inext)-ma(inextp)
      if (mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      bb=mj*fac
cdudu
      idum=-int(bb*1d6)
      ranr=bb
      return
      end
      
