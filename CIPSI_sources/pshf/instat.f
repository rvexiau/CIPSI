      subroutine instat(istat)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      common/open/nc,noc(3),oc(3),alpha(3,3),beta(3,3),rland
      common/pseudo/qn(doa)
      common/infoa/nat,ich,mul,n,nx,ne,na
      common/iofile/ir,iw
      oc(1)=2.0d0
      alpha(1,1)=2.0d0
      beta(1,1)=1.0d0
      if(istat.eq.1) go to 700
      oc(2)=1.0d0
      alpha(1,2)=1.0d0
      beta(1,2)=0.5d0
 700  nn=0
      do 750 i=1,na
      if(qn(i).ne.2) go to 800
 750  nn=nn+1
 800  noc(1)=nn
      imul=na-nn+1
      if(istat.ne.0) go to (1,2,3,4,5,6,7,8,8,10,10,12,13,13),istat
      write(iw,1000)
      if(imul.eq.1) go to 1
      if(imul.eq.2) go to 3
      if(imul.eq.(mul+1)) then
         write(iw,1100)
         go to 4
         else
         if(mul.eq.0.and.imul.eq.3) go to 2
         write(iw,1200)
         stop
      endif
c
c *** closed shell
c
 1    nc=1
      write(iw,2001) nc
      write(iw,2999)
      i=1
      j=na
      write(iw,3001) oc(1),i,j,alpha(1,1),beta(1,1)
      go to 500
c
c *** excited singlet state
c
 2    nc=3
      alpha(1,3)=1.0d0
      alpha(2,3)=0.5d0
      beta(1,3)=0.5d0
      beta(2,3)=-0.5d0
      oc(3)=1.0d0
      noc(2)=1
      noc(3)=1
      write(iw,2002) nc
      go to 100
c
c *** doublet state (non degenerate mos)
c
 3    nc=2
      noc(2)=1
      write(iw,2003) nc
      go to 100
c
c *** triplet state and other highest multiplicity half filled open
c *** shell configurations
c
 4    nc=2
      noc(2)=na-noc(1)
      alpha(2,2)=0.5d0
      beta(2,2)=0.5d0
      write(iw,2004) imul,nc
      go to 100
c
c *** two degenerate mos, one electron. doublet state
c
 5    nc=2
      oc(2)=0.5d0
      alpha(1,2)=0.5d0
      beta(1,2)=0.25d0
      noc(2)=2
      na=noc(1)+2
      write(iw,2005) nc
      go to 100
c
c *** two degenerate mos, three electrons. doublet state
c
 6    nc=2
      oc(2)=1.5d0
      noc(2)=2
      na=noc(1)+2
      alpha(1,2)=1.5d0
      alpha(2,2)=1.d0
      beta(1,2)=0.75d0
      beta(2,2)=0.5d0
      write(iw,2006) nc
      go to 100
c
c *** three degenerate mos, one electron
c
 7    nc=2
      oc(2)=0.333333333d0
      alpha(1,2)=0.333333333d0
      beta(1,2)=0.166666667d0
      noc(2)=3
      na=noc(1)+3
      write(iw,2007) nc
      go to 100
c
c *** one electron in a non degenerate mo, one in two degenerate mos
c *** of higher energy. singlet or triplet states
c
 8    nc=3
      oc(3)=0.5d0
      alpha(1,3)=0.5d0
      alpha(2,3)=0.25d0
      noc(2)=1
      noc(3)=2
      na=noc(1)+3
      beta(1,3)=0.25d0
c
c *** singlet state
c
      if(istat.eq.8) then
         beta(2,3)=-0.25d0
         write(iw,2008) nc
c
c *** triplet state
c
         else
         beta(2,3)=0.25d0
         write(iw,2009) nc
      endif
      go to 100
c
c *** three electrons in two degenerate mos, one electron in a
c *** non degenerate mo of higher energy. singlet or triplet states
c
 10   nc=3
      oc(2)=1.5d0
      oc(3)=1.0d0
      alpha(1,2)=1.50d0
      alpha(1,3)=1.0d0
      alpha(2,2)=1.0d0
      alpha(2,3)=0.75d0
      beta(1,2)=0.75d0
      beta(1,3)=0.5d0
      beta(2,2)=0.5d0
      noc(2)=2
      noc(3)=1
      na=noc(1)+3
c
c *** singlet
c
      if(istat.eq.10) then
         write(iw,2010) nc
         beta(2,3)=0.0d0
         else
c
c *** triplet
c
         write(iw,2011) nc
         beta(2,3)=0.5d0
      endif
      go to 100
c
c *** three degenerate mos, two electrons. triplet state
c
 12   nc=2
      oc(2)=0.666666667d0
      alpha(1,2)=0.666666667d0
      alpha(2,2)=0.166666667d0
      beta(1,2)=0.333333333d0
      beta(2,2)=0.166666667d0
      noc(2)=3
      na=noc(1)+3
      write(iw,2012) nc
      go to 100
c
c *** one electron in a non degenerate mo, one in three degenerate mos
c *** of higher energy. singlet or triplet state
c
 13   nc=3
      oc(2)=1.d0
      oc(3)=0.333333333d0
      noc(2)=1
      noc(3)=3
      na=noc(1)+4
      alpha(1,3)=0.333333333d0
      alpha(2,3)=0.166666667d0
      beta(1,3)=0.166666667d0
      if(istat.eq.13) then
c
c *** singlet
c
         beta(2,3)=-.166666667d0
         write(iw,2013) nc
         else
c
c *** triplet
c
         beta(2,3)=0.166666667d0
         write(iw,2014) nc
      endif
 100  do 600 i=1,nc
      do 600 j=i+1,nc
      alpha(j,i)=alpha(i,j)
      beta(j,i)=beta(i,j)
 600  continue
      write(iw,2999)
      i=1
      j=noc(1)
      nn=nc-1
      go to (200,300),nn
 200  do 250 k=1,2
      write(iw,3002)oc(k),i,j,(alpha(k,l),l=1,2),(beta(k,l),l=1,2)
      i=j+1
 250  j=na
      go to 500
 300  do 350 k=1,3
      write(iw,3003)oc(k),i,j,(alpha(k,l),l=1,3),(beta(k,l),l=1,3)
      i=j+1
 350  j=j+noc(k+1)
 500  continue
      return
 1000 format(1x,'the electronic state is not defined')
 1100 format(1x,'half filled open shell is assumed')
 1200 format(1x,'spin multiplicity is not compatible with half filled',
     1      ' open shell. program stops')
 2001 format(1x,'closed shell calculation. nc =',i2,/)
 2002 format(1x,'excited singlet state. nc =',i2)
 2003 format(1x,'doublet state. non degenerate mo. nc =',i3,/)
 2004 format(1x,'half filled configuration. spin multiplicity =',i2,
     1           '. nc =',i2,/)
 2005 format(1x,'two degenerate mos, one electron. doublet state. nc =',
     1           i2,/)
 2006 format(1x,'two degenerate mos, three electrons. doublet state.',
     1           ' nc =',i2,/)
 2007 format(1x,'three degenerate mos, one electron. doublet state.',
     1           ' nc =',i2,/)
 2008 format(1x,'one electron in a non degenerate mo, one in two',
     1       ' degenerate mos of higher energy. singlet state. nc =',i2
     2       ,/)
 2009 format(1x,'one electron in a non degenerate mo, one in two',
     1       ' degenerate mos of higher energy. triplet state. nc =',i2
     2       ,/)
 2010 format(1x,'three electrons in two degenerate mos, one electron',
     1       ' in a non degenerate mo of higher energy. singlet state.',
     2            ' nc =',i2,/)
 2011 format(1x,'three electrons in two degenerate mos, one electron',
     1       ' in a non degenerate mo of higher energy. triplet state.',
     2            ' nc =',i2,/)
 2012 format(1x,'three degenerate mos, two electrons. triplet state.',
     1            ' nc =',i2,/)
 2013 format(1x,'one electron in a non degenerate mo, one in three',
     1       ' degenerate mos of higher energy. singlet state. nc =',i2,
     2         /)
 2014 format(1x,'one electron in a non degenerate mo, one in three',
     1       ' degenerate mos of higher energy. triplet state. nc =',i2
     2         ,/)
 2999 format(3x,'oc',7x,'om',20x,'alpha',30x,'beta',/)
 3001 format(1x,f6.4,2x,i3,' -',i3,14x,f12.7,23x,f12.7)
 3002 format(1x,f6.4,2x,i3,' -',i3,6x,2f12.7,11x,2f12.7)
 3003 format(1x,f6.4,2x,i3,' -',i3,2x,3f12.7,3x,3f12.7)
      end
