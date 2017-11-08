      subroutine lecnew(ndbl,nts,nps,nexs,nds,yret)
      implicit real*8(a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      integer*4 nts,nps,nexs
      common  /cipsi/ norb,mnorb,noca,nocb,
     *igela,igelb,ncf,ncper,nsym,isym,isz,ntrsy,metat,
     *yprt,yion
      common/e0cip/ebmp(metz),eben(metz),tau
      dimension idot(16*lungo),necst(lungo),vs(lungo)
      dimension nts(*),nps(*),nexs(*),nds(*)
c
      ndbl=0
      nds(1)=0
      read(60,end=3487,err=9000) mw,nw,(idot(i),i=1,mw),(necst(i),i=1,
     & nw),(vs(i),i=1,nw)
c       write(6,*)'lec', mw,nw,(idot(i),i=1,mw),(necst(i),i=1,
c     & nw),(vs(i),i=1,nw)
      myy=0
      nyy=0
      do  lstf=1,nw
        nec=necst(lstf)
        if(nec.gt.8) then
          write(6,*) 'erreur degre d excitation',nec,' depasse les ',
     *    'possibilites du programme'
          stop
        end if
        if(vs(lstf).lt.tau) then
	  myy=myy+2*nec
	  else
          ndbl=ndbl+1
          nexs(ndbl)=nec
          if(nec.ne.0)then
	    do  kstf=1,nec
	      myy=myy+1
              nts(kstf+nyy)=idot(myy)
	      myy=myy+1
              nps(kstf+nyy)=idot(myy)
            end do
	    nyy=nyy+nec
          end if
        end if
      end do
      return
 3487 continue
      yret=.true.
      return
 9000 write(6,9010)
 9010 format(//,' ************ erreur sur la file 60 **************
     y  programme stoppe')
      stop
      end
