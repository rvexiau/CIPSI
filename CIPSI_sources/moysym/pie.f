      subroutine pie(nec,nts,nps,ndet)
      implicit real*8(a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      integer*4 trou,part,nexst,nd
      common/dinde/trou(8*ndimh),part(8*ndimh),nexst(ndimh),
     * nd(ndimh)
      dimension nts(*),nps(*)
      if(nec.ne.0)then 
	ndd=nd(ndet)
	do i=1,nec
	  trou(ndd+i)=nts(i)
	  part(ndd+i)=nps(i)
	end do
      end if
      nexst(ndet)=nec
      nd(ndet+1)=nd(ndet)+nec
      return
      end
