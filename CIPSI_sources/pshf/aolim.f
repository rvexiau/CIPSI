      subroutine aolim
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc)
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell,kmaz(ds)
      common/atlim/limlow(dc),limsup(dc)
      limlow(1)=1
      lat=1
      j=1
      do 10 i=1,nshell
      iat=katom(i)
      if(lat.eq.iat) go to 10
      lat=iat
      limsup(j)=kloc(i)-1
      j=j+1
      limlow(j)=kloc(i)
   10 continue
      limsup(j)=num
      return
      end
