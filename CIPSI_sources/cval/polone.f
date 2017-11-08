      subroutine polone
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical out,ymono
      parameter (EPSI=1.D-7)
      common/infnp/nprint
      common/iofile/ir,iw,ip,is,iq,ih
      common/polar/calfa(dc)
      common/infoa/nat,ich,mul,num,nnp,ne,na,nb,zan(dc),c(3,dc),ymono
      common/eigen2/f(doas),h(doas),s(doas)
      common/elchp/elec(3)
Ccc   RG : modified on 01/10/08 : Add operator alpha/R^4 in the MO basis
      common/polxx/polz
      call reada(h,nnp,3)
cCc ===========================
cCc = Read the overlap matrix =
cCc ===========================
      call reada(s,nnp,2) 
      rewind 17
      out=nprint.eq.1
      do 1000 ic=1,nat
      alfad=calfa(ic)
      if(alfad.le.EPSI) go to 1000
      do 30 k=1,3
      el=elec(k)
      read(17) (f(i),i=1,nnp)
      do 10 i=1,nat
      if(i.eq.ic) go to 10
      ri=(c(1,ic)-c(1,i))**2+(c(2,ic)-c(2,i))**2+(c(3,ic)-c(3,i))**2
      ri=ri**(-1.5d0)
      el=el+(c(k,ic)-c(k,i))*zan(i)*ri
   10 continue
      el=-el*alfad
      do 20 i=1,nnp
c      print *,'i, k, f(i),el ',i,k,f(i),-el/alfad
   20 h(i)=h(i)+f(i)*el
   30 continue
      write(6,*) 'matrice h modifiee par x,y,z'
      if(out)then
      jf=0
      do i=1,num
        ji=jf+1
        jf=jf+i
        write(iw,9999) i,(h(j),j=ji,jf)
      end do
      end if
      read(17) (f(i),i=1,nnp)
      do 40 i=1,nnp
c      print *,'i, f(i) ',i,f(i) 
   40 h(i)=h(i)-f(i)*alfad*0.5d0
      write(6,*) 'matrice h modifiee par 1/r**4'
cCc ======================
cCc = Here, I add 1/R**4 =
cCc ======================
cCc       if(polz/=0.d0) write(*,*) 'I add the 1/R**4 operator'
cCc       do i=1,nnp
cCc         h(i)=h(i)+0.25d0*polz*s(i)
cCc       enddo
      if (out) then
      jf =0
      do i=1,num
        ji=jf+1
        jf=jf+i
        write(iw,9999) i,(h(j),j=ji,jf)
      end do
      end if
c
 1000 continue
      call wrtda(h,nnp,3)
      if(out)then
      jf=0
      do 50 i=1,num
      ji=jf+1
      jf=jf+i
   50 write(iw,9999) i,(h(j),j=ji,jf)
 9999 format(i3,(12f10.6))
      end if
      return
      end
