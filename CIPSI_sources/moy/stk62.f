      subroutine stk62(ncf)                                   
c     version pour plus de 128 orbitales
c     stockage des indices en entier 2
c     version jpd juin 93
c                                                                       
c                                                                       
c        subroutine che compatta i determinanti e li scrive su tape 60  
c                                                                       
c                                                                       
      implicit real*8(a-h,o-x,z),logical*1(y)                           
      include 'pshf.prm'
      integer*4 nexst,trou,part,nd
      common/det/m,n,idet(16*lungo),nex(lungo),v(lungo) 
      common/dinde/trou(8*ndimh),part(8*ndimh),nexst(ndimh),
     * nd(ndimh)
      m=0
      n=0
      do icf=1,ncf
      nec=nexst(icf)
      if (nec.gt.8) then
      write(6,*) ' degre d''excitation superieure a 8',
     *' programme stoppe'
      end if
      if(nec.eq.0) go to 10                                             
      do 2 i=1,nec                                                      
      m=m+1
      idet(m)=trou(nd(icf)+i)
      m=m+1
      idet(m)=part(nd(icf)+i)
 2    continue                                                          
 10   continue                                                          
      n=n+1                                                             
      nex(n)=nec                                                        
      v(n)=1.d0
      if(n.eq.lungo) then                                             
      write(62)m,n,(idet(i),i=1,m),(nex(i),i=1,n),(v(i),i=1,n)          
      m=0                                                               
      n=0                                                               
      end if                                                            
      end do
      if(n.eq.0) return                                                 
      write(62) m,n,(idet(i),i=1,m),(nex(i),i=1,n),(v(i),i=1,n)         
      return                                                            
      end                                                               
