      subroutine pie(nec,nts,nps,xam)                                   
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
      common/det/m,n,idet(16*lungo),nex(lungo),v(lungo)
      dimension nts(8),nps(8)                            
      if (nec.gt.8) then
      write(6,*) ' degre d''excitation superieure a 8',
     *' programme stoppe',nec
      end if
      if(nec.eq.0) go to 10                                             
      do 2 i=1,nec                                                      
      m=m+1
      idet(m)=nts(i)
      m=m+1
      idet(m)=nps(i)
 2    continue                                                          
 10   continue                                                          
      n=n+1                                                             
      nex(n)=nec                                                        
      v(n)=xam                                                          
      if(n.lt.lungo) return                                             
      write(60)m,n,(idet(i),i=1,m),(nex(i),i=1,n),(v(i),i=1,n)          
c     write(6,*)'lec60',m,n,(idet(i),i=1,m),(nex(i),i=1,n),(v(i),i=1,n)          
      m=0                                                               
      n=0                                                               
      return                                                            
      entry fine                                                        
      if(n.eq.0) return                                                 
      write(60) m,n,(idet(i),i=1,m),(nex(i),i=1,n),(v(i),i=1,n)         
c     write(6,*) 'lec60',m,n,(idet(i),i=1,m),(nex(i),i=1,n),(v(i),i=1,n)         
      m=0                                                               
      n=0                                                               
      return                                                            
      end                                                               
