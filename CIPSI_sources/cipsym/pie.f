      subroutine pie(nec,nts,nps,xam,venvp,metat)                                   
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
      common/det/my,ny,idet(16*lungo),nexst(lungo),vst(lungo)
     * ,evp(lungo*metz),ntot
      dimension nts(8),nps(8),venvp(metz)                 
      if (nec.gt.8) then
      write(6,*) ' degre d''excitation superieure a 8',
     *' programme stoppe',nec
      end if
      if(nec.eq.0) go to 10                                             
      do 2 i=1,nec                                                      
      my=my+1
      idet(my)=nts(i)
      my=my+1
      idet(my)=nps(i)
 2    continue                                                          
 10   continue                                                          
      ny=ny+1            
      ntot=ntot+1
      nexst(ny)=nec                                                        
      vst(ny)=xam                
      do i=1,metat
        evp((ny-1)*metat+i)=venvp(i)
      enddo
      if(ny.lt.lungo) return                                             
      write(60)my,ny,(idet(i),i=1,my),(nexst(i),i=1,ny),(vst(i),i=1,ny)
     * ,(evp(i),i=1,ny*metat)          
c     write(6,*)'lec60',m,n,(idet(i),i=1,m),(nex(i),i=1,n),(v(i),i=1,n)          
      my=0                                                               
      ny=0                                                               
      return                                                            
      entry fine(metat)                                                        
      if(ny.eq.0) return                                                 
      write(60) my,ny,(idet(i),i=1,my),(nexst(i),i=1,ny),(vst(i),i=1,ny)        
     * ,(evp(i),i=1,ny*metat) 
c     write(6,*) 'lec60',m,n,(idet(i),i=1,m),(nex(i),i=1,n),(v(i),i=1,n)         
      my=0                                                               
      ny=0                                                               
      return                                                            
      end                                                               
