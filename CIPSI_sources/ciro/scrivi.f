      subroutine scrivi(ind,norb)                                       
      implicit real*8(a-h,o-x,z),logical*1(y)                           
      include 'pshf.prm'
      common/cd/r(ncouz*(doa**2)),ndd,ndt,ncf1,ncf2,metat,icf2,ncou,  
     *ibeg(ncouz),ietat(ncouz),jetat(ncouz),ycou(metz,metz),      
     &yiet(metz),yjet(metz)                                             
9969  format(i5,12f10.6)                                                
9971  format(/,5x,12(4x,i3,3x))                                         
      imax=0                                                            
      max=12                                                            
1601  imin=imax+1                                                       
      imax=imax+max                                                     
      if(imax.gt.norb) imax=norb                                        
      write(6,9971)(i,i=imin,imax)                                      
      do 1602 i=1,norb                                                  
1602  write(6,9969) i,(r(ind+(i-1)*norb+j),j=imin,imax)                 
      if(imax.eq.norb) go to 1603                                       
      go to 1601                                                        
1603  return                                                            
      end                                                               
