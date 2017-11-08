      subroutine dvint                                                  
c                                                                       
c     ----- gauss-hermite quadrature using minimum point formula -----  
c                                                                       
      implicit  double precision  (a-h,o-z)                             
      dimension min(6),max(6)                                           
      common/dstv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj     
     1 ,cx,cy,cz                                                        
      common/hermip/ h(21)                                              
      common/wermip/ w(21)                                              
      common /inout/ lec,imp                                            
c                                                                       
      data min /1,2,4,7,11,16/                                          
      data max /1,3,6,10,15,21/                                         
      data zero /0.0d+00/                                               
      xint=zero                                                         
      yint=zero                                                         
      zint=zero                                                         
      npts=(ni+nj+1)/2+1                                                
      imin=min(npts)                                                    
      imax=max(npts)                                                    
      do 13 i=imin,imax                                                 
      dum=h(i)/t                                                        
      ptx=dum+x0                                                        
      pty=dum+y0                                                        
      ptz=dum+z0                                                        
      px=ptx-cx                                                         
      py=pty-cy                                                         
      pz=ptz-cz                                                         
      ax=ptx-xi                                                         
      ay=pty-yi                                                         
      az=ptz-zi                                                         
      bx=ptx-xj                                                         
      by=pty-yj                                                         
      bz=ptz-zj                                                         
      go to (5,4,3,2,1),ni                                              
    1 px=px*ax                                                          
      py=py*ay                                                          
      pz=pz*az                                                          
    2 px=px*ax                                                          
      py=py*ay                                                          
      pz=pz*az                                                          
    3 px=px*ax                                                          
      py=py*ay                                                          
      pz=pz*az                                                          
    4 px=px*ax                                                          
      py=py*ay                                                          
      pz=pz*az                                                          
    5 go to (12,11,10,9,8,7,6),nj                                       
    6 px=px*bx                                                          
      py=py*by                                                          
      pz=pz*bz                                                          
    7 px=px*bx                                                          
      py=py*by                                                          
      pz=pz*bz                                                          
    8 px=px*bx                                                          
      py=py*by                                                          
      pz=pz*bz                                                          
    9 px=px*bx                                                          
      py=py*by                                                          
      pz=pz*bz                                                          
   10 px=px*bx                                                          
      py=py*by                                                          
      pz=pz*bz                                                          
   11 px=px*bx                                                          
      py=py*by                                                          
      pz=pz*bz                                                          
   12 dum=w(i)                                                          
      xint=xint+dum*px                                                  
      yint=yint+dum*py                                                  
      zint=zint+dum*pz                                                  
   13 continue                                                          
      return                                                            
      end                                                               
