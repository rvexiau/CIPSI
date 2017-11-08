      subroutine isoout(nt)                                             
      implicit real*8 (a-h,o-z)                                         
      integer shiftr                                                    
      common/machin/isingl,nbits                                        
      common/isopac/indin(48),indout(12)                                
c                                                                       
      mask1(iarg)=2**iarg-1                                             
cray
c     shiftr(iarg,ibit)= rshft(iarg, ibit)                              
c     land(iarg1,iarg2)= and(iarg1,iarg2)                               
cray
chp9000
      shiftr(iarg,ibit)= ishft(iarg,-ibit)                              
      land(iarg1,iarg2)= iand(iarg1,iarg2)                              
chp9000
      isopck(iarg,ibit)=land(shiftr(iarg,ibit),mask1(8))                
      ibit1=8                                                           
c                                                                       
      do 10 it=1,nt                                                     
   10 indin(it)=0                                                       
      it  =0                                                            
      nout=0                                                            
      max =0                                                            
   20 min =max+1                                                        
      minp=min                                                          
      max =max+(8/isingl)                                               
      maxp=max                                                          
      if(max.gt.nt) max=nt                                              
      nout=nout+1                                                       
      ipack=indout(nout)                                                
      do 30 i=minp,maxp                                                 
      ibit=ibit1*(maxp-i)                                               
      idum=isopck(ipack,ibit)                                           
      if(idum.eq.0) go to 30                                            
      it=it+1                                                           
      indin(it)=idum                                                    
   30 continue                                                          
c                                                                       
      if(max.lt.nt) go to 20                                            
      return                                                            
      end                                                               
