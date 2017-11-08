      subroutine isoin(nt)                                              
      implicit real*8 (a-h,o-z)                                         
      integer shiftl                                                    
      common/machin/isingl,nbits                                        
      common/isopac/indin(48),indout(12)                                
c                                                                       
cray
c     shiftl(iarg,ibit)= lshft(iarg, ibit)                              
c     lor(iarg1,iarg2) = or(iarg1,iarg2)                                
cray
chp9000
      shiftl(iarg,ibit)=ishft(iarg, ibit)                               
      lor(iarg1,iarg2) = ior(iarg1,iarg2)                               
chp9000
      isopck(iarg1,iarg2)=lor(shiftl(iarg1,8),shiftl(iarg2,0))          
c                                                                       
      nout=0                                                            
      max =0                                                            
   10 min =max+1                                                        
      max =max+(8/isingl)                                               
      if(max.gt.nt) max=nt                                              
      ipack=0                                                           
      do 20 it=min,max                                                  
      ipack=isopck(ipack,indin(it))                                     
   20 continue                                                          
      nout=nout+1                                                       
      indout(nout)=ipack                                                
      if(max.lt.nt) go to 10                                            
      return                                                            
      end                                                               
