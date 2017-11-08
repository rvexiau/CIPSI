      subroutine comb(lkt,ljt)                                          
      implicit  double precision  (a-h,o-z)                             
      dimension md1(4),md2(4)                                           
      dimension t(100)                                                  
      common/recomb/s(225,20),itpmax                                    
      common /inout/lec,imp                                             
      data c1,c2/0.8660254037844d0,0.5d0/                               
      data md1/1,3,6,10/,md2/1,3,5,7/                                   
      kd1=md1(lkt)                                                      
      kd2=md2(lkt)                                                      
      jd1=md1(ljt)                                                      
      jd2=md2(ljt)                                                      
      jk1=jd1*kd1                                                       
      jk2=jd1*kd2                                                       
c                                                                       
      do 9000 ig=1,itpmax                                               
      do 10 i=1,jk1                                                     
   10 t(i)=s(i,ig)                                                      
      if(lkt-3) 100,20,40                                               
   20 continue                                                          
      do 30 j1=1,jd1                                                    
      t1=t(j1)                                                          
      t2=t(j1+jd1)                                                      
      t3=t(j1+jd1+jd1)                                                  
      s(j1,ig)=c1*(t1-t2)                                               
      s(j1+jd1,ig)=t3-c2*(t1+t2)                                        
      do 25 l=3,5                                                       
   25 s(j1+(l-1)*jd1,ig)=t(j1+l*jd1)                                    
   30 continue                                                          
      go to 100                                                         
   40 write(imp,1000)                                                   
 1000 format(/' erreur pas de sym f')                                   
      stop                                                              
  100 continue                                                          
      if(ljt.lt.3) go to 9000                                           
      do 110 i=1,jk2                                                    
  110 t(i)=s(i,ig)                                                      
      do 130 k=1,kd2                                                    
      k1=(k-1)*jd1+1                                                    
      k2=(k-1)*jd2+1                                                    
      t1=t(k1)                                                          
      t2=t(k1+1)                                                        
      t3=t(k1+2)                                                        
      s(k2,ig)=c1*(t1-t2)                                               
      s(k2+1,ig)=t3-c2*(t1+t2)                                          
      do 125 l=3,5                                                      
  125 s(k2+l-1,ig)=t(k1+l)                                              
  130 continue                                                          
 9000 continue                                                          
      return                                                            
      end                                                               
