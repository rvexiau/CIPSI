c-----------------------------------------------------                          
      SUBROUTINE SHELL(A,IND,N)                                                 
      IMPLICIT REAL*8(A-H,O-Z)                                                  
      include 'pshf.prm'
      PARAMETER(ALN2I=1./0.69314748,TINY=1.E-5)                      
      DIMENSION A(NDETZ),IND(NDETZ)                                             
      LOGNB2=INT(ALOG(FLOAT(N))*ALN2I+TINY)                                     
      M=N                                                                       
      DO NN=1,LOGNB2                                                            
         M=M/2                                                                  
         K=N-M                                                                  
         DO J=1,K                                                               
            I=J                                                                 
   3        CONTINUE                                                            
            L=I+M                                                               
            IF(A(L).GT.A(I)) THEN                                               
                 T=A(I)                                                         
                 A(I)=A(L)                                                      
                 A(L)=T                                                         
                 INDI=IND(I)                                                    
                 IND(I)=IND(L)                                                  
                 IND(L)=INDI                                                    
                 I=I-M                                                          
                 IF(I.GE.1)GO TO 3                                              
            END IF                                                              
         END DO                                                                 
      END DO                                                                    
      RETURN                                                                    
      END                                                                       
