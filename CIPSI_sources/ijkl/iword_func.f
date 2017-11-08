c     RV 01/16 rewrite iword function in fortran 
      function iword1(iarg)
        implicit none
        integer :: iword1
        integer :: iarg,mask1,ibit,shiftr
        
        ibit=24
        mask1=2**8-1          
        
        shiftr = ishft(iarg,-ibit)                              
        iword1= iand(shiftr,mask1)                              
        return                                         
      end function
      
      function iword2(iarg)
        implicit none
        integer :: iword2
        integer :: iarg,mask1,ibit,shiftr
        
        ibit=16
        mask1=2**8-1          
        
        shiftr = ishft(iarg,-ibit)                              
        iword2= iand(shiftr,mask1)                              
        return                                         
      end function
      
      function iword3(iarg)
        implicit none
        integer :: iword3
        integer :: iarg,mask1,ibit,shiftr
        
        ibit=8
        mask1=2**8-1          
        
        shiftr = ishft(iarg,-ibit)                              
        iword3= iand(shiftr,mask1)                              
        return                                         
      end function
      
      function iword4(iarg)
        implicit none
        integer :: iword4
        integer :: iarg,mask1,ibit,shiftr
        
        ibit=0
        mask1=2**8-1          
        
        shiftr = ishft(iarg,-ibit)                              
        iword4= iand(shiftr,mask1)                              
        return                                         
      end function