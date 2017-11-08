subroutine wrtempt(nelac,ii,it,jt,jtbuf,ni,nj,ncf,numtab,isymat)
  implicit none
  integer,parameter :: smax=6  
 
  integer,dimension(ncf,smax),intent(inout) :: numtab  
  integer,intent(inout) :: it(smax),jt(smax),jtbuf(3)
  integer,intent(in) :: nelac,ni,nj,ii,isymat(smax),ncf
  integer :: ss,combi
  
  select case(nelac)
        
!     ******************************** 3e- ****************************************    
  case(3)
    do ss=2,4,2
      if(isymat(ss).eq.0) cycle      

      combi=ss+10*nj+100*ni
      select case(combi)
      
! *****  NELAC = 3 ; SS = 2 ; Ni=1 ; Nj=1 ******************************************      
      case(112)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
       
        
! *****  NELAC = 3 ; SS = 2 ; Ni=1 ; Nj=3 ******************************************            
      case(132)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+2

! *****  NELAC = 3 ; SS = 2 ; Ni=3 ; Nj=1 ******************************************    
      case(312)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
!       second determinant on line : store in buffer        
        if(jtbuf(1).eq.0) numtab(it(2)+1,2)=ii
        jtbuf(1)=jtbuf(1)+1
          
! *****  NELAC = 3 ; SS = 2 ; Ni=3 ; Nj=3 ******************************************    
      case(332)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+2
        
!       second determinant on line : store in buffer        
        if(jtbuf(1).eq.0) numtab(it(2)+1,2)=ii
        jtbuf(1)=jtbuf(1)+2
        
! *****  NELAC = 3 ; SS = 4 ; Ni=3 ; Nj=3 ******************************************         
      case(334)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1

      case(114,134)
!        no quartet possible with this determinant      
        cycle
! *****  NELAC = 3 ; end select on SS;Ni;Nj ******************************************        
      case default
        cycle
      end select               
    enddo
    
!     ******************************** 4e- ****************************************    
  case(4)
    do ss=1,5,2
      if(isymat(ss).eq.0) cycle      

      combi=ss+10*nj+100*ni
      select case(combi)
      
      case(111)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
      case(121)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1 
        
      case(211)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1 
        
      case(221)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1  
        
      case(161)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+2

      case(611)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        !       second determinant on line : store in buffer        
        if(jtbuf(1).eq.0) numtab(it(ss)+1,ss)=ii
        jtbuf(1)=jtbuf(1)+1

      case(661)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+2
        
        !       second determinant on line : store in buffer        
        if(jtbuf(1).eq.0) numtab(it(ss)+1,ss)=ii
        jtbuf(1)=jtbuf(1)+2
        
      case(621)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        !       second determinant on line : store in buffer        
        if(jtbuf(1).eq.0) numtab(it(ss)+1,ss)=ii
        jtbuf(1)=jtbuf(1)+1
        
      case(261)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+2
        
      case(663)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+3
        
        !  determinant in buffer        
        if(jtbuf(2).eq.0) numtab(it(ss)+1,ss)=ii
        jtbuf(2)=jtbuf(2)+3
        if(jtbuf(3).eq.0) numtab(it(ss)+2,ss)=ii
        jtbuf(3)=jtbuf(3)+3        
        
      case(623)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        !       second determinant on line : store in buffer        
        if(jtbuf(2).eq.0) numtab(it(ss)+1,ss)=ii
        jtbuf(2)=jtbuf(2)+1
        if(jtbuf(3).eq.0) numtab(it(ss)+2,ss)=ii
        jtbuf(3)=jtbuf(3)+3          
        
      case(263)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+3       
        
      case(665)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1                  
       
! *****  NELAC = 4 ; end select on SS;Ni;Nj ******************************************        
      case default
        cycle
      end select               
    enddo     
  end select
  
  return
     
end subroutine