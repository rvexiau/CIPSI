subroutine newdet(nelac,ncf,ordertab,alpha,deg,ygood,ydouble,yquad)
  implicit none
!     assign order rank depending on the spin 
!      2 determinants : +- ; -+
!      3 determinants : +-- ; -+- ; --+
!      6 determinants : ++-- ; +-+- ; +--+ ; -++- ; -+-+ ; --++  

!    add any missing determinants
  integer :: nelac,ncf
  logical,dimension(*) :: ygood,ydouble,yquad
  integer :: alpha(ncf,2),deg(ncf,10+1),ordertab(ncf,10)
  
  integer :: i,j,iwarn,na1,na2,na3,na4,na5,na6,nb1,nb2,nb3,nb4,nb5,nb6

  iwarn=0
  ordertab=0      
  
  select case(nelac)
 
! ***** 2 electron  ******************************************        
  case(2)
    do i=1,ncf
      ygood(i)=.true.
      select case(deg(i,1))
      case(0)
        ygood(i)=.false.    
      case(1)
        if(ydouble(i)) then
          cycle
        else
          ygood(i)=.false.
          iwarn=iwarn+1
        endif      
      case(2)
        if(alpha(deg(i,2),1).lt.alpha(deg(i,3),1)) then
          ordertab(i,1)=1
          ordertab(i,2)=2
        else  
          ordertab(i,1)=2
          ordertab(i,2)=1  
        endif  
      end select    
    enddo

! ***** 3 electron  ******************************************   
  case(3)
    do i=1,ncf
      ygood(i)=.true.
      select case(deg(i,1))
      ! ***** Number of combinaison : 0 *****************************************   
      case(0)
        ygood(i)=.false.
      !***** Number of combinaison : 1 ******************************************            
      case(1)
        if(ydouble(i)) then
          cycle
        else
          ygood(i)=.false.
          iwarn=iwarn+1
          cycle
        endif  
     ! ***** Number of combinaison : 2 ******************************************            
      case(2)        
        ygood(i)=.false.
        iwarn=iwarn+2
        cycle
        
     ! ***** Number of combinaison : 3 ******************************************            
      case(3)
        na1=alpha(deg(i,2),1)
        na2=alpha(deg(i,3),1)
        na3=alpha(deg(i,4),1)
      
        if(na1.lt.na2) then
          if(na1.lt.na3) then        
            ordertab(i,1)=1
            if(na2.lt.na3) then 
              ordertab(i,2)=2
              ordertab(i,3)=3
            else
              ordertab(i,2)=3
              ordertab(i,3)=2
            endif
          else
            ordertab(i,1)=2
            ordertab(i,2)=3
            ordertab(i,3)=1
          endif  
        else
          if(na1.lt.na3) then
            ordertab(i,1)=2
            ordertab(i,2)=1
            ordertab(i,3)=3
          else
            if(na2.lt.na3) then
              ordertab(i,1)=3
              ordertab(i,2)=1
              ordertab(i,3)=2
            else
              ordertab(i,1)=3
              ordertab(i,2)=2
              ordertab(i,3)=1          
            endif
          endif
        endif  
      
      case DEFAULT
        iwarn=iwarn+1
        write(*,*) 'warning, duplicate elements '
      end select
    ! end of loop on ncf    
    enddo

! ***** 4 electron  ******************************************     
  case(4)
    do i=1,ncf
      ygood(i)=.true.
      select case(deg(i,1))
      ! ***** Number of combinaison : 0 *****************************************   
      case(0)
        ygood(i)=.false.
      !***** Number of combinaison : 1 ******************************************            
      case(1)
        if(yquad(i)) then
          cycle
        else  
          ygood(i)=.false.
          iwarn=iwarn+1
        endif  
     ! ***** Number of combinaison : 2 ******************************************            
      case(2)        
        na1=alpha(deg(i,2),1)
        na2=alpha(deg(i,3),1)
        if(na1.lt.na2) then
          ordertab(i,1)=1
          ordertab(i,2)=2
        else  
          ordertab(i,1)=2
          ordertab(i,2)=1  
        endif  
     ! ***** Number of combinaison : 3 ******************************************            
      case(3)
        ygood(i)=.false.
        iwarn=iwarn+3
      case(4)
        ygood(i)=.false.
        iwarn=iwarn+4
      case(5)
        ygood(i)=.false.
        iwarn=iwarn+5
      case(6)  

      
      case DEFAULT
        iwarn=iwarn+1
        write(*,*) 'warning, duplicate elements '
      end select
    ! end of loop on ncf    
    enddo
  end select
  
  write(*,*) 'number of discarded elements ',iwarn
  
  return
end subroutine