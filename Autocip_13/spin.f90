subroutine sym_vec(netat,f_spin,nelac,ysort,group,isymat,is,i_fail)
  ! determination of the spin of each state by looking at *.spin file
  ! check with the symmetry of eigenvectors, by looking at the relative sign of
  !  the first pair of components which have the same absolute value
  !  is= 2S+1
  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: i,is(netat),netat,i_fail(netat),istatus,iline,l_str,m,ii,isb,nelac
  integer :: isymat  
  logical :: ysort,ygood
  real(dp) :: dnormes,dnormep,dnormed,dnormef
  character(len=180) :: f_spin
  character(len=*) :: group
  
  i_fail=0
  is=1
  if(.not.ysort) return ! no sorting done
  open(12,file=f_spin)
  
  call find_str( 12, ' metat     2S+1           e_var', istatus, iline, l_str )
    if(istatus==0) then
      write( 6, * )  'Missing meta card in SPIN file'
      i_fail=1
      return
    endif
  do i=1,netat
    read(12,*) m,is(i)
  ! for 2 electrons is=1/-1 for singlet/triplet, for 3 electrons is=1/-1 for doublet/quartet
    select case(is(i))
      case (1)
        select case (nelac)
          case (1)
            i_fail(i)=1
          case (2)
            i_fail(i)=0
          case (3)
            i_fail(i)=1
          case(4)
            i_fail(i)=0
          case default
            i_fail(i)=1
            write(6,*) 'Warning spin with more than 4 nelac not implemented'
        end select    
      case (2)
        select case (nelac)
          case (1)
            i_fail(i)=0
          case (2)
            i_fail(i)=1
          case (3)
            i_fail(i)=0
          case(4)
            i_fail(i)=1            
          case default
            i_fail(i)=1
            write(6,*) 'Warning spin with more than 4 nelac not implemented'
        end select    
      case (3)
        select case (nelac)
          case (1)
            i_fail(i)=1
          case (2)
            i_fail(i)=0
          case (3)
            i_fail(i)=1
          case(4)
            i_fail(i)=0            
          case default
            i_fail(i)=1
            write(6,*) 'Warning spin with more than 4 nelac not implemented'
        end select    
      case (4)
        select case (nelac)
          case (1)
            i_fail(i)=1
          case (2)
            i_fail(i)=1
          case (3)
            i_fail(i)=0
          case(4)
            i_fail(i)=1            
          case default
            i_fail(i)=1
            write(6,*) 'Warning spin with more than 4 nelac not implemented'
        end select    
      case (5)
        select case (nelac)
          case (1)
            i_fail(i)=1
          case (2)
            i_fail(i)=1
          case (3)
            i_fail(i)=1
          case(4)
            i_fail(i)=0            
          case default
            i_fail(i)=1
            write(6,*) 'Warning spin with more than 4 nelac not implemented'
        end select          
      case default
    end select  
  enddo
  
  if(group.eq.'CINFV'.or.group.eq.'cinfv'.or.group.eq.'DINFH'.or.group.eq.'dinfh') then ! sort parasite states
    rewind(12)
    call find_str( 12, 'Contribution of monoexcited determinant', istatus, iline, l_str )
    if(istatus==0) then
      write( 6, * )  'Missing norme card in SPIN file'
      i_fail=1
      return
    endif
    read(12,*)
    do i=1,netat
      read(12,*) ii,ygood,dnormes,dnormep,dnormed,dnormef
      if(.not.ygood) i_fail(i)=1    
    enddo
  endif
  
  close(12)
  return
end subroutine sym_vec
