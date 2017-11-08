  real*8 function hntd2(nelac,orbtab1,spintab1,orbtab2,spintab2)
    implicit none
    integer :: spintab(nelac),orbtab(nelac)
    integer :: i,j
    
    real*8,external :: ai
    
    do i=1,nelac
      do j=1,nelac
        n1=orbtab1(i)
        n2=orbtab2(j)
        if(n2.lt.n1) then
          cycle
        elseif(n2.eq.n1) then
          if(spintab(i).eq.spintab(j)) exit
          cycle
        else
          
        endif

      enddo
    enddo
    
    hntd2=0d0
  
  end function