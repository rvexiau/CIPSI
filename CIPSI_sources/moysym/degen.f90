subroutine degen(nelac,trou,part,nd,nexst,ygood,deg,ordertab,ncf,norb)
  implicit none
  integer,parameter :: smax=6
!     order for the spin eigenfucntions :
!      2 determinants : +- ; -+
!      3 determinants : +-- ; -+- ; --+
!      6 determinants : ++-- ; +-+- ; +--+ ; -++- ; -+-+ ; --++ 

!    order for the alpha array :
!    sort "part" first then if tied sort "trou"
  
  integer :: nelac,ncf,norb
  integer,dimension(*) :: trou,part,nexst,nd
  integer :: deg(ncf,11),ordertab(ncf,10)
  logical :: ygood(ncf)
   
  integer :: i,j,k,ne1,np,np2,np3,np4,nt2,iadd
  integer,dimension(:),allocatable :: alpha
  integer,dimension(:,:),allocatable :: orbtab
  logical :: ydeg
  logical,dimension(:),allocatable :: ydouble,yquad
  
  allocate(orbtab(ncf,nelac),alpha(ncf),ydouble(ncf),yquad(ncf))
  
  ydouble=.false.
  yquad=.false.
  orbtab=-1
  alpha=1
  deg=0
   
  select case(nelac)
! ************************* 2 electrons **********************************  
  case(2)    
    do i=1,ncf 
      ydeg=.false.
      ne1=nexst(i)
 !    find occupied orbitals (store the position of the up e- in alpha) 
      select case(ne1)
      case(0)
        ydouble(i)=.true.
      case(1)
        if(trou(nd(i)+1).eq.1) then
          np=part(nd(i)+1)
          alpha(i)=1
          orbtab(i,1)=1
          orbtab(i,2)=np
        else
          alpha(i)=2
          orbtab(i,1)=1
          orbtab(i,2)=part(nd(i)+1)-norb
        endif
      case(2)
        np=part(nd(i)+1)
        np2=part(nd(i)+2)-norb
        if(np.lt.np2)then
          alpha(i)=1
          orbtab(i,1)=np
          orbtab(i,2)=np2                    
        elseif(np.eq.np2) then
          ydouble(i)=.true.
        else
          alpha(i)=2
          orbtab(i,1)=np2
          orbtab(i,2)=np
        endif  
      end select
      
      do j=i-1,1,-1
        if(deg(j,1).eq.0)cycle
        if(ydouble(j)) cycle
        if(.not.(ne1.eq.nexst(j))) cycle          
        ydeg=.true.
        do k=1,nelac
          if(.not.(orbtab(i,k).eq.orbtab(j,k))) then
            ydeg=.false.
            exit
          endif  
        enddo
        if(ydeg) then
          deg(j,1)=deg(j,1)+1
          deg(i,1)=0
          ygood(i)=.false.
          deg(j,deg(j,1)+1) = i   
          ordertab(j,deg(j,1))=alpha(i)
          exit
        endif 
      enddo   
      if(.not.ydeg) then
        deg(i,1)=1
        deg(i,2)=i
        ordertab(i,1)=alpha(i)
      endif  
    enddo 

! ************************* 3 electrons **********************************  
  case(3)     
    do i=1,ncf 
      ydeg=.false.
      ne1=nexst(i)
 !    find occupied orbitals 
      select case(ne1)
      case(1)
        ydouble(i)=.true.
      case(2)
        if(trou(nd(i)+1).eq.1) then
          nt2=trou(nd(i)+2)
          if(nt2.eq.norb+1) then
            ydouble(i)=.true.
          elseif(nt2.eq.2) then
            np=part(nd(i)+1)
            orbtab(i,1)=1
            orbtab(i,2)=2
            orbtab(i,3)=np
            alpha(i)=3
          else
            orbtab(i,1)=1
            orbtab(i,2)=2
            orbtab(i,3)=part(nd(i)+1)-norb
            alpha(i)=2
          endif
        else
          if(trou(nd(i)+2).eq.norb+2) then
            ydouble(i)=.true.
          else
            orbtab(i,1)=1
            orbtab(i,2)=2
            orbtab(i,3)=part(nd(i)+1)-norb
            alpha(i)=1
          endif
        endif
      case(3)
        if(trou(nd(i)+1).eq.1) then
          nt2=trou(nd(i)+2)
          if(nt2.eq.2) then
            if(trou(nd(i)+3).eq.norb+1) then              
              np=part(nd(i)+1)
              np2=part(nd(i)+2)-norb
              if(np.lt.np2) then
                alpha(i)=1
                orbtab(i,1)=2                
                orbtab(i,2)=np
                orbtab(i,3)=np2
              elseif(np.eq.np2) then
                ydouble(i)=.true.
              else
                alpha(i)=2
                orbtab(i,1)=2                      
                orbtab(i,2)=np2
                orbtab(i,3)=np        
              endif     
            else
              np=part(nd(i)+1)
              np2=part(nd(i)+2)-norb
              if(np.lt.np2) then
                alpha(i)=1
                orbtab(i,1)=1                
                orbtab(i,2)=np
                orbtab(i,3)=np2
              elseif(np.eq.np2) then
                ydouble(i)=.true.
              else
                alpha(i)=2
                orbtab(i,1)=1                      
                orbtab(i,2)=np2            
                orbtab(i,3)=np        
              endif    
            endif
          else
            alpha(i)=3
            orbtab(i,1)=2
            orbtab(i,2)=part(nd(i)+1)-norb
            orbtab(i,3)=part(nd(i)+2)-norb
          endif
        else
          alpha(i)=3
          orbtab(i,1)=1
          orbtab(i,2)=part(nd(i)+1)-norb
          orbtab(i,3)=part(nd(i)+2)-norb          
        endif        
      case(4)
        np=part(nd(i)+1)
        np2=part(nd(i)+2)-norb
        np3=part(nd(i)+3)-norb
        if(np.lt.np2) then
          alpha(i)=1
          orbtab(i,1)=np
          orbtab(i,2)=np2
          orbtab(i,3)=np3
        elseif(np.eq.np2) then
          ydouble(i)=.true.
        elseif(np.lt.np3) then
          alpha(i)=2
          orbtab(i,1)=np2
          orbtab(i,2)=np
          orbtab(i,3)=np3    
        elseif(np.eq.np3) then
          ydouble(i)=.true.
        else
          alpha(i)=3
          orbtab(i,1)=np2
          orbtab(i,2)=np3
          orbtab(i,3)=np   
        endif
      end select
      
      if(ydouble(i)) then
        deg(i,1)=1
        deg(i,2)=i
      else
        do j=i-1,1,-1
          if(deg(j,1).eq.0.or.deg(j,1).ge.3)cycle
          if(ydouble(j)) cycle
          if(.not.(ne1.eq.nexst(j))) cycle          
          ydeg=.true.
          do k=1,nelac
            if(.not.(orbtab(i,k).eq.orbtab(j,k))) then
              ydeg=.false.
              exit
            endif  
          enddo
          if(ydeg) then
            deg(j,1)=deg(j,1)+1
            deg(i,1)=0
            ygood(i)=.false.
            deg(j,deg(j,1)+1) = i   
            ordertab(j,deg(j,1))=alpha(i)
            exit
          endif 
        enddo   
        if(.not.ydeg) then
          deg(i,1)=1
          deg(i,2)=i
          ordertab(i,1)=alpha(i)
        endif
      endif  
    enddo 

! ************************* 4 electrons **********************************      
  case(4)     
    do i=1,ncf 
      ydeg=.false.
      ne1=nexst(i)
 !    find occupied orbitals 
 
      select case(ne1)
      case(1)
        ydouble(i)=.true.
        np=part(nd(i)+1)
        if(trou(nd(i)+1).eq.1) then
          alpha(i)=1
          orbtab(i,1)=1
          orbtab(i,2)=2
          orbtab(i,3)=2
          orbtab(i,4)=np
        elseif(trou(nd(i)+1).eq.2) then
          alpha(i)=1
          orbtab(i,1)=1
          orbtab(i,2)=1
          orbtab(i,3)=2
          orbtab(i,4)=np        
        elseif(trou(nd(i)+1).eq.norb+1) then
          alpha(i)=2
          orbtab(i,1)=1
          orbtab(i,2)=2
          orbtab(i,3)=2
          orbtab(i,4)=np-norb         
        elseif(trou(nd(i)+1).eq.norb+2) then  
          alpha(i)=2
          orbtab(i,1)=1
          orbtab(i,2)=1
          orbtab(i,3)=2
          orbtab(i,4)=np-norb             
        endif
      case(2)
        np=part(nd(i)+1)
        np2=part(nd(i)+2)
        if(trou(nd(i)+1).eq.1) then
          if(trou(nd(i)+2).eq.2) then
            alpha(i)=1
            orbtab(i,1)=1
            orbtab(i,2)=2
            orbtab(i,3)=np
            orbtab(i,4)=np2    
          elseif(trou(nd(i)+2).eq.norb+1) then
            np2=np2-norb
            if(np.eq.np2) then
              yquad(i)=.true.
            elseif(np.lt.np2) then  
              ydouble(i)=.true.
              alpha(i)=1
              orbtab(i,1)=2
              orbtab(i,2)=2
              orbtab(i,3)=np
              orbtab(i,4)=np2
            else  
              ydouble(i)=.true.
              alpha(i)=2
              orbtab(i,1)=2
              orbtab(i,2)=2
              orbtab(i,3)=np2
              orbtab(i,4)=np
            endif  
          else
            np2=np2-norb
            if(np.lt.np2) then
              alpha(i)=2
              orbtab(i,1)=1
              orbtab(i,2)=2
              orbtab(i,3)=np
              orbtab(i,4)=np2   
            elseif(np.eq.np2) then
              ydouble(i)=.true.
              alpha(i)=1
              orbtab(i,1)=1
              orbtab(i,2)=2
              orbtab(i,3)=np
              orbtab(i,4)=np2   
            else
              alpha(i)=4
              orbtab(i,1)=1
              orbtab(i,2)=2
              orbtab(i,3)=np2
              orbtab(i,4)=np   
            endif  
          endif
        elseif(trou(nd(i)+1).eq.2) then
          if(trou(nd(i)+2).eq.norb+1) then
            np2=np2-norb
            if(np.lt.np2) then
              alpha(i)=3
              orbtab(i,1)=1
              orbtab(i,2)=2
              orbtab(i,3)=np
              orbtab(i,4)=np2
            elseif(np.eq.np2) then
              ydouble(i)=.true.
              alpha(i)=2
              orbtab(i,1)=1
              orbtab(i,2)=2
              orbtab(i,3)=np
              orbtab(i,4)=np 
            else
              alpha(i)=5
              orbtab(i,1)=1
              orbtab(i,2)=2
              orbtab(i,3)=np2
              orbtab(i,4)=np 
            endif
          else
            np2=np2-norb
            if(np.eq.np2) then
              yquad(i)=.true.
            elseif(np.lt.np2) then  
              ydouble(i)=.true.
              alpha(i)=1
              orbtab(i,1)=1
              orbtab(i,2)=1
              orbtab(i,3)=np
              orbtab(i,4)=np2
            else 
              ydouble(i)=.true.
              alpha(i)=2
              orbtab(i,1)=1
              orbtab(i,2)=1
              orbtab(i,3)=np2
              orbtab(i,4)=np
            endif  
          endif
        else
         alpha(i)=6
         orbtab(i,1)=1
         orbtab(i,2)=2
         orbtab(i,3)=np-norb
         orbtab(i,4)=np2-norb
        endif      
      case(3)
        np=part(nd(i)+1)
        np2=part(nd(i)+2)
        np3=part(nd(i)+3)-norb
        if(trou(nd(i)+1).eq.1) then
          if(trou(nd(i)+2).eq.2) then
            if(trou(nd(i)+3).eq.norb+1) then
              if(np3.lt.np) then
                alpha(i)=3
                orbtab(i,1)=2
                orbtab(i,2)=np3
                orbtab(i,3)=np
                orbtab(i,4)=np2
              elseif(np3.eq.np) then
                alpha(i)=1
                ydouble(i)=.true.
                orbtab(i,1)=2
                orbtab(i,2)=np3
                orbtab(i,3)=np
                orbtab(i,4)=np2
              elseif(np3.lt.np2) then
                alpha(i)=2
                orbtab(i,1)=2
                orbtab(i,2)=np
                orbtab(i,3)=np3
                orbtab(i,4)=np2
              elseif(np3.eq.np2) then
                alpha(i)=1
                ydouble(i)=.true.
                orbtab(i,1)=2
                orbtab(i,2)=np
                orbtab(i,3)=np3
                orbtab(i,4)=np2
              else
                alpha(i)=1
                orbtab(i,1)=2
                orbtab(i,2)=np
                orbtab(i,3)=np2
                orbtab(i,4)=np3
              endif
            else
              if(np3.lt.np) then
                alpha(i)=3
                orbtab(i,1)=1
                orbtab(i,2)=np3
                orbtab(i,3)=np
                orbtab(i,4)=np2
              elseif(np3.eq.np) then
                alpha(i)=1
                ydouble(i)=.true.
                orbtab(i,1)=1
                orbtab(i,2)=np3
                orbtab(i,3)=np
                orbtab(i,4)=np2
              elseif(np3.lt.np2) then
                alpha(i)=2
                orbtab(i,1)=1
                orbtab(i,2)=np
                orbtab(i,3)=np3
                orbtab(i,4)=np2
              elseif(np3.eq.np2) then
                alpha(i)=1
                ydouble(i)=.true.
                orbtab(i,1)=1
                orbtab(i,2)=np
                orbtab(i,3)=np3
                orbtab(i,4)=np2
              else
                alpha(i)=1
                orbtab(i,1)=1
                orbtab(i,2)=np
                orbtab(i,3)=np2
                orbtab(i,4)=np3
              endif
            endif
          else
            np2=np2-norb
            if(np.lt.np2) then
              alpha(i)=4
              orbtab(i,1)=2
              orbtab(i,2)=np
              orbtab(i,3)=np2
              orbtab(i,4)=np3
            elseif(np.eq.np2) then
              alpha(i)=2
              ydouble(i)=.true.
              orbtab(i,1)=2
              orbtab(i,2)=np
              orbtab(i,3)=np2
              orbtab(i,4)=np3
            elseif(np.lt.np3) then
              alpha(i)=5
              orbtab(i,1)=2
              orbtab(i,2)=np2
              orbtab(i,3)=np
              orbtab(i,4)=np3
            elseif(np.eq.np3) then
              alpha(i)=2
              ydouble(i)=.true.
              orbtab(i,1)=2
              orbtab(i,2)=np2
              orbtab(i,3)=np
              orbtab(i,4)=np3
            else
              alpha(i)=6
              orbtab(i,1)=2
              orbtab(i,2)=np2
              orbtab(i,3)=np3
              orbtab(i,4)=np
            endif
          endif
        else
          np2=np2-norb
          if(np.lt.np2) then
            alpha(i)=4
            orbtab(i,1)=1
            orbtab(i,2)=np
            orbtab(i,3)=np2
            orbtab(i,4)=np3
          elseif(np.eq.np2) then
            alpha(i)=2
            ydouble(i)=.true.
            orbtab(i,1)=1
            orbtab(i,2)=np
            orbtab(i,3)=np2
            orbtab(i,4)=np3
          elseif(np.lt.np3) then
            alpha(i)=5
            orbtab(i,1)=1
            orbtab(i,2)=np2
            orbtab(i,3)=np
            orbtab(i,4)=np3
          elseif(np.eq.np3) then
            alpha(i)=2
            ydouble(i)=.true.
            orbtab(i,1)=1
            orbtab(i,2)=np2
            orbtab(i,3)=np
            orbtab(i,4)=np3
          else
            alpha(i)=6
            orbtab(i,1)=1
            orbtab(i,2)=np2
            orbtab(i,3)=np3
            orbtab(i,4)=np
          endif
        endif      
      case(4)
        np=part(nd(i)+1)
        np2=part(nd(i)+2)
        np3=part(nd(i)+3)-norb
        np4=part(nd(i)+4)-norb
        if(np.lt.np3) then
          if(np2.lt.np3) then
            alpha(i)=1
            orbtab(i,1)=np
            orbtab(i,2)=np2
            orbtab(i,3)=np3
            orbtab(i,4)=np4   
          elseif(np2.eq.np3) then
            alpha(i)=1
            ydouble(i)=.true.
            orbtab(i,1)=np
            orbtab(i,2)=np2
            orbtab(i,3)=np3
            orbtab(i,4)=np4 
          elseif(np2.lt.np4) then
            alpha(i)=2
            orbtab(i,1)=np
            orbtab(i,2)=np3
            orbtab(i,3)=np2
            orbtab(i,4)=np4 
          elseif(np2.eq.np4) then
            alpha(i)=1
            ydouble(i)=.true.
            orbtab(i,1)=np
            orbtab(i,2)=np3
            orbtab(i,3)=np2
            orbtab(i,4)=np4 
          else
            alpha(i)=3
            orbtab(i,1)=np
            orbtab(i,2)=np3
            orbtab(i,3)=np4
            orbtab(i,4)=np2 
          endif
        elseif(np.eq.np3) then
          ydouble(i)=.true.
          if(np2.lt.np4) then
            alpha(i)=1
            orbtab(i,1)=np
            orbtab(i,2)=np3
            orbtab(i,3)=np2
            orbtab(i,4)=np4 
          elseif(np2.eq.np4) then
            yquad(i)=.true.
          else
            alpha(i)=2
            orbtab(i,1)=np
            orbtab(i,2)=np3
            orbtab(i,3)=np4
            orbtab(i,4)=np2
          endif
        elseif(np.lt.np4) then
          if(np2.lt.np4) then
            alpha(i)=4
            orbtab(i,1)=np3
            orbtab(i,2)=np
            orbtab(i,3)=np2
            orbtab(i,4)=np4 
          elseif(np2.eq.np4) then
            alpha(i)=2
            ydouble(i)=.true.
            orbtab(i,1)=np3
            orbtab(i,2)=np
            orbtab(i,3)=np2
            orbtab(i,4)=np4
          else
            alpha(i)=5
            orbtab(i,1)=np3
            orbtab(i,2)=np
            orbtab(i,3)=np4
            orbtab(i,4)=np2 
          endif
        elseif(np.eq.np4) then
          alpha(i)=2
          ydouble(i)=.true.
          orbtab(i,1)=np3
          orbtab(i,2)=np
          orbtab(i,3)=np4
          orbtab(i,4)=np2
        else
          alpha(i)=6
          orbtab(i,1)=np3
          orbtab(i,2)=np4
          orbtab(i,3)=np
          orbtab(i,4)=np2 
        endif
      end select
     
       ! associate equivalent determinant
     
      if(yquad(i)) then
        deg(i,1)=1
        deg(i,2)=i
      else
        do j=i-1,1,-1
          if(yquad(j)) cycle
          if(deg(j,1).eq.0.or.deg(j,1).ge.6) cycle
          if(ydouble(j).and.(.not.ydouble(i)).or.ydouble(i).and.(.not.ydouble(j))) cycle
          if(ydouble(j).and.deg(j,1).ge.2) cycle
          if(.not.(ne1.eq.nexst(j))) cycle          
          ydeg=.true.
          do k=1,nelac
            if(.not.(orbtab(i,k).eq.orbtab(j,k))) then
              ydeg=.false.
              exit
            endif  
          enddo
          if(ydeg) then
            deg(j,1)=deg(j,1)+1
            deg(i,1)=0
            ygood(i)=.false.
            deg(j,deg(j,1)+1) = i
            ordertab(j,deg(j,1))=alpha(i)
            exit
          endif 
        enddo   
        if(.not.ydeg) then
          deg(i,1)=1
          deg(i,2)=i
          ordertab(i,1)=alpha(i)
        endif
      endif  
      
    enddo  
  
  end select
  
  deallocate(orbtab,ydouble,yquad,alpha)
  
  return
end  