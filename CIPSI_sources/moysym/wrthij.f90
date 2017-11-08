subroutine wrthij(nbuf,itab,jtab,vbuf,nbufbuf,itabbuf,jtabbuf,vbufbuf,nhij,nelac,ncf,ii,jj,it,jt,jtbuf,ni,nj,hijsym,ordertab,numtab,isymat)
  implicit none
  integer,parameter :: lsize=8192
  integer,parameter :: smax=6  
  integer,parameter :: ndimh=10000000

! nhij size come from max value of 2S+1  
  integer,dimension(ncf,10),intent(in) :: ordertab
  integer,dimension(ncf,smax),intent(inout) :: numtab  
  integer,intent(inout) :: nbuf(smax),itab(lsize,smax),jtab(lsize,smax),nhij(smax),it(smax),jt(smax),jtbuf(3)
  integer,intent(inout) :: nbufbuf(3),itabbuf(ndimh,3),jtabbuf(ndimh,3)  
  real*8,intent(inout) :: vbuf(lsize,smax)
  real*8,intent(inout) :: vbufbuf(ndimh,3)  
  integer,intent(in) :: nelac,ni,nj,ii,jj,ncf,isymat(smax)
  real*8,intent(in) :: hijsym(10,10)
  real*8 :: fac(2,10),hij
  integer :: ss,i,j,ndiff,ij,combi
  real*8 :: val2d(2)
  real*8 :: val3d1(3),val3d2(3),val3q1(3)
  real*8 :: val4dsinglet1(6),val4dsinglet2(6),val4dtriplet1(6),val4dtriplet2(6),val4dtriplet3(6),val4dquintet(6)
  
  fac=0d0
  select case(nelac)
  
!      ******************************** 2 e- ****************************************  
  case (2)
    do ss=1,3,2
      if(isymat(ss).eq.0) cycle
      combi=ss+10*nj+100*ni
      select case(combi)
      
      case(111)
        hij=0d0   
        if(jt(1).eq.0) then
          it(1)=it(1)+1
          numtab(it(1),1)=ii
        endif
        jt(1)=jt(1)+1        
        hij = hijsym(1,1)
        if(ii.eq.jj) then
          hij=0.5d0*hij
        else
          if(dabs(hij).lt.1.d-10) cycle      
        endif
        nhij(1)=nhij(1)+1          
        nbuf(1)=nbuf(1)+1
        itab(nbuf(1),1)=it(1)
        jtab(nbuf(1),1)=jt(1)
        vbuf(nbuf(1),1)=hij 
        if(nbuf(1).eq.lsize) then
          write(70+1) itab(1:lsize,1),jtab(1:lsize,1),vbuf(1:lsize,1)
          nbuf(1)=0
        endif
        
      case(211)
        hij=0d0   
        if(jt(1).eq.0) then
          it(1)=it(1)+1
          numtab(it(1),1)=ii
        endif
        jt(1)=jt(1)+1    
        do i=1,ni
          hij = hij + hijsym(i,1)/sqrt(2.)
        enddo
        if(dabs(hij).lt.1.d-10) cycle      
        nhij(1)=nhij(1)+1          
        nbuf(1)=nbuf(1)+1
        itab(nbuf(1),1)=it(1)
        jtab(nbuf(1),1)=jt(1)
        vbuf(nbuf(1),1)=hij 
        if(nbuf(1).eq.lsize) then
          write(70+1) itab(1:lsize,1),jtab(1:lsize,1),vbuf(1:lsize,1)
          nbuf(1)=0
        endif
        
      case(121)
        hij=0d0   
        if(jt(1).eq.0) then
          it(1)=it(1)+1
          numtab(it(1),1)=ii
        endif
        jt(1)=jt(1)+1    
        do j=1,nj
          hij = hij + hijsym(1,j)/sqrt(2.)
        enddo
        if(dabs(hij).lt.1.d-10) cycle      
        nhij(1)=nhij(1)+1          
        nbuf(1)=nbuf(1)+1
        itab(nbuf(1),1)=it(1)
        jtab(nbuf(1),1)=jt(1)
        vbuf(nbuf(1),1)=hij 
        if(nbuf(1).eq.lsize) then
          write(70+1) itab(1:lsize,1),jtab(1:lsize,1),vbuf(1:lsize,1)
          nbuf(1)=0
        endif
        
      case(221)
        hij=0d0   
        if(jt(1).eq.0) then
          it(1)=it(1)+1
          numtab(it(1),1)=ii
        endif
        jt(1)=jt(1)+1    
        do i=1,ni
          do j=1,nj
            hij = hij + hijsym(i,j)/2.
          enddo
        enddo  
        if(ii.eq.jj) then
          hij=0.5d0*hij
        else
          if(dabs(hij).lt.1.d-10) cycle      
        endif   
        nhij(1)=nhij(1)+1          
        nbuf(1)=nbuf(1)+1
        itab(nbuf(1),1)=it(1)
        jtab(nbuf(1),1)=jt(1)
        vbuf(nbuf(1),1)=hij 
        if(nbuf(1).eq.lsize) then
          write(70+1) itab(1:lsize,1),jtab(1:lsize,1),vbuf(1:lsize,1)
          nbuf(1)=0
        endif
      
      case(223)
       hij=0d0   
        if(jt(3).eq.0) then
          it(3)=it(3)+1
          numtab(it(3),3)=ii
        endif
        jt(3)=jt(3)+1    
        if(ordertab(ii,1).eq.1) then
          fac(1,1)=1.
          fac(1,2)=-1.
        else
          fac(1,1)=-1.
          fac(1,2)=1.
        endif
        if(ordertab(jj,1).eq.1) then
          fac(2,1)=1.
          fac(2,2)=-1.
        else
          fac(2,1)=-1.
          fac(2,2)=1.
        endif
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)/2.
          enddo
        enddo  
        if(ii.eq.jj) then
          hij=0.5d0*hij
        else
          if(dabs(hij).lt.1.d-10) cycle      
        endif   
        nhij(3)=nhij(3)+1          
        nbuf(3)=nbuf(3)+1
        itab(nbuf(3),3)=it(3)
        jtab(nbuf(3),3)=jt(3)
        vbuf(nbuf(3),3)=hij 
        if(nbuf(3).eq.lsize) then
          write(70+3) itab(1:lsize,3),jtab(1:lsize,3),vbuf(1:lsize,3)
          nbuf(3)=0
        endif
        
      case default
        cycle
      end select  
    enddo  
    
    if(ii.eq.jj) then
      jt=0
    endif  
        
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
        
        hij = hijsym(1,1)
        if(ii.eq.jj) then
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij  
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif    
          
        else
        
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif      
            
          endif  
        endif
        
! *****  NELAC = 3 ; SS = 2 ; Ni=1 ; Nj=3 ******************************************            
      case(132)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0   
        val3d1=(/2.,1.,-1./)
        val3d1=val3d1/sqrt(6d0)
        val3d2=(/0.,1.,1./)
        val3d2=val3d2/sqrt(2d0)
        
        fac(2,1)=val3d1(ordertab(jj,1))
        fac(2,2)=val3d1(ordertab(jj,2))
        fac(2,3)=val3d1(ordertab(jj,3))
          
        do j=1,nj
          hij = hij + fac(2,j)*hijsym(1,j)
        enddo  
        
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        endif  
      
!       second determinant on the column
        jt(ss)=jt(ss)+1 
        fac(2,1)=val3d2(ordertab(jj,1))
        fac(2,2)=val3d2(ordertab(jj,2))
        fac(2,3)=val3d2(ordertab(jj,3))
          
        hij=0d0
        do j=1,nj
          hij = hij + fac(2,j)*hijsym(1,j)
        enddo
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij 
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif
        endif  

! *****  NELAC = 3 ; SS = 2 ; Ni=3 ; Nj=1 ******************************************    
      case(312)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0   
        val3d1=(/2.,1.,-1./)
        val3d1=val3d1/sqrt(6d0)
        val3d2=(/0.,1.,1./)
        val3d2=val3d2/sqrt(2d0)
          
        fac(1,1)=val3d1(ordertab(ii,1))
        fac(1,2)=val3d1(ordertab(ii,2))
        fac(1,3)=val3d1(ordertab(ii,3))      
      
        do i=1,ni
          hij = hij + fac(1,i)*hijsym(i,1)
        enddo  

        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        endif  
        
!       second determinant on line : store in buffer        
        if(jtbuf(1).eq.0) numtab(it(2)+1,2)=ii
        jtbuf(1)=jtbuf(1)+1
        fac(1,1)=val3d2(ordertab(ii,1))
        fac(1,2)=val3d2(ordertab(ii,2))
        fac(1,3)=val3d2(ordertab(ii,3))
      
        hij=0d0
        do i=1,ni
          hij = hij + fac(1,i)*hijsym(i,1)
        enddo
        if(dabs(hij).ge.1.d-10) then
           nhij(ss)=nhij(ss)+1          
           nbufbuf(1)=nbufbuf(1)+1
           itabbuf(nbufbuf(1),1)=it(ss)+1
           jtabbuf(nbufbuf(1),1)=jtbuf(1)
           vbufbuf(nbufbuf(1),1)=hij
         endif  
          
! *****  NELAC = 3 ; SS = 2 ; Ni=3 ; Nj=3 ******************************************    
      case(332)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0   
        val3d1=(/2.,1.,-1./)
        val3d1=val3d1/sqrt(6d0)
        val3d2=(/0.,1.,1./)
        val3d2=val3d2/sqrt(2d0)
    
        fac(1,1)=val3d1(ordertab(ii,1))
        fac(1,2)=val3d1(ordertab(ii,2))
        fac(1,3)=val3d1(ordertab(ii,3))      
        fac(2,1)=val3d1(ordertab(jj,1))
        fac(2,2)=val3d1(ordertab(jj,2))
        fac(2,3)=val3d1(ordertab(jj,3))
      
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(ii.eq.jj) then    
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij  
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif        
        else
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif          
          endif  
        endif
      
!       second determinant on the column
        if(jt(ss).lt.it(ss)) then
          jt(ss)=jt(ss)+1 
          fac(2,1)=val3d2(ordertab(jj,1))
          fac(2,2)=val3d2(ordertab(jj,2))
          fac(2,3)=val3d2(ordertab(jj,3))
          
          hij=0d0
          do i=1,ni
            do j=1,nj
              hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
            enddo
          enddo  
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij 
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif
          endif 
        endif  
        
!       second determinant on line : store in buffer        
        if(jtbuf(1).eq.0) numtab(it(2)+1,2)=ii
        jtbuf(1)=jtbuf(1)+1
        fac(1,1)=val3d2(ordertab(ii,1))
        fac(1,2)=val3d2(ordertab(ii,2))
        fac(1,3)=val3d2(ordertab(ii,3))
        fac(2,1)=val3d1(ordertab(jj,1))
        fac(2,2)=val3d1(ordertab(jj,2))
        fac(2,3)=val3d1(ordertab(jj,3))
      
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbufbuf(1)=nbufbuf(1)+1
          itabbuf(nbufbuf(1),1)=it(ss)+1
          jtabbuf(nbufbuf(1),1)=jtbuf(1)
          vbufbuf(nbufbuf(1),1)=hij
        endif  
     
        jtbuf(1)=jtbuf(1)+1
        fac(2,1)=val3d2(ordertab(jj,1))
        fac(2,2)=val3d2(ordertab(jj,2))
        fac(2,3)=val3d2(ordertab(jj,3))
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(ii.eq.jj) then   
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbufbuf(1)=nbufbuf(1)+1
          itabbuf(nbufbuf(1),1)=it(ss)+1
          jtabbuf(nbufbuf(1),1)=jtbuf(1)
          vbufbuf(nbufbuf(1),1)=hij       
        else
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbufbuf(1)=nbufbuf(1)+1
            itabbuf(nbufbuf(1),1)=it(ss)+1
            jtabbuf(nbufbuf(1),1)=jtbuf(1)
            vbufbuf(nbufbuf(1),1)=hij 
          endif
        endif
        
! *****  NELAC = 3 ; SS = 4 ; Ni=3 ; Nj=3 ******************************************         
      case(334)
      !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
       
        val3q1=(/1d0,-1d0,1d0/)
        fac(1,1)=val3q1(ordertab(ii,1))
        fac(1,2)=val3q1(ordertab(ii,2))
        fac(1,3)=val3q1(ordertab(ii,3))
        fac(2,1)=val3q1(ordertab(jj,1))
        fac(2,2)=val3q1(ordertab(jj,2))
        fac(2,3)=val3q1(ordertab(jj,3))
        
        hij=0d0        
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)* hijsym(i,j)/3d0
          enddo
        enddo  
        if(ii.eq.jj) then
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij  
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif        
        else
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif          
          endif  
        endif
      
! *****  NELAC = 3 ; end select on SS;Ni;Nj ******************************************        
      case default
        cycle
      end select               
    enddo

    if(ii.eq.jj) then
!     end of line, clear the 'doublet' buffer    
      jt=0
      jtbuf=0
      if(nbufbuf(1).gt.0) then
        ij=0
        do 
          ndiff=min(lsize-nbuf(2),nbufbuf(1))
          do i=1,ndiff
            ij=ij+1
            nbuf(2)=nbuf(2)+1
            itab(nbuf(2),2)=itabbuf(ij,1)
            jtab(nbuf(2),2)=jtabbuf(ij,1)
            vbuf(nbuf(2),2)=vbufbuf(ij,1)
            nbufbuf(1)=nbufbuf(1)-1
          enddo  
          if(nbuf(2).eq.lsize) then
            write(70+2) itab(1:lsize,2),jtab(1:lsize,2),vbuf(1:lsize,2)
            nbuf(2)=0
          endif
          if(nbufbuf(1).le.0)exit
        enddo          
        nbufbuf(1)=0
        it(2)=it(2)+1
      endif
    endif  
!     ******************************** 4e- ****************************************    
  case(4)
    do ss=1,5,2
      if(isymat(ss).eq.0) cycle      

      combi=ss+10*nj+100*ni
      select case(combi)
      
! *****  NELAC = 4 ; SS = 1 ; Ni=1 ; Nj=1 ******************************************      
      case(111)
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij = hijsym(1,1)
        if(ii.eq.jj) then
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij  
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        else  
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif               
          endif  
        endif
! *****  NELAC = 4 ; SS = 1 ; Ni=1 ; Nj=2 ******************************************      
      case(121)
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0     
        hij = (hijsym(1,1)+hijsym(1,2))/sqrt(2.0)  
        
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        endif  
      
! *****  NELAC = 4 ; SS = 1 ; Ni=2 ; Nj=1 ******************************************      
      case(211)
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0     
        hij = (hijsym(1,1)+hijsym(2,1))/sqrt(2.0)  
        
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        endif  
! *****  NELAC = 4 ; SS = 1 ; Ni=2 ; Nj=2 ******************************************      
      case(221)
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0     
        hij = (hijsym(1,1)+hijsym(1,2)+hijsym(2,1)+hijsym(2,2))/2.0  
        
        if(ii.eq.jj) then
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij  
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        else  
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif               
          endif  
        endif
! *****  NELAC = 4 ; SS = 1 ; Ni=1 ; Nj=6 ******************************************      
      case(161)
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0     
        val4dsinglet1=(/1d0,1d0,0d0,0d0,1d0,1d0/)
        val4dsinglet1=val4dsinglet1/2d0
        val4dsinglet2=(/-1d0,1d0,2d0,2d0,1d0,-1d0/)
        val4dsinglet2=val4dsinglet2/sqrt(12d0)     
        
        fac(2,1)=val4dsinglet1(ordertab(jj,1))
        fac(2,2)=val4dsinglet1(ordertab(jj,2))
        fac(2,3)=val4dsinglet1(ordertab(jj,3))
        fac(2,4)=val4dsinglet1(ordertab(jj,4))
        fac(2,5)=val4dsinglet1(ordertab(jj,5))
        fac(2,6)=val4dsinglet1(ordertab(jj,6))        
          
        do j=1,nj
          hij = hij + fac(2,j)*hijsym(1,j)
        enddo  
        
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        endif  
      
!       second determinant on the column
        jt(ss)=jt(ss)+1 
        fac(2,1)=val4dsinglet2(ordertab(jj,1))
        fac(2,2)=val4dsinglet2(ordertab(jj,2))
        fac(2,3)=val4dsinglet2(ordertab(jj,3))
        fac(2,4)=val4dsinglet2(ordertab(jj,4))
        fac(2,5)=val4dsinglet2(ordertab(jj,5))
        fac(2,6)=val4dsinglet2(ordertab(jj,6))        
          
        hij=0d0
        do j=1,nj
          hij = hij + fac(2,j)*hijsym(1,j)
        enddo
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij 
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif
        endif  
! *****  NELAC = 4 ; SS = 1 ; Ni=6 ; Nj=1 ******************************************      
      case(611)
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0                             
        val4dsinglet1=(/1d0,1d0,0d0,0d0,1d0,1d0/)
        val4dsinglet1=val4dsinglet1/2d0
        val4dsinglet2=(/-1d0,1d0,2d0,2d0,1d0,-1d0/)
        val4dsinglet2=val4dsinglet2/sqrt(12d0)    
          
        fac(1,1)=val4dsinglet1(ordertab(ii,1))
        fac(1,2)=val4dsinglet1(ordertab(ii,2))
        fac(1,3)=val4dsinglet1(ordertab(ii,3))   
        fac(1,4)=val4dsinglet1(ordertab(ii,4))
        fac(1,5)=val4dsinglet1(ordertab(ii,5))
        fac(1,6)=val4dsinglet1(ordertab(ii,6))        
      
        do i=1,ni
          hij = hij + fac(1,i)*hijsym(i,1)
        enddo  

        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        endif  
        
!       second determinant on line : store in buffer        
        if(jtbuf(1).eq.0) numtab(it(ss)+1,ss)=ii
        jtbuf(1)=jtbuf(1)+1
        fac(1,1)=val4dsinglet2(ordertab(ii,1))
        fac(1,2)=val4dsinglet2(ordertab(ii,2))
        fac(1,3)=val4dsinglet2(ordertab(ii,3))
        fac(1,4)=val4dsinglet2(ordertab(ii,4))
        fac(1,5)=val4dsinglet2(ordertab(ii,5))
        fac(1,6)=val4dsinglet2(ordertab(ii,6))        
      
        hij=0d0
        do i=1,ni
          hij = hij + fac(1,i)*hijsym(i,1)
        enddo
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbufbuf(1)=nbufbuf(1)+1
          itabbuf(nbufbuf(1),1)=it(ss)+1
          jtabbuf(nbufbuf(1),1)=jtbuf(1)
          vbufbuf(nbufbuf(1),1)=hij
        endif  
! *****  NELAC = 4 ; SS = 1 ; Ni=6 ; Nj=6 ******************************************      
      case(661)
            !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0                   
        val4dsinglet1=(/1d0,1d0,0d0,0d0,1d0,1d0/)
        val4dsinglet1=val4dsinglet1/2d0
        val4dsinglet2=(/-1d0,1d0,2d0,2d0,1d0,-1d0/)
        val4dsinglet2=val4dsinglet2/sqrt(12d0)   
    
        fac(1,1)=val4dsinglet1(ordertab(ii,1))
        fac(1,2)=val4dsinglet1(ordertab(ii,2))
        fac(1,3)=val4dsinglet1(ordertab(ii,3))
        fac(1,4)=val4dsinglet1(ordertab(ii,4))
        fac(1,5)=val4dsinglet1(ordertab(ii,5))
        fac(1,6)=val4dsinglet1(ordertab(ii,6))            
        fac(2,1)=val4dsinglet1(ordertab(jj,1))
        fac(2,2)=val4dsinglet1(ordertab(jj,2))
        fac(2,3)=val4dsinglet1(ordertab(jj,3))
        fac(2,4)=val4dsinglet1(ordertab(jj,4))
        fac(2,5)=val4dsinglet1(ordertab(jj,5))
        fac(2,6)=val4dsinglet1(ordertab(jj,6))        
      
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(ii.eq.jj) then    
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij  
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif        
        else
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif          
          endif  
        endif
      
!       second determinant on the column
        if(jt(ss).lt.it(ss)) then
          jt(ss)=jt(ss)+1 
          fac(2,1)=val4dsinglet2(ordertab(jj,1))
          fac(2,2)=val4dsinglet2(ordertab(jj,2))
          fac(2,3)=val4dsinglet2(ordertab(jj,3))
          fac(2,4)=val4dsinglet2(ordertab(jj,4))
          fac(2,5)=val4dsinglet2(ordertab(jj,5))
          fac(2,6)=val4dsinglet2(ordertab(jj,6))          
          
          hij=0d0
          do i=1,ni
            do j=1,nj
              hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
            enddo
          enddo  
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij 
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif
          endif 
        endif  
        
!       second determinant on line : store in buffer        
        if(jtbuf(1).eq.0) numtab(it(ss)+1,ss)=ii
        jtbuf(1)=jtbuf(1)+1
        fac(1,1)=val4dsinglet2(ordertab(ii,1))
        fac(1,2)=val4dsinglet2(ordertab(ii,2))
        fac(1,3)=val4dsinglet2(ordertab(ii,3))
        fac(1,4)=val4dsinglet2(ordertab(ii,4))
        fac(1,5)=val4dsinglet2(ordertab(ii,5))
        fac(1,6)=val4dsinglet2(ordertab(ii,6))            
        fac(2,1)=val4dsinglet1(ordertab(jj,1))
        fac(2,2)=val4dsinglet1(ordertab(jj,2))
        fac(2,3)=val4dsinglet1(ordertab(jj,3))
        fac(2,4)=val4dsinglet1(ordertab(jj,4))
        fac(2,5)=val4dsinglet1(ordertab(jj,5))
        fac(2,6)=val4dsinglet1(ordertab(jj,6))
      
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbufbuf(1)=nbufbuf(1)+1
          itabbuf(nbufbuf(1),1)=it(ss)+1
          jtabbuf(nbufbuf(1),1)=jtbuf(1)
          vbufbuf(nbufbuf(1),1)=hij
        endif  
     
        jtbuf(1)=jtbuf(1)+1
        fac(2,1)=val4dsinglet2(ordertab(jj,1))
        fac(2,2)=val4dsinglet2(ordertab(jj,2))
        fac(2,3)=val4dsinglet2(ordertab(jj,3))
        fac(2,4)=val4dsinglet2(ordertab(jj,4))
        fac(2,5)=val4dsinglet2(ordertab(jj,5))
        fac(2,6)=val4dsinglet2(ordertab(jj,6))        
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(ii.eq.jj) then   
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbufbuf(1)=nbufbuf(1)+1
          itabbuf(nbufbuf(1),1)=it(ss)+1
          jtabbuf(nbufbuf(1),1)=jtbuf(1)
          vbufbuf(nbufbuf(1),1)=hij       
        else
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbufbuf(1)=nbufbuf(1)+1
            itabbuf(nbufbuf(1),1)=it(ss)+1
            jtabbuf(nbufbuf(1),1)=jtbuf(1)
            vbufbuf(nbufbuf(1),1)=hij 
          endif
        endif
! *****  NELAC = 4 ; SS = 1 ; Ni=2 ; Nj=6 ******************************************
      case(261)  
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0                     
        val4dsinglet1=(/1d0,1d0,0d0,0d0,1d0,1d0/)
        val4dsinglet1=val4dsinglet1/2d0
        val4dsinglet2=(/-1d0,1d0,2d0,2d0,1d0,-1d0/)
        val4dsinglet2=val4dsinglet2/sqrt(12d0)          
        
        fac(2,1)=val4dsinglet1(ordertab(jj,1))
        fac(2,2)=val4dsinglet1(ordertab(jj,2))
        fac(2,3)=val4dsinglet1(ordertab(jj,3))
        fac(2,4)=val4dsinglet1(ordertab(jj,4))
        fac(2,5)=val4dsinglet1(ordertab(jj,5))
        fac(2,6)=val4dsinglet1(ordertab(jj,6))        
        
        do i=1,ni
          do j=1,nj
            hij = hij + fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        hij=hij/sqrt(2d0)
        
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        endif  
      
!       second determinant on the column
        jt(ss)=jt(ss)+1 
        fac(2,1)=val4dsinglet2(ordertab(jj,1))
        fac(2,2)=val4dsinglet2(ordertab(jj,2))
        fac(2,3)=val4dsinglet2(ordertab(jj,3))
        fac(2,4)=val4dsinglet2(ordertab(jj,4))
        fac(2,5)=val4dsinglet2(ordertab(jj,5))
        fac(2,6)=val4dsinglet2(ordertab(jj,6))        
          
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(2,j)*hijsym(i,j)
          enddo
        enddo
        hij=hij/sqrt(2d0)
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij 
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif
        endif  
! *****  NELAC = 4 ; SS = 1 ; Ni=6 ; Nj=2 ****************************************** 
      case(621) 
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0                
        val4dsinglet1=(/1d0,1d0,0d0,0d0,1d0,1d0/)
        val4dsinglet1=val4dsinglet1/2d0
        val4dsinglet2=(/-1d0,1d0,2d0,2d0,1d0,-1d0/)
        val4dsinglet2=val4dsinglet2/sqrt(12d0)      
          
        fac(1,1)=val4dsinglet1(ordertab(ii,1))
        fac(1,2)=val4dsinglet1(ordertab(ii,2))
        fac(1,3)=val4dsinglet1(ordertab(ii,3))   
        fac(1,4)=val4dsinglet1(ordertab(ii,4))
        fac(1,5)=val4dsinglet1(ordertab(ii,5))
        fac(1,6)=val4dsinglet1(ordertab(ii,6))        
      
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*hijsym(i,j)
          enddo
        enddo  
        hij=hij/sqrt(2d0)

        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        endif  
        
!       second determinant on line : store in buffer        
        if(jtbuf(1).eq.0) numtab(it(ss)+1,ss)=ii
        jtbuf(1)=jtbuf(1)+1
        fac(1,1)=val4dsinglet2(ordertab(ii,1))
        fac(1,2)=val4dsinglet2(ordertab(ii,2))
        fac(1,3)=val4dsinglet2(ordertab(ii,3))
        fac(1,4)=val4dsinglet2(ordertab(ii,4))
        fac(1,5)=val4dsinglet2(ordertab(ii,5))
        fac(1,6)=val4dsinglet2(ordertab(ii,6))        
      
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*hijsym(i,j)
          enddo  
        enddo
        hij=hij/sqrt(2d0)
        
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbufbuf(1)=nbufbuf(1)+1
          itabbuf(nbufbuf(1),1)=it(ss)+1
          jtabbuf(nbufbuf(1),1)=jtbuf(1)
          vbufbuf(nbufbuf(1),1)=hij
        endif  
! *****  NELAC = 4 ; SS = 3 ; Ni=2 ; Nj=2 ******************************************      
      case(223)
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        val2d=(/1.0,-1.0/)
        
        fac(1,1)=val2d(ordertab(ii,1))
        fac(1,2)=val2d(ordertab(ii,2))
        fac(2,1)=val2d(ordertab(jj,1))
        fac(2,2)=val2d(ordertab(jj,2))
        
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)  
          enddo
        enddo
        hij=0.5d0*hij
        
        if(ii.eq.jj) then
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij  
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        else  
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif               
          endif  
        endif
! *****  NELAC = 4 ; SS = 3 ; Ni=2 ; Nj=6 ******************************************      
      case(263)
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0   
        val4dtriplet1=(/-2d0,-1d0,1d0,-1d0,1d0,2d0/)
        val4dtriplet1=val4dtriplet1/sqrt(12d0)    
        val4dtriplet2=(/-1d0,1d0,-1d0,1d0,-1d0,1d0/)
        val4dtriplet2=val4dtriplet2/sqrt(6d0)         
        val4dtriplet3=(/0d0,-1d0,-1d0,1d0,1d0,0d0/)
        val4dtriplet3=val4dtriplet3/2d0

!        val4dtriplet1=(/2d0,-1d0,-1d0,1d0,1d0,-2d0/)
!        val4dtriplet1=val4dtriplet1/sqrt(12d0)     
!        val4dtriplet2=(/-1d0,-1d0,-1d0,1d0,1d0,1d0/)
!        val4dtriplet2=val4dtriplet2/sqrt(6d0)         
!        val4dtriplet3=(/0d0,-1d0,1d0,-1d0,1d0,0d0/)
!        val4dtriplet3=val4dtriplet3/sqrt(4d0)      
        val2d=(/1.0,-1.0/)
        
        fac(1,1)=val2d(ordertab(ii,1))
        fac(1,2)=val2d(ordertab(ii,2))    
        fac(2,1)=val4dtriplet1(ordertab(jj,1))
        fac(2,2)=val4dtriplet1(ordertab(jj,2))
        fac(2,3)=val4dtriplet1(ordertab(jj,3))
        fac(2,4)=val4dtriplet1(ordertab(jj,4))
        fac(2,5)=val4dtriplet1(ordertab(jj,5))
        fac(2,6)=val4dtriplet1(ordertab(jj,6))        
        
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        endif  
      
!       second determinant on the column
        jt(ss)=jt(ss)+1 
        fac(2,1)=val4dtriplet2(ordertab(jj,1))
        fac(2,2)=val4dtriplet2(ordertab(jj,2))
        fac(2,3)=val4dtriplet2(ordertab(jj,3))
        fac(2,4)=val4dtriplet2(ordertab(jj,4))
        fac(2,5)=val4dtriplet2(ordertab(jj,5))
        fac(2,6)=val4dtriplet2(ordertab(jj,6))        
          
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo

        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij 
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif
        endif  
        
!       third determinant on the column
        jt(ss)=jt(ss)+1 
        fac(2,1)=val4dtriplet3(ordertab(jj,1))
        fac(2,2)=val4dtriplet3(ordertab(jj,2))
        fac(2,3)=val4dtriplet3(ordertab(jj,3))
        fac(2,4)=val4dtriplet3(ordertab(jj,4))
        fac(2,5)=val4dtriplet3(ordertab(jj,5))
        fac(2,6)=val4dtriplet3(ordertab(jj,6))        
          
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo

        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij 
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif
        endif          
! *****  NELAC = 4 ; SS = 3 ; Ni=6 ; Nj=2 ******************************************      
      case(623)
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        hij=0d0     
        val4dtriplet1=(/-2d0,-1d0,1d0,-1d0,1d0,2d0/)
        val4dtriplet1=val4dtriplet1/sqrt(12d0)     
        val4dtriplet2=(/-1d0,1d0,-1d0,1d0,-1d0,1d0/)
        val4dtriplet2=val4dtriplet2/sqrt(6d0)         
        val4dtriplet3=(/0d0,-1d0,-1d0,1d0,1d0,0d0/)
        val4dtriplet3=val4dtriplet3/2d0
!        val4dtriplet1=(/2d0,-1d0,-1d0,1d0,1d0,-2d0/)
!        val4dtriplet1=val4dtriplet1/sqrt(12d0)     
!        val4dtriplet2=(/-1d0,-1d0,-1d0,1d0,1d0,1d0/)
!        val4dtriplet2=val4dtriplet2/sqrt(6d0)         
!        val4dtriplet3=(/0d0,-1d0,1d0,-1d0,1d0,0d0/)
!        val4dtriplet3=val4dtriplet3/sqrt(4d0)     
        val2d=(/1.0,-1.0/)
        
        fac(2,1)=val2d(ordertab(jj,1))
        fac(2,2)=val2d(ordertab(jj,2))    
        fac(1,1)=val4dtriplet1(ordertab(ii,1))
        fac(1,2)=val4dtriplet1(ordertab(ii,2))
        fac(1,3)=val4dtriplet1(ordertab(ii,3))
        fac(1,4)=val4dtriplet1(ordertab(ii,4))
        fac(1,5)=val4dtriplet1(ordertab(ii,5))
        fac(1,6)=val4dtriplet1(ordertab(ii,6))           
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif          
        endif  

!       second determinant on line : store in buffer        
        if(jtbuf(2).eq.0) numtab(it(ss)+1,ss)=ii
        jtbuf(2)=jtbuf(2)+1
        fac(1,1)=val4dtriplet2(ordertab(ii,1))
        fac(1,2)=val4dtriplet2(ordertab(ii,2))
        fac(1,3)=val4dtriplet2(ordertab(ii,3))
        fac(1,4)=val4dtriplet2(ordertab(ii,4))
        fac(1,5)=val4dtriplet2(ordertab(ii,5))
        fac(1,6)=val4dtriplet2(ordertab(ii,6))        
      
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo  
        enddo

        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbufbuf(2)=nbufbuf(2)+1
          itabbuf(nbufbuf(2),2)=it(ss)+1
          jtabbuf(nbufbuf(2),2)=jtbuf(2)
          vbufbuf(nbufbuf(2),2)=hij
        endif  
        
!       third determinant on line : store in buffer        
        if(jtbuf(3).eq.0) numtab(it(ss)+2,ss)=ii
        jtbuf(3)=jtbuf(3)+1
        fac(1,1)=val4dtriplet3(ordertab(ii,1))
        fac(1,2)=val4dtriplet3(ordertab(ii,2))
        fac(1,3)=val4dtriplet3(ordertab(ii,3))
        fac(1,4)=val4dtriplet3(ordertab(ii,4))
        fac(1,5)=val4dtriplet3(ordertab(ii,5))
        fac(1,6)=val4dtriplet3(ordertab(ii,6))        

        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo  
        enddo

        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbufbuf(3)=nbufbuf(3)+1
          itabbuf(nbufbuf(3),3)=it(ss)+2
          jtabbuf(nbufbuf(3),3)=jtbuf(3)
          vbufbuf(nbufbuf(3),3)=hij
        endif  

! *****  NELAC = 4 ; SS = 3 ; Ni=6 ; Nj=6 ******************************************      
      case(663)
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
              write(*,*) '***************************************************'
      write(*,*) ii,jj
        write(*,'(20i20)') (ordertab(ii,j),j=1,6)
        write(*,'(20i20)') (ordertab(jj,j),j=1,6)        
      do i=1,6
        write(*,'(20f20.12)') (hijsym(i,j),j=1,6)
      enddo
      write(*,*) '***************************************************'
      
        hij=0d0     
        val4dtriplet1=(/-2d0,-1d0,1d0,-1d0,1d0,2d0/)
        val4dtriplet1=val4dtriplet1/sqrt(12d0)     
        val4dtriplet2=(/-1d0,1d0,-1d0,1d0,-1d0,1d0/)
        val4dtriplet2=val4dtriplet2/sqrt(6d0)         
        val4dtriplet3=(/0d0,-1d0,-1d0,1d0,1d0,0d0/)
        val4dtriplet3=val4dtriplet3/2d0
!        val4dtriplet1=(/2d0,-1d0,-1d0,1d0,1d0,-2d0/)
!        val4dtriplet1=val4dtriplet1/sqrt(12d0)     
!        val4dtriplet2=(/-1d0,-1d0,-1d0,1d0,1d0,1d0/)
!        val4dtriplet2=val4dtriplet2/sqrt(6d0)         
!        val4dtriplet3=(/0d0,-1d0,1d0,-1d0,1d0,0d0/)
!        val4dtriplet3=val4dtriplet3/sqrt(4d0)    
    
        fac(1,1)=val4dtriplet1(ordertab(ii,1))
        fac(1,2)=val4dtriplet1(ordertab(ii,2))
        fac(1,3)=val4dtriplet1(ordertab(ii,3))
        fac(1,4)=val4dtriplet1(ordertab(ii,4))
        fac(1,5)=val4dtriplet1(ordertab(ii,5))
        fac(1,6)=val4dtriplet1(ordertab(ii,6))            
        fac(2,1)=val4dtriplet1(ordertab(jj,1))
        fac(2,2)=val4dtriplet1(ordertab(jj,2))
        fac(2,3)=val4dtriplet1(ordertab(jj,3))
        fac(2,4)=val4dtriplet1(ordertab(jj,4))
        fac(2,5)=val4dtriplet1(ordertab(jj,5))
        fac(2,6)=val4dtriplet1(ordertab(jj,6))        
      
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(ii.eq.jj) then    
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij  
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif        
        else
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif          
          endif  
        endif
      
!       second determinant on the column
        if(jt(ss).lt.it(ss)) then
          jt(ss)=jt(ss)+1 
          fac(2,1)=val4dtriplet2(ordertab(jj,1))
          fac(2,2)=val4dtriplet2(ordertab(jj,2))
          fac(2,3)=val4dtriplet2(ordertab(jj,3))
          fac(2,4)=val4dtriplet2(ordertab(jj,4))
          fac(2,5)=val4dtriplet2(ordertab(jj,5))
          fac(2,6)=val4dtriplet2(ordertab(jj,6))          
          
          hij=0d0
          do i=1,ni
            do j=1,nj
              hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
            enddo
          enddo  
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij 
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif
          endif 
        endif  
        
!       third determinant on the column
        if(jt(ss).lt.it(ss)) then
          jt(ss)=jt(ss)+1 
          fac(2,1)=val4dtriplet3(ordertab(jj,1))
          fac(2,2)=val4dtriplet3(ordertab(jj,2))
          fac(2,3)=val4dtriplet3(ordertab(jj,3))
          fac(2,4)=val4dtriplet3(ordertab(jj,4))
          fac(2,5)=val4dtriplet3(ordertab(jj,5))
          fac(2,6)=val4dtriplet3(ordertab(jj,6))          
          
          hij=0d0
          do i=1,ni
            do j=1,nj
              hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
            enddo
          enddo  
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij 
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif
          endif 
        endif  
        
!       second determinant on line : store in buffer        
        if(jtbuf(2).eq.0) numtab(it(ss)+1,ss)=ii
        jtbuf(2)=jtbuf(2)+1
        fac(1,1)=val4dtriplet2(ordertab(ii,1))
        fac(1,2)=val4dtriplet2(ordertab(ii,2))
        fac(1,3)=val4dtriplet2(ordertab(ii,3))
        fac(1,4)=val4dtriplet2(ordertab(ii,4))
        fac(1,5)=val4dtriplet2(ordertab(ii,5))
        fac(1,6)=val4dtriplet2(ordertab(ii,6))            
        fac(2,1)=val4dtriplet1(ordertab(jj,1))
        fac(2,2)=val4dtriplet1(ordertab(jj,2))
        fac(2,3)=val4dtriplet1(ordertab(jj,3))
        fac(2,4)=val4dtriplet1(ordertab(jj,4))
        fac(2,5)=val4dtriplet1(ordertab(jj,5))
        fac(2,6)=val4dtriplet1(ordertab(jj,6))
      
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbufbuf(2)=nbufbuf(2)+1
          itabbuf(nbufbuf(2),2)=it(ss)+1
          jtabbuf(nbufbuf(2),2)=jtbuf(2)
          vbufbuf(nbufbuf(2),2)=hij
        endif  
     
        jtbuf(2)=jtbuf(2)+1
        fac(2,1)=val4dtriplet2(ordertab(jj,1))
        fac(2,2)=val4dtriplet2(ordertab(jj,2))
        fac(2,3)=val4dtriplet2(ordertab(jj,3))
        fac(2,4)=val4dtriplet2(ordertab(jj,4))
        fac(2,5)=val4dtriplet2(ordertab(jj,5))
        fac(2,6)=val4dtriplet2(ordertab(jj,6))        
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(ii.eq.jj) then   
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbufbuf(2)=nbufbuf(2)+1
          itabbuf(nbufbuf(2),2)=it(ss)+1
          jtabbuf(nbufbuf(2),2)=jtbuf(2)
          vbufbuf(nbufbuf(2),2)=hij       
        else
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbufbuf(2)=nbufbuf(2)+1
            itabbuf(nbufbuf(2),2)=it(ss)+1
            jtabbuf(nbufbuf(2),2)=jtbuf(2)
            vbufbuf(nbufbuf(2),2)=hij 
          endif
        endif
        
        if(.not.(ii.eq.jj)) then
          jtbuf(2)=jtbuf(2)+1
          fac(2,1)=val4dtriplet3(ordertab(jj,1))
          fac(2,2)=val4dtriplet3(ordertab(jj,2))
          fac(2,3)=val4dtriplet3(ordertab(jj,3))
          fac(2,4)=val4dtriplet3(ordertab(jj,4))
          fac(2,5)=val4dtriplet3(ordertab(jj,5))
          fac(2,6)=val4dtriplet3(ordertab(jj,6))        
          hij=0d0
          do i=1,ni
            do j=1,nj
              hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
            enddo
          enddo  
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbufbuf(2)=nbufbuf(2)+1
            itabbuf(nbufbuf(2),2)=it(ss)+1
            jtabbuf(nbufbuf(2),2)=jtbuf(2)
            vbufbuf(nbufbuf(2),2)=hij 
          endif
        endif  
        
!       third determinant on line : store in buffer        
        if(jtbuf(3).eq.0) numtab(it(ss)+2,ss)=ii
        jtbuf(3)=jtbuf(3)+1
        fac(1,1)=val4dtriplet3(ordertab(ii,1))
        fac(1,2)=val4dtriplet3(ordertab(ii,2))
        fac(1,3)=val4dtriplet3(ordertab(ii,3))
        fac(1,4)=val4dtriplet3(ordertab(ii,4))
        fac(1,5)=val4dtriplet3(ordertab(ii,5))
        fac(1,6)=val4dtriplet3(ordertab(ii,6))            
        fac(2,1)=val4dtriplet1(ordertab(jj,1))
        fac(2,2)=val4dtriplet1(ordertab(jj,2))
        fac(2,3)=val4dtriplet1(ordertab(jj,3))
        fac(2,4)=val4dtriplet1(ordertab(jj,4))
        fac(2,5)=val4dtriplet1(ordertab(jj,5))
        fac(2,6)=val4dtriplet1(ordertab(jj,6))
      
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbufbuf(3)=nbufbuf(3)+1
          itabbuf(nbufbuf(3),3)=it(ss)+2
          jtabbuf(nbufbuf(3),3)=jtbuf(3)
          vbufbuf(nbufbuf(3),3)=hij
        endif  
     
        jtbuf(3)=jtbuf(3)+1
        fac(2,1)=val4dtriplet2(ordertab(jj,1))
        fac(2,2)=val4dtriplet2(ordertab(jj,2))
        fac(2,3)=val4dtriplet2(ordertab(jj,3))
        fac(2,4)=val4dtriplet2(ordertab(jj,4))
        fac(2,5)=val4dtriplet2(ordertab(jj,5))
        fac(2,6)=val4dtriplet2(ordertab(jj,6))        
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(dabs(hij).ge.1.d-10) then
          nhij(ss)=nhij(ss)+1          
          nbufbuf(3)=nbufbuf(3)+1
          itabbuf(nbufbuf(3),3)=it(ss)+2
          jtabbuf(nbufbuf(3),3)=jtbuf(3)
          vbufbuf(nbufbuf(3),3)=hij 
        endif
        
        jtbuf(3)=jtbuf(3)+1
        fac(2,1)=val4dtriplet3(ordertab(jj,1))
        fac(2,2)=val4dtriplet3(ordertab(jj,2))
        fac(2,3)=val4dtriplet3(ordertab(jj,3))
        fac(2,4)=val4dtriplet3(ordertab(jj,4))
        fac(2,5)=val4dtriplet3(ordertab(jj,5))
        fac(2,6)=val4dtriplet3(ordertab(jj,6))        
        hij=0d0
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        if(ii.eq.jj) then    
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbufbuf(3)=nbufbuf(3)+1
          itabbuf(nbufbuf(3),3)=it(ss)+2
          jtabbuf(nbufbuf(3),3)=jtbuf(3)
          vbufbuf(nbufbuf(3),3)=hij    
        else
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbufbuf(3)=nbufbuf(3)+1
            itabbuf(nbufbuf(3),3)=it(ss)+2
            jtabbuf(nbufbuf(3),3)=jtbuf(3)
            vbufbuf(nbufbuf(3),3)=hij 
          endif
        endif  
        

! *****  NELAC = 4 ; SS = 5 ; Ni=6 ; Nj=6 ******************************************      
      case(665)
        !      new line    
        if(jt(ss).eq.0) then
          it(ss)=it(ss)+1
          numtab(it(ss),ss)=ii        
        endif
        jt(ss)=jt(ss)+1
        
        hij=0d0         
        val4dquintet=(/1d0,-1d0,1d0,1d0,-1d0,1d0/)
        
        fac(1,1)=val4dquintet(ordertab(ii,1))
        fac(1,2)=val4dquintet(ordertab(ii,2))
        fac(1,3)=val4dquintet(ordertab(ii,3))
        fac(1,4)=val4dquintet(ordertab(ii,4))
        fac(1,5)=val4dquintet(ordertab(ii,5))
        fac(1,6)=val4dquintet(ordertab(ii,6))            
        fac(2,1)=val4dquintet(ordertab(jj,1))
        fac(2,2)=val4dquintet(ordertab(jj,2))
        fac(2,3)=val4dquintet(ordertab(jj,3))
        fac(2,4)=val4dquintet(ordertab(jj,4))
        fac(2,5)=val4dquintet(ordertab(jj,5))
        fac(2,6)=val4dquintet(ordertab(jj,6))        
      
        do i=1,ni
          do j=1,nj
            hij = hij + fac(1,i)*fac(2,j)*hijsym(i,j)
          enddo
        enddo  
        hij=hij/6d0
        if(ii.eq.jj) then    
          hij=0.5d0*hij
          nhij(ss)=nhij(ss)+1          
          nbuf(ss)=nbuf(ss)+1
          itab(nbuf(ss),ss)=it(ss)
          jtab(nbuf(ss),ss)=jt(ss)
          vbuf(nbuf(ss),ss)=hij  
          if(nbuf(ss).eq.lsize) then
            write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
            nbuf(ss)=0
          endif        
        else
          if(dabs(hij).ge.1.d-10) then
            nhij(ss)=nhij(ss)+1          
            nbuf(ss)=nbuf(ss)+1
            itab(nbuf(ss),ss)=it(ss)
            jtab(nbuf(ss),ss)=jt(ss)
            vbuf(nbuf(ss),ss)=hij
            if(nbuf(ss).eq.lsize) then
              write(70+ss) itab(1:lsize,ss),jtab(1:lsize,ss),vbuf(1:lsize,ss)
              nbuf(ss)=0
            endif          
          endif  
        endif
! *****  NELAC = 4 ; end select on SS;Ni;Nj ******************************************            
      end select
    enddo
      
    if(ii.eq.jj) then
!     end of line, clear the 'singlet' and 'triplet' buffers    
      jt=0
      jtbuf=0
      if(nbufbuf(1).gt.0) then
        ij=0
        do 
          ndiff=min(lsize-nbuf(1),nbufbuf(1))
          do i=1,ndiff
            ij=ij+1
            nbuf(1)=nbuf(1)+1
            itab(nbuf(1),1)=itabbuf(ij,1)
            jtab(nbuf(1),1)=jtabbuf(ij,1)
            vbuf(nbuf(1),1)=vbufbuf(ij,1)
            nbufbuf(1)=nbufbuf(1)-1
          enddo  
          if(nbuf(1).eq.lsize) then
            write(70+1) itab(1:lsize,1),jtab(1:lsize,1),vbuf(1:lsize,1)
            nbuf(1)=0
          endif
          if(nbufbuf(1).le.0)exit
        enddo          
        nbufbuf(1)=0
        it(1)=it(1)+1
      endif
      if(nbufbuf(2).gt.0) then
        ij=0
        do 
          ndiff=min(lsize-nbuf(3),nbufbuf(2))
          do i=1,ndiff
            ij=ij+1
            nbuf(3)=nbuf(3)+1
            itab(nbuf(3),3)=itabbuf(ij,2)
            jtab(nbuf(3),3)=jtabbuf(ij,2)
            vbuf(nbuf(3),3)=vbufbuf(ij,2)
            nbufbuf(2)=nbufbuf(2)-1
          enddo  
          if(nbuf(3).eq.lsize) then
            write(70+3) itab(1:lsize,3),jtab(1:lsize,3),vbuf(1:lsize,3)
            nbuf(3)=0
          endif
          if(nbufbuf(2).le.0)exit
        enddo          
        nbufbuf(2)=0
        it(3)=it(3)+1
      endif 
      if(nbufbuf(3).gt.0) then
        ij=0
        do 
          ndiff=min(lsize-nbuf(3),nbufbuf(3))
          do i=1,ndiff
            ij=ij+1
            nbuf(3)=nbuf(3)+1
            itab(nbuf(3),3)=itabbuf(ij,3)
            jtab(nbuf(3),3)=jtabbuf(ij,3)
            vbuf(nbuf(3),3)=vbufbuf(ij,3)
            nbufbuf(3)=nbufbuf(3)-1
          enddo  
          if(nbuf(3).eq.lsize) then
            write(70+3) itab(1:lsize,3),jtab(1:lsize,3),vbuf(1:lsize,3)
            nbuf(3)=0
          endif
          if(nbufbuf(3).le.0)exit
        enddo          
        nbufbuf(3)=0
        it(3)=it(3)+1
      endif          
    endif  
      
  case DEFAULT
       
     write(*,*) 'isymat work only for 2 to 4 e-'
     
  end select
  
  return
     
end subroutine