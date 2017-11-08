module grid_data
  integer, parameter :: dp=kind(1.d0)
  integer :: nz,ncoord
  real(dp),dimension(:),allocatable :: gridmin,gridmax,gridstep
  real(dp),dimension(:,:),allocatable :: coord,coord_dist
  real(dp),dimension(:,:,:),allocatable :: coord_cart
  save
end module

subroutine makegrid(ncore,M,molecaxis,refdist,group,grid_method,coordinate,lower_bound,grid_file)
   use grid_data
   use io_unit
    
   implicit none
   integer :: i,j,k,ncore,sgn,ngrid(ncoord),numbering(ncoord),molecaxis,npoint,axis(3)
   integer,dimension(:,:),allocatable :: step_rank
   real(dp),dimension(:,:),allocatable :: temp_coord   
   real(dp) :: M(ncore),Masseff,Mtot,d1,d2,d3,dist(ncoord),angle1,angle2,pi,refdist(ncoord),lower_bound,temp,temp2,kmat(3,3),xmat(6),massfact(9),tmp_val(3)
   real(dp) :: mi,mj,ml,mk
   character(len=80) :: fmt_output
   character(len=*) :: group,grid_method,coordinate,grid_file 
   real(dp) :: theta1,phi1,theta2,phi2,phi3
   logical :: logic
   
   ! LHS declaration
   real(dp) :: opt,phip,p,distance
   real(dp),dimension(:),allocatable :: vec

    
   pi=4.D0*DATAN(1.D0)
   p=50d0

! Automatic placement of molecular axis
  if(group=='CS'.or.group=='cs'.or.group=='CNV'.or.group=='cnv'.or.group=='CNH'.or.group=='cnh'.or.group=='DNH'.or.group=='dnh')then 
     axis(1)=2
     axis(2)=3
     axis(3)=1
     write( 6, * ) 'Molecular axis placed in x-direction due to the CS/CNV symmetry'
  else      
     axis(1)=1
     axis(2)=2
     axis(3)=3
     write( 6, * ) 'Molecular axis placed in z-direction'	
  endif
   
   if((grid_method.eq.'grid').or.(grid_method.eq.'GRID')) then
     nz=1
     do i=1,ncoord
       ngrid(i)=abs(gridmax(i)-gridmin(i))/gridstep(i)+1
       nz=nz*ngrid(i)
     enddo

     write( 6, * ) nz,' valeurs de grille demandees ' 
     allocate(coord(ncoord,nz+1),temp_coord(ncoord,nz+1),coord_dist(ncoord,nz+1)) ! JD 05/06 +1 for a possible reference state
     allocate(coord_cart(ncore,3,nz+1))
     allocate(step_rank(nz,ncoord))
     
     numbering = 1
     step_rank=1
     do i= 2, nz
       numbering(1) = numbering(1) + 1
       do j=1,ncoord-1
         if (numbering(j) .gt.ngrid(j)) then
           numbering(j)=1
           numbering(j+1)=numbering(j+1)+ 1
         endif
         step_rank(i,j)=numbering(j)
       enddo
       step_rank(i,ncoord)=numbering(ncoord)
     enddo

   ! compute the coordinate for each point of the grid
     do j=1,ncoord
       do i=1,nz
         temp_coord(j,i) = gridmin(j)+((step_rank(i,j)-1)*gridstep(j))
       enddo
       temp_coord(j,nz+1) = refdist(j)
     enddo
     deallocate(step_rank)
     
   elseif((grid_method.eq.'list').or.(grid_method.eq.'LIST')) then  
     open(10,file=grid_file)
     read(10,*) nz
     write( 6, * ) nz,' valeurs de grille demandees '
     allocate(coord(ncoord,nz+1),temp_coord(ncoord,nz+1),coord_dist(ncoord,nz+1)) ! JD 05/06 +1 for a possible reference state
     allocate(coord_cart(ncore,3,nz+1))
     do i=1,nz
       read(10,*) (temp_coord(j,i),j=1,ncoord)
     enddo
     temp_coord(:,nz+1) = refdist(:)
     close(10)
     
   elseif((grid_method.eq.'LHS').or.(grid_method.eq.'lhs')) then
     allocate(coord(ncoord,nz+1),coord_dist(ncoord,nz+1))
     allocate(coord_cart(ncore,3,nz+1))
     allocate(vec(ncoord),temp_coord(ncoord,nz+1))
     
     ! initialize with a seed based on the system clock 
     call random_initialize(0)
     
     do k=1,200
       ! generate a LH sample
       call latin_center(ncoord,nz,temp_coord)
       
       ! compute the phi_p criterion
       phip=0
       do i=1,nz-1
         do j=i+1,nz
           vec(1:ncoord) = coord(1:ncoord,i) - coord(1:ncoord,j)
           distance = sqrt(dot_product(vec(1:ncoord),vec(1:ncoord)))
           distance=distance**p
           phip=phip + 1/distance 
         enddo
       enddo
       phip=phip**(1/p)
       ! pick the best sample (minimal phi_p) 
       
       if((phip.lt.opt).or.(k.eq.1)) then
         opt=phip
         do j=1,nz
           do i=1,ncoord
             temp_coord(i,j)=coord(i,j)
           enddo  
         enddo
       endif  
     enddo

     ! rescale values of the LHS
     do i=1,ncoord
       do j=1,nz
         temp_coord(i,j)=temp_coord(i,j)*(gridmax(i)-gridmin(i))+gridmin(i)
       enddo  
       temp_coord(i,nz+1) = refdist(i) 
     enddo
     deallocate(vec)
   endif
    
   ! check if reference configuration is needed (l_diab=true) 
   if(refdist(1)==-1d0) then
     npoint=nz
   else
     npoint=nz+1
   endif
   
   ! rescaled the mass
   Mtot=0d0
   do i=1,ncore
     Mtot=Mtot+M(i)
   enddo
   do i=1,ncore
     M(i)=M(i)/Mtot
   enddo  
   
   coord_cart=0d0
   ! Définition du tableau contenant toutes les valeurs de distance internucléaire
   open(io_output,file='index.dat')
   write(io_output,*) '# Constructing the grid with '//coordinate//' coordinate'
   i=0
   do k = 1, npoint      
     select case (ncore)
     case (1)
       i=i+1
       coord(:,i)=temp_coord(:,k)
       coord_dist(:,i)=coord(:,i)
       coord_cart(1,axis(1),i) = 0d0
       coord_cart(1,axis(2),i) = 0d0
       coord_cart(1,axis(3),i) = 0d0 
     case (2)
       i=i+1
       coord(:,i)=temp_coord(:,k)
       coord_dist(:,i)=coord(:,i)
       coord_cart(1,axis(1),i) = 0d0
       coord_cart(1,axis(2),i) = 0d0
       coord_cart(1,axis(3),i) = M(2)*coord(1,i)/(M(1)+M(2))
       coord_cart(2,axis(1),i) = 0d0
       coord_cart(2,axis(2),i) = 0d0
       coord_cart(2,axis(3),i) = -M(1)*coord(1,i)/(M(1)+M(2))        
     case (3)
       if((coordinate(1:1).eq.'j').or.(coordinate(1:1).eq.'J')) then        
          ! Jacobi : R,r,theta
          ! see Ragni2007
         angle1=temp_coord(3,k)*pi/180d0
         Masseff=M(3)+M(2)
         dist(3) = temp_coord(2,k)
         dist(1) = dsqrt((temp_coord(2,k)*M(3)/Masseff)**2-2.*(M(3)/Masseff)*temp_coord(2,k)*temp_coord(1,k)*dcos(angle1)+temp_coord(1,k)**2)
         dist(2) = dsqrt((temp_coord(2,k)*M(2)/Masseff)**2+2.*(M(2)/Masseff)*temp_coord(2,k)*temp_coord(1,k)*dcos(angle1)+temp_coord(1,k)**2)
       elseif((coordinate(1:1).eq.'i').or.(coordinate(1:1).eq.'I')) then  
         dist(:)=temp_coord(:,k)
       elseif((coordinate(1:1).eq.'h').or.(coordinate(1:1).eq.'H')) then 
         ! hyperspherical : rho,Theta,Phi
         ! see Johnson1980
         Masseff=dsqrt(M(1)*M(2)*M(3)/(M(1)+M(2)+M(3)))
         d1=dsqrt((M(1)/Masseff)*(1-M(1)/(M(1)+M(2)+M(3))))    
         d2=dsqrt((M(2)/Masseff)*(1-M(2)/(M(1)+M(2)+M(3)))) 
         d3=dsqrt((M(3)/Masseff)*(1-M(3)/(M(1)+M(2)+M(3)))) 
         angle1=temp_coord(2,k)*pi/180d0
         angle2=temp_coord(3,k)*pi/180d0
         dist(1)=(temp_coord(1,k)*d3)*dsqrt(abs((1+dcos(2*angle1)*dcos(2*angle2-2*atan(M(2)/Masseff)))/2))
         dist(2)=(temp_coord(1,k)*d2)*dsqrt(abs((1+dcos(2*angle1)*dcos(2*angle2+2*atan(M(3)/Masseff)))/2))
         dist(3)=(temp_coord(1,k)*d1)*dsqrt(abs((1+dcos(2*angle1)*dcos(2*angle2))/2))
       else
         write(*,*) 'Error in coordinate definition'
         stop
       endif  

     ! discard invalid geometries   
       if((dist(1).lt.lower_bound).or.(dist(2).lt.lower_bound).or.(dist(3).lt.lower_bound)) then
         nz=nz-1
         cycle
       elseif(2*maxval(dist).gt.sum(dist)+0.0001) then
         nz=nz-1
         cycle
       else
         i=i+1
         coord(:,i)=temp_coord(:,k)
         coord_dist(:,i)=dist(:)
       endif
       
       ! internuclear distance definition : r1=AB ; r2=AC ; r3=BC 
       ! convert to cartesian. definition : atoms on the XY plane, atom A on the x-axis.
       Masseff=M(1)+M(2)+M(3)
       dist(1)=dist(1)**2
       dist(2)=dist(2)**2
       dist(3)=dist(3)**2
       
       coord_cart(:,3,i) = 0d0
       temp=dsqrt(abs(dist(1)*M(2)**2+dist(2)*M(3)**2+M(2)*M(3)*(dist(1)+dist(2)-dist(3))))
       if(temp.eq.0d0) then
         coord_cart(:,1,i) = 0d0
         
         coord_cart(1,2,i) = 0d0    
         coord_cart(2,2,i) = -dist(1)
         coord_cart(3,2,i) = dist(2)             
       else       
         coord_cart(1,1,i) = - temp/Masseff
         coord_cart(1,2,i) = 0d0     
         coord_cart(2,1,i) = (M(1)*(2d0*M(2)*dist(1)+M(3)*(dist(1)+dist(2)-dist(3))) +M(3)*(M(3)*(dist(1)-dist(2)-dist(3)) + M(2)*(dist(1)-dist(2)+dist(3))))/(2d0*temp*Masseff)
         coord_cart(2,2,i) = -M(3)*(dsqrt(abs(dist(1)**2+(dist(2)-dist(3))**2-2d0*dist(1)*(dist(2)+dist(3)))))/(2d0*temp)
         coord_cart(3,1,i) = (M(1)*(2d0*M(3)*dist(2)+M(2)*(dist(1)+dist(2)-dist(3))) -M(2)*(M(3)*(dist(1)-dist(2)-dist(3)) + M(2)*(dist(1)-dist(2)+dist(3))))/(2d0*temp*Masseff)
         coord_cart(3,2,i) = -M(2)*coord_cart(2,2,i)/M(3)
       endif
       
     case (4)
       if((coordinate(1:1).eq.'j').or.(coordinate(1:1).eq.'J')) then  
        ! Jacobi adapted to the system AB + CD
        ! Jacobi : R, r12, r34 , theta1,theta2,phi2
        ! vector r12 and r34 defined by spherical coord. (See Ragni2007b, H-vector for picture/definition)
      
         theta1=temp_coord(4,k)*pi/180d0
         theta2=temp_coord(5,k)*pi/180d0
         phi2=temp_coord(6,k)*pi/180d0
         mi=M(1)/(M(1)+M(2))
         mj=M(2)/(M(1)+M(2))
         mk=M(3)/(M(3)+M(4))
         ml=M(4)/(M(3)+M(4))
         temp=Cos(phi2)*Sin(theta1)*Sin(theta2)-Cos(theta1)*Cos(theta2)

         dist(1) = temp_coord(2,k)
         dist(2) = sqrt(abs((mj*temp_coord(2,k))**2+(ml*temp_coord(3,k))**2+temp_coord(1,k)**2+2*temp_coord(1,k)*(mj*temp_coord(2,k)*Cos(theta1)+ml*temp_coord(3,k)*Cos(theta2))-2*mj*ml*temp_coord(2,k)*temp_coord(3,k)*temp))
         dist(3) = sqrt(abs((mj*temp_coord(2,k))**2+(mk*temp_coord(3,k))**2+temp_coord(1,k)**2+2*temp_coord(1,k)*(mj*temp_coord(2,k)*Cos(theta1)-mk*temp_coord(3,k)*Cos(theta2))+2*mj*mk*temp_coord(2,k)*temp_coord(3,k)*temp))
         dist(4) = sqrt(abs((mi*temp_coord(2,k))**2+(ml*temp_coord(3,k))**2+temp_coord(1,k)**2-2*temp_coord(1,k)*(mi*temp_coord(2,k)*Cos(theta1)-ml*temp_coord(3,k)*Cos(theta2))+2*mi*ml*temp_coord(2,k)*temp_coord(3,k)*temp))
         dist(5) = sqrt(abs((mi*temp_coord(2,k))**2+(mk*temp_coord(3,k))**2+temp_coord(1,k)**2-2*temp_coord(1,k)*(mi*temp_coord(2,k)*Cos(theta1)+mk*temp_coord(3,k)*Cos(theta2))-2*mi*mk*temp_coord(2,k)*temp_coord(3,k)*temp))
         dist(6) = temp_coord(3,k)
       elseif((coordinate(1:1).eq.'b').or.(coordinate(1:1).eq.'B')) then  
        ! Jacobi used for DNH point group : same as above but the jacobi are adapted to the system AC + BD
        ! Jacobi : R, r13, r24 , theta1,theta2,phi2

         theta1=temp_coord(4,k)*pi/180d0
         theta2=temp_coord(5,k)*pi/180d0
         phi2=temp_coord(6,k)*pi/180d0
         mi=M(1)/(M(1)+M(3))
         mj=M(3)/(M(1)+M(3))
         mk=M(2)/(M(2)+M(4))
         ml=M(4)/(M(2)+M(4))
         temp=Cos(phi2)*Sin(theta1)*Sin(theta2)-Cos(theta1)*Cos(theta2)        
      
         dist(2) = temp_coord(2,k)
         dist(1) = sqrt(abs((mj*temp_coord(2,k))**2+(ml*temp_coord(3,k))**2+temp_coord(1,k)**2+2*temp_coord(1,k)*(mj*temp_coord(2,k)*Cos(theta1)+ml*temp_coord(3,k)*Cos(theta2))-2*mj*ml*temp_coord(2,k)*temp_coord(3,k)*temp))
         dist(3) = sqrt(abs((mj*temp_coord(2,k))**2+(mk*temp_coord(3,k))**2+temp_coord(1,k)**2+2*temp_coord(1,k)*(mj*temp_coord(2,k)*Cos(theta1)-mk*temp_coord(3,k)*Cos(theta2))+2*mj*mk*temp_coord(2,k)*temp_coord(3,k)*temp))
         dist(4) = sqrt(abs((mi*temp_coord(2,k))**2+(ml*temp_coord(3,k))**2+temp_coord(1,k)**2-2*temp_coord(1,k)*(mi*temp_coord(2,k)*Cos(theta1)-ml*temp_coord(3,k)*Cos(theta2))+2*mi*ml*temp_coord(2,k)*temp_coord(3,k)*temp))
         dist(6) = sqrt(abs((mi*temp_coord(2,k))**2+(mk*temp_coord(3,k))**2+temp_coord(1,k)**2-2*temp_coord(1,k)*(mi*temp_coord(2,k)*Cos(theta1)+mk*temp_coord(3,k)*Cos(theta2))-2*mi*mk*temp_coord(2,k)*temp_coord(3,k)*temp))
         dist(5) = temp_coord(3,k)    
       elseif((coordinate(1:1).eq.'i').or.(coordinate(1:1).eq.'I')) then  
         dist(:)=temp_coord(:,k)
       elseif((coordinate(1:1).eq.'h').or.(coordinate(1:1).eq.'H')) then  
         ! Hyperspherical coordinate, see Aquilanti1997
         theta1=temp_coord(2,k)*pi/180d0
         theta2=temp_coord(3,k)*pi/180d0
         phi1=temp_coord(4,k)*pi/180d0         
         phi2=temp_coord(5,k)*pi/180d0
         phi3=temp_coord(6,k)*pi/180d0         
         
         kmat(1,1) = cos(phi1)*cos(phi2)*cos(phi3)-sin(phi1)*sin(phi3)
         kmat(1,2) = -cos(phi3)*sin(phi1)-cos(phi1)*cos(phi2)*sin(phi3)
         kmat(1,3) = cos(phi1)*sin(phi2)
         kmat(2,1) = sin(phi1)*cos(phi2)*cos(phi3)+cos(phi1)*sin(phi3)
         kmat(2,2) = -cos(phi3)*sin(phi2)
         kmat(2,3) = sin(phi1)*sin(phi2)
         kmat(3,1) = cos(phi3)*cos(phi1)-sin(phi1)*cos(phi2)*sin(phi3)
         kmat(3,2) = sin(phi2)*sin(phi3)
         kmat(3,3) = cos(phi2)
         
         xmat(1) = temp_coord(1,k)**2*(cos(theta1)**2+3*sin(theta1)**2*(sin(theta2)**2*kmat(2,1)**2+cos(theta2)**2*kmat(3,1)**2))/3d0
         xmat(2) = temp_coord(1,k)**2*(cos(theta1)**2+3*sin(theta1)**2*(sin(theta2)**2*kmat(2,2)**2+cos(theta2)**2*kmat(3,2)**2))/3d0       
         xmat(3) = temp_coord(1,k)**2*(cos(theta1)**2+3*sin(theta1)**2*(sin(theta2)**2*kmat(2,3)**2+cos(theta2)**2*kmat(3,3)**2))/3d0   
         xmat(4) = (temp_coord(1,k)*sin(theta1))**2*(sin(theta2)**2*kmat(2,1)*kmat(2,2)+cos(theta2)**2*kmat(3,1)*kmat(3,2))
         xmat(5) = (temp_coord(1,k)*sin(theta1))**2*(sin(theta2)**2*kmat(2,1)*kmat(2,3)+cos(theta2)**2*kmat(3,1)*kmat(3,3))
         xmat(6) = (temp_coord(1,k)*sin(theta1))**2*(sin(theta2)**2*kmat(2,2)*kmat(2,3)+cos(theta2)**2*kmat(3,2)*kmat(3,3))   
         
         massfact(1) = sqrt(M(1)*M(3)/(M(2)*(M(1)+M(2)+M(3))))
         massfact(2) = sqrt(M(2)*M(3)/(M(1)*(M(1)+M(2)+M(3))))   
         massfact(3) = sqrt(M(4)*(M(1)+M(2))/(M(3)*(M(1)+M(2)+M(3)+M(4))))    
         massfact(4) = sqrt(M(3)*M(4)/((M(1)+M(2))*(M(1)+M(2)+M(3)+M(4))))     
         massfact(5) = sqrt(M(2)*M(4)*(M(1)+M(2)+M(3))/(M(1)*(M(1)+M(2))*(M(1)+M(2)+M(3)+M(4)))) 
         massfact(6) = sqrt(M(1)*M(4)*(M(1)+M(2)+M(3))/(M(2)*(M(1)+M(2))*(M(1)+M(2)+M(3)+M(4)))) 
         massfact(7) = M(1)*M(2)/(M(1)+M(2))
         massfact(8) = M(3)*(M(1)+M(2))/(M(1)+M(2)+M(3))
         massfact(9) = M(4)*(M(1)+M(2)+M(3))/(M(1)+M(2)+M(3)+M(4))
         Masseff = (M(1)*M(2)*M(3)*M(4)/(M(1)+M(2)+M(3)+M(4)))**(1d0/3d0)
         
         dist(1) = sqrt(abs(xmat(1)*Masseff/massfact(7)))
         dist(2) = sqrt(abs((Masseff/massfact(8))*(xmat(2)+xmat(1)*massfact(2)**2+2*massfact(2)*xmat(4))))
         dist(3) = sqrt(abs((Masseff/massfact(9))*(xmat(3)+xmat(2)*massfact(4)**2+xmat(1)*massfact(5)**2+2*massfact(4)*xmat(6)+2*massfact(5)*xmat(5)+2*massfact(4)*massfact(5)*xmat(4))))  
         dist(4) = sqrt(abs((Masseff/massfact(8))*(xmat(2)+xmat(1)*massfact(1)**2-2*massfact(1)*xmat(4))))
         dist(5) = sqrt(abs((Masseff/massfact(9))*(xmat(3)+xmat(2)*massfact(4)**2+xmat(1)*massfact(6)**2+2*massfact(4)*xmat(6)-2*massfact(6)*xmat(5)-2*massfact(6)*massfact(4)*xmat(4))))  
         dist(6) = sqrt(abs((Masseff/massfact(9))*(xmat(3)+xmat(2)*massfact(3)**2-2*massfact(3)*xmat(6))))
       else
         write(*,*) 'Error in coordinate definition'
         stop
       endif
       
       ! check validity of coordinate : triangle inequality for each trimer combinaison
       logic=.false.
       do j=1,4
         select case(j)
         case(1)
           tmp_val(1)=dist(1)
           tmp_val(2)=dist(2)
           tmp_val(3)=dist(4)      
         case(2)
           tmp_val(1)=dist(1)
           tmp_val(2)=dist(3)
           tmp_val(3)=dist(5)     
         case(3)
           tmp_val(1)=dist(2)
           tmp_val(2)=dist(3)
           tmp_val(3)=dist(6)     
         case(4)
           tmp_val(1)=dist(4)
           tmp_val(2)=dist(5)
           tmp_val(3)=dist(6)     
         end select  
         if(2*maxval(tmp_val).gt.sum(tmp_val)+0.0001) logic=.true.
       enddo  

      ! discard invalid geometries   
       if((dist(1).lt.lower_bound).or.(dist(2).lt.lower_bound).or.(dist(3).lt.lower_bound).or.(dist(4).lt.lower_bound).or.(dist(5).lt.lower_bound).or.(dist(6).lt.lower_bound)) then
         nz=nz-1
         cycle
       elseif(logic) then
         nz=nz-1
         cycle
       else
         i=i+1
         coord(:,i)=temp_coord(:,k)
         coord_dist(:,i)=dist(:)
       endif

       ! internuclear distance definition : r1=AB ; r2=AC ; r3=AD ; r4=BC ; r5=BD ; r6= CD
       ! convert to cartesian. definition : atoms A and B on the XY plane, center-of-mass of AB (and CD) on the x-axis.
       Masseff=M(1)+M(2)+M(3)+M(4)
       dist(1)=dist(1)**2
       dist(2)=dist(2)**2
       dist(3)=dist(3)**2
       dist(4)=dist(4)**2
       dist(5)=dist(5)**2
       dist(6)=dist(6)**2       
!       temp=(M(2)*dist(1))**2+(M(3)*dist(2))**2+(M(4)*dist(3))**2+M(2)*M(3)*(dist(1)**2+dist(2)**2-dist(4)**2)+M(2)*M(4)*(dist(1)**2+dist(3)**2-dist(5)**2)+M(3)*M(4)*(dist(2)**2+dist(3)**2-dist(6)**2)
!temp2=-(M(3)**2*(dist(1)**4 + (dist(2)**2 - dist(4)**2)**2 - 2*dist(1)**2*(dist(2)**2 + dist(4)**2))) - M(4)**2*(dist(1)**4 + (dist(3)**2 - dist(5)**2)**2 - 2*dist(1)**2*(dist(3)**2 + dist(5)**2)) + 2*M(3)*M(4)*(-dist(1)**4 - (dist(2)**2 - dist(4)**2)*(dist(3)**2 - dist(5)**2) + dist(1)**2*(dist(2)**2 + dist(3)**2 + dist(4)**2 + dist(5)**2 - 2*dist(6)**2))
       temp= sqrt(abs(M(1)*M(2)*(M(3)**2*(-dist(1)+dist(2)+dist(4))+M(4)**2*(-dist(1)+dist(3)+dist(5))+M(3)*M(4)*(-2*dist(1)+dist(2)+dist(3)+dist(4)+dist(5)-2*dist(6)))+M(1)**2*(M(3)**2*dist(2)+M(4)**2*dist(3)+M(3)*M(4)*(dist(2)+dist(3)-dist(6)))+M(2)**2*(M(3)**2*dist(4)+M(4)**2*dist(5)+M(3)*M(4)*(dist(4)+dist(5)-dist(6))) ))
       temp2=sqrt(abs( M(3)**2*(dist(1)**2+(dist(2)-dist(4))**2-2*dist(1)*(dist(2)+dist(4)))+M(4)**2*(dist(1)**2+(dist(3)-dist(5))**2-2*dist(1)*(dist(3)+dist(5)))+2*M(3)*M(4)*(dist(1)**2+(dist(2)-dist(4))*(dist(3)-dist(5))-dist(1)*(dist(2)+dist(3)+dist(4)+dist(5)-2*dist(6)))))
       coord_cart(1,1,i) = (M(2)*(M(3)**2*(dist(1)-dist(2)-dist(4))+M(4)**2*(dist(1)-dist(3)-dist(5))+M(2)*(M(3)*(dist(1)-dist(2)+dist(4))+M(4)*(dist(1)-dist(3)+dist(5)))+M(3)*M(4)*(2*dist(1)-dist(2)-dist(3)-dist(4)-dist(5)+2*dist(6)))-M(1)*(M(2)*(M(3)*(dist(1)+dist(2)-dist(4))+M(4)*(dist(1)+dist(3)-dist(5)))+2*(M(3)**2*dist(2)+M(4)**2*dist(3)+M(3)*M(4)*(dist(2)+dist(3)-dist(6)))))/(2*temp*Masseff)
       coord_cart(1,2,i) = M(2)*temp2/(2*temp)
       coord_cart(1,3,i) = 0d0
       coord_cart(2,1,i) =  (M(1)**2*(M(3)*(dist(1)+dist(2)-dist(4))+M(4)*(dist(1)+dist(3)-dist(5)))-M(1)*(M(3)**2*(-dist(1)+dist(2)+dist(4))+M(4)**2*(-dist(1)+dist(3)+dist(5))+M(2)*(M(3)*(dist(1)-dist(2)+dist(4))+M(4)*(dist(1)-dist(3)+dist(5)))+M(3)*M(4)*(-2*dist(1)+dist(2)+dist(3)+dist(4)+dist(5)-2*dist(6)))-2*M(2)*(M(3)**2*dist(4)+M(4)**2*dist(5)+M(3)*M(4)*(dist(4)+dist(5)-dist(6))))/(2*temp*Masseff)  
       coord_cart(2,2,i) = -M(1)*coord_cart(1,2,i)/M(2)  
       coord_cart(2,3,i) = 0d0
       coord_cart(3,1,i) = (M(1)**2*(2*M(3)*dist(2)+M(4)*(dist(2)+dist(3)-dist(6)))+M(1)*(M(2)*(2*M(3)*(-dist(1)+dist(2)+dist(4))+M(4)*(-2*dist(1)+dist(2)+dist(3)+dist(4)+dist(5)-2*dist(6)))+M(4)*(M(4)*(dist(2)-dist(3)-dist(6))+M(3)*(dist(2)-dist(3)+dist(6))))+M(2)*(M(2)*(2*M(3)*dist(4)+M(4)*(dist(4)+dist(5)-dist(6)))+M(4)*(M(4)*(dist(4)-dist(5)-dist(6))+M(3)*(dist(4)-dist(5)+dist(6))))) /(2*temp*Masseff)                   
       
       coord_cart(4,1,i) = (M(1)**2*(2*M(4)*dist(3)+M(3)*(dist(2)+dist(3)-dist(6)))+M(1)*(M(2)*(2*M(4)*(-dist(1)+dist(3)+dist(5))+M(3)*(-2*dist(1)+dist(2)+dist(3)+dist(4)+dist(5)-2*dist(6)))+M(3)*(-M(3)*(dist(2)-dist(3)+dist(6))+M(4)*(-dist(2)+dist(3)+dist(6))))+M(2)*(M(2)*(2*M(4)*dist(5)+M(3)*(dist(4)+dist(5)-dist(6)))+M(3)*(-M(3)*(dist(4)-dist(5)+dist(6))+M(4)*(-dist(4)+dist(5)+dist(6)))))/(2*temp*Masseff)
       if(abs(temp2).lt.1d-5) then
         coord_cart(3,2,i) = 0d0
         coord_cart(3,3,i) = 0d0
         coord_cart(4,2,i) = 0d0
         coord_cart(4,3,i) = 0d0
       else
         coord_cart(3,2,i) =M(4)*(M(1)*(M(4)*(dist(3)**2+2*dist(3)*dist(4)-dist(3)*dist(5)-dist(2)*(dist(3)+dist(5))-dist(3)*dist(6)+dist(5)*dist(6)+dist(1)*(dist(2)-dist(3)-dist(6)))+M(3)*(-dist(2)**2+dist(4)*(dist(3)-dist(6))+dist(1)*(dist(2)-dist(3)+dist(6))+dist(2)*(dist(3)+dist(4)-2*dist(5)+dist(6)))) +M(2)*(M(4)*(dist(3)*(dist(4)+dist(5)-dist(6))+dist(5)*(-2*dist(2)+dist(4)-dist(5)+dist(6))+dist(1)*(-dist(4)+dist(5)+dist(6)))-M(3)*(dist(2)*(dist(4)+dist(5)-dist(6))+dist(1)*(dist(4)-dist(5)+dist(6))+dist(4)*(-2*dist(3)-dist(4)+dist(5)+dist(6)))))/(2*temp*temp2)                   
         coord_cart(3,3,i) =M(4)*sqrt(abs(dist(2)**2*dist(5)+dist(1)**2*dist(6)+dist(4)*(dist(3)**2+dist(5)*dist(6)+dist(3)*(dist(4)-dist(5)-dist(6)))+dist(1)*(dist(2)*(dist(4)-dist(5)-dist(6))+dist(6)*(-dist(4)-dist(5)+dist(6))-dist(3)*(dist(4)-dist(5)+dist(6)))-dist(2)*(dist(3)*(dist(4)+dist(5)-dist(6))+dist(5)*(dist(4)-dist(5)+dist(6)))))/temp2
         coord_cart(4,2,i) = -M(3)*coord_cart(3,2,i)/M(4)
         coord_cart(4,3,i) = -M(3)*coord_cart(3,3,i)/M(4)
       endif
!  write(*,*) 'test',sqrt(abs(dist(2)**2*dist(5)+dist(1)**2*dist(6)+dist(4)*(dist(3)**2+dist(5)*dist(6)+dist(3)*(dist(4)-dist(5)-dist(6)))+dist(1)*(dist(2)*(dist(4)-dist(5)-dist(6))+dist(6)*(-dist(4)-dist(5)+dist(6))-dist(3)*(dist(4)-dist(5)+dist(6)))-dist(2)*(dist(3)*(dist(4)+dist(5)-dist(6))+dist(5)*(dist(4)-dist(5)+dist(6))))),abs(dist(2)**2*dist(5)+dist(1)**2*dist(6)+dist(4)*(dist(3)**2+dist(5)*dist(6)+dist(3)*(dist(4)-dist(5)-dist(6)))+dist(1)*(dist(2)*(dist(4)-dist(5)-dist(6))+dist(6)*(-dist(4)-dist(5)+dist(6))-dist(3)*(dist(4)-dist(5)+dist(6)))-dist(2)*(dist(3)*(dist(4)+dist(5)-dist(6))+dist(5)*(dist(4)-dist(5)+dist(6)))),dist(2)**2*dist(5)+dist(1)**2*dist(6)+dist(4)*(dist(3)**2+dist(5)*dist(6)+dist(3)*(dist(4)-dist(5)-dist(6)))+dist(1)*(dist(2)*(dist(4)-dist(5)-dist(6))+dist(6)*(-dist(4)-dist(5)+dist(6))-dist(3)*(dist(4)-dist(5)+dist(6)))-dist(2)*(dist(3)*(dist(4)+dist(5)-dist(6))+dist(5)*(dist(4)-dist(5)+dist(6))),coord_cart(3,3,i),temp2
!        coord_cart(1,1,i) = - sqrt(abs(temp))/Masseff
!        coord_cart(1,2,i) = 0d0
!        coord_cart(1,3,i) = 0d0
!        coord_cart(2,1,i) = ((M(3)+M(4))*(M(3)*(dist(1)**2-dist(2)**2-dist(4)**2)+M(4)*(dist(1)**2-dist(3)**2-dist(5)**2))+M(1)*(2*M(2)*dist(1)**2+M(3)*(dist(1)**2+dist(2)**2-dist(4)**2)+M(4)*(dist(1)**2+dist(3)**2-dist(5)**2))+M(2)*(M(3)*(dist(1)**2-dist(2)**2+dist(4)**2)+M(4)*(dist(1)**2-dist(3)**2+dist(5)**2))+2*M(3)*M(4)*dist(6)**2)/(2*Masseff*sqrt(abs(temp)))
!      !  coord_cart(2,2,i) =-(sqrt(abs(-M(3)**2*(dist(1)**4+(dist(2)**2-dist(4)**2)**2-2*dist(1)**2*(dist(2)**2+dist(4)**2))-M(4)**2*(dist(1)**4+(dist(3)**2-dist(5)**2)**2-2*dist(1)**2*(dist(3)**2+dist(5)**2))+M(3)*M(4)*(-dist(1)**4-(dist(2)**2-dist(4)**2)*(dist(3)**2-dist(5)**2)+dist(1)**2*(dist(2)**2+dist(3)**2+dist(4)**2+dist(5)**2-2*dist(6)**2)))))/(2*sqrt(abs(temp)))
!        coord_cart(2,2,i) =-sqrt(abs(temp2))/(2*sqrt(abs(temp)))    
!        coord_cart(2,3,i) = 0d0
!        if(abs(temp2).eq.0d0) then
!          coord_cart(3,1,i) =(-M(2)**2*(dist(1)**2-dist(2)**2+dist(4)**2)+M(1)*(2*M(3)*dist(2)**2+M(2)*(dist(1)**2+dist(2)**2-dist(4)**2)+M(4)*(dist(2)**2+dist(3)**2-dist(6)**2))+M(4)*(M(4)*(dist(2)**2-dist(3)**2-dist(6)**2)+M(3)*(dist(2)**2-dist(3)**2+dist(6)**2))-M(2)*(M(3)*(dist(1)**2-dist(2)**2-dist(4)**2)+M(4)*(dist(1)**2-2*dist(2)**2+dist(3)**2+dist(4)**2-2*dist(5)**2+dist(6)**2)))/(2*Masseff*sqrt(abs(temp)))
!          coord_cart(3,2,i) = 0d0
!          coord_cart(3,3,i) = 0d0
!          coord_cart(4,1,i) =(-M(2)**2*(dist(1)**2-dist(3)**2+dist(5)**2)+M(1)*(2*M(4)*dist(3)**2+M(2)*(dist(1)**2+dist(3)**2-dist(5)**2)+M(3)*(dist(2)**2+dist(3)**2-dist(6)**2))+M(3)*(-M(3)*(dist(2)**2-dist(3)**2+dist(6)**2)+M(4)*(-dist(2)**2+dist(3)**2+dist(6)**2))-M(2)*(M(4)*(dist(1)**2-dist(3)**2-dist(5)**2)+M(3)*(dist(1)**2+dist(2)**2-2*dist(3)**2-2*dist(4)**2+dist(5)**2+dist(6)**2)))/(2*Masseff*sqrt(abs(temp))) 
!          coord_cart(4,2,i) = 0d0
!          coord_cart(4,3,i) = 0d0
!        else
!          coord_cart(3,1,i) =(-M(2)**2*(dist(1)**2-dist(2)**2+dist(4)**2)+M(1)*(2*M(3)*dist(2)**2+M(2)*(dist(1)**2+dist(2)**2-dist(4)**2)+M(4)*(dist(2)**2+dist(3)**2-dist(6)**2))+M(4)*(M(4)*(dist(2)**2-dist(3)**2-dist(6)**2)+M(3)*(dist(2)**2-dist(3)**2+dist(6)**2))-M(2)*(M(3)*(dist(1)**2-dist(2)**2-dist(4)**2)+M(4)*(dist(1)**2-2*dist(2)**2+dist(3)**2+dist(4)**2-2*dist(5)**2+dist(6)**2)))/(2*Masseff*sqrt(abs(temp)))
!          coord_cart(3,2,i) =-(M(2)*(M(3)*(dist(1)**4 + (dist(2)**2 - dist(4)**2)**2 - 2*dist(1)**2*(dist(2)**2 + dist(4)**2)) + M(4)*(dist(1)**4 + (dist(2)**2 - dist(4)**2)*(dist(3)**2 - dist(5)**2)-dist(1)**2*(dist(2)**2 + dist(3)**2 + dist(4)**2 + dist(5)**2 - 2*dist(6)**2))) + M(4)*(M(4)*(-dist(3)**4 - 2*dist(3)**2*dist(4)**2 + dist(3)**2*dist(5)**2 + dist(2)**2*(dist(3)**2 + dist(5)**2) + dist(3)**2*dist(6)**2 - dist(5)**2*dist(6)**2 + dist(1)**2*(-dist(2)**2 + dist(3)**2 + dist(6)**2)) - M(3)*(-dist(2)**4 + dist(4)**2*(dist(3)**2 - dist(6)**2) + dist(1)**2*(dist(2)**2 - dist(3)**2 + dist(6)**2) + dist(2)**2*(dist(3)**2 + dist(4)**2 - 2*dist(5)**2 + dist(6)**2))))/(2*sqrt(abs(temp*temp2)))
!          coord_cart(3,3,i) = - M(4)*sqrt(abs(dist(2)**4*dist(5)**2 + dist(1)**4*dist(6)**2 + dist(4)**2*(dist(3)**4 + dist(5)**2*dist(6)**2 + dist(3)**2*(dist(4)**2 - dist(5)**2 - dist(6)**2)) + dist(1)**2*(dist(2)**2*(dist(4)**2 - dist(5)**2 - dist(6)**2) + dist(6)**2*(-dist(4)**2 - dist(5)**2 + dist(6)**2) - dist(3)**2*(dist(4)**2 - dist(5)**2 + dist(6)**2)) - dist(2)**2*(dist(3)**2*(dist(4)**2 + dist(5)**2 - dist(6)**2) + dist(5)**2*(dist(4)**2 - dist(5)**2 + dist(6)**2))))/(sqrt(abs(temp2)))
!          coord_cart(4,1,i) =(-M(2)**2*(dist(1)**2-dist(3)**2+dist(5)**2)+M(1)*(2*M(4)*dist(3)**2+M(2)*(dist(1)**2+dist(3)**2-dist(5)**2)+M(3)*(dist(2)**2+dist(3)**2-dist(6)**2))+M(3)*(-M(3)*(dist(2)**2-dist(3)**2+dist(6)**2)+M(4)*(-dist(2)**2+dist(3)**2+dist(6)**2))-M(2)*(M(4)*(dist(1)**2-dist(3)**2-dist(5)**2)+M(3)*(dist(1)**2+dist(2)**2-2*dist(3)**2-2*dist(4)**2+dist(5)**2+dist(6)**2)))/(2*Masseff*sqrt(abs(temp))) 
!          coord_cart(4,2,i) = -(M(2)*(M(4)*(dist(1)**4 + (dist(3)**2 - dist(5)**2)**2 - 2*dist(1)**2*(dist(3)**2 + dist(5)**2)) + M(3)*(dist(1)**4 + (dist(2)**2 - dist(4)**2)*(dist(3)**2 - dist(5)**2) - dist(1)**2*(dist(2)**2 + dist(3)**2 + dist(4)**2 + dist(5)**2 - 2*dist(6)**2))) + M(3)*(M(4)*(dist(3)**4 + 2*dist(3)**2*dist(4)**2 - dist(3)**2*dist(5)**2 - dist(2)**2*(dist(3)**2 + dist(5)**2) - dist(3)**2*dist(6)**2 + dist(5)**2*dist(6)**2 + dist(1)**2*(dist(2)**2 - dist(3)**2 - dist(6)**2)) + M(3)*(-dist(2)**4 + dist(4)**2*(dist(3)**2 - dist(6)**2) + dist(1)**2*(dist(2)**2 - dist(3)**2 + dist(6)**2) + dist(2)**2*(dist(3)**2 + dist(4)**2 - 2*dist(5)**2 + dist(6)**2)))) /(2*sqrt(abs(temp*temp2)))
!          coord_cart(4,3,i) = -M(3)*coord_cart(3,axis(3),i)/M(4)
!        endif       
        
        if(group=='DNH'.or.group=='dnh')then 
          ! apply a rotation to relocate atom A on the x axis
          theta1=atan(-coord_cart(1,2,i)/coord_cart(1,1,i))
          temp=coord_cart(1,1,i)
          temp2=coord_cart(1,2,i)
          coord_cart(1,1,i) = temp*cos(theta1)-temp2*sin(theta1)
          coord_cart(1,2,i) = temp*sin(theta1)+temp2*cos(theta1)
          temp=coord_cart(2,1,i)
          temp2=coord_cart(2,2,i)          
          coord_cart(2,1,i) = temp*cos(theta1)-temp2*sin(theta1)
          coord_cart(2,2,i) = temp*sin(theta1)+temp2*cos(theta1)
          temp=coord_cart(3,1,i)
          temp2=coord_cart(3,2,i)          
          coord_cart(3,1,i) = temp*cos(theta1)-temp2*sin(theta1)
          coord_cart(3,2,i) = temp*sin(theta1)+temp2*cos(theta1)
          temp=coord_cart(4,1,i)
          temp2=coord_cart(4,2,i)          
          coord_cart(4,1,i) = temp*cos(theta1)-temp2*sin(theta1)
          coord_cart(4,2,i) = temp*sin(theta1)+temp2*cos(theta1)          
        endif   
        
     case DEFAULT
       write(*,*) 'invalid number of center'
       stop
     end select  
     write( fmt_output, '(a6,i0,a11)' ) '(i6,2x',2*ncoord,'(3x,f15.4))'
     write(io_output,fmt_output) i, (coord(j, i ),j=1,ncoord), (coord_dist(j, i ),j=1,ncoord)
   enddo
   close(io_output)
   deallocate(temp_coord)

  return
end subroutine  
