C ***************************************************************
C *  performs wigner rotation of matrix elements stored in tab
C    uses precalculated wigner matrix in matrix  mat(l,i,j)
C    l=la,lb     j=1,2*l+1
C    rmat(i,j)  3x3 cartesian rotation matrix

      subroutine wigrot(la,lb,tab,mat,rmat)
      implicit none

      integer la,lb
      real*8 tab(4,10,10)
      real*8 mat(2,7,7)
      real*8 rmat(3,3)
 
      integer ivec,i,j,k
      real*8 temp1,temp2
      real*8 tac(7,7)

C transform first index

      do ivec=1,4
        do i=1,2*la+1
          do j=1,2*lb+1
            tac(i,j) = mat(1,i,1)*tab(ivec,1,j)
            do k=2,2*la+1
              tac(i,j) = tac(i,j)+mat(1,i,k)*tab(ivec,k,j)
            end do
          end do
        end do

c transform second index
        do i=1,2*la+1
          do j=1,2*lb+1
            tab(ivec,i,j) = mat(2,j,1)*tac(i,1)
            do k=2,2*lb+1
              tab(ivec,i,j) = tab(ivec,i,j)+mat(2,j,k)*tac(i,k)
c           tab(ivec,i,j) = mat(2,1,j)*tac(i,1)
c           do k=2,2*lb+1
c             tab(ivec,i,j) = tab(ivec,i,j)+mat(2,k,j)*tac(i,k)
            end do
          end do
        end do

C exchange p orbitals if necessary: canonical order (z,x,y)->PSHF order (x,y,z)
        if (la.eq.1) then
          do j=1,2*lb+1
            temp1=tab(ivec,1,j)
            tab(ivec,1,j) = tab(ivec,2,j)
            tab(ivec,2,j) = tab(ivec,3,j)
            tab(ivec,3,j) = temp1
          end do 
        end if

        if (lb.eq.1) then
          do i=1,2*la+1
            temp1=tab(ivec,i,1)
            tab(ivec,i,1) = tab(ivec,i,2)
            tab(ivec,i,2) = tab(ivec,i,3)
            tab(ivec,i,3) = temp1
          end do 
        end if

      end do


C back rotation of electrostatic field
      do i=1,2*la+1
         do j=1,2*lb+1
           temp1      = rmat(1,1)*tab(2,i,j)+rmat(2,1)*tab(3,i,j)
     &                 +rmat(3,1)*tab(4,i,j)

           temp2      = rmat(1,2)*tab(2,i,j)+rmat(2,2)*tab(3,i,j)
     &                 +rmat(3,2)*tab(4,i,j)

           tab(4,i,j) = rmat(1,3)*tab(2,i,j)+rmat(2,3)*tab(3,i,j)
     &                 +rmat(3,3)*tab(4,i,j)
           tab(2,i,j) = temp1
           tab(3,i,j) = temp2
         end do
      end do

99    return
      end
