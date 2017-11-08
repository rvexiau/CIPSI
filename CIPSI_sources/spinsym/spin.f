      program spinsym
      implicit real*8(a-h,o-p,r-z),logical*4(q)
      character*26 timeday
      include 'pshf.prm'   
      logical*1 qion,yocd,yspin
      integer*4 ne,trou,part,iorb,ispin,nd
      character*8 group  
      dimension fmul(metz),c(metz*ndetz),e(ndetz),itmp(ndetz)
      integer,dimension(:),allocatable :: spin,norme   
      real*8,dimension(:),allocatable :: sol    
      integer :: dettab(6)            
      dimension itsyp(doa)
      namelist/sort/group,nusym,nelac,seuil
      common/spob/isz,noca,nvr,qion
      common/spo/ne(ndetz+1),nd(ndetz+1),
     1trou(8*ndetz),part(8*ndetz),
     2iorb(2*doa),ispin(2*doa),yocd(2*doa,ndetz)      
  
      call fdate(timeday)
      write(6,9999) timeday
 9999 format(/,80(1h*),/,8(10h   spin     ),/,  10x,a25,/,
     * 80(1h*))
       
      seuil=1d-8 
      read(5,sort)      
      rewind(5)
      call openf   
      read(70) ncf,nelac,nsmax,yspin
c     lecture des vecteurs propres
      if(yspin) then
      read(59)ncf,nvr,(c(k),k=1,ncf*nvr),(e(k),k=1,nvr),
     * (itmp(k),k=1,nvr)
      else
      read(59)ncf,nvr,(c(k),k=1,ncf*nvr),(e(k),k=1,nvr)
      endif
      allocate(spin(nvr),sol(nvr))      
c    
      if(.not.(yspin)) then
c     lecture des determinants 
      yocd=.true.
      read(64)ndetr,ii,mnorb,qion,(ne(i),nd(i),i=1,ndetr),
     *(trou(i),part(i),i=1,ii),
     *(iorb(i),i=1,mnorb+1),(ispin(i),i=1,mnorb+1)
     
      if(ndetr.lt.ncf) then
        write(6,*) 'error : ndetr too low'
        stop
      endif
      write(6,*) 'nombre de determinants :',ncf
      
      call mult(fmul,ncf,c)  
      do i=1,nvr
        del=1.+4.*fmul(i)                                                
        sol(i)=(-1.+dsqrt(del))/2.        
        spin(i)=1+int(0.5+2.*sol(i))    
      enddo
      
      else
        do i=1,nvr
        spin(i)=itmp(i)
        sol(i)=-5.
        enddo
      endif
      
      
      write(6,*) 'multiplicite de spin'
      write(6,*) 'metat     2S+1           e_var'
      do i=1,nvr
         write(6,37) i,spin(i),e(i)
37       format(i3,9x,i3,7x,f20.12)
      enddo
      
      write(6,*)
      write(6,*)
      write(6,*)
      
       if(group.eq.'CINFV'.or.group.eq.'cinfv'.or.group.eq.'DINFH'
     *     .or.group.eq.'dinfh') then
     
      write(*,'(a39)') 'Contribution of monoexcited determinant'
      write(*,'(3x,4x,10x,a5,10x,12x,a2,11x,10x,a5,10x,10x,a3)')
     *    'Sigma','Pi','Delta','Phi'
        rewind(64)
       read(64)ndetr,ii,mnorb,qion,(ne(i),nd(i),i=1,ndetr),
     *(trou(i),part(i),i=1,ii),
     *(iorb(i),i=1,mnorb+1),(ispin(i),i=1,mnorb+1)
      read(40) nsym,norb,noc,ntrans,                  
     1 (itsyp(i),i=1,norb)  
      do i=norb+1,2*norb
      itsyp(i)=itsyp(i-norb)
      enddo
      itsyp(mnorb+1)=1
      if(qion) then
        nelac=nelac+1
      endif  
      
      allocate(norme(ncf))
      norme=nelac
      
      do i=1,ncf
      
      ! construct reference determinant
      select case (nelac/2)
      case (1)
        dettab=(/1,norb+1,0,0,0,0/)
      case (2)
        dettab=(/1,norb+1,2,norb+2,0,0/)
      case (3)
        dettab=(/1,norb+1,2,norb+2,3,norb+3/)
      case DEFAULT 
          write(*,*) 'wrong number of electron' 
          stop
      end select
        
       ! construct virtual determinant by removing trou from ref
      do jj=1,ne(i)
        nt=trou(nd(i)+jj)
        np=part(nd(i)+jj)
        do k=1,nelac
          if(nt.eq.dettab(k)) then
            dettab(k)=np
            exit
          endif  
        enddo        
      enddo
        
      ! check number of sigma orbitals in the virtual determinant   
      do jj=1,nelac
        n1=itsyp(dettab(jj))
        if(n1.eq.1.or.
     *   ((group.eq.'DINFH'.or.group.eq.'dinfh')
     *     .and.n1.eq.2)) then
        norme(i)=norme(i)-1    
        else
          nsymdet=n1
        endif  
      enddo
      
      ! compute symmetry
      select case(norme(i))
      case(0)
        norme(i)=0
      case(1)
        if(group.eq.'CINFV'.or.group.eq.'cinfv') then
          select case(nsymdet)
          case(1)
            norme(i)=0
          case(2,3)
            norme(i)=1
          case(4,5)
            norme(i)=2
          case(6,7)
            norme(i)=3
          case DEFAULT
            norme(i)=-1
          end select
        else
          select case(nsymdet)
          case(1,2)
            norme(i)=0
          case(3,4,5,6)
            norme(i)=1
          case(7,8,9,10)
            norme(i)=2
          case(11,12,13,14)
            norme(i)=3
          case DEFAULT
            norme(i)=-1
          end select    
        endif
      case DEFAULT
        norme(i)=-1
      end select  
      
      ! end of determinant sorting
      enddo
      
      
      do j=1,nvr
      dnormes=0d0
      dnormep=0d0
      dnormed=0d0
      dnormef=0d0
      qgood=.true.
      
      do i=1,ncf             
        select case(norme(i))
        case(-1)
          cycle
        case(0)
          dnormes=dnormes+c(ncf*(j-1)+i)*c(ncf*(j-1)+i)
        case(1)
          dnormep=dnormep+c(ncf*(j-1)+i)*c(ncf*(j-1)+i)
        case(2)
          dnormed=dnormed+c(ncf*(j-1)+i)*c(ncf*(j-1)+i)
        case(3)
          dnormef=dnormef+c(ncf*(j-1)+i)*c(ncf*(j-1)+i)
        case DEFAULT
          cycle
        end select    
      enddo
      
      if(group.eq.'CINFV'.or.group.eq.'cinfv') then

      select case(nusym)
      case (1)
        dnorme=dnormep+dnormed+dnormef
        if(dnorme.gt.dnormes) qgood=.false.
      case (2,3)
        dnorme=dnormes+dnormed+dnormef
        if(dnorme.gt.dnormep) qgood=.false.      
      case (4,5)
        dnorme=dnormes+dnormep+dnormef
        if(dnorme.gt.dnormed) qgood=.false.      
      case (6,7)
        dnorme=dnormes+dnormep+dnormed
        if(dnorme.gt.dnormef) qgood=.false.      
      end select
      elseif(group.eq.'DINFH'.or.group.eq.'dinfh') then
      select case(nusym)
      case (1,2)
        dnorme=dnormep+dnormed+dnormef
        if(dnorme.gt.dnormes) qgood=.false.
      case (3,4,5,6)
        dnorme=dnormes+dnormed+dnormef
        if(dnorme.gt.dnormep) qgood=.false.      
      case (7,8,9,10)
        dnorme=dnormes+dnormep+dnormef
        if(dnorme.gt.dnormed) qgood=.false.      
      case (11,12,13,14)
        dnorme=dnormes+dnormep+dnormed
        if(dnorme.gt.dnormef) qgood=.false.      
      end select
      endif
      
      write(*,'(i2,1x,l1,3x,4f25.15)') j,qgood,dnormes,dnormep,dnormed,dnormef
      enddo
      
      deallocate(norme)
      endif
     
      deallocate(spin,sol)
      
      call fdate(timeday)
      write(6,9998) timeday
 9998 format(/,80(1h*),/,8(10h fin spin   ),/, 10x,a25,/,
     * 80(1h*))
      
      stop
99    write(6,*) 'error reading file f04'    
      
      end program