subroutine makedata_pshf(opt,workdir,f_pshf,f_temp,f_temp1,prefix,ncore,nelac,mcharg,iz,group,&
	&l_efield,l_efield_1,efield,trec,rang)
! modification d'un fichier f_pshf initial, pour chaque valeur de grille
! -----ATTENTION a la mise en page de f_pshf initial
! voir fichier modele: modele_tests
! cree le fichier f_temp (=f_pshf, avec valeur de grille modifiee)
! cree le fichier f_temp1 (donnees pour le 2eme passage pshf)
! par examen des differentes NAMELIST dans f_temp
! -----ATTENTION: le nom de chaque NAMELIST, et le END associe
!                 doivent etre ecrits SEULS sur une ligne
!                 les variables lues doivent figurer dans les lignes
!                 situees entre &NOM et &END
! -----ATTENTION: insertion dans f_temp1 d'une ligne blanche entre les
!                 NAMELIST, pour parer a des CLOSE ou REWIND manquants
!                 dans pshf.f
!
  use io_unit 
  use grid_data
  implicit none
  integer :: lwork,i,k,il,iz,j
  integer :: number_cent,ncore,rang,ncent,nelac,mcharg
  real(dp) :: trec
  character(len=180) :: workdir,cipsidir,f_pshf,f_temp,f_temp1
  character(len=30) :: coord_str,charg_str,fmt
  character(len=80) :: ligne
  character(len=*) :: group,opt
  character(len=10) :: prefix
  character(len=200) ::	command
  character(len=200),dimension(:),allocatable :: text
  logical :: l_efield,l_efield_1
  logical :: in_namelist
  real(dp),dimension(:) :: efield(3)
  real(dp),dimension(:) :: calfa(10)

! namelists for rcut and rpol
  integer plder,plmax,dc
  logical ysoft		
  parameter (plder=3,dc=500)	! these parameters (esp. 'dc') depend on the configuration pshf.prm during compile time of rcut
  real(dp),dimension(:,:) :: rcut1(dc,0:plder)
  real(dp),dimension(:) :: temp(dc)
  integer,dimension(:) :: no(2)
  namelist/rcutval/plmax,rcut1,ysoft,no

  real(dp),dimension(:) :: rcut(dc)
  namelist/rcval/rcut,ysoft
  
  logical :: ymono
  namelist/vpol/calfa,ymono

  call get_environment_variable("CIPSI_ROOT", cipsidir)  !repertoire de cipsi
  cipsidir=trim(cipsidir)
  lwork=len_trim(workdir)

  open(io_input+rang*100,file=workdir(1:lwork)//'/'//f_pshf)
  open(io_namelist+rang*100,file=workdir(1:lwork)//'/'//f_temp)
!   write( 6, * ) '*******************************************************'
!   write( 6, * ) 'Si plantage a la lecture des donnees, penser a verifier'
!   write( 6, * ) 'la mise en page du fichier .dat:'
!   write( 6, * ) '- nom du NAMELIST et &END sur des lignes separees'
!   write( 6, * ) '- les ncore (ou ncore-1 pour DINFH) NAMELIST &CENT a la suite'
!   write( 6, * ) '*******************************************************'

  rcut=0d0
  rcut1=0d0
  read(io_input+rang*100,rcutval)	! read namelist for rcut 
  read(io_input+rang*100,rcval)	! read namelist for rpol
  read(io_input+rang*100,vpol)	! read namelist for calfa  
  rewind (io_input+rang*100)
  
  if(group=='CNH'.or.group=='cnh'.or.group.eq.'DNH'.or.group.eq.'dnh'.or.group.eq.'DINFH'.or.group.eq.'dinfh'.or.group.eq.'CNV'.or.group.eq.'cnv') then ! inversion center, reorder namelist data
    temp=0d0
    temp(1:ncore)=calfa(1:ncore)
    do i=1,ncore/2
      calfa(2*i-1)=temp(i)
      calfa(2*i)=calfa(2*i-1)
    enddo  
    temp(1:ncore)=rcut(1:ncore)
    do i=1,ncore/2
      rcut(2*i-1)=temp(i)
      rcut(2*i)=rcut(2*i-1)
    enddo 
    do j=0,plder
      temp(1:ncore)=rcut1(1:ncore,j)
      do i=1,ncore/2
        rcut1(2*i-1,j)=temp(i)
        rcut1(2*i,j)=rcut1(2*i-1,j)
      enddo 
    enddo
  endif
  
  il=0
  do                                      !nbre de lignes dans F_PSHF
    read(io_input+rang*100,'(a200)',end=1)
    il=il+1
  enddo
1 rewind (io_input+rang*100)
  allocate(text(il))
  do i=1,il                              !lecture des lignes
    read(io_input+rang*100,'(a200)')text(i)
    write(io_namelist+rang*100,'(a200)')text(i)
  enddo
  close(io_input+rang*100)

  ! Add cartesian coordinate of each center 
  ncent=ncore
  if(group=='DINFH'.or.group=='dinfh'.or.group=='CNH'.or.group=='cnh'.or.group=='DNH'.or.group=='dnh'.or.group.eq.'CNV'.or.group.eq.'cnv')ncent=ncent/2
  do k=1,ncent
    rewind (io_namelist+rang*100)
    deallocate(text)
    il=0
    do                                            !nbre de lignes dans F_PSHF
      read(io_namelist+rang*100,'(a200)',end=2)
      il=il+1
    enddo
2   rewind (io_namelist+rang*100)
    allocate(text(il))
    do i=1,il                         	!read lines again
      read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
    enddo
    rewind (io_namelist+rang*100)				!go back and rewrite complete file

    in_namelist = .false.
    number_cent=0
    do i=1,il                         	!look for cent namelist
      if(.not.in_namelist) then
         if((text(i)(1:5)=='&cent').or.(text(i)(1:5)=='&CENT')) then
           number_cent=number_cent+1
           if(number_cent==k) in_namelist = .true.
         endif
         write(io_namelist+rang*100,'(a200)')text(i)
      else				!if found add GROUP variable at the end 
         if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
           in_namelist = .false.
           write(coord_str,'(f0.8)')coord_cart(k,1,iz)
	   write(io_namelist+rang*100,'(a40)')'X='//coord_str//','
           write(coord_str,'(f0.8)')coord_cart(k,2,iz)
	   write(io_namelist+rang*100,'(a40)')'Y='//coord_str//','
           write(coord_str,'(f0.8)')coord_cart(k,3,iz)
	   write(io_namelist+rang*100,'(a40)')'Z='//coord_str//','	   
           write(io_namelist+rang*100,'(a200)')text(i)
           if(number_cent==ncent) then
             write(io_namelist+rang*100,'(a5)')'&CENT' 
             write(io_namelist+rang*100,'(a4)')'&END' 
           endif
         else
           write(io_namelist+rang*100,'(a200)')text(i)
         endif
      endif
    enddo
  enddo
  rewind (io_namelist+rang*100)  

! JD 06/06 Write 'group' into hondo-namelist  (required for automatic change of sym group 
!		in the chase of an external electric field)

!  change SCFINP
  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist+rang*100,'(a200)',end=3)
    il=il+1
  enddo
3 rewind (io_namelist+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist+rang*100)				!go back and rewrite complete file

  in_namelist = .false.
  do i=1,il                         	!look for hondo namelist
    if(.not.in_namelist) then
       if((text(i)(1:7)=='&scfinp').or.(text(i)(1:7)=='&SCFINP')) then
         in_namelist = .true.
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    else				!if found add GROUP variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
	 select case(nelac)
	 case (1)
	   write(io_namelist+rang*100,*)' QN=1.,'	 
	 case (2)
	   write(io_namelist+rang*100,*)' QN=2.,'	 
	 case (3)
	   write(io_namelist+rang*100,*)' QN=2.,1.,'	 
	 case (4)
	   write(io_namelist+rang*100,*)' QN=2.,2.,'	 
	 case (5)
	   write(io_namelist+rang*100,*)' QN=2.,2.,1.,'	 
	 case (6)
	   write(io_namelist+rang*100,*)' QN=2.,2.,2.,'	 
	 case DEFAULT
	   write(*,*) 'too many electron'
	   stop
	 end select   
         write(io_namelist+rang*100,'(a200)')text(i)
       else
         write(io_namelist+rang*100,'(a200)')text(i)
       endif
    endif
  enddo
  rewind (io_namelist+rang*100)
  
! change HONDO  
  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist+rang*100,'(a200)',end=4)
    il=il+1
  enddo
4 rewind (io_namelist+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist+rang*100)				!go back and rewrite complete file

  in_namelist = .false.
  do i=1,il                         	!look for hondo namelist
    if(.not.in_namelist) then
       if((text(i)(1:6)=='&hondo').or.(text(i)(1:6)=='&HONDO')) then
         in_namelist = .true.
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    else				!if found add GROUP variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
         
	 if(group=='CNV'.or.group=='cnv')write(io_namelist+rang*100,*)' NAXIS=2'
	 if(group=='CNH'.or.group=='cnh')write(io_namelist+rang*100,*)' NAXIS=2'
	 if(group=='DNH'.or.group=='dnh')write(io_namelist+rang*100,*)' NAXIS=2'	 
         write(charg_str,'(i0)')mcharg	 
	 write(io_namelist+rang*100,*)' MCHARG='//charg_str//','
	 if(mod(nelac,2)==0) then
           write(io_namelist+rang*100,*)' SZ=0.,'
         else
           write(io_namelist+rang*100,*)' SZ=0.5,'
         endif  
	 write(io_namelist+rang*100,*)' GROUP='//group//','        
         write(io_namelist+rang*100,'(a200)')text(i)
       else
         write(io_namelist+rang*100,'(a200)')text(i)
       endif
    endif
  enddo
  rewind (io_namelist+rang*100)
  
  ! namelist &IJKLIN
  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist+rang*100,'(a200)',end=5)
    il=il+1
  enddo
5 rewind (io_namelist+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist+rang*100)				!go back and rewrite complete file
  
  in_namelist = .false.
  do i=1,il                         	!look for ijklin namelist
    if(.not.in_namelist) then
       if((text(i)(1:7)=='&ijklin').or.(text(i)(1:7)=='&IJKLIN')) then
         in_namelist = .true.
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    else				!if found add GROUP variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
         write(io_namelist+rang*100,*)' TREC=',trec,','
         write(io_namelist+rang*100,*)' NREC= 0,'         
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    endif
  enddo
  rewind (io_namelist+rang*100)
  
  ! namelist &VPOL
  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist+rang*100,'(a200)',end=6)
    il=il+1
  enddo
6 rewind (io_namelist+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist+rang*100)				!go back and rewrite complete file
  in_namelist = .false.
  do i=1,il                         	!look for vpol namelist
    if(.not.in_namelist) then
       if((text(i)(1:5)=='&vpol').or.(text(i)(1:5)=='&VPOL')) then
         in_namelist = .true.
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    else				!if found add MREC variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
         write(io_namelist+rang*100,*)' MREC= 0,'         
         write(fmt,'(i0)')ncore	
         write(io_namelist+rang*100,'(a7,'//fmt//'f20.12,a1)')' CALFA=',(calfa(j),j=1,ncore),','         
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    endif
  enddo
  rewind (io_namelist+rang*100)

! JD 04/06 Add electrical field to file temp.dat 
  if(l_efield.and.l_efield_1)then	!if l_efield is true and temp.dat should be 
					!modified, add line with field-components 
					
    deallocate(text)		! read 60 again
    rewind (io_namelist+rang*100)
    il=0
    do                                            !nbre de lignes dans F_PSHF
      read(io_namelist+rang*100,'(a200)',end=8)
      il=il+1
    enddo
8   rewind (io_namelist+rang*100)
    allocate(text(il))
    do i=1,il                         	!read lines again
      read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
    enddo
    rewind (io_namelist+rang*100)  					
    in_namelist = .false.
    do i=1,il                         	!look for hondo namelist
      if(.not.in_namelist) then
         if((text(i)(1:6)=='&hondo').or.(text(i)(1:6)=='&HONDO')) then
           in_namelist = .true.
         endif
         write(io_namelist+rang*100,'(a200)')text(i)
      else				!if found add ELEC variable at the end 
         if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
           in_namelist = .false.
           write(io_namelist+rang*100,'(a6,f14.10,a1,f14.10,a1,f14.10,a1)')' ELEC=',efield(1),',',efield(2),',',efield(3),','
	   write(io_namelist+rang*100,'(a200)')text(i)
         else
           write(io_namelist+rang*100,'(a200)')text(i)
         endif
      endif
    enddo
    rewind (io_namelist+rang*100)   
  endif
  
! *************************************************  
! add prefix tag to the *fil namelists 
  
  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist+rang*100,'(a200)',end=9)
    il=il+1
  enddo
9 rewind (io_namelist+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist+rang*100)				!go back and rewrite complete file
  in_namelist = .false.
  do i=1,il                         
    if(.not.in_namelist) then
       if((text(i)(1:7)=='&rpofil').or.(text(i)(1:7)=='&RPOFIL')) then
         in_namelist = .true.
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    else				!if found add prefix variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
         write(io_namelist+rang*100,'(a18)')"prefix='"//trim(prefix)//"',"             
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    endif
  enddo
  rewind (io_namelist+rang*100)  

  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist+rang*100,'(a200)',end=10)
    il=il+1
  enddo
10 rewind (io_namelist+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist+rang*100)				!go back and rewrite complete file
  in_namelist = .false.
  do i=1,il                         	
    if(.not.in_namelist) then
       if((text(i)(1:7)=='&rcufil').or.(text(i)(1:7)=='&RCUFIL')) then
         in_namelist = .true.
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    else				!if found add prefix variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
         write(io_namelist+rang*100,'(a18)')"prefix='"//trim(prefix)//"',"             
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    endif
  enddo
  rewind (io_namelist+rang*100)  
  
  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist+rang*100,'(a200)',end=11)
    il=il+1
  enddo
11 rewind (io_namelist+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist+rang*100)				!go back and rewrite complete file
  in_namelist = .false.
  do i=1,il                         	
    if(.not.in_namelist) then
       if((text(i)(1:7)=='&cvafil').or.(text(i)(1:7)=='&CVAFIL')) then
         in_namelist = .true.
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    else				!if found add prefix variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
         write(io_namelist+rang*100,'(a18)')"prefix='"//trim(prefix)//"',"             
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    endif
  enddo
  rewind (io_namelist+rang*100)    
  
  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist+rang*100,'(a200)',end=12)
    il=il+1
  enddo
12 rewind (io_namelist+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist+rang*100)				!go back and rewrite complete file
  in_namelist = .false.
  do i=1,il                         	
    if(.not.in_namelist) then
       if((text(i)(1:6)=='&ijfil').or.(text(i)(1:6)=='&IJFIL')) then
         in_namelist = .true.
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    else				!if found add prefix variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
         write(io_namelist+rang*100,'(a18)')"prefix='"//trim(prefix)//"',"             
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    endif
  enddo
  rewind (io_namelist+rang*100)    
  
  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist+rang*100,'(a200)',end=13)
    il=il+1
  enddo
13 rewind (io_namelist+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist+rang*100)				!go back and rewrite complete file
  in_namelist = .false.
  do i=1,il                         	
    if(.not.in_namelist) then
       if((text(i)(1:7)=='&fokfil').or.(text(i)(1:7)=='&FOKFIL')) then
         in_namelist = .true.
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    else				!if found add prefix variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
         write(io_namelist+rang*100,'(a18)')"prefix='"//trim(prefix)//"',"             
       endif
       write(io_namelist+rang*100,'(a200)')text(i)
    endif
  enddo
  rewind (io_namelist+rang*100)    
  
  
! JD 02/09 If new chain with rpol+rcut is used, prepare namelists automatically
  if((opt=='rpol').or.(opt=='RPOL'))then
    deallocate(text)			! read again
    il=0
    do                                            !nbre de lignes dans F_PSHF
      read(io_namelist+rang*100,'(a200)',end=7)
      il=il+1
    enddo
7   rewind (io_namelist+rang*100)
    allocate(text(il))
    do i=1,il                         	!read lines again
      read(io_namelist+rang*100,'(a200)')text(i)		!store file temporarily in memory
    enddo
    rewind (io_namelist+rang*100)				!go back and rewrite complete file
    in_namelist = .false.
    do i=1,il                         	!look for rcutval and rcval namelist
      if(.not.in_namelist) then
	if((text(i)(1:8)=='&rcutval').or.(text(i)(1:8)=='&RCUTVAL').or.(text(i)(1:6)=='&rcval').or.(text(i)(1:6)=='&RCVAL')) then
	  in_namelist = .true.
	else
	  write(io_namelist+rang*100,'(a200)')text(i)
	endif
      else 	! if found, erase namelist
	if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
	  in_namelist = .false.
	endif
      endif
    enddo
    
    do i=1,ncore
      rcut1(i,0)=-abs(rcut1(i,0))
      rcut(i)=abs(rcut1(i,0))
    enddo

    write(io_namelist+rang*100,rcutval)
    write(io_namelist+rang*100,rcval)
    
    rewind (io_namelist+rang*100)    
  endif

  
  deallocate(text)
  
! ***************************** Fichier temp1.dat
  ! JD 04/06
  open(io_namelist2+rang*100,file=workdir(1:lwork)//'/'//f_temp1)!creation fichier f_temp1
                                                !pour pshf+cv
  ! Copie complete
  rewind (io_namelist+rang*100)  
  il=0
  do                                      !nbre de lignes dans F_PSHF
    read(io_namelist+rang*100,'(a200)',end=14)
    il=il+1
  enddo
14 rewind (io_namelist+rang*100)
  allocate(text(il))
  do i=1,il                              !lecture des lignes
    read(io_namelist+rang*100,'(a200)')text(i)
    write(io_namelist2+rang*100,'(a200)')text(i)
  enddo
  close(io_namelist+rang*100)
  rewind (io_namelist2+rang*100)
  
  !namelist &PSFIL
  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist2+rang*100,'(a200)',end=15)
    il=il+1
  enddo
15 rewind (io_namelist2+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist2+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist2+rang*100)				!go back and rewrite complete file
  in_namelist = .false.
  do i=1,il                         	!look for vpol namelist
    if(.not.in_namelist) then
       if((text(i)(1:6)=='&psfil').or.(text(i)(1:6)=='&PSFIL')) then
         in_namelist = .true.
       endif
       write(io_namelist2+rang*100,'(a200)')text(i)
    else				!if found add pqrs and info variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
         write(io_namelist2+rang*100,'(a16)')' pqrs='//"'pqrs_cv'"//','
         write(io_namelist2+rang*100,'(a16)')' info='//"'info_cv'"//','        
       endif
       write(io_namelist2+rang*100,'(a200)')text(i)
    endif
  enddo
  rewind (io_namelist2+rang*100)
  
  !namelist &HONDO
  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist2+rang*100,'(a200)',end=16)
    il=il+1
  enddo
16 rewind (io_namelist2+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist2+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist2+rang*100)				!go back and rewrite complete file
  in_namelist = .false.
  do i=1,il                         	!look for vpol namelist
    if(.not.in_namelist) then
       if((text(i)(1:6)=='&hondo').or.(text(i)(1:6)=='&HONDO')) then
         in_namelist = .true.
       endif
       write(io_namelist2+rang*100,'(a200)')text(i)
    else				!if found add pqrs and info variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
         if(l_efield.and.(.not.l_efield_1))then
  	   write(io_namelist2+rang*100,'(a6,f14.10,a1,f14.10,a1,f14.10,a1)')' ELEC=',efield(1),',',efield(2),',',efield(3),','
         endif
         write(io_namelist2+rang*100,'(a9)')' IREST=7,'         !insertion automatique en derniere ligne      
       endif
       write(io_namelist2+rang*100,'(a200)')text(i)
    endif
  enddo
  rewind (io_namelist2+rang*100)  
  
  !namelist &SCFINP
  deallocate(text)
  il=0
  do                                            !nbre de lignes dans F_PSHF
    read(io_namelist2+rang*100,'(a200)',end=17)
    il=il+1
  enddo
17 rewind (io_namelist2+rang*100)
  allocate(text(il))
  do i=1,il                         	!read lines again
    read(io_namelist2+rang*100,'(a200)')text(i)		!store file temporarily in memory
  enddo
  rewind (io_namelist2+rang*100)				!go back and rewrite complete file
  in_namelist = .false.
  do i=1,il                         	!look for vpol namelist
    if(.not.in_namelist) then
       if((text(i)(1:7)=='&scfinp').or.(text(i)(1:7)=='&SCFINP')) then
         in_namelist = .true.
       endif
       write(io_namelist2+rang*100,'(a200)')text(i)
    else				!if found add pqrs and info variable at the end 
       if((text(i)(1:4)=='&end').or.(text(i)(1:4)=='&END')) then   
         in_namelist = .false.
         write(io_namelist2+rang*100,'(a9)')' MAXIT=0,'          
       endif
       write(io_namelist2+rang*100,'(a200)')text(i)
    endif
  enddo
  rewind (io_namelist2+rang*100)  
  close(io_namelist2+rang*100)  
 
  ! ajout de la racine de cipsi au chemin des repertoire
  command= 'sed -i "s# psnl_fil=# psnl_fil='''//trim(cipsidir)//'#" '//workdir(1:lwork)//'/'//f_temp
  call system(command)
  command= 'sed -i "s# ps_molcas_fil=# ps_molcas_fil='''//trim(cipsidir)//'#" '//workdir(1:lwork)//'/'//f_temp
  call system(command)
  command= 'sed -i "s# f10=# f10='''//trim(cipsidir)//'#" '//workdir(1:lwork)//'/'//f_temp
  call system(command)
  command= 'sed -i "s# psnl_fil=# psnl_fil='''//trim(cipsidir)//'#" '//workdir(1:lwork)//'/'//f_temp1
  call system(command)
  command= 'sed -i "s# ps_molcas_fil=# ps_molcas_fil='''//trim(cipsidir)//'#" '//workdir(1:lwork)//'/'//f_temp1
  call system(command)
  
  return
end