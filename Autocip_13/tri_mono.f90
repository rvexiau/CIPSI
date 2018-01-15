subroutine tri_mono(l_work,workdir,l_max,nvec,ncore,f_pshfcv,nf)
  ! programme de lecture des fichiers de resultats monoelectroniques
  !   'AB.pshfcv'
  ! ou AB est la molecule choisie.
  ! AB=molecule diatomique de type alcalin (ie, etat fondamental sigma)
  ! Execution: tri_mono<tri_mono.in
  ! Le programme lit une liste de fichiers a examiner dans tri_mono.in. Le nom de chaque
  ! fichier (FILE_IN) est specifie avec son chemin relatif a l'endroit de l'appel,
  ! et avec la distance (R_IN)correspondante (pour un test de securite).
  ! Le programme recherche l'endroit a lire a l'aide de la routine FIND_STR.
  ! Le nombre d'energies desirees dans chaque symetrie est precise (NVEC).
  ! Le programme analyse les vecteurs propres pour determiner la symetrie
  ! Fichiers de sortie:
  ! 	- energies SCF et CPP en fonction de R
  !	- energies monoelectroniques tries par symetrie
  ! INPUT:
  ! workdir: repertoire de travail
  ! at1,at2: atomes traites
  ! l_max: lambda max des etats
  ! nvec: nombre de vecteurs propres a traiter
  ! file_in: nom du fichier a traiter (de typeAB.pshfcv)
  ! r_in: valeurs de R considerees
  !
  ! OUTPUT: en fonction de r_in
  ! e_1sr: 1/R
  ! e_scf: energie SCF
  ! e_cpp: core-polarization term
  ! e_tri: energies monoelectroniques (incluant e_1sr et e_cpp)
  !
  ! VARIABLES INTERMEDIAIRES
  ! oa1,oa2: orbitales (de type oa_typ) sur chaque centre
  ! n_shell: nbre de couches; 
  ! sh_typ: nbre de couches de chaque type
  ! n_oa: nbre d'OA
  ! oa_pos: position de la premiere OA de chaque type (sym)
  ! n_basis: nbre d'OA total
  ! n_mono: nbre d'OM dans chaque symetrie
  !
  ! SUBROUTINES:
  ! find_str: recherche d'ne chaine dans un fichier
  ! longchain: longueur d'une chaine (sans blanc)
  ! cex2: conversion d'un entier (<99) en chaine de longueur 2
  ! divise: resultat de la division de deux entiers
  ! version 1: fevrier 2004, O.D.
  use io_unit
  implicit none

  integer, parameter :: dp=kind(1.d0)
  integer :: nsym,nvec,nf,istatus,iline,l_str,n_basis,l_max,n_ug,ns_tri,at_num,l_work,nprint,ncore
  integer :: n_shell(2),n_oa(2)
  integer :: ifile,ie,iv,ib,i,j,il,iug
  integer :: divise,longchain
  integer :: l_om(16),l_oa(4),sh_typ(4,2),oa_pos(4,2),n_mono(4,2)
  integer, dimension(:), allocatable :: centre,l_orb
  integer, dimension(:,:), allocatable :: n_tri
  real(dp) :: r,at_q,x1,y1,z1,x2,y2,z2
  real(dp), dimension(:), allocatable :: e_mono,e_1sr,e_cpp,e_scf,r_in
  real(dp), dimension(:,:), allocatable :: v_mono,int_dip_z,int_dip_x,int_dip_y
  real(dp), dimension(:,:,:,:,:), allocatable :: dip_mono_z,dip_mono_x,dip_mono_y
  real(dp), dimension(:,:,:,:), allocatable :: e_tri
  character(len=80) :: file_in,str_to_find,c,file
  character(len=*) :: f_pshfcv(nf),workdir
  character(len=22) :: fmt
  character(len=2) :: cex2,atom1,atom2,pseudo
  character(len=4) :: sym(16)
  character(len=1) :: oa_typ(4),ug(2)
  character(len=3), dimension(:), allocatable :: atom
  character(len=5), dimension(:), allocatable :: orb
  character(len=1), dimension(:), allocatable :: oa1,oa2
  logical :: exist

  data oa_typ/'s','p','d','f'/
  data l_oa/0,1,2,3/
  data sym/' s  ','x   ','y   ','z   ','dsig','dpix','dpiy','ddxx','ddxy',&
&          'fsig','fpix','fpiy','fdxx','fdxy','ffix','ffiy'/
  data l_om/0,1,1,0,0,1,1,2,2,0,1,1,2,2,3,3/
  data ug/'g','u'/

  n_ug=1
      
  open(io_output,file=workdir(1:l_work)//'/tri_e_mono.out',iostat=istatus)
  write(io_output,*)' ',nf,' fichiers '
  
  allocate(r_in(nf),e_1sr(nf),e_cpp(nf),e_scf(nf))
  allocate(e_tri(nvec,4,2,nf))
  allocate(n_tri(nvec,2))
  allocate(dip_mono_z(nvec,4,nvec,4,nf),dip_mono_y(nvec,4,nvec,4,nf),dip_mono_x(nvec,4,nvec,4,nf))
  nprint=-1  

  do ifile=1,nf                    ! boucle liste fichiers
    file_in=f_pshfcv(ifile)
    open(io_work,file=workdir(1:l_work)//'/'//file_in,iostat=istatus)
    write(io_output,*)'ouverture ',file_in(1:longchain(file_in))

    rewind(io_work)                 ! lecture de la base

    call find_str(io_work,'                     atom pseudo', &
  &               istatus,iline,l_str)
    if(istatus==0) then
      write(*,*) 'Missing # BASIS card in PSHFCV file'
      write(io_output,*)iline,' lignes lues dans ',file_in
      stop
    endif
    do i=1,4           !saute 4 lignes
      read(io_work,*)
    enddo
                              ! lecture oa centre 1
    oa_pos=0
    sh_typ=0
    n_oa=0
    read(io_work,'(20x,a2,4x,a2,9x,i2,6x,f4.1,3(f12.7,3x),i5)') &
 &            atom1,pseudo,at_num,at_q,x1,y1,z1,n_shell(1)
    write(io_output,'(20x,a2,4x,a2,9x,i2,6x,f4.1,3(f12.7,3x),i5)') &
 &            atom1,pseudo,at_num,at_q,x1,y1,z1,n_shell(1)
    if(ifile==1)allocate(oa1(n_shell(1)))
    do i=1,divise(n_shell(1),5)
      read(io_work,'(91x,5(x,a1))')(oa1(j),j=(i-1)*5+1,min(i*5,n_shell(1)))
    enddo
    write(io_output,*)oa1
    do j=1,n_shell(1)
      do i=1,4                  ! calcul nombre total d'oa et de chaque type
        if(oa1(j)==oa_typ(i)) then
          sh_typ(i,1)=sh_typ(i,1)+1
          n_oa(1)=n_oa(1)+2*l_oa(i)+1
        endif
      enddo
    enddo
    if(sh_typ(1,1)>0)oa_pos(1,1)=1  ! position de la premiere OA de chaque type
    do i=2,4
      if(sh_typ(i,1)>0)oa_pos(i,1)=(2*l_oa(i-1)+1)*sh_typ(i-1,1)+oa_pos(i-1,1)+(2*l_oa(i-1))
    enddo
    write(io_output,*)'# OA centre 1 =',n_oa(1)
    write(io_output,*)'# dans chaque symetrie'
    write(io_output,*)(sh_typ(i,1),i=1,4)
    write(io_output,*)'position de la premiere de chaque symetrie'
    write(io_output,*)(oa_pos(i,1),i=1,4)

    if(ncore.eq.2)then
     read(io_work,*)
                              ! lecture oa centre 2
     read(io_work,'(20x,a2,4x,a2,9x,i2,6x,f4.1,3(f12.7,3x),i5)') &
  &       atom2,pseudo,at_num,at_q,x2,y2,z2,n_shell(2)
     write(io_output,'(20x,a2,4x,a2,9x,i2,6x,f4.1,3(f12.7,3x),i5)') &
  &       atom2,pseudo,at_num,at_q,x2,y2,z2,n_shell(2)

     r_in(ifile)=abs(z1)+abs(z2)
     if(ifile==1)allocate(oa2(n_shell(2)))
     do i=1,divise(n_shell(2),5)
       read(io_work,'(91x,5(x,a1))')(oa2(j),j=(i-1)*5+1,min(i*5,n_shell(2)))
     enddo
     write(io_output,*)oa2
     do j=1,n_shell(2)
       do i=1,4                  ! calcul nombre total d'oa et de chaque type
         if(oa2(j)==oa_typ(i))then
           sh_typ(i,2)=sh_typ(i,2)+1
           n_oa(2)=n_oa(2)+2*l_oa(i)+1
         endif
       enddo
     enddo
     if(sh_typ(1,2)>0)oa_pos(1,2)=1  ! position de la premiere OA de chaque type
     do i=2,4
       if(sh_typ(i,2)>0)oa_pos(i,2)=(2*l_oa(i-1)+1)*sh_typ(i-1,2)+oa_pos(i-1,2)+(2*l_oa(i-1))
     enddo
     write(io_output,*)'# OA centre 2 =',n_oa(2)
     write(io_output,*)'# dans chaque symetrie'
     write(io_output,*)(sh_typ(i,2),i=1,4)
     write(io_output,*)'position de la premiere de chaque symetrie'
     write(io_output,*)(oa_pos(i,2),i=1,4)
    endif
!   if(.not.((atom1 /= at1).and.(atom2 /= at2)))then
!      write(6,*)'verifier les atomes fichier ',file_in
!      stop
!    endif

    rewind(io_work)                  !nbre total de fonctions de base
    call find_str(io_work,' total number of basis functions   =', &
  &               istatus,iline,l_str)
    if(istatus==0) then
      write(*,*) 'Missing # BASIS FUNCTIONS card in PSHFCV file'
      write(io_output,*)iline,' lignes lues dans ',file_in
      stop
    endif
    rewind(io_work)
    do il=1,iline-1
      read(io_work,*)
    enddo
    read(io_work,'(A'//cex2(l_str)//',i5)')str_to_find(1:l_str),n_basis
    write(io_output,*)n_basis,' fonctions de base'

    !if(n_basis/=n_oa(1)+n_oa(2))then
    !  write(6,*)'probleme nombre OA, fichier ',file_in
    !  write(6,*)'n_basis=',n_basis,' n-oa(1+2)=',n_oa(1)+n_oa(2)
    !endif

    rewind(io_work)                     !1/R
    call find_str(io_work,' energie nucleaire                        stockee sur la file 33', &
  &               istatus,iline,l_str)
    if(istatus==0) then
      write(*,*) 'Missing NUCLEAR ENERGY card in PSHFCV file'
      stop
    endif
    read(io_work,*)e_1sr(ifile)
    write(io_output,*)'1/R=',e_1sr(ifile)

    rewind(io_work)                      !energie SCF
    call find_str(io_work,'  energie SCF (couches fermees)           stockee sur la file 33', &
  &               istatus,iline,l_str)
    if(istatus==0) then
      write(*,*) 'Missing E SCF card in PSHFCV file'
      stop
    endif
    read(io_work,*)e_scf(ifile)
    write(io_output,*)'Energie SCF=',e_scf(ifile)

    rewind(io_work)                      !energie SCF+ polarisation
    call find_str(io_work,' energie SCF + CPP                        stockee sur la file 33', &
  &               istatus,iline,l_str)
    if(istatus==0) then
      write(*,*) 'Missing E SCF-POLARISATION card in PSHFCV file'
      stop
    endif
    read(io_work,*)e_cpp(ifile)
    write(io_output,*)'Energie SCF+CPP=',e_cpp(ifile)

    rewind(io_work)                       !energies
    call find_str(io_work,' energies monoelectroniques alpha         stockee sur la file 33', &
  &               istatus,iline,l_str)
    if(istatus==0) then
      write(*,*) 'Missing ENERGIES MONOELECTRONIQUES card in PSHFCV file'
      stop
    endif
    if(ifile==1)allocate(e_mono(n_basis))
    do iv=1,divise(n_basis,5)     ! lecture par blocs de 5
      read(io_work,*)(e_mono(ie),ie=(iv-1)*5+1,min(n_basis,iv*5))
      write(io_output,'(10f12.7)')(e_mono(ie),ie=(iv-1)*5+1,min(n_basis,iv*5))
    enddo

    rewind(io_work)                       !vecteurs propres
    call find_str(io_work,'          eigenvectors alpha',istatus,iline,l_str)
    if(istatus==0) then
      write(*,*) 'Missing EIGENVECTORS card in PSHFCV file'
      stop
    endif
    do il=1,9                ! saute 9 lignes
      read(io_work,*)
    enddo
    
    if(ifile==1) then
      allocate(v_mono(n_basis,n_basis),atom(n_basis),centre(n_basis),orb(n_basis),l_orb(n_basis))
      l_orb=0
    endif
    
    do ib=1,divise(nvec,12)     ! lecture par blocs de 12
      do iv=1,n_basis
        read(io_work,'(i4,a3,i4,x,a4,12f12.5)')i,atom(iv),centre(iv),orb(iv),&
 &          (v_mono(iv,ie),ie=(ib-1)*12+1,min((ib*12),nvec))
!        write(6,'(i4,a3,i4,x,a4,12f12.5)')i,atom(iv),centre(iv),orb(iv),&
! &          (v_mono(iv,ie),ie=(ib-1)*12+1,min((ib*12),nvec))
        do i=1,16
          if(orb(iv)==sym(i))l_orb(iv)=l_om(i)+1
        enddo
      enddo
      do il=1,8              ! saute 8 lignes
        read(io_work,*)
      enddo
    enddo
    
    ! RV2015 : ajout des moments dipolaires
    
    rewind(io_work)                       !integrales de moments dipolaires
    call find_str(io_work,'          z-dipole moment integrals',istatus,iline,l_str)
    if(istatus/=0) then
      do il=1,6                ! saute 6 lignes
        read(io_work,*)
      enddo
      if(ifile==1)allocate(int_dip_z(n_basis,n_basis))
      do ib=1,divise(n_basis,7)     ! lecture par blocs de 12
        do iv=1,n_basis
          read(io_work,'(i5,1x,3a4,7e15.8)')i,atom(iv),centre(iv),orb(iv),&
 &        (int_dip_z(iv,ie),ie=(ib-1)*7+1,min((ib*7),n_basis))
        enddo
        do il=1,5              ! saute 5 lignes
          read(io_work,*)
        enddo
      enddo
    endif  
    
    rewind(io_work)                       !integrales de moments dipolaires
    call find_str(io_work,'          y-dipole moment integrals',istatus,iline,l_str)
    if(istatus/=0) then
      nprint=1
      do il=1,6                ! saute 6 lignes
        read(io_work,*)
      enddo
      if(ifile==1)allocate(int_dip_y(n_basis,n_basis))
      do ib=1,divise(n_basis,7)     ! lecture par blocs de 12
        do iv=1,n_basis
          read(io_work,'(i5,1x,3a4,7e15.8)')i,atom(iv),centre(iv),orb(iv),&
 &        (int_dip_y(iv,ie),ie=(ib-1)*7+1,min((ib*7),n_basis))
        enddo
        do il=1,5              ! saute 5 lignes
          read(io_work,*)
        enddo
      enddo
    endif  
    
    rewind(io_work)                       !integrales de moments dipolaires
    call find_str(io_work,'          x-dipole moment integrals',istatus,iline,l_str)
    if(istatus/=0) then
      do il=1,6                ! saute 6 lignes
        read(io_work,*)
      enddo
      if(ifile==1)allocate(int_dip_x(n_basis,n_basis))
      do ib=1,divise(n_basis,7)     ! lecture par blocs de 12
        do iv=1,n_basis
          read(io_work,'(i5,1x,3a4,7e15.8)')i,atom(iv),centre(iv),orb(iv),&
 &        (int_dip_x(iv,ie),ie=(ib-1)*7+1,min((ib*7),n_basis))
        enddo
        do il=1,5              ! saute 5 lignes
          read(io_work,*)
        enddo
      enddo
    endif  
    
    close(io_work)
                             ! analyse des vecteurs propres
    n_ug=1               !hetero
    !if(at1==at2)n_ug=2   !homo
    write(io_output,*)'ATTENTION!!! EVERY SYSTEM TREATED AS HETERO-NUCLEAR'
                             !comptage des symetries
    n_mono=0
    if(n_ug==1)then      !heteronucleaire
      do ie =1,nvec      ! lambda de l'OM, par analyse de la 1ere composante non nulle sur le centre 1
        do iv=1,n_oa(1)  ! comptage du nombre de vecteurs de chaque sym.
          if(v_mono(iv,ie)/=0)then
            n_mono(l_orb(iv),1)=n_mono(l_orb(iv),1)+1
            n_tri(ie,1)=n_mono(l_orb(iv),1)
            n_tri(ie,2)=l_orb(iv)
            e_tri(n_mono(l_orb(iv),1),l_orb(iv),1,ifile)=e_mono(ie)+e_1sr(ifile)+(e_cpp(ifile)-e_scf(ifile))
            exit
          endif
        enddo
      enddo
    else                 !homonucleaire
      do ie =1,nvec
        do iv=1,n_oa(1)
          if(v_mono(iv,ie)/=0)then
!            write(6,*)'composante non nulle pour ie,iv=',ie,iv
!            write(6,*)v_mono(iv,ie),v_mono(iv+n_oa(1),ie)
            if(abs(v_mono(iv,ie)/v_mono(iv+n_oa(1),ie)-1.)<1.e-5)then       !g
!              write(6,*)'g'
              n_mono(l_orb(iv),1)=n_mono(l_orb(iv),1)+1
              n_tri(ie,1)=n_mono(l_orb(iv),1)
              e_tri(n_mono(l_orb(iv),1),l_orb(iv),1,ifile)=e_mono(ie)+e_1sr(ifile)+(e_cpp(ifile)-e_scf(ifile))
!              e_tri(n_tri,l_orb(iv),1,if)=e_mono(ie)
              exit
            elseif(abs(v_mono(iv,ie)/v_mono(iv+n_oa(1),ie)+1.)<1.e-5)then   !u
!              write(6,*)'u'
              n_mono(l_orb(iv),2)=n_mono(l_orb(iv),2)+1
              n_tri(ie,1)=n_mono(l_orb(iv),2)
              e_tri(n_mono(l_orb(iv),2),l_orb(iv),2,ifile)=e_mono(ie)+e_1sr(ifile)+(e_cpp(ifile)-e_scf(ifile))
!              e_tri(n_tri,l_orb(iv),2,if)=e_mono(ie)
              exit
	    else	! Just a temporary bugfix ! JD 04/06
!              write(6,*)'u'
              n_mono(l_orb(iv),2)=n_mono(l_orb(iv),2)+1
              n_tri(ie,1)=n_mono(l_orb(iv),2)
              e_tri(n_mono(l_orb(iv),2),l_orb(iv),2,ifile)=e_mono(ie)+e_1sr(ifile)+(e_cpp(ifile)-e_scf(ifile))
!              e_tri(n_tri,l_orb(iv),2,if)=e_mono(ie)
              exit
            endif
          endif
        enddo
      enddo
    endif
    write(io_output,*)'# de vecteurs propres dans chaque symetrie'
    write(io_output,'(4i4)')((n_mono(j,i),j=1,4),i=1,2)
    l_max=max(l_max,maxval(l_orb))   
    write(io_output,*)'l_max=',l_max
    
                                                  ! RV2015 : analyse des moments dipolaires
    if(nprint==1) then
      do i=1,nvec
        do j=1,nvec
         dip_mono_z(n_tri(i,1),n_tri(i,2),n_tri(j,1),n_tri(j,2),ifile)=0d0
         dip_mono_y(n_tri(i,1),n_tri(i,2),n_tri(j,1),n_tri(j,2),ifile)=0d0
         dip_mono_x(n_tri(i,1),n_tri(i,2),n_tri(j,1),n_tri(j,2),ifile)=0d0
           do il=1,n_basis
             do iv=1,n_basis
              dip_mono_z(n_tri(i,1),n_tri(i,2),n_tri(j,1),n_tri(j,2),ifile) = dip_mono_z(n_tri(i,1),n_tri(i,2),n_tri(j,1),n_tri(j,2),ifile)&
              + int_dip_z(il,iv)*v_mono(iv,j)*v_mono(il,i)
              dip_mono_y(n_tri(i,1),n_tri(i,2),n_tri(j,1),n_tri(j,2),ifile) = dip_mono_y(n_tri(i,1),n_tri(i,2),n_tri(j,1),n_tri(j,2),ifile)&
              + int_dip_y(il,iv)*v_mono(iv,j)*v_mono(il,i)
              dip_mono_x(n_tri(i,1),n_tri(i,2),n_tri(j,1),n_tri(j,2),ifile) = dip_mono_x(n_tri(i,1),n_tri(i,2),n_tri(j,1),n_tri(j,2),ifile)&
              + int_dip_x(il,iv)*v_mono(iv,j)*v_mono(il,i)
             enddo
           enddo  
        enddo
      enddo 
      
      ! degenerate state handling : valid only for s-(p/d/f...) transitions, heternuclear 
      do j=2,l_max
        do iv = 1,n_mono(j,1),2
           if((dip_mono_z(1,1,iv,j,ifile)==0d0).OR.(dip_mono_z(1,1,iv+1,j,ifile)==0d0)) then
              dip_mono_z(1,1,iv,j,ifile) = dip_mono_z(1,1,iv,j,ifile) + dip_mono_z(1,1,iv+1,j,ifile)
              dip_mono_z(1,1,iv+1,j,ifile) = 0d0
           endif  
           if((dip_mono_y(1,1,iv,j,ifile)==0d0).OR.(dip_mono_y(1,1,iv+1,j,ifile)==0d0)) then
              dip_mono_y(1,1,iv,j,ifile) = dip_mono_y(1,1,iv,j,ifile) + dip_mono_y(1,1,iv+1,j,ifile)
              dip_mono_y(1,1,iv+1,j,ifile) = 0d0
           endif   
           if((dip_mono_x(1,1,iv,j,ifile)==0d0).OR.(dip_mono_x(1,1,iv+1,j,ifile)==0d0)) then
              dip_mono_x(1,1,iv,j,ifile) = dip_mono_x(1,1,iv,j,ifile) + dip_mono_x(1,1,iv+1,j,ifile)
              dip_mono_x(1,1,iv+1,j,ifile) = 0d0
           endif   
        enddo
      enddo
      
      nprint=2 ! dipole integral found at at least one distance
    endif  
    
  enddo                          !fin liste fichiers
  deallocate(v_mono,atom,centre,orb,l_orb,oa1,e_mono)
  if(allocated(int_dip_z))deallocate(int_dip_z)
  if(allocated(int_dip_y))deallocate(int_dip_y)
  if(allocated(int_dip_x))deallocate(int_dip_x)
  if(allocated(oa2))deallocate(oa2)

  do iug=1,n_ug                  ! ouverture dans workdir des fichiers d'energie par symetrie
    do i=1,l_max
      if(n_ug==1)then
        file=workdir(1:l_work)//'/'//'e_'//oa_typ(i)
      else
        file=workdir(1:l_work)//'/'//'e_'//oa_typ(i)//ug(iug)
      endif
!        open(iug*10+i,file=workdir(1:l_work)//'/'//'e_'//oa_typ(i))
      inquire(file=file,EXIST=exist)
      if(exist)then
        write(io_output,*)'File ',file,' already exists'
        write(io_output,*)'opened as OLD'
        write(io_output,*)'new values will be appended and then sorted by R values'
        open(iug*10+i,file=file,status='old',position='append')
      else
        write(io_output,*)'File ',file,' opened as NEW'
        open(iug*10+i,file=file,status='new')
        write(iug*10+i,*)'#  R(au)   1/R(au)     E_scf(au)   E_cpp(au)   E_mono(au)'
      endif
    enddo
  enddo
  
  !fmt='(f9.3,'//cex2(nvec+3)//'f12.7)'
  ! DB 12/2013 : Pour éviter les "escaliers" quand la variation d'énergie devient faible d'un point à l'autre
  fmt='(f9.3,3f12.7,'//cex2(nvec)//'e19.10)'
  do ifile=1,nf
    do iug=1,n_ug
      do i=1,l_max
        write(iug*10+i,fmt)r_in(ifile),e_1sr(ifile),e_scf(ifile),e_cpp(ifile)-e_scf(ifile), &
 &           (e_tri(ie,i,iug,ifile),ie=1,n_mono(i,iug))
      enddo
    enddo
  enddo

  ! RV 2015 : ecriture des moments dipolaires avec le fondamental de la symetrie s
  !dz
  if(nprint/=-1) then
    do i=1,1
      do j=1,l_max
       file=workdir(1:l_work)//'/'//'dz_mono_1'//oa_typ(i)//oa_typ(j)
       inquire(file=file,EXIST=exist)
       if(exist)then
         open(50+j*10+i,file=file,status='old',position='append')
       else
        open(50+j*10+i,file=file,status='new')
        write(50+j*10+i,*)'#  R(au)   dip_mono(au)'
       endif
     enddo
    enddo   
 
    fmt='(f9.3,'//cex2(nvec)//'e19.10)'
    do ifile=1,nf
      do i=1,1
        do j=1,l_max
          write(50+j*10+i,fmt)r_in(ifile),(dip_mono_z(1,i,ie,j,ifile),ie=1,n_mono(j,1))
        enddo
      enddo
    enddo
  
  !dv
    do i=1,1
      do j=1,l_max
       file=workdir(1:l_work)//'/'//'dy_mono_1'//oa_typ(i)//oa_typ(j)
       inquire(file=file,EXIST=exist)
       if(exist)then
         open(50+j*10+i,file=file,status='old',position='append')
       else
        open(50+j*10+i,file=file,status='new')
        write(50+j*10+i,*)'#  R(au)   dip_mono(au)'
       endif
      enddo
    enddo   
 
    fmt='(f9.3,'//cex2(nvec)//'e19.10)'
    do ifile=1,nf
      do i=1,l_max
        do j=1,l_max
          write(50+j*10+i,fmt)r_in(ifile),(dip_mono_y(1,i,ie,j,ifile),ie=1,n_mono(j,1))
        enddo
      enddo
    enddo
  
  !dx
    do i=1,l_max
      do j=1,l_max
        file=workdir(1:l_work)//'/'//'dx_mono_1'//oa_typ(i)//oa_typ(j)
        inquire(file=file,EXIST=exist)
        if(exist)then
          open(50+j*10+i,file=file,status='old',position='append')
        else
         open(50+j*10+i,file=file,status='new')
         write(50+j*10+i,*)'#  R(au)   dip_mono(au)'
        endif
      enddo
    enddo   
 
    fmt='(f9.3,'//cex2(nvec)//'e19.10)'
    do ifile=1,nf
      do i=1,l_max
        do j=1,l_max
          write(50+j*10+i,fmt)r_in(ifile),(dip_mono_x(1,i,ie,j,ifile),ie=1,n_mono(j,1))
        enddo
      enddo
    enddo
  endif
  
  deallocate(r_in,e_1sr,e_cpp,e_scf)
  deallocate(e_tri,n_tri)
  
  close(io_input)
  close(io_output)

  return
end subroutine tri_mono

function longchain(fich)
  ! real length of a character chain
  implicit none
  integer :: i,longchain
  character(len=*),intent(in) :: fich
  do i = 1,120
    if(fich(i:i)==' ')exit
  enddo
  longchain = i-1
  return
end function longchain
