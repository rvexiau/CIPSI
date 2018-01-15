subroutine makedata_cip(workdir,f_sym,f_bd,prefix,group,isymat,sym,nusym,itsym,mult_tab,isz,nelac,metat,noac,noac_ref,nexmax,l_ref,rang)
!  creation du fichier de donnees de cip,moy,bd, pour chaque symetrie
!  a partir d'un fichier f_sym type, precise dans autocip.in
!  les informations relatives au spin (ISYMAT), nom de symetrie, numero de
!  symetrie (NUSYM) et nombre d'etats desires (NETAT)
!  sont traitees automatiquement ici.
!  Les autres donnees ne sont pas changees, et peuvent donc etre
!  ajustees a chaque calcul par l'utilisateur
!
!  CONVENTION sur la notation des symetries: voir programme principal
 use grid_data
  implicit none
  integer :: lwork,rang,i,j,k,isymat(6)
  integer :: l_sym,metat,nusym,noac,noac2,noac_ref,nelac,isz,numref,sym_det,nexmax,lecdet
  integer, dimension(noac) :: itsym,inum
  integer, dimension(:), allocatable :: numac
  character(len=*) :: workdir,f_bd,f_sym
  character(len=*) :: sym,group
  character(len=180) :: f_bd_loop  
  character(len=10) :: prefix
  logical :: l_ref,qictot
  integer,dimension(14,14) :: mult_tab
  
  lwork=len_trim(workdir)
   qictot=.true.
   lecdet=4
   l_sym=len_trim(f_sym)  
   f_bd=f_sym(1:l_sym)//sym     

   if (noac_ref.gt.noac) then
      noac_ref=-1
      write(*,*) 'warning, noac_ref too big : FCI method'
   endif

  if (noac_ref.ge.0) then
    select case(nelac)
    case (1)
       numref=1
       sym_det=1
    case (2)
       numref=1
       sym_det=itsym(1)
    case (3)
       numref=2
       sym_det=itsym(1)  
       sym_det=mult_tab(sym_det,itsym(1))       
    case (4)
       numref=2
       sym_det=itsym(1)
       sym_det=mult_tab(sym_det,itsym(1))
       sym_det=mult_tab(sym_det,itsym(2))       
    case DEFAULT
       write(*,*) '5+ electrons not implemented in makedata_cip'
    end select
    
    noac2=noac_ref+numref
    allocate(numac(noac2))    
    ! add occupied orbital
    do i=1,numref
      numac(i)=i
    enddo
    ! add virtual orbitals of correct symmetry
    k=1
    do i= numref+1,noac
      if(k.gt.noac_ref) exit
      if(mult_tab(sym_det,itsym(i))==nusym) then
        numac(k+numref)=i
        k=k+1
      endif
    enddo
    
    call makedata_cip_template(workdir,lwork,f_sym,f_bd,prefix,group,(/-1,0,0,0,0,0/),sym,nusym,numac,isz,nelac,metat,noac2,nexmax,qictot,l_ref,rang,lecdet)
  else
    allocate(numac(noac))
    do i=1,noac
      numac(i)=i
    enddo
    call makedata_cip_template(workdir,lwork,f_sym,f_bd,prefix,group,(/-1,0,0,0,0,0/),sym,nusym,numac,isz,nelac,metat,noac,nexmax,qictot,l_ref,rang,lecdet)
  endif
  
  lecdet=5
  l_sym=len_trim(f_sym)
  f_bd_loop=trim(f_bd)//'_loop'
  call makedata_cip_template(workdir,lwork,f_sym,f_bd_loop,prefix,group,(/-1,0,0,0,0,0/),sym,nusym,numac,isz,nelac,metat,-1,-1,qictot,l_ref,rang,lecdet)
  f_bd_loop=trim(f_bd)//'_end'
  qictot=.true.
  call makedata_cip_template(workdir,lwork,f_sym,f_bd_loop,prefix,group,isymat,sym,nusym,numac,isz,nelac,metat,-1,-1,qictot,l_ref,rang,lecdet)

  deallocate(numac)

  return
end subroutine


subroutine makedata_cip_template(workdir,lwork,f_sym,f_bd,prefix,group,isymat,sym,nusym,numac,isz,nelac,metat,noac,nexmax,qictot,l_ref,rang,lecdet)
!  creation du fichier de donnees de cip,moy,bd, pour chaque symetrie
!  a partir d'un fichier f_sym type, precise dans autocip.in
!  les informations relatives au spin (ISYMAT), nom de symetrie, numero de
!  symetrie (NUSYM), nombre d'orbitales actives (NOAC)
!  et nombre d'etats desires (NETAT)
!  sont traitees automatiquement ici.
!  Les autres donnees ne sont pas changees, et peuvent donc etre
!  ajustees a chaque calcul par l'utilisateur
!
!  CONVENTION sur la notation des symetries: voir programme principal
  use io_unit
  implicit none
  integer :: lwork,l_str,istatus,ist_1,ist_2,i,il0,il1,il2
  integer :: l_sym,metat,nusym,noac,nelac,isz,rang,isymat(6),nexmax,iselec,lecdet
  integer, dimension(abs(noac)) :: numac  
  character(len=*) :: workdir,f_bd,f_sym
  character(len=80) :: ligne
  character(len=*) :: sym,group
  character(len=7) :: n_list,n_list_c
  character(len=10) :: prefix
  character(len=3) :: cex3  
  logical :: ymoyen,l_ref,qhef,qictot,parity

  iselec=1
  ymoyen=.true.
  qhef=.true.
  open(io_input+rang*100,file=workdir(1:lwork)//'/'//f_sym)
  open(io_namelist+rang*100,file=workdir(1:lwork)//'/'//f_bd)

  n_list = '&cipfil'
  n_list_c= '&CIPFIL'

  call find_str(io_input+rang*100,n_list,istatus,il0,l_str)  !namelist &CIPFIL ou &FCIFIL
  if(istatus==0)then
    rewind io_input+rang*100
    call find_str(io_input+rang*100,n_list_c,istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing ',n_list,' card in ',f_sym
      stop
    endif
  endif
  !write( 6, * ) n_list,' debute en ligne',il0
  call find_str(io_input+rang*100,'&end',ist_1,il1,l_str)
  rewind io_input+rang*100
  do i=1,il0
    read(io_input+rang*100,*)
  enddo
  call find_str(io_input+rang*100,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist ',n_list,' not closed in ',f_sym
    stop
  endif
  !write( 6, * ) n_list,' finit en ligne',il0+min(il1,il2)
  rewind io_input+rang*100
  do i=1,il0-1
    read(io_input+rang*100,*)
  enddo
  do i=il0,il0+min(il1,il2)-1
    read(io_input+rang*100,'(a80)')ligne
    write(io_namelist+rang*100,'(a80)')ligne
  enddo
  write(io_namelist+rang*100,'(a18)')"prefix='"//trim(prefix)//"'," 
  write(io_namelist+rang*100,'(a14)')'f04='//"'f04"//sym//"',"  !insertion de f04=, en fin de NAMELIST
!  if(l_ref)then
!    write(io_namelist+rang*100,'(a17)')' det='//"'drcl"//sym//"',"  !insertion de detcl=, en fin de NAMELIST
!  else
!    write(io_namelist+rang*100,'(a16)')' det='//"'dcl"//sym//"',"  !insertion de detcl=, en fin de NAMELIST
!  endif
  if(l_ref)then
    write(io_namelist+rang*100,'(a19)')' bdvec='//"'bdr"//sym//"',"  !insertion de bdvec=, en fin de NAMELIST
  else
    write(io_namelist+rang*100,'(a18)')' bdvec='//"'bd"//sym//"',"  !insertion de bdvec=, en fin de NAMELIST
  endif
  write(io_namelist+rang*100,'(a4)')'&END'                   !
  rewind io_input+rang*100
  write(io_namelist+rang*100,'(a80)')

  call find_str(io_input+rang*100,'&icinp',istatus,il0,l_str)  !namelist &ICINP
  if(istatus==0)then
    rewind io_input+rang*100
    call find_str(io_input+rang*100,'&ICINP',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &ICINP card in ',f_sym
      stop
    endif
  endif
  !write( 6, * ) '&ICINP debute en ligne',il0
  call find_str(io_input+rang*100,'&end',ist_1,il1,l_str)
  rewind io_input+rang*100
  do i=1,il0
    read(io_input+rang*100,*)
  enddo
  call find_str(io_input+rang*100,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &ICINP not closed in ',f_sym
    stop
  endif
  !write( 6, * ) '&ICINP finit en ligne',il0+min(il1,il2)
  rewind (io_input+rang*100)
  do i=1,il0-1
    read(io_input+rang*100,*)
  enddo
  read(io_input+rang*100,'(a80)')ligne
  write(io_namelist+rang*100,'(a80)')ligne
  if(noac.ge.0) write(io_namelist+rang*100,'(a7,'//cex3(noac)//'(i0,a1))')'NUMAC=',(numac(i),",",i=1,noac)    !insertion de NUMAC=, en debut de NAMELIST
  do i=il0+1,il0+min(il1,il2)-1
    read(io_input+rang*100,'(a80)')ligne
    write(io_namelist+rang*100,'(a80)')ligne
  enddo
  write(io_namelist+rang*100,'(a5,i0,a1)')' ISZ=',isz,','      !insertion de ISZ=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a7,i0,a1)')' NUSYM=',nusym,','  !insertion de NUSYM=, en fin de NAMELIST
  if(noac.ge.0) write(io_namelist+rang*100,'(a6,i0,a1)')' NOAC=',noac,','    !insertion de NOAC=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a7,i0,a1)')' METAT=',metat,','  !insertion de METAT=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a8,i0,a1)')' LECDET=',lecdet,','  !insertion de LECDET=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a11)')' ITYPER=-1,' !insertion de ITYPER=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a8,l1,a1)')' YMOYEN=',ymoyen,','  !insertion de YMOYEN=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a7,i0,a1)')' NELAC=',nelac,','  !insertion de NELAC=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a8)')' YBRD=F,'  !insertion de YBRD=, en fin de NAMELIST  
  if(nexmax.gt.0) write(io_namelist+rang*100,'(a7,i0,a1)')' NEMAX=',nexmax,','  !insertion de NEMAX=, en fin de NAMELIST  
  parity = MOD(nelac,2)
  if(parity) then
    write(io_namelist+rang*100,'(a8)')' YION=T,'  !insertion de YION=, en fin de NAMELIST 
  else
    write(io_namelist+rang*100,'(a8)')' YION=F,'  !insertion de YION=, en fin de NAMELIST 
  endif
  write(io_namelist+rang*100,'(a8,i3,a1)')' ISELEC=',iselec,','  !insertion de ISELEC=, en fin de NAMELIST        
  write(io_namelist+rang*100,'(a10)')' TEST=0d0,'
  write(io_namelist+rang*100,'(a9)')' TAU=0d0,'
  write(io_namelist+rang*100,'(a4)')'&END'                   !
  write(io_namelist+rang*100,'(a80)')
  rewind (io_input+rang*100)

  call find_str(io_input+rang*100,'&moyfil',istatus,il0,l_str)  !namelist &MOYFIL
  if(istatus==0)then
    rewind (io_input+rang*100)
    call find_str(io_input+rang*100,'&MOYFIL',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &MOYFIL card in ',f_sym
      stop
    endif
  endif
  call find_str(io_input+rang*100,'&end',ist_1,il1,l_str)
  rewind (io_input+rang*100)
  do i=1,il0
    read(io_input+rang*100,*)
  enddo
  call find_str(io_input+rang*100,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &MOYFIL not closed in ',f_sym
    stop
  endif
  rewind (io_input+rang*100)
  do i=1,il0-1
    read(io_input+rang*100,*)
  enddo
  do i=il0,il0+min(il1,il2)-1
    read(io_input+rang*100,'(a80)')ligne
    write(io_namelist+rang*100,'(a80)')ligne
  enddo
  write(io_namelist+rang*100,'(a18)')"prefix='"//trim(prefix)//"'," 
  write(io_namelist+rang*100,'(a15)')'f04='//"'f04"//sym//"',"  !insertion de f04=, en fin de NAMELIST
  if(l_ref)then
    write(io_namelist+rang*100,'(a21)')' det_cl='//"'drcl"//sym//"',"  !insertion de detcl=, en fin de NAMELIST
    write(io_namelist+rang*100,'(a22)')' det_spin='//"'drspin"//sym//"',"  !insertion de detspin=, en fin de NAMELIST         
  else
    write(io_namelist+rang*100,'(a20)')' det_cl='//"'dcl"//sym//"',"  !insertion de detcl=, en fin de NAMELIST
    write(io_namelist+rang*100,'(a22)')' det_spin='//"'dspin"//sym//"',"  !insertion de detspin=, en fin de NAMELIST      
  endif
  write(io_namelist+rang*100,'(a4)')'&END'                   !
  write(io_namelist+rang*100,'(a80)')
  rewind (io_input+rang*100)

  call find_str(io_input+rang*100,'&moyen',istatus,il0,l_str)  !namelist &MOYEN
  if(istatus==0)then
    rewind (io_input+rang*100)
    call find_str(io_input+rang*100,'&MOYEN',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &MOYEN card in ',f_sym
      stop
    endif
  endif
  call find_str(io_input+rang*100,'&end',ist_1,il1,l_str)
  rewind (io_input+rang*100)
  do i=1,il0
    read(io_input+rang*100,*)
  enddo
  call find_str(io_input+rang*100,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &MOYEN not closed in ',f_sym
    stop
  endif
  rewind (io_input+rang*100)
  do i=1,il0-1
    read(io_input+rang*100,*)
  enddo
  do i=il0,il0+min(il1,il2)-1                     
    read(io_input+rang*100,'(a80)')ligne
    write(io_namelist+rang*100,'(a80)')ligne
  enddo
  write(io_namelist+rang*100,'(a9)')' TAU=0d0,'!insertion de TAU=, en fin de NAMELIST      
  write(io_namelist+rang*100,'(a11)')' yspipro=F,'    
  if(isymat(1).ge.0) then
    write(io_namelist+rang*100,'(a7,i0,a1)')' NELAC=',nelac,','  
    write(io_namelist+rang*100,'(a8,6(i2,a1))')' ISYMAT=',(isymat(i),',',i=1,6)  
    if((isymat(1).eq.1).or.(isymat(2).eq.1).or.(isymat(3).eq.1) .or.(isymat(4).eq.1) .or.(isymat(5).eq.1) .or.(isymat(6).eq.1))  then
      write(io_namelist+rang*100,'(a11)')' yspipro=T,'  
    endif  
  endif  
  write(io_namelist+rang*100,'(a4)')'&END'
  write(io_namelist+rang*100,'(a80)')
  rewind (io_input+rang*100)

  call find_str(io_input+rang*100,'&bdfil',istatus,il0,l_str)  !namelist &BDFIL
  if(istatus==0)then
    rewind (io_input+rang*100)
    call find_str(io_input+rang*100,'&BDFIL',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &BDFIL card in ',f_sym
      stop
    endif
  endif
  !write( 6, * ) '&BDFIL debute en ligne',il0
  call find_str(io_input+rang*100,'&end',ist_1,il1,l_str)
  rewind (io_input+rang*100)
  do i=1,il0
    read(io_input+rang*100,*)
  enddo
  call find_str(io_input+rang*100,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &BDFIL not closed in ',f_sym
    stop
  endif
  !write( 6, * ) '&BDFIL finit en ligne',il0+min(il1,il2)
  rewind (io_input+rang*100)
  do i=1,il0-1
    read(io_input+rang*100,*)
  enddo
  do i=il0,il0+min(il1,il2)-1
    read(io_input+rang*100,'(a80)')ligne
    write(io_namelist+rang*100,'(a80)')ligne
  enddo
  write(io_namelist+rang*100,'(a18)')"prefix='"//trim(prefix)//"'," 
  write(io_namelist+rang*100,'(a15)')'f04='//"'f04"//sym//"',"      !insertion de f04=, en fin de NAMELIST
  if(l_ref)then
    write(io_namelist+rang*100,'(a20)')' bdvec='//"'bdr"//sym//"',"  !insertion de bdvec=, en fin de NAMELIST
    write(io_namelist+rang*100,'(a21)')' det_cl='//"'drcl"//sym//"',"  !insertion de detcl=, en fin de NAMELIST  
    write(io_namelist+rang*100,'(a22)')' det_spin='//"'drspin"//sym//"',"  !insertion de detspin=, en fin de NAMELIST       
  else
    write(io_namelist+rang*100,'(a19)')' bdvec='//"'bd"//sym//"',"  !insertion de bdvec=, en fin de NAMELIST
    write(io_namelist+rang*100,'(a20)')' det_cl='//"'dcl"//sym//"',"  !insertion de detcl=, en fin de NAMELIST  
    write(io_namelist+rang*100,'(a22)')' det_spin='//"'dspin"//sym//"',"  !insertion de detspin=, en fin de NAMELIST      
  endif
  write(io_namelist+rang*100,'(a4)')'&END'                   !
  write(io_namelist+rang*100,'(a80)')
  rewind (io_input+rang*100)

  call find_str(io_input+rang*100,'&option',istatus,il0,l_str)  !namelist &OPTION
  if(istatus==0)then
    rewind (io_input+rang*100)
    call find_str(io_input+rang*100,'&OPTION',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &OPTION card in ',f_sym
      stop
    endif
  endif
  !write( 6, * ) '&OPTION debute en ligne',il0
  call find_str(io_input+rang*100,'&end',ist_1,il1,l_str)
  rewind (io_input+rang*100)
  do i=1,il0
    read(io_input+rang*100,*)
  enddo
  call find_str(io_input+rang*100,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &OPTION not closed in ',f_sym
    stop
  endif
  !write( 6, * ) '&OPTION finit en ligne',il0+min(il1,il2)
  rewind (io_input+rang*100)
  do i=1,il0-1
    read(io_input+rang*100,*)
  enddo
  do i=il0,il0+min(il1,il2)-1
    read(io_input+rang*100,'(a80)')ligne
    write(io_namelist+rang*100,'(a80)')ligne
  enddo
  write(io_namelist+rang*100,'(a7,i3,a1)')'METAT=',metat,','  !insertion de METAT=, en fin de NAMELIST 
  write(io_namelist+rang*100,'(a6,l1,a1)')' QHEF=',qhef,','  !insertion de QHEF=, en fin de NAMELIST 
  write(io_namelist+rang*100,'(a8,l1,a1)')' QICTOT=',qictot,','  !insertion de QICTOT=, en fin de NAMELIST  
  write(io_namelist+rang*100,'(a4)')'&END'                   !
  write(io_namelist+rang*100,'(a80)')
  rewind (io_input+rang*100)
  
  call find_str(io_input+rang*100,'&sort',istatus,il0,l_str)  !namelist &SORT
  if(istatus==0)then
    rewind (io_input+rang*100)
    call find_str(io_input+rang*100,'&SORT',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &SORT card in ',f_sym
      stop
    endif
  endif
  !write( 6, * ) '&OPTION debute en ligne',il0
  call find_str(io_input+rang*100,'&end',ist_1,il1,l_str)
  rewind (io_input+rang*100)
  do i=1,il0
    read(io_input+rang*100,*)
  enddo
  call find_str(io_input+rang*100,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &SORT not closed in ',f_sym
    stop
  endif
  !write( 6, * ) '&OPTION finit en ligne',il0+min(il1,il2)
  rewind (io_input+rang*100)
  do i=1,il0-1
    read(io_input+rang*100,*)
  enddo
  do i=il0,il0+min(il1,il2)-1
    read(io_input+rang*100,'(a80)')ligne
    write(io_namelist+rang*100,'(a80)')ligne
  enddo
  write(io_namelist+rang*100,*)'GROUP='//group//','  !insertion de GROUP=, en fin de NAMELIST 
  write(io_namelist+rang*100,'(a7,i0,a1)')' NUSYM=',nusym,','  !insertion de NUSYM=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a7,i0,a1)')' NELAC=',nelac,','  !insertion de NELAC=, en fin de NAMELIST   
  write(io_namelist+rang*100,'(a4)')'&END'                   !
  write(io_namelist+rang*100,'(a80)')
  rewind (io_input+rang*100)
  close(io_input+rang*100)
  close(io_namelist+rang*100)
  return
end
