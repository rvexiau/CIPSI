subroutine makedata_ciro(workdir,f_mdip,f_ciro,prefix,sym_1,sym_2,calfa,ncore,ntr,ietat,jetat,rang)
! Creation du fichier de donnees pour le calcul des moments de transiton
! d'une combinaison de symetries donnee, a partir du fichier type F_MDIP
  use io_unit
  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: lwork,l_str,istatus,ist_1,ist_2,i,il0,il1,il2,l1,l2,ncore
  integer :: l_mdip,metat,ntr,rang
  integer,dimension(ntr) :: ietat,jetat
  real(dp) :: calfa(ncore)
  character(len=180) :: workdir,f_ciro,f_mdip
  character(len=80) :: ligne
  character(len=*) :: sym_1,sym_2
  character(len=10) :: prefix 
  logical :: ymoyen
  

  ymoyen = .true.                !result from bdav
  lwork=len_trim(workdir)
  l_mdip=len_trim(f_mdip)
  l1=len_trim(sym_1)
  l2=len_trim(sym_2)
  f_ciro=f_mdip(1:l_mdip)//sym_1(1:l1)//sym_2(1:l2)
  open(io_input+rang*100,file=workdir(1:lwork)//'/'//f_mdip)
  open(io_namelist+rang*100,file=workdir(1:lwork)//'/'//f_ciro)

  call find_str(io_input+rang*100,'&rofil',istatus,il0,l_str)    !namelist &ROFIL
  if(istatus==0)then
    rewind io_input+rang*100
    call find_str(io_input+rang*100,'&ROFIL',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &ROFIL card in ',f_mdip
      stop
    endif
  endif
  !write( 6, * ) '&ROFIL debute en ligne',il0
  call find_str(io_input+rang*100,'&end',ist_1,il1,l_str)
  rewind io_input+rang*100
  do i=1,il0
    read(io_input+rang*100,*)
  enddo
  call find_str(io_input+rang*100,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &ROFIL not closed in ',f_mdip
    stop
  endif
  !write( 6, * ) '&ROFIL finit en ligne',il0+min(il1,il2)
  rewind io_input+rang*100
  do i=1,il0-1
    read(io_input+rang*100,*)
  enddo
  do i=il0,il0+min(il1,il2)-1
    read(io_input+rang*100,'(a80)')ligne
    write(io_namelist+rang*100,'(a80)')ligne
  enddo
  write(io_namelist+rang*100,'(a18)')"prefix='"//trim(prefix)//"'," 
  write(io_namelist+rang*100,'(a16)')'f04a='//"'f04"//sym_1//"',"      !insertion de f04a=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a16)')'f04b='//"'f04"//sym_2//"',"      !insertion de f04b=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a20)')'detcla='//"'dcl"//sym_1//"',"  !insertion de bdveca=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a20)')'detclb='//"'dcl"//sym_2//"',"  !insertion de bdvecb=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a25)')'detspina='//"'dspin"//sym_1//"',"  !insertion de bdveca=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a25)')'detspinb='//"'dspin"//sym_2//"',"  !insertion de bdvecb=, en fin de NAMELIST  
  write(io_namelist+rang*100,'(a20)')'bdveca='//"'bd"//sym_1//"',"  !insertion de bdveca=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a20)')'bdvecb='//"'bd"//sym_2//"',"  !insertion de bdvecb=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a4)')'&END'                   !
  write(io_namelist+rang*100,'(a80)')
  rewind io_input+rang*100

  call find_str(io_input+rang*100,'&roinp',istatus,il0,l_str)    !namelist &ROINP
  if(istatus==0)then
    rewind io_input+rang*100
    call find_str(io_input+rang*100,'&ROINP',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &ROINP card in ',f_mdip
      stop
    endif
  endif
  !write( 6, * ) '&ROINP debute en ligne',il0
  call find_str(io_input+rang*100,'&end',ist_1,il1,l_str)
  rewind io_input+rang*100
  do i=1,il0
    read(io_input+rang*100,*)
  enddo
  call find_str(io_input+rang*100,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &ROINP not closed in ',f_mdip
    stop
  endif
  !write( 6, * ) '&ROINP finit en ligne',il0+min(il1,il2)
  rewind io_input+rang*100
  do i=1,il0-1
    read(io_input+rang*100,*)
  enddo
  do i=il0,il0+min(il1,il2)-1
    read(io_input+rang*100,'(a80)')ligne
    write(io_namelist+rang*100,'(a80)')ligne
  enddo
  write(io_namelist+rang*100,'(a7,30(i2,a1))')' IETAT=',(ietat(i),',',i=1,ntr)
  write(io_namelist+rang*100,'(a7,30(i2,a1))')' JETAT=',(jetat(i),',',i=1,ntr)
  write(io_namelist+rang*100,'(a7,l1,a1)')'YMOYEN=',ymoyen,','  !insertion de YMOYEN=, en fin de NAMELIST
  write(io_namelist+rang*100,'(a7,10f10.4)')' CALFA=',(calfa(i),i=1,ncore)    !polarisabilites tirees de f_pshf  
  write(io_namelist+rang*100,'(a4)')'&END'                   !
  write(io_namelist+rang*100,'(a80)')
  close(io_input+rang*100)
  close(io_namelist+rang*100)

  return
end