subroutine makedata_som(workdir,prefix,f_diab,f_som,noac)

  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: lwork,l_diab,nusym,noac,istatus,il0,l_str,ist_1,ist_2
  integer :: il1,il2,i
  character(len=180) :: workdir,f_diab,f_som,ligne
  character(len=10) :: prefix

  lwork=len_trim(workdir)
  l_diab=len_trim(f_diab)
  f_som=f_diab(1:l_diab)//'som'

  open(60,file=workdir(1:lwork)//'/'//f_diab)
  open(70,file=workdir(1:lwork)//'/'//f_som)
  !write( 6, * ) 'lecture de  ',f_diab

  call find_str(60,'&somfil',istatus,il0,l_str)
  if(istatus==0)then
    rewind 60
    call find_str(60,'&SOMFIL',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &SOMFIL card in ',f_diab
      stop
    endif
  endif
  !write( 6, * ) '&SOMFIL debute en ligne',il0
  call find_str(60,'&end',ist_1,il1,l_str)
  rewind 60
  do i=1,il0
    read(60,*)
  enddo
  call find_str(60,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &SOMFIL not closed in ',f_diab
    stop
  endif
  !write( 6, * ) '&SOMFIL finit en ligne',il0+min(il1,il2)
  rewind 60
  do i=1,il0-1
    read(60,*)
  enddo 
  do i=il0+1,il0+min(il1,il2)-1
    read(60,'(a80)')ligne
    write(70,'(a80)')ligne
  enddo
  write(70,'(a18)')"prefix='"//trim(prefix)//"',"
  write(70,'(a4)')'&END'                   !
  write(70,'(a80)')

  rewind 60
  call find_str(60,'&datom',istatus,il0,l_str)
  if(istatus==0)then
    rewind 60
    call find_str(60,'&DATOM',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &DATOM card in ',f_diab
      stop
    endif
  endif
  !write( 6, * ) '&DATOM debute en ligne',il0
  call find_str(60,'&end',ist_1,il1,l_str)
  rewind 60
  do i=1,il0
    read(60,*)
  enddo
  call find_str(60,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &DATOM not closed in ',f_diab
    stop
  endif
  !write( 6, * ) '&DATOM finit en ligne',il0+min(il1,il2)
  rewind 60
  do i=1,il0-1
    read(60,*)
  enddo
  do i=il0,il0+min(il1,il2)-1
    read(60,'(a80)')ligne
    write(70,'(a80)')ligne
  enddo
  write(70,'(a7,i3,a1)')'  naor=',noac,','
  write(70,'(a7,i3,a1)')'  nomr=',noac,','
  write(70,'(a4)')'&END'                   !
  write(70,'(a80)')

  rewind 60
  close(60)
  close(70)
  return
end

subroutine makedata_hdiab(workdir,f_diab,f_hdiab,sym,noca,nrot,num,numr)

  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: lwork,l_diab,nusym,noac,istatus,il0,l_str,ls
  integer :: ist_1,ist_2,il1,il2,i
  integer :: noca,nrot
  integer, dimension(:) :: num(nrot),numr(nrot)
  character(len=180) :: workdir,f_diab,ligne
  character(len=*) :: sym
  character(len=*),intent(out) :: f_hdiab

  lwork=len_trim(workdir)  
  l_diab=len_trim(f_diab)
  ls=len_trim(sym)
  f_hdiab=f_diab(1:l_diab)//sym(1:ls)

  open(60,file=workdir(1:lwork)//'/'//f_diab)
  open(70,file=workdir(1:lwork)//'/'//f_hdiab)
  !write( 6, * ) 'lecture de  ',f_diab

  call find_str(60,'&hefil',istatus,il0,l_str)
  if(istatus==0)then
    rewind 60
    call find_str(60,'&HEFIL',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &HEFIL card in ',f_hdiab
      stop
    endif
  endif
  !write( 6, * ) '&HEFIL debute en ligne',il0
  call find_str(60,'&end',ist_1,il1,l_str)
  rewind 60
  do i=1,il0
    read(60,*)
  enddo
  call find_str(60,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &HEFIL not closed in ',f_hdiab
    stop
  endif
  !write( 6, * ) '&HEFIL finit en ligne',il0+min(il1,il2)
  rewind 60
  do i=1,il0-1
    read(60,*)
  enddo
  do i=il0,il0+min(il1,il2)-1
    read(60,'(a80)')ligne
    write(70,'(a80)')ligne
  enddo
  write(70,'(a13)')'det='//"'dcl"//sym//"'," 
  write(70,'(a15)')'detr='//"'drcl"//sym//"',"  
  write(70,'(a14)')'bdvec='//"'bd"//sym//"',"  !insertion de bdvec=, en fin de NAMELIST
  write(70,'(a16)')'bdvecr='//"'bdr"//sym//"',"  !insertion de bdvec=, en fin de NAMELIST
  write(70,'(a17)')'diavec='//"'dia"//sym//"',"  !insertion de bdvec=, en fin de NAMELIST
  write(70,'(a4)')'&END'                   !
  write(70,'(a80)')

  rewind 60
  call find_str(60,'&hefinp',istatus,il0,l_str)
  if(istatus==0)then
    rewind 60
    call find_str(60,'&HEFINP',istatus,il0,l_str)
    if(istatus==0) then
      write( 6, * )  'missing &HEFINP card in ',f_hdiab
      stop
    endif
  endif
  !write( 6, * ) '&HEFINP debute en ligne',il0
  call find_str(60,'&end',ist_1,il1,l_str)
  rewind 60
  do i=1,il0
    read(60,*)
  enddo
  call find_str(60,'&END',ist_2,il2,l_str)
  if((ist_1==0).and.(ist_2==0)) then
    write( 6, * )  'namelist &HEFINP not closed in ',f_hdiab
    stop
  endif
  !write( 6, * ) '&HEFINP finit en ligne',il0+min(il1,il2)
  rewind 60
  do i=1,il0-1
    read(60,*)
  enddo
  do i=il0,il0+min(il1,il2)-1
    read(60,'(a80)')ligne
    write(70,'(a80)')ligne
  enddo
  write(70,'(a6,i2)')' noca=',noca
  write(70,'(a6,i2)')' nrot=',nrot
  write(70,'(a5,30(i2,a1))')' num=',(num(i),',',i=1,nrot)
  write(70,'(a6,30(i2,a1))')' numr=',(numr(i),',',i=1,nrot)
  write(70,'(a4)')'&END'                   !
  rewind 60
  write(70,'(a80)')

  rewind 60
  close(60)
  close(70)
  return
end