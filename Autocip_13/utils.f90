subroutine fetch_refom(prefix,lpr,workdir,refdir,f_diab)
! procedure to copy the saved reference states into the working directory
  implicit none

  integer :: lpr,lw,lwr
  character(len=*) :: workdir,refdir,f_diab
  character(len=10) :: prefix
  character(len=400) :: command
  character(len=10) :: info,som,ijom,omref	! JD 05/06 namelist/somfil/ for SOM
  
  namelist/somfil/prefix,info,ijom,omref

  omref='omref       '
  lw=len_trim(workdir)
  lwr=len_trim(refdir)  
  open(10,file=workdir(1:lw)//'/'//f_diab)     !lecture du NAMELIST somfil
  read(10,somfil)
  close(10)

  command='cp -f '//refdir(1:lwr)//'/'//prefix(1:lpr)//'om '//workdir(1:lw)//'/'//prefix(1:lpr)//omref
  write(6,*) command
  call system(command)

end subroutine fetch_refom

subroutine fetch_refbd(prefix,lpr,workdir,refdir,sym)

! procedure to safe the output of ijkl as reference for HDIAB
  implicit none

  integer :: lpr,lw,lwr
  character(len=*) :: workdir,refdir,sym
  character(len=10) :: prefix
  character(len=400) :: command

  lw=len_trim(workdir)
  lwr=len_trim(refdir)
  command='cp -f '//refdir(1:lwr)//'/'//prefix(1:lpr)//"bdr"//sym//' '//workdir(1:lw)//'/'//prefix(1:lpr)//"bdr"//sym
  call system(command)
  command='cp -f '//refdir(1:lwr)//'/'//prefix(1:lpr)//"drcl"//sym//' '//workdir(1:lw)//'/'//prefix(1:lpr)//"drcl"//sym
  call system(command)

end subroutine fetch_refbd

subroutine copy_file(file,source,dest)

!  use iflport
  implicit none

  character(len=*) :: file,source,dest
  integer :: lf,ls,ld,i_sys
  character(len=400) :: command

  lf=len_trim(file)
  ls=len_trim(source)
  ld=len_trim(dest)

  command='cp -f '//source(1:ls)//file(1:lf)//' '//dest(1:ld)//file(1:lf)
  call system(command)
!  i_sys = systemqq(command)
!  write( 6, * ) i_sys,command

end subroutine copy_file

subroutine copy2_file(sourcefile,destfile)

!  use iflport
  implicit none

  character(len=*) :: sourcefile,destfile
  integer :: ls,ld,i_sys
  character(len=300) :: command

  ls=len_trim(sourcefile)
  ld=len_trim(destfile)

  command='cp -f '//sourcefile(1:ls)//' '//destfile(1:ld)
  call system(command)  
  ! i_sys = systemqq(command)
  ! write( 6, * ) i_sys,command

end subroutine copy2_file

function cex(nj)
  integer, intent(in) :: nj
  character(len=*) :: cex
  character(len=2) :: cex2
  character(len=3) :: cex3
  
  if(nj.gt.99)then
     cex = cex3(nj)
  else
     cex = cex2(nj)
  endif

  return
end function cex

function cex2(nj)
  ! Converts an integer into a string of length 2
  implicit none
  integer, intent(in) :: nj
  character(len=2) :: cex2
  integer :: dnj,unj

  if(nj.gt.99) then
     write(*,*) 'Argument too big in function cex2: ',nj
     stop
  endif
  dnj = nj/10
  unj = nj - dnj*10
  cex2 = achar(dnj+48)//achar(unj+48)
  return
end function cex2

function cex3(nj)
  ! Converts an integer into a string of length 2
  implicit none
  integer, intent(in) :: nj
  character(len=3) :: cex3
  integer :: cnj,dnj,unj

  if(nj.gt.999) then
     write(*,*) 'Argument too big in function cex3: ',nj
     stop
  endif
  cnj = nj/100
  dnj = (nj - 100*cnj)/10
  unj = nj - cnj*100 - dnj*10
  cex3 = achar(cnj+48)//achar(dnj+48)//achar(unj+48)
  return
end function cex3

function divise(n1,n2)
  ! donne le quotient de n1 par n2, augmente de 1 si le reste est different de 0
  ! utile pour calculer en combien de blocs n2 divise n1
  implicit none
  integer, intent(in) :: n1,n2
  integer :: divise
  if(mod(n1,n2)==0)then
    divise=n1/n2
  else
    divise=n1/n2+1
  endif
  return
end function divise

function longu(fich)
  ! donne la position du premier caractere non blanc d'une chaine
  implicit none
  integer :: i,longu
  character(len=*),intent(in) :: fich
  do i = 1,120
    if(fich(i:i)/=' ')exit
  enddo
  longu = i

  return
end function longu

function eip(atom)
  !
  ! Retruns the ionization potential of 'atom' in atomic units
  !
  ! atom = 'nAa' where n = mass number
  !                    Aa = symbol
  ! source: "appraoched" numbers from quantum chmistry calculations
  !
  ! (after C. Dion, 28/07/2000: O.D. avril 2005)
  !
  implicit none
  integer, parameter :: dp=kind(1.d0)
  character(len=*), intent(in) :: atom
  real(dp) :: eip

  select case (trim(atom))
  case ('6Li')
     eip=0.198147
  case ('7Li')
     eip=0.198147
  case ('Li')
     eip=0.198147
  case ('23Na')
     eip=0.188859
  case ('Na')
     eip=0.188859
  case ('39K')
     eip=0.159521
  case ('40K')
     eip=0.159521
  case ('41K')
     eip=0.159521
  case ('K')
     eip=0.159521
  case ('85Rb')
     eip=0.15303
  case ('87Rb')
     eip=0.15303
  case ('Rb')
     eip=0.15303
  case ('133Cs')
     eip=0.143098
  case ('Cs')
     eip=0.143098
  case ('Fr')
     eip=0.14967049
  case default
     write(*,*) 'Element ',trim(atom),' unknown'
     eip = 0.
     !stop
  end select

end function eip

function znuc(atom)
  !
  ! Retruns the ionization potential of 'atom' in atomic units
  !
  ! atom = 'nAa' where n = mass number
  !                    Aa = symbol
  ! source: "appraoched" numbers from quantum chmistry calculations
  !
  ! (after C. Dion, 28/07/2000: O.D. avril 2005)
  !
  implicit none
  integer, parameter :: dp=kind(1.d0)
  character(len=*), intent(in) :: atom
  real(dp) :: znuc

  select case (trim(atom))
  case ('1H')
     znuc=1
  case ('2H','D')
     znuc=1
  case ('3H','T')
     znuc=1
  case ('H')
     znuc=1
  case ('3He')
     znuc=2
  case ('4He')
     znuc=2
  case ('He')
     znuc=2
 case ('6Li')
     znuc=1
  case ('7Li')
     znuc=1
  case ('Li')
     znuc=1
  case ('9Be')
     znuc=2
  case ('Be')
     znuc=2
  case ('23Na')
     znuc=1
  case ('Na')
     znuc=1
  case ('24Mg')
     znuc=2
  case ('25Mg')
     znuc=2
  case ('26Mg')
     znuc=2
  case ('Mg')
     znuc=2
  case ('27Al')
     znuc=3
  case ('Al')
     znuc=3
  case ('39K')
     znuc=1
  case ('40K')
     znuc=1
  case ('41K')
     znuc=1
  case ('K')
     znuc=1
  case ('40Ca')
     znuc=2
  case ('42Ca')
     znuc=2
  case ('43Ca')
     znuc=2
  case ('44Ca')
     znuc=2
  case ('46Ca')
     znuc=2
  case ('48Ca')
     znuc=2
  case ('Ca')
     znuc=2
  case ('85Rb')
     znuc=1
  case ('87Rb')
     znuc=1
  case ('Rb')
     znuc=1
  case ('84Sr')
     znuc=2
  case ('86Sr')
     znuc=2
  case ('87Sr')
     znuc=2
  case ('88Sr')
     znuc=2
  case ('Sr')
     znuc=2
  case ('133Cs')
     znuc=1
  case ('Cs')
     znuc=1
  case ('130Ba')
     znuc=2
  case ('132Ba')
     znuc=2
  case ('134Ba')
     znuc=2
  case ('135Ba')
     znuc=2
  case ('136Ba')
     znuc=2
  case ('137Ba')
     znuc=2
  case ('138Ba')
     znuc=2
  case ('Ba')
     znuc=2
  case ('Fr')
     znuc=1
  case ('168Yb')
     znuc=2
  case ('170Yb')
     znuc=2
  case ('171Yb')
     znuc=2
  case ('172Yb')
     znuc=2
  case ('173Yb')
     znuc=2
  case ('174Yb')
     znuc=2
  case ('176Yb')
     znuc=2
  case ('Yb')
     znuc=2
  case default
     write(*,*) 'Element ',trim(atom),' unknown'
     znuc = 0.
     !stop
  end select

end function znuc

subroutine find_str(io,str_to_find,istatus,iline,ilen)
  !
  ! Given the number of a unit, find the next appearance of string str_to_find
  !
  ! C. Dion, 10/2000
  !
  ! Input:
  !    io = unit number
  !    str_to_find = string to search for
  !
  ! Output:
  !    istatus = 1 if string is found, 0 if string is not found
  !    iline : location of the string
  !
  ! (Note:  the calling routine is responsible of rewinding the unit if
  !  necessary)
  !
  implicit none
  integer, intent(in) :: io
  integer, intent(out) :: iline,ilen
  character(len=*), intent(in) :: str_to_find
  integer, intent(out) :: istatus
  character(len=len(str_to_find)) :: blockname
  character(len=5) :: fmt
  character(len=2) :: cex2
  character(len=3) :: cex3
  fmt='(A'//cex2(len(str_to_find))//')'
  ilen=len(str_to_find)
  istatus=0
  iline=0

!  write( 6, * ) str_to_find
!  write( 6, * ) len(str_to_find)
  do
     read(io,fmt,end=5) blockname
!  write(6,fmt)blockname
     iline = iline+1
     if (blockname==str_to_find) then
        istatus=1
        return
     end if
  end do

5 continue
  return

end subroutine find_str

subroutine bash_exec(command,dir,name,io,bash_command)
  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: io
  character(len=*) :: command,name,dir,bash_command
  character(len=400) :: command2
  integer*8 :: time_begin,time_end,ir
  real(dp)  :: total_time
  
  call system_clock(count=time_begin, count_rate=ir)
  
  command2=trim(bash_command)//' '//command
  call system(trim(dir)//command2)
      
  call system_clock(count=time_end, count_rate=ir)
  total_time=real(time_end - time_begin,kind=8)/real(ir,kind=8)
  write(io,*) '    exec de '//trim(name)//' en ',total_time,'s'
  
end subroutine bash_exec
