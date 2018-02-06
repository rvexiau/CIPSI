subroutine calcul_pshf(iz,opt,bindir,lb,workdir,threaddir,lw,f_tmp,f_tmp1,prefix,f_pshfcv,noac,rang,bash_command)
!  procedure PSHF,RPOL(RCUT),CVAL,PSHF
!  use iflport
  use io_unit
  implicit none

  integer, parameter :: dp=kind(1.d0)
  integer :: lw,lp,lt,lt1,lpr,lb,lz,kz,i,nrec
  integer :: longu,istatus,iline,l_str,noac,rang
  integer :: iz
  logical :: i_sys
  character(len=*) :: bindir,workdir,threaddir,f_tmp,f_tmp1,opt,f_pshfcv
  character(len=400) :: psnl_fil,ps_molcas_fil,path,command,dir
  character(len=200) :: bash_command
  character(len=10) :: cz,prefix,nrec_str
  
  namelist/psfil/prefix,psnl_fil,ps_molcas_fil

  open(io_input+rang*100,file=trim(threaddir)//'/'//f_tmp)
  read(io_input+rang*100,psfil)                         !lecture namelist
  close(io_input+rang*100)
  lpr=len_trim(prefix)
  lt=len_trim(f_tmp)
  lt1=len_trim(f_tmp1)
  path=workdir(1:lw)//'/'//prefix(1:lpr)
  lp=len_trim(path)
  write(cz,'(i0)')iz
  kz=longu(cz)
  lz=len_trim(cz)
  dir='cd '//trim(threaddir)//';'     !RV 01/16 needed for parallel code : execute binaries in the work directory of the current thread

  !efface les fichiers precedents
  command='rm -f '//prefix(1:lpr)//'*'
  call system(trim(dir)//command)
  command='rm -f BD_*'
  call system(trim(dir)//command)
  command='rm -f TMP*'
  call system(trim(dir)//command)
  command='rm -f fort*'
  call system(trim(dir)//command)
  command='rm -f F10*'
  call system(trim(dir)//command)
  command='rm -f COMB_VCPP'
  call system(trim(dir)//command)
                                        ! premier passage SCF
  command=bindir(1:lb)//'/pshf.exe<'//f_tmp(1:lt)//'>'//prefix(1:lpr-1)//'.pshf'
  call bash_exec(command,dir,'pshf',io_output+rang*100,bash_command)
      open(io_work+rang*100,file=trim(threaddir)//'/'//prefix(1:lpr-1)//'.pshf')     !extraction de NREC et mod des namelists
      call find_str(io_work+rang*100,'nombre de blocs sur la file'&
      &                  ,istatus,iline,l_str)                !du fichier de resultats
      read(io_work+rang*100,*)i,nrec
      close(io_work+rang*100)
      write(nrec_str,'(i0)') nrec
      command= 'sed -i "s# MREC=.*# MREC='//nrec_str//',#" '//f_tmp(1:lt)
      call system(trim(dir)//command)
      
      
                                       ! calcul integrales champ electrique local
  if((opt=='rpol').or.(opt=='RPOL'))then
     ! JD 01/09 new chain for correct long range behaviour following Romain Guerout
     command=bindir(1:lb)//'/rpol.exe<'//f_tmp(1:lt)//'>'//prefix(1:lpr-1)//'.rpol'
     call bash_exec(command,dir,'rpol',io_output+rang*100,bash_command)

     ! correction coeur valence
      command=bindir(1:lb)//'/cval.exe<'//f_tmp(1:lt)//'>'//prefix(1:lpr-1)//'.cval'
      call bash_exec(command,dir,'cval',io_output+rang*100,bash_command)
      
      open(io_work+rang*100,file=trim(threaddir)//'/'//prefix(1:lpr-1)//'.cval')     !extraction de NREC et mod des namelists
      call find_str(io_work+rang*100,'nombre de blocs sur la file'&
      &                  ,istatus,iline,l_str)                !du fichier de resultats
      read(io_work+rang*100,*)i,nrec
      close(io_work+rang*100)
      write(nrec_str,'(i0)') nrec
      command= 'sed -i "s# MREC=.*# MREC='//nrec_str//',#" '//f_tmp(1:lt)
      call system(trim(dir)//command)
      
     ! second passage SCF
      command=bindir(1:lb)//'/pshf.exe<'//f_tmp1(1:lt1)//'>'//prefix(1:lpr-1)//'.pshfcv0'//cz(kz:kz+lz)
      call bash_exec(command,dir,'pshf',io_output+rang*100,bash_command) 

      ! copy and rename files
      call copy_file('info_cv',trim(threaddir)//'/'//prefix(1:lpr),trim(threaddir)//'/'//prefix(1:lpr)//'_rpol_')
      call copy_file('pqrs_cv',trim(threaddir)//'/'//prefix(1:lpr),trim(threaddir)//'/'//prefix(1:lpr)//'_rpol_')
      call copy2_file(trim(threaddir)//'/'//prefix(1:lpr)//'info_cv',trim(threaddir)//'/'//prefix(1:lpr)//'info')
      call copy2_file(trim(threaddir)//'/'//prefix(1:lpr)//'pqrs_cv',trim(threaddir)//'/'//prefix(1:lpr)//'pqrs')

      ! RG 13 mai 2009
      ! save the vcpp created by rpol
      call copy2_file(trim(threaddir)//'/'//prefix(1:lpr)//'vcpp',trim(threaddir)//'/'//'VCPP_RPOL')      
      ! 2nd run
      command=bindir(1:lb)//'/rcut.exe<'//f_tmp(1:lt)//'>'//prefix(1:lpr-1)//'.rcut'
      call bash_exec(command,dir,'rcut',io_output+rang*100,bash_command)

  elseif((opt=='rcut').or.(opt=='RCUT'))then
    command=bindir(1:lb)//'/rcut.exe<'//f_tmp(1:lt)//'>'//prefix(1:lpr-1)//'.rcut'
    call bash_exec(command,dir,'rcut',io_output+rang*100,bash_command)
  endif

  ! 2nd correction coeur valence
  command=bindir(1:lb)//'/cval.exe<'//f_tmp(1:lt)//'>'//prefix(1:lpr-1)//'.cval2'
  call bash_exec(command,dir,'cval',io_output+rang*100,bash_command)
  
  ! last run of pshf
  command=bindir(1:lb)//'/pshf.exe<'//f_tmp1(1:lt1)//'>'//path(1:lp-1)//'.pshfcv'//cz(kz:kz+lz)
  call bash_exec(command,dir,'pshf',io_output+rang*100,bash_command)
  f_pshfcv=prefix(1:lpr-1)//'.pshfcv'//cz(kz:kz+lz)
 
  ! RG 13 mai 2009
  ! This passage in rcut has generated a combined (RPOL+RCUT) vcpp file
  call copy2_file(trim(threaddir)//'/'//'COMB_VCPP',trim(threaddir)//'/'//prefix(1:lpr)//'vcpp')

  open(io_work+rang*100,file=path(1:lp-1)//'.pshfcv'//cz(kz:kz+lz))     !extraction de NOAC (nbre orbitales actives)
  call find_str(io_work+rang*100,' total number of basis functions   ='&
  &                  ,istatus,iline,l_str)                !du fichier de resultats
  rewind io_work+rang*100
  do i=1,iline-1
    read(io_work+rang*100,*)
  enddo
  read(io_work+rang*100,'(36x,i5)')noac
  close(io_work+rang*100)
  
  open(io_work+rang*100,file=trim(threaddir)//'/'//prefix(1:lpr-1)//'.cval2')     !extraction de NREC et mod des namelists
  call find_str(io_work+rang*100,'nombre de blocs sur la file'&
  &                  ,istatus,iline,l_str)                !du fichier de resultats
  read(io_work+rang*100,*)i,nrec
  close(io_work+rang*100)
  write(nrec_str,'(i0)') nrec
  command= 'sed -i "s# NREC=.*# NREC='//nrec_str//',#" '//f_tmp(1:lt)
  call system(trim(dir)//command)
  
  command='rm '//prefix(1:lpr)//'vcpp'
return
end

subroutine calcul_ijkl(iz,bindir,lb,workdir,threaddir,lw,f_tmp,prefix,f_ijkl,itsym,trec,i_bug,rang,bash_command)
!  procedure pour IJKL,FOCK
!  use iflport
  use io_unit
  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: lw,lp,lt,lpr,lb,kz,lz,norb,i,temp,error,i_bug,iz
  integer :: longu,iline,l_str,istatus,rang
  integer, dimension(*) :: itsym
  character(len=*) :: bindir,workdir,threaddir,f_tmp,prefix,f_ijkl
  character(len=200) :: path,command,string,dir,bash_command
  character(len=20) :: trec_str  
  character(len=10) :: cz
  real(dp) :: trec,trec_temp
  
  lpr=len_trim(prefix)
  lt=len_trim(f_tmp)
  path=workdir(1:lw)//'/'//prefix(1:lpr)
  lp=len_trim(path)
  write(cz,'(i0)')iz
  kz=longu(cz)
  lz=len_trim(cz)
  f_ijkl=prefix(1:lpr-1)//'.ijkl'//cz(kz:kz+lz)
  i_bug=0
  trec_temp=trec
  dir='cd '//trim(threaddir)//';'
                                          ! integrales bielectroniques
  ! RV 01/16 loop to catch bug and try to automatically correct trec 
  do i=1,100

    command=bindir(1:lb)//'/ijkl.exe<'//f_tmp(1:lt)//'>'//path(1:lp-1)// &
          & '.ijkl'//cz(kz:kz+lz)
    call bash_exec(command,dir,'ijkl',io_output+rang*100,bash_command)

    open(io_work+rang*100,file=path(1:lp-1)//'.ijkl'//cz(kz:kz+lz))
    call find_str(io_work+rang*100,' !! ERROR in ijkl', istatus,iline,l_str)          
    if(istatus==0) then
      i_bug=1
      exit
    endif  
    read(io_work+rang*100,*) string
    read(io_work+rang*100,'(i2)') error
    close(io_work+rang*100)
    select case(error)
    case (1)
      ! increase trec
      trec_temp=trec_temp*2
      write(io_output+rang*100,*) 'TREC increased to :',trec_temp
      write(trec_str,'(f0.12)') trec_temp
      command= 'sed -i "s# TREC=.*# TREC='//trec_str//',#" '//f_tmp(1:lt)
      call system(trim(dir)//command)
      command= 'rm *f25 *scr *f50 *ijom'
      call system(trim(dir)//command)
      cycle
    case (2)
      ! decrease trec
      trec_temp=trec_temp/3
      write(io_output+rang*100,*) 'TREC decreased to :',trec_temp      
      write(trec_str,'(f0.12)') trec_temp
      command= 'sed -i "s# TREC=.*# TREC='//trec_str//',#" '//f_tmp(1:lt)
      call system(trim(dir)//command)
      command= 'rm *f25 *scr *f50 *ijom'
      call system(trim(dir)//command)
      cycle
    case DEFAULT
     write(io_output+rang*100,*) 'error in ijkl file'
     return
    end select
  enddo
  
  ! reput initial trec value
  write(trec_str,'(f0.12)') trec
  command= 'sed -i "s# TREC=.*# TREC='//trec_str//',#" '//f_tmp(1:lt)
  call system(trim(dir)//command)
                                          ! integrales de Fock
  command=bindir(1:lb)//'/fock.exe<'//f_tmp(1:lt)//'>'//path(1:lp-1)//'.fock'//cz(kz:kz+lz)
  call bash_exec(command,dir,'fock',io_output+rang*100,bash_command)
  
  open(io_work+rang*100,form='unformatted',status='unknown       ',file=trim(threaddir)//'/'//prefix(1:lpr)//'ijkl')
  read(io_work+rang*100,IOSTAT=istatus) temp,norb,temp,temp,(itsym(i),i=1,norb)
  if(istatus==1) then
     write(io_output+rang*100,*) 'error in ijkl file'
     i_bug=0
     return
  endif
  close(io_work+rang*100)
  
  return
end

subroutine calcul_som(iz,bindir,lb,workdir,lw,threaddir,f_som,prefix,rang,bash_command)
!  procedure pour SOM
!  use iflport
  use io_unit
  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: lw,lp,lb,lpr,lsom,kz,lz,i_sys,iz,rang
  integer :: longu
  character(len=*) :: bindir,workdir,prefix,f_som,threaddir
  character(len=200) :: path,command,dir,bash_command
  character(len=10) :: cz

  lpr=len_trim(prefix)
  path=workdir(1:lw)//'/'//prefix(1:lpr)
  lp=len_trim(path)
  lsom=len_trim(f_som)
  dir='cd '//trim(threaddir)//';'

  write(cz,'(i0)')iz
  kz=longu(cz)
  lz=len_trim(cz)
  command='rm -f *_som'
  call system(trim(dir)//command)

  command=bindir(1:lb)//'/som.exe<'//f_som(1:lsom)//'>'//path(1:lp-1)//'.som'//cz(kz:kz+lz)
  call bash_exec(command,dir,'som',io_output+rang*100,bash_command)

  return
end

subroutine calcul_cip(iz,bindir,lb,workdir,threaddir,lw,f_bd,prefix,isymat,sym,f_bdf,metat,netat_pertu,envp_max,det_bounds,imax_s,tau_init,tau_step,e_corr,ncf,nbase,rang,bash_command)
! procedure pour CIP,MOY,BD,
! RV 01/16 rewritten for CIPSI algorithm
! le fichier de donnees F_BD depend maintenant de la symetrie SYM
  use io_unit
  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: lw,lb,lpr,lp,lbd,kz,lz,ls,temp,rang,iz,isymat(6),mw,nw,ntot,ideb,ideb2,imax_s,ind,nbase,counter
  integer :: longu,iline,l_str,istatus,i,j,k,ncf,mcf,netat,metat,metat2,det_bounds(4),netat_pertu
  integer,dimension(:),allocatable :: numtab
  real(dp) :: envp,tau,tau_moy
  real(dp) ::  envp_max,tau_init,tau_step
  real(dp),dimension(metat) :: e_corr
  real(dp),dimension(:),allocatable :: energ2,energ2old,envp_lim,xam
  character(len=*) :: bindir,workdir,threaddir,f_bd,prefix,sym
  character(len=*),intent(out) :: f_bdf
  character(len=100) :: path, f_bd_calc,dir
  character(len=200) :: command,bash_command
  character(len=20) :: tau_str,moy_str,bd_str
  character(len=10) :: cz
  
  integer*8 :: time_begin,time_end,ir
  real(dp)  :: total_time
  
  allocate(energ2(metat))

  f_bd_calc=f_bd
  lpr=len_trim(prefix)
  lbd=len_trim(f_bd_calc)
  path=workdir(1:lw)//'/'//prefix(1:lpr)
  lp=len_trim(path)
  write(cz,'(i0)') iz   !ecriture de z, pour transformation en caractere
  kz=longu(cz)
  lz=len_trim(cz)
  ls=len_trim(sym)
  f_bdf=prefix(1:lpr-1)//'.bd'//sym(1:ls)//cz(kz:kz+lz)
  tau=tau_init
  dir='cd '//trim(threaddir)//';'
  if(isymat(1).ge.0) then  
    moy_str='moysym.exe'
    bd_str='bdsym.exe'
  else
    moy_str='moy.exe'
    bd_str='bd.exe'
  endif

  if(imax_s.gt.0) then
    write(io_output+rang*100,*) 'CIPSI method'
  
    ! generate a reference space     
    do i=1,imax_s
      if(i.gt.1) then
        ! get the number of reference determinant 
        open(io_work+rang*100,file=path(1:lp-1)//'.cip'//sym(1:ls)//cz(kz:kz+lz))
        call find_str(io_work+rang*100,'nombre de determinants sur la file 60', istatus,iline,l_str)                
        if(istatus==0) then
          write(io_output+rang*100,*) 'error in cip file'
          return
        endif
        read(io_work+rang*100,*,IOSTAT=istatus) ncf,metat2
        close(io_work+rang*100)
      endif

      command= 'sed -i "s# ITYPER=.*# ITYPER=-1,#" '//f_bd_calc(1:lbd)
      call system(trim(dir)//command)
      command=bindir(1:lb)//'/cipsym.exe<'//f_bd_calc(1:lbd)//'>'//path(1:lp-1)//'.cip'//sym(1:ls)//cz(kz:kz+lz)
      call bash_exec(command,dir,'cip',io_output+rang*100,bash_command) 
      f_bd_calc=trim(f_bd)//'_loop'
      lbd=len_trim(f_bd_calc)
      
      command= 'sed -i "s# TAU=.*# TAU=0d0,#" '//f_bd_calc(1:lbd)
      call system(trim(dir)//command)
      command=bindir(1:lb)//'/moy.exe<'//f_bd_calc(1:lbd)//'>'//path(1:lp-1)//'.moy'//sym(1:ls)//cz(kz:kz+lz)
      call bash_exec(command,dir,'moy',io_output+rang*100,bash_command)
      command=bindir(1:lb)//'/bd.exe<'//f_bd_calc(1:lbd)//'>'//path(1:lp-1)//'.bd'//sym(1:ls)//cz(kz:kz+lz)
      call bash_exec(command,dir,'bd',io_output+rang*100,bash_command) 
      
   ! compute perturbation   
      if(i.eq.imax_s.or.((ncf.ge.metat).and.(ncf.ge.det_bounds(1)))) then
        tau=0d0
        nbase=ncf
        write(io_output+rang*100,*) '     Final S size : ',ncf        
      else
        write(io_output+rang*100,*) '     Determinants S : ',ncf
      endif
      
      if(ncf.ge.det_bounds(2)) then
        write(*,*) ' S size too big for calculation (restart with different inputs), index : ',iz
        return
      endif
      
      write(tau_str,'(f20.12)') tau
      command= 'sed -i "s# TEST=.*# TEST='//tau_str//',#" '//f_bd_calc(1:lbd)
      call system(trim(dir)//command)
      command= 'sed -i "s# TAU=.*# TAU='//tau_str//',#" '//f_bd_calc(1:lbd)
      call system(trim(dir)//command)
      command= 'sed -i "s# ITYPER=.*# ITYPER=2,#" '//f_bd_calc(1:lbd)
      call system(trim(dir)//command)
      command=bindir(1:lb)//'/cipsym.exe<'//f_bd_calc(1:lbd)//'>'//path(1:lp-1)//'.cip'//sym(1:ls)//cz(kz:kz+lz)
      call bash_exec(command,dir,'cip',io_output+rang*100,bash_command)

      ! Reference matrix with at least minmal number (det_bounds(1)) of determinant 
      if ((ncf.ge.metat).and.(ncf.ge.det_bounds(1))) exit
      tau = tau*tau_step
    enddo
    
    ! get the number of total determinant 
    open(io_work+rang*100,file=path(1:lp-1)//'.cip'//sym(1:ls)//cz(kz:kz+lz))
    call find_str(io_work+rang*100,'nombre de determinants sur la file 60', istatus,iline,l_str)                
    if(istatus==0) then
      write(io_output+rang*100,*) 'error in cip file'
      return
    endif
    read(io_work+rang*100,*,IOSTAT=istatus) ncf,metat2
    close(io_work+rang*100)
    
   ! get the perturbation envp 
    energ2=0d0
    open(io_work+rang*100,file=path(1:lp-1)//'.cip'//sym(1:ls)//cz(kz:kz+lz))
    call find_str(io_work+rang*100,' energies CIPSI envp', istatus,iline,l_str)                
    if(istatus==0) then
      write(io_output+rang*100,*) 'error in cip file'
      return
    endif
    read(io_work+rang*100,*,IOSTAT=istatus) (energ2(j),j=1,metat2)
    close(io_work+rang*100)
    
      
    open(io_work+rang*100,form='unformatted',status='unknown       ',file=trim(threaddir)//'/'//prefix(1:lpr)//'det')
    envp=0d0
    allocate(envp_lim(ncf*metat2),xam(ncf),numtab(ncf))
    ideb=0
    ideb2=0
    do  
    read(io_work+rang*100,end=3487)  mw,nw,(temp,i=1,mw),(temp,i=1,nw),(xam(ideb+i),i=1,nw),(envp_lim(ideb2+i),i=1,nw*metat2)
    ideb2=ideb2+nw*metat2
    ideb=ideb+nw
    enddo 
    3487 continue    

    call system_clock(count=time_begin, count_rate=ir)
    

    call mrgrnk(xam,numtab,ncf)

    mcf=0
    tau=1d0
    counter=0
    do i=ncf,1,-1
      counter=counter+1
      ind=numtab(i)
      do j=1,metat2
        energ2(j)=energ2(j)-envp_lim((ind-1)*metat2+j)
      enddo
      envp=0d0
      do j=1,netat_pertu
        envp=max(envp,dabs(energ2(j)))
      enddo
      if((envp.lt.envp_max).and.(counter.ge.det_bounds(3))) exit
      if(counter.ge.det_bounds(4)) then
        write(*,*) 'too many determinants for calculation (restart with different inputs), index : ',iz,' ; energ. precision reached : ',envp
        return
      endif
    enddo
    tau=xam(ind)
    deallocate(envp_lim,xam)
    close(io_work+rang*100) 
    
    call system_clock(count=time_end, count_rate=ir)
    total_time=real(time_end - time_begin,kind=8)/real(ir,kind=8)
    write(io_output+rang*100,*) '  tri en  ',total_time,'s'
      
    ! Select the size of partition M
    f_bd_calc=trim(f_bd)//'_end'
    lbd=len_trim(f_bd_calc)
      
    write(tau_str,'(f20.12)') tau
    command= 'sed -i "s# TAU=.*# TAU='//tau_str//',#" '//f_bd_calc(1:lbd)
    call system(trim(dir)//command)

    command=bindir(1:lb)//'/'//trim(moy_str)//'<'//f_bd_calc(1:lbd)//'>'//path(1:lp-1)//'.moy'//sym(1:ls)//cz(kz:kz+lz)
    call bash_exec(command,dir,'moy',io_output+rang*100,bash_command)
    command=bindir(1:lb)//'/'//trim(bd_str)//'<'//f_bd_calc(1:lbd)//'>'//path(1:lp-1)//'.bd'//sym(1:ls)//cz(kz:kz+lz)
    call bash_exec(command,dir,'bd',io_output+rang*100,bash_command)
 
  else
    write(io_output+rang*100,*) 'FCI method'
    energ2=0d0
    command=bindir(1:lb)//'/cipsym.exe<'//f_bd_calc(1:lbd)//'>'//path(1:lp-1)//'.cip'//sym(1:ls)//cz(kz:kz+lz)
    call bash_exec(command,dir,'cip',io_output+rang*100,bash_command)
    f_bd_calc=trim(f_bd)//'_end'
    lbd=len_trim(f_bd_calc)
    command=bindir(1:lb)//'/'//trim(moy_str)//'<'//f_bd_calc(1:lbd)//'>'//path(1:lp-1)//'.moy'//sym(1:ls)//cz(kz:kz+lz)
    call bash_exec(command,dir,'moy',io_output+rang*100,bash_command)
    command=bindir(1:lb)//'/'//trim(bd_str)//'<'//f_bd_calc(1:lbd)//'>'//path(1:lp-1)//'.bd'//sym(1:ls)//cz(kz:kz+lz)
    call bash_exec(command,dir,'bd',io_output+rang*100,bash_command)
  endif 
   
   ! write down the number of determinant used
  open(io_work+rang*100,file=path(1:lp-1)//'.bd'//sym(1:ls)//cz(kz:kz+lz))
  call find_str(io_work+rang*100,' ncf,netat,ncper,iget', istatus,iline,l_str)                
  if(istatus==0) then
    write(io_output+rang*100,*) 'error in bd file'
    return
  endif
  read(io_work+rang*100,*) ncf
  close(io_work+rang*100)
  write(io_output+rang*100,*) '     Final computation, determinant :',ncf
  
  if(isymat(1).ge.0) then
    command=bindir(1:lb)//'/spinsym.exe<'//f_bd_calc(1:lbd)//'>'//path(1:lp-1)//'.spin'//sym(1:ls)//cz(kz:kz+lz)
    call bash_exec(command,dir,'spin',io_output+rang*100,bash_command)
  elseif(isymat(1).eq.-99) then
    write(*,*) 'skip spin computation (the sorting will crash)'
  else
    command=bindir(1:lb)//'/spin.exe<'//f_bd_calc(1:lbd)//'>'//path(1:lp-1)//'.spin'//sym(1:ls)//cz(kz:kz+lz)
    call bash_exec(command,dir,'spin',io_output+rang*100,bash_command)
  endif
  
  do j=1,metat
    e_corr(j)=abs(energ2(j))
  enddo
    
  deallocate(energ2)    
  return
end

subroutine calcul_hdiab(iz,bindir,lb,workdir,lw,threaddir,f_hdiab,sym,prefix,rang,bash_command)
!  procedure pour SOM
!  use iflport
  use io_unit
  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: lw,lp,lpr,ldiab,ls,kz,lz,i_sys,lb,iz,rang
  integer :: longu
  character(len=*) :: bindir,workdir,prefix,f_hdiab,sym,threaddir
  character(len=200) :: path,command,dir,bash_command
  character(len=10) :: cz  

  write(cz,'(i0)')iz
  kz=longu(cz)
  lz=len_trim(cz)

  lpr=len_trim(prefix)
  path=workdir(1:lw)//'/'//prefix(1:lpr)
  lp=len_trim(path)
  ls=len_trim(sym)
  ldiab=len_trim(f_hdiab)
  dir='cd '//trim(threaddir)//';'  

  command=bindir(1:lb)//'/hdiab.exe<'//f_hdiab(1:ldiab)//'>'//path(1:lp-1)//'.hdiab'//sym(1:ls)//cz(kz:kz+lz)
  call bash_exec(command,dir,'hdiab',io_output+rang*100,bash_command)

  return
end

subroutine calcul_ciro(iz,bindir,lb,workdir,threaddir,lw,f_ciro,prefix,sym_1,sym_2,f_dip,rang,bash_command)
! procedure pour CIRO
! le fichier de donnees F_CIRO depend maintenant des combinaisons de symetries COM
  use io_unit
  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: lw,lb,lpr,lp,lcir,l1,l2,kz,lz,iz
  integer :: longu,rang
  character(len=*) :: bindir,workdir,threaddir,f_ciro,prefix,f_dip
  character(len=200) :: path,command,dir,bash_command
  character(len=*) :: sym_1,sym_2
  character(len=10) :: cz

  l1=len_trim(sym_1)
  l2=len_trim(sym_2)
  lpr=len_trim(prefix)
  lcir=len_trim(f_ciro)
  path=workdir(1:lw)//'/'//prefix(1:lpr)
  lp=len_trim(path)
  dir='cd '//trim(threaddir)//';'
  
  write(cz,'(i0)')iz
  kz=longu(cz)
  lz=len_trim(cz)
  f_dip=prefix(1:lpr-1)//'.ciro'//sym_1(1:l1)//sym_2(1:l2)//cz(kz:kz+lz)

  command=bindir(1:lb)//'/ciro.exe<'//f_ciro(1:lcir)//'>'//path(1:lp-1)//&
           &'.ciro'//sym_1(1:l1)//sym_2(1:l2)//cz(kz:kz+lz)
  call bash_exec(command,dir,'ciro',io_output+rang*100,bash_command)

  return
end
