! Program to exec the CIPSI chain of Toulouse for atoms and molecules :
! see User's manual for details
!
! traduction en FORTRAN90 du programme fourni par F. Spiegelman en janvier 2003
! O.D., fevrier 2004
! O.D., mars 2005: option de tri singlet/triplet (CASMAT)
! O.D., avril 2005: tri en distance des fichiers sortants
! O.D., juin 2005: tri en symetrie des etats et de moments dipolaires (a partir de CIPMOY)
! O.D., juillet 2005: tri ci-dessus, en donnat les seuils de tri en donnees (seuil_abs,seuil_rel)
!                     elimination des symetries parasites par examen des n_bug premieres composantes
!                     par comparaison a seuil_bug
!
! J. Deiglmayr, april/juin 2006: 
!
! major changes in the format of autocip.in:
!	1. Additional switches
!	   - after the line with logical variables l_cip, l_dip, ...
!	     two lines have been added:
!            l_efield,l_efield_firstpshf	
!	     l_diab,l_calcref
!
!       2. The following lines have been split into "configurations blocks"
!	   labeled by a marker of the form ":name:" where name=symmetries, dipolemoments, 
!          diabatization and electricfield.
!	   For symmetries and dipolemoments the old ordering of autocip7 has been kept
!	   
!          :electricfield:
!	      'direction'         'parallel' or 'orthogonal'
!             nval		  number of field values
!	      E1, E2, ...	  list of nval field-values
!
!	ATTENTION:
!	    If l_efield == t the output is not sorted by symmetries but simply 
!                    written to files ending in _allsym
! O. Dulieu, juin 2011
!
!  important correction of the core-core polarizarion term
!  creation of a function ZNUC to define the charge of the core involved in V_cc
!
! R. Vexiau, april 2016
!   major modifications to allow polyatomic calculations
!   parallelisation done with OpenMP routines
!   allow the use of a perturbative selection of the actve space (CIPSI)
!
! R. Vexiau, may 2017
!   added determination of the sign of the tdm (not 100% accurate)

program auto
 !$ use OMP_LIB
  use io_unit
  use grid_data
  implicit none

  integer :: lwork,nsym,ncom,nzR,lbin,l_m,l_max,lpr
  integer :: nvec,ncf,netat,netat_pertu,ncper,ncore,rang,ntr,det_bounds(4)
  integer :: i,j,k,ki,kj,iz,iz1,iz2,jz,kz,ls,ie,je,if,iv,jv,isz
  integer :: nelac,metat,noac,noac_ref,i_step,loop
  integer :: istatus,iline,l_str,l1,l2,i_sys
  integer :: longu,divise
  integer :: nexmax,imax_s       				! RV 02/16
  integer :: ns,nsmax							! RV 02/16
  integer :: iunt,mcharg,nprint,irest,naxis
  integer :: efield_dirn,efield_nval,iefield   			! JD 04/06 
  integer :: lworkroot,lefield_chars,inerr,err_aloc		! JD 04/06
  integer :: diabnoca,diabnrot					! JD 05/06
  integer :: molecaxis						! JD 05/06
  integer :: i_bug						! RV 02/16
  integer :: axis
  integer,dimension(:),allocatable :: isymat 
  integer,dimension(:),allocatable :: nusym,itsym,jz_s,i_fail,ndet,nbase
  integer,dimension(:),allocatable :: diabnum,diabnumr		!JD 05/06
  integer,dimension(:,:),allocatable :: ietat,jetat,i_sym,iz_f,is_com,is_field
  integer,dimension(:,:,:),allocatable :: is
  integer,dimension(:,:,:,:),allocatable :: ie_s		! RV 02/16 added spin dimension
  integer,dimension(14,14) :: mult_tab   			! RV 01/16  
  real(dp) :: sz
  real(dp) :: m1overm2					! JD 04/06
  real(dp) :: eip,znuc,e_tmp,e_ref0,mass
  real(dp),dimension(:),allocatable :: z_nuc,calfa,e_ip,Ma
  real(dp) :: a_rep,b_rep,c_rep
  real(dp) :: envp_max,tau_init,tau_step,trec			! RV 12/15
  real(dp) :: M1,M2,w                                           ! OD 02/10; reduced mass for center of mass calc.
  real(dp) :: lower_bound
  real(dp) :: phase_test,phase_variation			! RV 05/17
  real(dp),allocatable :: dip_value(:,:,:),phase(:,:)		! RV 05/17
  real(dp),dimension(:),allocatable :: refdist						! JD 05/06
  real(dp),dimension(:) :: efield(3),elec(3)			! JD 04/06
  real(dp),dimension(:),allocatable :: efield_vals		! JD 04/06
  real(dp),dimension(:,:),allocatable :: zval_mat		! JD 05/06
  real(dp),dimension(:),allocatable :: e_ref,e_bd,e_bdtot,e_ccind,e_ccdis,e_ccrep
  real(dp),dimension(:,:),allocatable :: e_bdref,e_bdind,e_bdall,e_thresh
  real(dp),dimension(:,:,:),allocatable :: dx,dy,dz,e_corr
  character(len=180) :: workdir,tmpdir,bindir,cipsidir,f_pshf,f_pshf1,f_tmp,f_tmp1,threaddir,refdir,grid_file,field_dir
  character(len=180) :: f_mono,f_mdip,f_diab,f_som		!JD 05/06
  character(len=180) :: f_sym
  character(len=180) :: workdirroot,efield_chars		      	!JD 04/06 
  character(len=180) :: fmt_bd,fmt_ref,fmt_ind,fmt_all,fmt_com,fmt_dip,fmt_ndet,f_spin
  character(len=1),dimension(:),allocatable :: spin_str
  character(len=3) :: cex3
  character(len=2) :: cex2
  character(len=2),dimension(:),allocatable :: atom
  character(len=4) :: opt_rc, rang_str
  character(len=5) :: opt_mdip
  character(len=5) :: group
  character(len=20) :: grid_method,coordinate 				! RV 03/16
  character(len=10) :: efield_dir				! JD 05/06
  character(len=4),dimension(:),allocatable :: sym
  character(len=4),dimension(:,:),allocatable :: com
  character(len=180),dimension(:),allocatable :: f_bd,f_pshfcv,f_ciro,f_ijkl,f_hdiab
  character(len=180),dimension(:,:),allocatable ::f_bdf,f_dip
  character(len=10) :: prefix
  character(len=10) :: cz
  character(len=80) :: ligne
  character(len=180) :: psnl_fil,ps_molcas_fil,path
  character(len=200) ::	command,dir,bash_command
  character(len=180),dimension(:),allocatable :: f_bd_ref		! JD 05/06
  character(len=1) :: axisname(3)
  logical :: l_mono,l_cip,l_dip,l_rc,l_vcc,l_calc,l_fci
  logical :: ymono,exist
  logical :: l_efield, l_efield_1,l_ref,l_ffield	! JD 04/06
  logical :: l_diab,l_calcref				! JD 05/06
  logical :: l_com                                      ! OD 02/10
  logical :: ysort					! RV 02/16
  logical :: phase_sort				! RV 05/17
  logical :: parity
  logical,dimension(:,:,:,:),allocatable :: l_accept		! RV 05/16
  
   ! timing variables RV 12/15
  integer*8 :: time_begin,time_end,time_begin_iz,time_end_iz,ir
  real(dp)  :: total_time,total_time_iz
  
  namelist/psfil/prefix,psnl_fil,ps_molcas_fil
  namelist/hondo/group,sz,iunt,mcharg,nprint,irest,elec,naxis
  namelist/vpol/calfa,ymono
 
  call timestamp()
  
  nsmax=6 ! max value of 2S+1, if changed isymat and the binary moy/bd need to be modified
  allocate(isymat(nsmax),spin_str(nsmax))
  spin_str(1)='s'
  spin_str(2)='d'
  spin_str(3)='t'
  spin_str(4)='q'
  spin_str(5)='p'
  spin_str(6)='h'  

  data axisname/'x','y','z'/  ! Character Strings for the direction of the electric field
  
  write(*,*) ' BEGIN AUTOCIP **********'
  call system_clock(count=time_begin, count_rate=ir) !progam computation time
  call get_environment_variable("CIPSI_ROOT", cipsidir)  !repertoire de cipsi
  call get_environment_variable("CIPSI_WORK", tmpdir,i,istatus)   ! temporary work folder  
  if(istatus.eq.1) tmpdir='work'                   ! set a default value if the env variable is empty
  call getcwd(workdirroot)				! RV 11/16 : current directory ; results are put in the exec folder   
  lworkroot = len_trim( workdirroot )  
  
  l_com=.false.
  l_ref=.false. 	! JD 05/06 	this flag changes the output format for BD and MOY
			!		in order to save the reference states for SOM and HDIAB 
  ! Unité numéro 5 = entrée standard
  ! Unité numéro 6 = sortie standard
  read( 5, * ) l_calc                !true: calcul +tri; false tri sur fichiers existants
  read( 5, * ) bindir                !repertoire des binaires
  bindir=trim(cipsidir)//bindir  
  lbin = len_trim( bindir )

  read(5,*) bash_command			!RV 01/16 : additional bash command to execute valgrind/time/ etc ...
  read( 5, * )                              !repertoire de travail, can be used by scripts calling autocip
  
  read(5,*) ncore,nelac				! RV 01/16 multicore grid definition
  select case(ncore)
  case (1)
    ncoord=1
  case (2)
    ncoord=1
  case (3)
    ncoord=3
  case (4)
    ncoord=6
  case DEFAULT
    write(*,*) 'invalid number of center'
  end select

  if(mod(nelac,2)==0) then
    isz=0
  else
    isz=1
  endif  
  
  allocate(atom(ncore),z_nuc(ncore),Ma(ncore),calfa(ncore),e_ip(ncore))
  read( 5, * ) (atom(i),i=1,ncore)            !atomes a traiter
  do i=1,ncore
    z_nuc(i)=znuc(atom(i))
    Ma(i)= mass(atom(i), 'a.u.')
  enddo
  
  read( 5, * ) f_pshf                 !fichier de donnees de PSHF
  read( 5, * ) opt_rc                 !option RPOL/RCUT
  read( 5, * ) l_mono,l_max,nvec      		!tri monoelectronique
  read( 5, * ) l_cip,l_dip,l_vcc,a_rep,b_rep,c_rep     !options de calcul, et parmateres de v_rep
  l_efield_1=.true.
  read( 5, * ) l_efield,l_efield_1  	!external electric field
  read( 5, * ) l_diab,l_calcref	   	!l_diab call som and hdiab
					!l_calcref  calculate reference states
					
  open(io_input,file=workdirroot(1:lworkroot)//'/'//f_pshf)    !lecture du NAMELIST HONDO, pour
  read(io_input,psfil)                                  	 !et le prefixe
  read(io_input,vpol)                                   	 !et la polarisabilite
  rewind(io_input)
  close(io_input)
  lpr = len_trim(prefix)
  
  mcharg=-nelac
  do i=1,ncore
    mcharg=mcharg+z_nuc(i)
  enddo

! now look for the different configuration blocks in the input file

  ! reworked by RV april 2016: definition of the coordinate grid
  rewind 5
  call find_str(5,':grid:',istatus,l1,l_str)
  if(istatus==0)then
    write( 6, * ) 'Problem: Configuration-block :grid: not specified'
    stop
  endif
  read(5,*) group,l_com
  if(l_com) then
    read(5,*)
  else
    read(5,*) (Ma(i),i=1,ncore)
  endif  
  read(5,*) grid_method
  read(5,*) coordinate
  allocate(gridmin(ncoord),gridmax(ncoord),gridstep(ncoord),refdist(ncoord))
  read(5,*) (gridmin(i),i=1,ncoord)
  read(5,*) (gridmax(i),i=1,ncoord)
  if((grid_method.eq.'grid').or.(grid_method.eq.'GRID')) then 
    read(5,*) (gridstep(i),i=1,ncoord)
  elseif((grid_method.eq.'LHS').or.(grid_method.eq.'lhs')) then
    read(5,*) nz
  elseif((grid_method.eq.'list').or.(grid_method.eq.'LIST')) then
    read(5,*) grid_file
  else
    write(*,*) 'error in grid method; choose grid, lhs or list'
    stop
  endif
  read(5,*) lower_bound
  
  refdist=-1d0
  ! JD 05/06 read parameters for diabatization procedure 
  if(l_cip.and.l_diab)then
     rewind 5
     call find_str(5,':diabatization:',istatus,l1,l_str)
     if(istatus==0)then
        write( 6, * ) 'Problem: Configuration-block :diabatization: not specified although l_diab=true'
        stop
     endif
     read( 5, * ) f_diab
     read( 5, * ) (refdist(j),j=1,ncoord)
     read( 5, * ) diabnoca
     read( 5, * ) diabnrot
     allocate(diabnum(diabnrot),diabnumr(diabnrot))
     read( 5, * ) diabnum
     read( 5, * ) diabnumr
  endif

! generate the coordinate grid
  call makegrid(ncore,Ma,molecaxis,refdist,group,grid_method,coordinate,lower_bound,grid_file)

  allocate(f_pshfcv(nz+1)) ! JD 05/06 +1 for a possible reference state
  
  if(l_efield)then
     rewind 5
     call find_str(5,':electricfield:',istatus,l1,l_str)
     if(istatus==0)then
	write( 6, * ) 'Problem: Configuration-block :electricfield: not specified although l_efield=true'
	stop
     endif
     if((group.eq.'DINFH').or.(group.eq.'dinfh')) then
       write(6,*) 'Problem : electric field breaks u/g symmetry, change the group '
       stop
     endif
     read( 5, * ) efield_dir
     if((efield_dir(1:1).eq.'x').or.(efield_dir(1:1).eq.'X'))then 
       efield_dirn = 1					
     elseif((efield_dir(1:1).eq.'y').or.(efield_dir(1:1).eq.'Y'))then 
       efield_dirn = 2	
     elseif((efield_dir(1:1).eq.'z').or.(efield_dir(1:1).eq.'Z'))then 	
       efield_dirn = 3	
     else      
	write( 6, * ) "efield_dir must be along one of the cartesian axis "
      	stop  
     endif
     field_dir=workdirroot(1:lworkroot)//'/efield_'//axisname(efield_dirn)//'_'

     read( 5, * ) efield_nval
     if(.not.(efield_nval==0))then 
        allocate(efield_vals(efield_nval))
        read( 5, * ) ,efield_vals
	write( 6, * ) 'electric field-values passed: ',efield_vals
	read(5,*) l_ffield
     else
	read( 5, * ) 
	l_efield = .false.
        efield_nval = 1
	workdir = workdirroot
	lwork = lworkroot
     endif
  else
     l_ffield=.false.
     efield_nval = 1
     workdir = workdirroot
     lwork = lworkroot
  endif
  ! end JD 04/06
  
  ! Condition : si calcul d'interactions de configurations demandé  
  if(l_cip) then
  
    rewind 5
    call find_str(5,':symmetries:',istatus,l1,l_str)
    if(istatus==0)then
      write( 6, * ) 'Problem: Configuration-block :symmetries: not specified'
      stop
    endif
 
    read( 5, * ) f_sym                !fichier de donnees type de CIP,MOY,BD. Ex. : carbsym.dat
    read( 5, * ) nsym,metat       !nbre de symetries a traiter,nbre d'etats voulus
    
    if(l_diab) allocate(f_hdiab(nsym),f_bd_ref(nsym))    
    allocate(sym(nsym),STAT=err_aloc)
    allocate(ndet(nsym),nbase(nsym))
    
    read( 5, * ) sym                  !symetries a traiter
    read( 5, * ) (isymat(i),i=1,6) !spin to compute (RV 02/16)    
    if(isymat(1).eq.-99) ysort=.false.  ! spin/sort binary will not be called
    allocate(f_bd(nsym),nusym(nsym),f_bdf(nsym,nz+1),f_ijkl(nz+1))
    allocate(e_corr(metat,nsym,nz+1))
    
    rewind 5
    call find_str(5,':method:',istatus,l1,l_str)
    if(istatus==0)then
      write( 6, * ) 'Problem: Configuration-block :method: not specified'
      stop
    endif
    read(5,*) l_fci
    if(l_fci) then
      imax_s=0
    else
      imax_s=1000		! maximum numbre of iteration
    endif
    read(5,*) noac_ref,nexmax                                  ! number of active orbital for the initial reference CAS
    parity = mod(nelac,2)
    if(parity.and.nexmax.ge.0) nexmax=nexmax+1		! take into account the virtual electron added by the method
 
    read( 5, * ) netat_pertu,envp_max ! n_etat_envp, envp maximum apres calcul
    read( 5, * ) tau_init,tau_step ! nbre de det. de reference minimun; tau initial, tau step
    read( 5, * ) (det_bounds(i),i=1,4)			! minimal and maximal number of determinants for reference space ; min and max number for total space
    read( 5, * ) ysort					! autosort parasite state by looking at monoexcited determinant
    read( 5, * ) trec
    
    
    ! Condition : calcul des moments dipolaires de transition demandé
    if(l_dip)then                 !etape 3
    
      rewind 5
      call find_str(5,':dipolemoments:',istatus,l1,l_str)
      if(istatus==0)then
         write( 6, * ) 'Problem: Configuration-block :dipolemoments: not specified although l_dip=true'
         stop
      endif
      
      read( 5, * ) f_mdip             !fichier de donnees type de CIRO
      read( 5, * ) phase_sort,phase_variation ! attempt to correct the sign of tdm      
      read( 5, * ) ncom               !nombre de combinaisons de transition
      if(ncore.ne.2) phase_sort=.false.		! only work for dimers

      allocate(f_ciro(ncom),com(ncom,2),is_com(ncom,2),f_dip(ncom,nz))       
      do i=1,ncom                 !combinaisons de symetries        
        read( 5, * ) com(i,1),com(i,2)  !nbre d'etats dans chaque symetrie
        ntr=metat*metat             !nbre total de transitions
      enddo
        
      allocate(ietat(ntr,ncom),jetat(ntr,ncom))
        
      do i=1,ncom 
        j=0
        do ki=1,metat
          do kj=1,metat
            j=j+1
            ietat(j,i)=ki                   !liste des combinaisons pour CIRO
            jetat(j,i)=kj
          enddo
        enddo
      enddo

      do i=1,ncom                 !verification des combinaisons par rapport
        do j=1,2                  !aux symetries donnees plus haut
          do k=1,nsym
            if(com(i,j)==sym(k))then
              is_com(i,j)=k       !numero de la symetrie identifiee, dans la liste calculee
              istatus=1
            endif
          enddo
          if(istatus==0)then
            write( 6, * ) 'erreur de symetrie, combinaison ',i
            stop
          endif
        enddo
      enddo

    ! Fin condition demande de calcul des moments de transition  
    endif
  ! Fin condition demande de calcul d'interactions de configurations
  else
  ! unused values will be set to zero
    nsym=1
    allocate(ndet(nsym),nbase(nsym))
  endif
  
  ! array that will indicate if bd was successful
  if(l_ffield) then
    allocate(is_field(nsmax,nsym))
    is_field=0
    allocate(l_accept(nz+1,efield_nval,nsmax,nsym))
    l_accept=.false.
  endif  

  write(*,*) "Starting computing main loop *****"
  
! ************************************************************************************************************************  
! ************************************************************************************************************************  
!    START COMPUTATION
! ************************************************************************************************************************  
! ************************************************************************************************************************  

! JD 04/06
! Boucle champ électrique. Englobe TOUT le calcul.
  do iefield=1,efield_nval
    efield(1)=0.0
    efield(2)=0.0
    efield(3)=0.0

    if(l_efield)then 
	! Create a new directory for every electric field value. 
	! Therefore the changes in the rest of the code are minimal

	write(efield_chars,'(E9.2)')efield_vals(iefield)
	efield_chars=trim(adjustl(efield_chars))
  	lefield_chars=len_trim(efield_chars)

	write( 6, * ) 
	write( 6, * ) '##################################'
	write( 6, * ) 'Calcul pour e-field=',efield_chars(1:lefield_chars)
	write( 6, * ) '##################################'
	write( 6, * ) 
	
	efield(efield_dirn)=efield_vals(iefield)

	workdir = trim(field_dir)//efield_chars(1:lefield_chars)
	lwork = len_trim(workdir)

  	command='mkdir -p '//workdir(1:lwork)  ! create new directory
  	call system(command)
  	command='cp '//workdirroot(1:lworkroot)//'/'//prefix(1:lpr-1)//'*.dat '//workdir(1:lwork) ! copy date files
  	call system(command)
    endif
! end JD 04/06

! Attention: INSIDE LOOP for e-field until END!!!

  rang=0
  command='mkdir -p '//workdir(1:lwork)//'/work'
  call system(command)
  
! if l_calc=False, set the name for the output files and skip the calculation loop  
  if (.not.l_calc) then
    do iz=1,nz+1
      write(cz,'(i0)')iz
      f_pshfcv(iz)=prefix(1:lpr-1)//'.pshfcv'//trim(adjustl(cz))
    enddo
  else
  
! Calculate reference states for diabatization if requested
  if(l_diab)then

      write( 6, * ) 'calculation of reference states at : ',refdist

      f_tmp='temp.dat'
      f_tmp1='temp1.dat'
      refdir = trim(tmpdir)//'/workref'
      command='mkdir -p '//trim(refdir)
      call system(command)
      open(io_output+rang*100, file=workdir(1:lwork)//'/work/work_ref.stdout')
      
      command='cp '//trim(f_pshf)//' '//trim(refdir)//'/'
      call system(command)
      if(l_cip) then
        command='cp '//trim(f_sym)//' '//trim(refdir)//'/'
        call system(command)
      endif      
      call makedata_pshf(opt_rc,refdir,f_pshf,f_tmp,f_tmp1,prefix,ncore,nelac,mcharg,nz+1,&
		&group,l_efield,l_efield_1,efield,trec,rang) 
      if(l_cip)call sym_check(group,nsym,sym,nusym,mult_tab)

      call calcul_pshf(nz+1,opt_rc,bindir,lbin,workdir,refdir,lwork,f_tmp,f_tmp1,prefix,f_pshfcv(nz+1),noac,rang,bash_command )
      allocate(itsym(noac))
      call calcul_ijkl(nz+1,bindir,lbin,workdir,refdir,lwork,f_tmp,prefix,f_ijkl(nz+1),itsym,trec,i_bug,rang,bash_command)

      l_ref=.true.
      do i=1,nsym
        call makedata_cip(refdir,f_sym,f_bd_ref(i),prefix,group,isymat,sym(i),nusym(i),itsym,mult_tab,isz,nelac,metat,noac,noac_ref,nexmax,l_ref,rang) 
        call calcul_cip(nz+1,bindir,lbin,workdir,refdir,lwork,f_bd_ref(i),prefix,isymat,sym(i),f_bdf(i,nz+1),metat,netat_pertu,envp_max,det_bounds,imax_s,tau_init,tau_step,e_corr(:,i,nz+1),ndet(i),nbase(i),rang,bash_command)
      enddo
      deallocate(itsym)
      l_ref=.false.
      close(io_output+rang*100)
  endif

  ! RV 12/15 time of each geometry


  !RV 01/16 OpenMP 
  !$OMP PARALLEL DEFAULT(PRIVATE), &
  !$OMP& FIRSTPRIVATE(f_pshf,f_sym,f_diab,f_mdip,opt_rc,prefix,lpr,group,l_efield,l_efield_1,efield,molecaxis,trec,isz), &
  !$OMP& FIRSTPRIVATE(ncore,isymat,nsym,nusym,bindir,lbin,workdir,tmpdir,field_dir,lwork,nelac,mcharg,metat,netat_pertu,envp_max,det_bounds,tau_init,tau_step,noac_ref,nexmax,imax_s), &
  !$OMP& FIRSTPRIVATE(ncom,com,ntr,ietat,jetat,calfa), &
  !$OMP& FIRSTPRIVATE(diabnoca,diabnrot,diabnum,diabnumr,refdir), &
  !$OMP& FIRSTPRIVATE(l_ref,l_cip,l_dip,l_diab,fmt_ndet,bash_command), &
  !$OMP& SHARED(e_corr,f_pshfcv,coord,coord_dist,coord_cart,Ma,nz,ncoord,sym,l_accept,nsmax), &
  !$OMP& SHARED(io_ndet,io_input,io_namelist,io_namelist2,io_output,io_work)
  
  !$OMP SINGLE
    !$ write(*,*) 'Computation with OMP_THREAD =',OMP_GET_NUM_THREADS()
    open(io_ndet, file = workdir( 1 : lwork )//'/n_time.dat' )
    write(io_ndet,*) '#  number of symmetry : ',nsym
    write(io_ndet,'(a100)') '# distance / cpu num. / number of det (space S and total) for each sym / computation time (s)'
  !$OMP END SINGLE
    write( fmt_ndet, '(a16,i0,a21)' )  '(i5,2x,a3,i0,2x,', nsym,'(i8,2x,i8,2x),f20.2)'
  
!  create a working directory for each thread
  !$ rang= OMP_GET_THREAD_NUM()
   write(rang_str,'(i0)') rang
   threaddir = trim(tmpdir)//'/work'//trim(rang_str)  
   command='mkdir -p '//trim(threaddir)
   call system(command)
   command='cp '//trim(f_pshf)//' '//trim(threaddir)//'/'
   call system(command)
   if(l_cip) then
     command='cp '//trim(f_sym)//' '//trim(threaddir)//'/'
     call system(command)
     if(l_dip) then
       command='cp '//trim(f_mdip)//' '//trim(threaddir) //'/'
       call system(command)
     endif
   endif
   if(l_diab) then
     command='cp '//trim(f_diab)//' '//trim(threaddir) //'/'
     call system(command)
   endif   
  open(io_output+rang*100, file=workdir(1:lwork)//'/work/work'//trim(rang_str)//'.stdout')
  write( fmt_all, '(a5,i0,a11)' ) '(a13,', ncoord,'(3x,f15.4))'
  
  ! #### Boucle sur les distances internucléaires ####
  loop=0
  !$OMP DO SCHEDULE(DYNAMIC)
  do iz = 1, nz

    call system_clock(count=time_begin_iz, count_rate=ir)
    loop=loop+1
    ndet = 0 ! initialisation
    nbase = -1

    write( 6, fmt_all ) 'gridpoint = ', (coord(j,iz),j=1,ncoord)
    write(io_output+rang*100, fmt_all ) 'gridpoint = ', (coord(j,iz),j=1,ncoord)
 
    f_tmp='temp.dat'
    f_tmp1='temp1.dat'
    call makedata_pshf(opt_rc,threaddir,f_pshf,f_tmp,f_tmp1,prefix,ncore,nelac,mcharg,iz,&
		&group,l_efield,l_efield_1,efield,trec,rang) !  donnees pshf -- JD 04/06 parameters extended
    if(l_cip) call sym_check(group,nsym,sym,nusym,mult_tab)

    call calcul_pshf( iz,opt_rc,bindir,lbin,workdir,threaddir,lwork,f_tmp,f_tmp1,prefix,f_pshfcv(iz),noac,rang,bash_command )

    ! Si demande de calcul d'interaction de configurations
    ! Calcul bielectroniques + IC totale
    if(l_cip) then
      if( allocated(itsym) ) deallocate( itsym )				! RV 12/15 orbital symmetries computed by ijkl, needed by makedata_cip
      allocate(itsym(noac))
      call calcul_ijkl(iz,bindir,lbin,workdir,threaddir,lwork,f_tmp,prefix,f_ijkl(iz),itsym,trec,i_bug,rang,bash_command)
      if (i_bug==0) then
        write(*,*) 'error in ijkl, geometry ',iz,' skipped'
        deallocate(itsym)
        loop=loop-1
        cycle
      endif  
         
      if( loop == 1 ) then
        if(l_diab) call makedata_som(threaddir,prefix,f_diab,f_som,noac)
        do i=1, nsym
          call makedata_cip(threaddir,f_sym,f_bd(i),prefix,group,isymat,sym(i),nusym(i),itsym,mult_tab,isz,nelac,metat,noac,noac_ref,nexmax,l_ref,rang) 
          if(l_diab) call makedata_hdiab(threaddir,f_diab,f_hdiab(i),sym(i),diabnoca,diabnrot,diabnum,diabnumr)
        enddo
      endif  
      deallocate(itsym)

      do i = 1, nsym
        call calcul_cip(iz,bindir,lbin,workdir,threaddir,lwork,f_bd(i),prefix,isymat,sym(i),f_bdf(i,iz),metat,netat_pertu,envp_max,det_bounds,imax_s,tau_init,tau_step,e_corr(:,i,iz),ndet(i),nbase(i),rang,bash_command)
        if(l_diab) then
          call fetch_refom(prefix,lpr,threaddir,refdir,f_diab)
          call fetch_refbd(prefix,lpr,threaddir,refdir,sym(i))        
          call calcul_som(iz,bindir,lbin,workdir,lwork,threaddir,f_som,prefix,rang,bash_command)
          call calcul_hdiab(iz,bindir,lbin,workdir,lwork,threaddir,f_hdiab(i),sym(i),prefix,rang,bash_command)	
        endif
      enddo

      ! Si demande de calcul de moment dipolaire
      if(l_dip) then
        do i=1,ncom
          if(loop==1) call makedata_ciro(threaddir,f_mdip,f_ciro(i),prefix,com(i,1),com(i,2),calfa,ncore,ntr,ietat(:,i),jetat(:,i),rang)
          call calcul_ciro(iz,bindir,lbin,workdir,threaddir,lwork,f_ciro(i),prefix,com(i,1),com(i,2),f_dip(i,iz),rang,bash_command)
        enddo
      endif
    ! Fin de condition sur le calcul d'interaction de configuration
    endif
    
    ! time of the step
    call system_clock(count=time_end_iz, count_rate=ir)
    total_time_iz=real(time_end_iz - time_begin_iz,kind=8)/real(ir,kind=8)
    if (.not.l_cip) then 
      ndet=0
      nbase=-1
    endif  
    !$OMP CRITICAL
    write(io_ndet,fmt_ndet) iz,'cpu',rang, (nbase(i),ndet(i),i=1,nsym),total_time_iz
    !$OMP END CRITICAL
  enddo
  !$OMP END DO
  ! #### Fin boucle sur les distances internucléaires ####
  
  !$OMP SINGLE
  close(io_ndet)
  !$OMP END SINGLE  
  
  dir='cd '//trim(threaddir)//';'
  lpr = len_trim( prefix )
  command='rm -f '//prefix(1:lpr)//'*'  !efface les fichiers binaires
  command=trim(dir)//command  
  call system(command)
  command='rm -f BD_*'  !efface les fichiers binaires
  command=trim(dir)//command  
  call system(command)
  command='rm -f TMP*'  !efface les fichiers binaires
  command=trim(dir)//command
  call system(command)
  command='rm -f fort*'  !efface les fichiers binaires
  command=trim(dir)//command
  call system(command)
  command='rm -f F10*'  !efface les fichiers binaires
  command=trim(dir)//command
  call system(command)
  command='rm -f COMB_VCPP'  ! RG 13 mai 2009
  command=trim(dir)//command  
  call system(command)
  
  close(io_output+rang*100)
  !$OMP END PARALLEL
  endif
  
  ! start sorting the results
  write(*,*) "Main loop done, sorting the results"
  
  if(l_diab)then
    nzR=nz+1	! sort also reference state
  else
    nzR=nz
  endif

  !sauve les energies monoelectroniques
  if( l_mono ) then
    lpr=len_trim(prefix)
    call tri_mono(lwork,workdir(1:lwork),l_max,nvec,ncore,f_pshfcv(1:nzR),nzR)
    !end JD 04/06
  endif


  ! Si l'interaction de configurations n'est pas demandée, 
  ! on passe à la valeur de champ électrique suivante
  if(.not.l_cip) then 
    command='gzip -f '//workdir(1:lwork)//'/*pshfcv*'
    call system(command)
    command='mkdir -p '//workdir(1:lwork)//'/Results/'
    call system(command)
    command='mv -f '//workdir(1:lwork)//'/*.gz '//workdir(1:lwork)//'/Results/'
    call system(command)
    cycle
  endif

  
   !energies totales
  ! ############ 
  ! ############ Début tri, génération des fichiers de potentiel, moments dipolaires, etc.
  ! ############
  
  if(allocated(e_ref)) deallocate(e_ref,e_ccind)
  allocate(e_ref(nzR),e_ccind(nzR))

  ! Si effet de volume (coeur-coeur) et énergie de dispersion (coeur-coeur) demandés
  if(l_vcc) then
     if(allocated(e_ccdis)) deallocate(e_ccdis,e_ccrep)
     allocate(e_ccdis(nzR),e_ccrep(nzR))
  endif

  
  if(allocated(is)) deallocate(is,iz_f,jz_s)
  allocate(is(nzR,nsmax,nsym))
  allocate(iz_f(nzR,nsym))
  allocate(jz_s(nsym))
  is = 0

  if(.not.l_calc) call system('rm -f '//workdir(1:lwork)//'/v_*')
  if(.not.l_calc) call system('rm -f '//workdir(1:lwork)//'/vpur_*')
  if(.not.l_calc) call system('rm -f '//workdir(1:lwork)//'/vccind_*')
  if(.not.l_calc) call system('rm -f '//workdir(1:lwork)//'/vccall_*')

  ! #### (tri) Boucle sur les symétries ####
  sym_energy:do i = 1, nsym

    !bare potential file (original) 
    inquire( file = workdir( 1 : lwork )//'/vpur_'//sym(i), EXIST = exist )
    if(exist)then
      open( 40+i, file = workdir( 1 : lwork )//'/vpur_'//sym(i), status = 'old', position = 'append' )
    else
      open( 40+i, file = workdir( 1 : lwork )//'/vpur_'//sym(i), status = 'new' )
    endif

    jz=0
    jz_s(i)=0

    ! #### (tri) Boucle sur les distances ####
    z_energy:do iz = 1, nzR
      ! DB 04/2014 : réécriture plus légère des lignes précédentes en tirant profit des fonctions intrinsèques du fortran90
      write(cz,'(i0)') iz
      f_bdf(i,iz)=prefix(1:lpr-1)//'.bd'//trim(sym(i))//trim(adjustl(cz))

      ! Ouverture du fichier "bd" dont le nom vient d'être écrit
      open(50+iz,file=workdir(1:lwork)//'/'//f_bdf(i,iz))
      
      ! Où est l'énergie du determinant de reference dans le fichier "bd" ?
      call find_str( 50+iz, '  energie de reference', istatus, iline, l_str )
      if( istatus == 0 ) then
        write( 6, * )  'Missing ENERGIES REF card in BD file'
        go to 11
      endif
      read( 50+iz, * ) e_ref0

      ! Quel est le nombre d'etats calculés dans le fichier "bd" ?
      call find_str(50+iz,' ncf,netat,ncper,iget', istatus,iline,l_str)                
      if(istatus==0) then
        write( 6, * )  'Missing NETAT card in BD file'
        go to 11
      endif
      read( 50+iz , * ) ncf, netat, ncper
      
      ! DB 04/2014 : En fortran 95, le descripteur "i0" écrit un nombre entier sur le nombre de caractères qui convient.
      ! Ex. Le nombre 22 est écrit sur 2 caractères, le nombre 100 est écrit sur 3 caractères.
      write(fmt_bd,'(a1,i0,a6,i0,a7)')  '(',ncoord,'f12.4,', netat+1,'f20.12)'
      write(fmt_ref,'(a1,i0,a6,i0,a7)') '(',ncoord,'f12.4,', netat+2,'f20.12)'
      write(fmt_ind,'(a1,i0,a6,i0,a7)') '(',ncoord,'f12.4,', netat+2,'f20.12)'
      write(fmt_all,'(a1,i0,a6,i0,a7)') '(',ncoord,'f12.4,', netat+3,'f20.12)'

      if(.not.allocated(e_bd)) allocate(e_bd(netat))
      if(.not.allocated(e_bdtot)) allocate(e_bdtot(netat))
      if(.not.allocated(e_bdref)) allocate(e_bdref(netat,nzR))
      if(.not.allocated(e_bdind)) allocate(e_bdind(netat,nzR))
      if(.not.allocated(e_thresh)) allocate(e_thresh(netat,nzR))
      if((l_vcc).and.(.not.allocated(e_bdall))) allocate(e_bdall(netat,nzR))

      ! RV 2015 : lecture direct du résultat final dans bd
      call find_str(50+iz,'     var          var + epolnuc', &
      &               istatus,iline,l_str)
      if(istatus==0) then
        write( 6, * )  'Missing ENERGIES VAR card in BD file'
        go to 11
      endif

      do iv = 1,netat                      
        read( 50+iz , * ) e_bd(iv), e_bdtot(iv)
      enddo 
      ! Ecriture de e_bd (énergies propres telles qu'écrites dans le fichier "bd") dans 'vpur_[...]'
      write(40+i,fmt_bd) (coord(je,iz),je=1,ncoord), (e_bd(ie),ie=1,netat)

      rewind(50+iz)
      call find_str( 50+iz, ' vecteurs propres', istatus, iline, l_str )
      if(istatus==0) then
        write( 6, * )  'Missing VECTEURS PROPRES card in BD file'
        go to 11
      endif

      if(.not.allocated(ie_s)) then 
        allocate(ie_s(netat,nzR,nsmax,nsym))
        ie_s=0
      endif  
      if(.not.allocated(i_sym)) allocate(i_sym(netat,nzR))

      ! Si on n'a pas trouvé les énergies ou les vecteurs propres, on ne passe pas par ici, cf. les "go to 11" un peu plus haut.
      ! C'est pourquoi jz compte les bons fichiers "bd".
      jz = jz + 1
      jz_s( i ) = jz                          !nbre de valeurs de z ou BD est valide
      iz_f( jz, i ) = iz                      ! Pour chaque symétrie (indice i), on enregistre les valeurs de z (indice iz) ou BD est valide
      kz = iz_f( jz, i )
      
      e_ref( jz ) = e_ref0
      e_ccind( jz ) = -0.5 * ( z_nuc(2)*z_nuc(2)*calfa(1)+z_nuc(1)*z_nuc(1)*calfa(2) ) / coord(1, kz )**4
      
      if( l_vcc ) then
        e_ccdis(jz)=-1.5*calfa(1)*calfa(2)*e_ip(1)*e_ip(2)/(e_ip(1)+e_ip(2))/coord(1,kz)**6
        e_ccrep(jz)=a_rep*(coord(1,kz)**b_rep)*exp(-c_rep*coord(1,kz))
      endif


! RV 01/16 spin sorted by reading the output of the new spin binary
      allocate(i_fail(netat))
      write(cz,'(i0)') iz
      f_spin=workdir(1:lwork)//'/'//prefix(1:lpr-1)//'.spin'//trim(sym(i))//trim(adjustl(cz))
      call sym_vec( netat,f_spin,nelac,ysort,group,isymat(1),i_sym(:,jz),i_fail)

      do ie=1,netat                     
        if(i_fail(ie)==0)then                            
          ns=i_sym(ie,jz)
          is(jz,ns,i)=is(jz,ns,i)+1
          ie_s(is(jz,ns,i),jz,ns,i)=ie
          e_bdref(ie,jz)=e_bd(ie)
          e_bdind(ie,jz)=e_bdtot(ie)
          e_thresh(ie,jz)=e_corr(ie,i,jz)
          if(l_vcc)e_bdall(ie,jz)=e_bdind(ie,jz)+e_ccdis(jz)+e_ccrep(jz)
        endif
      enddo
      deallocate(i_fail)
11    close(50+iz)
      if(allocated(e_bd))deallocate(e_bd)
      if(allocated(e_bdtot))deallocate(e_bdtot)
    enddo z_energy   !index iz

    do ns=1,nsmax       
      if(maxval(is(:,ns,i)).eq.0) cycle ! no results for this spin sym

      ! Open output files

      inquire(file=workdir(1:lwork)//'/v_'//spin_str(ns)//sym(i),EXIST=exist)          !bare potential file (filtered)
      if(exist)then
        open(10+i,file=workdir(1:lwork)//'/v_'//spin_str(ns)//sym(i),status='old',position='append')
      else
        open(10+i,file=workdir(1:lwork)//'/v_'//spin_str(ns)//sym(i),status='new')
      endif
      inquire(file=workdir(1:lwork)//'/vccind_'//spin_str(ns)//sym(i),EXIST=exist)      !potential file with induction term
      if(exist)then
        open(20+i,file=workdir(1:lwork)//'/vccind_'//spin_str(ns)//sym(i),status='old',position='append')
      else
        open(20+i,file=workdir(1:lwork)//'/vccind_'//spin_str(ns)//sym(i),status='new')
      endif
      inquire(file=workdir(1:lwork)//'/vthresh_'//spin_str(ns)//sym(i),EXIST=exist)      !potential file with induction term
      if(exist)then
        open(80+i,file=workdir(1:lwork)//'/vthresh_'//spin_str(ns)//sym(i),status='old',position='append')
      else
        open(80+i,file=workdir(1:lwork)//'/vthresh_'//spin_str(ns)//sym(i),status='new')
      endif
      if(l_vcc)then
        inquire(file=workdir(1:lwork)//'/vccall_'//spin_str(ns)//sym(i),EXIST=exist)      !potential file with all core-core terms
        if(exist)then
          open(30+i,file=workdir(1:lwork)//'/vccall_'//spin_str(ns)//sym(i),status='old',position='append')
        else
          open(30+i,file=workdir(1:lwork)//'/vccall_'//spin_str(ns)//sym(i),status='new')
        endif
        do j=1,ncore
          e_ip(j)=eip(atom(j))
        enddo
      endif

    z_energy_w: do iz=1,jz_s(i)
        write( 10+i, fmt_ref ) (coord(je,iz_f(iz,i)),je=1,ncoord), e_ref(iz), (e_bdref(ie_s(je,iz,ns,i),iz),je=1,is(iz,ns,i)),&
       &                                          (0.,je=1,maxval(is(:,ns,i))-is(iz,ns,i))
        write(20+i,fmt_ind)(coord(je,iz_f(iz,i)),je=1,ncoord),(e_bdind(ie_s(je,iz,ns,i),iz),je=1,is(iz,ns,i)),&
       &     (0.,je=1,maxval(is(:,ns,i))-is(iz,ns,i)) 
        if(l_vcc)&
       &      write(30+i,fmt_all)(coord(je,iz_f(iz,i)),je=1,ncoord),e_ccrep(iz),(e_bdall(ie_s(je,iz,ns,i),iz),je=1,is(iz,ns,i)),&
       &     (0.,je=1,maxval(is(:,ns,i))-is(iz,ns,i))
        write(80+i,fmt_ind)(coord(je,iz_f(iz,i)),je=1,ncoord),envp_max,(e_thresh(ie_s(je,iz,ns,i),iz),je=1,is(iz,ns,i)),&
       &     (0.,je=1,maxval(is(:,ns,i))-is(iz,ns,i))
       
        if(l_ffield) l_accept(iz_f(iz,i),iefield,ns,i)=.true.        
        
      enddo z_energy_w
      
      if(l_ffield) then
        if(iefield.eq.1) then
          is_field(ns,i)=maxval(is(:,ns,i))
        else        
          is_field(ns,i)=min(is_field(ns,i),maxval(is(:,ns,i)))
        endif
      endif
      
      close(10+i)
      close(20+i)
      close(80+i)
      if(l_vcc) close(30+i)  
    enddo  !spin ns=(2*S+1)
    
    if(allocated(e_bdref))deallocate(e_bdref)
    if(allocated(e_bdind))deallocate(e_bdind)
    if(allocated(e_thresh))deallocate(e_thresh)
    if(l_vcc.and.allocated(e_bdall))deallocate(e_bdall)
    
    close(40+i)
  enddo  sym_energy  !index i

  command='gzip -f '//workdir(1:lwork)//'/*.bd*'
  call system(command)
  command='gzip -f '//workdir(1:lwork)//'/*.cip*'
  call system(command)
  command='gzip -f '//workdir(1:lwork)//'/*.moy*'
  call system(command)
  command='gzip -f '//workdir(1:lwork)//'/*ijkl*'
  call system(command)
  command='gzip -f '//workdir(1:lwork)//'/*pshfcv*'
  call system(command)
  command='gzip -f '//workdir(1:lwork)//'/*spin*'
  call system(command)
  command='gzip -f '//workdir(1:lwork)//'/*fock*'
  call system(command)
                                    !moments dipolaires
  if(.not.l_dip) then 
    command='mkdir -p '//workdir(1:lwork)//'/Results/'
    call system(command)
    command='mv -f '//workdir(1:lwork)//'/*.gz '//workdir(1:lwork)//'/Results/'
    call system(command)
    cycle
  endif
 ! if(l_efield)  then
 !   command='mkdir -p '//workdir(1:lwork)//'/Results/'
 !   call system(command)
 !   command='mv -f '//workdir(1:lwork)//'/*.gz '//workdir(1:lwork)//'Results/'
 !   call system(command)
 !   cycle
 ! endif	! JD05/06 the remaining procedure was not adapted to calculations within an
			! electric field. Dipole moments have to be extracted directly from the ciro files

  if(allocated(dx))deallocate(dx,dy,dz)

  allocate(dx(0:metat,metat,nz))
  allocate(dy(0:metat,metat,nz))
  allocate(dz(0:metat,metat,nz))
  write(fmt_com,'(a4,i0,a5,i0,a12)')  '(a1,',ncoord,'(12x),', metat,'(4x,2i3,4x))'
  write(fmt_dip,'(a1,i0,a6,i0,a6)')  '(',ncoord,'f12.4,', metat,'e14.6)'

  if(.not.l_calc)call system('rm -f '//workdir(1:lwork)//'/dx*')
  if(.not.l_calc)call system('rm -f '//workdir(1:lwork)//'/dy*')
  if(.not.l_calc)call system('rm -f '//workdir(1:lwork)//'/dz*')

  com_dip:do i=1,ncom    
    l1=len_trim(com(i,1))
    l2=len_trim(com(i,2))

    z_dip_r:do iz=1,nz
      write(cz,'(i0)') iz
      f_dip(i,iz)=prefix(1:lpr-1)//'.ciro'//com(i,1)(1:l1)//com(i,2)(1:l2)//trim(adjustl(cz))
      open(70+iz,file=workdir(1:lwork)//'/'//f_dip(i,iz))
      call find_str(70+iz,&
  &        ' etats         dmx       dmy       dmz     te(cm-1)   te(eV)     fe', &
  &        istatus,iline,l_str)
      if(istatus==0) then
        do ki=0,metat
          do kj=1,metat
            dx(ki,kj,iz)=0d0
            dy(ki,kj,iz)=0d0
            dz(ki,kj,iz)=0d0
          enddo
        enddo
      else
        do ki=1,metat
          do kj=1,metat
            read(70+iz,*)k,k,dx(ki,kj,iz),dy(ki,kj,iz),dz(ki,kj,iz)
          enddo
        enddo
      endif
      close(70+iz)
    enddo z_dip_r

    spin:do ns=1,nsmax
      if(maxval(is(:,ns,:)).eq.0) cycle 
      
      ! open files
      inquire(file=workdir(1:lwork)//'/dx_'//spin_str(ns)//com(i,1)(1:l1)//com(i,2)(1:l2),EXIST=exist)
      if(exist)then
        open(10+i,file=workdir(1:lwork)//'/dx_'//spin_str(ns)//com(i,1)(1:l1)//com(i,2)(1:l2),status='old',position='append')
      else
        open(10+i,file=workdir(1:lwork)//'/dx_'//spin_str(ns)//com(i,1)(1:l1)//com(i,2)(1:l2),status='new')
      endif
      inquire(file=workdir(1:lwork)//'/dy_'//spin_str(ns)//com(i,1)(1:l1)//com(i,2)(1:l2),EXIST=exist)
      if(exist)then
        open(30+i,file=workdir(1:lwork)//'/dy_'//spin_str(ns)//com(i,1)(1:l1)//com(i,2)(1:l2),status='old',position='append')
      else
        open(30+i,file=workdir(1:lwork)//'/dy_'//spin_str(ns)//com(i,1)(1:l1)//com(i,2)(1:l2),status='new')
      endif
      inquire(file=workdir(1:lwork)//'/dz_'//spin_str(ns)//com(i,1)(1:l1)//com(i,2)(1:l2),EXIST=exist)
      if(exist)then
        open(50+i,file=workdir(1:lwork)//'/dz_'//spin_str(ns)//com(i,1)(1:l1)//com(i,2)(1:l2),status='old',position='append')
      else
        open(50+i,file=workdir(1:lwork)//'/dz_'//spin_str(ns)//com(i,1)(1:l1)//com(i,2)(1:l2),status='new')
      endif

      write(10+i,*)'# *****************************************'
      write(10+i,*)'# Spin symmetry : '//spin_str(ns)
      write(10+i,*)'# *****************************************'
      write(30+i,*)'# *****************************************'      
      write(30+i,*)'# Spin symmetry : '//spin_str(ns)
      write(30+i,*)'# *****************************************'
      write(50+i,*)'# *****************************************'      
      write(50+i,*)'# Spin symmetry : '//spin_str(ns)
      write(50+i,*)'# *****************************************'      
      do ki=1,is(1,ns,is_com(i,1))
        write(10+i,fmt_com)'#',(ki,kj,kj=1,is(1,ns,is_com(i,2)))
        write(10+i,fmt_com)'#',(ie_s(ki,1,ns,is_com(i,1)),ie_s(kj,1,ns,is_com(i,2)),kj=1,is(1,ns,is_com(i,2)))
        write(30+i,fmt_com)'#',(ki,kj,kj=1,is(1,ns,is_com(i,2)))
        write(30+i,fmt_com)'#',(ie_s(ki,1,ns,is_com(i,1)),ie_s(kj,1,ns,is_com(i,2)),kj=1,is(1,ns,is_com(i,2)))
        write(50+i,fmt_com)'#',(ki,kj,kj=1,is(1,ns,is_com(i,2)))
        write(50+i,fmt_com)'#',(ie_s(ki,1,ns,is_com(i,1)),ie_s(kj,1,ns,is_com(i,2)),kj=1,is(1,ns,is_com(i,2)))
        
        allocate(dip_value(3,maxval(is(:,ns,is_com(i,2))),3),phase(maxval(is(:,ns,is_com(i,2))),3))        
        dip_value=0d0
        phase=1d0
        
        do iz1=1,jz_s(is_com(i,1))                             !z valides pour la 1ere sym de la transition
          do iz2=1,jz_s(is_com(i,2))                 !recherche du z identique valides pour la 2eme sym de la transition
            if(iz_f(iz2,is_com(i,2))==iz_f(iz1,is_com(i,1)))then
              if((ie_s(ki,iz1,ns,is_com(i,1)).lt.1))cycle
              
            ! attempt to obtain the sign of the transition dipole moment :
            ! the sign change value when the change of derivative between two adjacent R value is greater than phase_variation, only work for dimer  
              do kj=1,is(iz2,ns,is_com(i,2))
                dip_value(3,kj,1)=abs(dx(ie_s(ki,iz1,ns,is_com(i,1)),ie_s(kj,iz2,ns,is_com(i,2)),iz_f(iz1,is_com(i,1))))
                dip_value(3,kj,2)=abs(dy(ie_s(ki,iz1,ns,is_com(i,1)),ie_s(kj,iz2,ns,is_com(i,2)),iz_f(iz1,is_com(i,1))))
                dip_value(3,kj,3)=abs(dz(ie_s(ki,iz1,ns,is_com(i,1)),ie_s(kj,iz2,ns,is_com(i,2)),iz_f(iz1,is_com(i,1))))              
                
                if((iz1.lt.3).or.((com(i,1)(1:l1).eq.com(i,2)(1:l2)).and.(ki.eq.kj)).or.(.not.phase_sort)) then
                  phase(kj,1)=sign(1d0,dx(ie_s(ki,iz1,ns,is_com(i,1)),ie_s(kj,iz2,ns,is_com(i,2)),iz_f(iz1,is_com(i,1))))
                  phase(kj,2)=sign(1d0,dy(ie_s(ki,iz1,ns,is_com(i,1)),ie_s(kj,iz2,ns,is_com(i,2)),iz_f(iz1,is_com(i,1))))
                  phase(kj,3)=sign(1d0,dz(ie_s(ki,iz1,ns,is_com(i,1)),ie_s(kj,iz2,ns,is_com(i,2)),iz_f(iz1,is_com(i,1))))    
                else
                  do axis=1,3
                    phase_test=((dip_value(3,kj,axis)-dip_value(2,kj,axis))/(dip_value(3,kj,axis)+1d-20))-((dip_value(2,kj,axis)-dip_value(1,kj,axis))/(dip_value(2,kj,axis)+1d-20))
                    if(phase_test.gt.phase_variation) phase(kj,axis)=-phase(kj,axis)
                  enddo 
                endif  
                
                dip_value(1,:,:)=dip_value(2,:,:)
                dip_value(2,:,:)=dip_value(3,:,:)                
              enddo
              
              
              write(10+i,fmt_dip)(coord(je,iz_f(iz1,is_com(i,1))),je=1,ncoord),(phase(kj,1)*abs(dx(ie_s(ki,iz1,ns,is_com(i,1)),ie_s(kj,iz2,ns,is_com(i,2)),iz_f(iz1,is_com(i,1)))),&
             &                   kj=1,is(iz2,ns,is_com(i,2))),(0.,je=1,maxval(is(:,ns,is_com(i,2)))-is(iz2,ns,is_com(i,2)))
              write(30+i,fmt_dip)(coord(je,iz_f(iz1,is_com(i,1))),je=1,ncoord),(phase(kj,2)*abs(dy(ie_s(ki,iz1,ns,is_com(i,1)),ie_s(kj,iz2,ns,is_com(i,2)),iz_f(iz1,is_com(i,1)))),&
             &                   kj=1,is(iz2,ns,is_com(i,2))),(0.,je=1,maxval(is(:,ns,is_com(i,2)))-is(iz2,ns,is_com(i,2)))
              write(50+i,fmt_dip)(coord(je,iz_f(iz1,is_com(i,1))),je=1,ncoord),(phase(kj,3)*abs(dz(ie_s(ki,iz1,ns,is_com(i,1)),ie_s(kj,iz2,ns,is_com(i,2)),iz_f(iz1,is_com(i,1)))),&
             &                   kj=1,is(iz2,ns,is_com(i,2))),(0.,je=1,maxval(is(:,ns,is_com(i,2)))-is(iz2,ns,is_com(i,2)))
            endif
          enddo
        enddo
        
        deallocate(dip_value,phase)
      enddo
      
      close(10+i)
      close(30+i)
      close(50+i)
      
     ! specific sorting for permanent dipole moment 
      if(com(i,1)(1:l1).eq.com(i,2)(1:l2)) then
      
        select case(ncoord)
        
        case (1)		! dimer : permanent dipole along z axis
          inquire(file=workdir(1:lwork)//'/dperm_'//spin_str(ns)//com(i,1)(1:l1),EXIST=exist)
          if(exist)then
            open(10+i,file=workdir(1:lwork)//'/dperm_'//spin_str(ns)//com(i,1)(1:l1),status='old',position='append')
          else
            open(10+i,file=workdir(1:lwork)//'/dperm_'//spin_str(ns)//com(i,1)(1:l1),status='new')
          endif   
        
          do iz=1,jz_s(is_com(i,1)) 
!            if((ie_s(ki,iz1,ns,is_com(i,1)).lt.1))cycle
            write(10+i,fmt_dip)(coord(je,iz_f(iz,is_com(i,1))),je=1,ncoord),(dz(ie_s(ki,iz,ns,is_com(i,1)),ie_s(ki,iz,ns,is_com(i,1)),iz_f(iz,is_com(i,1))),&
           &                   ki=1,is(iz,ns,is_com(i,1))),(0.,je=1,maxval(is(:,ns,is_com(i,1)))-is(iz,ns,is_com(i,1)))        
          enddo      
          close(10+i)
          
        case (3) 		! trimer : permanent dipole along x and y axis
          inquire(file=workdir(1:lwork)//'/dpermx_'//spin_str(ns)//com(i,1)(1:l1),EXIST=exist)
          if(exist)then
            open(10+i,file=workdir(1:lwork)//'/dpermx_'//spin_str(ns)//com(i,1)(1:l1),status='old',position='append')
          else
            open(10+i,file=workdir(1:lwork)//'/dpermx_'//spin_str(ns)//com(i,1)(1:l1),status='new')
          endif   
        
          do iz=1,jz_s(is_com(i,1)) 
!            if((ie_s(ki,iz1,ns,is_com(i,1)).lt.1))cycle
            write(10+i,fmt_dip)(coord(je,iz_f(iz,is_com(i,1))),je=1,ncoord),(dx(ie_s(ki,iz,ns,is_com(i,1)),ie_s(ki,iz,ns,is_com(i,1)),iz_f(iz,is_com(i,1))),&
           &                   ki=1,is(iz,ns,is_com(i,1))),(0.,je=1,maxval(is(:,ns,is_com(i,1)))-is(iz,ns,is_com(i,1)))        
          enddo      
          close(10+i) 
          
          inquire(file=workdir(1:lwork)//'/dpermy_'//spin_str(ns)//com(i,1)(1:l1),EXIST=exist)
          if(exist)then
            open(10+i,file=workdir(1:lwork)//'/dpermy_'//spin_str(ns)//com(i,1)(1:l1),status='old',position='append')
          else
            open(10+i,file=workdir(1:lwork)//'/dpermy_'//spin_str(ns)//com(i,1)(1:l1),status='new')
          endif   
        
          do iz=1,jz_s(is_com(i,1)) 
!            if((ie_s(ki,iz1,ns,is_com(i,1)).lt.1))cycle
            write(10+i,fmt_dip)(coord(je,iz_f(iz,is_com(i,1))),je=1,ncoord),(dy(ie_s(ki,iz,ns,is_com(i,1)),ie_s(ki,iz,ns,is_com(i,1)),iz_f(iz,is_com(i,1))),&
           &                   ki=1,is(iz,ns,is_com(i,1))),(0.,je=1,maxval(is(:,ns,is_com(i,1)))-is(iz,ns,is_com(i,1)))        
          enddo      
          close(10+i)          
        end select
        
      endif
      
    enddo spin 
  enddo com_dip  !index i

  command='gzip -f '//workdir(1:lwork)//'/*ciro*'
  call system(command)
 
  command='mkdir -p '//workdir(1:lwork)//'/Results/'
  call system(command)
  command='mv -f '//workdir(1:lwork)//'/*.gz '//workdir(1:lwork)//'/Results/'
  call system(command)
  
  enddo ! JD 04/06 End of loop over e-field values !
  
  ! finite field computation
  if(l_ffield.and.l_cip) then
    command='mkdir -p '//workdirroot(1:lworkroot)//'/ffield'
    call system(command)
    
    do i=1,nsym
      do ns=1,nsmax
        if(is_field(ns,i).eq.0) cycle 
        call finitefield(field_dir,spin_str(ns)//sym(i),l_accept(:,:,ns,i),efield_nval,efield_vals,is_field(ns,i))   
      enddo
    enddo
  endif
  
  
  deallocate(isymat,spin_str)
  if(allocated(l_accept)) deallocate(l_accept)
  
  !RV 15/12 write total time
  call system_clock(count=time_end, count_rate=ir)
  total_time=real(time_end - time_begin,kind=8)/real(ir,kind=8)
  write(*,*) 'Computation time : ',total_time,' seconds '

end program auto


