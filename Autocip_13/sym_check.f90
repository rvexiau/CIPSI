subroutine sym_check(group,nsym,sym,nusym,mult_tab)
! RV 01/2016 : also give multiplication table of the group
!              the tables have to be the same than the one define in BLOCKDATA1.f of ijkl
!  CONVENTION la notation des symetries:
!       - groupes DINFH,CINFV,CNV,CNH,CS,C1 seuls admis ici
!  Les symetries sont identifiees par leur numero, selon la table
!  de correspondance de HONDO:
!  DINFH:
!  1: spg(sigma_plus_g),d2g(delta_(x^2-y^2)g);  3:pxg(pi_x_g);  5:pyg(pi_y_g);
!  2: spu(sigma_plus_u),d2u(delta_(x^2-y^2)u);  4:pxu(pi_x_u);  6:pyu(pi_y_u);
!  7=1 ; 9:smg(sigma_moins_g),dcg(delta_(xy)g);
!  8=2 ;10:smu(sigma_moins_u),dcg(delta_(xy)u);
!  CINFV
!  1: sp(sigma_plus),d2(delta_(x^2-y^2));
!  2: px(pi_x)
!  3: py(pi_y)
!  4=1
!  5: sm(sigma_moins),dc(delta_(xy))
!
!  JD Juin 06: added support for C2V (group=CNV and naxis=2) and CS to support calculations of polarizability
!  With orthorgonal electric field :
!     C2V: 1, 2, 3, 4  
!      CS: a1p, a2p    
!  RV Janvier 16 : added support for C2H (group=CNH and naxis=2) to support calculations of tetramer
!  		symmetries labeled by
!  With orthorgonal electric field :
!   - 1 number for the symmetry:
!     C2H : ag, bg, au, bu  
    

  implicit none
  integer, parameter :: dp=kind(1.d0)
  integer :: nsym,nusym(nsym),i,j
  character(len=3) :: sym_dih(10),sym_dih_c(10)
  character(len=2) :: sym_civ(5),sym_civ_c(5)
  character(len=2) :: sym_c2v(4)
  character(len=2) :: sym_c2h(4)  
  character(len=3) :: sym_d2h(8)    
  character(len=3) :: sym_cs(2)  
  character(len=1) :: sym_c1(1)    
  character(len=5) :: sym_grp(7),sym_grp_c(7)
  character(len=*) :: sym(nsym),group
  integer :: itsdh(14,14),itscv(7,7),its8(8,8),its4(4,4),its2(2,2),its1(1,1)
  integer,dimension(14,14) :: mult_tab   

  data sym_dih/'spg','spu','pxg','pxu','pyg','pyu','d2g','d2u','smg','smu'/
  data sym_dih_c/'SPG','SPU','PXG','PXU','PYG','PYU','D2G','D2U','SMG','SMU'/
  data sym_civ/'sp','px','py','d2','sm'/
  data sym_civ_c/'SP','PX','PY','D2','SM'/
  data sym_c2v/'a1','a2','b1','b2'/
  data sym_c2h/'ag','bu','au','bg'/
  data sym_d2h/'a1g','a1u','b3u','b3g','b1g','b1u','b2u','b2g'/  
  data sym_cs/'a1p','a2p'/
  data sym_c1/'a'/  
  data sym_grp/'dinfh','cinfv','cnv','cnh','cs','c1','dnh'/
  data sym_grp_c/'DINFH','CINFV','CNV','CNH','CS','C1','DNH'/
  data itsdh / 1, 2, 3, 4, 5, 6, 1, 2, 9,10, 3, 4, 5, 6,&
               2, 1, 4, 3, 6, 5, 2, 1,10, 9, 4, 3, 6, 5,&
               3, 4, 1, 2, 9,10, 3, 4, 5, 6, 1, 2, 9,10,&
               4, 3, 2, 1,10, 9, 4, 3, 6, 5, 2, 1,10, 9,&
               5, 6, 9,10, 1, 2, 5, 6, 3, 4, 9,10, 1, 2,&
               6, 5,10, 9, 2, 1, 6, 5, 4, 3,10, 9, 2, 1,&
               1, 2, 3, 4, 5, 6, 1, 2, 9,10, 3, 4, 5, 6,&
               2, 1, 4, 3, 6, 5, 2, 1,10, 9, 4, 3, 6, 5,&
               9,10, 5, 6, 3, 4, 9,10, 1, 2, 5, 6, 3, 4,&
              10, 9, 6, 5, 4, 3,10, 9, 2, 1, 6, 5, 4, 3,&
               3, 4, 1, 2, 9,10, 3, 4, 5, 6, 1, 2, 9,10,&
               4, 3, 2, 1,10, 9, 4, 3, 6, 5, 2, 1,10, 9,&
               5, 6, 9,10, 1, 2, 5, 6, 3, 4, 9,10, 1, 2,&
               6, 5,10, 9, 2, 1, 6, 5, 4, 3,10, 9, 2, 1/
  data itscv/ 1,2,3,1,5,2,3,&
              2,1,5,2,3,1,5,&
              3,5,1,3,2,5,1,&
              1,2,3,1,5,2,3,&
              5,3,2,5,1,3,2,&
              2,1,5,2,3,1,5,&
              3,5,1,3,2,5,1/
  data its2 / 1,2,&
              2,1 /
  data its4 / 1,2,3,4,&
              2,1,4,3,&
              3,4,1,2,&
              4,3,2,1       /
  data its8 / 1,2,3,4,5,6,7,8,&
              2,1,4,3,6,5,8,7,&
              3,4,1,2,7,8,5,6,&
              4,3,2,1,8,7,6,5,&       
              5,6,7,8,1,2,3,4,&
              6,5,8,7,2,1,4,3,&
              7,8,5,6,3,4,1,2,&
              8,7,6,5,4,3,2,1       /
  data its1 /1/

  mult_tab=0
  nusym=0

  if((group==sym_grp(1)).or.(group==sym_grp_c(1)))then
    do i=1,nsym
      do j=1,10
        if((sym(i)(1:3)==sym_dih(j)).or.(sym(i)(1:3)==sym_dih_c(j)))nusym(i)=j
      enddo
      if(nusym(i)==0)then
        write( 6, * ) 'erreur symetrie ',i
        stop
      endif
    enddo
    do i=1,14
      do j=1,14
        mult_tab(j,i)=itsdh(j,i)
      enddo
    enddo  
  elseif((group==sym_grp(2)).or.(group==sym_grp_c(2)))then
    do i=1,nsym
      do j=1,5
        if((sym(i)(1:2)==sym_civ(j)).or.(sym(i)(1:2)==sym_civ_c(j)))nusym(i)=j
      enddo
      if(nusym(i)==0)then
        write( 6, * ) 'erreur symetrie ',i
        stop
      endif
    enddo
    do i=1,7
      do j=1,7
        mult_tab(j,i)=itscv(j,i)
      enddo
    enddo  
  elseif((group==sym_grp(3)).or.(group==sym_grp_c(3)))then
    do i=1,nsym
      do j=1,4
        if(sym(i)(1:2)==sym_c2v(j))nusym(i)=j
      enddo
      if(nusym(i)==0)then
        write( 6, * ) 'erreur symetrie ',i
        stop
      endif
      sym(i)=sym(i)(1:2)//'p'
    enddo
    do i=1,4
      do j=1,4
        mult_tab(j,i)=its4(j,i)
      enddo
    enddo    
  elseif((group==sym_grp(4)).or.(group==sym_grp_c(4)))then
    do i=1,nsym
      do j=1,4
        if(sym(i)(1:2)==sym_c2h(j))nusym(i)=j
      enddo
      if(nusym(i)==0)then
        write( 6, * ) 'erreur symetrie ',i
        stop
      endif
    enddo
    do i=1,4
      do j=1,4
        mult_tab(j,i)=its4(j,i)
      enddo
    enddo     
  elseif((group==sym_grp(5)).or.(group==sym_grp_c(5)))then
    do i=1,nsym
      do j=1,2
        if(sym(i)(1:3)==sym_cs(j))nusym(i)=j
      enddo
      if(nusym(i)==0)then
        write( 6, * ) 'erreur symetrie ',i
        stop
      endif
    enddo
    do i=1,2
      do j=1,2
        mult_tab(j,i)=its2(j,i)
      enddo
    enddo      
  elseif((group==sym_grp(6)).or.(group==sym_grp_c(6)))then
    do i=1,nsym
      if(sym(i)(1:1)==sym_c1(1))nusym(i)=1
      if(nusym(i)==0)then
        write( 6, * ) 'erreur symetrie ',i
        stop
      endif
    enddo   
    mult_tab(1,1) = its1(1,1)
  elseif((group==sym_grp(7)).or.(group==sym_grp_c(7)))then
    do i=1,nsym
      do j=1,8
        if(sym(i)(1:3)==sym_d2h(j))nusym(i)=j
      enddo
      if(nusym(i)==0)then
        write( 6, * ) 'erreur symetrie ',i
        stop
      endif
    enddo
    do i=1,8
      do j=1,8
        mult_tab(j,i)=its8(j,i)
      enddo
    enddo         
  else
    write( 6, * ) 'erreur groupe de symetrie'
    stop
  endif

  return
end