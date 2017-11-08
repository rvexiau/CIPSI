!    Set of functions defining the fundamental constants to be used 
!    in any program
!                  LAST MODIFICATION: 23.10.2000, C. Dion
!
!    Source: NIST Web page http://physics.nist.gov

module convmod
  implicit none
  integer, parameter :: dp=kind(1.d0)
  real(dp), parameter :: amu2kg=1.66053873d-27, & 
       au2kg=9.10938188d-31, &
       au2A=0.5291772083d0, &
       au2s=2.418884326500d-17, &
       au2m=0.5291772083d-10, &
       au2C=1.602176462d-19, &
       au2C_m=8.47835267d-30, & ! dipole moment (C m)
       au2C2m2J=1.648777251d-41, & ! polarizability (C2 m2 J-1)
       au2J=4.35974381d-18, & ! Energies
       au2Hz=6.579683920735d15, &
       au2cm=2.194746313710d5, &
       au2K=3.1577465d5, &
       au2eV=27.2113834d0, &
       au2V=27.2113834d0, &
       J2cm=5.03411762d22, &
       J2eV=6.24150974d18, &
       J2Hz=1.50919050d33, &
       J2K=7.242964d22, &
       eV2cm=8.06554477d3, &
       eV2Hz=2.417989491d14, &
       eV2K=1.1604506d4
end module convmod



function conv(u_type,u_in,u_out)
  !
  ! This function returns the conversion factor needed to convert 
  ! from unit 'u_in' to unit 'u_out'
  !
  ! Usage:
  ! y=x*conv('u_type','u_in','u_out')
  !
  ! where 'u_in' is the unit of x and 'u_out' is the unit of y
  ! and 'u_type' the type of of unit
  !
  ! u_type = E (energy)
  !        = L (length)
  !        = A (area)
  !        = V (volume)
  !        = T (time)
  !        = M (mass)
  !        = EP (electric potential)
  ! 
  ! All conversions proceed from 'u_in' units to SI units, then 
  ! from SI units to 'u_out' units.  This may cause a small discrepancy
  ! between this result and the one obtained using the approved value,
  ! but should in any case be less than the uncertainty on the 
  ! approved values themselves
  !
  use convmod
  implicit none
  character(len=*), intent(in) :: u_type,u_in,u_out
  real(dp) :: conv
  
  type_of_u: select case (trim(u_type))
     
  case ('E','e')
     ! Energy

     select case (trim(u_in))
     case ('J','j','SI')
        conv=1.d0
     case ('au','a.u.','hartree','Hartree','Eh')
        conv=1.d0*au2J
     case ('eV','ev','EV')
        conv=1.d0/J2eV
     case ('cm-1','1/cm','cm')
        conv=1.d0/J2cm
     case ('Hz','HZ','hz')
        conv=1.d0/J2Hz
     case default
        write(*,*) 'Unit ',trim(u_in),' for type energy unknown in ',&
             & 'function conv'
        stop
     end select

     select case (trim(u_out))
     case ('J','j','SI')
        ! do nothing
     case ('au','a.u.','hartree','Hartree','Eh')
        conv=conv/au2J
     case ('eV','ev','EV')
        conv=conv*J2eV
     case ('cm-1','1/cm','cm')
        conv=conv*J2cm
     case ('Hz','HZ','hz')
        conv=conv*J2Hz
     case default
        write(*,*) 'Unit ',trim(u_out),' for type energy unknown in ',&
             & 'function conv'
        stop
     end select
     
  case ('L','l')
     ! Length
     select case (trim(u_in))
     case ('m','M','SI')
        conv=1.d0
     case ('cm','CM')
        conv=1.d-2
     case ('mm','MM')
        conv=1.d-3
     case ('km','KM')
        conv=1.d3
     case ('micron','Micron')
        conv=1.d-6
     case ('A','angstrom','Angstrom')
        conv=1.d-10
     case ('au','a.u.','bohr','Bohr','a0')
        conv=1.d0*au2m
     case default
        write(*,*) 'Unit ',trim(u_in),' for type length unknown in ',&
             & 'function conv'
        stop
     end select

     select case (trim(u_out))
     case ('m','M','SI')
        ! do nothing
     case ('cm','CM')
        conv=conv*1.d2
     case ('mm','MM')
        conv=conv*1.d3
     case ('km','KM')
        conv=conv*1.d-3
     case ('A','angstrom','Angstrom')
        conv=conv*1.d10
     case ('micron','Micron')
        conv=conv*1.d6
     case ('au','a.u.','bohr','Bohr','a0')
        conv=conv/au2m
     case default
        write(*,*) 'Unit ',trim(u_out),' for type length unknown in ',&
             & 'function conv'
        stop
     end select


  case ('A','a')
     ! Area
     write(*,*) 'Unit conversion for type area not yet implemented ', &
          & 'in function conv'
     stop

  case ('V','v')
     ! Volume
     write(*,*) 'Unit conversion for type volume not yet implemented ', &
          & 'in function conv'
     stop

  case ('T','t')
     ! Time

     select case (trim(u_in))
     case ('s','S','SI')
        conv=1.d0
     case ('ms','MS')
        conv=1.d-3
     case ('ns','NS')
        conv=1.d-9
     case ('ps','PS')
        conv=1.d-12
     case ('fs','FS')
        conv=1.d-15
     case ('au','a.u.')
        conv=au2s
     case default
        write(*,*) 'Unit ',trim(u_in),' for type time unknown in ',&
             & 'function conv'
        stop
     end select

     select case (trim(u_out))
     case ('s','S','SI')
        ! do nothing
     case ('ms','MS')
        conv=conv*1.d3
     case ('ns','NS')
        conv=conv*1.d9
     case ('ps','PS')
        conv=conv*1.d12
     case ('fs','FS')
        conv=conv*1.d15
     case ('au','a.u.')
        conv=conv/au2s
     case default
        write(*,*) 'Unit ',trim(u_in),' for type time unknown in ',&
             & 'function conv'
        stop
     end select

  case ('M','m')
     ! Mass

     select case (trim(u_in))
     case ('kg','KG','SI')
        conv=1.d0
     case ('g','G')
        conv=1.d-3
     case ('amu','AMU','a.m.u.','A.M.U.')
        conv=amu2kg
     case ('au','a.u.')
        conv=au2kg
     case default
        write(*,*) 'Unit ',trim(u_in),' for type mass unknown in ',&
             & 'function conv'
        stop
     end select

     select case (trim(u_out))
     case ('kg','KG','SI')
        ! do nothing
     case ('g','G')
        conv=conv*1.d3
     case ('amu','AMU','a.m.u.','A.M.U.')
        conv=conv/amu2kg
     case ('au','a.u.')
        conv=conv/au2kg
     case default
        write(*,*) 'Unit ',trim(u_in),' for type mass unknown in ',&
             & 'function conv'
        stop
     end select

  case ('EP','ep')
     ! Electrical potential

     select case (trim(u_in))
     case ('V','v','SI')
        conv=1.d0
     case ('au','a.u.')
        conv=au2V
     case default
        write(*,*) 'Unit ',trim(u_in),' for type electric potential ',&
             & 'unknown in function conv'
        stop
     end select

     select case (trim(u_out))
     case ('V','v','SI')
        ! do nothing
     case ('au','a.u.')
        conv=conv/au2V
     case default
        write(*,*) 'Unit ',trim(u_in),' for type electric potential ',&
             & 'unknown in function conv'
        stop
     end select


  case default
     write(*,*) 'Unit type ',u_type,' unknown in function conv'
     stop
     
  end select type_of_u

  return

end function conv



recursive function fundcst(name,unit) result(cst)
  use convmod
  implicit none
  character(len=*), intent(in) :: name,unit
  real(dp) :: cst
  real(dp), external :: conv
  
  ! All base units are SI

  select case (trim(name))

  case ('pi','Pi','PI')
     cst=3.14159265358979323846d0

  case ('c','C') ! Speed of light in vacuum
     cst=299792458.d0 ! m/s

     select case (trim(unit))
     case ('ms-1','m s-1','m/s','SI')
     case ('cms-1','cm s-1','cm/s')
        cst=cst*1.d2
     case ('au','a.u.')
        cst=cst/au2m*au2s
     case default
        write(*,*) 'Unit ',trim(unit),' not defined for c'
        stop
     end select

  case ('h','Planck','planck') ! Planck's constant
     cst=6.62606876d-34
     
     select case (trim(unit))
     case ('SI')
     case ('eVs','eV s')
        cst=cst*J2eV
     case ('au','a.u.')
        cst=cst/au2J/au2s
     case default
        write(*,*) 'Unit ',trim(unit),' not defined for h'
        stop
     end select

  case ('hbar') ! Planck's constant / 2 pi
     cst=1.054571596d-34 
     
     select case (trim(unit))
     case ('Js','J s','SI')
     case ('eVs','eV s')
        cst=cst*J2eV
     case ('au','a.u.')
        cst = 1.d0
     case default
        write(*,*) 'Unit ',trim(unit),' not defined for hbar'
        stop
     end select

  case ('Na','NA','Avogadro','avogadro') ! Avogadro's number
     cst=6.02214199d23

  case ('k','kB','Boltzmann','boltzmann')
     cst=1.3806503d-23 ! J/K

     select case (trim(unit))
     case ('JK-1','J K-1','J/K','SI')
     case ('eVK-1','eV K-1','eV/K')
        cst=cst*J2eV
     case ('au','a.u.')
        cst=cst/au2J
     case default
        write(*,*) 'Unit ',trim(unit),' not defined for kB'
        stop
     end select
     
  case ('R') ! molar gas constant
     cst=8.314472d0

     select case (trim(unit))
     case ('Jmol-1K-1','JK-1mol-1','J mol-1 K-1','J K-1 mol-1','J/K/mol','SI')
     case ('eVmol-1K-1','eVK-1mol-1','eV mol-1 K-1','eV K-1 mol-1','eV/K/mol')
        cst=cst*J2eV
     case ('au','a.u.')
        cst=cst/au2J
     case default
        write(*,*) 'Unit ',trim(unit),' not defined for molar gas constant'
        stop
     end select
     
  case ('me') ! electron mass
     cst=9.10938188d-31*conv('M','kg',unit) 

  case ('e') ! electron charge
     cst=1.602176462d-19 ! C

     select case (trim(unit))
     case ('C','SI')
     case ('au','a.u.')
        cst=1.d0
     case default
        write(*,*) 'Unit ',trim(unit),' not defined for electron charge'
        stop
     end select

  case ('mp') ! proton mass
     cst=1.67262158d-27*conv('M','kg',unit)

  case ('F','Faraday','faraday') ! Faraday constant
     cst=96485.3415d0
     select case (trim(unit))
     case ('C','Cmol-1','C mol-1','C/mol','SI')
     case ('au','a.u.')
        cst=cst/au2C
     case default
        write(*,*) 'Unit ',trim(unit),' not defined for Faraday constant'
        stop
     end select

  case ('Rinf','Rydberg','rydberg','Ryd') ! Rydberg constant
     cst=10973731.568549d0 ! m-1

     select case (trim(unit))
     case ('m-1','1/m','SI')
     case ('cm-1','1/cm')
        cst=cst*1.d-2
     case ('J')
        cst=cst*1.d-2/J2cm
     case ('eV','ev','EV')
        cst=cst*1.d-2/eV2cm
     case default
        write(*,*) 'Unit ',trim(unit),' not defined for Rydberg constant'
        stop
     end select
     
  case ('mu0','MU0','magnetic') ! magnetic constant
     cst=4.d-7*fundcst('pi','') ! N A^-2
     
     select case (trim(unit))
     case ('NA-2','N/A2','SI')
     case default
        write(*,*) 'Unit ',trim(unit),' not defined for magnetic constant'
        stop
     end select

  case ('e0','E0','epsilon','electric') ! electric constant
     cst=1.d0/(fundcst('mu0','SI')*fundcst('c','SI')**2) ! F/m

     select case (trim(unit))
     case ('F/m','Fm-1','SI')
     case default
        write(*,*) 'Unit ',trim(unit),' not defined for electric constant'
        stop
     end select

  case ('alpha','fs','fine-structure','FS') ! fine-structure constant
     cst=7.297352533d-3

  case default
     write(*,*) 'Fundamental constant ',trim(name),' not defined'
     stop
  end select

  return

end function fundcst



function rydberg()
  ! Rydberg constant in cm**(-1)
  !
  ! Obsolete!  To be replaced by function fundcst
  !
  implicit none
  integer, parameter :: dp=kind(1.d0)
  real(dp) :: rydberg,fundcst
  rydberg = fundcst('rydberg','cm-1')
  return
end function rydberg



function mass(atom,unit)
  !
  ! Retruns the atomic mass of 'atom'
  ! 
  ! atom = 'nAa' where n = mass number
  !                    Aa = symbol
  !
  ! unit = any mass unit understood by function conv
  !
  ! If no mass number is given, atomic weigth is retruned
  !
  ! Source:  http://physics.nist.gov/cgi-bin/Compositions/index.html
  ! (C. Dion, 28/07/2000)
  !
  implicit none
  integer, parameter :: dp=kind(1.d0)
  character(len=*), intent(in) :: atom,unit
  real(dp) :: mass,conv

  select case (trim(atom))
  case('NO')
     mass=30
  case ('1H')
     mass=1.0078250321d0
  case ('2H','D')
     mass=2.0141017780d0
  case ('3H','T')
     mass=3.0160492675d0
  case ('H')
     mass=1.00794d0
  case ('3He')
     mass=3.0160293097d0
  case ('4He')
     mass=4.0026032497d0
  case ('He')
     mass=4.002602d0
  case ('6Li')
     mass=6.0151223d0
  case ('7Li')
     mass=7.0160040d0
  case ('Li')
     mass=6.941d0
  case ('9Be')
     mass=9.0121821d0
  case ('Be')
     mass=9.012182d0
  case ('10B')
     mass=10.0129370d0
  case ('11B')
     mass=11.0093055d0
  case ('B')
     mass=10.811d0
  case ('12C')
     mass=12.d0
  case ('13C')
     mass=13.0033548378d0
  case ('C')
     mass=12.0107d0
  case ('14N')
     mass=14.0030740052d0
  case ('15N')
     mass=15.0001088984d0
  case ('N')
     mass=14.00674d0
  case ('16O')
     mass=15.9949146221d0
  case ('17O')
     mass=16.99913150d0
  case ('18O')
     mass=17.9991604d0
  case ('O')
     mass=15.9994d0
  case ('19F','F')
     mass=18.9984032d0
  case ('20Ne')
     mass=19.9924401759d0
  case ('21Ne')
     mass=20.99384674d0
  case ('22Ne')
     mass=21.99138551d0
  case ('Ne')
     mass=20.1797d0
  case ('23Na')
     mass=22.98976967d0
  case ('Na')
     mass=22.989770d0
  case ('24Mg')
     mass=23.98504190d0
  case ('25Mg')
     mass=24.98583702d0
  case ('26Mg')
     mass=25.98259304d0
  case ('Mg')
     mass=24.3050d0
  case ('27Al')
     mass=26.98153844d0
  case ('Al')
     mass=26.981538d0
  case ('28Si')
     mass=27.9769265327d0
  case ('29Si')
     mass=28.97649472d0
  case ('30Si')
     mass=29.97377022d0
  case ('Si')
     mass=28.0855d0
  case ('31P')
     mass=30.97376151d0
  case ('P')
     mass=30.973761d0
  case ('32S')
     mass=31.97207069d0
  case ('33S')
     mass=32.97145850d0
  case ('34S')
     mass=33.96786683d0
  case ('36S')
     mass=35.96708088d0
  case ('S')
     mass=32.066d0
  case ('35Cl')
     mass=34.96885271d0
  case ('37Cl')
     mass=36.96590260d0
  case ('Cl')
     mass=35.4527d0
  case ('36Ar')
     mass=35.96754628d0
  case ('38Ar')
     mass=37.9627322d0
  case ('40Ar')
     mass=39.962383123d0
  case ('Ar')
     mass=39.948d0
  case ('39K')
     mass=38.9637069d0
  case ('40K')
     mass=39.96399867d0
  case ('41K')
     mass=40.96182597d0
  case ('K')
     mass=39.0983d0
  case ('40Ca')
     mass=39.9625912d0
  case ('42Ca')
     mass=41.9586183d0
  case ('43Ca')
     mass=42.9587668d0
  case ('44Ca')
     mass=43.9554811d0
  case ('46Ca')
     mass=45.9536928d0
  case ('48Ca')
     mass=47.952534d0
  case ('Ca')
     mass=40.078d0
  case ('45Sc')
     mass=44.9559102d0
  case ('Sc')
     mass=44.955910d0
  case ('46Ti')
     mass=45.9526295d0
  case ('47Ti')
     mass=46.9517638d0
  case ('48Ti')
     mass=47.9479471d0
  case ('49Ti')
     mass=48.9478708d0
  case ('50Ti')
     mass=49.9447921d0
  case ('Ti')
     mass=47.867d0
  case ('50V')
     mass=49.9471628d0
  case ('51V')
     mass=50.9439637d0
  case ('V')
     mass=50.9415d0
  case ('50Cr')
     mass=49.9460496d0
  case ('52Cr')
     mass=51.9405119d0
  case ('53Cr')
     mass=52.9406538d0
  case ('54Cr')
     mass=53.9388849d0
  case ('Cr')
     mass=51.9961d0
  case ('55Mn')
     mass=54.9380496d0
  case ('Mn')
     mass=54.938049d0
  case ('54Fe')
     mass=53.9396148d0
  case ('56Fe')
     mass=55.9349421d0
  case ('57Fe')
     mass=56.9353987d0
  case ('58Fe')
     mass=57.9332805d0
  case ('Fe')
     mass=55.845d0
  case ('59Co')
     mass=58.9332002d0
  case ('Co')
     mass=58.933200d0
!_____________________________________________________
!Ni  58   57.9353479(15)     68.0769(89)   58.6934(2) 
!    60   59.9307906(15)     26.2231(77)              
!    61   60.9310604(15)     1.1399(6)                
!    62   61.9283488(15)     3.6345(17)               
!    64   63.9279696(16)     0.9256(9)                
!_____________________________________________________
!Cu  63   62.9296011(15)     69.17(3)      63.546(3)  
!    65   64.9277937(19)     30.83(3)                 
!_____________________________________________________
!Zn  64   63.9291466(18)     48.63(60)     65.39(2)   
!    66   65.9260368(16)     27.90(27)                
!    67   66.9271309(17)     4.10(13)                 
!    68   67.9248476(17)     18.75(51)                
!    70   69.925325(4)       0.62(3)                  
!_____________________________________________________
!Ga  69   68.925581(3)       60.108(9)     69.723(1)  
!    71   70.9247050(19)     39.892(9)                
!_____________________________________________________
!Ge  70   69.9242504(19)     20.84(87)     72.61(2)   
!    72   71.9220762(16)     27.54(34)                
!    73   72.9234594(16)     7.73(5)                  
!    74   73.9211782(16)     36.28(73)                
!    76   75.9214027(16)     7.61(38)                 
!_____________________________________________________
!As  75   74.9215964(18)     100           74.92160(2)
!_____________________________________________________
!Se  74   73.9224766(16)     0.89(4)       78.96(3)   
!    76   75.9192141(16)     9.37(29)                 
!    77   76.9199146(16)     7.63(16)                 
!    78   77.9173095(16)     23.77(28)                
!    80   79.9165218(20)     49.61(41)                
!    82   81.9167000(22)     8.73(22)                 
!
  case ('79Br')
     mass=78.9183376d0
  case ('81Br')
     mass=80.916291d0
  case ('Br')
     mass=79.904d0
!
!Kr  78   77.920386(7)       0.35(1)       83.80(1)   
!    80   79.916378(4)       2.28(6)                  
!    82   81.9134846(28)     11.58(14)                
!    83   82.914136(3)       11.49(6)                 
!    84   83.911507(3)       57.00(4)                 
!    86   85.9106103(12)     17.30(22)                
!
  case ('85Rb')
     mass=84.9117893d0
  case ('87Rb')
     mass=86.9091835d0
  case ('Rb')
     mass=85.4678d0
  case ('84Sr')
     mass=83.913425d0
  case ('86Sr')
     mass=85.9092624d0
  case ('87Sr')
     mass=86.9088793d0
  case ('88Sr')
     mass=87.9056143d0
  case ('Sr')
     mass=87.62d0
  case ('89Y')
     mass=88.9058479d0
  case ('Y')
     mass=88.90585d0
!_____________________________________________________
!Zr  90   89.9047037(23)     51.45(40)     91.224(2)  
!    91   90.9056450(23)     11.22(5)                 
!    92   91.9050401(23)     17.15(8)                 
!    94   93.9063158(25)     17.38(28)                
!    96   95.908276(3)       2.80(9)                  
!_____________________________________________________
  case ('93Nb')
     mass=92.9063775d0
  case ('Nb')
     mass=92.90638d0
!_____________________________________________________
!Mo  92   91.906810(4)       14.84(35)     95.94(1)   
!    94   93.9050876(20)     9.25(12)                 
!    95   94.9058415(20)     15.92(13)                
!    96   95.9046789(20)     16.68(2)                 
!    97   96.9060210(20)     9.55(8)                  
!    98   97.9054078(20)     24.13(31)                
!    100  99.907477(6)       9.63(23)                 
!_____________________________________________________
!Tc  97   96.906365(5)                     [98]       
!    98   97.907216(4)                                
!    99   98.9062546(21)                              
!_____________________________________________________
!Ru  96   95.907598(8)       5.54(14)      101.07(2)  
!    98   97.905287(7)       1.87(3)                  
!    99   98.9059393(21)     12.76(14)                
!    100  99.9042197(22)     12.60(7)                 
!    101  100.9055822(22)    17.06(2)                 
!    102  101.9043495(22)    31.55(14)                
!    104  103.905430(4)      18.62(27)                
!_____________________________________________________
  case ('103Rh')
     mass=102.905504d0
  case('Rh')
     mass=102.90550d0
!_____________________________________________________
!Pd  102  101.905608(3)      1.02(1)       106.42(1)  
!    104  103.904035(5)      11.14(8)                 
!    105  104.905084(5)      22.33(8)                 
!    106  105.903483(5)      27.33(3)                 
!    108  107.903894(4)      26.46(9)                 
!    110  109.905152(12)     11.72(9)                 
!_____________________________________________________
!Ag  107  106.905093(6)      51.839(8)     107.8682(2)
!    109  108.904756(3)      48.161(8)                
!_____________________________________________________
!Cd  106  105.906458(6)      1.25(6)       112.411(8) 
!    108  107.904183(6)      0.89(3)                  
!    110  109.903006(3)      12.49(18)                
!    111  110.904182(3)      12.80(12)                
!    112  111.9027572(30)    24.13(21)                
!    113  112.9044009(30)    12.22(12)                
!    114  113.9033581(30)    28.73(42)                
!    116  115.904755(3)      7.49(18)                 
!_____________________________________________________
!In  113  112.904061(4)      4.29(5)       114.818(3) 
!    115  114.903878(5)      95.71(5)                 
!_____________________________________________________
!Sn  112  111.904821(5)      0.97(1)       118.710(7) 
!    114  113.902782(3)      0.66(1)                  
!    115  114.903346(3)      0.34(1)                  
!    116  115.901744(3)      14.54(9)                 
!    117  116.902954(3)      7.68(7)                  
!    118  117.901606(3)      24.22(9)                 
!    119  118.903309(3)      8.59(4)                  
!    120  119.9021966(27)    32.58(9)                 
!    122  121.9034401(29)    4.63(3)                  
!    124  123.9052746(15)    5.79(5)                  
!_____________________________________________________
!Sb  121  120.9038180(24)    57.21(5)      121.760(1) 
!    123  122.9042157(22)    42.79(5)                 
!_____________________________________________________
!Te  120  119.904020(11)     0.09(1)       127.60(3)  
!    122  121.9030471(20)    2.55(12)                 
!    123  122.9042730(19)    0.89(3)                  
!    124  123.9028195(16)    4.74(14)                 
!    125  124.9044247(20)    7.07(15)                 
!    126  125.9033055(20)    18.84(25)                
!    128  127.9044614(19)    31.74(8)                 
!    130  129.9062228(21)    34.08(62)                
!
  case ('127I')
     mass=126.904468d0
  case ('I')
     mass=126.90447d0
!
!Xe  124  123.9058958(21)    0.09(1)       131.29(2)  
!    126  125.904269(7)      0.09(1)                  
!    128  127.9035304(15)    1.92(3)                  
!    129  128.9047795(9)     26.44(24)                
!    130  129.9035079(10)    4.08(2)                  
!    131  130.9050819(10)    21.18(3)                 
!    132  131.9041545(12)    26.89(6)                 
!    134  133.9053945(9)     10.44(10)                
!    136  135.907220(8)      8.87(16)                 
!
  case ('133Cs')
     mass=132.905447d0
  case ('Cs')
     mass=132.90545d0
  case ('130Ba')
     mass=129.906310d0
  case ('132Ba')
     mass=131.905056d0
  case ('134Ba')
     mass=133.904503d0
  case ('135Ba')
     mass=134.905683d0
  case ('136Ba')
     mass=135.904570d0
  case ('137Ba')
     mass=136.905821d0
  case ('138Ba')
     mass=137.905241d0
  case ('Ba')
     mass=137.327d0
!_____________________________________________________
!La  138  137.907107(4)      0.090(1)      138.9055(2)
!    139  138.906348(3)      99.910(1)                
!_____________________________________________________
!Ce  136  135.907140(50)     0.185(2)      140.116(1) 
!    138  137.905986(11)     0.251(2)                 
!    140  139.905434(3)      88.450(51)               
!    142  141.909240(4)      11.114(51)               
!_____________________________________________________
  case ('141Pr')
     mass=140.907648d0
  case ('Pr')
     mass=140.90765d0
!_____________________________________________________
!Nd  142  141.907719(3)      27.2(5)       144.24(3)  
!    143  142.909810(3)      12.2(2)                  
!    144  143.910083(3)      23.8(3)                  
!    145  144.912569(3)      8.3(1)                   
!    146  145.913112(3)      17.2(3)                  
!    148  147.916889(3)      5.7(1)                   
!    150  149.920887(4)      5.6(2)                   
!_____________________________________________________
!Pm  145  144.912744(4)                    [145]      
!    147  146.915134(3)                               
!_____________________________________________________
!Sm  144  143.911995(4)      3.07(7)       150.36(3)  
!    147  146.914893(3)      14.99(18)                
!    148  147.914818(3)      11.24(10)                
!    149  148.917180(3)      13.82(7)                 
!    150  149.917271(3)      7.38(1)                  
!    152  151.919728(3)      26.75(16)                
!    154  153.922205(3)      22.75(29)                
!_____________________________________________________
!Eu  151  150.919846(3)      47.81(3)      151.964(1) 
!    153  152.921226(3)      52.19(3)                 
!_____________________________________________________
!Gd  152  151.919788(3)      0.20(1)       157.25(3)  
!    154  153.920862(3)      2.18(3)                  
!    155  154.922619(3)      14.80(12)                
!    156  155.922120(3)      20.47(9)                 
!    157  156.923957(3)      15.65(2)                 
!    158  157.924101(3)      24.84(7)                 
!    160  159.927051(3)      21.86(19)                
!_____________________________________________________
  case ('159Tb')
     mass=158.925343d0
  case ('Tb')
     mass=158.92534d0
!_____________________________________________________
!Dy  156  155.924278(7)      0.06(1)       162.50(3)  
!    158  157.924405(4)      0.10(1)                  
!    160  159.925194(3)      2.34(8)                  
!    161  160.926930(3)      18.91(24)                
!    162  161.926795(3)      25.51(26)                
!    163  162.928728(3)      24.90(16)                
!    164  163.929171(3)      28.18(37)                
!_____________________________________________________
  case ('165Ho')
     mass=164.930319d0
  case ('Ho')
     mass=164.93032d0
!_____________________________________________________
!Er  162  161.928775(4)      0.14(1)       167.26(3)  
!    164  163.929197(4)      1.61(3)                  
!    166  165.930290(3)      33.61(35)                
!    167  166.932045(3)      22.93(17)                
!    168  167.932368(3)      26.78(26)                
!    170  169.935460(3)      14.93(27)                
!_____________________________________________________
  case ('169Tm')
     mass=168.934211d0
  case ('Tm')
     mass=168.93421d0
!_____________________________________________________
  case ('168Yb')
     mass=167.933894d0
  case ('170Yb')
     mass=169.934759d0
  case ('171Yb')
     mass=170.936322d0
  case ('172Yb')
     mass=171.9363777d0
  case ('173Yb')
     mass=172.9382068d0
  case ('174Yb')
     mass=173.9388581d0
  case ('176Yb')
     mass=175.942568d0
  case ('Yb')
     mass=173.04d0
!Yb  168  167.933894(5)      0.13(1)       173.04(3)  
!    170  169.934759(3)      3.04(15)                 
!    171  170.936322(3)      14.28(57)                
!    172  171.9363777(30)    21.83(67)                
!    173  172.9382068(30)    16.13(27)                
!    174  173.9388581(30)    31.83(92)                
!    176  175.942568(3)      12.76(41)                
!_____________________________________________________
!Lu  175  174.9407679(28)    97.41(2)      174.967(1) 
!    176  175.9426824(28)    2.59(2)                  
!_____________________________________________________
!Hf  174  173.940040(3)      0.16(1)       178.49(2)  
!    176  175.9414018(29)    5.26(7)                  
!    177  176.9432200(27)    18.60(9)                 
!    178  177.9436977(27)    27.28(7)                 
!    179  178.9458151(27)    13.62(2)                 
!    180  179.9465488(27)    35.08(16)                
!_____________________________________________________
!Ta  180  179.947466(3)      0.012(2)      180.9479(1)
!    181  180.947996(3)      99.988(2)                
!_____________________________________________________
!W   180  179.946706(5)      0.12(1)       183.84(1)  
!    182  181.948206(3)      26.50(16)                
!    183  182.9502245(29)    14.31(4)                 
!    184  183.9509326(29)    30.64(2)                 
!    186  185.954362(3)      28.43(19)                
!_____________________________________________________
!Re  185  184.9529557(30)    37.40(2)      186.207(1) 
!    187  186.9557508(30)    62.60(2)                 
!_____________________________________________________
!Os  184  183.952491(3)      0.02(1)       190.23(3)  
!    186  185.953838(3)      1.59(3)                  
!    187  186.9557479(30)    1.96(2)                  
!    188  187.9558360(30)    13.24(8)                 
!    189  188.9581449(30)    16.15(5)                 
!    190  189.958445(3)      26.26(2)                 
!    192  191.961479(4)      40.78(19)                
!_____________________________________________________
!Ir  191  190.960591(3)      37.3(2)       192.217(3) 
!    193  192.962924(3)      62.7(2)                  
!_____________________________________________________
!Pt  190  189.959930(7)      0.014(1)      195.078(2) 
!    192  191.961035(4)      0.782(7)                 
!    194  193.962664(3)      32.967(99)               
!    195  194.964774(3)      33.832(10)               
!    196  195.964935(3)      25.242(41)               
!    198  197.967876(4)      7.163(55)                
!
  case ('197Au')
     mass=196.966552d0
  case ('Au')
     mass=196.96655d0
!_____________________________________________________
!Hg  196  195.965815(4)      0.15(1)       200.59(2)  
!    198  197.966752(3)      9.97(20)                 
!    199  198.968262(3)      16.87(22)                
!    200  199.968309(3)      23.10(19)                
!    201  200.970285(3)      13.18(9)                 
!    202  201.970626(3)      29.86(26)                
!    204  203.973476(3)      6.87(15)                 
!_____________________________________________________
!Tl  203  202.972329(3)      29.524(14)    204.3833(2)
!    205  204.974412(3)      70.476(14)               
!_____________________________________________________
!Pb  204  203.973029(3)      1.4(1)        207.2(1)   
!    206  205.974449(3)      24.1(1)                  
!    207  206.975881(3)      22.1(1)                  
!    208  207.976636(3)      52.4(1)                  
!_____________________________________________________
  case ('209Bi')
     mass=208.980383d0
  case ('Bi')
     mass=208.98038d0
!_____________________________________________________
!Po  209  208.982416(3)                    [209]      
!    210  209.982857(3)                               
!_____________________________________________________
!At  210  209.987131(9)                    [210]      
!    211  210.987481(4)                               
!_____________________________________________________
!Rn  211  210.990585(8)                    [222]      
!    220  220.0113841(29)                             
!    222  222.0175705(27)                             
!_____________________________________________________
!Fr  223  223.0197307(29)                  [223]      
!_____________________________________________________
!Ra  223  223.018497(3)                    [226]      
!    224  224.0202020(29)                             
!    226  226.0254026(27)                             
!    228  228.0310641(27)                             
!_____________________________________________________
!Ac  227  227.0277470(29)                  [227]      
!_____________________________________________________
!Th  230  230.0331266(22)                  232.0381(1)
!    232  232.0380504(22)    100                      
!_____________________________________________________
!Pa  231  231.0358789(28)    100           231.03588(2
!_____________________________________________________
!U   233  233.039628(3)                    238.0289(1)
!    234  234.0409456(21)    0.0055(2)                
!    235  235.0439231(21)    0.7200(51)               
!    236  236.0455619(21)                             
!    238  238.0507826(21)    99.2745(106)             
!_____________________________________________________

  case default
     write(*,*) 'Element ',trim(atom),' unknown'
     stop
  end select

  select case (trim(unit))

  case ('amu','a.m.u.','AMU','A.M.U.')

  case default
     mass=mass*conv('M','amu',unit)

  end select

end function mass
