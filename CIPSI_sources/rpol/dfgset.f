      subroutine dfgset
c
c***********************************************************************
c                     ----- set pure D F G coeff -----
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
      parameter (ZERO=0.0D+00)
c
      common /dfgcof/ dcoef(6,6),fcoef(10,10),gcoef(15,15)
      common /cofdfg/ dcart(6,6),fcart(10,10),gcart(15,15)
c
      do j=1,6
         do i=1,6
            dcoef(i,j) = ZERO
            dcart(i,j) = ZERO
         enddo
      enddo
      dcoef( 1, 1) =-0.500000000000020D+00
      dcoef( 2, 1) =-0.500000000000020D+00
      dcoef( 3, 1) = 0.100000000000000D+01
      dcoef( 5, 2) = 0.100000000000000D+01
      dcoef( 6, 3) = 0.100000000000000D+01
      dcoef( 1, 4) = 0.866025403784470D+00
      dcoef( 2, 4) =-0.866025403784470D+00
      dcoef( 4, 5) = 0.100000000000000D+01
      dcoef( 1, 6) = 0.447213595499940D+00
      dcoef( 2, 6) = 0.447213595499940D+00
      dcoef( 3, 6) = 0.447213595499940D+00
      dcart( 1, 1) =-0.333333333333329D+00
      dcart( 4, 1) = 0.577350269189605D+00
      dcart( 6, 1) = 0.745355992499950D+00
      dcart( 1, 2) =-0.333333333333329D+00
      dcart( 4, 2) =-0.577350269189605D+00
      dcart( 6, 2) = 0.745355992499950D+00
      dcart( 1, 3) = 0.666666666666658D+00
      dcart( 6, 3) = 0.745355992499980D+00
      dcart( 5, 4) = 0.100000000000000D+01
      dcart( 2, 5) = 0.100000000000000D+01
      dcart( 3, 6) = 0.100000000000000D+01
      do j=1,10
         do i=1,10
            fcoef(i,j) = ZERO
         enddo
      enddo
      fcoef( 3, 1) = 0.999999999999980D+00
      fcoef( 5, 1) =-0.670820393249900D+00
      fcoef( 7, 1) =-0.670820393249900D+00
      fcoef( 1, 2) =-0.612372435695820D+00
      fcoef( 6, 2) =-0.273861278752580D+00
      fcoef( 8, 2) = 0.109544511501030D+01
      fcoef( 2, 3) =-0.612372435695820D+00
      fcoef( 4, 3) =-0.273861278752580D+00
      fcoef( 9, 3) = 0.109544511501030D+01
      fcoef( 5, 4) = 0.866025403784470D+00
      fcoef( 7, 4) =-0.866025403784470D+00
      fcoef(10, 5) = 0.100000000000000D+01
      fcoef( 1, 6) = 0.790569415042100D+00
      fcoef( 6, 6) =-0.106066017177980D+01
      fcoef( 2, 7) =-0.790569415042100D+00
      fcoef( 4, 7) = 0.106066017177980D+01
      fcoef( 1, 8) = 0.654653670707980D+00
      fcoef( 6, 8) = 0.292770021884560D+00
      fcoef( 8, 8) = 0.292770021884560D+00
      fcoef( 2, 9) = 0.654653670707980D+00
      fcoef( 4, 9) = 0.292770021884560D+00
      fcoef( 9, 9) = 0.292770021884560D+00
      fcoef( 3,10) = 0.654653670707980D+00
      fcoef( 5,10) = 0.292770021884560D+00
      fcoef( 7,10) = 0.292770021884560D+00
      fcart( 2, 1) =-0.244948974278320D+00
      fcart( 6, 1) = 0.316227766016839D+00
      fcart( 8, 1) = 0.916515138991149D+00
      fcart( 3, 2) =-0.244948974278320D+00
      fcart( 7, 2) =-0.316227766016839D+00
      fcart( 9, 2) = 0.916515138991149D+00
      fcart( 1, 3) = 0.400000000000015D+00
      fcart(10, 3) = 0.916515138991153D+00
      fcart( 3, 4) =-0.182574185835062D+00
      fcart( 7, 4) = 0.707106781186559D+00
      fcart( 9, 4) = 0.683130051063977D+00
      fcart( 1, 5) =-0.447213595499977D+00
      fcart( 4, 5) = 0.577350269189605D+00
      fcart(10, 5) = 0.683130051063986D+00
      fcart( 2, 6) =-0.182574185835062D+00
      fcart( 6, 6) =-0.707106781186559D+00
      fcart( 8, 6) = 0.683130051063977D+00
      fcart( 1, 7) =-0.447213595499977D+00
      fcart( 4, 7) =-0.577350269189605D+00
      fcart(10, 7) = 0.683130051063986D+00
      fcart( 2, 8) = 0.730296743340235D+00
      fcart( 8, 8) = 0.683130051064003D+00
      fcart( 3, 9) = 0.730296743340235D+00
      fcart( 9, 9) = 0.683130051064003D+00
      fcart( 5,10) = 0.100000000000000D+01
      do j=1,15
         do i=1,15
            gcoef(i,j) = ZERO
         enddo
      enddo
      gcoef( 1, 1) = 0.375000000000070D+00
      gcoef( 2, 1) = 0.375000000000070D+00
      gcoef( 3, 1) = 0.100000000000020D+01
      gcoef(10, 1) = 0.219577516413420D+00
      gcoef(11, 1) =-0.878310065653690D+00
      gcoef(12, 1) =-0.878310065653690D+00
      gcoef( 5, 2) =-0.896421457000670D+00
      gcoef( 8, 2) = 0.119522860933420D+01
      gcoef(14, 2) =-0.400891862868570D+00
      gcoef( 7, 3) =-0.896421457000670D+00
      gcoef( 9, 3) = 0.119522860933420D+01
      gcoef(13, 3) =-0.400891862868570D+00
      gcoef( 1, 4) =-0.559016994375070D+00
      gcoef( 2, 4) = 0.559016994375070D+00
      gcoef(11, 4) = 0.981980506062010D+00
      gcoef(12, 4) =-0.981980506062010D+00
      gcoef( 4, 5) =-0.422577127364270D+00
      gcoef( 6, 5) =-0.422577127364270D+00
      gcoef(15, 5) = 0.113389341902770D+01
      gcoef( 5, 6) = 0.790569415042100D+00
      gcoef(14, 6) =-0.106066017177980D+01
      gcoef( 7, 7) =-0.790569415042100D+00
      gcoef(13, 7) = 0.106066017177980D+01
      gcoef( 1, 8) = 0.739509972887590D+00
      gcoef( 2, 8) = 0.739509972887590D+00
      gcoef(10, 8) =-0.129903810567670D+01
      gcoef( 4, 9) = 0.111803398874970D+01
      gcoef( 6, 9) =-0.111803398874970D+01
      gcoef( 1,10) = 0.333333333333350D+00
      gcoef( 2,10) = 0.333333333333350D+00
      gcoef( 3,10) = 0.333333333333350D+00
      gcoef(10,10) = 0.195180014589690D+00
      gcoef(11,10) = 0.195180014589690D+00
      gcoef(12,10) = 0.195180014589690D+00
      gcoef( 1,11) =-0.372677996249970D+00
      gcoef( 2,11) =-0.372677996249970D+00
      gcoef( 3,11) = 0.745355992499940D+00
      gcoef(10,11) =-0.218217890235960D+00
      gcoef(11,11) = 0.109108945117980D+00
      gcoef(12,11) = 0.109108945117980D+00
      gcoef( 7,12) = 0.487950036474280D+00
      gcoef( 9,12) = 0.487950036474280D+00
      gcoef(13,12) = 0.218217890235990D+00
      gcoef( 5,13) = 0.487950036474280D+00
      gcoef( 8,13) = 0.487950036474280D+00
      gcoef(14,13) = 0.218217890235990D+00
      gcoef( 4,14) = 0.487950036474280D+00
      gcoef( 6,14) = 0.487950036474280D+00
      gcoef(15,14) = 0.218217890235990D+00
      gcoef( 1,15) = 0.645497224367910D+00
      gcoef( 2,15) =-0.645497224367910D+00
      gcoef(11,15) = 0.188982236504580D+00
      gcoef(12,15) =-0.188982236504580D+00
      gcart( 1, 1) = 0.857142857142726D-01
      gcart( 4, 1) =-0.127775312999958D+00
      gcart( 8, 1) = 0.169030850945672D+00
      gcart(10, 1) = 0.599999999999963D+00
      gcart(11, 1) =-0.383325938999961D+00
      gcart(15, 1) = 0.663940002206980D+00
      gcart( 1, 2) = 0.857142857142726D-01
      gcart( 4, 2) = 0.127775312999958D+00
      gcart( 8, 2) = 0.169030850945672D+00
      gcart(10, 2) = 0.599999999999963D+00
      gcart(11, 2) =-0.383325938999961D+00
      gcart(15, 2) =-0.663940002206980D+00
      gcart( 1, 3) = 0.228571428571389D+00
      gcart(10, 3) = 0.599999999999959D+00
      gcart(11, 3) = 0.766651877999913D+00
      gcart( 5, 4) =-0.169030850945695D+00
      gcart( 9, 4) = 0.447213595500036D+00
      gcart(14, 4) = 0.878310065653659D+00
      gcart( 2, 5) =-0.358568582800373D+00
      gcart( 6, 5) = 0.316227766016834D+00
      gcart(13, 5) = 0.878310065653649D+00
      gcart( 5, 6) =-0.169030850945695D+00
      gcart( 9, 6) =-0.447213595500036D+00
      gcart(14, 6) = 0.878310065653659D+00
      gcart( 3, 7) =-0.358568582800373D+00
      gcart( 7, 7) =-0.316227766016834D+00
      gcart(12, 7) = 0.878310065653649D+00
      gcart( 2, 8) = 0.478091443733830D+00
      gcart(13, 8) = 0.878310065653668D+00
      gcart( 3, 9) = 0.478091443733830D+00
      gcart(12, 9) = 0.878310065653668D+00
      gcart( 1,10) = 0.975900072948535D-01
      gcart( 8,10) =-0.577350269189606D+00
      gcart(10,10) = 0.683130051064037D+00
      gcart(11,10) =-0.436435780472049D+00
      gcart( 1,11) =-0.390360029179404D+00
      gcart( 4,11) = 0.436435780471966D+00
      gcart(10,11) = 0.683130051064050D+00
      gcart(11,11) = 0.218217890236032D+00
      gcart(15,11) = 0.377964473009290D+00
      gcart( 1,12) =-0.390360029179404D+00
      gcart( 4,12) =-0.436435780471966D+00
      gcart(10,12) = 0.683130051064050D+00
      gcart(11,12) = 0.218217890236032D+00
      gcart(15,12) =-0.377964473009290D+00
      gcart( 3,13) =-0.267261241912473D+00
      gcart( 7,13) = 0.707106781186563D+00
      gcart(12,13) = 0.654653670707971D+00
      gcart( 2,14) =-0.267261241912473D+00
      gcart( 6,14) =-0.707106781186563D+00
      gcart(13,14) = 0.654653670707971D+00
      gcart( 5,15) = 0.755928946018445D+00
      gcart(14,15) = 0.654653670707969D+00
      return
      end