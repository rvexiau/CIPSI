      subroutine openf
      character*80 prefix
      character*10 info,pqrs,pqrscv,infocv,ener,vcpp,om,
     &             vcpp0,vcpp0c
      character*90 blan90,info1,pqrs1,pqrscv1,infocv1,ener1,vcpp1,
     &             om1,vcpp01,vcpp0c1
      logical mixte
      common/cnl/mixte
      namelist/cvafil/prefix,info,pqrs,infocv,pqrscv,ener,vcpp,
     < om,info1,pqrs1,infocv1,pqrscv1,ener1,vcpp1,om1,mixte
      mixte=.false.
      blan90=' '
      pqrs=  'pqrs      '
      info=  'info      '
      infocv='info_cv   '
      pqrscv='pqrs_cv   '
      vcpp=  'vcpp      '
      vcpp0= 'vcpp0     '
      vcpp0c= 'vcpp0c   '
      ener=  'ener      '
      om  =  'om        '
      prefix=' '
      pqrs1=blan90
      info1=blan90
      pqrscv1=blan90
      infocv1=blan90
      vcpp1=blan90
      vcpp01=blan90
      vcpp0c1=blan90
      ener1=blan90
      om1=blan90
      read(5,cvafil)
      write(6,cvafil)
      call nomfil(2,info,info1,prefix,'unknown   ')
      call nomfil(3,infocv,infocv1,prefix,'unknown   ')
      call nomfil(8,pqrs,pqrs1,prefix,'unknown   ')
      call nomfil(9,pqrscv,pqrscv1,prefix,'unknown   ')
      call nomfil(11,om,om1,prefix,'unknown   ')
      call nomfil(17,vcpp,vcpp1,prefix,'unknown   ')
      call nomfil(33,ener,ener1,prefix,'unknown   ')
      if(mixte) then
      call nomfil(18,vcpp0,vcpp01,prefix,'unknown   ')
      call nomfil(19,vcpp0c,vcpp0c1,prefix,'unknown   ')
      endif
      return
      end
