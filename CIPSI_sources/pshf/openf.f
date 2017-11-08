      subroutine openf
      character*10 info,om,om_essai,pqrs,ener
      character*90 blan90,info1,om1,om_essai1,pqrs1,ener1
      character*80 psnl_fil,ps_molcas_fil,egp_fil,prefix
C Changed 26/11/97 Malte GROSS
      include 'pshf.prm'
      namelist /psfil/ prefix,info,om,om_essai,pqrs,ener,
     *info1,om1,om_essai1,pqrs1,ener1,psnl_fil,egp_fil,
     *ps_molcas_fil
      blan90=' '
      psnl_fil=psnl_fil_def
      ps_molcas_fil=ps_molcas_fil_def
      egp_fil=egp_fil_def
      info1=blan90
      pqrs1=blan90
      om1=blan90
      om_essai1=blan90
      ener1=blan90
      prefix=   '          '
      info=     'info      '
      pqrs=     'pqrs      '
      om =      'om        '
      om_essai= 'om_essai  '
      ener=     'ener      '
      read(5,psfil)
      write(6,psfil)
      call nomfil(2,info,info1,prefix,'unknown   ')
      call nomfil(8,pqrs,pqrs1,prefix,'unknown   ')
      call nomfil(11,om,om1,prefix,'unknown   ')
      call nomfil(12,om_essai,om_essai1,prefix,'unknown   ')
      call nomfil(33,ener,ener1,prefix,'unknown   ')
      open (20,form='unformatted',status='old',file=psnl_fil)
      open (23,form='unformatted',status='old',file=ps_molcas_fil)
c     open (21,form='unformatted',status='old',file=egp_fil)
      write(6,*)'fin ouverture des files'
      return
      end
