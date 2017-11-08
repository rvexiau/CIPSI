      subroutine openf
Ccc   modified on 24/09/07 : nomfil() wants a C*10 as a last parameter…
      character*10 info,vcpp
      character*90 blan90,info1,vcpp1
      character*80 prefix
      character*80 f10
      namelist /rcufil/ prefix,info,vcpp,info1,vcpp1,f10
      blan90=' '
      info1=blan90
      vcpp1=blan90
      f10='DERPL03'
      prefix=   '          '
      info=     'info      '
      vcpp=     'vcpp      '
      read(5,rcufil)
      write(6,rcufil)
      call nomfil(2,info,info1,prefix,'old       ')   ! RG 24 sept 2007
      call nomfil(17,vcpp,vcpp1,prefix,'unknown   ')   ! RG 24 sept 2007
!      open (10,file=f10,shared,form='unformatted',
      open (10,file=f10,form='unformatted',
     * status='old')
      return
      end
