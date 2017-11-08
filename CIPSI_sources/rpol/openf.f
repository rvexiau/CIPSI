      subroutine openf
      include 'pshf.prm'
      character*10 info,vcpp
      character*90 blan90,info1,vcpp1
      character*80 prefix
      character*80 f10
      namelist /rpofil/ prefix,info,vcpp,info1,vcpp1
      blan90=''
      info1=blan90
      vcpp1=blan90
      prefix=''
      info=  'info        '
      vcpp=  'vcpp        '
      read(5,rpofil)
      write(6,rpofil)
      call nomfil(2,info,info1,prefix,'unknown   ')
      call nomfil(17,vcpp,vcpp1,prefix,'unknown    ')
      return
      end
