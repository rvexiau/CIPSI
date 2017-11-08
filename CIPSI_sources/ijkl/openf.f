      subroutine openf
      character*80 prefix
      character*10 f25,f50,info,pqrs,ijom,om,scr
      character*90 pqrs1,info1,ijom1,om1,f251,f501,scr1,blan90
      common /op/prefix,om,om1
      namelist /ijfil/ prefix,f25,f50,info,pqrs,ijom,om,scr
     *,pqrs1,info1,ijom1,om1,f251,f501,scr1
      blan90=' '
      pqrs1=blan90
      info1=blan90
      ijom1=blan90
      om1=blan90
      f251=blan90
      f501=blan90
      scr1=blan90
      om='om        '
      f25='f25       '
      f50='f50       '
      info= 'info      '
      pqrs='pqrs      '
      ijom='ijom      '
      scr ='scr      '
      prefix= '          '
      read(5,ijfil)
      write(6,ijfil)
      call nomfil(25,f25,f251,prefix,'new       ')
C Changed by Malte Gross 26/9/96 - didn't run on pqthp5
C     open (unit=40,file=form='unformatted',status='scratch')
      call nomfil(40,scr,scr1,prefix,'new       ')
      call nomfil(50,f50,f501,prefix,'new       ')
      call nomfil(2,info,info1,prefix,'old       ')
      call nomfil(8,pqrs,pqrs1,prefix,'old       ')
      call nomfil(10,ijom,ijom1,prefix,'new       ')
      call nomfil(11,om,om1,prefix,'unknown   ')
      return
      end
