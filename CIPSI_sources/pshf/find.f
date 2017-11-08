      subroutine find(a,lng,symb,symbl,nsymb,n)
      character*1 a(20),symb(100,8)
      integer*2 symbl(100)
c      write(6,*) 'symbl',symbl
c      write(6,*) 'symb',symb
c      write(6,*) 'nsymb',nsymb
      do i=1,nsymb
c         write(6,'(1x,20a1)')(symb(i,j),j=1,symbl(i))
      enddo
      n=0
   10 n=n+1
      if (n.gt.nsymb) go to 50
      if (lng.ne.symbl(n)) go to 10
      do 20 i=1,lng
      if (a(i).ne.symb(n,i)) go to 10
   20 continue
      return
   50 write(6,60) (a(i),i=1,lng)
   60 format (//5x,'symbol printed is unresolved'/5x,20a1//)
      stop
      end
