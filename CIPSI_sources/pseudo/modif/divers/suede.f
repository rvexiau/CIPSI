        character*80 ident,psname
         character*2 rep
         real*8 work(10000)
         open(unit=9,form='unformatted',file='pot')
         open(unit=10,form='formatted',file='potclair')

10         write(6,*)' nom du fichier de non locaux?, "  " to stop'
        psname='                                          '
         read(5,'(a)') psname
         if(psname(1:4).eq.'    ')then
             write(9)psname,0,0,0
             stop
         endif
        open(unit=8,form='unformatted',file=psname)
        rep='   '
1       read(8,end=2,err=4)ident,l,nm1,nm2,
     *           (work(i),i=1,nm1+nm1+nm2+nm2)
         write(6,*)' label: ',ident
         write(6,*)' nproj ',l,' m1=',nm1,' m2=',nm2
        if(l.lt.0)goto10
       write(6,*)' exp',
     *           (work(i),i=1,nm1+nm2)
       write(6,*)' coeff',
     *       (work(i),i=nm1+nm2+1,2*(nm1+nm2))
         if(index(rep,'A').eq.0.and.index(rep,'a').eq.0)then
         write(6,*)' selection Yes or Not . or All ?'
         read(5,'(a)') rep
         if(index(rep,'z').ne.0)stop
         if(index(rep,'A').ne.0.or.index(rep,'a').ne.0)rep='AY'
         endif

          if(index(rep,'Y').ne.0.or.index(rep,'y').ne.0)then
                 write(6,*)' selected: ',ident
            if(ident(1:1).ne.' ')then
            do i=80,2,-1
            ident(i:i)=ident(i-1:i-1)
             enddo
            ident(1:1)=' '
            endif
    
       write(9)ident,l,nm1,nm2,
     *           (work(i),i=1,nm1+nm2),
     *       (work(i),i=nm1+nm2+1,2*(nm1+nm2))
       write(10,*)ident,l,nm1,nm2,
     *           (work(i),i=1,nm1+nm2),
     *       (work(i),i=nm1+nm2+1,2*(nm1+nm2))
           do j=0,l
          read(8,end=3)k,nz,nb,
     *   (work(i),i=1,nb+nz+nz*nb)
           write(6,*) ' ang ',j,k,' nzeta',nz,' ncont',nb
          write(9)k,nz,nb,
     *   (work(i),i=1,nb+nz+nz*nb)
          write(10,*)k,nz,nb
          write(10,112)(work(i),i=1,nb)
          write(10,112)(work(i),i=nb+1,nb+nz)
	    do iiprim=0,nz-1
       write(10,112)(work(i),i=nb+nz+1+iiprim,nb+nz+nb*nz,nz)
             enddo
	    enddo
 111      format(2x,1e20.14)
 112      format(3d20.13)
         else

           do j=0,l
          read(8,end=3,err=3)k,nz,nb,
     *   (work(i),i=1,nb+nz+nz*nb)
           write(6,*) ' ang ',j,k,' nzeta',nz,' ncont',nb
           enddo
         endif
         goto1
2        continue
        goto10
3        write(6,*)'error  end of file for : '
           write(6,*) ' ang ',j,k,' nzeta',nz,' ncont',nb
         close (unit=9)
         close(unit=10)
         close(unit=8)
         stop
4         write(6,*)'error  end of file for : '
           write(6,*) ident
         end
