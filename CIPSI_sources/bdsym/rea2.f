      subroutine  rea2(d,m,qend)
      implicitreal*8(a-h,o-p,r-z),logical*4(q)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple
      dimension d(m)
       nn=nn+1
        read(2,end=20)d
       return
20       write(6,*)' fin fichier sur 2,mbuf',m,nn
          if(m-1.gt.ncf)then
          write(6,*)' ncf trouve ',m-1,' > ncf suppose ',ncf
          stop
          else
          ncf=m-1
          write(6,*)' taille matrice ',ncf
            qend=.true.
            m=m-1
          endif
          return
            end
