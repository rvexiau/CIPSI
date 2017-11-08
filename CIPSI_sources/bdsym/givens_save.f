c$nostandard system
      subroutine givens (h,e,v,n,nev,nvec)
      implicit real*8 (a-h,o-x,z)
      include 'bd.prm'
c
c                  h    matrice a diagonaliser rangee sous la forme d'un
c                       demi matrice superieure
c                  e    valeurs propres en ordre croissant
c                       rangees dans un vecteur de dimension  n
c
c                  v    matrices des vecteurs propres
c                        de dimension nvec
c
c                  n    dimension de la matrice
c                  nev   nombre de valaurs propres demandees
c                  nvec  nombre de vecteurs  propres demandes
c
      integer ag,rn,alf
      real*8 norm,lambd,l,mult
      dimension h(1),e(1),v(1),b(ndetz),c(ndetz),p(ndetz),q(ndetz),
     1r(ndetz),w(ndetz),vy(ndetz+2),in(ndetz),z(ndetz)
      dimension hs1(ndetz),hs2(ndetz)
c***********************************************************************
c     nstart(i)=indice du premier element de la ligne i d'une demi-matri
c     superieure
      nstart(i)=n*(n+1)/2-(n-i+1)*(n-i+2)/2+1
c***********************************************************************
c  soit i  le nombre de decimales pour la machine utilisee. on prend
c  epsi=10**(-(i+3))   et   epsdg=epsv= 10**(-(i-2))
      i=15
      epsi=10.0**(-(i+3))
      epsv=10.0**(-(i-2))
      epsdg=epsv
c     epsi=1.d-10
c     epsdg=1.0d-5
c     epsv=1.d-5
      itevm=3
      joym=15
c       modm1 = 2**15 - 1
      modm1 = 32767
c dlgmax= puissance de 10 maximale
      dlgmax=36.d0
      rn=modm1 - 2
      alf=259
      agnes=2.0**(-14)
      np1=n+1
      nm1=n-1
      nm2=n-2
      ng=1
      if(nm2)99,9,200
  200 continue
c transfert de la premiere ligne de h dans hs1
      do 2 i=1,n
    2 hs1(i)=h(i)
c
c     debut de la tridiagonalisation
c
c     hs1 contient la ligne ir  a tridiagonaliser
c
      do 801 ir=1,nm2
      c(ir)=hs1(ir)
      ip1=ir+1
      ip2=ir+2
c
c     determination des rotations annulant les elements ir+2,..,n
c     le sinus de l'angle de rotation est stocke a la place de l'element
c     annule
      do 802 i=ip2,n
      t=hs1(i)
      if(dabs(t).lt.1.d-8) then
      hs1(i)=1.
      hs2(i-ip1)=0.
      else
      s=hs1(ip1)
      ss=dsqrt(s*s+t*t)
      hs1(ip1)=ss
      s=s/ss
      t=t/ss
      hs1(i)=s
      hs2(i-ip1)=t
      endif
  802 continue
c
c
c     ish3 : indice-ip1 du premier element de la ligne ip1
c
      ish3=nstart(ip1)-ip1
c
      do 805 k=ip2,n
      kp1=k+1
c
c     ish4 :indice-k du premier element de la ligne k
      ish4=nstart(k)-k
c
c
c     on fait les rotations entre les lignes ir+1 et k k=ir+2,...,n
c
      s=hs1(k)
      t=hs2(k-ip1)
      if(t) 820,821,820
  820 continue
      ih4k=ish4+k
      ih3k=ish3+k
      ih3p=ish3+ip1
      u=h(ih3p)
      l=h(ih3k)
      tt=h(ih4k)
c
      h(ih3p)= (u*s+l*t)*s+(l*s+tt*t)*t
      h(ih3k) =-(u*s+l*t)*t+(l*s+tt*t)*s
      h(ih4k) = (u*t-l*s)*t+ (-l*t+tt*s)*s
c
c
      if(kp1-n) 840,840,841
  840 continue
c
      do 807 ll=kp1,n
      ih3l=ish3+ll
      ih4l=ish4+ll
      u=h(ih3l)
      l=h(ih4l)
      h(ih3l)= u*s+l*t
      h(ih4l)=-u*t+l*s
  807 continue
c
  821 continue
c
      if(kp1-n) 823,823,841
  823 continue
      ih3k=ish3+k
      do 808 ll=kp1,n
      s=hs1(ll)
      t=hs2(ll-ip1)
      if(t) 822,808,822
  822 continue
      ih4l=ish4+ll
      u=h(ih3k)
      l=h(ih4l)
      h(ih3k)= u*s+l*t
      h(ih4l)= -u*t+l*s
  808 continue
  841 continue
c
c     fin de rotation de la ligne k
c
  805 continue
      b(ir)=hs1(ip1)
c
c     on stocke dans hs2 le cosinus de la rotation ou le cosinus plus
c     3. (si le sinus est negatif)
      do 3410 k=ip2,n
      s=hs1(k)
      if(s.ge.0.d0) go to 3410
      hs2(k-ip1)=hs2(k-ip1)+3.0d00
 3410 continue
c
c     on transfere la ligne suivante dans hs1
      do 3050 k=ip1,n
 3050 hs1(k)=h(ish3+k)
c   les cosinus des rotations stockes dans hs2 sont transferes dans h
      ish1=nstart(ir)-ir
      do 3020 k=ip2,n
 3020 h(ish1+k)=hs2(k-ip1)
  801 continue
c
      c(nm1)=h(nstart(nm1))
      b(nm1)=h(nstart(nm1)+1)
      c(n)=h(nstart(n))
      b(n)=0.d0
      go to 850
    9 continue
      iiii=0
      c(nm1)=h(iiii+1)
      b(nm1)=h(iiii+2)
      c(n)=h(iiii+3)
      b(n)=0.d0
c
  850 continue
      norm=dabs(c(1))+dabs(b(1))
      do 10 i=2,n
      t=dabs(c(i))+dabs(b(i))+dabs(b(i-1))
      if(norm-t) 104,10,10
  104 norm=t
   10 continue
      do 11 i=1,n
   11 w(i)=b(i)*b(i)
      k=1
      l=-norm
      mult=norm*epsi
      do 12 i=1,nev
   12 e(i)=norm
   13 u=e(k)
   14 lambd=0.5*(l+u)
      if(lambd-l-mult) 30,30,106
  106 if(u-lambd-mult) 30,30,107
  107 ag=0
      i=1
   16 s=c(i)-lambd
   18 if(dabs(s)-1.d-21) 20,108,108
  108 if(s) 110,110,109
  109 ag=ag+1
  110 i=i+1
      if(i-n) 111,111,22
  111 s=c(i)-lambd-w(i-1)/s
      go to 18
   20 ag=ag+1
      i=i+2
      if(i-n) 16,16,22
   22 if(ag-n+k) 24,24,112
  112 l=lambd
      go to 14
   24 u=lambd
      if(n-ag-nev) 113,113,114
  113 m=n-ag
      go to 115
  114 m=nev
  115 do 26 i=k,m
   26 e(i)=lambd
      go to 14
   30 e(k)=lambd
      k=k+1
      if(k-nev) 13,13,116
  116 if(nvec) 40,999,40
   40 ii=0
      inn=-n
      epsin=norm*epsdg
      do 90 i=1,nvec
      inn=inn+n
      ieps=0
      vy(n+1)=0.
      vy(n+2)=0.
  402 t=e(i)
      ieps=ieps+1
      go to (430,404),ieps
  404 j=i-ii
      if(j) 406,406,407
  406 s=norm
      go to 408
  407 s=t-e(j)
  408 l=t
      do 414 k=i,nev
      u=e(k)
      if(u-l - epsin) 410,416,416
  410 l=u
  414 continue
  416 tt=u-t
      if(tt-s) 418,420,420
  418 s=tt
  420 u=dabs(t)
      if(u-s) 422,424,424
  422 u=s
  424 t=t + u*epsdg*0.1
      ir=0
  430 continue
      do 44 j=1,n
      p(j)=0.
      q(j)=b(j)
      r(j)=c(j)-t
      go to (42,44),ieps
   42 vy(j)=1.0
   44 continue
      do 50 j=1,nm1
      if(dabs(r(j))+dabs(b(j))) 152,152,154
  152 in(j)=0
      w(j)=0.
      r(j)=1.d-30
      go to 50
  154 continue
      if(dabs(r(j))-dabs(b(j)))49,117,117
  117 mult=b(j)/r(j)
      in(j)=0
      go to 48
   49 mult=r(j)/b(j)
      in(j)=1
      r(j)=b(j)
      t=r(j+1)
      r(j+1)=q(j)
      q(j)=t
      p(j)=q(j+1)
      q(j+1)=0.
   48 w(j)=mult
      q(j+1)=q(j+1)-mult*p(j)
      r(j+1)=r(j+1)-mult*q(j)
      if(dabs(r(j))-1.d-30) 118,118,50
  118 r(j)=1.d-30
   50 continue
      if(dabs(r(n))-1.d-30)  119,119,120
  119 r(n)=1.d-30
  120 go to (1202,145),ieps
 1202 if(i-nvec) 155,121,121
  155 if(dabs(e(i+1)-e(i)) - epsin) 53,121,121
  121 if(ii) 122,55,122
  122 kk=1
      go to 54
   53 kk=2
   54 ii=ii+1
      if(ii-1) 123,55,123
  123 joy=0
      ir=0
   51 joy=joy+1
      if(joy-joym) 124,95,95
  124 do 52 j=1,n
      vy(j)=ran(alf)
   52 continue
   55 itev=0
   56 itev=itev+1
      do 66 ji=1,n
      k=n-ji+1
   62 t=vy(k)
      tnum= t-vy(k+1)*q(k)-vy(k+2)*p(k)
      tden=r(k)
      if(dabs(tnum).gt.1.d-30) then
      dltnum=dlog(dabs(tnum))
      else
      dltnum=-30.d0
      end if
      if((dltnum-dlog(dabs(tden))).gt.dlgmax) then
      do 64 j=1,n
   64 vy(j)=vy(j)*1.d-5
      go to 62
      end if
      vy(k)=tnum/tden
   66 continue
      if(itev-1) 145,145,131
  131 s=e(i)
      tt=dabs((c(1)-s)*vy(1) + b(1)*vy(2))
      t=dabs(vy(1))
      do 135 j=2,n
      u =dabs(b(j-1)*vy(j-1) + (c(j)-s)*vy(j) + b(j)*vy(j+1))
      if(tt-u ) 132,133,133
  132 tt=u
  133 u =dabs(vy(j))
      if(t-u ) 134,135,135
  134 t=u
  135 continue
      s=tt/t
      if(s-epsv) 136,136,139
  136 if(itev-2) 69,69,137
  137 write(6,100) i,itev
      go to 69
  139 write(6,101) i,epsv,itev ,tt,t,s
      if(itev-itevm) 145,69,69
  145 do 68 j=1,nm1
      if(in(j)) 144,144,67
  144 vy(j+1)=vy(j+1)-w(j)*vy(j)
      go to 68
   67 t=vy(j)
      vy(j)=vy(j+1)
      vy(j+1)=t-w(j)*vy(j+1)
   68 continue
      go to 56
   69 continue
      t=1./t
      do 98 j=1,n
   98 vy(j)=vy(j)*t
      if(ii-1) 77,147,72
   72 ji=i-ii+1
      m=i-1
      t=0.
      do 70 j=1,n
      u=vy(j)
   70 t=t + u*u
      t=1./t
      ag=1
      ikk=n*(ji-2)
      do 75 k=ji,m
      ikk=ikk+n
      s=0.
      do 73 j=1,n
      jvk=ikk+j
   73 s=s+vy(j)*v(jvk)
      tt=s*z(k)
      do 74 j=1,n
      jvk=ikk+j
   74 vy(j) = vy(j)-tt*v(jvk)
      u=s*tt*t
      t=t/(1.-u)
      if(u-0.75) 75,75,76
   76 ag=2
   75 continue
      go to (160,143),ag
  143 itev=1
      ir=ir+1
      go to (145,1432),ir
 1432 go to (402,51),ieps
  160 if(joy-2) 147,147,146
  146 write(6,100) i,ii,joy,ir
  147 s=0.
      do 162 j=1,n
      jvi=inn+j
      u=vy(j)
      s=s + u*u
  162 v(jvi)=u
      z(i)=1./s
c
c
      go to (78,90),kk
   77 do 777 j=1,n
      jvi=inn+j
  777 v(jvi)=vy(j)
   78 continue
      ii=0
   90 continue
      nbl1=0
      if(nm2)85,85,81
   81 do 84 j=1,nm2
      k=n-j-1
      ks=j+2
      kp=k+1
      m=k+2
      ir=n+1
      is=j+1
      ish1=nstart(k)-k
c
      do 83 kk=m,n
      is=is-1
      ir=ir-1
      t=h(ish1+ir)
      if(t.gt.1.5d0) go to 3460
      if(t.eq.0.d0) go to 83
      s=dsqrt(1.d0-t*t)
      go to 3470
 3460 t=t-3.0d0
      s=-dsqrt(1.0d0-t*t)
 3470 continue
      inn=-n
      do 89 joy=1,nvec
      inn=inn+n
      kpjoy=inn+kp
      irjoy=inn+ir
      u=v(kp joy)
      l=v(ir joy)
      v(kp joy)=u*s-l*t
      v(ir joy)=u*t+l*s
   89 continue
   83 continue
   84 continue
   85 continue
      inn=-n
      do 87 joy=1,nvec
      inn=inn+n
      s=0.
      do 86 j=1,n
      jjoy=inn+j
   86 s=s+v(j joy)*v(j joy)
      s=1./dsqrt(s)
      do 88 j=1,n
      jjoy=inn+j
   88 v(j joy)=v(j joy)*s
   87 continue
  999 return
   95 write(6,102)i,ii,joym
      call exit
  100 format(1h0,10i5)
  101 format(1h0,i4,e16.8,i3,3e18.8)
  102 format(1h0,i4,i3)
      return
 99   continue
      e(1)=h(1)
      v(1)=1.0
      return
      end
