      subroutine jandk                                                  
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 iandj,kandl,same                                        
      common/iofile/ir,iw,ipr,is                                        
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),     
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell                                                            
      common/misc/maxll,iandj,kandl,same                                
      common/shlnos/qq4,inew,jnew,knew,lnew,lit,ljt,lkt,llt,            
     1 mini,minj,mink,minl,maxi,maxj,maxk,maxl,                         
     1 loci,locj,lock,locl,ij,kl,ijkl,nij                               
     1,ni,nj,nk,nl                                                      
      common/times/ti,tx,tim                                            
      common/shlt/ tol,cutoff,icount                                    
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst            
c                                                                       
c     ----- size of gout _                                              
c           1 if s or k shells                                          
c          81 if p      shells                                          
c         256 if      l shells                                          
c        1296 if d or m shells                                          
c       10000 if f or n shells                                          
c       50625 if g or o shells                                          
c                                                                       
      common/gout/gout(10000)                                           
      common/symtry/invt(48),iso(ds,12),ptr(3,144),dtr(6,288)           
     1,ftr(10,480),nt                                                   
      common/wrtc/idum1,idum2,idum3,idum4,indx(4,5),q4,vdum,q4x(5),     
     1 valx(5),nwrt                                                     
      common/output/ nprint                                             
      common/fil2/a(dc),ns(dc),ks(dc),newsh(ds,48)
      dimension m0(48),m1(48),m2(48),m3(48)                             
      dimension mm0(48),mm1(48),mm2(48),mm3(48)                         
 9998 format(1x,'ist,jst,kst,lst',4i4,3x,'nrec,intloc',2i5,3x,'  total t
     1ime',f8.3)                                                        
      call secnd(tim0)                                                  
      ntwd=(nt+3)/4                                                     
      nwrt=0                                                            
c                                                                       
c     ----- initialize parameters                                       
c                                                                       
      call debut                                                        
c
c     --------boucle quadruple sur les elements des petites listes p1
c
      do 13000 in=ist,nshell
      do 1020 it=1,nt
      mm0(it)=newsh(in,it)
      if (mm0(it).gt.in) go to 13000
 1020 continue
c
      do 12000 jn=1,in
      do 1030 it=1,nt
      mm1(it)=newsh(jn,it)
      if (mm1(it).gt.jn) go to 12000
 1030 continue
c
      do 11000 kn=1,in
      do 1040 it=1,nt
      mm2(it)=newsh(kn,it)
      if (mm2(it).gt.kn) go to 11000
 1040 continue
c
      maxln=kn
      if (in.eq.kn) maxln=jn
c
      do 10000 ln=1,maxln
      do 1050 it=1,nt
      mm3(it)=newsh(ln,it)
      if (mm3(it).gt.ln) go to 10000
 1050 continue
c
c      ......boucle sur les quadruplets de transf. de sym.
c
c     ----- ishell -----                                                
c                                                                       
      do 9000 iop=1,nt                                                  
      iii=mm0(iop)
      ip=iop
 1100 ip=ip+1
      if (ip.gt.nt) go to 1110
      if (iii.eq.mm0(ip)) go to 9000
      go to 1100
c                                                                       
c     ----- jshell -----                                                
c                                                                       
 1110 do 8000 jop=1,nt                                                  
      jjj=mm1(jop)
      jp=jop
 1120 jp=jp+1
      if (jp.gt.nt) go to 1130
      if (jjj.eq.mm1(jp)) go to 8000
      go to 1120
c
 1130 iii=mm0(iop)
      if (iii.ge.jjj) go to 1140
      if (in.eq.jn) go to 8000
      nd=iii
      iii=jjj
      jjj=nd
c                                                                       
c     ----- kshell -----                                                
c                                                                       
 1140 do 7000 kop=1,nt                                                  
      kk=mm2(kop)
      kp=kop
 1150 kp=kp+1
      if (kp.gt.nt) go to 1160
      if (kk.eq.mm2(kp)) go to 7000
      go to 1150
c                                                                       
c     ----- lshell -----                                                
c                                                                       
 1160 do 6000 lop=1,nt                                                  
      ll=mm3(lop)
      lp=lop
 1170 lp=lp+1
      if (lp.gt.nt) go to 1200
      if (ll.eq.mm3(lp)) go to 6000
      go to 1170
 1200 ii=iii
c
      jj=jjj
      kk=mm2(kop)
c
      if (kk.ge.ll) go to 1210
      if (kn.eq.ln) go to 6000
      nd=ll
      ll=kk
      kk=nd
 1210 if (ii.gt.kk) go to 1240
      if (ii.eq.kk) go to 1220
      nd=ii
      ii=kk
      kk=nd
      go to 1230
 1220 if (jj.ge.ll) go to 1240
 1230 if (in.eq.kn.and.jn.eq.ln) go to 6000
      nd=jj
      jj=ll
      ll=nd
c
 1240 continue
      do  it=1,nt
      m0(it)=newsh(ii,it)
      m1(it)=newsh(jj,it)
      m2(it)=newsh(kk,it)
      m3(it)=newsh(ll,it)
      end do
c
      n4=0
      do 1570 it=1,nt
      ld=m3(it)
      kd=m2(it)
      if (kd.ge.ld) go to 1520
      nd=kd
      kd=ld
      ld=nd
 1520 id=m0(it)
      jd=m1(it)
      if (id.ge.jd) go to 1530
      nd=id
      id=jd
      jd=nd
 1530 if (id.gt.kd) go to 1560
      if (id.eq.kd) go to 1540
      nd=id
      id=kd
      kd=nd
      go to 1550
 1540 if (jd.ge.ld) go to 1560
 1550 nd=jd
      jd=ld
      ld=nd
 1560 if (id.lt.ii) go to 1570
      if (id.gt.ii) go to 6000
      if (jd.lt.jj) go to 1570
      if (jd.gt.jj) go to 6000
      if (kd.lt.kk) go to 1570
      if (kd.gt.kk) go to 6000
      if (ld.lt.ll) go to 1570
      if (ld.gt.ll) go to 6000
      n4=n4+1
 1570 continue
c                                                                       
c     ----- calculate q4 factor for this group of shells                
c                                                                       
      qq4=dble(float(nt))/dble(float(n4))
c
      inew=ii
      jnew=jj
      knew=kk
      lnew=ll
c                                                                       
c     ----- get information about i-shell and j-shell                   
c                                                                       
      call shells(1)                                                    
c                                                                       
c     ----- form pairs of primitives from i-shell and j-shell           
c                                                                       
      call ijprim                                                       
      if(nij.eq.0)go to 6000 
c                                                                       
c     ----- get information about k-shell and l-shell                   
c                                                                       
      call shells(2)                                                    
      do 80 i=1,ijkl                                                    
   80 gout(i)=0.0d+00                                                   
c                                                                       
c     ----- compute two-electron integrals and                          
c     ----- write them out on tape (is)                                 
c                                                                       
      call genral                                                       
      call timit(0)                                                     
      if(tim.lt.timlim) go to 6000                                      
c                                                                       
c     ----- call final prior to exit from this overlay                  
c                                                                       
      if(nprint.eq.1) call wrt(1)                                       
      call final(0)                                                     
      go to 14000                                                       
 6000 continue                                                          
 7000 continue                                                          
 8000 continue                                                          
 9000 continue                                                          
10000 continue
11000 continue
12000 continue
cjhc  write(iw,*) 'in,nrec,intloc',in,nrec,intloc
      write(iw,*) 'in, nrec, intloc  ',in,nrec,icount-1
13000 continue
      if(nprint.eq.1) call wrt(1)                                       
      call final(1)                                                     
14000 continue                                                          
      call timit(0)  
      call secnd(t)
      t=t-tim0                                                        
      write(iw,9999) t,tim                                              
 9999 format(10x,'elapsed time = ',f10.3,5x,'total time = ',f10.3)      
      return                                                            
      end                                                               
