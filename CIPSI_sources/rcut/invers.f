      subroutine invers
      implicit integer (a-z)
      double precision temp,xi,yi,zi,xj,yj,zj,cx,cy,cz
      common/parinv/xi,yi,zi,xj,yj,zj,cx,cy,cz,lit,ljt,i1,j1,i2,j2,
     1              mini,maxi,mazi,minj,maxj,mazj,loci,locj
c
c
c
                temp=xj
                xj=xi
                xi=temp
                temp=yj
                yj=yi
                yi=temp
                temp=zj
                zj=zi
                zi=temp
c
                itemp=ljt
                ljt=lit
                lit=itemp
c
                itemp=j1
                j1=i1
                i1=itemp
                itemp=j2
                j2=i2
                i2=itemp
c
                itemp=minj
                minj=mini
                mini=itemp
                itemp=maxj
                maxj=maxi
                maxi=itemp
                itemp=mazj
                mazj=mazi
                mazi=itemp
c
                itemp=locj
                locj=loci
                loci=itemp
c
c
c
      return
      end
