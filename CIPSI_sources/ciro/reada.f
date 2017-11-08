      subroutine reada (a,n,idaf)                                       
      implicit real*8 (a-h,o-z)                                         
      common/iofile/ir,iw,ip,is,ipk,idf,nw,ioda(20)                     
      dimension a(n)                                                    
      nav=ioda(idaf)                                                    
      read(idf,rec=nav)a                                                
      return                                                            
      end                                                               
