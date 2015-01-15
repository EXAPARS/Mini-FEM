      subroutine lap_comm_mpi(nnt,prec,nblk,pbuf,nit,lni,npara,inni,nni)
      
      implicit   none
      integer    nnt,nit,lni,npara,inni,nni,nblk
      real*8     prec,pbuf
      dimension  npara(nit,3),inni(nit+1),nni(lni)
      dimension  prec(nnt),pbuf(*)

      if(nblk > 1) then
        call SPdOPE(nit,lni,npara,inni,nni,+1,nnt,nnt,prec,pbuf,'sum')
      endif      

      return
      end
