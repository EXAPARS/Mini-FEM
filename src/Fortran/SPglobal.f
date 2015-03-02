!-----------------------------------------------------------------------
!
! Global operations for double
!
!-----------------------------------------------------------------------
      subroutine SPGdSUM(value,nbval,buffer)
      implicit   none
      include   'para.h'
      integer    nbval
      double precision value(nbval),buffer(nbval)
!
      integer    i
!
      if(nblk.le.1) return
!
      do i=1,nbval
        buffer(i) = value(i)
      enddo
      sendcount = nbval
      call MPI_Allreduce(buffer,value,sendcount,MPI_double_precision,   &
     &                   MPI_sum,MPI_comm_myw,MPI_err)
      if(MPI_err.ne.0) call SP_stop('SPGdSUM/Allreduce',MPI_err)
!
      return
      end
!-----------------------------------------------------------------------
!
! Global operations for Logical
!
!-----------------------------------------------------------------------
      subroutine SPGlOR(value,nbval,buffer)
      implicit   none
      include   'para.h'
      integer    nbval
      logical    value(nbval),buffer(nbval)
!
      integer    i
!
      if(nblk.le.1) return
!
      do i=1,nbval
        buffer(i) = value(i)
      enddo
      sendcount = nbval
      
!      call MPI_Allreduce(buffer,value,sendcount,MPI_logical,MPI_lor,    &
!     &                   MPI_comm_myw,MPI_err)
      call MPI_Allreduce(buffer,value,sendcount,MPI_logical,MPI_lor,    &
     &                   MPI_COMM_WORLD,MPI_err)

      if(MPI_err.ne.0) call SP_stop('SPGlOR/Allreduce',MPI_err)
!
      return
      end
