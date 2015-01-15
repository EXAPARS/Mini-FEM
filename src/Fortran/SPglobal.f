!-----------------------------------------------------------------------
!
! Global operations for Double Real
!
!-----------------------------------------------------------------------
      subroutine SPGdMAX(value,nbval,buffer)
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
     &                   MPI_max,MPI_comm_myw,MPI_err)
      if(MPI_err.ne.0) call SP_stop('SPGdMAX/Allreduce',MPI_err)
!
      return
      end
!-----------------------------------------------------------------------
      subroutine SPGdMIN(value,nbval,buffer)
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
     &                   MPI_min,MPI_comm_myw,MPI_err)
      if(MPI_err.ne.0) call SP_stop('SPGdMIN/Allreduce',MPI_err)
!
      return
      end
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
! Global operations for Integer
!
!-----------------------------------------------------------------------
      subroutine SPGiMAX(value,nbval,buffer)
      implicit   none
      include   'para.h'
      integer    nbval
      integer    value(nbval),buffer(nbval)
!
      integer    i
!
      if(nblk.le.1) return
!
      do i=1,nbval
        buffer(i) = value(i)
      enddo
      sendcount = nbval
      call MPI_Allreduce(buffer,value,sendcount,MPI_integer,MPI_max,    &
     &                   MPI_comm_myw,MPI_err)
      if(MPI_err.ne.0) call SP_stop('SPGiMAX/Allreduce',MPI_err)
!
      return
      end
!-----------------------------------------------------------------------
      subroutine SPGiSUM(value,nbval,buffer)
      implicit   none
      include   'para.h'
      integer    nbval
      integer    value(nbval),buffer(nbval)
!
      integer    i
!
      if(nblk.le.1) return
!
      do i=1,nbval
        buffer(i) = value(i)
      enddo
      sendcount = nbval
      call MPI_Allreduce(buffer,value,sendcount,MPI_integer,MPI_sum,    &
     &                   MPI_comm_myw,MPI_err)
      if(MPI_err.ne.0) call SP_stop('SPGiSUM/Allreduce',MPI_err)
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
