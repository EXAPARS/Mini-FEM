      include  'mpif.h'
!
      integer   MPI_status
      common /SPWsta/ MPI_status(MPI_status_size)
!
!     MPI_comm_myw : communicator for my World
!     nblk         : total number of blocks
!     iblk         : rank of local block
!     MPI_out      : output file number
!                    > 0 : local outputs => [PATH/]output[_iblk]
!                    <=0 : local outputs => STDOUT.
!
      integer   MPI_comm_myw,nblk,iblk,MPI_out
      common /SPWmyw/ MPI_comm_myw,nblk,iblk,MPI_out
!
      integer   MPI_err,sendcount,recvcount
