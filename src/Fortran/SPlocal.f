!-----------------------------------------------------------------------
!
! Local communication for Double Real Data
!
!-----------------------------------------------------------------------
      subroutine SPdOPE(nit,lni,npara,inni,nni,ndof,nnt,ldv,v,b,ope)
      implicit   none
      include   'para.h'
      character  ope*3
      integer    nit,lni,ndof,nnt,ldv
      integer    npara(nit,3),inni(nit+1),nni(lni)
      double precision v(ldv,*),b(ndof,lni,2)
      integer    dtype
      parameter (dtype=MPI_double_precision)
!
      intrinsic  MIN,MAX
      integer    iop,intf,ni,ni1,ni2,idof,i,nvproc,itype
      logical    lerr(2),lbuf(2),fv
!
      integer    give,get
      data       give,get                                          /1,2/
! NM
      integer ierror
      call MPI_COMM_RANK(MPI_COMM_WORLD, iblk, ierror)
      iblk = iblk + 1
! end NM
!
! 2.  Initializing reception from adjacent blocks
!     *) Type of message from block K = 100 + K
!
      ni2    = 0
      do intf=1,nit
        ni1    = ni2 + 1
        ni2    = inni(intf+1)
        recvcount = (ni2 - ni1 + 1)*ndof
        nvproc = npara(intf,1) - 1
        itype  = 100 + npara(intf,1)
        call MPI_Irecv(b(1,ni1,get),recvcount,dtype,nvproc,itype,       &
     &                 MPI_COMM_WORLD,npara(intf,3),MPI_err)
      enddo
!
! 3.  Buffering local data
!
      ni2    = 0
      do intf=1,nit
        ni1    = ni2 + 1
        ni2    = inni(intf+1)
        do idof=1,ndof
          do ni=ni1,ni2
            i      = nni(ni)
            b(idof,ni,give) = v(idof,i)
          enddo
        enddo
      enddo
!
! 4.  Sending local data to adjacent blocks
!     *) Type of all messages = 100 + local block number
!
      itype  = 100 + iblk
      ni2    = 0
      do intf=1,nit
        ni1    = ni2 + 1
        ni2    = inni(intf+1)
        sendcount = (ni2 - ni1 + 1)*ndof
        nvproc = npara(intf,1)-1
        call MPI_Send(b(1,ni1,give),sendcount,dtype,nvproc,itype,       &
     &                MPI_COMM_WORLD,MPI_err)
      enddo
!
! 5.  Waiting/Checking incoming data
!
      do intf=1,nit
        call MPI_Wait(npara(intf,3),MPI_status,MPI_err)
      enddo
!
! 5.  Assembling local and incoming data
!
      ni2    = 0
      do intf=1,nit
        ni1    = ni2 + 1
        ni2    = inni(intf+1)
        do idof=1,ndof
          do ni=ni1,ni2
            i      = nni(ni)
            v(idof,i) = v(idof,i) + b(idof,ni,get)
          enddo
        enddo
      enddo
!
      return
80001 format(i3.3,': ',a,5i14)
      end
