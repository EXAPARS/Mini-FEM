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
      integer    iop,intf,ni,ni1,ni2,idof,i,nvproc,itype,nerr
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
! 1.  Checking inputs
!
      if(ndof.le.0) call SP_stop('SPdOPE/[ndof]...',-1)
!
      lerr(1) = nnt.le.ndof
      lerr(2) = ldv.lt.ndof
!      call SPGlOR(lerr,2,lbuf)
!      if(lerr(1).or.lerr(2)) call SP_stop('SPdOPE/[ndof,nnt,ldv]...',-1)
      fv     = ldv.lt.nnt
!
      if     (ope.eq.'MIN' .or. ope.eq.'min') then
        iop    =  1
      else if(ope.eq.'MAX' .or. ope.eq.'max') then
        iop    =  2
      else if(ope.eq.'SUM' .or. ope.eq.'sum') then
        iop    =  3
      else
        call SP_stop('SPdOPE/[ope]...',-1)
      endif
!
! 2.  Initializing reception from adjacent blocks
!     *) Type of message from block K = 100 + K
!
      nerr   = 0
      ni2    = 0
      print *, "Fortran"
      do intf=1,nit
        ni1    = ni2 + 1
        ni2    = inni(intf+1)
        recvcount = (ni2 - ni1 + 1)*ndof
        nvproc = npara(intf,1) - 1
        itype  = 100 + npara(intf,1)
        call MPI_Irecv(b(1,ni1,get),recvcount,dtype,nvproc,itype,       &
     &                 MPI_COMM_WORLD,npara(intf,3),MPI_err)
        if(MPI_err.ne.0) nerr = nerr + 1
      enddo
      lerr(1) = nerr.gt.0
!      call SPGlOR(lerr,1,lbuf)
!      if(lerr(1)) call SP_stop('SPdOPE/MPI_Irecv...',-2)
!
! 3.  Buffering local data
!
      ni2    = 0
      if(fv) then
        do intf=1,nit
          ni1    = ni2 + 1
          ni2    = inni(intf+1)
          do idof=1,ndof
            do ni=ni1,ni2
              i      = nni(ni)
              b(idof,ni,give) = v(idof,i)
              if ((intf == 1) .AND. (idof == 1) .AND. (ni == ni1)) then
                print*, iblk-1,": prec(", idof, ",",  i, ") =",v(idof,i)
              endif
            enddo
          enddo
        enddo
      else
        do intf=1,nit
          ni1    = ni2 + 1
          ni2    = inni(intf+1)
          do idof=1,ndof
            do ni=ni1,ni2
              i      = nni(ni)
              b(idof,ni,give) = v(i,idof)
            enddo
          enddo
        enddo
      endif
!
! 4.  Sending local data to adjacent blocks
!     *) Type of all messages = 100 + local block number
!
      itype  = 100 + iblk
      nerr   = 0
      ni2    = 0
      do intf=1,nit
        ni1    = ni2 + 1
        ni2    = inni(intf+1)
        sendcount = (ni2 - ni1 + 1)*ndof
        nvproc = npara(intf,1)-1
!        print *, iblk, ":", ni1, ni2, sendcount, nvproc, itype
        call MPI_Send(b(1,ni1,give),sendcount,dtype,nvproc,itype,       &
     &                MPI_COMM_WORLD,MPI_err)
        if(MPI_err.ne.0) nerr = nerr + 1
      enddo
      lerr(1) = nerr.gt.0
!      call SPGlOR(lerr,1,lbuf)
!      if(lerr(1)) call SP_stop('SPdOPE/MPI_Send...',-4)
!
! 5.  Waiting/Checking incoming data
!
      nerr   = 0
      ni2    = 0
      do intf=1,nit
        ni1    = ni2 + 1
        ni2    = inni(intf+1)
        sendcount = (ni2 - ni1 + 1)*ndof
        
        call MPI_Wait(npara(intf,3),MPI_status,MPI_err)
        
!        if(MPI_err.eq.0) then
!          call MPI_Get_count(MPI_status,dtype,recvcount,MPI_err)
!          if     (MPI_err.ne.0) then
!            nerr   = nerr + 1
!          else if(recvcount.ne.sendcount) then
!            nerr   = nerr + 1
!          endif
!        else
!          nerr   = nerr + 1
!        endif
      enddo
      lerr(1) = nerr.gt.0
!      call SPGlOR(lerr,1,lbuf)
!      if(lerr(1)) call SP_stop('SPdOPE/MPI_Wait & MPI_Get_count...',-5)
!
! 5.  Assembling local and incoming data
!

      ni2    = 0
      if(fv) then
        if     (iop.eq.1) then
          do intf=1,nit
            ni1    = ni2 + 1
            ni2    = inni(intf+1)
            do idof=1,ndof
              do ni=ni1,ni2
                i      = nni(ni)
                v(idof,i) = MIN(v(idof,i),b(idof,ni,get))
              enddo
            enddo
          enddo
        else if(iop.eq.2) then
          do intf=1,nit
            ni1    = ni2 + 1
            ni2    = inni(intf+1)
            do idof=1,ndof
              do ni=ni1,ni2
                i      = nni(ni)
                v(idof,i) = MAX(v(idof,i),b(idof,ni,get))
              enddo
            enddo
          enddo
        else if(iop.eq.3) then
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
        endif
      else
        if     (iop.eq.1) then
          do intf=1,nit
            ni1    = ni2 + 1
            ni2    = inni(intf+1)
            do idof=1,ndof
              do ni=ni1,ni2
                i      = nni(ni)
                v(i,idof) = MIN(v(i,idof),b(idof,ni,get))
              enddo
            enddo
          enddo
        else if(iop.eq.2) then
          do intf=1,nit
            ni1    = ni2 + 1
            ni2    = inni(intf+1)
            do idof=1,ndof
              do ni=ni1,ni2
                i      = nni(ni)
                v(i,idof) = MAX(v(i,idof),b(idof,ni,get))
              enddo
            enddo
          enddo
        else if(iop.eq.3) then
          do intf=1,nit
            ni1    = ni2 + 1
            ni2    = inni(intf+1)
            do idof=1,ndof
              do ni=ni1,ni2
                i      = nni(ni)
                v(i,idof) = v(i,idof) + b(idof,ni,get)
              enddo
            enddo
          enddo
        endif
      endif      
!
      return
80001 format(i3.3,': ',a,5i14)
      end
