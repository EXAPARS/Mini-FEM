!-----------------------------------------------------------------------
      subroutine SP_stop(texte,MPI_err)
!-----------------------------------------------------------------------
      implicit   none
      include   'para.h'
      character*(*) texte
!
      integer    ier
      logical    lMPI_out
!
      lMPI_out = .false.
      if(MPI_out.gt.0) then
        INQUIRE(unit=MPI_out,opened=lMPI_out)
      endif
!
      if(MPI_err.eq.0)  then
        if(lMPI_out) then
          write(MPI_out,80001) texte
          close(MPI_out)
        else
          print 80002,iblk,texte
        endif
        call MPI_finalize(ier)
      else
        if(lMPI_out) then
          write(MPI_out,90001) texte
          close(MPI_out)
        else
          print 90002,iblk,texte
        endif
        call MPI_abort(MPI_COMM_WORLD, MPI_err, ier)
      endif
!
      stop
80001 format(/'--> SP_stop : ',a)
80002 format(/i3.3,': --> SP_stop : ',a)
90001 format(/'!!! SP_stop : ',a)
90002 format(/i3.3,': !!! SP_stop : ',a)
      end
