      subroutine VecDot(nnt,ndof,idof,n1,n2,vec1,vec2,
#ifdef MPI
     &                  iwsp,buf,
#ifdef NOGVAR
     &                  nblk,nfout,
#endif
#endif
     &                  dot)
      implicit   none
      integer    nnt,ndof,idof,n1,n2
      real*8     vec1(nnt*ndof,n1),vec2(nnt*ndof,n2)
      real*8     dot(n1,n2)
#ifdef MPI
      integer    iwsp(nnt)
      real*8     buf(n1*n2)
#ifdef NOGVAR
      integer    nblk,nfout
#else
      include   'para.h'
#endif
#endif
c
      integer    i,i1,i2,ndl,j
      real*8     dDOT
c
      real*8     zero,one
      data       zero,one                                    /0.d0,1.d0/
c

	  
#ifdef MPI
      if(nblk.gt.1) then
        do i1=1,n1
          do i2=1,n2
            dot(i1,i2) = zero
          enddo
        enddo
        if(idof.gt.0) then
          do i=1,nnt
          if(iwsp(i).gt.0) then
            do i1=1,n1
            do i2=1,n2
              dot(i1,i2) = dot(i1,i2)
     &                   + dDOT(ndof,vec1(i,i1),nnt,vec2(i,i2),nnt)
            enddo
            enddo
          endif
          enddo
        else
          j      = 1
          do i=1,nnt
          if(iwsp(i).gt.0) then
            do i1=1,n1
            do i2=1,n2
              dot(i1,i2) = dot(i1,i2)
     &                   + dDOT(ndof,vec1(j,i1),1,vec2(j,i2),1)
            enddo
            enddo
          endif
          j      = j + ndof
          enddo
        endif

#ifdef NOGVAR
        call SPGdSUM(dot,n1*n2,buf,nfout)
#else
        call SPGdSUM(dot,n1*n2,buf)
#endif
      else
#endif
        ndl    = nnt*ndof
        call dGEMM('T','N',n1,n2,ndl,one,vec1,ndl,vec2,ndl,zero,dot,n1)
#ifdef MPI
      endif
#endif
c
      return
      end
