      subroutine DQMRD4(nnt,jlog,nnd,nlim,ier)
      implicit   none
      integer    nnt,jlog,nnd,nlim,ier
      dimension  jlog(nnt),nlim(*)
c
c.... variables locales
c
      integer    nndx,i,k
c
      ier    = 0
      nndx   = nnd
      nnd    = 0
      do i=1,nnt
        if(jlog(i).ne.0) then
          nnd    = nnd + 1
          if(nnd.le.nndx) nlim(nnd) = i
        endif
      enddo
c
      if(nnd.gt.nndx) ier = -1
c
      return
      end
