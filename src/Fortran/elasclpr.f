! Recursive version of ELASCLPR for Cilk recursion (LT)
      recursive subroutine DC_ELASCLPR(nspa,nsv,nnt,nbe,nne,coor,       &
     &                     innv,nnv,stoc,firstElem,lastElem,ierr,ier)
!-----------------------------------------------------------------------
!   Calcul de la matrice :
!    stoc = -div(sigma)
!    sigma = lambda Tr(Epsilon) Id + 2 mu Epsilon
!    Epsilon = 1/2 (grad(coor)+grad(coor)t)
!    3 lambda + 2 mu = 0
!    lambda = 1 / volume
!-----------------------------------------------------------------------
      implicit   none
      integer    nspa,nsv,nnt,nbe,nne,innv,nnv,ierr,ier
      integer    firstElem,lastElem
      real*8     coor,stoc
      dimension  nne(nsv,nbe),coor(nspa,nnt),innv(nnt+1),nnv(*)
      dimension  stoc(nspa,nspa,*)
!
      integer    i,j,k,l,m,ni,nj,nk,ki,kj
      real*8     xs(3,4),gs(3,4),vol,gsij
      real*8     xlam,xmu
      real*8     xk,xe,xnu
      real*8     a(3,3)
      logical    lerr(1)
!
      real*8     vzero
      DATA       vzero         /1.d-38/
!
      do l=firstElem,lastElem
        do k=1,nsv
          nk     = nne(k,l)
          xs(1,k) = coor(1,nk)
          xs(2,k) = coor(2,nk)
          xs(3,k) = coor(3,nk)
        enddo
        call M3Dcoef(xs,nne(1,l),vzero,ierr,ier,vol,gs,+1)
        lerr(1) = (ier .ne. 0)
        if(lerr(1)) ier=1
        if(ier.ne.0) return
!
        xmu=1.d0
        xlam=0.25d0*xmu
!
        do i=1,nsv
          ni     = nne(i,l)
          do j=1,nsv
            nj     = nne(j,l)
            do m=innv(ni)+1,innv(ni+1)
              if(nnv(m).eq.nj) goto 21
            enddo
            goto 20
 21         continue
c
            a(1,1)=(xlam+2.d0*xmu)*gs(1,j)*gs(1,i)+
     &             xmu*(gs(2,j)*gs(2,i)+gs(3,j)*gs(3,i))
            a(1,2)=(xlam+xmu)*gs(1,j)*gs(2,i)
            a(1,3)=(xlam+xmu)*gs(1,j)*gs(3,i)
c
            a(2,1)=(xlam+xmu)*gs(2,j)*gs(1,i)
            a(2,2)=(xlam+2.d0*xmu)*gs(2,j)*gs(2,i)+
     &             xmu*(gs(1,j)*gs(1,i)+gs(3,j)*gs(3,i))
            a(2,3)=(xlam+xmu)*gs(2,j)*gs(3,i)
c
            a(3,1)=(xlam+xmu)*gs(3,j)*gs(1,i)
            a(3,2)=(xlam+xmu)*gs(3,j)*gs(2,i)
            a(3,3)=(xlam+2.d0*xmu)*gs(3,j)*gs(3,i)+
     &             xmu*(gs(1,j)*gs(1,i)+gs(2,j)*gs(2,i))
c
            do ki=1,nspa
              do kj=1,nspa
                stoc(ki,kj,m) = stoc(ki,kj,m) + a(ki,kj)
              enddo
            enddo
 20       enddo
        enddo
      enddo
!
      return
      end

! Separation of MPI communications from ELASCLPR (LT)
      subroutine ela_comm_mpi(nspa,nnt,prec,nblk,pbuf,nit,lni,npara,    &
     &                        inni,nni)

      implicit   none
      integer    nspa,nnt,nit,lni,npara,inni,nni,ndof,nblk
      real*8     prec,pbuf
      dimension  npara(nit,3),inni(nit+1),nni(lni)
      dimension  prec(nspa,nspa,nnt),pbuf(*)

      if(nblk .gt. 1) then
       ndof=nspa*nspa
       call SPdOPE(nit,lni,npara,inni,nni,ndof,nnt,ndof,prec,pbuf,'sum')
      endif

      return
      end

! Separation of prec inversion from ELASCLPR (LT)
      recursive subroutine ela_invert_prec(nspa,nnt,innv,nnv,prec,ier,  &
     &                                     markf,i)

      implicit   none
      integer    nspa,nnt,innv,nnv,ier,markf
      real*8     prec
      dimension  markf(nnt,nspa),innv(nnt+1),nnv(*)
      dimension  prec(nspa,nspa,nnt)
!
      integer    i,m,ki,kj
      real*8     a(3,3),work(3)
      integer    ipvt(3),info,lwork,lda
!
      lda = 3
      lwork = 3
!
! Correspond to "esymprec"
      do ki=1,nspa
        if(markf(i,ki).ne.0) then
          do kj=1,nspa
            prec(ki,kj,i) = 0.d0
            prec(kj,ki,i) = 0.d0
          end do
          prec(ki,ki,i) = 1.d0
        end if
      end do
!
      do m=innv(i)+1,innv(i+1)
        if(nnv(m).eq.i) goto 61
      enddo
      goto 60
 61   continue
      do ki=1,nspa
        do kj=1,nspa
          a(ki,kj) = prec(ki,kj,i)
        enddo
      enddo
      call DGETRF(nspa,nspa,a,lda,IPVT,INFO)
      if(info.ne.0) then
        print*,"!!! info[DGETRF]=",info
        ier=info
      end if
      call DGETRI(nspa,a,lda,IPVT,WORK,LWORK,INFO)
      if(info.ne.0) then
        print*,"!!! info[DGETRI]=",info
        ier=info
      end if
      do ki=1,nspa
        do kj=1,nspa
          prec(ki,kj,i) = a(ki,kj)
        enddo
      enddo
!
 60   return
      end
