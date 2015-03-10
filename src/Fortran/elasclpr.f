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
