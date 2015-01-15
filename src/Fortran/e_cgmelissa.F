      SUBROUTINE E_CGMelissa(NIVIMP,NDIM,N,INVN,NVN,mat,STOCK,PREC,
     &                     NND,NLIM,JLOG,B,X,LDW,WORK,ITER,RESIDU,
     &                     Me,NbProcs,Maitre,NbIntf,BVoisin,
     &                     INN_Intf,L_NN_Intf,NN_Intf,Lwsp,wsp,pBUF,
     &                     iout,ierr,icvg,IER)
C-----------------------------------------------------------------------
C CGMelissa  résoud le système linéaire A.X = B  en utilisant la méthode
C itérative du Gradient Conjugué Préconditionné (PCG) sur l'ensemble des
C NDIM composantes de X.
C Rappel : PCG applique le préconditionnement en modifiant  la  métrique
C  du gradient.
C
C A est stockée  Morse  dans STOCK avec comme pointeurs  INVN et NVN.
C Le préconditionneur bloc diagonal est stocké dans PREC, i.e. PREC = D^(-1).
C-----------------------------------------------------------------------
C INPUT :
C -------
C NIVIMP        : I    NIVEAU D INFORMATION
C NDIM          : I    NOMBRE D'INCONNUES PAR POINT
C N             : I    DIMENSION DE LA MATRICE (NOMBRE DE NOEUDS)
C INVN(N)       : I    INDEX POINTEUR MORSE
C NVN(*)        : I    POINTEUR MORSE
c mat           : I    type de stockage CSR
c                      == -1  demi-inferieure
c                      == 0   complete
C STOCK(NDIM,NDIM,*) : R8   MATRICE
C PREC(NDIM,NDIM,*)  : R8   MATRICE DE PRECONDITIONNEMENT
C NND           : I    NOMBRE DE NOEUDS DIRICHLET
C NLIM(NND)     : I    LISTE DES NOEUDS DIRICHLET
C JLOG(N)       : I    CODES AUX NOEUDS
C B(N,NDIM)     : R8   SECOND MEMBRE
C X(N,NDIM)     : R8   SOLUTION INITIALE
C LDW           : I    PREMIERE DIMENSION DU TABLEAU WORK
C                      (= NDIM * N en general)
C WORK(LDW,4)   : R8   VECTEUR DE TRAVAIL POUR :
C                      RESIDU, VECTEUR DE DIRECTION, etc...
C ITER          : I    NOMBRE MAXIMUM D'ITERATIONS A EFFECTUER
C RESIDU        : R8   LA VALEUR A ATTEINDRE POUR LA CONVERGENCE :
C                      norm(b - A*x) / norm(b)
C
C OUTPUT :
C --------
C ITER          : I    NOMBRE D'ITERATIONS EFFECTUEES
C X(N,NDIM)     : R8   SOLUTION
C RESIDU        : R8   VALEUR DE FINALE DE : norm(b - A*x) / norm(b)
C IER           : I    ERREUR
C                      ==  0: Successful !!
C                      ==  1: Nombre d'iterations maxi effectue.
C                      ==  2: Arret avec volumes positifs.
C                      == -1: N <= 0 ou NDIM=/2,3
C                      == -2: LDW  < NDIM*N
C                      == -3: ITER <=  0
C-----------------------------------------------------------------------
C Auteurs : Olivier ROQUES                         Date : 23 Avr 1999
C           Quang Vinh DINH                        Date > novembre 2006
C           Gilbert Roge                           Date > 2008-2012
C-----------------------------------------------------------------------
      IMPLICIT   NONE
C
C.... Variables INPUT/OUTPUT
C
      INTEGER    NIVIMP,NDIM,N,INVN,NVN,mat,NND,NLIM,JLOG
      INTEGER    LDW,ITER,IER
      REAL*8     STOCK,PREC,B,X,WORK,RESIDU
      real*8     DDOT
      DIMENSION  INVN(N+1),NVN(*),NLIM(NND),JLOG(N)
      DIMENSION  STOCK(NDIM,NDIM,*)
      DIMENSION  PREC(NDIM,NDIM,*),B(N,NDIM),X(N,NDIM),WORK(LDW,4)
C
C.... Variables ajoutées pour les interfaces
C
      integer    Me,NbProcs,Maitre,NbIntf,L_NN_Intf,iout,icvg,ierr
      integer    BVoisin(NbIntf,3),INN_Intf(*),NN_Intf(*)
      logical    Lwsp(N)
      real*8     wsp(N),pBUF(*)
C
C.... Variables locales
C
      INTEGER    I,J,K,M,M1,M2,MAXIT,NSX
      REAL*8     SCLR1,ALPHA,BETA,RHO,RHO1,TOL,BNRM2
c
      intrinsic  SQRT,MOD
      real*8     E_pDDOT,res2,res,sca
      integer    Fout,Ferr,Fcvg,Ncvg
      logical    msz,Close_cvg,lcvg
      common /CGmeli1/ Fout,Ferr,Fcvg,Ncvg,Close_cvg
c
      real*8     MPI_Wtime
      real*8     mpiwt0,dmpiwt
      common /CGmeli2/ mpiwt0
C
C.... Aliases sur les vecteurs de travail.
C.... Constantes réelles
C
      integer    R,Z,P,Q
      real*8     DZERO,DONE,ZERMAC
      data       R,Z,P,Q /1,2,3,4/
CGR   data       DZERO,DONE,ZERMAC /0.0D+0,1.0D+0,1.d-17/
C sur Tom02 le 23 juin 2009
C EPSILON= 0.222044604925031308E-15
C TINY= 0.890029543402880553E-305
C HUGE= 0.112355820928894747E+306
C on prend donc zermac=2.d-16
c et TOL ne devra pas etre inferieur a zermac
      data       DZERO,DONE,ZERMAC /0.0D+0,1.0D+0,2.d-16/
C
C.... Test des paramètres en entrée
C
      IER    = 0
      NSX    = N*NDIM
      IF     (N.LE.0 .or. (NDIM.ne.2 .and. NDIM.ne.3)) THEN
        IER    = -1
        goto 999
      ELSEIF (LDW .LT. NSX) THEN
        IER    = -2
        goto 999
      ELSEIF (ITER .LE. 0) THEN
        IER    = -3
        goto 999
      ENDIF
c
      if(iout.le.0) then
        Fout   = 6
      else
        Fout   = iout
      endif
      if(ierr.le.0) then
        Ferr   = 6
      else
        Ferr   = ierr
      endif
      if     (icvg.lt.0) then
        Fcvg   = -icvg
        Close_cvg = .false.
      else if(icvg.gt.0) then
        Fcvg   = icvg
        Close_cvg = .true.
      else
        Fcvg   = 88
        Close_cvg = .true.
      endif
      if     (NIVIMP.eq.0) then
        Ncvg   = ITER
      else if(NIVIMP.eq.1) then
        Ncvg   = 100
      else if(NIVIMP.eq.2) then
        Ncvg   =  50
      else if(NIVIMP.eq.3) then
        Ncvg   =  20
      else if(NIVIMP.eq.4) then
        Ncvg   =  10
      else if(NIVIMP.eq.5) then
        Ncvg   =   1
      endif
      lcvg   = Me.eq.Maitre .and. NIVIMP.gt.0
      MAXIT  = ITER
      TOL    = RESIDU
      if(TOL.LT.ZERMAC) THEN
        print*," TOL, ZERMAC:",TOL,ZERMAC
        print*,
     &" on prend donc TOL=ZERMAC pour eviter les instabilites num."
        TOL=ZERMAC
      end if
c     print*," b=",sqrt(ddot(nsx,b,1,b,1))
c     print*," x=",sqrt(ddot(nsx,x,1,x,1))
C
C.... Calcul norm(b) = norm(b-A.Xd)   (d pour Dirichlet)
C
      msz    = .true. ! mise à zéro de Z en dehors des Dirichlet
      call E_ESSBCD(NDIM,N,NND,NLIM,JLOG,X,WORK(1,Z),msz)
c     print*," z=",sqrt(ddot(nsx,WORK(1,Z),1,WORK(1,Z),1))
      call DCOPY(NSX,B,1,WORK(1,R),1)
      call E_PRODMV (NDIM,N,INVN,NVN,mat,STOCK,WORK(1,Z),WORK(1,Q),
     &             NND,NLIM,JLOG,NbProcs,NbIntf,L_NN_Intf,
     &             BVoisin,INN_Intf,NN_Intf,pBUF)
      SCLR1  = -DONE
      call DAXPY(NSX,SCLR1,WORK(1,Q),1,WORK(1,R),1)
      res2   = E_pDDOT(NDIM,N,WORK(1,R),WORK(1,R),NbProcs,Lwsp)
      res    = SQRT(res2)
      if(lcvg) then
        write(Fcvg,'(a,i1,a,es21.12)')
     & '# ',ndim,'-components : R=',res
        mpiwt0 = MPI_Wtime()
      endif
      BNRM2  = res
      if(BNRM2.lt.ZERMAC) then
        if(lcvg) write(Fout,'(a,e10.3,a)')
     & 'Initial residual <',ZERMAC,' : no iteration needed.'
        goto 99
      endif
C
      ITER   = 0
C
C.... Définition du résidu initial.
C.... R = B
C
      call DCOPY(NSX,B,1,WORK(1,R),1)
      res2   = E_pDDOT(NDIM,N,X,X,NbProcs,Lwsp)
c     print*," x=",sqrt(ddot(nsx,x,1,x,1))
C
C...... R = R - A.X
C
      IF(SQRT(res2).GT.DZERO) THEN
        call E_PRODMV (NDIM,N,INVN,NVN,mat,STOCK,X,WORK(1,Q),
     &               NND,NLIM,JLOG,NbProcs,NbIntf,L_NN_Intf,
     &               BVoisin,INN_Intf,NN_Intf,pBUF)
c     print*," q=",sqrt(ddot(nsx,WORK(1,q),1,WORK(1,q),1))
        SCLR1  = -DONE
        call DAXPY(NSX,SCLR1,WORK(1,Q),1,WORK(1,R),1)
        res2   = E_pDDOT(NDIM,N,WORK(1,R),WORK(1,R),NbProcs,Lwsp)
        res    = SQRT(res2)
        RESIDU = SQRT(res2) / BNRM2
        IF(RESIDU.LT.TOL) GOTO 20
      ENDIF
c     print*," r=",sqrt(ddot(nsx,WORK(1,r),1,WORK(1,r),1))
c
C
C....... Perform Preconditioned Conjugate Gradient iteration.
C
 10   ITER   = ITER + 1
C
C....... Résolution avec la matrice de préconditionnement.
C....... Z = PREC^(-1).R
C
         call E_SOLVEP (NDIM,N,INVN,NVN,PREC,WORK(1,R),
     &                NND,NLIM,JLOG,WORK(1,Z))
C
C....... RHO = (R,Z)
C
         RHO    = E_pDDOT(NDIM,N,WORK(1,R),WORK(1,Z),NbProcs,Lwsp)
C
C....... Calcul le vecteur de direction P.
C
         IF(ITER.GT.1) THEN
C
C......... Z = Z + (RHO / RHO1) * P
C
           BETA   = RHO / RHO1
           call DAXPY(NSX,BETA,WORK(1,P),1,WORK(1,Z),1)
         ENDIF
C
C....... P = Z
C
         call DCOPY(NSX,WORK(1,Z),1,WORK(1,P),1)
C
C....... Q = A.P
C
         call E_PRODMV (NDIM,N,INVN,NVN,mat,STOCK,WORK(1,P),WORK(1,Q),
     &                NND,NLIM,JLOG,NbProcs,NbIntf,L_NN_Intf,
     &                BVoisin,INN_Intf,NN_Intf,pBUF)
C
C....... ALPHA = RHO / (P,Q)
C
         sca    = E_pDDOT(NDIM,N,WORK(1,P),WORK(1,Q),NbProcs,Lwsp)
         ALPHA  = RHO / sca
C
C....... Calcul la solution courante X.
C....... X = X + ALPHA * P
C
         call DAXPY(NSX,ALPHA,WORK(1,P),1,X,1)
C
C....... Calcul le vecteur résidu
C
C....... R = R - ALPHA * Q
C
         SCLR1  = -ALPHA
         call DAXPY(NSX,SCLR1,WORK(1,Q),1,WORK(1,R),1)
C
C....... Calcul la norme du résidu et compare avec la tolérance.
C
         res2   = E_pDDOT(NDIM,N,WORK(1,R),WORK(1,R),NbProcs,Lwsp)
         res    = SQRT(res2)
         RESIDU = res / BNRM2
C
C....... La solution a convergé.
C
         IF(RESIDU.LE.TOL) THEN
           IER    = 0
           GOTO 20
         ENDIF
C
C....... Affichage.
C
         if(lcvg .and. MOD(ITER,Ncvg).eq.0) then
           dmpiwt = MPI_Wtime() - mpiwt0
           write(Fcvg,'(i5,es21.12,1pe10.3)') ITER, res, dmpiwt
         endif
C
C....... Nombre itérations max atteint.
C
         IF(ITER.EQ.MAXIT) THEN
           IER    = 1
           GOTO 20
         ENDIF
C
         RHO1   = RHO
C
C....... Retour itération suivante.
C
      GOTO 10
C
C.... Fin des itérations
C
 20   if(lcvg) then
        if     (ier.eq.0) then
          write(Fout,'(a,e10.3,a,i6,a)')
     & 'Convergence test (<=',TOL,') satisfied with',ITER,' iterations.'
        else if(ier.eq.1) then
          write(Fout,'(a,i6,a,e10.3,a)')
     & 'Max number of iterations (=',MAXIT,
     & ') not sufficient to satisfy convergence test (<=',TOL,
     & ') but continue.'
        endif
        write(Fcvg,*) ' '
      endif
C
 99   if(lcvg .and. Close_cvg) close(Fcvg)
 999  RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine E_PRODMV (NDIM,NBP,INVN,NVN,mat,STOCK,X,Y,NND,NLIM,
     &                     JLOG,nblk,nit,lni,npara,inni,nni,pbuf)
C-----------------------------------------------------------------------
C Calcul de Y = STOCK.X
C-----------------------------------------------------------------------
C INPUT :
C -------
C NDIM          : I    NOMBRE D'INCONNUES PAR POINT
C NBP           : I    NOMBRE DE POINTS DU MAILLAGE
C INVN(NBP+1)   : I    INDEX POINTEUR MORSE
C NVN(*)        : I    POINTEUR MORSE
c mat           : I    -1(demi-inferieure) ou 0(complete)
C STOCK(NDIM,NDIM,*) : R8   MATRICE
C X(NBP,NDIM)   : R8   VECTEUR
C NND           : I    NOMBRE DE NOEUDS DIRICHLET
C NLIM(NND)     : I    LISTE DES NOEUDS DIRICHLET
C JLOG(N)       : I    LISTE DES CODES AUX NOEUDS
C
C OUTPUT :
C --------
C Y(NBP,NDIM)   : R8   VECTEUR
C-----------------------------------------------------------------------
C Auteurs : Olivier ROQUES                         Date : 22 Fev 1999
C           Quang Vinh DINH (version parallèle)    Date > novembre 2006
C-----------------------------------------------------------------------
      IMPLICIT   NONE
C
C.... Paramètres INPUT/OUTPUT
C
      INTEGER    NDIM,NBP,INVN,NVN,mat,NND,NLIM,JLOG
      REAL*8     STOCK,X,Y
      DIMENSION  INVN(NBP+1),NVN(*),NLIM(NND),JLOG(NBP)
      DIMENSION  STOCK(NDIM,NDIM,*),X(NBP,NDIM),Y(NBP,NDIM)
c
      integer    nblk,nit,lni
      integer    npara(nit,3),inni(nit+1),nni(lni)
      real*8     pbuf(ndim,lni,2)
C
C.... Variables locales
C
      INTEGER   I,J,KI,KJ,M,M1,M2
C
      DO KI=1,NDIM
        DO I=1,NBP
          Y(I,KI) = 0.d0
        ENDDO
      ENDDO
c
      if(mat.eq.0) then
c
c$OMP PARALLEL DO DEFAULT(NONE)
c$OMP& SHARED(nbp,ndim,invn,nvn,stock,x,y)
c$OMP& PRIVATE(i,j,ki,kj,m,m1,m2)
      DO KI=1,NDIM
      DO KJ=1,NDIM
      DO I=1,NBP
        M1     = INVN(I) + 1
        M2     = INVN(I + 1)
        DO M=M1,M2
          J      = NVN(M)
          Y(I,KI) = Y(I,KI) + STOCK(KI,KJ,M) * X(J,KJ)
        ENDDO
      ENDDO
      ENDDO
      ENDDO
c$OMP END PARALLEL DO
c
      else
c
#ifdef OPENMP
c$OMP PARALLEL DO DEFAULT(NONE)
c$OMP& SHARED(nbp,ndim,invn,nvn,stock,x,y)
c$OMP& PRIVATE(i,j,ki,kj,m,m1,m2)
      DO KI=1,NDIM
      DO KJ=1,NDIM
      DO I=1,NBP
        M1     = INVN(I) + 1
        M2     = INVN(I + 1)
        Y(I,KI) = Y(I,KI)+STOCK(KI,KJ,M2) * X(I,KJ)
        DO M=M1,M2-1
          J      = NVN(M)
          Y(I,KI) = Y(I,KI) + STOCK(KI,KJ,M) * X(J,KJ)
        ENDDO
      ENDDO
      ENDDO
      ENDDO
c$OMP END PARALLEL DO
      DO KI=1,NDIM
      DO KJ=1,NDIM
      DO I=1,NBP
        M1     = INVN(I) + 1
        M2     = INVN(I + 1)
        DO M=M1,M2-1
          J      = NVN(M)
          Y(J,KJ) = Y(J,KJ) + STOCK(KI,KJ,M) * X(I,KI)
        ENDDO
      ENDDO
      ENDDO
      ENDDO
#else
c
      DO KI=1,NDIM
      DO KJ=1,NDIM
      DO I=1,NBP
        M1     = INVN(I) + 1
        M2     = INVN(I + 1)
        Y(I,KI) = Y(I,KI)+STOCK(KI,KJ,M2) * X(I,KJ)
        DO M=M1,M2-1
          J      = NVN(M)
          Y(I,KI) = Y(I,KI) + STOCK(KI,KJ,M) * X(J,KJ)
          Y(J,KJ) = Y(J,KJ) + STOCK(KI,KJ,M) * X(I,KI)
        ENDDO
      ENDDO
      ENDDO
      ENDDO
c
#endif
      endif
c
      if(nblk.gt.1) then
        call SPdOPE(nit,lni,npara,inni,nni,NDIM,NBP,NBP,Y,pbuf,'sum')
      endif
c
      call E_ESSBCH(NDIM,NBP,NND,NLIM,JLOG,Y)
C
C.... Fin
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE E_SOLVEP (NDIM,NBP,INVN,NVN,PREC,B,
     &                   NND,NLIM,JLOG,X)
C-----------------------------------------------------------------------
C Calcul de X = P^(-1).B où P est la matrice de préconditionnement bloc diag.
C    PREC = P^(-1) = D^(-1)
C-----------------------------------------------------------------------
C INPUT :
C -------
C NDIM          : I    NOMBRE D'INCONNUES PAR POINT
C NBP           : I    NOMBRE DE POINTS DU MAILLAGE
C INVN(NBP+1)   : I    INDEX POINTEUR MORSE
C NVN(*)        : I    POINTEUR MORSE
C PREC(NDIM,NDIM,*)  : R8   MATRICE DE PRECONDITIONNEMENT
C B(NBP,NDIM)   : R8   SECOND MEMBRE
C NND           : I    NOMBRE DE NOEUDS DIRICHLET
C NLIM(NND)     : I    LISTE DES NOEUDS DIRICHLET
C JLOG(N)       : I    LISTE DES CODES AUX NOEUDS
C
C OUTPUT :
C --------
C X(NBP,NDIM)   : R8   SOLUTION
C-----------------------------------------------------------------------
C Auteurs : Olivier ROQUES                         Date : 22 Fev 1999
C           Quang Vinh DINH                        Date > novembre 2006
C-----------------------------------------------------------------------
C
      IMPLICIT   NONE
C
C.... Paramètres INPUT/OUTPUT
C
      INTEGER    NDIM,NBP,INVN,NVN,NND,NLIM,JLOG
      REAL*8     PREC,X,B
      DIMENSION  INVN(NBP+1),NVN(*),NLIM(NND),JLOG(NBP)
      DIMENSION PREC(NDIM,NDIM,*),X(NBP,NDIM),B(NBP,NDIM)
C
C.... Variables locales
C
      INTEGER    I,KI,KJ
C
C.... Initialisation du second membre pour les noeuds Dirichlet
C
      call E_ESSBCH(NDIM,NBP,NND,NLIM,JLOG,B)
C
C.... Préconditionnement diagonal
C
      DO KI=1,NDIM
        DO I=1,NBP
          X(I,KI) = 0.d0
        ENDDO
      ENDDO
C
      DO KI=1,NDIM
      DO KJ=1,NDIM
      DO I=1,NBP
        X(I,KI) = X(I,KI)+PREC(KI,KJ,I) * B(I,KJ)
      ENDDO
      ENDDO
      ENDDO
C
C.... Fin
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine E_ESSBCD(NDIM,NBP,NND,NLIM,JLOG,D,X,msz)
      implicit   none
      include   'DefMesh_jlog.h'
C
C.... Imposition des C.L. non-homogènes aux noeuds Dirichlet
C
      INTEGER    NDIM,NBP,NND,NLIM,JLOG
      REAL*8     D,X
      logical    msz
      DIMENSION  NLIM(NND),JLOG(NBP)
      DIMENSION  D(NBP,NDIM),X(NBP,NDIM)
C
      integer    i,j,nl
C
      if(msz) then
        do j=1,ndim
        do i=1,nbp
          x(i,j) = 0.d0
        enddo
        enddo
      endif
C
      if     (NDIM.EQ.2) then
        DO I=1,NND
          NL = NLIM(I)
          J  = JLOG(NL)
          IF     (J.EQ.jlsymx) THEN
            X(NL,1) = D(NL,1)
          ELSEIF (J.EQ.jlsymz) THEN
            X(NL,2) = D(NL,2)
          ELSE
            X(NL,1) = D(NL,1)
            X(NL,2) = D(NL,2)
          ENDIF
        ENDDO
      else if(NDIM.EQ.3) then
        DO I=1,NND
          NL = NLIM(I)
          J  = JLOG(NL)
          IF     (J.EQ.jlsymx) THEN
            X(NL,1) = D(NL,1)
          ELSEIF (J.EQ.jlsymy) THEN
            X(NL,2) = D(NL,2)
          ELSEIF (J.EQ.jlsymz) THEN
            X(NL,3) = D(NL,3)
          ELSE
            X(NL,1) = D(NL,1)
            X(NL,2) = D(NL,2)
            X(NL,3) = D(NL,3)
          ENDIF
        ENDDO
      endif
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine E_ESSBCH(NDIM,NBP,NND,NLIM,JLOG,X)
      implicit   none
      include   'DefMesh_jlog.h'
C
C.... Imposition des C.L. homogènes aux noeuds Dirichlet
C
      INTEGER    NDIM,NBP,NND,NLIM,JLOG
      REAL*8     X
      DIMENSION  NLIM(NND),JLOG(NBP)
      DIMENSION  X(NBP,NDIM)
C
      integer    i,j,nl
C
      if     (NDIM.EQ.2) then
        do I=1,NND
          NL = NLIM(I)
          J  = JLOG(NL)
          IF     (J.EQ.jlsymx) THEN
            X(NL,1) = 0.D0
          ELSEIF (J.EQ.jlsymz) THEN
            X(NL,2) = 0.D0
          ELSE
            X(NL,1) = 0.D0
            X(NL,2) = 0.D0
          ENDIF
        ENDDO
      else if(NDIM.EQ.3) then
        do I=1,NND
          NL = NLIM(I)
          J  = JLOG(NL)
          IF     (J.EQ.jlsymx) THEN
            X(NL,1) = 0.D0
          ELSEIF (J.EQ.jlsymy) THEN
            X(NL,2) = 0.D0
          ELSEIF (J.EQ.jlsymz) THEN
            X(NL,3) = 0.D0
          ELSE
            X(NL,1) = 0.D0
            X(NL,2) = 0.D0
            X(NL,3) = 0.D0
          ENDIF
        ENDDO
      endif
C
      RETURN
      END
C-----------------------------------------------------------------------
      subroutine E_ESSBCM(NDIM,NBP,NND,NLIM,JLOG,MARKF)
      implicit   none
      include   'DefMesh_jlog.h'
C
C.... Construction du vecteur cache pour les C.L. aux noeuds Dirichlet
C
      INTEGER    NDIM,NBP,NND,NLIM,JLOG,MARKF
      DIMENSION  NLIM(NND),JLOG(NBP)
      DIMENSION  MARKF(NBP,NDIM)
C
      integer    i,j,nl
C
      do j=1,NDIM
      do i=1,NBP
        MARKF(i,j) = 0
      enddo
      enddo
C
      if     (NDIM.EQ.2) then
        do I=1,NND
          NL = NLIM(I)
          J  = JLOG(NL)
          IF     (J.EQ.jlsymx) THEN
            MARKF(NL,1) = 1
          ELSEIF (J.EQ.jlsymz) THEN
            MARKF(NL,2) = 1
          ELSE
            MARKF(NL,1) = 1
            MARKF(NL,2) = 1
          ENDIF
        ENDDO
      else if(NDIM.EQ.3) then
        do I=1,NND
          NL = NLIM(I)
          J  = JLOG(NL)
          IF     (J.EQ.jlsymx) THEN
            MARKF(NL,1) = 1
          ELSEIF (J.EQ.jlsymy) THEN
            MARKF(NL,2) = 1
          ELSEIF (J.EQ.jlsymz) THEN
            MARKF(NL,3) = 1
          ELSE
            MARKF(NL,1) = 1
            MARKF(NL,2) = 1
            MARKF(NL,3) = 1
          ENDIF
        ENDDO
      endif
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function E_pDDOT(NDIM,NBP,X,Y,NbProcs,Lwsp)
      implicit  none
      integer   NDIM,NBP,NbProcs
      logical   Lwsp(NBP)
      real*8    X(NBP,NDIM),Y(NBP,NDIM)
c
      integer   i,j
      real*8    DDOT,sca,buf
c
      if(NbProcs.gt.1) then
        sca    = 0.d0
        do j=1,ndim
        do i=1,nbp
        if(Lwsp(i)) then
          sca    = sca + X(i,j)*Y(i,j)
        endif
        enddo
        enddo
        call SPGdSUM(sca,1,buf)
      else
        sca    = DDOT(NBP*NDIM,X,1,Y,1)
      endif
c
      E_pDDOT  = sca
      return
      end
