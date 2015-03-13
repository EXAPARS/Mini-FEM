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
