      subroutine M3Dcoef(x,n,vmin,ierr,ier,v,b,job)
      implicit   none
c-----------------------------------------------------------------------
c- Calcul des coefficients geometriques par element tetraedrique
c         base (x,y,z) ---> base (a,b,c)
c-----------------------------------------------------------------------
c input :
c -----
c  x(i,j) : i eme coordonnee du j eme sommet de l'element
c  n(j)   : numero global du j eme sommet
c  vmin   : volume minimale admise
c  ierr   : numero de l'unite de sortie pour messages
c output :
c ------
c  ier    : indicateur de sortie
c  v      : 6*(volume de l'element)
c  b(i,j) : v*(derivee en la i eme coordonnee de la j eme fct de base)
c> si job = 1
c  v      : volume de l'element
c  b(i,j) : derivee en la i eme coordonnee de la j eme fct de base
c-----------------------------------------------------------------------
      integer    n(4), ierr, ier, job
      real*8     x(3,4), vmin, v, b(3,4)
c
      real*8     xa, xb, xc, ya, yb, yc, za, zb, zc, dv
      real*8     yzbc, yzca, yzab, zxbc, zxca, zxab, xybc, xyca, xyab
c
      ier    = 0
c
      xa     = x(1,1) - x(1,4)
      xb     = x(1,3) - x(1,4)
      xc     = x(1,2) - x(1,4)
      ya     = x(2,1) - x(2,4)
      yb     = x(2,3) - x(2,4)
      yc     = x(2,2) - x(2,4)
      za     = x(3,1) - x(3,4)
      zb     = x(3,3) - x(3,4)
      zc     = x(3,2) - x(3,4)
c
      yzbc   = yb*zc - yc*zb
      yzca   = yc*za - ya*zc
      yzab   = ya*zb - yb*za
      zxbc   = zb*xc - zc*xb
      zxca   = zc*xa - za*xc
      zxab   = za*xb - zb*xa
      xybc   = xb*yc - xc*yb
      xyca   = xc*ya - xa*yc
      xyab   = xa*yb - xb*ya
      v      = xa*yzbc + ya*zxbc + za*xybc
      if(v.le.6.d0*vmin) then
        ier    = -1
        if(ierr.gt.0) write(ierr,90001) v,n(1),n(2),n(3),n(4)
        return
      endif
c
      b(1,1) = yzbc
      b(2,1) = zxbc
      b(3,1) = xybc
      b(1,2) = yzab
      b(2,2) = zxab
      b(3,2) = xyab
      b(1,3) = yzca
      b(2,3) = zxca
      b(3,3) = xyca
      b(1,4) = - (b(1,1)+b(1,2)+b(1,3))
      b(2,4) = - (b(2,1)+b(2,2)+b(2,3))
      b(3,4) = - (b(3,1)+b(3,2)+b(3,3))
c
      if(job.eq.1) then
        dv     = 1.d0/v
        b(1,1) = b(1,1)*dv
        b(2,1) = b(2,1)*dv
        b(3,1) = b(3,1)*dv
        b(1,2) = b(1,2)*dv
        b(2,2) = b(2,2)*dv
        b(3,2) = b(3,2)*dv
        b(1,3) = b(1,3)*dv
        b(2,3) = b(2,3)*dv
        b(3,3) = b(3,3)*dv
        b(1,4) = b(1,4)*dv
        b(2,4) = b(2,4)*dv
        b(3,4) = b(3,4)*dv
        v      = v/6.d0
      endif
c
      return
90001 format('!!! M3Dcoef : 6*volume =',e10.3,
     &       ' pour l''element de sommets (',i8,3(',',i8),')')
      end
