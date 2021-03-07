c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine bounds  --  check periodic boundary conditions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "bounds" finds the center of mass of each molecule and
c     translates any stray molecules back into the periodic box
c
c
      subroutine bounds
      use atomid
      use atoms
      use boxes
      use math
      use molcul
      implicit none
      integer i,j,k
      integer init,stop
      real*8 weigh,corr
      real*8 xmid,ymid,zmid
      real*8 xfrac,yfrac,zfrac
      real*8 xcom,ycom,zcom
c
c
c     locate the center of mass of each molecule
c
      do i = 1, nmol
         init = imol(1,i)
         stop = imol(2,i)
         xmid = 0.0d0
         ymid = 0.0d0
         zmid = 0.0d0
         do j = init, stop
            k = kmol(j)
            weigh = mass(k)
            xmid = xmid + x(k)*weigh
            ymid = ymid + y(k)*weigh
            zmid = zmid + z(k)*weigh
         end do
         weigh = molmass(i)
         xmid = xmid / weigh
         ymid = ymid / weigh
         zmid = zmid / weigh
c
c     get fractional coordinates of center of mass
c
         if (triclinic) then
            zfrac = zmid / gamma_term
            yfrac = (ymid - zfrac*beta_term) / gamma_sin
            xfrac = xmid - yfrac*gamma_cos - zfrac*beta_cos
         else if (monoclinic) then
            zfrac = zmid / beta_sin
            yfrac = ymid
            xfrac = xmid - zfrac*beta_cos
         else
            zfrac = zmid
            yfrac = ymid
            xfrac = xmid
         end if
c
c     translate center of mass into the periodic box
c
         if (dodecadron) then
            xfrac = xfrac - xbox*nint(xfrac/xbox)
            yfrac = yfrac - ybox*nint(yfrac/ybox)
            zfrac = zfrac - root2*zbox*nint(zfrac/(zbox*root2))
            corr = xbox2 * int(abs(xfrac/xbox)+abs(yfrac/ybox)
     &                               +abs(root2*zfrac/zbox))
            xfrac = xfrac - sign(corr,xfrac)
            yfrac = yfrac - sign(corr,yfrac)
            zfrac = zfrac - sign(corr,zfrac)*root2
         else if (octahedron) then
            xfrac = xfrac - xbox*nint(xfrac/xbox)
            yfrac = yfrac - ybox*nint(yfrac/ybox)
            zfrac = zfrac - zbox*nint(zfrac/zbox)
            corr = box23 * int(abs(xfrac/xbox)+abs(yfrac/ybox)
     &                                +abs(zfrac/zbox))
            xfrac = xfrac - sign(corr,xfrac)
            yfrac = yfrac - sign(corr,yfrac)
            zfrac = zfrac - sign(corr,zfrac)
         else
            xfrac = xfrac - xbox*nint(xfrac/xbox)
            yfrac = yfrac - ybox*nint(yfrac/ybox)
            zfrac = zfrac - zbox*nint(zfrac/zbox)
         end if
c
c     convert translated fractional center of mass to Cartesian
c
         if (triclinic) then
            xcom = xfrac + yfrac*gamma_cos + zfrac*beta_cos
            ycom = yfrac*gamma_sin + zfrac*beta_term
            zcom = zfrac * gamma_term
         else if (monoclinic) then
            xcom = xfrac + zfrac*beta_cos
            ycom = yfrac
            zcom = zfrac * beta_sin
         else
            xcom = xfrac
            ycom = yfrac
            zcom = zfrac
         end if
c
c     translate coordinates via offset from center of mass
c
         do j = init, stop
            k = kmol(j)
            x(k) = x(k) - xmid + xcom
            y(k) = y(k) - ymid + ycom
            z(k) = z(k) - zmid + zcom
         end do
      end do
      return
      end
