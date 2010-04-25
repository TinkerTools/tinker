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
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'molcul.i'
      integer i,j,k
      integer init,stop
      real*8 weigh
      real*8 xmid,ymid,zmid
      real*8 xfrac,yfrac,zfrac
      real*8 xoff(maxatm)
      real*8 yoff(maxatm)
      real*8 zoff(maxatm)
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
c     save atomic coordinates relative to center of mass
c
         do j = init, stop
            k = kmol(j)
            xoff(k) = x(k) - xmid
            yoff(k) = y(k) - ymid
            zoff(k) = z(k) - zmid
         end do
c
c     get fractional coordinates of center of mass
c
         if (orthogonal .or. octahedron) then
            zfrac = zmid
            yfrac = ymid
            xfrac = xmid
         else if (monoclinic) then
            zfrac = zmid / beta_sin
            yfrac = ymid
            xfrac = xmid - zfrac*beta_cos
         else if (triclinic) then
            zfrac = zmid / gamma_term
            yfrac = (ymid - zfrac*beta_term) / gamma_sin
            xfrac = xmid - yfrac*gamma_cos - zfrac*beta_cos
         end if
c
c     translate center of mass into the periodic box
c
         dowhile (xfrac .gt. xbox2)
            xfrac = xfrac - xbox
         end do
         dowhile (xfrac .lt. -xbox2)
            xfrac = xfrac + xbox
         end do
         dowhile (yfrac .gt. ybox2)
            yfrac = yfrac - ybox
         end do
         dowhile (yfrac .lt. -ybox2)
            yfrac = yfrac + ybox
         end do
         dowhile (zfrac .gt. zbox2)
            zfrac = zfrac - zbox
         end do
         dowhile (zfrac .lt. -zbox2)
            zfrac = zfrac + zbox
         end do
c
c     truncated octahedron needs to have corners removed
c
         if (octahedron) then
            if (abs(xfrac)+abs(yfrac)+abs(zfrac) .gt. box34) then
               xfrac = xfrac - sign(xbox2,xfrac)
               yfrac = yfrac - sign(ybox2,yfrac)
               zfrac = zfrac - sign(zbox2,zfrac)
            end if
         end if
c
c     convert fractional center of mass back to Cartesian
c
         if (orthogonal .or. octahedron) then
            xmid = xfrac
            ymid = yfrac
            zmid = zfrac
         else if (monoclinic) then
            xmid = xfrac + zfrac*beta_cos
            ymid = yfrac
            zmid = zfrac * beta_sin
         else if (triclinic) then
            xmid = xfrac + yfrac*gamma_cos + zfrac*beta_cos
            ymid = yfrac*gamma_sin + zfrac*beta_term
            zmid = zfrac * gamma_term
         end if
c
c     translate coordinates via offset from center of mass
c
         do j = init, stop
            k = kmol(j)
            x(k) = xoff(k) + xmid
            y(k) = yoff(k) + ymid
            z(k) = zoff(k) + zmid
         end do
      end do
      return
      end
