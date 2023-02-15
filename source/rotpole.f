c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotpole  --  rotate multipoles to global frame  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rotpole" constructs the global atomic multipoles by applying
c     a rotation matrix to convert from local to global frame
c
c
      subroutine rotpole (poltype)
      use atoms
      use mpole
      use repel
      implicit none
      integer i
      real*8 a(3,3)
      logical planar
      character*5 poltype
c
c
c     rotate local multipoles to global frame at each site
c
      call upcase (poltype)
      if (poltype .eq. 'MPOLE') then
         do i = 1, n
            if (pollist(i) .ne. 0) then
               call rotmat (i,a,planar)
               call rotsite (i,a,planar,pole,rpole)
            end if
         end do
      else if (poltype .eq. 'REPEL') then
         do i = 1, n
            if (replist(i) .ne. 0) then
               call rotmat (i,a,planar)
               call rotsite (i,a,planar,repole,rrepole)
            end if
         end do
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotrpole  --  rotate multipoles to local frame  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rotrpole" constructs the local atomic multipoles by applying
c     a rotation matrix to convert from global to local frame
c
c
      subroutine rotrpole (poltype)
      use atoms
      use mpole
      use repel
      implicit none
      integer i
      real*8 a(3,3)
      logical planar
      character*5 poltype
c
c
c     rotate global multipoles to local frame at each site
c
      call upcase (poltype)
      if (poltype .eq. 'MPOLE') then
         do i = 1, n
            if (pollist(i) .ne. 0) then
               call rotmat (i,a,planar)
               call invert (3,a)
               planar = .false.
               call rotsite (i,a,planar,rpole,pole)
            end if
         end do
      else if (poltype .eq. 'REPEL') then
         do i = 1, n
            if (replist(i) .ne. 0) then
               call rotmat (i,a,planar)
               call invert (3,a)
               planar = .false.
               call rotsite (i,a,planar,rrepole,repole)
            end if
         end do
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine rotmat  --  local-to-global rotation matrix  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "rotmat" finds the rotation matrix that rotates the local
c     coordinate system into the global frame at a specified atom
c
c
      subroutine rotmat (i,a,planar)
      use atoms
      use math
      use mpole
      implicit none
      integer i,ix,iy,iz
      real*8 r,dot
      real*8 eps,angle
      real*8 xi,yi,zi
      real*8 dx,dy,dz
      real*8 dx1,dy1,dz1
      real*8 dx2,dy2,dz2
      real*8 dx3,dy3,dz3
      real*8 dx4,dy4,dz4
      real*8 a(3,3)
      logical planar
      character*8 axetyp
c
c
c     get coordinates and frame definition for multipole site
c
      xi = x(i)
      yi = y(i)
      zi = z(i)
      iz = zaxis(i)
      ix = xaxis(i)
      iy = abs(yaxis(i))
      axetyp = polaxe(i)
      planar = .false.
c
c     use the identity matrix as the default rotation matrix
c
      a(1,1) = 1.0d0
      a(2,1) = 0.0d0
      a(3,1) = 0.0d0
      a(1,2) = 0.0d0
      a(2,2) = 1.0d0
      a(3,2) = 0.0d0
      a(1,3) = 0.0d0
      a(2,3) = 0.0d0
      a(3,3) = 1.0d0
c
c     get Z-Only rotation matrix elements for z-axis only
c
      if (axetyp .eq. 'Z-Only') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = 1.0d0
         dy = 0.0d0
         dz = 0.0d0
         dot = a(1,3)
         eps = 0.707d0
         if (abs(dot) .gt. eps) then
            dx = 0.0d0
            dy = 1.0d0
            dot = a(2,3)
         end if
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     get Z-then-X rotation matrix elements for z- and x-axes
c
      else if (axetyp .eq. 'Z-then-X') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     get Bisector rotation matrix elements for z- and x-axes
c
      else if (axetyp .eq. 'Bisector') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dot = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dot*a(1,3)
         dy = dy2 - dot*a(2,3)
         dz = dz2 - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     get Z-Bisect rotation matrix elements for z- and x-axes;
c     use alternate x-axis if central atom is close to planar
c
      else if (axetyp .eq. 'Z-Bisect') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(iy) - xi
         dy = y(iy) - yi
         dz = z(iy) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx = dx / r
         dy = dy / r
         dz = dz / r
         dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         angle = 180.0d0 - radian*acos(dot)
c        eps = 15.0d0
         eps = 0.0d0
         if (angle .lt. eps) then
            planar = .true.
            dx = dy1*dz2 - dz1*dy2
            dy = dz1*dx2 - dx1*dz2
            dz = dx1*dy2 - dy1*dx2
            dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
            if (dot .lt. 0.0d0) then
               dx = -dx
               dy = -dy
               dz = -dz
               dot = -dot
            end if
         end if
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     get 3-Fold rotation matrix elements for z- and x-axes;
c     use alternate z-axis if central atom is close to planar
c
      else if (axetyp .eq. '3-Fold') then
         dx = x(iz) - xi
         dy = y(iz) - yi
         dz = z(iz) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(ix) - xi
         dy = y(ix) - yi
         dz = z(ix) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = x(iy) - xi
         dy = y(iy) - yi
         dz = z(iy) - zi
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx3 = dx / r
         dy3 = dy / r
         dz3 = dz / r
         dx = dx1 + dx2 + dx3
         dy = dy1 + dy2 + dy3
         dz = dz1 + dz2 + dz3
         r = sqrt(dx*dx + dy*dy + dz*dz)
c        eps = 0.15d0
         eps = 0.0d0
         if (r .lt. eps) then
            planar = .true.
            dx2 = x(ix) - x(iz)
            dy2 = y(ix) - y(iz)
            dz2 = z(ix) - z(iz)
            dx3 = x(iy) - x(iz)
            dy3 = y(iy) - y(iz)
            dz3 = z(iy) - z(iz)
            dx4 = dy2*dz3 - dz2*dy3
            dy4 = dz2*dx3 - dx2*dz3
            dz4 = dx2*dy3 - dy2*dx3
            dot = dx4*dx + dy4*dy + dz4*dz
            if (dot .gt. 0.0d0) then
               dx = dx4
               dy = dy4
               dz = dz4
            else
               dx = -dx4
               dy = -dy4
               dz = -dz4
            end if
            r = sqrt(dx*dx + dy*dy + dz*dz)
         end if
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dot = dx1*a(1,3) + dy1*a(2,3) + dz1*a(3,3)
         dx = dx1 - dot*a(1,3)
         dy = dy1 - dot*a(2,3)
         dz = dz1 - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
      end if
c
c     finally, find rotation matrix elements for the y-axis
c
      a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
      a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rotsite  --  rotate input multipoles to final  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotsite" rotates atomic multipoles from the input to final
c     frame at a specified atom by applying a rotation matrix
c
c
      subroutine rotsite (ii,a,planar,inpole,outpole)
      use mpole
      implicit none
      integer i,j,k,m,ii
      real*8 spole(maxpole)
      real*8 a(3,3)
      real*8 mp(3,3)
      real*8 rp(3,3)
      real*8 inpole(maxpole,*)
      real*8 outpole(maxpole,*)
      logical planar
      character*8 axetyp
c
c
c     copy input multipoles and modify at planar sites
c
      do i = 1, maxpole
         spole(i) = inpole(i,ii)
      end do
      if (planar) then
         axetyp = polaxe(ii)
         if (axetyp .eq. 'Z-Bisect') then
            spole(2) = 0.0d0
            spole(7) = 0.0d0
            spole(11) = 0.0d0
            spole(5) = 0.5d0 * (spole(5)+spole(9))
            spole(9) = spole(5)
         else if (axetyp .eq. '3-Fold') then
            do i = 2, maxpole
               spole(i) = 0.0d0
            end do
         end if
      end if
c
c     monopoles are the same in any coordinate frame
c
      outpole(1,ii) = spole(1)
c
c     rotate input dipoles to final coordinate frame
c
      do i = 2, 4
         outpole(i,ii) = 0.0d0
         do j = 2, 4
            outpole(i,ii) = outpole(i,ii) + spole(j)*a(i-1,j-1)
         end do
      end do
c
c     rotate input quadrupoles to final coordinate frame
c
      k = 5
      do i = 1, 3
         do j = 1, 3
            mp(i,j) = spole(k)
            rp(i,j) = 0.0d0
            k = k + 1
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            if (j .lt. i) then
               rp(i,j) = rp(j,i)
            else
               do k = 1, 3
                  do m = 1, 3
                     rp(i,j) = rp(i,j) + a(i,k)*a(j,m)*mp(k,m)
                  end do
               end do
            end if
         end do
      end do
      k = 5
      do i = 1, 3
         do j = 1, 3
            outpole(k,ii) = rp(i,j)
            k = k + 1
         end do
      end do
      return
      end
