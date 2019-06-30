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
c     "rotpole" constructs the set of atomic multipoles in the global
c     frame by applying the correct rotation matrix for each site
c
c
      subroutine rotpole
      use mpole
      implicit none
      integer i
      real*8 a(3,3)
c     real*8 d1(3,3)
c     real*8 d2(5,5)
c     logical use_qi
c
c
c     rotate the atomic multipoles at each site in turn
c
      do i = 1, npole
         call rotmat (i,a)
         call rotsite (i,a)
      end do
c
c     set QI dipole rotation matrix by permuting the Cartesian
c     rotation matrix to adhere to spherical harmonic ordering
c
c     do i = 1, npole
c        d1(1,1) = a(3,3)
c        d1(2,1) = a(3,1)
c        d1(3,1) = a(3,2)
c        d1(1,2) = a(1,3)
c        d1(2,2) = a(1,1)
c        d1(3,2) = a(1,2)
c        d1(1,3) = a(2,3)
c        d1(2,3) = a(2,1)
c        d1(3,3) = a(2,2)
c        call shrotmat (d1,d2)
c        call shrotsite (i,d1,d2)
c     end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rotmat  --  find global frame rotation matrix  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotmat" finds the rotation matrix that rotates the local
c     coordinate system into the global frame at a multipole site
c
c
      subroutine rotmat (i,a)
      use atoms
      use mpole
      implicit none
      integer i,ii
      integer ix,iy,iz
      real*8 r,dot
      real*8 xi,yi,zi
      real*8 dx,dy,dz
      real*8 dx1,dy1,dz1
      real*8 dx2,dy2,dz2
      real*8 dx3,dy3,dz3
      real*8 a(3,3)
c
c
c     get coordinates and frame definition for the multipole site
c
      ii = ipole(i)
      xi = x(ii)
      yi = y(ii)
      zi = z(ii)
      iz = zaxis(i)
      ix = xaxis(i)
      iy = abs(yaxis(i))
c
c     use the identity matrix as the default rotation matrix
c
      a(1,1) = 1.0d0
      a(2,1) = 0.0d0
      a(3,1) = 0.0d0
      a(1,3) = 0.0d0
      a(2,3) = 0.0d0
      a(3,3) = 1.0d0
c
c     z-only frame rotation matrix elements for z-axis only
c
      if (polaxe(i) .eq. 'Z-Only') then
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
         if (abs(dot) .gt. 0.866d0) then
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
c     z-then-x frame rotation matrix elements for z- and x-axes
c
      else if (polaxe(i) .eq. 'Z-then-X') then
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
c     bisector frame rotation matrix elements for z- and x-axes
c
      else if (polaxe(i) .eq. 'Bisector') then
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
c     z-bisect frame rotation matrix elements for z- and x-axes
c
      else if (polaxe(i) .eq. 'Z-Bisect') then
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
         dx = dx - dot*a(1,3)
         dy = dy - dot*a(2,3)
         dz = dz - dot*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     3-fold frame rotation matrix elements for z- and x-axes
c
      else if (polaxe(i) .eq. '3-Fold') then
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
c     ##  subroutine rotsite  --  rotate multipoles at single site  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotsite" rotates the local frame atomic multipoles at a
c     specified site into the global coordinate frame by applying
c     a rotation matrix
c
c
      subroutine rotsite (isite,a)
      use atoms
      use mpole
      implicit none
      integer i,j,k,m
      integer isite
      real*8 a(3,3)
      real*8 mp(3,3)
      real*8 rp(3,3)
c
c
c     monopoles have the same value in any coordinate frame
c
      rpole(1,isite) = pole(1,isite)
c
c     rotate the dipoles to the global coordinate frame
c
      do i = 2, 4
         rpole(i,isite) = 0.0d0
         do j = 2, 4
            rpole(i,isite) = rpole(i,isite) + pole(j,isite)*a(i-1,j-1)
         end do
      end do
c
c     rotate the quadrupoles to the global coordinate frame
c
      k = 5
      do i = 1, 3
         do j = 1, 3
            mp(i,j) = pole(k,isite)
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
            rpole(k,isite) = rp(i,j)
            k = k + 1
         end do
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine shrotmat  --  QI quadrupole rotation matrix  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "shrotmat" finds the rotation matrix that converts spherical
c     harmonic quadrupoles from the local to the global frame given
c     the required dipole rotation matrix
c
c
      subroutine shrotmat (d1,d2)
      use math
      implicit none
      real*8 d1(3,3)
      real*8 d2(5,5)
c
c
c     build the quadrupole rotation matrix
c
      d2(1,1) = 0.5d0*(3.0d0*d1(1,1)*d1(1,1) - 1.0d0)
      d2(2,1) = sqrtthree*d1(1,1)*d1(2,1)
      d2(3,1) = sqrtthree*d1(1,1)*d1(3,1)
      d2(4,1) = 0.5d0*sqrtthree*(d1(2,1)*d1(2,1) - d1(3,1)*d1(3,1))
      d2(5,1) = sqrtthree*d1(2,1)*d1(3,1)
      d2(1,2) = sqrtthree*d1(1,1)*d1(1,2)
      d2(2,2) = d1(2,1)*d1(1,2) + d1(1,1)*d1(2,2)
      d2(3,2) = d1(3,1)*d1(1,2) + d1(1,1)*d1(3,2)
      d2(4,2) = d1(2,1)*d1(2,2) - d1(3,1)*d1(3,2)
      d2(5,2) = d1(3,1)*d1(2,2) + d1(2,1)*d1(3,2)
      d2(1,3) = sqrtthree*d1(1,1)*d1(1,3)
      d2(2,3) = d1(2,1)*d1(1,3) + d1(1,1)*d1(2,3)
      d2(3,3) = d1(3,1)*d1(1,3) + d1(1,1)*d1(3,3)
      d2(4,3) = d1(2,1)*d1(2,3) - d1(3,1)*d1(3,3)
      d2(5,3) = d1(3,1)*d1(2,3) + d1(2,1)*d1(3,3)
      d2(1,4) = 0.5d0*sqrtthree*(d1(1,2)*d1(1,2) - d1(1,3)*d1(1,3))
      d2(2,4) = d1(1,2)*d1(2,2) - d1(1,3)*d1(2,3)
      d2(3,4) = d1(1,2)*d1(3,2) - d1(1,3)*d1(3,3)
      d2(4,4) = 0.5d0*(d1(2,2)*d1(2,2) - d1(3,2)*d1(3,2)
     &             - d1(2,3)*d1(2,3) + d1(3,3)*d1(3,3))
      d2(5,4) = d1(2,2)*d1(3,2) - d1(2,3)*d1(3,3)
      d2(1,5) = sqrtthree*d1(1,2)*d1(1,3)
      d2(2,5) = d1(2,2)*d1(1,3) + d1(1,2)*d1(2,3)
      d2(3,5) = d1(3,2)*d1(1,3) + d1(1,2)*d1(3,3)
      d2(4,5) = d1(2,2)*d1(2,3) - d1(3,2)*d1(3,3)
      d2(5,5) = d1(3,2)*d1(2,3) + d1(2,2)*d1(3,3)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine shrotsite  --  rotate SH mpoles local to global  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "shrotsite" converts spherical harmonic multipoles from the
c     local to the global frame given required rotation matrices
c
c
      subroutine shrotsite (i,d1,d2)
      use mpole
      implicit none
      integer i
      real*8 d1(3,3)
      real*8 d2(5,5)
c
c
c     find the global frame spherical harmonic multipoles
c
      srpole(1,i) = spole(1,i)
      srpole(2,i) = spole(2,i)*d1(1,1) + spole(3,i)*d1(2,1)
     &                 + spole(4,i)*d1(3,1)
      srpole(3,i) = spole(2,i)*d1(1,2) + spole(3,i)*d1(2,2)
     &                 + spole(4,i)*d1(3,2)
      srpole(4,i) = spole(2,i)*d1(1,3) + spole(3,i)*d1(2,3)
     &                 + spole(4,i)*d1(3,3)
      srpole(5,i) = spole(5,i)*d2(1,1) + spole(6,i)*d2(2,1)
     &                 + spole(7,i)*d2(3,1) + spole(8,i)*d2(4,1)
     &                 + spole(9,i)*d2(5,1)
      srpole(6,i) = spole(5,i)*d2(1,2) + spole(6,i)*d2(2,2)
     &                 + spole(7,i)*d2(3,2) + spole(8,i)*d2(4,2)
     &                 + spole(9,i)*d2(5,2)
      srpole(7,i) = spole(5,i)*d2(1,3) + spole(6,i)*d2(2,3)
     &                 + spole(7,i)*d2(3,3) + spole(8,i)*d2(4,3)
     &                 + spole(9,i)*d2(5,3)
      srpole(8,i) = spole(5,i)*d2(1,4) + spole(6,i)*d2(2,4)
     &                 + spole(7,i)*d2(3,4) + spole(8,i)*d2(4,4)
     &                 + spole(9,i)*d2(5,4)
      srpole(9,i) = spole(5,i)*d2(1,5) + spole(6,i)*d2(2,5)
     &                 + spole(7,i)*d2(3,5) + spole(8,i)*d2(4,5)
     &                 + spole(9,i)*d2(5,5)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine qirotmat  --  QI interatomic rotation matrix  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "qirotmat" finds a rotation matrix that describes the
c     interatomic vector
c
c
      subroutine qirotmat (i,k,rinv,d1)
      use atoms
      use mpole
      implicit none
      integer i,k
      real*8 rinv,r,dot
      real*8 dx,dy,dz
      real*8 a(3,3)
      real*8 d1(3,3)
c
c     set the z-axis to be the interatomic vector
c
      dx = x(k) - x(i)
      dy = y(k) - y(i)
      dz = z(k) - z(i)
      a(1,3) = dx * rinv
      a(2,3) = dy * rinv
      a(3,3) = dz * rinv
c
c     find an x-axis that is orthogonal to the z-axis
c
      if (y(i).ne.y(k) .or. z(i).ne.z(k)) then
        dx = dx + 1.0d0 
      else
        dy = dy + 1.0d0 
      end if
      dot = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
      dx = dx - dot*a(1,3)
      dy = dy - dot*a(2,3)
      dz = dz - dot*a(3,3)
      r = 1.d0 / sqrt(dx*dx + dy*dy + dz*dz)
      a(1,1) = dx * r
      a(2,1) = dy * r
      a(3,1) = dz * r
c
c     get rotation matrix elements for the y-axis
c
      a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
      a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
c
c     reorder to account for spherical harmonics
c
      d1(1,1) = a(3,3)
      d1(2,1) = a(1,3)
      d1(3,1) = a(2,3)
      d1(1,2) = a(3,1)
      d1(2,2) = a(1,1)
      d1(3,2) = a(2,1)
      d1(1,3) = a(3,2)
      d1(2,3) = a(1,2)
      d1(3,3) = a(2,2)
      return
      end
