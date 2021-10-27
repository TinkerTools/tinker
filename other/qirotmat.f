c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2021  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotpole  --  rotate multipoles to the QI frame  ##
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
      real*8 d1(3,3)
      real*8 d2(5,5)
      logical use_qi
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
      do i = 1, npole
         d1(1,1) = a(3,3)
         d1(2,1) = a(3,1)
         d1(3,1) = a(3,2)
         d1(1,2) = a(1,3)
         d1(2,2) = a(1,1)
         d1(3,2) = a(1,2)
         d1(1,3) = a(2,3)
         d1(2,3) = a(2,1)
         d1(3,3) = a(2,2)
         call shrotmat (d1,d2)
         call shrotsite (i,d1,d2)
      end do
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
c     literature reference:
c
c     A. C. Simmonett, F. C. Pickard IV, H. F. Schaefer III and
c     B. R. Brooks, "An Efficient Algorithm for Multipole Energies
c     and Derivatives Based on Spherical Harmonics and Extensions
c     to Particle Mesh Ewald", Journal of Chemical Physics, 140,
c     184101 (2014)
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
      d2(2,1) = root3*d1(1,1)*d1(2,1)
      d2(3,1) = root3*d1(1,1)*d1(3,1)
      d2(4,1) = 0.5d0*root3*(d1(2,1)*d1(2,1) - d1(3,1)*d1(3,1))
      d2(5,1) = root3*d1(2,1)*d1(3,1)
      d2(1,2) = root3*d1(1,1)*d1(1,2)
      d2(2,2) = d1(2,1)*d1(1,2) + d1(1,1)*d1(2,2)
      d2(3,2) = d1(3,1)*d1(1,2) + d1(1,1)*d1(3,2)
      d2(4,2) = d1(2,1)*d1(2,2) - d1(3,1)*d1(3,2)
      d2(5,2) = d1(3,1)*d1(2,2) + d1(2,1)*d1(3,2)
      d2(1,3) = root3*d1(1,1)*d1(1,3)
      d2(2,3) = d1(2,1)*d1(1,3) + d1(1,1)*d1(2,3)
      d2(3,3) = d1(3,1)*d1(1,3) + d1(1,1)*d1(3,3)
      d2(4,3) = d1(2,1)*d1(2,3) - d1(3,1)*d1(3,3)
      d2(5,3) = d1(3,1)*d1(2,3) + d1(2,1)*d1(3,3)
      d2(1,4) = 0.5d0*root3*(d1(1,2)*d1(1,2) - d1(1,3)*d1(1,3))
      d2(2,4) = d1(1,2)*d1(2,2) - d1(1,3)*d1(2,3)
      d2(3,4) = d1(1,2)*d1(3,2) - d1(1,3)*d1(3,3)
      d2(4,4) = 0.5d0*(d1(2,2)*d1(2,2) - d1(3,2)*d1(3,2)
     &             - d1(2,3)*d1(2,3) + d1(3,3)*d1(3,3))
      d2(5,4) = d1(2,2)*d1(3,2) - d1(2,3)*d1(3,3)
      d2(1,5) = root3*d1(1,2)*d1(1,3)
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
