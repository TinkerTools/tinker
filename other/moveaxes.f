c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program moveaxes  --  switch local axes for multipoles  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "moveaxes" converts atomic multipoles from the original local
c     axis definitions to the bisector method; in the original mode,
c     the "zaxis" atom defines the z-direction and the "xaxis" atom
c     lies in the positive-x half of the z,x-plane; in the alternate
c     "bisector" the z-direction is taken as the bisector of the
c     angle (zaxis atom-central atom-xaxis atom), the x-direction is
c     then perpendicular to the z-axis and in the plane
c
c
      program moveaxes
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'units.i'
      integer i,k,ia
      real*8 a(3,3)
c
c
c     get the molecule and set up the energy parameters
c
      call initial
      call getxyz
      call mechanic
c
c     convert dipole and quadrupole moments back to atomic units
c
      do i = 1, npole
         do k = 2, 4
            pole(k,i) = pole(k,i) / bohr
         end do
         do k = 5, 13
            pole(k,i) = 3.0d0 * pole(k,i) / bohr**2
         end do
      end do
c
c     rotate multipoles to interconvert local coordinate systems
c
      do i = 1, npole
         ia = ipole(i)
         write (iout,10)  i,ia,zaxis(i),xaxis(i),polaxe(i),
     &                    (pole(k,i),k=1,5),pole(8,i),
     &                    pole(9,i),(pole(k,i),k=11,13)
   10    format (/,i6,4x,i6,6x,i6,1x,i6,3x,a8,4x,f9.5,/,50x,3f9.5,
     &              /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
         if (polaxe(i) .eq. 'Z-then-X') then
            polaxe(i) = 'Bisector'
            call rotmatx (i,a)
            call rotsite (i,a)
            write (iout,20)  i,ia,zaxis(i),xaxis(i),polaxe(i),
     &                       (rpole(k,i),k=1,5),rpole(8,i),
     &                       rpole(9,i),(rpole(k,i),k=11,13)
   20       format (/,i6,4x,i6,6x,i6,1x,i6,3x,a8,4x,f9.5,/,50x,3f9.5,
     &                 /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
         end if
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##########################
c     ##                      ##
c     ##  subroutine rotmatx  ##
c     ##                      ##
c     ##########################
c
c
c     "rotmatx" find the rotation matrix that converts from the local
c     coordinate system at each multipole site from the "Z-then-X"
c     system to the "bisector" system
c
c
      subroutine rotmatx (i,a)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i
      real*8 dx,dy,dz,r
      real*8 dxz1,dyz1,dzz1
      real*8 dxz2,dyz2,dzz2
      real*8 dxz,dyz,dzz
      real*8 dxx1,dyx1,dzx1
      real*8 dxx2,dyx2,dzx2
      real*8 dotxz,a(3,3)
c
c
c     z1 vector
c
      dx = x(zaxis(i)) - x(ipole(i))
      dy = y(zaxis(i)) - y(ipole(i))
      dz = z(zaxis(i)) - z(ipole(i))
      r = sqrt(dx*dx + dy*dy + dz*dz)
      dxz1 = dx / r
      dyz1 = dy / r
      dzz1 = dz / r
c
c     z2 vector
c
      dx = x(xaxis(i)) - x(ipole(i))
      dy = y(xaxis(i)) - y(ipole(i))
      dz = z(xaxis(i)) - z(ipole(i))
      r = sqrt(dx*dx + dy*dy + dz*dz)
      dxz2 = dx / r
      dyz2 = dy / r
      dzz2 = dz / r
c
c     bisector of z1 and z2
c
      dx = dxz1 + dxz2
      dy = dyz1 + dyz2
      dz = dzz1 + dzz2
      r = sqrt(dx*dx + dy*dy + dz*dz)
      dxz = dx / r
      dyz = dy / r
      dzz = dz / r
c
c     x vector in (z-x) system
c
      dotxz = dxz2*dxz1 + dyz2*dyz1 + dzz2*dzz1
      dxx1 = dxz2 - dotxz*dxz1
      dyx1 = dyz2 - dotxz*dyz1
      dzx1 = dzz2 - dotxz*dzz1
      r = sqrt(dxx1*dxx1 + dyx1*dyx1 + dzx1*dzx1)
      dxx1 = dxx1 / r
      dyx1 = dyx1 / r
      dzx1 = dzx1 / r
c
c     x vector in (zz-x) system
c
      dotxz = dxz2*dxz + dyz2*dyz + dzz2*dzz
      dxx2 = dxz2 - dotxz*dxz
      dyx2 = dyz2 - dotxz*dyz
      dzx2 = dzz2 - dotxz*dzz
      r = sqrt(dxx2*dxx2 + dyx2*dyx2 + dzx2*dzx2)
      dxx2 = dxx2 / r
      dyx2 = dyx2 / r
      dzx2 = dzx2 / r
      a(1,1) = dxx2*dxx1 + dyx2*dyx1 + dzx2*dzx1
      a(1,2) = 0.0d0
      a(1,3) = dxx2*dxz1 + dyx2*dyz1 + dzx2*dzz1
      a(2,1) = 0.0d0
      a(2,2) = 1.0d0
      a(2,3) = 0.0d0
      a(3,1) = dxz*dxx1 + dyz*dyx1 + dzz*dzx1
      a(3,2) = 0.0d0
      a(3,3) = dxz*dxz1 + dy *dyz1 + dzz*dzz1
      return
      end
