c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine kinetic  --  compute kinetic energy components  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "kinetic" computes the total kinetic energy and kinetic energy
c     contributions to the pressure tensor by summing over velocities
c
c
      subroutine kinetic (eksum,ekin)
      use sizes
      use atomid
      use atoms
      use bath
      use group
      use mdstuf
      use moldyn
      use rgddyn
      use units
      use usage
      implicit none
      integer i,j,k
      integer start,stop
      real*8 eksum,weigh
      real*8 term,value
      real*8 xr,yr,zr
      real*8 x2,y2,z2
      real*8 xcm,ycm,zcm
      real*8 ekin(3,3)
      real*8 inert(3,3)
c
c
c     zero out the total kinetic energy and its outer product
c
      eksum = 0.0d0
      do i = 1, 3
         do j = 1, 3
            ekin(j,i) = 0.0d0
         end do
      end do
c
c     get the total kinetic energy and tensor for atomic sites
c
      if (integrate .ne. 'RIGIDBODY') then
         do i = 1, n
            if (use(i)) then
               term = 0.5d0 * mass(i) / convert
               do j = 1, 3
                  do k = 1, 3
                     value = term * v(j,i) * v(k,i)
                     ekin(k,j) = ekin(k,j) + value
                  end do
               end do
            end if
         end do
         eksum = ekin(1,1) + ekin(2,2) + ekin(3,3)
c
c     get the total kinetic energy and tensor for rigid bodies
c
      else
         do i = 1, ngrp
            start = igrp(1,i)
            stop = igrp(2,i)
            xcm = 0.0d0
            ycm = 0.0d0
            zcm = 0.0d0
            do j = start, stop
               k = kgrp(j)
               weigh = mass(k)
               xcm = xcm + x(k)*weigh
               ycm = ycm + y(k)*weigh
               zcm = zcm + z(k)*weigh
            end do
            xcm = xcm / grpmass(i)
            ycm = ycm / grpmass(i)
            zcm = zcm / grpmass(i)
c
c     find the inertial tensor relative to the center of mass
c
            do j = 1, 3
               do k = 1, 3
                  inert(k,j) = 0.0d0
               end do
            end do
            do j = start, stop
               k = kgrp(j)
               xr = x(k) - xcm
               yr = y(k) - ycm
               zr = z(k) - zcm
               x2 = xr * xr
               y2 = yr * yr
               z2 = zr * zr
               weigh = mass(k)
               inert(1,1) = inert(1,1) + weigh*(y2+z2)
               inert(2,1) = inert(2,1) - weigh*xr*yr
               inert(3,1) = inert(3,1) - weigh*xr*zr
               inert(2,2) = inert(2,2) + weigh*(x2+z2)
               inert(3,2) = inert(3,2) - weigh*yr*zr
               inert(3,3) = inert(3,3) + weigh*(x2+y2)
            end do
            inert(1,2) = inert(2,1)
            inert(1,3) = inert(3,1)
            inert(2,3) = inert(3,2)
c
c     increment the kinetic energy due to translational motion
c
            term = 0.5d0 * grpmass(i) / convert
            do j = 1, 3
               do k = 1, 3
                  value = term * vc(j,i) * vc(k,i)
                  ekin(k,j) = ekin(k,j) + value
                  if (j .eq. k)  eksum = eksum + value
               end do
            end do
c
c     increment the kinetic energy due to rotational motion
c
            term = 0.5d0 / convert
            do j = 1, 3
               do k = 1, 3
                  value = term * inert(k,j) * wc(j,i) * wc(k,i)
                  eksum = eksum + value
               end do
            end do
         end do
      end if
c
c     get the kinetic energy for Bussi-Parrinello barostat
c
      if (isobaric .and. barostat.eq.'BUSSI') then
         term = dble(nfree) * gasconst * kelvin * taupres * taupres
         value = 0.5d0 * term * eta * eta
         do j = 1, 3
            ekin(j,j) = ekin(j,j) + value/3.0d0
         end do
         eksum = eksum + value
      end if
      return
      end
