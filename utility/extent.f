c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2007  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program extent  --  find maximum nonperiodic distance  ##
c     ##                                                         ##
c     #############################################################
c
c
      program extent
      use sizes
      use atoms
      use iounit
      implicit none
      integer i,k
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 r2,rmax
c
c
c     get the structure to be analyzed for its extent
c
      call initial
      call getxyz
c
c     find the largest pairwise distance in the system
c
      rmax = 0.0d0
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do k = i+1, n
            xk = x(k)
            yk = y(k)
            zk = z(k)
            r2 = (xi-xk)**2 + (yi-yk)**2 + (zi-zk)**2
            rmax = max(r2,rmax)
         end do
      end do
      rmax = sqrt(rmax)
c
c     print out the maximum distance in the system
c
      write (iout,10)  rmax
   10 format (/,' Maximum Distance :  ',f12.4)
c
c     perform any final tasks before program exit
c
      call final
      end
