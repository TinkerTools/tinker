c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program euler  --  Euler angle-rotation matrix conversion  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "euler" tests an Euler angle --> rotation matrix --> Euler angle
c     sequence of conversions for consistency; this is a check on the
c     method used by the "roteuler" routine in the main TINKER code
c
c
      program euler
      implicit none
      real*8 radian,pi
      parameter (radian=57.29577951308232088d0)
      parameter (pi=3.141592653589793238d0)
      integer i,k
      real*8 phi,theta,psi,eps
      real*8 cphi,ctheta,cpsi
      real*8 sphi,stheta,spsi
      real*8 a(3,3),b(3,3)
      logical err(3)
      character*240 record
c
c
c     get initial values of the Euler angles
c
      write (*,10)
   10 format (/,' Enter Euler Angle Values :  ',$)
      read (*,20)  record
   20 format (a240)
      read (record,*,err=30,end=30)  phi,theta,psi
   30 continue
c
c     write out the initial Euler angles
c
      write (*,40)  phi,theta,psi
   40 format (/,' Input Angles :',5x,3f12.4)
c
c     convert Euler angles from degrees to radians
c
      phi = phi / radian
      theta = theta / radian
      psi = psi / radian
c
c     construct the rotation matrix from Euler angle values
c
      cphi = cos(phi)
      sphi = sin(phi)
      ctheta = cos(theta)
      stheta = sin(theta)
      cpsi = cos(psi)
      spsi = sin(psi)
      a(1,1) = ctheta * cphi
      a(2,1) = spsi*stheta*cphi - cpsi*sphi
      a(3,1) = cpsi*stheta*cphi + spsi*sphi
      a(1,2) = ctheta * sphi
      a(2,2) = spsi*stheta*sphi + cpsi*cphi
      a(3,2) = cpsi*stheta*sphi - spsi*cphi
      a(1,3) = -stheta
      a(2,3) = ctheta * spsi
      a(3,3) = ctheta * cpsi
c
c     write out the initial rotation matrix
c
      write (*,50)
   50 format (/,' Initial Rotation Matrix :',/)
      do i = 1, 3
         write (*,60)  (a(i,k),k=1,3)
   60    format (20x,3f12.4)
      end do
c
c     set the tolerance for Euler angles and rotation elements
c
      eps = 1.0d-8
c
c     get a trial value of theta from a single rotation element
c
      theta = -asin(min(1.0d0,max(-1.0d0,a(1,3))))
      ctheta = cos(theta)
      stheta = -a(1,3)
c
c     set the phi/psi difference when theta is either 90 or -90
c
      if (abs(ctheta) .le. eps) then
         phi = 0.0d0
         if (abs(a(3,1)) .lt. eps) then
            psi = asin(min(1.0d0,max(-1.0d0,-a(2,1)/a(1,3))))
         else if (abs(a(2,1)) .lt. eps) then
            psi = acos(min(1.0d0,max(-1.0d0,-a(3,1)/a(1,3))))
         else
            psi = atan(a(2,1)/a(3,1))
         end if
c
c     set the phi and psi values for all other theta values
c
      else
         if (abs(a(1,1)) .lt. eps) then
            phi = asin(min(1.0d0,max(-1.0d0,a(1,2)/ctheta)))
         else if (abs(a(1,2)) .lt. eps) then
            phi = acos(min(1.0d0,max(-1.0d0,a(1,1)/ctheta)))
         else
            phi = atan(a(1,2)/a(1,1))
         end if
         if (abs(a(3,3)) .lt. eps) then
            psi = asin(min(1.0d0,max(-1.0d0,a(2,3)/ctheta)))
         else if (abs(a(2,3)) .lt. eps) then
            psi = acos(min(1.0d0,max(-1.0d0,a(3,3)/ctheta)))
         else
            psi = atan(a(2,3)/a(3,3))
         end if
      end if
c
c     find sine and cosine of the trial phi and psi values
c
      cphi = cos(phi)
      sphi = sin(phi)
      cpsi = cos(psi)
      spsi = sin(psi)
c
c     reconstruct the diagonal of the rotation matrix
c
      b(1,1) = ctheta * cphi
      b(2,2) = spsi*stheta*sphi + cpsi*cphi
      b(3,3) = ctheta * cpsi
c
c     compare the correct matrix diagonal to rebuilt diagonal
c
      do i = 1, 3
         err(i) = .false.
         if (abs(a(i,i)-b(i,i)) .gt. eps)  err(i) = .true.
      end do
c
c     alter Euler angles to get correct rotation matrix values
c
      if (err(1) .and. err(2))  phi = phi - sign(pi,phi)
      if (err(1) .and. err(3))  theta = -theta + sign(pi,theta)
      if (err(2) .and. err(3))  psi = psi - sign(pi,psi)
c
c     construct the rotation matrix from Euler angle values
c
      cphi = cos(phi)
      sphi = sin(phi)
      ctheta = cos(theta)
      stheta = sin(theta)
      cpsi = cos(psi)
      spsi = sin(psi)
      b(1,1) = ctheta * cphi
      b(2,1) = spsi*stheta*cphi - cpsi*sphi
      b(3,1) = cpsi*stheta*cphi + spsi*sphi
      b(1,2) = ctheta * sphi
      b(2,2) = spsi*stheta*sphi + cpsi*cphi
      b(3,2) = cpsi*stheta*sphi - spsi*cphi
      b(1,3) = -stheta
      b(2,3) = ctheta * spsi
      b(3,3) = ctheta * cpsi
c
c     write out the final Euler angles
c
      write (*,70)  phi*radian,theta*radian,psi*radian
   70 format (/,' Final Angles :',5x,3f12.4)
c
c     write out the final rotation matrix
c
      write (*,80)
   80 format (/,' Final Rotation Matrix :',/)
      do i = 1, 3
         write (*,90)  (b(i,k),k=1,3)
   90    format (20x,3f12.4)
      end do
      end
