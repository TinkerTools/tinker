c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine bspline  --  get B-spline coefficients  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "bspline" calculates the coefficients for an n-th order
c     B-spline approximation
c
c
      subroutine bspline (x,n,c)
      implicit none
      integer i,k,n
      real*8 x,denom
      real*8 c(n)
c
c
c     initialize the B-spline as the linear case
c
      c(1) = 1.0d0 - x
      c(2) = x
c
c     compute standard B-spline recursion to n-th order
c
      do k = 3, n
         denom = 1.0d0 / dble(k-1)
         c(k) = x * c(k-1) * denom
         do i = 1, k-2
            c(k-i) = ((x+dble(i))*c(k-i-1)
     &                  + (dble(k-i)-x)*c(k-i)) * denom
         end do
         c(1) = (1.0d0-x) * c(1) * denom
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine bspline1  --  get B-spline coeffs & derivs  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "bspline1" calculates the coefficients and derivative
c     coefficients for an n-th order B-spline approximation
c
c
      subroutine bspline1 (x,n,c,d)
      implicit none
      integer i,k,n
      real*8 x,denom
      real*8 c(n),d(n)
c
c
c     initialize the B-spline as the linear case
c
      c(1) = 1.0d0 - x
      c(2) = x
c
c     compute standard B-spline recursion to n-1-th order
c
      do k = 3, n-1
         denom = 1.0d0 / dble(k-1)
         c(k) = x * c(k-1) * denom
         do i = 1, k-2
            c(k-i) = ((x+dble(i))*c(k-i-1)
     &                  + (dble(k-i)-x)*c(k-i)) * denom
         end do
         c(1) = (1.0d0-x) * c(1) * denom
      end do
c
c     get the derivative from n-1-th order coefficients
c
      d(1) = -c(1)
      do i = 2, n-1
         d(i) = c(i-1) - c(i)
      end do
      d(n) = c(n-1)
c
c     use one final recursion to get n-th order coefficients
c
      denom = 1.0d0 / dble(n-1)
      c(n) = x * c(n-1) * denom
      do i = 1, n-2
         c(n-i) = ((x+dble(i))*c(n-i-1)
     &               + (dble(n-i)-x)*c(n-i)) * denom
      end do
      c(1) = (1.0d0-x) * c(1) * denom
      return
      end
