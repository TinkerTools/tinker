c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine numgrad  --  numerical gradient of a function  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "numgrad" computes the gradient of the objective function
c     "evalue" with respect to Cartesian coordinates of the atoms
c     via a one-sided or two-sided numerical differentiation
c
c
      subroutine numgrad (evalue,g,eps)
      use atoms
      implicit none
      integer i
      real*8 evalue,eps
      real*8 e,e0,old
      real*8 g(3,*)
      logical twosided
      external evalue
c
c
c     chose between use of one-sided or two-sided gradient
c
      twosided = .true.
      if (.not. twosided)  e0 = evalue ()
c
c     compute the numerical gradient from function values
c
      do i = 1, n
         old = x(i)
         if (twosided) then
            x(i) = x(i) - 0.5d0*eps
            e0 = evalue ()
         end if
         x(i) = x(i) + eps
         e = evalue ()
         x(i) = old
         g(1,i) = (e - e0) / eps
         old = y(i)
         if (twosided) then
            y(i) = y(i) - 0.5d0*eps
            e0 = evalue ()
         end if
         y(i) = y(i) + eps
         e = evalue ()
         y(i) = old
         g(2,i) = (e - e0) / eps
         old = z(i)
         if (twosided) then
            z(i) = z(i) - 0.5d0*eps
            e0 = evalue ()
         end if
         z(i) = z(i) + eps
         e = evalue ()
         z(i) = old
         g(3,i) = (e - e0) / eps
      end do
      return
      end
