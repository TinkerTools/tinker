c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2020  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine qrsolver  --  QR factorization as linear solver  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "qrsolver" uses a QR factorization method to solve the linear
c     system Ax = b, returning "x" in "b"; "A" is the upper triangle
c     of a symmetric matrix including the diagonal stored by rows
c
c     literature reference:
c
c     W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling,
c     "Numerical Recipes: The Art of Scientific Computing, 2nd Edition",
c     Cambridge University Press, 1992, Section 2.10
c
c
      subroutine qrsolver (nvar,a,b)
      use iounit
      implicit none
      integer i,j,k,nvar
      real*8 amax,sigma
      real*8 sum,tau
      real*8 a(*)
      real*8 b(*)
      real*8, allocatable :: c(:)
      real*8, allocatable :: d(:)
      real*8, allocatable :: af(:,:)
      logical singular
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (c(nvar))
      allocate (d(nvar))
      allocate (af(nvar,nvar))
c
c     copy input upper triangle into the full matrix
c
      k = 0
      do i = 1, nvar
         do j = i, nvar
            k = k + 1
            af(j,i) = a(k)
            af(i,j) = af(j,i)
         end do
      end do
c
c     perform QR factorization of the input matrix
c
      singular = .false.
      do k = 1, nvar-1
         amax = 0.0d0
         do i = k, nvar
            amax = max(amax,abs(af(i,k)))
         end do
         if (amax .eq. 0.0d0) then
            singular = .true.
            c(k) = 0.0d0
            d(k) = 0.0d0
         else
            do i = k, nvar
               af(i,k) = af(i,k) / amax
            end do
            sum = 0.0d0
            do i = k, nvar
               sum = sum + af(i,k)**2
            end do
            sigma = sign(sqrt(sum),af(k,k))
            af(k,k) = af(k,k) + sigma
            c(k) = sigma * af(k,k)
            d(k) = -amax * sigma
            do j = k+1, nvar
               sum = 0.0d0
               do i = k, nvar
                  sum = sum + af(i,k)*af(i,j)
               end do
               tau = sum / c(k)
               do i = k, nvar
                  af(i,j) = af(i,j) - tau*af(i,k)
               end do
            end do
         end if
      end do
      d(nvar) = af(nvar,nvar)
      if (d(nvar) .eq. 0.0d0)  singular = .true.
      if (singular) then
         write (iout,10)
   10    format (/,' QRSOLVER  --  Input Matrix Singular during',
     &              ' QR Factorization')
         call fatal
      end if
c
c     use factored matrix to solve the linear equations
c
      do j = 1, nvar-1
         sum = 0.0d0
         do i = j, nvar
            sum = sum + af(i,j)*b(i)
         end do
         tau = sum / c(j)
         do i = j, nvar
            b(i) = b(i) - tau*af(i,j)
         end do
      end do
      b(nvar) = b(nvar) / d(nvar)
      do i = nvar-1, 1, -1
         sum = 0.0d0
         do j = i+1, nvar
            sum = sum + af(i,j)*b(j)
         end do
         b(i) = (b(i)-sum) / d(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (c)
      deallocate (d)
      deallocate (af)
      return
      end
