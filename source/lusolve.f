c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2020  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine lusolve  --  LU factorization as linear solver  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "lusolve" uses a LU factorization with partial pivoting to solve
c     the linear system Ax = b, returning "x" in "b"; "A" is the upper
c     triangle of a symmetric matrix and the diagonal stored by rows
c
c     literature reference:
c
c     W. H. Press, B. P. Flannery, S. A. Teukolsky and W. T. Vetterling,
c     "Numerical Recipes: The Art of Scientific Computing, 2nd Edition",
c     Cambridge University Press, 1992, Section 2.3
c
c
      subroutine lusolve (nvar,a,b)
      use iounit
      implicit none
      integer i,j,k,m
      integer nvar,imax
      integer, allocatable :: indx(:)
      real*8 amax,sum
      real*8 eps,temp
      real*8 a(*)
      real*8 b(*)
      real*8, allocatable :: vv(:)
      real*8, allocatable :: af(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (indx(nvar))
      allocate (vv(nvar))
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
c     perform LU factorization of the input matrix
c
      do i = 1, nvar
         amax = 0.0d0
         do j = 1, nvar
            if (abs(af(i,j)) .gt. amax)  amax = abs(af(i,j))
         end do
         if (amax .eq. 0.0d0) then
            write (iout,10)
   10       format (/,' LUSOLVE  --  Input Matrix Singular during',
     &                 ' LU Factorization')
            call fatal
         end if
         vv(i) = 1.0d0 / amax
      end do
      eps = 1.0d-10
      do j = 1, nvar
         do i = 1, j-1
            sum = af(i,j)
            do k = 1, i-1
               sum = sum - af(i,k)*af(k,j)
            end do
            af(i,j) = sum
         end do
         amax = 0.0d0
         do i = j, nvar
            sum = af(i,j)
            do k = 1, j-1
               sum = sum - af(i,k)*af(k,j)
            end do
            af(i,j) = sum
            temp = vv(i) * abs(sum)
            if (temp .ge. amax) then
               imax = i
               amax = temp
            end if
         end do
         if (j .ne. imax) then
            do k = 1, nvar
               temp = af(imax,k)
               af(imax,k) = af(j,k)
               af(j,k) = temp
            end do
            vv(imax) = vv(j)
         end if
         indx(j) = imax
         if (af(j,j) .eq. 0.0d0)  af(j,j) = eps
         if (j .ne. nvar) then
            temp = 1.0d0 / af(j,j)
            do i = j+1, nvar
               af(i,j) = af(i,j) * temp
            end do
         end if
      end do
c
c     use factored matrix to solve the linear equations
c
      m = 0
      do i = 1, nvar
         k = indx(i)
         sum = b(k)
         b(k) = b(i)
         if (m .ne. 0) then
            do j = m, i-1
               sum = sum - af(i,j)*b(j)
            end do
         else if (sum .ne. 0.0d0) then
            m = i
         end if
         b(i) = sum
      end do
      do i = nvar, 1, -1
         sum = b(i)
         do j = i+1, nvar
            sum = sum - af(i,j)*b(j)
         end do
         b(i) = sum / af(i,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (indx)
      deallocate (vv)
      deallocate (af)
      return
      end
