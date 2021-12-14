c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2020  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine simplex  --  Nelder-Mead simplex optimization  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "simplex" is a general multidimensional Nelder-Mead simplex
c     optimization routine requiring only repeated evaluations of
c     the objective function
c
c     literature reference:
c
c     R. O'Neill, "Algorithm AS 47: Function Minimization Using a
c     Simplex Procedure", Applied Statistics, 20, 338-345 (1971)
c
c
      subroutine simplex (nvar,iter,ntest,x0,y0,step,toler,fvalue)
      use inform
      use iounit
      use keys
      use minima
      implicit none
      real*8 ccoeff,ecoeff
      real*8 rcoeff,eps
      parameter (ccoeff=0.5d0)
      parameter (ecoeff=2.0d0)
      parameter (rcoeff=1.0d0)
      parameter (eps=0.001d0)
      integer i,j,k,nvar
      integer iter,next
      integer ihi,ilo
      integer ntest,jtest
      real*8 toler,tol
      real*8 fvalue,step
      real*8 x,z,del
      real*8 y0,ylo
      real*8 ystar,y2star
      real*8 x0(*)
      real*8, allocatable :: xmin(:)
      real*8, allocatable :: pbar(:)
      real*8, allocatable :: pstar(:)
      real*8, allocatable :: p2star(:)
      real*8, allocatable :: y(:)
      real*8, allocatable :: p(:,:)
      logical done
      character*20 keyword
      character*240 record
      character*240 string
      external fvalue
c
c
c     set default parameters for the optimization
c
      if (maxiter .eq. 0)  maxiter = 1000000
      if (iprint .lt. 0)  iprint = 1000
      if (iwrite .lt. 0)  iwrite = 1000
c
c     search the keywords for optimization parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:8) .eq. 'MAXITER ') then
            read (string,*,err=10,end=10)  maxiter
         end if
   10    continue
      end do
c
c     initialization of various counters and variables
c
      done = .false.
      iter = 0
      jtest = ntest
      del = 1.0d0
      tol = toler * dble(nvar)
c
c     print header information about the optimization method
c
      if (iprint .gt. 0) then
         write (iout,20)
   20    format (/,' Nelder-Mead Simplex Optimization :')
         write (iout,30)
   30    format (/,' NM Iter     F Value      G RMS      F Move',
     &              '   X Move    Comment')
         flush (iout)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (xmin(nvar))
      allocate (pbar(nvar))
      allocate (pstar(nvar))
      allocate (p2star(nvar))
      allocate (y(nvar+1))
      allocate (p(nvar,nvar+1))
c
c     initialize or restart with the base function value
c
      dowhile (.not. done)
         do i = 1, nvar
            p(i,nvar+1) = x0(i)
         end do
         y(nvar+1) = fvalue (x0)
         iter = iter + 1
c
c     define the initial simplex as an "nvar+1" polytope
c
         do j = 1, nvar
            x = x0(j)
            x0(j) = x0(j) + step*del
            do i = 1, nvar
               p(i,j) = x0(i)
            end do
            y(j) = fvalue (x0)
            iter = iter + 1
            x0(j) = x
         end do
c
c     find highest and lowest values; highest will be replaced
c
         ilo = nvar + 1
         ylo = y(ilo)
         do i = 1, nvar
            if (y(i) .le. ylo) then
               ilo = i
               ylo = y(ilo)
            end if
         end do
c
c     set "y0" to be the current highest function value
c
         dowhile (iter .lt. maxiter)
            ihi = nvar + 1
            y0 = y(ihi)
            do i = 1, nvar
               if (y(i) .ge. y0) then
                  ihi = i
                  y0 = y(ihi)
               end if
            end do
c
c     calculate "pbar", the centroid of the simplex vertices
c     excepting the vertex with the highest function value
c
            do i = 1, nvar
               pbar(i) = 0.0d0
               do j = 1, nvar+1
                  pbar(i) = pbar(i) + p(i,j)
               end do
               pbar(i) = (pbar(i)-p(i,ihi)) / dble(nvar)
            end do
c
c     reflection through the centroid of the vertices
c
            do i = 1, nvar
               pstar(i) = pbar(i) + rcoeff*(pbar(i)-p(i,ihi))
            end do
            ystar = fvalue (pstar)
            iter = iter + 1
c
c     successful reflection, so try simplex extension
c
            if (ystar .lt. ylo) then
               do i = 1, nvar
                  p2star(i) = pbar(i) + ecoeff*(pstar(i)-pbar(i))
               end do
               y2star = fvalue (p2star)
               iter = iter + 1
c
c     retain extension or contraction of the simplex
c
               if (ystar .lt. y2star) then
                  do i = 1, nvar
                     p(i,ihi) = pstar(i)
                  end do
                  y(ihi) = ystar
               else
                  do i = 1, nvar
                     p(i,ihi) = p2star(i)
                  end do
                  y(ihi) = y2star
               end if
c
c     no extension of the simplex will be used
c
            else
               k = 0
               do i = 1, nvar+1
                  if (ystar .lt. y(i)) then
                     k = k + 1
                  end if
               end do
               if (1 .lt. k) then
                  do i = 1, nvar
                     p(i,ihi) = pstar(i)
                  end do
                  y(ihi) = ystar
c
c     contraction on the "ihi" side of the centroid
c
               else if (k .eq. 0) then
                  do i = 1, nvar
                     p2star(i) = pbar(i) + ccoeff*(p(i,ihi)-pbar(i))
                  end do
                  y2star = fvalue (p2star)
                  iter = iter + 1
c
c     perform contraction of the whole simplex
c
                  if (y(ihi) .lt. y2star) then
                     do j = 1, nvar+1
                        do i = 1, nvar
                           p(i,j) = 0.5d0 * (p(i,j)+p(i,ilo))
                        end do
                        do i = 1, nvar
                           xmin(i) = p(i,j)
                        end do
                        y(j) = fvalue (xmin)
                        iter = iter + 1
                     end do
                     ilo = nvar + 1
                     ylo = y(ilo)
                     do i = 1, nvar
                        if (y(i) .le. ylo) then
                           ilo = i
                           ylo = y(ilo)
                        end if
                     end do
                     goto 40
c
c     retain the contraction of the simplex
c
                  else
                     do i = 1, nvar
                        p(i,ihi) = p2star(i)
                     end do
                     y(ihi) = y2star
                  end if
c
c     contraction on the reflection side of the centroid
c
               else if (k .eq. 1) then
                  do i = 1, nvar
                     p2star(i) = pbar(i) + ccoeff*(pstar(i)-pbar(i))
                  end do
                  y2star = fvalue (p2star)
                  iter = iter + 1
c
c     check whether to retain reflection of the simplex
c
                  if (y2star .le. ystar) then
                     do i = 1, nvar
                        p(i,ihi) = p2star(i)
                     end do
                     y(ihi) = y2star
                  else
                     do i = 1, nvar
                        p(i,ihi) = pstar(i)
                     end do
                     y(ihi) = ystar
                  end if
               end if
            end if
c
c     check to see if the "ylo" value has improved
c
            if (y(ihi) .lt. ylo) then
               ylo = y(ihi)
               ilo = ihi
            end if
c
c     check to see if the desired minimum has been reached
c
            jtest = jtest -1
            if (jtest .eq. 0) then
               if (iter .le. maxiter) then
                  jtest = ntest
                  x = 0.0d0
                  do i = 1, nvar+1
                     x = x + y(i)
                  end do
                  x = x / dble(nvar+1)
                  z = 0.0d0
                  do i = 1, nvar+1
                     z = z + (y(i)-x)**2
                  end do
                  if (z .le. tol) then
                     goto 50
                  end if
               end if
            end if
   40       continue
         end do
   50    continue
c
c     factorial tests to check if "y0" is a local minimum
c
         do i = 1, nvar
            xmin(i) = p(i,ilo)
         end do
         y0 = y(ilo)
         done = .true.
         if (iter .ge. maxiter) then
            write (iout,60)
   60       format (/,' SIMPLEX  --  Maximum Number of Iterations',
     &                 ' Exceeded')
         else
            do i = 1, nvar
               del = step * eps
               xmin(i) = xmin(i) + del
               z = fvalue (xmin)
               iter = iter + 1
               if (z .lt. y0) then
                  done = .false.
                  goto 70
               end if
               xmin(i) = xmin(i) - del - del
               z = fvalue (xmin)
               iter = iter + 1
               if (z .lt. y0) then
                  done = .false.
                  goto 70
               end if
               xmin(i) = xmin(i) + del
            end do
         end if
   70    continue
c
c     set return to current minimum, restart if warranted
c
         do i = 1, nvar
            x0(i) = xmin(i)
         end do
         del = eps
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xmin)
      deallocate (pbar)
      deallocate (pstar)
      deallocate (p2star)
      deallocate (y)
      deallocate (p)
      return
      end
