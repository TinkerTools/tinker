c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2018  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine shake  --  SHAKE distance constraint method  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "shake" implements the SHAKE algorithm by correcting atomic
c     positions to maintain interatomic distance and absolute spatial
c     constraints
c
c     literature reference:
c
c     J. P. Ryckaert, G. Ciccotti and H. J. C. Berendsen, "Numerical
c     Integration of the Cartesian Equations of Motion of a System
c     with Constraints: Molecular Dynamics of n-Alkanes", Journal of
c     Computational Physics, 23, 327-341 (1977)
c
c
      subroutine shake (xold,yold,zold)
      use atomid
      use atoms
      use freeze
      use group
      use inform
      use iounit
      use usage
      implicit none
      integer i,j,k
      integer ia,ib,mode
      integer niter,maxiter
      integer start,stop
      real*8 eps,sor,dt
      real*8 xr,yr,zr
      real*8 xo,yo,zo
      real*8 dot,rma,rmb
      real*8 weigh,dist2
      real*8 delta,term
      real*8 xterm,yterm,zterm
      real*8 xold(*)
      real*8 yold(*)
      real*8 zold(*)
      logical done
      logical, allocatable :: moved(:)
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (moved(n))
      allocate (update(n))
c
c     initialize the lists of atoms previously corrected
c
      do i = 1, n
         if (use(i)) then
            moved(i) = .true.
         else
            moved(i) = .false.
         end if
         update(i) = .false.
      end do
c
c     set the iteration counter, termination and tolerance
c
      maxiter = 500
      sor = 1.25d0
      eps = rateps
c
c     apply SHAKE to adjust distances to constraint values
c
      niter = 0
      done = .false.
      do while (.not.done .and. niter.lt.maxiter)
         niter = niter + 1
         done = .true.
         do i = 1, nrat
            ia = irat(1,i)
            ib = irat(2,i)
            if (moved(ia) .or. moved(ib)) then
               xr = x(ib) - x(ia)
               yr = y(ib) - y(ia)
               zr = z(ib) - z(ia)
               if (frzimage(i))  call image (xr,yr,zr)
               dist2 = xr*xr + yr*yr + zr*zr
               delta = krat(i)**2 - dist2
               if (abs(delta) .gt. eps) then
                  done = .false.
                  update(ia) = .true.
                  update(ib) = .true.
                  xo = xold(ib) - xold(ia)
                  yo = yold(ib) - yold(ia)
                  zo = zold(ib) - zold(ia)
                  if (frzimage(i))  call image (xo,yo,zo)
                  dot = xr*xo + yr*yo + zr*zo
                  rma = 1.0d0 / mass(ia)
                  rmb = 1.0d0 / mass(ib)
                  term = 0.5d0 * sor * delta / ((rma+rmb) * dot)
                  xterm = xo * term
                  yterm = yo * term
                  zterm = zo * term
                  x(ia) = x(ia) - xterm*rma
                  y(ia) = y(ia) - yterm*rma
                  z(ia) = z(ia) - zterm*rma
                  x(ib) = x(ib) + xterm*rmb
                  y(ib) = y(ib) + yterm*rmb
                  z(ib) = z(ib) + zterm*rmb
               end if
            end if
         end do
         do i = 1, n
            moved(i) = update(i)
            update(i) = .false.
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (moved)
      deallocate (update)
c
c     write information on the number of iterations needed
c
      if (niter .eq. maxiter) then
         write (iout,10)
   10    format (/,' SHAKE  --  Warning, Distance Constraints',
     &              ' not Satisfied')
         call prterr
         call fatal
      else if (debug) then
         write (iout,20)  niter
   20    format (' SHAKE   --  Distance Constraints met at',i6,
     &              ' Iterations')
      end if
c
c     any rigid planar water molecules are handled separately
c
      dt = 0.0d0
      call settle (dt,xold,yold,zold)
c
c     apply any group position constraints via exact reset
c
      do i = 1, nratx
         ia = iratx(i)
         mode = kratx(i)
         xr = 0.0d0
         yr = 0.0d0
         zr = 0.0d0
         start = igrp(1,ia)
         stop = igrp(2,ia)
         do j = start, stop
            k = kgrp(j)
            weigh = mass(k) / grpmass(ia)
            if (mode .gt. 2) then
               xr = xr + x(k)*weigh
            end if
            if (mode .gt. 1) then
               yr = yr + y(k)*weigh
            end if
            zr = zr + z(k)*weigh
         end do
         do j = start, stop
            k = kgrp(j)
            x(k) = x(k) - xr
            y(k) = y(k) - yr
            z(k) = z(k) - zr
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine shakeg  --  SHAKE gradient vector correction  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "shakeg" modifies the gradient to remove components along any
c     holonomic distance contraints using a variant of SHAKE
c
c     literature reference:
c
c     Y. Duan, S. Kumar, J. M. Rosenberg and P. A. Kollman, "Gradient
c     SHAKE: An Improved Method for Constrained Energy Minimization in
c     Macromolecular Simulations", Journal of Computational Chemistry,
c     16, 1351-1356 (1995)
c
c
      subroutine shakeg (derivs)
      use atoms
      use freeze
      use inform
      use iounit
      use usage
      implicit none
      integer i,ia,ib
      integer niter,maxiter
      real*8 eps,sor
      real*8 xr,yr,zr
      real*8 xf,yf,zf
      real*8 dist2,delta,term
      real*8 xterm,yterm,zterm
      real*8 derivs(3,*)
      logical done
      logical, allocatable :: moved(:)
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (moved(n))
      allocate (update(n))
c
c     initialize the lists of atoms previously corrected
c
      do i = 1, n
         if (use(i)) then
            moved(i) = .true.
         else
            moved(i) = .false.
         end if
         update(i) = .false.
      end do
c
c     set the iteration counter, termination and tolerance
c
      maxiter = 500
      sor = 1.15d0
      eps = rateps
c
c     adjust the gradient to remove constraint components
c
      niter = 0
      done = .false.
      do while (.not.done .and. niter.lt.maxiter)
         niter = niter + 1
         done = .true.
         do i = 1, nrat
            ia = irat(1,i)
            ib = irat(2,i)
            if (moved(ia) .or. moved(ib)) then
               xr = x(ib) - x(ia)
               yr = y(ib) - y(ia)
               zr = z(ib) - z(ia)
               if (frzimage(i))  call image (xr,yr,zr)
               dist2 = xr*xr + yr*yr + zr*zr
               xf = derivs(1,ib) - derivs(1,ia)
               yf = derivs(2,ib) - derivs(2,ia)
               zf = derivs(3,ib) - derivs(3,ia)
               delta = xr*xf + yr*yf + zr*zf
               if (abs(delta) .gt. eps) then
                  done = .false.
                  update(ia) = .true.
                  update(ib) = .true.
                  term = 0.5d0 * sor * delta / dist2
                  xterm = xr * term
                  yterm = yr * term
                  zterm = zr * term
                  derivs(1,ia) = derivs(1,ia) + xterm
                  derivs(2,ia) = derivs(2,ia) + yterm
                  derivs(3,ia) = derivs(3,ia) + zterm
                  derivs(1,ib) = derivs(1,ib) - xterm
                  derivs(2,ib) = derivs(2,ib) - yterm
                  derivs(3,ib) = derivs(3,ib) - zterm
               end if
            end if
         end do
         do i = 1, n
            moved(i) = update(i)
            update(i) = .false.
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (moved)
      deallocate (update)
c
c     write information on the number of iterations needed
c
      if (niter .eq. maxiter) then
         write (iout,10)
   10    format (/,' SHAKEG  --  Warning, Gradient Constraints',
     &              ' not Satisfied')
         call prterr
         call fatal
      else if (debug) then
         write (iout,20)  niter
   20    format (' SHAKEG  --  Gradient Constraints met at',i6,
     &              ' Iterations')
      end if
c
c     any rigid planar water molecules are handled separately
c
      call settleg (derivs)
      return
      end
