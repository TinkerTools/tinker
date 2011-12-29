c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine shake  --  SHAKE distance constraint algorithm  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "shake" implements the SHAKE algorithm by correcting atomic
c     positions to maintain interatomic distance and absolute spatial
c     constraints
c
c     literature reference:
c
c     J.-P. Ryckaert, G. Cicotti, H. J. C. Berendsen, "Numerical
c     Integration of the Cartesian Equations of Motion of a System
c     with Constraints: Molecular Dynamics of n-Alkanes", Journal
c     of Computational Physics, 23, 327-341 (1977)
c
c
      subroutine shake (dt,xold,yold,zold)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'freeze.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'moldyn.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k
      integer ia,ib,mode
      integer niter,maxiter
      integer start,stop
      real*8 dt,eps,sor
      real*8 xr,yr,zr
      real*8 xo,yo,zo
      real*8 dot,rma,rmb
      real*8 delta,dist2
      real*8 weigh,vterm,term
      real*8 xterm,yterm,zterm
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 xold(maxatm)
      real*8 yold(maxatm)
      real*8 zold(maxatm)
      logical done
      logical moved(maxatm)
      logical update(maxatm)
c
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
      vterm = 1.0d0 / (dt*dt*convert)
c
c     apply shake to the list of pairwise distance values
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
               if (ratimage(i))  call image (xr,yr,zr)
               dist2 = xr**2 + yr**2 + zr**2
               delta = krat(i)**2 - dist2
               if (abs(delta) .gt. eps) then
                  done = .false.
                  update(ia) = .true.
                  update(ib) = .true.
                  xo = xold(ib) - xold(ia)
                  yo = yold(ib) - yold(ia)
                  zo = zold(ib) - zold(ia)
                  if (ratimage(i))  call image (xo,yo,zo)
                  dot = xr*xo + yr*yo + zr*zo
                  rma = 1.0d0 / mass(ia)
                  rmb = 1.0d0 / mass(ib)
                  term = sor * delta / (2.0d0 * (rma+rmb) * dot)
                  xterm = xo * term
                  yterm = yo * term
                  zterm = zo * term
                  x(ia) =  x(ia) - xterm*rma
                  y(ia) =  y(ia) - yterm*rma
                  z(ia) =  z(ia) - zterm*rma
                  x(ib) =  x(ib) + xterm*rmb
                  y(ib) =  y(ib) + yterm*rmb
                  z(ib) =  z(ib) + zterm*rmb
c
c     increment the internal virial tensor components
c     (Note: this virial code seems to not be correct....)
c
                  term = term * vterm
                  xterm = xr * term
                  yterm = yr * term
                  zterm = zr * term
                  vxx = xr * xterm
                  vyx = yr * xterm
                  vzx = zr * xterm
                  vyy = yr * yterm
                  vzy = zr * yterm
                  vzz = zr * zterm
                  vir(1,1) = vir(1,1) - vxx
                  vir(2,1) = vir(2,1) - vyx
                  vir(3,1) = vir(3,1) - vzx
                  vir(1,2) = vir(1,2) - vyx
                  vir(2,2) = vir(2,2) - vyy
                  vir(3,2) = vir(3,2) - vzy
                  vir(1,3) = vir(1,3) - vzx
                  vir(2,3) = vir(2,3) - vzy
                  vir(3,3) = vir(3,3) - vzz
               end if
            end if
         end do
         do i = 1, n
            moved(i) = update(i)
            update(i) = .false.
         end do
      end do
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
c     apply group position constraints via exact reset
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
            if (mode .gt. 2)  xr = xr + x(k)*weigh
            if (mode .gt. 1)  yr = yr + y(k)*weigh
            zr = zr + z(k)*weigh
         end do
         do j = start, stop
            k = kgrp(j)
            x(k) =  x(k) - xr
            y(k) =  y(k) - yr
            z(k) =  z(k) - zr
         end do
      end do
      return
      end
