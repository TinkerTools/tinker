c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine palmer  --  shake/rattle constraint algorithm  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "palmer" maintains constrained interatomic distances during
c     molecular dynamics by correcting positions and velocities
c     via application of Palmer's modified "shake/rattle" algorithm
c
c     literature reference:
c
c     B. J. Palmer, "Direct Application of SHAKE to the Velocity
c     Verlet Algorithm", Journal of Computational Physics, 104,
c     470-472 (1993)
c
c
      subroutine palmer (dt,xold,yold,zold)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'moldyn.i'
      include 'shake.i'
      include 'usage.i'
      integer i,ia,ib,niter,maxiter
      real*8 dt,xr,yr,zr,xo,yo,zo
      real*8 eps,dist2,delta,rma,rmb
      real*8 dot,term,xterm,yterm,zterm
      real*8 xold(maxatm),yold(maxatm),zold(maxatm)
      logical done,moved(maxatm),update(maxatm)
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
      maxiter = 100
      niter = 0
      done = .false.
      eps = 0.000001d0
c
c     apply iterative correction to positions and velocities
c
      dowhile (.not.done .and. niter.lt.maxiter)
         niter = niter + 1
         done = .true.
         do i = 1, nrat
            ia = irat(1,i)
            ib = irat(2,i)
            if (moved(ia) .or. moved(ib)) then
               xr = x(ib) - x(ia)
               yr = y(ib) - y(ia)
               zr = z(ib) - z(ia)
               dist2 = xr**2 + yr**2 + zr**2
               delta = krat(i)**2 - dist2
               if (abs(delta) .gt. eps) then
                  done = .false.
                  update(ia) = .true.
                  update(ib) = .true.
                  xo = xold(ib) - xold(ia)
                  yo = yold(ib) - yold(ia)
                  zo = zold(ib) - zold(ia)
                  dot = xr*xo + yr*yo + zr*zo
                  rma = 1.0d0 / mass(ia)
                  rmb = 1.0d0 / mass(ib)
                  term = delta / (2.0d0 * (rma+rmb) * dot)
                  xterm = xo * term
                  yterm = yo * term
                  zterm = zo * term
                  x(ia) = x(ia) - xterm*rma
                  y(ia) = y(ia) - yterm*rma
                  z(ia) = z(ia) - zterm*rma
                  x(ib) = x(ib) + xterm*rmb
                  y(ib) = y(ib) + yterm*rmb
                  z(ib) = z(ib) + zterm*rmb
                  rma = rma / dt
                  rmb = rmb / dt
                  v(1,ia) = v(1,ia) - xterm*rma
                  v(2,ia) = v(2,ia) - yterm*rma
                  v(3,ia) = v(3,ia) - zterm*rma
                  v(1,ib) = v(1,ib) + xterm*rmb
                  v(2,ib) = v(2,ib) + yterm*rmb
                  v(3,ib) = v(3,ib) + zterm*rmb
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
   10    format (' PALMER  --  Warning, Position Constraints',
     &              ' not Satisfied')
         call fatal
      else if (debug) then
         write (iout,20)  niter
   20    format (' PALMER  --  Position Constraints met at',
     &              i4,' Iterations')
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine verlet  --  Verlet molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "verlet" performs a single molecular dynamics time step
c     by means of the velocity Verlet multistep recursion formula
c
c
      subroutine verlet (istep,dt,dt_2,dt2_2,ndump)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'moldyn.i'
      include 'shake.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,istep,ndump
      real*8 dt,dt_2,dt2_2,e_tot,e_kin,e_pot
      real*8 temp,pres,derivs(3,maxatm)
      real*8 x_old(maxatm),y_old(maxatm),z_old(maxatm)
c
c
c     store the current atom positions, then find new atom
c     positions and half-step velocities via Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            x_old(i) = x(i)
            y_old(i) = y(i)
            z_old(i) = z(i)
            x(i) = x(i) + v(1,i)*dt + a(1,i)*dt2_2
            y(i) = y(i) + v(2,i)*dt + a(2,i)*dt2_2
            z(i) = z(i) + v(3,i)*dt + a(3,i)*dt2_2
            do j = 1, 3
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
         end if
      end do
c
c     apply "rattle" to correct atom positions and velocities
c
      if (use_rattle)  call palmer (dt,x_old,y_old,z_old)
c
c     get the potential energy and atomic forces
c
      call gradient (e_pot,derivs)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               a(j,i) = -convert * derivs(j,i) / mass(i)
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
         end if
      end do
c
c     sum the total kinetic energy over all atoms
c
      e_kin = 0.0d0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               e_kin = e_kin + mass(i)*v(j,i)**2
            end do
         end if
      end do
      e_kin = 0.5d0 * e_kin / convert
c
c     determine system temperature and total energy
c
      temp = 2.0d0 * e_kin / (dble(3*nuse-nrat-6) * gasconst)
      e_tot = e_kin + e_pot
c
c     control temperature and pressure via external bath coupling
c
      if (isothermal)  call temper (dt,temp)
      if (isobaric)  call pressure (dt,pres,e_kin)
c
c     compute any averages or statistics for this step
c
      call mdstat (istep,dt,e_tot,e_pot,e_kin,temp,pres,ndump)
      return
      end
