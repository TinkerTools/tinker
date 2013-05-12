c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2011  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine respa  --  r-RESPA molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "respa" performs a single multiple time step molecular dynamics
c     step using the reversible reference system propagation algorithm
c     (r-RESPA) via a velocity Verlet core with a position Verlet inner
c     loop and potential split into fast- and slow-evolving portions
c
c     literature references:
c
c     D. D. Humphreys, R. A. Friesner and B. J. Berne, "A Multiple-
c     Time-Step Molecular Dynamics Algorithm for Macromolecules",
c     Journal of Physical Chemistry, 98, 6885-6892 (1994)
c
c     X. Qian and T. Schlick, "Efficient Multiple-Time-Step Integrators
c     with Distance-Based Force Splitting for Particle-Mesh-Ewald
c     Molecular Dynamics Simulations", Journal of Chemical Physics,
c     115, 4019-4029 (2001)
c
c
      subroutine respa (istep,dt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'freeze.i'
      include 'moldyn.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k
      integer istep
      integer nalt
      real*8 dt,dt_2
      real*8 dta,dta_2
      real*8 epot,etot
      real*8 eksum,eps
      real*8 temp,pres
      real*8 ealt,dalt
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 viralt(3,3)
      real*8 xold(maxatm)
      real*8 yold(maxatm)
      real*8 zold(maxatm)
      real*8 xolda(maxatm)
      real*8 yolda(maxatm)
      real*8 zolda(maxatm)
      real*8 derivs(3,maxatm)
c
c
c     set some time values for the dynamics integration
c
      eps =  0.00000001d0
      dalt = 0.00025d0
      nalt = int(dt/(dalt+eps)) + 1
      dalt = dble(nalt)
      dt_2 = 0.5d0 * dt
      dta = dt / dalt
      dta_2 = 0.5d0 * dta
c
c     make half-step temperature and pressure corrections
c
      call temper (dt)
c
c     store the current atom positions, then find half-step
c     velocities via velocity Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
         end if
      end do
c
c     initialize virial from fast-evolving potential energy terms
c
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = 0.0d0
         end do
      end do
c
c     find fast-evolving positions via position Verlet recursion
c
      do k = 1, nalt
         do i = 1, n
            if (use(i)) then
               xolda(i) = x(i)
               yolda(i) = y(i)
               zolda(i) = z(i)
               x(i) = x(i) + v(1,i)*dta_2
               y(i) = y(i) + v(2,i)*dta_2
               z(i) = z(i) + v(3,i)*dta_2
            end if
         end do
c
c     determine constraint-corrected half-step positions
c
         do i = 1, 3
            do j = 1, 3
               vir(j,i) = 0.0d0
            end do
         end do
         if (use_rattle)  call shake (dta_2,xolda,yolda,zolda)
         do i = 1, 3
            do j = 1, 3
               viralt(j,i) = viralt(j,i) + vir(j,i)/dalt
            end do
         end do
c
c     get the fast-evolving potential energy and atomic forces
c
         call gradfast (ealt,derivs)
c
c     use Newton's second law to get fast-evolving accelerations;
c     then velocities and positions via position Verlet recursion
c
         do i = 1, n
            if (use(i)) then
               do j = 1, 3
                  aalt(j,i) = -convert * derivs(j,i) / mass(i)
                  v(j,i) = v(j,i) + aalt(j,i)*dta
               end do
               xolda(i) = x(i)
               yolda(i) = y(i)
               zolda(i) = z(i)
               x(i) = x(i) + v(1,i)*dta_2
               y(i) = y(i) + v(2,i)*dta_2
               z(i) = z(i) + v(3,i)*dta_2
            end if
         end do
         if (use_rattle)  call shake (dta_2,xolda,yolda,zolda)
c
c     increment average virial from fast-evolving potential terms
c
         do i = 1, 3
            do j = 1, 3
               viralt(j,i) = viralt(j,i) + vir(j,i)/dalt
            end do
         end do
      end do
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt,xold,yold,zold)
c
c     get the slow-evolving potential energy and atomic forces
c
      call gradslow (epot,derivs)
c
c     use Newton's second law to get the slow accelerations;
c     find full-step velocities using velocity Verlet recursion
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
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     total potential and virial from sum of fast and slow parts
c
      epot = epot + ealt
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do
      end do
c
c     make full-step temperature and pressure corrections
c
      call temper2 (dt,eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
      call mdrest (istep)
      return
      end
