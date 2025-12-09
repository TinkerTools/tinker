c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 2011 by John Chodera & Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ghmcstep  --  generalized hybrid Monte Carlo  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ghmcstep" performs a stochastic dynamics step via a generalized
c     hybrid Monte Carlo (GHMC) algorithm that ensures exact sampling
c     from the Boltzmann density
c
c     literature references:
c
c     T. Lelievre, M. Rousset and G. Stoltz, "Free Energy Computations:
c     A Mathematical Perspective", Imperial College Press, London, 2010
c     [Algorithm 2.11]
c
c     T. Lelievre, M. Rousset and G. Stoltz, "Langevin Dynamics
c     with Constraints and Computation of Free Energy Differences",
c     Mathematics of Computation, 81, 2071-2125 (2012) [eq 3.16-3.18]
c
c     original version written by John D. Chodera, University of
c     California, Berkeley, November 2010
c
c
      subroutine ghmcstep (istep,dt)
      use atoms
      use atomid
      use bath
      use freeze
      use inform
      use iounit
      use moldyn
      use units
      use usage
      use virial
      implicit none
      integer i,j,k
      integer istep
      integer nreject
      real*8 dt,dt_2
      real*8 energy
      real*8 epot,etot
      real*8 epold,etold
      real*8 eksum,de
      real*8 temp,pres
      real*8 random,ratio
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: vold(:,:)
      real*8, allocatable :: aold(:,:)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: alpha(:,:)
      real*8, allocatable :: beta(:,:)
      logical first
      external energy
      external random
      save nreject
      save epot
      save first
      data first  / .true. /
c
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
c
c     use current energy as previous value for initial step
c
      if (first) then
         first = .false.
         nreject = 0
         epot = energy ()
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (vold(3,n))
      allocate (aold(3,n))
      allocate (derivs(3,n))
      allocate (alpha(3,n))
      allocate (beta(3,n))
c
c     evolve velocities according to midpoint Euler for half-step
c
      call ghmcterm (istep,dt,alpha,beta)
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            v(j,k) = v(j,k)*alpha(j,k) + beta(j,k)
         end do
      end do
c
c     find constraint-corrected velocities prior to Verlet step
c
      if (use_freeze)  call rattle2 (dt)
c
c     find the kinetic energy and store the energy values
c
      call kinetic (eksum,ekin,temp)
      epold = epot
      etold = eksum + epot
c
c     store the current positions and derivatives, find half-step
c     velocities and full-step positions via Verlet recursion
c
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            vold(j,k) = v(j,k)
            aold(j,k) = a(j,k)
            v(j,k) = v(j,k) + a(j,k)*dt_2
         end do
         xold(k) = x(k)
         yold(k) = y(k)
         zold(k) = z(k)
         x(k) = x(k) + v(1,k)*dt
         y(k) = y(k) + v(2,k)*dt
         z(k) = z(k) + v(3,k)*dt
      end do
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_freeze)  call rattle (dt,xold,yold,zold)
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Verlet recursion
c
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            a(j,k) = -ekcal * derivs(j,k) / mass(k)
            v(j,k) = v(j,k) + a(j,k)*dt_2
         end do
      end do
c
c     find the constraint-corrected full-step velocities
c
      if (use_freeze)  call rattle2 (dt)
c
c     compute kinetic energy, temperature and total energy
c
      call kinetic (eksum,ekin,temp)
      etot = eksum + epot
c
c     accept or reject according to the Metropolis criterion;
c     note velocities have flipped sign upon rejection
c
      de = (etot-etold) / (gasconst*kelvin)
      if (de.gt.0.0d0 .and. random().gt.exp(-de)) then
         nreject = nreject + 1
         epot = epold
         do i = 1, nuse
            k = iuse(i)
            x(k) = xold(k)
            y(k) = yold(k)
            z(k) = zold(k)
            do j = 1, 3
               v(j,k) = -vold(j,k)
               a(j,k) = aold(j,k)
            end do
         end do
      end if
      if (mod(istep,1000) .eq. 0) then
         ratio = 1.0d0 - dble(nreject)/1000.0d0
         nreject = 0
         write (iout,10)  ratio
   10    format (/,' GHMC Acceptance Ratio',6x,f8.3,
     &              ' for the Last 1000 Steps')
      end if
c
c     update velocities using midpoint Euler for half-step
c
      call ghmcterm (istep,dt,alpha,beta)
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            v(j,k) = v(j,k)*alpha(j,k) + beta(j,k)
         end do
      end do
c
c     update the constraint-corrected full-step velocities
c
      if (use_freeze) then
         call rattle2 (dt)
         do i = 1, nuse
            k = iuse(i)
            xold(k) = x(k)
            yold(k) = y(k)
            zold(k) = z(k)
         end do
      end if
c
c     compute full-step kinetic energy and pressure correction
c
      call kinetic (eksum,ekin,temp)
      call pressure (dt,ekin,pres,stress)
      call pressure2 (epot,temp)
c
c     final constraint step to enforce position convergence
c
      if (use_freeze)  call shake (xold,yold,zold)
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (vold)
      deallocate (aold)
      deallocate (derivs)
      deallocate (alpha)
      deallocate (beta)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot,eksum)
      call mdrest (istep)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ghmcterm  --  GHMC friction & fluctuation terms  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ghmcterm" finds the friction and fluctuation terms needed
c     to update velocities during GHMC stochastic dynamics
c
c
      subroutine ghmcterm (istep,dt,alpha,beta)
      use atoms
      use atomid
      use bath
      use stodyn
      use units
      use usage
      implicit none
      integer i,j,k
      integer istep
      real*8 dt,dt_4
      real*8 gamma,sigma
      real*8 term,normal
      real*8 alpha(3,*)
      real*8 beta(3,*)
      logical first
      external normal
      save first
      data first  / .true. /
c
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(fgamma))  allocate (fgamma(n))
c
c     set atomic friction coefficients to the global value
c
         do i = 1, n
            fgamma(i) = friction
         end do
      end if
c
c     find friction coefficients scaled by accessibility
c
      if (use_sdarea)  call sdarea (istep)
c
c     compute the viscous friction and fluctuation terms
c
      dt_4 = 0.25d0 * dt
      term = sqrt(boltzmann*kelvin*dt)
      do i = 1, nuse
         k = iuse(i)
         gamma = dt_4 * fgamma(k)
         sigma = term * sqrt(fgamma(k)/mass(k))
         do j = 1, 3
            alpha(j,k) = (1.0d0-gamma) / (1.0d0+gamma)
            beta(j,k) = normal() * sigma / (1.0d0+gamma)
         end do
      end do
      return
      end
