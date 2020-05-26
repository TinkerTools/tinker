c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 2011 by John Chodera & Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ghmcstep  --  generalized hybrid MC time step  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ghmcstep" performs a single stochastic dynamics time step via
c     the generalized hybrid Monte Carlo (GHMC) algorithm to ensure
c     exact sampling from the Boltzmann density
c
c     literature references:
c
c     T. Lelievre, M. Rousset and G. Stoltz, "Free Energy Computations:
c     A Mathematical Perspective", Imperial College Press, London, 2010,
c     Algorithm 2.11
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
      use iounit
      use mdstuf
      use moldyn
      use units
      use usage
      use virial
      implicit none
      integer i,j,k
      integer istep,nrej
      real*8 dt,dt_2
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
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: alpha(:,:)
      real*8, allocatable :: beta(:,:)
      external random
      save epot,nrej
c
c
c     compute the half time step value
c
      dt_2 = 0.5d0 * dt
c
c     perform dynamic allocation of some local arrays
c
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
c     accumulate the kinetic energy and store the energy values
c
      call kinetic (eksum,ekin,temp)
      epold = epot
      etold = eksum + epot
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (vold(3,n))
      allocate (derivs(3,n))
c
c     store the current positions and velocities, find half-step
c     velocities and full-step positions via Verlet recursion
c
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            vold(j,k) = v(j,k)
            aalt(j,k) = a(j,k)
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
      if (use_rattle)  call rattle (dt,xold,yold,zold)
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
c
c     use current values as previous energies for first step
c
      if (istep .eq. 1) then
         nrej = 0
         epold = epot
         etold = eksum + epot
      end if
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
      if (use_rattle)  call rattle2 (dt)
c
c     determine the kinetic energy, temperature and total energy
c
      call kinetic (eksum,ekin,temp)
      etot = eksum + epot
c
c     accept or reject according to Metropolis scheme;
c     note that velocities are reversed upon rejection
c
      de = (etot-etold) / (gasconst*kelvin)
      if (de.gt.0.0d0 .and. random().gt.exp(-de)) then
         nrej = nrej + 1
         ratio = 1.0d0 - dble(nrej)/dble(istep)
         write (iout,10)  ratio
   10    format (' GHMC Step Rejected',6x,'Acceptance Ratio',f8.3)
         epot = epold
         do i = 1, nuse
            k = iuse(i)
            x(k) = xold(k)
            y(k) = yold(k)
            z(k) = zold(k)
            do j = 1, 3
               v(j,k) = -vold(j,k)
               a(j,k) = aalt(j,k)
            end do
         end do
      end if
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
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (vold)
      deallocate (derivs)
      deallocate (alpha)
      deallocate (beta)
c
c     update the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     compute and control the temperature and pressure
c
      call kinetic (eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot,eksum)
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
      real*8 dt,dt_2,dt_4
      real*8 gamma,sigma
      real*8 normal
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
c     set the atomic friction coefficients to the global value
c
         do i = 1, n
            fgamma(i) = friction * mass(i)
         end do
      end if
c
c     set the value of the friction coefficient for each atom
c
      if (use_sdarea)  call sdarea (istep)
c
c     get the viscous friction and fluctuation terms for GHMC
c
      dt_2 = 0.5d0 * dt
      dt_4 = 0.25d0 * dt
      do i = 1, nuse
         k = iuse(i)
         gamma = dt_4 * fgamma(k) / mass(k)
         sigma = sqrt(2.0d0*boltzmann*kelvin*fgamma(k))
         do j = 1, 3
            alpha(j,k) = (1.0d0-gamma) / (1.0d0+gamma)
            beta(j,k) = normal() * sqrt(dt_2) * sigma
     &                     / ((1.0d0+gamma)*mass(k))
         end do
      end do
      return
      end
