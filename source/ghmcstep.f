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
c     T. Lelievre, M. Rousset and G. Stoltz, "Langevin Dynamics with
c     Constraints and Computation of Free Energy Differences", arXiv:
c     1006.4914v2, 18 Apr 2011, Equations 3.16-3.18
c
c     original version written by John D. Chodera, University of
c     California, Berkeley, November 2010
c
c
      subroutine ghmcstep (istep,dt)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmtyp.i'
      include 'bath.i'
      include 'iounit.i'
      include 'mdstuf.i'
      include 'moldyn.i'
      include 'shake.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j
      integer istep,nrej
      real*8 dt,dt_2
      real*8 etot,eksum,epot
      real*8 epold,etold,de
      real*8 temp,pres,ratio
      real*8 random,term
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 xold(maxatm)
      real*8 yold(maxatm)
      real*8 zold(maxatm)
      real*8 vold(3,maxatm)
      real*8 derivs(3,maxatm)
      real*8 alpha(3,maxatm)
      real*8 beta(3,maxatm)
      save nrej
      save epot
c
c
c     compute the half time step value
c
      dt_2 = 0.5d0 * dt
c
c     evolve momenta according to midpoint Euler for half timestep
c
      call ghmcterm (istep,dt,alpha,beta)
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = v(j,i)*alpha(j,i) + beta(j,i)
            end do
         end if
      end do      
c
c     accumulate the kinetic energy and store the energy values
c
      call kinetic (eksum,ekin)
      epold = epot
      etold = eksum + epot
c
c     store the current positions and velocities, find half-step
c     velocities and full-step positions via Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               vold(j,i) = v(j,i)
               aold(j,i) = a(j,i)
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
            x(i) = x(i) + v(1,i)*dt 
            y(i) = y(i) + v(2,i)*dt
            z(i) = z(i) + v(3,i)*dt
         end if
      end do
c
c     get constraint-corrected positions (TODO: Check to make
c     sure velocities are supposed to be corrected here)
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
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               a(j,i) = -convert * derivs(j,i) / mass(i)
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
         end if
      end do
c
c     find the constraint-corrected velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     accumulate the kinetic energy, temperature and total energy
c
      call kinetic (eksum,ekin)
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
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
         do i = 1, n
            if (use(i)) then
               x(i) = xold(i)
               y(i) = yold(i)
               z(i) = zold(i)
               do j = 1, 3
                  v(j,i) = -vold(j,i)
                  a(j,i) = aold(j,i)
               end do
            end if
         end do         
      end if
c
c     evolve momenta according to midpoint Euler for half timestep
c
      call ghmcterm (istep,dt,alpha,beta)
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = v(j,i)*alpha(j,i) + beta(j,i)
            end do
         end if
      end do      
c
c     find the constraint-corrected velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     compute and control the temperature and pressure
c
      call kinetic (eksum,ekin)
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
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
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmtyp.i'
      include 'bath.i'
      include 'stodyn.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,istep
      real*8 dt,dt_2,dt_4
      real*8 normal,gamma,sigma
      real*8 alpha(3,maxatm)
      real*8 beta(3,maxatm)
c
c
c     set the atomic friction coefficients to the global value
c
      if (istep .eq. 1) then
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
      do i = 1, n
         if (use(i)) then
            gamma = dt_4 * fgamma(i) / mass(i)
            sigma = sqrt(2.0d0*boltzmann*kelvin*fgamma(i))
            do j = 1, 3
               alpha(j,i) = (1.0d0-gamma) / (1.0d0+gamma)
               beta(j,i) = normal() * sqrt(dt_2) * sigma 
     &                        / ((1.0d0+gamma)*mass(i))
            end do
         end if
      end do
      return
      end
