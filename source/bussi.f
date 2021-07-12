c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2010 by Teresa Head-Gordon & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine bussi  --  Bussi NPT molecular dynamics step  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "bussi" performs a single molecular dynamics time step via
c     the Bussi-Parrinello isothermal-isobaric algorithm
c
c     literature reference:
c
c     G. Bussi, T. Zykova-Timan and M. Parrinello, "Isothermal-Isobaric
c     Molecular Dynamics using Stochastic Velocity Rescaling", Journal
c     of Chemical Physics, 130, 074101 (2009)
c
c     original version written by Teresa Head-Gordon, October 2010
c
c
      subroutine bussi (istep,dt)
      use atomid
      use atoms
      use bath
      use boxes
      use freeze
      use ielscf
      use mdstuf
      use moldyn
      use polar
      use units
      use usage
      implicit none
      integer i,j,k
      integer istep
      real*8 dt,dt_2,dt_x
      real*8 dt2_2,dt3_2
      real*8 epot,etot,eksum
      real*8 expterm,sinhterm
      real*8 kt,w,temp,pres
      real*8 part1,part2
      real*8 factor,term
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
c
c
c     set some time values, constants and barostat mass
c
      dt_2 = 0.5d0 * dt
      dt2_2 = dt_2 * dt_2
      dt3_2 = dt2_2 * dt_2
      kt = boltzmann * kelvin
      w = dble(nfree) * kt * taupres * taupres
c
c     get Beeman integration coefficients for velocity updates
c
      factor = dble(bmnmix)
      dt_x = dt / factor
      part1 = 0.5d0*factor + 1.0d0
      part2 = part1 - 2.0d0
c
c     make half-step temperature correction and get pressure
c
      call temper (dt_2,eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
c
c     get half-step Beeman velocities and update barostat velocity
c
      eta = eta + 3.0d0*(volbox*(pres-atmsph)*ekcal/prescon
     &                         + 2.0*kt)*dt_2/w
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            eta = eta + mass(k)*a(j,k)*v(j,k)*dt2_2/w
     &               + mass(k)*a(j,k)*a(j,k)*dt3_2/(3.0d0*w)
            v(j,k) = v(j,k) + (part1*a(j,k)-aalt(j,k))*dt_x
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (derivs(3,n))
c
c     store the current atom positions, then alter positions
c     and velocities via coupling to the barostat
c
      term = eta * dt
      expterm = exp(term)
      sinhterm = sinh(term)
      do i = 1, nuse
         k = iuse(i)
         xold(k) = x(k)
         yold(k) = y(k)
         zold(k) = z(k)
         x(k) = x(k)*expterm + v(1,k)*sinhterm/eta
         y(k) = y(k)*expterm + v(2,k)*sinhterm/eta
         z(k) = z(k)*expterm + v(3,k)*sinhterm/eta
         do j = 1, 3
            v(j,k) = v(j,k) / expterm
         end do
      end do
c
c     set the new box dimensions and other lattice values;
c     current version assumes isotropic pressure
c
      xbox = xbox * expterm
      ybox = ybox * expterm
      zbox = zbox * expterm
      call lattice
c
c     apply Verlet half-step updates for any auxiliary dipoles
c
      if (use_ielscf) then
         do i = 1, nuse
            k = iuse(i)
            do j = 1, 3
               vaux(j,k) = vaux(j,k) + aaux(j,k)*dt_2
               vpaux(j,k) = vpaux(j,k) + apaux(j,k)*dt_2
               uaux(j,k) = uaux(j,k) + vaux(j,k)*dt
               upaux(j,k) = upaux(j,k) + vpaux(j,k)*dt
            end do
         end do
         call temper2 (dt,temp)
      end if
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt,xold,yold,zold)
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
c
c     use Newton's second law to get the next accelerations
c
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            aalt(j,k) = a(j,k)
            a(j,k) = -ekcal * derivs(j,k) / mass(k)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (derivs)
c
c     get full-step Beeman velocities and update barostat velocity
c
      eta = eta + 3.0d0*(volbox*(pres-atmsph)*ekcal/prescon
     &                         + 2.0*kt)*dt_2/w
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            eta = eta + mass(k)*a(j,k)*v(j,k)*dt2_2/w
     &               + mass(k)*a(j,k)*a(j,k)*dt3_2/(3.0d0*w)
            v(j,k) = v(j,k) + (part2*a(j,k)+aalt(j,k))*dt_x
         end do
      end do
c
c     apply Verlet full-step updates for any auxiliary dipoles
c
      if (use_ielscf) then
         term = 2.0d0 / (dt*dt)
         do i = 1, nuse
            k = iuse(i)
            do j = 1, 3
               aaux(j,k) = term * (uind(j,k)-uaux(j,k))
               apaux(j,k) = term * (uinp(j,k)-upaux(j,k))
               vaux(j,k) = vaux(j,k) + aaux(j,k)*dt_2
               vpaux(j,k) = vpaux(j,k) + apaux(j,k)*dt_2
            end do
         end do
      end if
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     make full-step temperature correction and get pressure
c
      call temper (dt_2,eksum,ekin,temp)
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
      call mdrest (istep)
      return
      end
