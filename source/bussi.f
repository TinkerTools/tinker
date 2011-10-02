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
c     "bussi" performs a single molecular dynamics time step
c     via the Bussi-Parrinello isothermal-isobaric algorithm
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
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'boxes.i'
      include 'mdstuf.i'
      include 'moldyn.i'
      include 'shake.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,istep
      real*8 dt,dt_2
      real*8 etot,epot
      real*8 eksum,ekbaro
      real*8 temp,pres
      real*8 factor,term
      real*8 random,normal
      real*8 kt,scale,sum
      real*8 w,c,d,r,s,si
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 xold(maxatm)
      real*8 yold(maxatm)
      real*8 zold(maxatm)
      real*8 derivs(3,maxatm)
c
c
c     set some time values, constants and barostat mass
c
      dt_2 = 0.5d0 * dt
      kt = boltzmann * kelvin
      w = dble(nfree) * kt * taupres * taupres
c
c     get kinetic energy and temperature, including barostat
c
      call kinetic (eksum,ekin)
      ekbaro = 0.5d0 * w * eta * eta / convert
      eksum = eksum + ekbaro
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
c
c     make thermostat updates and couple to atomic velocities
c
      c = exp(-dt_2/tautemp)
      d = (1.0d0-c) * (kelvin/temp) / dble(nfree)
      r = normal ()
      s = 0.0d0
      do i = 1, nfree-1
         si = normal ()
         s = s + si*si
      end do
      scale = c + (s+r*r)*d + 2.0d0*r*sqrt(c*d)
      scale = sqrt(scale)
      if ((r+sqrt(c/d)) .lt. 0.0d0)  scale = -scale
      eta = scale * eta
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = scale * v(j,i)
            end do
         end if
      end do
c
c     calculate the stress tensor for anisotropic systems, and
c     set isotropic pressure to the average of tensor diagonal
c
      factor = prescon / volbox
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = factor * (2.0d0*ekin(j,i)-vir(j,i))
         end do
      end do
      pres = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0
c
c     get half-step velocities and update barostat velocity variable
c
      eta = eta + 3.0d0*(volbox*(pres-atmsph)*convert/prescon
     &                         + 2.0*kt)*dt_2/w
      do i = 1, n
         if (use(i)) then
            sum = 0.0d0
            do j = 1, 3
               eta = eta + mass(i)*a(j,i)*v(j,i)*dt_2*dt_2/w
               sum = sum + a(j,i)*a(j,i)
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
            eta = eta + mass(i)*abs(sum)*dt_2*dt_2*dt_2/(3.0d0*w)
        end if
      end do
c
c     store the current atom positions, then alter positions
c     and velocities via coupling to the barostat 
c
      term = eta * dt
      do i = 1, n
         if (use(i)) then
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
            x(i) = x(i)*exp(term) + v(1,i)*sinh(term)/eta
            y(i) = y(i)*exp(term) + v(2,i)*sinh(term)/eta
            z(i) = z(i)*exp(term) + v(3,i)*sinh(term)/eta
            do j = 1, 3
               v(j,i) = v(j,i) * exp(-term)
            end do
         end if
      end do
c
c     set the new box dimensions and other lattice values;
c     current version assumes isotropic pressure
c
      xbox = xbox * exp(term)
      ybox = ybox * exp(term)
      zbox = zbox * exp(term)
      call lattice
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
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               a(j,i) = -convert * derivs(j,i) / mass(i)
            end do
         end if
      end do
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     find new half-step velocities via Verlet recursion and
c     update barostat velocity variable
c
      eta = eta + 3.0d0*(volbox*(pres-atmsph)*convert/prescon
     &                         + 2.0*kt)*dt_2/w
      do i = 1, n
         if (use(i)) then
            sum = 0.0d0
            do j = 1, 3
               eta = eta + mass(i)*a(j,i)*v(j,i)*dt_2*dt_2/w
               sum = sum + a(j,i)*a(j,i)
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
            eta = eta + mass(i)*abs(sum)*dt_2*dt_2*dt_2/(3.0d0*w)
        end if
      end do
c
c     update thermostat and couple to atomic velocities 
c
      c = exp(-dt_2/tautemp)
      d = (1.0d0-c) * (kelvin/temp) / dble(nfree)
      r = normal ()
      s = 0.0d0
      do i = 1, nfree-1
         si = normal ()
         s = s + si*si
      end do
      scale = c + (s+r*r)*d + 2.0d0*r*sqrt(c*d)
      scale = sqrt(scale)
      if (r+sqrt(c/d) .lt. 0.0d0)  scale = -scale
      eta = scale * eta
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = scale * v(j,i)
            end do
         end if
      end do
c
c     get kinetic energy and temperature, including barostat
c
      call kinetic (eksum,ekin)
      ekbaro = 0.5d0 * w * eta * eta / convert
      eksum = eksum + ekbaro
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
c
c     calculate the stress tensor for anisotropic systems, and
c     set isotropic pressure to the average of tensor diagonal
c
      factor = prescon / volbox
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = factor * (2.0d0*ekin(j,i)-vir(j,i))
         end do
      end do
      pres = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0
c
c     system energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
      return
      end
