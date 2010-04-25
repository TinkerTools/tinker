c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine beeman  --  Beeman molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "beeman" performs a single molecular dynamics time step
c     by means of a Beeman multistep recursion formula; the
c     actual coefficients are Brooks' "Better Beeman" values
c
c     literature references:
c
c     D. Beeman, "Some Multistep Methods for Use in Molecular
c     Dynamics Calculations", Journal of Computational Physics,
c     20, 130-139 (1976)
c
c     B. R. Brooks, "Algorithms for Molecular Dynamics at Constant
c     Temperature and Pressure", DCRT Report, NIH, April 1988
c
c
      subroutine beeman (istep,dt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'moldyn.i'
      include 'shake.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,istep
      real*8 dt,etot
      real*8 dt_8,dt2_8
      real*8 eksum,epot
      real*8 temp,pres
      real*8 xterm,yterm,zterm
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 xold(maxatm)
      real*8 yold(maxatm)
      real*8 zold(maxatm)
      real*8 derivs(3,maxatm)
c
c
c     set some time values for the dynamics integration
c
      dt_8 = 0.125d0 * dt
      dt2_8 = dt * dt_8
c
c     make half-step temperature and pressure corrections
c
      call temper (dt)
c
c     store the current atom positions, then find new atom
c     positions and half-step velocities via Beeman recursion
c
      do i = 1, n
         if (use(i)) then
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
            xterm = 5.0d0*a(1,i) - aold(1,i)
            yterm = 5.0d0*a(2,i) - aold(2,i)
            zterm = 5.0d0*a(3,i) - aold(3,i)
            x(i) = x(i) + v(1,i)*dt + xterm*dt2_8
            y(i) = y(i) + v(2,i)*dt + yterm*dt2_8
            z(i) = z(i) + v(3,i)*dt + zterm*dt2_8
            v(1,i) = v(1,i) + xterm*dt_8
            v(2,i) = v(2,i) + yterm*dt_8
            v(3,i) = v(3,i) + zterm*dt_8
         end if
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
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Beeman recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               aold(j,i) = a(j,i)
               a(j,i) = -convert * derivs(j,i) / mass(i)
               v(j,i) = v(j,i) + (3.0d0*a(j,i)+aold(j,i))*dt_8
            end do
         end if
      end do
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     accumulate the kinetic energy and its outer product
c
      call kinetic (eksum,ekin)
c
c     make full-step temperature and pressure corrections
c
      call temper2 (dt,eksum,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
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
