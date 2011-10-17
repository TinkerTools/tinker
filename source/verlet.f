c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine verlet  --  Verlet molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "verlet" performs a single molecular dynamics time step
c     via the velocity Verlet multistep recursion formula
c
c
      subroutine verlet (istep,dt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'moldyn.i'
      include 'shake.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,istep
      real*8 dt,dt_2,etot
      real*8 eksum,epot
      real*8 temp,pres
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
      dt_2 = 0.5d0 * dt
c
c     make half-step temperature and pressure corrections
c
      call temper (dt)
c
c     store the current atom positions, then find half-step
c     velocities and full-step positions via Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
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
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt,xold,yold,zold)
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
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
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
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
