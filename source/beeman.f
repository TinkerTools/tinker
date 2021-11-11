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
c     via the Beeman multistep recursion formula; uses original
c     coefficients or Bernie Brooks' "Better Beeman" values
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
      use atomid
      use atoms
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
      real*8 dt,dt_x,factor
      real*8 etot,eksum,epot
      real*8 temp,pres
      real*8 part1,part2
      real*8 dt_2,term
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
c
c
c     set time values and coefficients for Beeman integration
c
      factor = dble(bmnmix)
      dt_x = dt / factor
      part1 = 0.5d0*factor + 1.0d0
      part2 = part1 - 2.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (derivs(3,n))
c
c     store the current atom positions, then find half-step
c     velocities and full-step positions via Beeman recursion
c
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            v(j,k) = v(j,k) + (part1*a(j,k)-aalt(j,k))*dt_x
         end do
         xold(k) = x(k)
         yold(k) = y(k)
         zold(k) = z(k)
         x(k) = x(k) + v(1,k)*dt
         y(k) = y(k) + v(2,k)*dt
         z(k) = z(k) + v(3,k)*dt
      end do
c
c     apply Verlet half-step updates for any auxiliary dipoles
c
      if (use_ielscf) then
         dt_2 = 0.5d0 * dt
         do i = 1, nuse
            k = iuse(i)
            do j = 1, 3
               vaux(j,k) = vaux(j,k) + aaux(j,k)*dt_2
               vpaux(j,k) = vpaux(j,k) + apaux(j,k)*dt_2
               uaux(j,k) = uaux(j,k) + vaux(j,k)*dt
               upaux(j,k) = upaux(j,k) + vpaux(j,k)*dt
            end do
         end do
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
c     make half-step temperature and pressure corrections
c
      call temper2 (dt,temp)
      call pressure2 (epot,temp)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Beeman recursion
c
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            aalt(j,k) = a(j,k)
            a(j,k) = -ekcal * derivs(j,k) / mass(k)
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
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (derivs)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     make full-step temperature and pressure corrections
c
      call temper (dt,eksum,ekin,temp)
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
