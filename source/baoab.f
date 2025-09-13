c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 2016 by Charles Matthews & Ben Leimkuhler  ##
c     ##        COPYRIGHT (C)  2017  by  Jay William Ponder        ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine baoab  --  BAOAB stochastic dynamics step  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "baoab" implements a stochastic dynamics time step using
c     a geodesic BAOAB scheme with optional holonomic constraints
c
c     literature reference:
c
c     B. Leimkuhler and C. Matthews, "Efficient Molecular Dynamics
c     Using Geodesic Integration and Solvent-Solute Splitting",
c     Proceedings of the Royal Society A, 472, 20160138 (2016)
c
c     J. Jung and Y. Sugita, "Langevin Integration for Isothermal-
c     Isobaric Condition with a Large Time Step", Journal of Chemical
c     Physics, 162, 104108 (2025)
c
c
      subroutine baoab (istep,dt)
      use atomid
      use atoms
      use freeze
      use moldyn
      use units
      use usage
      use virial
      implicit none
      integer i,j,k
      integer istep
      integer nrattle
      real*8 dt,dt_2,dtr
      real*8 etot,epot
      real*8 eksum
      real*8 temp,pres
      real*8 drattle
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 virrat(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: vfric(:)
      real*8, allocatable :: vrand(:,:)
      real*8, allocatable :: derivs(:,:)
c
c
c     set some time values for the dynamics integration
c
      nrattle = 1
      if (use_rattle)  nrattle = 5
      drattle = dble(nrattle)
      dt_2 = 0.5d0 * dt
      dtr = dt_2 / drattle
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (vfric(n))
      allocate (vrand(3,n))
      allocate (derivs(3,n))
c
c     use a first B step to find the half-step velocities
c 
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            v(j,k) = v(j,k) + a(j,k)*dt_2
         end do
      end do
c
c     take the first A step to get the half-step positions
c
      do j = 1, nrattle
         do i = 1, nuse
            k = iuse(i)
            xold(k) = x(k)
            yold(k) = y(k)
            zold(k) = z(k)
            x(k) = x(k) + v(1,k)*dtr
            y(k) = y(k) + v(2,k)*dtr
            z(k) = z(k) + v(3,k)*dtr
         end do
         if (use_rattle) then
            call rattle (dtr,xold,yold,zold)
            call rattle2 (dtr)
            do i = 1, 3
               do k = 1, 3
                  vir(k,i) = 0.0d0
               end do
            end do
         end if
      end do
c
c     use an O step to get frictional and random components
c
      call oprep (dt,vfric,vrand)
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3 
            v(j,k) = v(j,k)*vfric(k) + vrand(j,k)
         end do
      end do
      if (use_rattle) then
         call rattle2 (dt)
         do i = 1, 3
            do j = 1, 3
               virrat(j,i) = vir(j,i)
               vir(j,i) = 0.0d0
            end do
         end do
      end if
c
c     take a second A step to get the full-step positions
c
      do j = 1, nrattle
         do i = 1, nuse
            k = iuse(i)
            xold(k) = x(k)
            yold(k) = y(k)
            zold(k) = z(k)
            x(k) = x(k) + v(1,k)*dtr
            y(k) = y(k) + v(2,k)*dtr
            z(k) = z(k) + v(3,k)*dtr
         end do
         if (use_rattle) then
            call rattle (dtr,xold,yold,zold)
            call rattle2 (dtr)
            do i = 1, 3
               do k = 1, 3
                  virrat(k,i) = virrat(k,i) + vir(k,i)/drattle
                  vir(k,i) = 0.0d0
               end do
            end do
         end if
      end do
c
c     get the potential energy and atomic forces
c 
      call gradient (epot,derivs)
c
c     compute the kinetic energy from half-step velocities
c
      call kinetic (eksum,ekin,temp)
c
c     second B step for accelerations and full-step velocities
c
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            a(j,k) = -ekcal * derivs(j,k) / mass(k)
            v(j,k) = v(j,k) + a(j,k)*dt_2
         end do
      end do 
      if (use_rattle) then
         call rattle2 (dt)
         do i = 1, 3
            do j = 1, 3
               vir(j,i) = vir(j,i) + virrat(j,i)
            end do
         end do
         do i = 1, nuse
            k = iuse(i)
            xold(k) = x(k)
            yold(k) = y(k)
            zold(k) = z(k)
         end do
      end if
c
c     compute full-step kinetic energy and pressure correction;
c     half-step kinetic energy gives better pressure control
c
c     call kinetic (eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
      call pressure2 (epot,temp)
c
c     final constraint step to enforce position convergence
c
      if (use_rattle)  call shake (xold,yold,zold)
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (vfric)
      deallocate (vrand)
      deallocate (derivs)
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
c     ############################################################
c     ##                                                        ##
c     ##  subroutine obabo  --  OBABO stochastic dynamics step  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "obabo" implements a stochastic dynamics time step using
c     a geodesic OBABO scheme with optional holonomic constraints
c
c     literature reference:
c
c     B. Leimkuhler and C. Matthews, "Efficient Molecular Dynamics
c     Using Geodesic Integration and Solvent-Solute Splitting",
c     Proceedings of the Royal Society A, 472, 20160138 (2016)
c
c     G. Bussi and M. Parrinello, "Accurate Sampling Using Langevin
c     Dynamics", Physical Review E, 75, 056707 (2007)
c
c
      subroutine obabo (istep,dt)
      use atomid
      use atoms
      use freeze
      use moldyn
      use units
      use usage
      use virial
      implicit none
      integer i,j,k
      integer istep
      integer nrattle
      real*8 dt,dt_2,dtr
      real*8 etot,epot
      real*8 eksum
      real*8 temp,pres
      real*8 drattle
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 virrat(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: vfric(:)
      real*8, allocatable :: vrand(:,:)
      real*8, allocatable :: derivs(:,:)
c
c
c     set some time values for the dynamics integration
c
      nrattle = 1
      if (use_rattle)  nrattle = 5
      drattle = dble(nrattle)
      dt_2 = 0.5d0 * dt
      dtr = dt / drattle
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (vfric(n)) 
      allocate (vrand(3,n))
      allocate (derivs(3,n))
c
c     take first O step for frictional and random components
c
      call oprep (dt_2,vfric,vrand)
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            v(j,k) = v(j,k)*vfric(k) + vrand(j,k)
         end do
      end do
      if (use_rattle) then
         do i = 1, 3
            do j = 1, 3
               vir(j,i) = 0.0d0
            end do
         end do
         call rattle2 (dt_2)
         do i = 1, 3
            do j = 1, 3
               virrat(j,i) = vir(j,i)
               vir(j,i) = 0.0d0
            end do
         end do
      end if
c
c     use a first B step to find the half-step velocities
c 
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            v(j,k) = v(j,k) + a(j,k)*dt_2
         end do
      end do
c
c     take full-step A step according to the BAOAB sequence
c
      do j = 1, nrattle
         do i = 1, nuse
            k = iuse(i)
            xold(k) = x(k)
            yold(k) = y(k)
            zold(k) = z(k)
            x(k) = x(k) + v(1,k)*dtr
            y(k) = y(k) + v(2,k)*dtr
            z(k) = z(k) + v(3,k)*dtr
         end do
         if (use_rattle) then
            call rattle (dtr,xold,yold,zold)
            call rattle2 (dtr)
            do i = 1, 3
               do k = 1, 3
                  virrat(k,i) = virrat(k,i) + vir(k,i)/drattle
                  vir(k,i) = 0.0d0
               end do
            end do
         end if
      end do
c
c     get the potential energy and atomic forces
c 
      call gradient (epot,derivs)
c
c     compute the kinetic energy from half-step velocities
c
c     call kinetic (eksum,ekin,temp)
c
c     second B step for accelerations and full-step velocities
c
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            a(j,k) = -ekcal * derivs(j,k) / mass(k)
            v(j,k) = v(j,k) + a(j,k)*dt_2
         end do
      end do 
      if (use_rattle)  call rattle2 (dt)
c
c     use second O step for frictional and random components
c
      call oprep (dt_2,vfric,vrand)
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3 
            v(j,k) = v(j,k)*vfric(k) + vrand(j,k)
         end do
      end do
      if (use_rattle) then
         call rattle2 (dt_2)
         do i = 1, 3
            do j = 1, 3
               vir(j,i) = vir(j,i) + virrat(j,i)
            end do
         end do
         do i = 1, nuse
            k = iuse(i)
            xold(k) = x(k)
            yold(k) = y(k)
            zold(k) = z(k)
         end do
      end if
c
c     compute the kinetic energy and control the pressure
c
      call kinetic (eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
      call pressure2 (epot,temp)
c
c     final constraint step to enforce position convergence
c
      if (use_rattle)  call shake (xold,yold,zold)
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (vfric) 
      deallocate (vrand)
      deallocate (derivs)
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine oprep  --  frictional & random terms for BAOAB  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "oprep" sets up frictional and random terms needed to update
c     positions and velocities for the BAOAB and OBABO integrators
c
c
      subroutine oprep (dt,vfric,vrand)
      use atoms
      use atomid
      use bath
      use stodyn
      use units
      use usage
      implicit none
      integer i,j,k
      real*8 dt,ktm
      real*8 egdt,normal
      real*8 vsig,vnorm
      real*8 vfric(*)
      real*8 vrand(3,*)
      logical first
      external normal
      save first
      data first  / .true. /
c
c
c     set the atomic friction coefficients to the global value
c
      if (first) then
         first = .false.
         if (.not. allocated(fgamma))  allocate (fgamma(n))
         do i = 1, n
            fgamma(i) = friction
         end do
      end if
c
c     get the frictional and random terms for a BAOAB step
c
      egdt = exp(-friction*dt)
      do i = 1, nuse
         k = iuse(i)
         vfric(k) = egdt
         ktm = boltzmann * kelvin / mass(k)
         vsig = sqrt(ktm*(1.0d0-egdt*egdt))
         do j = 1, 3
            vnorm = normal ()
            vrand(j,k) = vsig * vnorm
         end do
      end do
      return
      end
