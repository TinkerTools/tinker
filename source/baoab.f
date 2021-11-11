c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 2016 by Charles Matthews & Ben Leimkuhler  ##
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
c     "baoab" implements a constrained stochastic dynamics time
c     step using the geodesic BAOAB scheme
c
c     literature reference:
c
c     B. Leimkuhler and C. Matthews, "Efficient Molecular Dynamics
c     Using Geodesic Integration and Solvent-Solute Splitting",
c     Proceedings of the Royal Society A, 472, 20160138 (2016)
c
c
      subroutine baoab (istep,dt)
      use atomid
      use atoms
      use freeze
      use mdstuf
      use moldyn
      use units
      use usage
      use virial  
      use limits
      use potent
      implicit none
      integer i,j,k
      integer istep
      integer nrattle
      real*8 dt,dtr
      real*8 dt_2,dt_4
      real*8 etot,epot
      real*8 eksum
      real*8 temp,pres
      real*8 ekin(3,3)
      real*8 stress(3,3)
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
      dt_2 = 0.5d0 * dt
      dt_4 = 0.25d0 * dt
      nrattle = 1
      dtr = dt_2 / dble(nrattle)
c
c     make half-step temperature and pressure corrections
c
      call temper (dt,eksum,ekin,temp)
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (derivs(3,n))
      allocate (vfric(n)) 
      allocate (vrand(3,n))
c
c     find half-step velocities via the Verlet recursion
c 
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            v(j,k) = v(j,k) + a(j,k)*dt_2
         end do 
      end do
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     take first A step according to the BAOAB sequence
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
         if (use_rattle)  call rattle (dtr,xold,yold,zold)
         if (use_rattle)  call rattle2 (dtr)   
      end do 
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     update velocities with frictional and random components
c
      call oprep (dt,vfric,vrand)
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3 
            v(j,k) = v(j,k)*vfric(k) + vrand(j,k)
         end do
      end do
      if (use_rattle)  call rattle2 (dt)
c
c     take second A step according to the BAOAB sequence
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
         if (use_rattle)  call rattle (dtr,xold,yold,zold)
         if (use_rattle)  call rattle2 (dtr)   
      end do 
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
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (derivs)
      deallocate (vfric) 
      deallocate (vrand)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     compute and control the temperature and pressure
c
      call kinetic (eksum,ekin,temp)
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
c     "oprep" sets up the frictional and random terms needed to
c     update positions and velocities for the BAOAB integrator
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
      egdt = exp(-(friction * dt))
      do i = 1, nuse
         k = iuse(i)
         vfric(k) = egdt 
         ktm = boltzmann * kelvin / mass(k) 
         vsig = sqrt(ktm*(1.0-egdt*egdt))
         do j = 1, 3 
            vnorm = normal () 
            vrand(j,k) = vsig * vnorm
         end do 
      end do
      return
      end
