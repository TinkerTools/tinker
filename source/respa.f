c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2011  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine vrespa  --  Verlet r-RESPA molecular dynamics  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "vrespa" performs a multiple time step (MTS) molecular dynamics
c     step using the reversible reference system propagation algorithm
c     (r-RESPA) via a velocity Verlet recursion with the potential
c     split into fast- and slow-evolving components
c
c     literature references:
c
c     D. D. Humphreys, R. A. Friesner and B. J. Berne, "A Multiple-
c     Time-Step Molecular Dynamics Algorithm for Macromolecules",
c     Journal of Physical Chemistry, 98, 6885-6892 (1994)
c
c     X. Qian and T. Schlick, "Efficient Multiple-Time-Step Integrators
c     with Distance-Based Force Splitting for Particle-Mesh-Ewald
c     Molecular Dynamics Simulations", Journal of Chemical Physics,
c     115, 4019-4029 (2001)
c
c
      subroutine vrespa (istep,dt)
      use atomid
      use atoms
      use freeze
      use ielscf
      use mdstuf
      use moldyn
      use polar
      use units
      use usage
      use virial
      implicit none
      integer i,j,k,m
      integer istep
      real*8 dt,dt_2
      real*8 dta,dta_2
      real*8 epot,etot
      real*8 eksum,term
      real*8 temp,pres
      real*8 drespa,efast
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 virfast(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
c
c
c     set some time values for the dynamics integration
c
      drespa = dble(nrespa)
      dta = dt / drespa
      dt_2 = 0.5d0 * dt
      dta_2 = 0.5d0 * dta
c
c     find half-step velocities via velocity Verlet recursion
c
      do i = 1, nuse
         m = iuse(i)
         do j = 1, 3
            v(j,m) = v(j,m) + a(j,m)*dt_2
         end do
      end do
c
c     initialize virial from fast-evolving potential energy terms
c
      do i = 1, 3
         do j = 1, 3
            virfast(j,i) = 0.0d0
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
c     get fast-evolving velocities and positions via Verlet recursion
c
      do k = 1, nrespa
         do i = 1, nuse
            m = iuse(i)
            do j = 1, 3
               v(j,m) = v(j,m) + aalt(j,m)*dta_2
            end do
            xold(m) = x(m)
            yold(m) = y(m)
            zold(m) = z(m)
            x(m) = x(m) + v(1,m)*dta
            y(m) = y(m) + v(2,m)*dta
            z(m) = z(m) + v(3,m)*dta
         end do
         if (use_rattle)  call rattle (dta,xold,yold,zold)
c
c     find the fast-evolving potential energy and atomic forces
c
         call gradfast (efast,derivs)
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using Verlet recursion
c
         do i = 1, nuse
            m = iuse(i)
            do j = 1, 3
               aalt(j,m) = -ekcal * derivs(j,m) / mass(m)
               v(j,m) = v(j,m) + aalt(j,m)*dta_2
            end do
         end do
         if (use_rattle)  call rattle2 (dta)
c
c     find average virial from fast-evolving potential terms
c
         do i = 1, 3
            do j = 1, 3
               virfast(j,i) = virfast(j,i) + vir(j,i)/drespa
            end do
         end do
      end do
c
c     apply Verlet half-step updates for any auxiliary dipoles
c
      if (use_ielscf) then
         do i = 1, nuse
            m = iuse(i)
            do j = 1, 3
               vaux(j,m) = vaux(j,m) + aaux(j,m)*dt_2
               vpaux(j,m) = vpaux(j,m) + apaux(j,m)*dt_2
               uaux(j,m) = uaux(j,m) + vaux(j,m)*dt
               upaux(j,m) = upaux(j,m) + vpaux(j,m)*dt
            end do
         end do
      end if
c
c     get the slow-evolving potential energy and atomic forces
c
      call gradslow (epot,derivs)
      epot = epot + efast
c
c     compute and make the half-step temperature correction
c
      call temper2 (dt,temp)
c
c     use Newton's second law to get the slow accelerations;
c     find full-step velocities using Verlet recursion
c
      do i = 1, nuse
         m = iuse(i)
         do j = 1, 3
            a(j,m) = -ekcal * derivs(j,m) / mass(m)
            v(j,m) = v(j,m) + a(j,m)*dt_2
         end do
      end do
c
c     apply Verlet full-step updates for any auxiliary dipoles
c
      if (use_ielscf) then
         term = 2.0d0 / (dt*dt)
         do i = 1, nuse
            m = iuse(i)
            do j = 1, 3
               aaux(j,m) = term * (uind(j,m)-uaux(j,m))
               apaux(j,m) = term * (uinp(j,m)-upaux(j,m))
               vaux(j,m) = vaux(j,m) + aaux(j,m)*dt_2
               vpaux(j,m) = vpaux(j,m) + apaux(j,m)*dt_2
            end do
         end do
      end if
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle) then
         do i = 1, nuse
            m = iuse(i)
            xold(m) = x(m)
            yold(m) = y(m)
            zold(m) = z(m)
         end do
         call rattle2 (dt)
      end if
c
c     increment total virial from sum of fast and slow parts
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + virfast(j,i)
         end do
      end do
c
c     compute full-step temperature and pressure corrections
c
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,ekin,temp,pres,stress)
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine brespa  --  Beeman r-RESPA molecular dynamics  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "brespa" performs a multiple time step (MTS) molecular dynamics
c     step using the reversible reference system propagation algorithm
c     (r-RESPA) via a Beeman recursion with the potential split into
c     fast- and slow-evolving components
c
c     literature references:
c
c     D. D. Humphreys, R. A. Friesner and B. J. Berne, "A Multiple-
c     Time-Step Molecular Dynamics Algorithm for Macromolecules",
c     Journal of Physical Chemistry, 98, 6885-6892 (1994)
c
c     X. Qian and T. Schlick, "Efficient Multiple-Time-Step Integrators
c     with Distance-Based Force Splitting for Particle-Mesh-Ewald
c     Molecular Dynamics Simulations", Journal of Chemical Physics,
c     115, 4019-4029 (2001)
c
c     D. Beeman, "Some Multistep Methods for Use in Molecular
c     Dynamics Calculations", Journal of Computational Physics,
c     20, 130-139 (1976)
c
c
      subroutine brespa (istep,dt)
      use atomid
      use atoms
      use freeze
      use ielscf
      use mdstuf
      use moldyn
      use polar
      use units
      use usage
      use virial
      implicit none
      integer i,j,k,m
      integer istep
      real*8 dt,dt_2
      real*8 dmix,dta
      real*8 dtx,dtax
      real*8 epot,etot
      real*8 eksum,term
      real*8 temp,pres
      real*8 part1,part2
      real*8 drespa,efast
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 virfast(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
c
c
c     set some time values for the dynamics integration
c
      drespa = dble(nrespa)
      dmix = dble(bmnmix)
      part1 = 0.5d0*dmix + 1.0d0
      part2 = part1 - 2.0d0
      dtx = dt / dmix
      dta = dt / drespa
      dtax = dta / dmix
      dt_2 = 0.5d0 * dt
c
c     find half-step velocities via the Beeman recursion
c
      do i = 1, nuse
         m = iuse(i)
         do j = 1, 3
            v(j,m) = v(j,m) + (part1*a(j,m)-aslow(j,m))*dtx
         end do
      end do
c
c     initialize virial from fast-evolving potential energy terms
c
      do i = 1, 3
         do j = 1, 3
            virfast(j,i) = 0.0d0
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
c     get fast-evolving velocities and positions via Beeman recursion
c
      do k = 1, nrespa
         do i = 1, nuse
            m = iuse(i)
            do j = 1, 3
               v(j,m) = v(j,m) + (part1*aalt(j,m)-afast(j,m))*dtax
            end do
            xold(m) = x(m)
            yold(m) = y(m)
            zold(m) = z(m)
            x(m) = x(m) + v(1,m)*dta
            y(m) = y(m) + v(2,m)*dta
            z(m) = z(m) + v(3,m)*dta
         end do
         if (use_rattle)  call rattle (dta,xold,yold,zold)
c
c     find the fast-evolving potential energy and atomic forces
c
         call gradfast (efast,derivs)
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using Beeman recursion
c
         do i = 1, nuse
            m = iuse(i)
            do j = 1, 3
               afast(j,m) = aalt(j,m)
               aalt(j,m) = -ekcal * derivs(j,m) / mass(m)
               v(j,m) = v(j,m) + (part2*aalt(j,m)+afast(j,m))*dtax
            end do
         end do
         if (use_rattle)  call rattle2 (dta)
c
c     find average virial from fast-evolving potential terms
c
         do i = 1, 3
            do j = 1, 3
               virfast(j,i) = virfast(j,i) + vir(j,i)/drespa
            end do
         end do
      end do
c
c     apply Verlet half-step updates for any auxiliary dipoles
c
      if (use_ielscf) then
         do i = 1, nuse
            m = iuse(i)
            do j = 1, 3
               vaux(j,m) = vaux(j,m) + aaux(j,m)*dt_2
               vpaux(j,m) = vpaux(j,m) + apaux(j,m)*dt_2
               uaux(j,m) = uaux(j,m) + vaux(j,m)*dt
               upaux(j,m) = upaux(j,m) + vpaux(j,m)*dt
            end do
         end do
      end if
c
c     get the slow-evolving potential energy and atomic forces
c
      call gradslow (epot,derivs)
      epot = epot + efast
c
c     compute and make the half-step temperature correction
c
      call temper2 (dt,temp)
c
c     use Newton's second law to get the slow accelerations;
c     find full-step velocities using Beeman recursion
c
      do i = 1, nuse
         m = iuse(i)
         do j = 1, 3
            aslow(j,m) = a(j,m)
            a(j,m) = -ekcal * derivs(j,m) / mass(m)
            v(j,m) = v(j,m) + (part2*a(j,m)+aslow(j,m))*dtx
         end do
      end do
c
c     apply Verlet full-step updates for any auxiliary dipoles
c
      if (use_ielscf) then
         term = 2.0d0 / (dt*dt)
         do i = 1, nuse
            m = iuse(i)
            do j = 1, 3
               aaux(j,m) = term * (uind(j,m)-uaux(j,m))
               apaux(j,m) = term * (uinp(j,m)-upaux(j,m))
               vaux(j,m) = vaux(j,m) + aaux(j,m)*dt_2
               vpaux(j,m) = vpaux(j,m) + apaux(j,m)*dt_2
            end do
         end do
      end if
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle) then
         do i = 1, nuse
            m = iuse(i)
            xold(m) = x(m)
            yold(m) = y(m)
            zold(m) = z(m)
         end do
         call rattle2 (dt)
      end if
c
c     increment total virial from sum of fast and slow parts
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + virfast(j,i)
         end do
      end do
c
c     compute full-step temperature and pressure corrections
c
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,ekin,temp,pres,stress)
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine srespa  --  BAOAB r-RESPA stochastic dynamics  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "srespa" performs a multiple time step (MTS) stochastic dynamics
c     step using the reversible reference system propagation algorithm
c     (r-RESPA) via a BAOAB recursion with the potential split into
c     fast- and slow-evolving components
c
c     literature reference:
c
c     D. D. Humphreys, R. A. Friesner and B. J. Berne, "A Multiple-
c     Time-Step Molecular Dynamics Algorithm for Macromolecules",
c     Journal of Physical Chemistry, 98, 6885-6892 (1994)
c
c     B. Leimkuhler and C. Matthews, "Efficient Molecular Dynamics
c     Using Geodesic Integration and Solvent-Solute Splitting",
c     Proceedings of the Royal Society A, 472, 20160138 (2016)
c
c
      subroutine srespa (istep,dt)
      use atomid
      use atoms
      use freeze
      use mdstuf
      use moldyn
      use units
      use usage
      use virial
      implicit none
      integer i,j,k,m
      integer istep
      integer nrattle
      real*8 dt,dt_2
      real*8 dta,dta_2
      real*8 dtar,dtar_2
      real*8 etot,epot
      real*8 eksum,eksave
      real*8 temp,tave,pres
      real*8 drespa,efast
      real*8 drattle
      real*8 ekin(3,3)
      real*8 ekave(3,3)
      real*8 stress(3,3)
      real*8 virrat(3,3)
      real*8 virfast(3,3)
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
      drespa = dble(nrespa)
      nrattle = 1
      if (use_rattle)  nrattle = 3
      drattle = dble(nrattle)
      dta = dt / drespa
      dtar = dta / drattle
      dt_2 = 0.5d0 * dt
      dta_2 = 0.5d0 * dta
      dtar_2 = 0.5d0 * dtar
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
c     use outer B step to find half-step slow velocities
c
      do i = 1, nuse
         m = iuse(i)
         do j = 1, 3
            v(j,m) = v(j,m) + a(j,m)*dt_2
         end do
      end do
c
c     initialize kinetic energy and virial from fast potentials
c
      eksave = 0.0d0
      do i = 1, 3
         do j = 1, 3
            ekave(j,i) = 0.0d0
            virfast(j,i) = 0.0d0
         end do
      end do
      tave = 0.0d0
c
c     use inner B step to find half-step fast velocities
c
      do k = 1, nrespa
         do i = 1, nuse
            m = iuse(i)
            do j = 1, 3
               v(j,m) = v(j,m) + aalt(j,m)*dta_2
            end do
         end do
c
c     take an inner A step to get fast half-step positions
c
         do j = 1, nrattle
            do i = 1, nuse
               m = iuse(i)
               xold(m) = x(m)
               yold(m) = y(m)
               zold(m) = z(m)
               x(m) = x(m) + v(1,m)*dtar_2
               y(m) = y(m) + v(2,m)*dtar_2
               z(m) = z(m) + v(3,m)*dtar_2
            end do
            if (use_rattle) then
               call rattle (dtar_2,xold,yold,zold)
               call rattle2 (dtar_2)
               do i = 1, 3
                  do m = 1, 3
                     vir(m,i) = 0.0d0
                  end do
               end do
            end if
         end do
c
c     use inner O step to get frictional and random components
c
         call oprep (dta,vfric,vrand)
         do i = 1, nuse
            m = iuse(i)
            do j = 1, 3 
               v(j,m) = v(j,m)*vfric(m) + vrand(j,m)
            end do
         end do
         if (use_rattle) then
            call rattle2 (dta)
            do i = 1, 3
               do j = 1, 3
                  virrat(j,i) = vir(j,i)
                  vir(j,i) = 0.0d0
               end do
            end do
         end if
c
c     take second inner A step to get the full-step positions
c
         do j = 1, nrattle
            do i = 1, nuse
               m = iuse(i)
               xold(m) = x(m)
               yold(m) = y(m)
               zold(m) = z(m)
               x(m) = x(m) + v(1,m)*dtar_2
               y(m) = y(m) + v(2,m)*dtar_2
               z(m) = z(m) + v(3,m)*dtar_2
            end do
            if (use_rattle) then
               call rattle (dtar_2,xold,yold,zold)
               call rattle2 (dtar_2)
               do i = 1, 3
                  do m = 1, 3
                     virrat(m,i) = virrat(m,i) + vir(m,i)/drattle
                     vir(m,i) = 0.0d0
                  end do
               end do
            end if
         end do
c
c     get the fast-evolving potential energy and atomic forces
c 
         call gradfast (efast,derivs)
c
c     find average kinetic energy from fast-evolving potentials
c
         call kinetic (eksum,ekin,temp)
         eksave = eksave + eksum/drespa
         do i = 1, 3
            do j = 1, 3
               ekave(j,i) = ekave(j,i) + ekin(j,i)/drespa
            end do
         end do
         tave = tave + temp/drespa
c
c     inner B step for fast accelerations and full-step velocities
c
         do i = 1, nuse
            m = iuse(i)
            do j = 1, 3
               aalt(j,m) = -ekcal * derivs(j,m) / mass(m)
               v(j,m) = v(j,m) + aalt(j,m)*dta_2
            end do
         end do 
         if (use_rattle) then
            call rattle2 (dta)
            do i = 1, 3
               do j = 1, 3
                  vir(j,i) = vir(j,i) + virrat(j,i)
               end do
            end do
         end if
c
c     average the virial from fast-evolving potential terms
c
         do i = 1, 3
            do j = 1, 3
               virfast(j,i) = virfast(j,i) + vir(j,i)/drespa
            end do
         end do
      end do
c
c     transfer average kinetic energy and temperature values
c
      eksum = eksave
      do i = 1, 3
         do j = 1, 3
            ekin(j,i) = ekave(j,i)
         end do
      end do
      temp = tave
c
c     get the slow-evolving potential energy and atomic forces
c 
      call gradslow (epot,derivs)
      epot = epot + efast
c
c     outer B step for slow accelerations and full-step velocities
c
      do i = 1, nuse
         m = iuse(i)
         do j = 1, 3
            a(j,m) = -ekcal * derivs(j,m) / mass(m)
            v(j,m) = v(j,m) + a(j,m)*dt_2
         end do
      end do 
      if (use_rattle) then
         call rattle2 (dt)
         do i = 1, nuse
            m = iuse(i)
            xold(m) = x(m)
            yold(m) = y(m)
            zold(m) = z(m)
         end do
      end if
c
c     increment total virial from sum of fast and slow parts
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + virfast(j,i)
         end do
      end do
c
c     compute full-step kinetic energy and pressure correction;
c     prior kinetic energy gives better pressure control
c
c     call kinetic (eksum,ekin,temp)
      call pressure (dt,ekin,temp,pres,stress)
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradfast  --  fast energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradfast" calculates the potential energy and first derivatives
c     for the fast-evolving local valence potential energy terms
c
c
      subroutine gradfast (energy,derivs)
      use limits
      use potent
      implicit none
      real*8 energy
      real*8 derivs(3,*)
      logical save_vdw,save_repel
      logical save_disp,save_charge
      logical save_chgdpl,save_dipole
      logical save_mpole,save_polar
      logical save_chgtrn,save_rxnfld
      logical save_solv,save_list
c
c
c     save the original state of slow-evolving potentials
c
      save_vdw = use_vdw
      save_repel = use_repel
      save_disp = use_disp
      save_charge = use_charge
      save_chgdpl = use_chgdpl
      save_dipole = use_dipole
      save_mpole = use_mpole
      save_polar = use_polar
      save_chgtrn = use_chgtrn
      save_rxnfld = use_rxnfld
      save_solv = use_solv
      save_list = use_list
c
c     turn off slow-evolving nonbonded potential energy terms
c
      use_vdw = .false.
      use_repel = .false.
      use_disp = .false.
      use_charge = .false.
      use_chgdpl = .false.
      use_dipole = .false.
      use_mpole = .false.
      use_polar = .false.
      use_chgtrn = .false.
      use_rxnfld = .false.
      use_solv = .false.
      use_list = .false.
c
c     get energy and gradient for fast-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of slow-evolving potentials
c
      use_vdw = save_vdw
      use_repel = save_repel
      use_disp = save_disp
      use_charge = save_charge
      use_chgdpl = save_chgdpl
      use_dipole = save_dipole
      use_mpole = save_mpole
      use_polar = save_polar
      use_chgtrn = save_chgtrn
      use_rxnfld = save_rxnfld
      use_solv = save_solv
      use_list = save_list
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradslow  --  slow energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradslow" calculates the potential energy and first derivatives
c     for the slow-evolving nonbonded potential energy terms
c
c
      subroutine gradslow (energy,derivs)
      use potent
      implicit none
      real*8 energy
      real*8 derivs(3,*)
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_angtor,save_tortor
      logical save_geom,save_metal
      logical save_extra
c
c
c     save the original state of fast-evolving potentials
c
      save_bond = use_bond
      save_angle = use_angle
      save_strbnd = use_strbnd
      save_urey = use_urey
      save_angang = use_angang
      save_opbend = use_opbend
      save_opdist = use_opdist
      save_improp = use_improp
      save_imptor = use_imptor
      save_tors = use_tors
      save_pitors = use_pitors
      save_strtor = use_strtor
      save_angtor = use_angtor
      save_tortor = use_tortor
      save_geom = use_geom
      save_metal = use_metal
      save_extra = use_extra
c
c     turn off fast-evolving valence potential energy terms
c
      use_bond = .false.
      use_angle = .false.
      use_strbnd = .false.
      use_urey = .false.
      use_angang = .false.
      use_opbend = .false.
      use_opdist = .false.
      use_improp = .false.
      use_imptor = .false.
      use_tors = .false.
      use_pitors = .false.
      use_strtor = .false.
      use_angtor = .false.
      use_tortor = .false.
      use_geom = .false.
      use_metal = .false.
      use_extra = .false.
c
c     get energy and gradient for slow-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of fast-evolving potentials
c
      use_bond = save_bond
      use_angle = save_angle
      use_strbnd = save_strbnd
      use_urey = save_urey
      use_angang = save_angang
      use_opbend = save_opbend
      use_opdist = save_opdist
      use_improp = save_improp
      use_imptor = save_imptor
      use_tors = save_tors
      use_pitors = save_pitors
      use_strtor = save_strtor
      use_angtor = save_angtor
      use_tortor = save_tortor
      use_geom = save_geom
      use_metal = save_metal
      use_extra = save_extra
      return
      end
