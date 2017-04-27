!ALBAUGH
      subroutine auxrespa (istep,dt)
      use sizes
      use atomid
      use atoms
      use freeze
      use ielscf
      use moldyn
      use polar
      use units
      use usage
      use virial
      use potent
      implicit none
      integer i,j,k,l
      integer istep
      integer nalt
      integer nmed
      real*8 dt,dt_2
      real*8 dta,dta_2
      real*8 dtm,dtm_2
      real*8 epot,etot
      real*8 eksum,eps
      real*8 temp,pres
      real*8 ealt,dalt
      real*8 dmed
      real*8 term
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 viralt(3,3)
      real*8 ep,e1,e2,e3
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: derivs1(:,:)
      real*8, allocatable :: derivs2(:,:)
      real*8, allocatable :: derivs3(:,:)
c
c
c     set some time values for the dynamics integration
c
      eps =  0.00000001d0
      dmed = 0.001d0
      nmed = int(dt/(dmed+eps)) + 1
      dmed = dble(nmed)
      dt_2 = 0.5d0 * dt
      dtm = dt / dmed
      dtm_2 = 0.5d0 * dtm
      dalt = 0.00025d0
      nalt = int(dtm/(dalt+eps)) + 1
      dalt = dble(nalt)
      dta = dtm / dalt
      dta_2 = 0.5d0 * dta
      
      allocate(derivs(3,n))
      allocate(derivs1(3,n))
      allocate(derivs2(3,n))
      allocate(derivs3(3,n))
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
c
c     store the current atom positions, then find half-step
c     velocities via velocity Verlet recursion
c
      
      call temper2 (dt,temp)
      
      do l = 1, nmed
         do k = 1, nalt
            do i = 1, n
               if (use(i)) then
                  do j = 1, 3
                     vaux(j,i) = vaux(j,i) + aaux(j,i)*dta_2
                     vpaux(j,i) = vpaux(j,i) + apaux(j,i)*dta_2
                     v(j,i) = v(j,i) + a(j,i)*dta_2
                     uaux(j,i) = uaux(j,i) + vaux(j,i)*dta
                     upaux(j,i) = upaux(j,i) + vpaux(j,i)*dta
                  end do
                  xold(i) = x(i)
                  yold(i) = y(i)
                  zold(i) = z(i)
                  x(i) = x(i) + v(1,i)*dta
                  y(i) = y(i) + v(2,i)*dta
                  z(i) = z(i) + v(3,i)*dta
               end if
            end do
            if (use_rattle)  call rattle (dta,xold,yold,zold)
            call gradfastaux(e1,derivs1)
            do i = 1, n
               if (use(i)) then
                  do j = 1, 3
                     derivs(j,i) = derivs1(j,i)
                     aaux(j,i) = 0.0d0
                     apaux(j,i) = 0.0d0
                  end do
               end if
            end do
            
            do i = 1, 3
               do j = 1, 3
                  viralt(j,i) = viralt(j,i) + 
     &                          vir(j,i)/(dble(nalt)*dble(nmed))
               end do
            end do
            
            if (k .eq. nalt) then
               call gradintaux(e2,derivs2)
               term = 2.0d0 / (dtm*dtm)
               if (use_ielscf) then
                  do i = 1, n
                     if (use(i)) then
                        do j = 1, 3
                           aaux(j,i)=term*(uind(j,i)
     &                                   -uaux(j,i))*dble(nalt)
                           apaux(j,i)=term*(uinp(j,i)
     &                                   -upaux(j,i))*dble(nalt)
                           derivs(j,i)=derivs(j,i)
     &                                   +dble(nalt)*derivs2(j,i)
                        end do
                     end if
                  end do
               else
                  do i = 1, n
                     if (use(i)) then
                        do j = 1, 3
                           aaux(j,i)=gamma_aux*term*(uind(j,i)
     &                                   -uaux(j,i))*dble(nalt)
                           apaux(j,i)=gamma_aux*term*(uinp(j,i)
     &                                   -upaux(j,i))*dble(nalt)
                           derivs(j,i)=derivs(j,i)
     &                                   +dble(nalt)*derivs2(j,i)
                        end do
                     end if
                  end do
               end if
               do i = 1, 3
                  do j = 1, 3
                     viralt(j,i) = viralt(j,i) + vir(j,i)/dble(nmed)
                  end do
               end do
               if(l .eq. nmed) then
                  call gradslowaux(e3,derivs3)
                  do i = 1, n
                     if (use(i)) then
                        do j = 1, 3
                           derivs(j,i)=derivs(j,i)
     &                        +derivs3(j,i)*dble(nalt)*dble(nmed)
                        end do
                     end if
                  end do
                  do i = 1, 3
                     do j = 1, 3
                        viralt(j,i) = viralt(j,i) + vir(j,i)
                     end do
                  end do
               end if
            end if
            do i = 1, n
               if (use(i)) then
                  do j = 1, 3
                     a(j,i) = -convert * derivs(j,i) / mass(i)
                     v(j,i) = v(j,i) + a(j,i)*dta_2
                     vaux(j,i) = vaux(j,i) + aaux(j,i)*dta_2
                     vpaux(j,i) = vpaux(j,i) + apaux(j,i)*dta_2
                  end do
               end if
            end do
            if (use_rattle)  call rattle2 (dta)
         end do
      end do
      
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = viralt(j,i)
         end do
      end do
      
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
      epot = e1 + e3 + e2
      etot = eksum + epot
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot,eksum)
      call mdrest (istep)
      
      deallocate (derivs)
      deallocate (derivs1)
      deallocate (derivs2)
      deallocate (derivs3)
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      return
      end
      
      


      subroutine gradfastaux (energy,derivs)
      use limits
      use potent
      implicit none
      real*8 energy
      real*8 derivs(3,*)
      logical save_vdw,save_charge
      logical save_chgdpl,save_dipole
      logical save_mpole,save_polar
      logical save_rxnfld,save_solv
      logical save_list
c
c
c     save the original state of slow-evolving potentials
c
      save_vdw = use_vdw
      save_charge = use_charge
      save_chgdpl = use_chgdpl
      save_dipole = use_dipole
      save_mpole = use_mpole
      save_polar = use_polar
      save_rxnfld = use_rxnfld
      save_solv = use_solv
      save_list = use_list
c
c     turn off slow-evolving nonbonded potential energy terms
c
      use_vdw = .false.
      use_charge = .false.
      use_chgdpl = .false.
      use_dipole = .false.
      use_mpole = .false.
      use_polar = .false.
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
      use_charge = save_charge
      use_chgdpl = save_chgdpl
      use_dipole = save_dipole
      use_mpole = save_mpole
      use_polar = save_polar
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
      subroutine gradslowaux (energy,derivs)
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
      logical save_tortor,save_geom
      logical save_metal,save_extra
      logical save_polar
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
      save_tortor = use_tortor
      save_geom = use_geom
      save_metal = use_metal
      save_extra = use_extra
      save_polar = use_polar
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
      use_tortor = .false.
      use_geom = .false.
      use_metal = .false.
      use_extra = .false.
      use_polar = .false.
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
      use_tortor = save_tortor
      use_geom = save_geom
      use_metal = save_metal
      use_extra = save_extra
      use_polar = save_polar
      return
      end
      
      
      
      
      
      
      

      subroutine gradintaux (energy,derivs)
      use limits
      use potent
      implicit none
      real*8 energy
      real*8 derivs(3,*)
      logical save_vdw,save_charge
      logical save_chgdpl,save_dipole
      logical save_mpole,save_polar
      logical save_rxnfld,save_solv
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_tortor,save_geom
      logical save_metal,save_extra
c
c
c     save the original state of slow-evolving potentials
c
      save_vdw = use_vdw
      save_charge = use_charge
      save_chgdpl = use_chgdpl
      save_dipole = use_dipole
      save_mpole = use_mpole
      save_rxnfld = use_rxnfld
      save_solv = use_solv
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
      save_tortor = use_tortor
      save_geom = use_geom
      save_metal = use_metal
      save_extra = use_extra
c
c     turn off slow-evolving nonbonded potential energy terms
c
      use_vdw = .false.
      use_charge = .false.
      use_chgdpl = .false.
      use_dipole = .false.
      use_mpole = .false.
      use_rxnfld = .false.
      use_solv = .false.
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
      use_tortor = .false.
      use_geom = .false.
      use_metal = .false.
      use_extra = .false.
c
c     get energy and gradient for fast-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of slow-evolving potentials
c
      use_vdw = save_vdw
      use_charge = save_charge
      use_chgdpl = save_chgdpl
      use_dipole = save_dipole
      use_mpole = save_mpole
      use_rxnfld = save_rxnfld
      use_solv = save_solv
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
      use_tortor = save_tortor
      use_geom = save_geom
      use_metal = save_metal
      use_extra = save_extra
      return
      end
