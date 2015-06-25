c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2011  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ########################################################################
c     ##                                                                    ##
c     ##  subroutine stoch_respa2  --  r-RESPA with four time step levels  ##
c     ##                                                                    ##
c     ########################################################################
c
c
c  
c
      subroutine stoch_respa2 (istep,dt)
      use sizes
      use atomid
      use freeze
      use moldyn
      use units
      use virial
      use atoms
      use bath
      use bound
      use deriv
      use inform
      use iounit
      use keys
      use mdstuf
      use potent
      use solute
      use stodyn
      use usage
      implicit none      
      integer i,j,k
      integer ires_bond,ires_short,ires_tors
      integer istep
      real*8 eone,etwo,ethree,efour
      real*8 dt,dt_2
      real*8 dti,dtm
      real*8 epot,etot
      real*8 eksum,eps
      real*8 temp,pres
      real*8 ealt,dalt
      real*8 kBT,LkT,tempstor,len_fac
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 viralt(3,3)
      real*8 force_sum(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: derivs2(:,:)
      logical write_isok
     


c
c
c     set some values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
      dtm = dt / (dble(nres_short))
      dti = dtm / (dble(nres_tors*nres_bond))

      kBT = boltzmann*kelvin
      LkT = dble(len_nhc)*kBT
      len_fac=dble(len_nhc)/dble(len_nhc+1)
      write_isok = .false.


c
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
      allocate (derivs2(3,n))
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))

      if (XO_RESPA) call applyisok (dt)

c
c
c     begin inner RESPA loop
c

      do ires_short = 1, nres_short

         if (XM_RESPA) call applyisok (dtm)

         do ires_tors = 1, nres_tors
            do ires_bond = 1, nres_bond
c
c
c     nose-like isokinetic operator

               if (XI_RESPA) call applyisok (dti)

c
c
c     stochastic noise operator
c      
               call v2random (dti)      


c
c     force-like isokinetic operator
c

               call applyisok2 (dti)

            
c
c    Evolve positions
c
               do i = 1, n
                  if (use(i)) then         
                     xold(i) = x(i)
                     yold(i) = y(i)
                     zold(i) = z(i)
                     x(i) = x(i) + v(1,i)*dti
                     y(i) = y(i) + v(2,i)*dti
                     z(i) = z(i) + v(3,i)*dti
                  end if
               end do
               if (use_rattle) call rattle (dti,xold,yold,zold)



c
c     get the potential energy and atomic forces
c
               call grad1 (eone,derivs)
                  if (ires_bond.eq.nres_bond) then
                     call grad2 (etwo,derivs2)
                     do i = 1, n
                        do j = 1, 3
                           derivs(j,i) = derivs(j,i) 
     &                            + (mult_bond*derivs2(j,i))
                        end do
                     end do
                     if (ires_tors.eq.nres_tors) then
                        call grad3 (ethree,derivs2)                 
                        do i = 1, n
                           do j = 1, 3
                              derivs(j,i) = derivs(j,i) 
     &                              + (mult_tors*derivs2(j,i)) 
                           end do
                        end do                
                        if (ires_short.eq.nres_short) then
                           call grad4 (efour,derivs2)
                           do i = 1, n
                              do j = 1, 3
                                 derivs(j,i) = derivs(j,i)
     &                                 + (mult_short*derivs2(j,i))   
                              end do
                           end do           
                        end if
                     end if
                  end if   

                  do i = 1, n
                     do j = 1, 3
                        a(j,i) = -convert * derivs(j,i) / mass(i)
                     end do
                  end do

                  if (use_rattle) call rattle2 (dti)

c
c
c     force-like isokinetic operator
c
               call applyisok2 (dti)
c
c
c     nose-like isokinetic operator
c
               if (XI_RESPA) call applyisok (dti)

            end do
         end do

         if (XM_RESPA) call applyisok (dtm)

      end do

      if (XO_RESPA) call applyisok (dt)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)    
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)  
c
c     total energy is sum of kinetic and potential energies
c
      call kinetic (eksum,ekin)
      epot = eone + etwo + ethree + efour
      etot = eksum + epot
      call pressure (dt,epot,ekin,temp,pres,stress)      
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
c
c     if desired, report isokinetic constraint (actual/target) 
c     
      if (write_isok) then
        tempstor=0.0d0
        do i = 1, n
          do j = 1, 3
            tempstor=tempstor+(mass(i)*v(j,i)*v(j,i))
            do k = 1, len_nhc
                  tempstor=tempstor+len_fac*q_iso1(k,j,i)
     &                      *v_iso1(k,j,i)*v_iso1(k,j,i)      
            end do
         end do 
       end do            
      write (6,*) tempstor/(3.0d0*n*LkT)
      end if
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot,eksum)
      call mdrest (istep)

      return
      end

c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine applyisok                                        ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "applyisok" represents the Liouville operator iL_N in Leimkuhler, Margul, and Tuckerman (2013)
c
c
      subroutine applyisok (dt)
      use sizes
      use atomid
      use atoms
      use bath
      use freeze
      use moldyn
      use units
      use usage
      use virial
      implicit none
      integer ns,i,j,k,iatom,inhc,ixyz
      real*8 dt,dtc,dts
      real*8 dt2,dt4,dt8 
      real*8 tempstor,kBT,LkT
      real*8 w(3)
      real*8 len_fac

      len_fac = dble(len_nhc)/dble(len_nhc+1)

      kBT = boltzmann*kelvin
      LkT = (dble(len_nhc))*kBT    

      ns = 3
      dtc = dt / dble(respa_therm_nc)
      w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
      w(2) = 1.0d0 - 2.0d0*w(1)
      w(3) = w(1)
    
c
c     use multiple time steps and Suzuki-Yoshida decomposition
c
      do k = 1, respa_therm_nc
         do j = 1, ns

            dts = w(j) * dtc
            dt2 = 0.5d0 * dts
            dt4 = 0.25d0 * dts
            dt8 = 0.125d0 * dts
         
            do iatom = 1, n
               do inhc = 1, len_nhc
                  do ixyz = 1, 3
                     tempstor = q_iso1(inhc,ixyz,iatom) 
     &                *v_iso1(inhc,ixyz,iatom)*v_iso1(inhc,ixyz,iatom) 
     &                   - kBT
                     v_iso2(inhc,ixyz,iatom) = v_iso2(inhc,ixyz,iatom)
     &                   + dt4*tempstor/q_iso2(inhc,ixyz,iatom)
                  end do
               end do
            end do 
         
            do iatom = 1, n
               do ixyz = 1, 3
                  tempstor = mass(iatom)*v(ixyz,iatom)*v(ixyz,iatom) 
                  do inhc = 1, len_nhc
                     tempstor = tempstor + (len_fac * 
     &                   q_iso1(inhc,ixyz,iatom) * 
     &                   v_iso1(inhc,ixyz,iatom) * 
     &                   v_iso1(inhc,ixyz,iatom) * 
     &                   exp(-2.0 * dt2 * v_iso2(inhc,ixyz,iatom)))
                  end do
                  tempstor = sqrt(LkT/tempstor)
                  v(ixyz,iatom) = v(ixyz,iatom) * tempstor
                  do inhc = 1, len_nhc
                     v_iso1(inhc,ixyz,iatom)=v_iso1(inhc,ixyz,iatom) *
     &                   tempstor*exp(-dt2*v_iso2(inhc,ixyz,iatom))   
                  end do
               end do
            end do

            do iatom = 1, n
               do inhc = 1, len_nhc
                  do ixyz = 1, 3
                     tempstor = q_iso1(inhc,ixyz,iatom) 
     &                *v_iso1(inhc,ixyz,iatom)*v_iso1(inhc,ixyz,iatom)
     &                   - kBT
                     v_iso2(inhc,ixyz,iatom) = v_iso2(inhc,ixyz,iatom)
     &                   + dt4*tempstor/q_iso2(inhc,ixyz,iatom)
                  end do
               end do
            end do 

         end do
      end do
      return
      end
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine applyisok2                                       ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "applyisok2" represents the Liouville operator iL_v in Leimkuhler, Margul, and Tuckerman (2013)
c
c
      subroutine applyisok2 (dt)
      use sizes
      use atomid
      use atoms
      use bath
      use freeze
      use moldyn
      use units
      use usage
      use virial
      implicit none
      real*8 kBT,LkT,isoA,isoB
      real*8 dt,dt_2,rootB,arg
      real*8 iso_s,iso_sdot
      real*8 epot,etot
      integer iatom,inhc,ixyz,i,j,k

      dt_2 = 0.5d0*dt
      kBT = boltzmann*kelvin
      LkT = (dble(len_nhc))*kBT


      do iatom = 1, n
         do ixyz = 1, 3
            isoA = mass(iatom)*a(ixyz,iatom)*v(ixyz,iatom)/LkT
            isoB = mass(iatom)*a(ixyz,iatom)*a(ixyz,iatom)/LkT
            rootB = sqrt(isoB)
            arg = dt_2*rootB
            if (arg > 0.00001) then
               iso_s = (1.0/rootB)*sinh(arg) + 
     &             (isoA/isoB)*(cosh(arg)-1.0)
               iso_sdot = cosh(arg) + (isoA/rootB)*sinh(arg)
            else 
               iso_s = ((((isoB*isoA/24.0)*dt_2 +isoB/6.0)*dt_2 + 
     &             0.5*isoA)*dt_2 + 1.0)*dt_2
               iso_sdot = (((isoB*isoA/6.0)*dt_2 +0.5*isoB)*dt_2 + 
     &             isoA)*dt_2 + 1.0
            end if
            v(ixyz,iatom) = (v(ixyz,iatom) + (a(ixyz,iatom))*iso_s)
     &       /iso_sdot
            do inhc = 1, len_nhc
               v_iso1(inhc,ixyz,iatom)=v_iso1(inhc,ixyz,iatom)/iso_sdot
            end do
         end do
      end do


      return
      end
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine v2random                                         ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "v2random" represents the Liouville operator iL_OU in Leimkuhler, Margul, and Tuckerman (2013)
c
c
      subroutine v2random (dt)
      use sizes
      use atomid
      use atoms
      use bath
      use freeze
      use moldyn
      use units
      use usage
      use virial
      implicit none
      real*8 dt,dt_2,kBT
      real*8 stoch_e,sigma,OU_term
      real*8 normal
      integer iatom,inhc,ixyz

      kBT = boltzmann*kelvin
      sigma=sqrt(kBT*(1.0d0-exp(-2.0d0*stoch_gamma*dt))/q_iso2(1,1,1))
      stoch_e = exp(-stoch_gamma*dt)

      do iatom = 1, n
         do ixyz = 1, 3
            do inhc = 1, len_nhc
               v_iso2(inhc,ixyz,iatom) = v_iso2(inhc,ixyz,iatom) *
     &             stoch_e
               v_iso2(inhc,ixyz,iatom) = v_iso2(inhc,ixyz,iatom) +
     &             (sigma * normal () )
            end do
         end do
      end do

      return
      end


c
c
c     #########################################################################
c     ##                                                                     ##
c     ##  subroutine grad1  -                                                ##
c     ##                                                                     ##
c     #########################################################################
c
c
c
c
      subroutine grad1 (eone,derivs)
      use limits
      use mdstuf
      use potent
      implicit none
      real*8 eone
      real*8 derivs(3,*)
      logical save_vdw,save_charge
      logical save_chgdpl,save_dipole
      logical save_mpole,save_polar
      logical save_rxnfld,save_solv
      logical save_list
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_tortor,save_geom
      logical save_metal,save_extra
      logical save_vdwSR,save_vdwLR
      logical save_chgSR,save_chgLR
      logical save_mpoleSR,save_mpoleLR
c
c
c     save the original state of slow-evolving potentials
c
      save_vdw = use_vdw
      save_vdwSR = use_vdwSR
      save_vdwLR = use_vdwLR
      save_chgSR = use_chgSR
      save_chgLR = use_chgLR
      save_mpoleSR = use_mpoleSR
      save_mpoleLR = use_mpoleLR      
      save_charge = use_charge
      save_chgdpl = use_chgdpl
      save_dipole = use_dipole
      save_mpole = use_mpole
      save_polar = use_polar
      save_rxnfld = use_rxnfld
      save_solv = use_solv
      save_list = use_list
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
      use_vdwLR = .false.
      use_vdwSR = .false.
      use_chgLR = .false.
      use_chgSR = .false.
      use_mpoleLR = .false.
      use_mpoleSR = .false.
      use_charge = .false.
      use_chgdpl = .false.
      use_dipole = .false.
      use_mpole = .false.
      use_polar = .false.
      use_rxnfld = .false.
      use_solv = .false.
      use_list = .false.
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
      call grad_respa (eone,derivs)
         
c
c     restore the original state of slow-evolving potentials
c
      use_vdw = save_vdw
      use_vdwSR = save_vdwSR
      use_vdwLR = save_vdwLR
      use_chgSR = save_chgSR
      use_chgLR = save_chgLR  
      use_mpoleSR = save_mpoleSR
      use_mpoleLR = save_mpoleLR    
      use_charge = save_charge
      use_chgdpl = save_chgdpl
      use_dipole = save_dipole
      use_mpole = save_mpole
      use_polar = save_polar
      use_rxnfld = save_rxnfld
      use_solv = save_solv
      use_list = save_list      
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
c
c
c
c
c     #########################################################################
c     ##                                                                     ##
c     ##  subroutine grad2  -                                                ##
c     ##                                                                     ##
c     #########################################################################
c
c
c
c
      subroutine grad2 (etwo,derivs)
      use limits
      use mdstuf
      use potent
      implicit none
      real*8 etwo
      real*8 derivs(3,*)
      logical save_vdw,save_charge
      logical save_chgdpl,save_dipole
      logical save_mpole,save_polar
      logical save_rxnfld,save_solv
      logical save_list
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_tortor,save_geom
      logical save_metal,save_extra
      logical save_vdwSR,save_vdwLR
      logical save_chgSR,save_chgLR
      logical save_mpoleSR,save_mpoleLR

c
c
c     save the original state of slow-evolving potentials
c
      save_vdw = use_vdw
      save_vdwSR = use_vdwSR
      save_vdwLR = use_vdwLR
      save_chgSR = use_chgSR
      save_chgLR = use_chgLR
      save_mpoleSR = use_mpoleSR
      save_mpoleLR = use_mpoleLR
      save_charge = use_charge
      save_chgdpl = use_chgdpl
      save_dipole = use_dipole
      save_mpole = use_mpole
      save_polar = use_polar
      save_rxnfld = use_rxnfld
      save_solv = use_solv
      save_list = use_list
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


c
c     turn off fast-evolving bonded potential energy terms
c
      use_bond = .false.
      use_angle = .false.
      use_strbnd = .false.
      use_urey = .false.
      use_angang = .false.
      use_opbend = .false.
      use_opdist = .false.
      use_vdw = .false.      
      use_vdwLR = .false.
      use_vdwSR = .false.
      use_chgLR = .false.
      use_chgSR = .false.
      use_mpoleLR = .false.
      use_mpoleSR = .false.
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
      call grad_respa (etwo,derivs)
         
c
c     restore the original state of slow-evolving potentials
c
      use_vdw = save_vdw
      use_vdwSR = save_vdwSR
      use_vdwLR = save_vdwLR
      use_chgSR = save_chgSR
      use_chgLR = save_chgLR
      use_mpoleSR = save_mpoleSR
      use_mpoleLR = save_mpoleLR      
      use_charge = save_charge
      use_chgdpl = save_chgdpl
      use_dipole = save_dipole
      use_mpole = save_mpole
      use_polar = save_polar
      use_rxnfld = save_rxnfld
      use_solv = save_solv
      use_list = save_list      
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
c
c 
c
c
c
c     #########################################################################
c     ##                                                                     ##
c     ##  subroutine grad3  -                                                ##
c     ##                                                                     ##
c     #########################################################################
c
c
c
c
      subroutine grad3 (ethree,derivs)
      use limits
      use mdstuf
      use potent
      implicit none
      real*8 ethree
      real*8 derivs(3,*)
      integer ires_short
      logical save_vdw,save_charge
      logical save_chgdpl,save_dipole
      logical save_mpole,save_polar
      logical save_rxnfld,save_solv
      logical save_list
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_tortor,save_geom
      logical save_metal,save_extra
      logical save_vdwSR,save_vdwLR
      logical save_chgSR,save_chgLR
      logical save_mpoleSR,save_mpoleLR

c
c
c     save the original state of slow-evolving potentials
c
      save_vdw = use_vdw
      save_vdwSR = use_vdwSR
      save_vdwLR = use_vdwLR
      save_chgSR = use_chgSR
      save_chgLR = use_chgLR
      save_mpoleSR = use_mpoleSR
      save_mpoleLR = use_mpoleLR      
      save_charge = use_charge
      save_chgdpl = use_chgdpl
      save_dipole = use_dipole
      save_mpole = use_mpole
      save_polar = use_polar
      save_rxnfld = use_rxnfld
      save_solv = use_solv
      save_list = use_list
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
c     turn off fast-evolving bonded potential energy terms
c
      use_bond = .false.
      use_angle = .false.
      use_strbnd = .false.
      use_urey = .false.
      use_angang = .false.
      use_opbend = .false.
      use_opdist = .false.
      use_vdwLR = .false.
      use_chgLR = .false.
      use_mpoleLR = .false.
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
      call grad_respa (ethree,derivs)
         
c
c     restore the original state of slow-evolving potentials
c
      use_vdw = save_vdw
      use_vdwSR = save_vdwSR
      use_vdwLR = save_vdwLR
      use_chgSR = save_chgSR
      use_chgLR = save_chgLR 
      use_mpoleSR = save_mpoleSR
      use_mpoleLR = save_mpoleLR     
      use_charge = save_charge
      use_chgdpl = save_chgdpl
      use_dipole = save_dipole
      use_mpole = save_mpole
      use_polar = save_polar
      use_rxnfld = save_rxnfld
      use_solv = save_solv
      use_list = save_list      
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
c
c 
c
c 
c
c
c
c     #########################################################################
c     ##                                                                     ##
c     ##  subroutine grad4  -                                                ##
c     ##                                                                     ##
c     #########################################################################
c
c
c
c
      subroutine grad4 (efour,derivs)
      use limits
      use mdstuf
      use potent
      implicit none
      real*8 efour
      real*8 derivs(3,*)
      logical save_vdw,save_charge
      logical save_chgdpl,save_dipole
      logical save_mpole,save_polar
      logical save_rxnfld,save_solv
      logical save_list
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_tortor,save_geom
      logical save_metal,save_extra
      logical save_vdwSR,save_vdwLR
      logical save_chgSR,save_chgLR
      logical save_mpoleSR,save_mpoleLR

c
c
c     save the original state of slow-evolving potentials
c
      save_vdw = use_vdw
      save_vdwSR = use_vdwSR
      save_vdwLR = use_vdwLR
      save_chgSR = use_chgSR
      save_chgLR = use_chgLR    
      save_mpoleSR = use_mpoleSR
      save_mpoleLR = use_mpoleLR  
      save_charge = use_charge
      save_chgdpl = use_chgdpl
      save_dipole = use_dipole
      save_mpole = use_mpole
      save_polar = use_polar
      save_rxnfld = use_rxnfld
      save_solv = use_solv
      save_list = use_list
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
c     turn off fast-evolving bonded potential energy terms
c
      use_vdw = .false.
      use_bond = .false.
      use_angle = .false.
      use_strbnd = .false.
      use_urey = .false.
      use_angang = .false.
      use_opbend = .false.
      use_opdist = .false.
      use_vdwSR = .false.
      use_chgSR = .false.
      use_chgdpl = .false.
      use_dipole = .false.
      use_rxnfld = .false.
      use_solv = .false.
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
      call grad_respa (efour,derivs)
         
c
c     restore the original state of slow-evolving potentials
c
      use_vdw = save_vdw
      use_vdwSR = save_vdwSR
      use_vdwLR = save_vdwLR
      use_chgSR = save_chgSR
      use_chgLR = save_chgLR
      use_mpoleSR = save_mpoleSR
      use_mpoleLR = save_mpoleLR
      use_charge = save_charge
      use_chgdpl = save_chgdpl
      use_dipole = save_dipole
      use_mpole = save_mpole
      use_polar = save_polar
      use_rxnfld = save_rxnfld
      use_solv = save_solv
      use_list = save_list      
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
c
c 







