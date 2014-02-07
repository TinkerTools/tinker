c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2011  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine stochisok                                   ##
c     ##                                                         ##
c     #############################################################
c
c
c
      subroutine stochisok (istep,dt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'freeze.i'
      include 'moldyn.i'
      include 'units.i'
      include 'virial.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bond.i'
      include 'bound.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'mdstuf.i'
      include 'potent.i'
      include 'solute.i'
      include 'stodyn.i'
      include 'usage.i'

      integer i,j,k
      integer istep
      integer nalt
      real*8 dt,dt_2,tempstor
      real*8 dta,dta_2,kBT,LkT
      real*8 epot,etot
      real*8 eksum,eps
      real*8 temp,pres
      real*8 ealt,dalt
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 viralt(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
c
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
      kBT = boltzmann*kelvin
      LkT = len_nhc*kBT
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (derivs(3,n))

c
c     save atom positions
c
      do i = 1, n
         if (use(i)) then
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
         end if
      end do


c
c    {v1-v2} evolution and Nose-like isokinetic term
c
      call applyisok (dt)
c
c    Random noise
c
      call v2random (dt)
c
c    Force-like isokinetic term
c
      call applyisok2 (dt)
c
c    Evolve positions
c

      if (use_rattle) then
         do i = 1, n
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
            x(i) = x(i) + v(1,i)*dt
            y(i) = y(i) + v(2,i)*dt
            z(i) = z(i) + v(3,i)*dt
         end do
         call rattle (dta,xold,yold,zold)
      else
         do i = 1, n
            if (use(i)) then
               x(i) = x(i) + v(1,i)*dt
               y(i) = y(i) + v(2,i)*dt
               z(i) = z(i) + v(3,i)*dt
            end if
         end do
      end if

c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)

c
c     use Newton's second law to get the next accelerations    
c
      do i = 1, n
         do j = 1, 3
            a(j,i) = -convert * derivs(j,i) / mass(i)
         end do
      end do
c
c    Force-like isokinetic term
c
      call applyisok2 (dt)
c
c    {v1-v2} evolution and Nose-like isokinetic term
c
      call applyisok (dt)

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

      call kinetic (eksum,ekin)
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
      etot = eksum + epot  
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot)
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
c     "applyisok" represents the Liuville operator iL_N in Leimkuhler, Margul, and Tuckerman (2013)
c
c
      subroutine applyisok (dt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'freeze.i'
      include 'moldyn.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer nc,ns,i,j,k,iatom,inhc,ixyz
      real*8 dt,dtc,dts
      real*8 dt2,dt4,dt8 
      real*8 tempstor,kBT,LkT
      real*8 w(3)
      real*8 len_fac

      len_fac = dble(len_nhc)/dble(len_nhc+1)

      kBT = boltzmann*kelvin
      LkT = len_nhc*kBT


      nc = 5
      ns = 3
      dtc = dt / dble(nc)
      w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
      w(2) = 1.0d0 - 2.0d0*w(1)
      w(3) = w(1)
    
c
c     use multiple time steps and Suzuki-Yoshida decomposition
c
      do k = 1, nc
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
c     "applyisok2" represents the Liuville operator iL_v in Leimkuhler, Margul, and Tuckerman (2013)
c
c
      subroutine applyisok2 (dt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'freeze.i'
      include 'moldyn.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      real*8 kBT,LkT,isoA,isoB
      real*8 dt,dt_2,rootB,arg
      real*8 iso_s,iso_sdot
      real*8 epot,etot
      integer iatom,inhc,ixyz,i,j,k

      dt_2 = 0.5d0*dt
      kBT = boltzmann*kelvin
      LkT = len_nhc*kBT


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
c     "v2random" represents the Liuville operator iL_OU in Leimkuhler, Margul, and Tuckerman (2013)
c
c
      subroutine v2random (dt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'freeze.i'
      include 'moldyn.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      real*8 dt,dt_2,kBT
      real*8 gamma,stoch_e,sigma
      real*8 normal
      integer iatom,inhc,ixyz

      kBT = boltzmann*kelvin
      gamma = 1.0 / (200.0*dt)
      sigma = sqrt(kBT* (1.0-exp(-2.0*gamma*dt) )/q_iso2(1,1,1) )
      stoch_e = exp(-gamma*dt)

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
