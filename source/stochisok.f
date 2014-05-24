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
      use sizes
      use atomid
      use freeze
      use moldyn
      use units
      use virial
      use atoms
      use bath
      use bound
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
      integer istep
      real*8 dt,dt_2
      real*8 dta,dta_2
      real*8 epot,etot
      real*8 eksum,eps
      real*8 temp,pres
      real*8 ealt,dalt
      real*8 kBT,LkT,tempstor,temp2,len_fac
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 viralt(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)



c
c
c     set some values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
      dta = dt / (dble(nrespa))
      dta_2 = 0.5d0 * dta

      kBT = boltzmann*kelvin
      LkT = len_nhc*kBT
      len_fac=dble(len_nhc)/dble(len_nhc+1)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))

c
c
c     initialize virial from fast-evolving potential energy terms
c
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = 0.0d0
         end do
      end do

c
c
c     begin inner RESPA loop
c

      do k = 1, nrespa
            do i = 1, n
                  if (use(i)) then         
                        xold(i) = x(i)
                        yold(i) = y(i) 
                        zold(i) = z(i) 
                  end if
            end do

c
c
c     nose-like isokinetic operator
c
            call applyisok (dta)
c
c
c     stochastic noise operator
c      
            call v2random (dta)
c
c
c     force-like isokinetic operator
c

            call applyisok2 (dta)
c
c    Evolve positions
c
            if (use_rattle) then
                  do i = 1, n
                    xold(i) = x(i)
                    yold(i) = y(i)
                    zold(i) = z(i)
                    x(i) = x(i) + v(1,i)*dta
                    y(i) = y(i) + v(2,i)*dta
                    z(i) = z(i) + v(3,i)*dta
                  end do
                  call rattle (dta,xold,yold,zold)
            else
                  do i = 1, n
                    if (use(i)) then
                        x(i) = x(i) + v(1,i)*dta
                        y(i) = y(i) + v(2,i)*dta
                        z(i) = z(i) + v(3,i)*dta
                    end if
                  end do
            end if
c
c     get the potential energy and atomic forces
c
            if (k.eq.nrespa) then
               call gradfast (ealt,derivs)

               if (use_rattle) then
                  do i = 1, n
                     if (use(i)) then
                        do j = 1, 3
                           a(j,i) = -convert * derivs(j,i) / mass(i)
                        end do
                     end if
                  end do
                  call gradslow (epot,derivs)
                  do i = 1, n
                     do j = 1, 3
                        a(j,i) = a(j,i) - ( nrespa * convert 
     &                        * derivs(j,i) / mass(i) )
                     end do
                  end do
                  call rattle2 (dta)
               else
                  do i = 1, n
                     if (use(i)) then
                        do j = 1, 3
                           a(j,i) = -convert * derivs(j,i) / mass(i)
                        end do
                     end if
                  end do
                  call gradslow (epot,derivs)
                  do i = 1, n
                     do j = 1, 3
                        a(j,i) = a(j,i) - ( nrespa * convert 
     &                        * derivs(j,i) / mass(i) )
                     end do
                  end do                  
               end if
            
            else 
               call gradfast (ealt,derivs)
               if (use_rattle) then
                  do i = 1, n
                     if (use(i)) then
                        do j = 1, 3
                           a(j,i) = -convert * derivs(j,i) / mass(i)
                        end do
                     end if
                  end do
                  call rattle2 (dta)
               else
                  do i = 1, n
                     if (use(i)) then
                        do j = 1, 3
                           a(j,i) = -convert * derivs(j,i) / mass(i)
                        end do
                     end if
                  end do
               end if               
            end if

c
c
c     force-like isokinetic operator
c
            call applyisok2 (dta)
c
c
c     nose-like isokinetic operator
c
            call applyisok (dta)

      end do
 
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     total potential and virial from sum of fast and slow parts
c
      epot = epot + ealt
      do i = 1, 3
            do j = 1, 3
                  vir(j,i) = vir(j,i) + viralt(j,i)
            end do
      end do


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
      etot = eksum + epot

      call kinetic (eksum,ekin)
      call pressure (dt,epot,ekin,temp,pres,stress)      
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)


c      tempstor=0.0d0
c       do i = 1, n
c         do j = 1, 3
c            tempstor=tempstor+(mass(i)*v(j,i)*v(j,i))
c            do k = 1, len_nhc
c                  tempstor=tempstor+len_fac*q_iso1(k,j,i)
c     &                      *v_iso1(k,j,i)*v_iso1(k,j,i)      
c            end do
c         end do 
c       end do            
c
c     print (actual/target) isokinetic constraint
c
c      write (iout,10) tempstor/(3.0d0*n*LkT)
c 10    format (3f16.14) 



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
      integer nc,ns,i,j,k,iatom,inhc,ixyz
      real*8 dt,dtc,dts
      real*8 dt2,dt4,dt8 
      real*8 tempstor,kBT,LkT
      real*8 w(3)
      real*8 len_fac

      len_fac = dble(len_nhc)/dble(len_nhc+1)

      kBT = boltzmann*kelvin
      LkT = (dble(len_nhc))*kBT    

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
c      do i = 1, 1
c         do j = 1, 1
c            write (*,10) v(j,i)
c            do k = 1, len_nhc
c               write (*,10) v_iso1(k,j,i)
c            end do
c         end do
c      end do
10    format (3f16.6)
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
c            write (*,10) isoA
c            write (*,10) isoB
c            write (*,10) LkT
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
      real*8 gamma,stoch_e,sigma
      real*8 normal
      integer iatom,inhc,ixyz

      kBT = boltzmann*kelvin
      gamma = 1.0d0 / (200.0d0*dt)
      sigma = sqrt(kBT* (1.0d0-exp(-2.0d0*gamma*dt) )/q_iso2(1,1,1) )
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