c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2011  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine isok_respa2                                 ##
c     ##                                                         ##
c     #############################################################
c
c
c
      subroutine isok_respa2 (istep,dt)
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
      integer ires_tors,ires_short,ires_bond
      integer istep
      real*8 eone,etwo,ethree,efour
      integer nalt
      real*8 dt,dt_2,tempstor,temp2
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
      real*8, allocatable :: derivs2(:,:)
      logical write_isok

c
c     set some values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
      dta = dt / (dble(nres_short*nres_tors*nres_bond))
      dta_2 = 0.5d0 * dta

      kBT = boltzmann*kelvin
      LkT = len_nhc*kBT
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (derivs(3,n))
      allocate (derivs2(3,n))

c
c     initialize virial from fast-evolving potential energy terms
c
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = 0.0d0
         end do
      end do

      if (XO_RESPA) call appisok (dt)

c
c
c     begin inner RESPA loop
c

      do ires_short = 1, nres_short
         do ires_tors = 1, nres_tors
            do ires_bond = 1, nres_bond


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
      if (XI_RESPA) call appisok (dta)

c
c
c     force-like isokinetic operator
c

            call appisok2 (dta)
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
          call grad1 (eone,derivs)
            if (ires_bond.eq.nres_bond) then
               call grad2 (etwo,derivs2)
               do i = 1, n
                  do j = 1, 3
                     derivs(j,i) = derivs(j,i) 
     &                      + (mult_bond*derivs2(j,i))
                  end do
               end do
               if (ires_tors.eq.nres_tors) then
                  call grad3 (ethree,derivs2,ires_short)                 
                  do i = 1, n
                     do j = 1, 3
                        derivs(j,i) = derivs(j,i) 
     &                        + (mult_tors*derivs2(j,i)) 
                     end do
                  end do                
                  if (ires_short.eq.nres_short) then
                     call grad4 (efour,derivs2)
                     do i = 1, n
                        do j = 1, 3
                           derivs(j,i) = derivs(j,i)
     &                           + (mult_short*derivs2(j,i))
                           
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

         if (use_rattle) call rattle2 (dta)

c
c
c     force-like isokinetic operator
c
            call appisok2 (dta)

c
c     nose-like isokinetic operator
c
      if (XI_RESPA) call appisok (dta)

            end do
         end do 
      end do

      if (XO_RESPA) call appisok (dt)

 
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     total potential and virial from sum of fast and slow parts
c


c      do i = 1, 3
c            do j = 1, 3
c                  vir(j,i) = vir(j,i) + viralt(j,i)
c            end do
c      end do


c
c     perform deallocation of some local arrays
c

      deallocate (derivs)    
      deallocate (derivs2)    
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)  
c
c     total energy is sum of kinetic and potential energies
c
      epot = eone + etwo + ethree + efour
      etot = eksum + epot

      call kinetic (eksum,ekin)
      call pressure (dt,epot,ekin,temp,pres,stress)      
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)

      if (.true.) then
         tempstor=0.0d0
         do i = 1, n
            tempstor=tempstor+(mass(i)*v(1,i)*v(1,i))
            temp2=0.5*q_iso1(1,1,i)*v_iso1(1,1,i)*v_iso1(1,1,i)
            tempstor=tempstor+temp2      
         
            tempstor=tempstor+(mass(i)*v(2,i)*v(2,i))
            temp2=0.5*q_iso1(1,2,i)*v_iso1(1,2,i)*v_iso1(1,2,i)
            tempstor=tempstor+temp2   
         
            tempstor=tempstor+(mass(i)*v(3,i)*v(3,i))
            temp2=0.5*q_iso1(1,3,i)*v_iso1(1,3,i)*v_iso1(1,3,i)
            tempstor=tempstor+temp2   
         end do
         write (6,*) tempstor/(3.0*n*kBT)
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
c     "applyisok" represents the Liuville operator iL_N in Leimkuhler, Margul, and Tuckerman (2013)
c
c
      subroutine appisok (dt)
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
      integer nc,ns,i,j,k,iatom,inhc,ixyz,index
      integer len_m1,len_p1,len_m2
      real*8 dt,dtc,dts
      real*8 dt2,dt4,dt8 
      real*8 tempstor,kBT,LkT,temp2,temp3
      real*8 w(respa_therm_nsy)

      len_m1=len_nhc-1
      len_p1=len_nhc+1
      len_m2=len_nhc-2

      kBT = boltzmann*kelvin
      LkT = len_nhc*kBT

      dtc = dt / dble(respa_therm_nc)
      if (respa_therm_nsy .eq. 3) then
         w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
         w(2) = 1.0d0 - 2.0d0*w(1)
         w(3) = w(1)
      else if (respa_therm_nsy .eq. 5) then
         w(1) = 1.0d0 / (4.0d0-4.0d0**(1.0d0/3.0d0))
         w(2) = w(1)
         w(4) = w(1)
         w(5) = w(1)
         w(3) = 1 - (4.0d0 * w(1))
      else if (respa_therm_nsy .eq. 7) then
         w(1) = 0.784513610477560
         w(2) = 0.235573213359357
         w(3) = -1.17767998417887
         w(5) = w(3)
         w(6) = w(2)
         w(7) = w(1)
         w(4) = 1.0d0 - 2.0d0*(w(1)+w(2)+w(3))
      else if (respa_therm_nsy.eq.15) then
         w(1) =  0.102799849391985
         w(2) = -1.96061023297549
         w(3) =  1.93813913762276
         w(4) = -0.158240635368243
         w(5) = -1.44485223686048
         w(6) =  0.253693336566229
         w(7) =  0.914844246229740
         w(8) = 1.0d0 - 2.0d0*(w(1)+w(2)+w(3)+w(4)+w(5)+w(6)+w(7))
         w(9) = w(1)
         w(10) = w(2)
         w(11) = w(3)
         w(12) = w(4)
         w(13) = w(5)
         w(14) = w(6)
         w(15) = w(7)
      else
         write (iout,10)
   10    format (/,' Suzuki-Yoshida decomposition of selected order
     &          not supported. Try RESPA-THERM-SY = 3,5,7, or 15')
         call fatal
      end if

    
c
c     use multiple time steps and Suzuki-Yoshida decomposition
c
      do k = 1, respa_therm_nc
         do j = 1, respa_therm_nsy

            dts = w(j) * dtc
            dt2 = 0.5d0 * dts
            dt4 = 0.25d0 * dts
            dt8 = 0.125d0 * dts

            do iatom = 1, n
              do ixyz = 1, 3
                do inhc = 2, len_nhc
                  f_iso1(inhc,ixyz,iatom)=((q_iso1(inhc,ixyz,iatom)
     &                          *v_iso1(inhc-1,ixyz,iatom)
     &                          *v_iso1(inhc-1,ixyz,iatom))-kBT)
                end do
              end do
            end do

            do iatom = 1, n
              do ixyz = 1, 3
                v_iso1(len_nhc,ixyz,iatom)=v_iso1(len_nhc,ixyz,iatom)
     &                          +dt4*f_iso1(len_nhc,ixyz,iatom)
     &                          /q_iso1(len_nhc,ixyz,iatom)
              end do
            end do


            do iatom = 1, n
              do inhc = len_m1, 2, -1
                do ixyz = 1, 3
                  index=inhc+1
                  tempstor = exp(-dt8*v_iso1(index,ixyz,iatom))
                  v_iso1(inhc,ixyz,iatom)=v_iso1(inhc,ixyz,iatom)
     &                        *tempstor*tempstor              
                  v_iso1(inhc,ixyz,iatom)=v_iso1(inhc,ixyz,iatom)
     &                        +dt4*tempstor*f_iso1(inhc,ixyz,iatom)
     &                        /q_iso1(index,ixyz,iatom)
                end do
              end do
            end do


            do iatom = 1, n
              do ixyz = 1, 3
                temp2=v(ixyz,iatom)*v(ixyz,iatom)*mass(iatom)
     &           +(0.5*q_iso1(1,ixyz,iatom)*v_iso1(1,ixyz,iatom)
     &           *v_iso1(1,ixyz,iatom)
     &           *exp(-2.0*dt2*v_iso1(2,ixyz,iatom)))
                temp3=sqrt(kBT/temp2)
                v(ixyz,iatom)=v(ixyz,iatom)*temp3
                v_iso1(1,ixyz,iatom)=v_iso1(1,ixyz,iatom)
     &           *temp3*exp(-dt2*v_iso1(2,ixyz,iatom))
              end do
            end do

            do iatom = 1, n
              do ixyz = 1, 3
                do inhc = 2, len_nhc
                  f_iso1(inhc,ixyz,iatom)=((q_iso1(inhc,ixyz,iatom)
     &                          *v_iso1(inhc-1,ixyz,iatom)
     &                          *v_iso1(inhc-1,ixyz,iatom))-kBT)
                end do
              end do
            end do

            do iatom = 1, n
              do inhc = 2, len_m1
                do ixyz = 1, 3
                  index=inhc+1
                  tempstor = exp(-dt8*v_iso1(index,ixyz,iatom))
                  v_iso1(inhc,ixyz,iatom)=v_iso1(inhc,ixyz,iatom)
     &                        *tempstor*tempstor              
                  v_iso1(inhc,ixyz,iatom)=v_iso1(inhc,ixyz,iatom)
     &                        +dt4*tempstor*f_iso1(inhc,ixyz,iatom)
     &                        /q_iso1(index,ixyz,iatom)
                end do
              end do
            end do

            do iatom = 1, n
              do ixyz = 1, 3
                v_iso1(len_nhc,ixyz,iatom)=v_iso1(len_nhc,ixyz,iatom)
     &                          +dt4*f_iso1(len_nhc,ixyz,iatom)
     &                          /q_iso1(len_nhc,ixyz,iatom)
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
c
      subroutine appisok2 (dt)
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
            isoA = mass(iatom)*a(ixyz,iatom)*v(ixyz,iatom)/kBT
            isoB = mass(iatom)*a(ixyz,iatom)*a(ixyz,iatom)/kBT
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

            v(ixyz,iatom) = (v(ixyz,iatom) + ((a(ixyz,iatom))*iso_s))
     &       /iso_sdot
            v_iso1(1,ixyz,iatom)=v_iso1(1,ixyz,iatom)/iso_sdot
         end do
      end do
  

      return
      end
