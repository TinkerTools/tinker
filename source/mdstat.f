c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdstat  --  compute averages over a trajectory  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mdstat" is called at each molecular dynamics time step to
c     form statistics on various average values and fluctuations,
c     and to periodically save the state of the trajectory
c
c
      subroutine mdstat (istep,dt,etot,epot,ekin,temp,pres)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cutoff.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'mdstuf.i'
      include 'molcul.i'
      include 'units.i'
      include 'usage.i'
      include 'warp.i'
      integer istep,modstep
      real*8 dt,temp,pres
      real*8 etot,epot,ekin
      real*8 pico,dens
      real*8 fluctuate,fluctuate2
      real*8 intfluct,intfluct2
      real*8 potfluct,potfluct2
      real*8 kinfluct,kinfluct2
      real*8 tfluct,pfluct,dfluct
      real*8 tfluct2,pfluct2,dfluct2
      real*8 etot_sum,etot2_sum
      real*8 eint_sum,eint2_sum
      real*8 etot_ave,etot2_ave
      real*8 eint_ave,eint2_ave
      real*8 epot_sum,epot2_sum
      real*8 ekin_sum,ekin2_sum
      real*8 epot_ave,epot2_ave
      real*8 ekin_ave,ekin2_ave
      real*8 temp_sum,temp2_sum
      real*8 temp_ave,temp2_ave
      real*8 pres_sum,pres2_sum
      real*8 pres_ave,pres2_ave
      real*8 dens_sum,dens2_sum
      real*8 dens_ave,dens2_ave
      save etot_sum,etot2_sum
      save eint_sum,eint2_sum
      save epot_sum,epot2_sum
      save ekin_sum,ekin2_sum
      save temp_sum,temp2_sum
      save pres_sum,pres2_sum
      save dens_sum,dens2_sum
c
c
c     set number of steps for block averages of properties
c
      modstep = mod(istep,iprint)
c
c     zero out summation variables for new averaging period
c
      if (modstep.eq.1 .or. iprint.eq.1) then
         etot_sum = 0.0d0
         etot2_sum = 0.0d0
         epot_sum = 0.0d0
         epot2_sum = 0.0d0
         ekin_sum = 0.0d0
         ekin2_sum = 0.0d0
         eint_sum = 0.0d0
         eint2_sum = 0.0d0
         temp_sum = 0.0d0
         temp2_sum = 0.0d0
         pres_sum = 0.0d0
         pres2_sum = 0.0d0
         dens_sum = 0.0d0
         dens2_sum = 0.0d0
      end if
c
c     print energy, temperature and pressure for current step
c
      if (verbose) then
         if (modstep .eq. 1) then
            if (use_bounds .and. integrate.ne.'STOCHASTIC') then
               write (iout,10)
   10          format (/,4x,'MD Step',6x,'E Total',3x,'E Potential',
     &                    5x,'E Kinetic',7x,'Temp',7x,'Pres',/)
            else
               write (iout,20)
   20          format (/,4x,'MD Step',6x,'E Total',3x,'E Potential',
     &                    5x,'E Kinetic',7x,'Temp',/)
            end if
         end if
         if (use_bounds .and. integrate.ne.'STOCHASTIC') then
            write (iout,30)  istep,etot,epot,ekin,temp,pres
   30       format (i10,3f14.4,2f11.2)
         else
            write (iout,40)  istep,etot,epot,ekin,temp
   40       format (i10,3f14.4,f11.2)
         end if
      end if
c
c     print header for the averages over a group of recent steps
c
      if (modstep .eq. 0) then
         pico = dble(istep) * dt
         write (iout,50)  iprint,istep
   50    format (/,' Average Values for the last',i6,' out of',
     &              i9,' Dynamics Steps')
         if (digits .ge. 10) then
            write (iout,60)  pico
 60         format (/,' Simulation Time',3x,f21.10,' Picosecond')
         else if (digits .ge. 8) then
            write (iout,61)  pico
 61         format (/,' Simulation Time',3x,f19.8,' Picosecond')
         else if (digits .ge. 6) then
            write (iout,62)  pico
 62         format (/,' Simulation Time',3x,f17.6,' Picosecond')
         else 
            write (iout,63)  pico
 63         format (/,' Simulation Time',3x,f15.4,' Picosecond')
         end if
      end if
c
c     compute total energy and fluctuation for recent steps
c
      etot_sum = etot_sum + etot
      etot2_sum = etot2_sum + etot**2
      if (modstep .eq. 0) then
         etot_ave = etot_sum / dble(iprint)
         etot2_ave = etot2_sum / dble(iprint)
         fluctuate2 = etot2_ave - etot_ave**2
         if (fluctuate2 .gt. 0.0d0) then
            fluctuate = sqrt(fluctuate2)
         else
            fluctuate = 0.0d0
         end if
         if (digits .ge. 10) then
            write (iout,70)  etot_ave,fluctuate
 70         format (' Total Energy',6x,f21.10,' Kcal/mole',3x,
     &           '(+/-',f15.10,')')
         else if (digits .ge. 8) then
            write (iout,71)  etot_ave,fluctuate
 71         format (' Total Energy',6x,f19.8,' Kcal/mole',3x,
     &           '(+/-',f13.8,')')
         else if (digits .ge. 6) then
            write (iout,72)  etot_ave,fluctuate
 72         format (' Total Energy',6x,f17.6,' Kcal/mole',3x,
     &           '(+/-',f11.6,')')
         else 
            write (iout,73)  etot_ave,fluctuate
 73         format (' Total Energy',6x,f15.4,' Kcal/mole',3x,
     &           '(+/-',f9.4,')')
         end if
      end if
c
c     compute average potential energy and its fluctuation
c
      epot_sum = epot_sum + epot
      epot2_sum = epot2_sum + epot**2
      if (modstep .eq. 0) then
         epot_ave = epot_sum / dble(iprint)
         epot2_ave = epot2_sum / dble(iprint)
         potfluct2 = epot2_ave - epot_ave**2
         if (potfluct2 .gt. 0.0d0) then
            potfluct = sqrt(potfluct2)
         else
            potfluct = 0.0d0
         end if
         if (digits .ge. 10) then
            write (iout,80)  epot_ave,potfluct
 80         format (' Potential Energy',2x,f21.10,' Kcal/mole',3x,
     &           '(+/-',f15.10,')')
         else if (digits .ge. 8) then
            write (iout,81)  epot_ave,potfluct
 81         format (' Potential Energy',2x,f19.8,' Kcal/mole',3x,
     &           '(+/-',f13.8,')')
         else if (digits .ge. 6) then
            write (iout,82)  epot_ave,potfluct
 82         format (' Potential Energy',2x,f17.6,' Kcal/mole',3x,
     &           '(+/-',f11.6,')')
         else 
            write (iout,83)  epot_ave,potfluct
 83         format (' Potential Energy',2x,f15.4,' Kcal/mole',3x,
     &           '(+/-',f9.4,')')
         end if
      end if
c
c     compute average kinetic energy and its fluctuation
c
      ekin_sum = ekin_sum + ekin
      ekin2_sum = ekin2_sum + ekin**2
      if (modstep .eq. 0) then
         ekin_ave = ekin_sum / dble(iprint)
         ekin2_ave = ekin2_sum / dble(iprint)
         kinfluct2 = ekin2_ave - ekin_ave**2
         if (kinfluct2 .gt. 0.0d0) then
            kinfluct = sqrt(kinfluct2)
         else
            kinfluct = 0.0d0
         end if
         if (digits .ge. 10) then
            write (iout,90)  ekin_ave,kinfluct
 90         format (' Kinetic Energy',4x,f21.10,' Kcal/mole',3x,
     &           '(+/-',f15.10,')')
         else if (digits .ge. 8) then
            write (iout,91)  ekin_ave,kinfluct
 91         format (' Kinetic Energy',4x,f19.8,' Kcal/mole',3x,
     &           '(+/-',f13.8,')')
         else if (digits .ge. 6) then
            write (iout,92)  ekin_ave,kinfluct
 92         format (' Kinetic Energy',4x,f17.6,' Kcal/mole',3x,
     &           '(+/-',f11.6,')')
         else 
            write (iout,93)  ekin_ave,kinfluct
 93         format (' Kinetic Energy',4x,f15.4,' Kcal/mole',3x,
     &           '(+/-',f9.4,')')
         end if
      end if
c
c     compute average intermolecular energy and its fluctuation
c
      if (nmol.ne.1 .and. nmol.ne.n .and. .not.use_ewald) then
         eint_sum = eint_sum + einter
         eint2_sum = eint2_sum + einter**2
         if (modstep .eq. 0) then
            eint_ave = eint_sum / dble(iprint)
            eint2_ave = eint2_sum / dble(iprint)
            intfluct2 = eint2_ave - eint_ave**2
            if (intfluct2 .gt. 0.0d0) then
               intfluct = sqrt(intfluct2)
            else
               intfluct = 0.0d0
            end if
            if (digits .ge. 10) then
               write (iout,100)  eint_ave,intfluct
 100           format (' Intermolecular',4x,f21.10,' Kcal/mole',3x,
     &              '(+/-',f15.10,')')
            else if (digits .ge. 8) then
               write (iout,101)  eint_ave,intfluct
 101           format (' Intermolecular',4x,f19.8,' Kcal/mole',3x,
     &              '(+/-',f13.8,')')
            else if (digits .ge. 6) then
               write (iout,102)  eint_ave,intfluct
 102           format (' Intermolecular',4x,f17.6,' Kcal/mole',3x,
     &              '(+/-',f11.6,')')
            else 
               write (iout,103)  eint_ave,intfluct
 103           format (' Intermolecular',4x,f15.4,' Kcal/mole',3x,
     &              '(+/-',f9.4,')')
            end if
         end if
      end if
c
c     compute the average temperature and its fluctuation
c
      temp_sum = temp_sum + temp
      temp2_sum = temp2_sum + temp**2
      if (modstep .eq. 0) then
         temp_ave = temp_sum / dble(iprint)
         temp2_ave = temp2_sum / dble(iprint)
         tfluct2 = temp2_ave - temp_ave**2
         if (tfluct2 .gt. 0.0d0) then
            tfluct = sqrt(tfluct2)
         else
            tfluct = 0.0d0
         end if
         if (digits .ge. 10) then
            write (iout,110)  temp_ave,tfluct
 110        format (' Temperature',7x,f21.8,' Kelvin',6x,
     &           '(+/-',f15.8,')')
         else if (digits .ge. 8) then
            write (iout,111)  temp_ave,tfluct
 111        format (' Temperature',7x,f19.6,' Kelvin',6x,
     &           '(+/-',f13.6,')')
         else if (digits .ge. 6) then
            write (iout,112)  temp_ave,tfluct
 112        format (' Temperature',7x,f17.4,' Kelvin',6x,
     &           '(+/-',f11.4,')')
         else 
            write (iout,113)  temp_ave,tfluct
 113        format (' Temperature',7x,f15.2,' Kelvin',6x,
     &           '(+/-',f9.2,')')
         end if
      end if
c
c     compute the average pressure and its fluctuation
c
      if (use_bounds) then
         pres_sum = pres_sum + pres
         pres2_sum = pres2_sum + pres**2
         if (modstep .eq. 0) then
            pres_ave = pres_sum / dble(iprint)
            pres2_ave = pres2_sum / dble(iprint)
            pfluct2 = pres2_ave - pres_ave**2
            if (pfluct2 .gt. 0.0d0) then
               pfluct = sqrt(pfluct2)
            else
               pfluct = 0.0d0
            end if
            if (digits .ge. 10) then
               write (iout,120)  pres_ave,pfluct
 120           format (' Pressure',10x,f21.8,' Atmosphere',2x,
     &              '(+/-',f15.8,')')
            else if (digits .ge. 8) then
               write (iout,121)  pres_ave,pfluct
 121           format (' Pressure',10x,f19.6,' Atmosphere',2x,
     &              '(+/-',f13.6,')')
            else if (digits .ge. 6) then
               write (iout,122)  pres_ave,pfluct
 122           format (' Pressure',10x,f17.4,' Atmosphere',2x,
     &              '(+/-',f11.4,')')
            else 
               write (iout,123)  pres_ave,pfluct
 123           format (' Pressure',10x,f15.2,' Atmosphere',2x,
     &              '(+/-',f9.2,')')
            end if
         end if
c
c     compute the average density and its fluctuation
c
         dens = (1.0d24/volbox) * (totmass/avogadro)
         dens_sum = dens_sum + dens
         dens2_sum = dens2_sum + dens**2
         if (modstep .eq. 0) then
            dens_ave = dens_sum / dble(iprint)
            dens2_ave = dens2_sum / dble(iprint)
            dfluct2 = dens2_ave - dens_ave**2
            if (dfluct2 .gt. 0.0d0) then
               dfluct = sqrt(dfluct2)
            else
               dfluct = 0.0d0
            end if
            if (digits .ge. 10) then
               write (iout,130)  dens_ave,dfluct
 130           format (' Density',11x,f21.10,' Grams/cc',4x,
     &              '(+/-',f15.10,')')
            else if (digits .ge. 8) then
               write (iout,131)  dens_ave,dfluct
 131           format (' Density',11x,f19.8,' Grams/cc',4x,
     &              '(+/-',f13.8,')')
            else if (digits .ge. 6) then
               write (iout,132)  dens_ave,dfluct
 132           format (' Density',11x,f17.6,' Grams/cc',4x,
     &              '(+/-',f11.6,')')
            else 
               write (iout,133)  dens_ave,dfluct
 133           format (' Density',11x,f15.4,' Grams/cc',4x,
     &              '(+/-',f9.4,')')
            end if
         end if
      end if
c
c     note deformation value for potential energy smoothing
c
      if (use_smooth) then
         if (modstep .eq. 0) then
            if (digits .ge. 10) then
               write (iout,140)  deform
 140           format (' Deformation',7x,f21.9,' Sqr Angs')
            else if (digits .ge. 8) then
               write (iout,141)  deform
 141           format (' Deformation',7x,f19.7,' Sqr Angs')
            else if (digits .ge. 6) then
               write (iout,142)  deform
 142           format (' Deformation',7x,f17.5,' Sqr Angs')
            else 
               write (iout,143)  deform
 143           format (' Deformation',7x,f15.3,' Sqr Angs')
            end if
         end if
      end if
      return
      end
