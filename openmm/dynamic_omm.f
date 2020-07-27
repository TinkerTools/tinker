c
c
c     ############################################################
c     ##                  COPYRIGHT (C) 2015                    ##
c     ##     by Mark Friedrichs, Lee-Ping Wang & Jay Ponder     ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program dynamic_omm  --  molecular dynamics via OpenMM API  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "dynamic_omm" computes a molecular or stochastic dynamics
c     trajectory via an interface to the OpenMM GPU code for the
c     computation of forces and dynamics integration steps
c
c
      program dynamic_omm
      use atoms
      use bath
      use bndstr
      use bound
      use boxes
      use inform
      use iounit
      use keys
      use mdstuf
      use openmm
      use openmp
      use potent
      use stodyn
      use usage
      implicit none
      integer i,istep,nstep
      integer mode,next
      integer nextStep
      integer nextUpdate
      integer updateCalls
      integer callMdStat
      integer callMdSave
      real*8 e,dt
      real*8 dtsave,speed
      real*8 elapsed,cpu
      real*8, allocatable :: derivs(:,:)
      logical exist
      logical updateEachStep
      character*1 answer
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set up the structure and molecular mechanics calculation
c
      call initial
      call getxyz
      call mechanic
      call kopenmm
c
c     multipole and polarization are computed together in OpenMM
c
      if (use_mpole .and. .not.use_polar .or.
     &       .not.use_mpole .and. use_polar) then
         use_mpole = .true.
         use_polar = .true.
         call kmpole
         call kpolar
         call mutate
      end if
c
c     initialize the temperature, pressure, integrator and GPU ID
c
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.
      integrate = 'VERLET'
c
c     check for keywords containing any altered parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         end if
      end do
c
c     initialize the simulation length as number of time steps
c
      nstep = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  nstep
   10 continue
      dowhile (nstep .lt. 0)
         write (iout,20)
   20    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,30,err=40)  nstep
   30    format (i10)
         if (nstep .lt. 0)  nstep = 0
   40    continue
      end do
c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  dt
   50 continue
      do while (dt .lt. 0.0d0)
         write (iout,60)
   60    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,70,err=80)  dt
   70    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
   80    continue
      end do
      dt = 0.001d0 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate saves
c
      dtsave = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  dtsave
   90 continue
      do while (dtsave .lt. 0.0d0)
         write (iout,100)
  100    format (/,' Enter Time between saves in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,110,err=120)  dtsave
  110    format (f20.0)
         if (dtsave .le. 0.0d0)  dtsave = 0.1d0
  120    continue
      end do
      iwrite = nint(dtsave/dt)
c
c     get choice of statistical ensemble for periodic system
c
      if (use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=130,end=130)  mode
  130    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write (iout,140)
  140       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,150,err=160)  mode
  150       format (i10)
            if (mode .le. 0)  mode = 1
  160       continue
         end do
         if (integrate.eq.'BUSSI' .or. integrate.eq.'NOSE-HOOVER'
     &                .or. integrate.eq.'GHMC') then
            if (mode .ne. 4) then
               mode = 4
               write (iout,170)
  170          format (/,' Switching to NPT Ensemble as Required',
     &                    ' by Chosen Integrator')
            end if
         end if
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=180,end=180)  kelvin
  180       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,190)
  190          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,200,err=210)  kelvin
  200          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  210          continue
            end do
         end if
         if (mode.eq.3 .or. mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=220,end=220)  atmsph
  220       continue
            do while (atmsph .lt. 0.0d0)
               write (iout,230)
  230          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,240,err=250)  atmsph
  240          format (f20.0)
               if (atmsph .le. 0.0d0)  atmsph = 1.0d0
  250          continue
            end do
         end if
      end if
c
c     use constant energy or temperature for nonperiodic system
c
      if (.not. use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=260,end=260)  mode
  260    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write (iout,270)
  270       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &              /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,280,err=290)  mode
  280       format (i10)
            if (mode .le. 0)  mode = 1
  290       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=300,end=300)  kelvin
  300       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,310)
  310          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,320,err=330)  kelvin
  320          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  330          continue
            end do
         end if
      end if
c
c     set the frequency with which data is returned from the GPU
c
c     if updateEachStep is true, take single steps on the GPU
c     and retrieve data (positions/velocities/energies) each step;
c     if updateEachStep is false, then take multiple time steps
c     on the GPU before sending data updates to the CPU with the
c     number of steps per update set by the value of "iwrite"
c
      updateEachStep = .false.
      call nextarg (string,exist)
      if (exist)  read (string,*,err=340,end=340)  answer
  340 continue
      if (.not. exist) then
         write (iout,350)
  350    format (/,' Return Data from the GPU at Every Time Step',
     &              ' [N] :  ',$)
         read (input,360)  record
  360    format (a240)
         next = 1
      end if
      call upcase (answer)
      if (answer .eq. 'Y')  updateEachStep = .true.
c
c     perform the setup functions needed to run dynamics
c
      call mdinit
c
c     get Tinker energy/gradient values for initial structure
c
      if (verbose) then
         allocate (derivs(3,n))
         call gradient (e,derivs)
         deallocate (derivs)
      end if
c
c     map Tinker data structures to OpenMM wrapper structures
c
      call ommdata ()
c
c     setup the required potential energy terms within OpenMM
c
      call openmm_init (ommHandle,dt)
c
c     compare the energy and gradient between Tinker and OpenMM
c
      if (verbose)  call openmm_test ()
c
c     print out a header line for the dynamics computation
c
      if (integrate .eq. 'VERLET') then
         write (iout,370)
  370    format (/,' Molecular Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'STOCHASTIC') then
         write (iout,380)
  380    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'BAOAB') then
         write (iout,390)
  390    format (/,' Constrained Stochastic Dynamics Trajectory',
     &              ' via BAOAB Algorithm')
      else if (integrate .eq. 'BUSSI') then
         write (iout,400)
  400    format (/,' Molecular Dynamics Trajectory via',
     &              ' Bussi-Parrinello NPT Algorithm')
      else if (integrate .eq. 'NOSE-HOOVER') then
         write (iout,410)
  410    format (/,' Molecular Dynamics Trajectory via',
     &              ' Nose-Hoover NPT Algorithm')
      else if (integrate .eq. 'GHMC') then
         write (iout,420)
  420    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Generalized Hybrid Monte Carlo')
      else if (integrate .eq. 'RIGIDBODY') then
         write (iout,430)
  430    format (/,' Molecular Dynamics Trajectory via',
     &              ' Rigid Body Algorithm')
      else if (integrate .eq. 'RESPA') then
         write (iout,440)
  440    format (/,' Molecular Dynamics Trajectory via',
     &              ' r-RESPA MTS Algorithm')
      else
         write (iout,450)
  450    format (/,' Molecular Dynamics Trajectory via',
     &              ' Modified Beeman Algorithm')
      end if
c
c     initialize some counters used during the MD steps
c
      istep = 0
      nextStep = 1
      updateCalls = 0
c
c     integrate equations of motion to take the MD time steps
c
      call settime
      if (updateEachStep) then
         callMdStat = 1
         callMdSave = 1
         do while (istep .lt. nstep)
            call openmm_take_steps (ommHandle,nextStep)
            istep = istep + nextStep
            updateCalls = updateCalls + 1
            call openmm_update (ommHandle,dt,istep,
     &                          callMdStat,callMdSave)
         end do
      else
         nextUpdate = iwrite
         callMdStat = 0
         callMdSave = 1
         do while (istep .lt. nstep)
            nextStep = nextUpdate - istep
            nextUpdate = nextUpdate + iwrite
            if (nextStep+istep .gt. nstep) then
               nextStep = nstep - istep
            end if
            call openmm_take_steps (ommHandle,nextStep)
            istep = istep + nextStep
            updateCalls = updateCalls + 1
            call openmm_update (ommHandle,dt,istep,
     &                          callMdStat,callMdSave)
         end do
      end if
      call gettime (elapsed,cpu)
c
c     print performance and timing information
c
      speed = 0.0d0
      if (elapsed .ne. 0.0d0)  speed = 86.4d0 * nstep * dt / elapsed
      write (iout,460)  speed,elapsed,nstep,updateCalls,
     &                  1000.0d0*dt,n,nthread
  460 format (/,' Performance:  ns/day',9x,f14.4,
     &        /,15x,'Wall Time',6x,f14.4,
     &        /,15x,'Steps',14x,i10,
     &        /,15x,'Updates',12x,i10,
     &        /,15x,'Time Step',6x,f14.4,
     &        /,15x,'Atoms',14x,i10,
     &        /,15x,'Threads',12x,i10)
c
c     perform any final tasks before program exit
c
      call openmm_cleanup (ommHandle)
      call final
      end
