c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program dynamic  --  run molecular or stochastic dynamics  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dynamic" computes a molecular or stochastic dynamics trajectory
c     in one of the standard statistical mechanical ensembles and using
c     any of several possible integration methods
c
c
      program dynamic
      use sizes
      use atoms
      use bath
      use bndstr
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
      integer i,istep,nstep
      integer mode,next
      real*8 dt,dtdump
      logical exist
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
c
c     initialize integrator and temperature/pressure coupling
c
      integrate = 'BEEMAN'
      kelvin = -1.0d0
      atmsph = -1.0d0
      isothermal = .false.
      isobaric = .false.
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
         else if (keyword(1:12) .eq. 'TEMPERATURE ') then
            read (string,*,err=10,end=10)  kelvin
         else if (keyword(1:9) .eq. 'PRESSURE ') then
            read (string,*,err=10,end=10)  atmsph
         end if
   10    continue
      end do
c
c     initialize the simulation length as number of time steps
c
      nstep = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=20,end=20)  nstep
   20 continue
      dowhile (nstep .lt. 0)
         write (iout,30)
   30    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,40,err=50)  nstep
   40    format (i10)
         if (nstep .lt. 0)  nstep = 0
   50    continue
      end do
c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  dt
   60 continue
      do while (dt .lt. 0.0d0)
         write (iout,70)
   70    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,80,err=90)  dt
   80    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
   90    continue
      end do
      dt = 0.001d0 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=100,end=100)  dtdump
  100 continue
      do while (dtdump .lt. 0.0d0)
         write (iout,110)
  110    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,120,err=130)  dtdump
  120    format (f20.0)
         if (dtdump .le. 0.0d0)  dtdump = 0.1d0
  130    continue
      end do
      iwrite = nint(dtdump/dt)
c
c     get choice of statistical ensemble for periodic system
c
      if (use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=140,end=140)  mode
  140    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write (iout,150)
  150       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,160,err=170)  mode
  160       format (i10)
            if (mode .le. 0)  mode = 1
  170       continue
         end do
         if (integrate.eq.'BUSSI' .or. integrate.eq.'NOSE-HOOVER'
     &                .or. integrate.eq.'GHMC') then
            if (mode .ne. 4) then
               mode = 4
               write (iout,180)
  180          format (/,' Switching to NPT Ensemble as Required',
     &                    ' by Chosen Integrator')
            end if
         end if
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            if (kelvin .lt. 0.0d0) then
               call nextarg (string,exist)
               if (exist)  read (string,*,err=190,end=190)  kelvin
  190          continue
               do while (kelvin .lt. 0.0d0)
                  write (iout,200)
  200             format (/,' Enter the Desired Temperature in Degrees',
     &                       ' K [298] :  ',$)
                  read (input,210,err=220)  kelvin
  210             format (f20.0)
                  if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  220             continue
               end do
            end if
         end if
         if (mode.eq.3 .or. mode.eq.4) then
            isobaric = .true.
            if (atmsph .lt. 0.0d0) then
               call nextarg (string,exist)
               if (exist)  read (string,*,err=230,end=230)  atmsph
  230          continue
               do while (atmsph .lt. 0.0d0)
                  write (iout,240)
  240             format (/,' Enter the Desired Pressure in Atmospheres'
     &                       ' [1.0] :  ',$)
                  read (input,250,err=260)  atmsph
  250             format (f20.0)
                  if (atmsph .le. 0.0d0)  atmsph = 1.0d0
  260             continue
               end do
            end if
         end if
      end if
c
c     use constant energy or temperature for nonperiodic system
c
      if (.not. use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=270,end=270)  mode
  270    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write (iout,280)
  280       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &              /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,290,err=300)  mode
  290       format (i10)
            if (mode .le. 0)  mode = 1
  300       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            if (kelvin .lt. 0.0d0) then
               call nextarg (string,exist)
               if (exist)  read (string,*,err=310,end=310)  kelvin
  310          continue
               do while (kelvin .lt. 0.0d0)
                  write (iout,320)
  320             format (/,' Enter the Desired Temperature in Degrees',
     &                       ' K [298] :  ',$)
                  read (input,330,err=340)  kelvin
  330             format (f20.0)
                  if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  340             continue
               end do
            end if
         end if
      end if
c
c     initialize any holonomic constraints and setup dynamics
c
      call shakeup
      call mdinit
c
c     print out a header line for the dynamics computation
c
      if (integrate .eq. 'VERLET') then
         write (iout,350)
  350    format (/,' Molecular Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'STOCHASTIC') then
         write (iout,360)
  360    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'BUSSI') then
         write (iout,370)
  370    format (/,' Molecular Dynamics Trajectory via',
     &              ' Bussi-Parrinello NPT Algorithm')
      else if (integrate .eq. 'NOSE-HOOVER') then
         write (iout,380)
  380    format (/,' Molecular Dynamics Trajectory via',
     &              ' Nose-Hoover NPT Algorithm')
      else if (integrate .eq. 'GHMC') then
         write (iout,390)
  390    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Generalized Hybrid Monte Carlo')
      else if (integrate .eq. 'RIGIDBODY') then
         write (iout,400)
  400    format (/,' Molecular Dynamics Trajectory via',
     &              ' Rigid Body Algorithm')
      else if (integrate .eq. 'RESPA') then
         write (iout,410)
  410    format (/,' Molecular Dynamics Trajectory via',
     &              ' r-RESPA MTS Algorithm')
      else
         write (iout,420)
  420    format (/,' Molecular Dynamics Trajectory via',
     &              ' Modified Beeman Algorithm')
      end if
      flush (iout)
c
c     integrate equations of motion to take a time step
c
      do istep = 1, nstep
         if (integrate .eq. 'VERLET') then
            call verlet (istep,dt)
         else if (integrate .eq. 'STOCHASTIC') then
            call sdstep (istep,dt)
         else if (integrate .eq. 'BUSSI') then
            call bussi (istep,dt)
         else if (integrate .eq. 'NOSE-HOOVER') then
            call nose (istep,dt)
         else if (integrate .eq. 'GHMC') then
            call ghmcstep (istep,dt)
         else if (integrate .eq. 'RIGIDBODY') then
            call rgdstep (istep,dt)
         else if (integrate .eq. 'RESPA') then
            call respa (istep,dt)
         else
            call beeman (istep,dt)
         end if
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
