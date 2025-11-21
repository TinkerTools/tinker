c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program mdavg  --  statistics from molecular dynamics log  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mdavg" is a simple utility to read the output from a Tinker
c     dynamics simulation and compute average and standard deviation
c     for quantities such as total energy, potential energy, kinetic
c     energy, temperature, pressure and density
c
c
      program mdavg
      use inform
      use iounit
      implicit none
      integer i,ilog
      integer nask,nread
      integer nstep,nblock
      integer maxblock
      integer start,stop
      integer freeunit
      real*8 block,time
      real*8 val,var
      real*8 etot,vtot
      real*8 epot,vpot
      real*8 ekin,vkin
      real*8 temp,vtemp
      real*8 pres,vpres
      real*8 dens,vdens
      logical exist,query
      logical doinfo,proceed
      character*240 logfile
      character*240 record
      character*240 string
c
c
c     default unit numbers for I/O and command line arguments
c
      input = 5
      iout = 6
      call command
c
c     zero out the individual average and variance values
c
      doinfo = .false.
      proceed = .false.
      maxblock = 1000000
      nread = 0
      nstep = 0
      nblock = 0
      time = 10000000.0d0
      etot = 0.0d0
      vtot = 0.0d0
      epot = 0.0d0
      vpot = 0.0d0
      ekin = 0.0d0
      vkin = 0.0d0
      temp = 0.0d0
      vtemp = 0.0d0
      pres = 0.0d0
      vpres = 0.0d0
      dens = 0.0d0
      vdens = 0.0d0
c
c     try to get a filename from the command line arguments
c
      call nextarg (logfile,exist)
      if (exist) then
         call basefile (logfile)
         call suffix (logfile,'log','old')
         inquire (file=logfile,exist=exist)
      end if
c
c     ask for the user specified dynamics log filename
c
      if (.not. exist)  call promo
      nask = 0
      do while (.not.exist .and. nask.lt.maxask)
         doinfo = .true.
         nask = nask + 1
         write (iout,10)
   10    format (/,' Enter Molecular Dynamics Log File Name :  ',$)
         read (input,20)  logfile
   20    format (a240)
         call basefile (logfile)
         call suffix (logfile,'log','old')
         inquire (file=logfile,exist=exist)
      end do
      if (.not. exist)  call fatal
c
c     find the first and last block to use in the analysis
c
      start = 1
      stop = maxblock
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=30,end=30)  start
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  stop
   30 continue
c
c     make interactive query for the range of blocks to use
c
      if (doinfo) then
         write (iout,40)
   40    format (/,' Numbers of First & Last Block of Steps :  ',$)
         read (input,50)  record
   50    format (a240)
         read (record,*,err=60,end=60)  start,stop
   60    continue
      end if
c
c     open the input molecular dynamics log file
c
      ilog = freeunit ()
      open (unit=ilog,file=logfile,status='old')
      rewind (unit=ilog)
c
c     get block average values from the dynamics log file
c
      do i = 1, maxblock
         read (ilog,70,err=80,end=80)  record
   70    format (a240)
         if (record(2:15) .eq. 'Average Values') then
            nread = nread + 1
            string = record(29:34)
            read (string,*)  nstep
            proceed = .false.
            if (nread.ge.start .and. nread.le.stop) then
               proceed = .true.
               nblock = nblock + 1
            end if
         else if (record(2:16) .eq. 'Simulation Time') then
            string = record(21:36)
            read (string,*)  val
            time = val / dble(nread)
         end if
         if (proceed) then
            if (record(2:13) .eq. 'Total Energy') then
               string = record(21:36)
               read (string,*)  val
               string = record(54:62)
               read (string,*)  var
               etot = etot + val
               vtot = vtot + var*var
            else if (record(2:17) .eq. 'Potential Energy') then
               string = record(21:36)
               read (string,*)  val
               string = record(54:62)
               read (string,*)  var
               epot = epot + val
               vpot = vpot + var*var
            else if (record(2:15) .eq. 'Kinetic Energy') then
               string = record(21:36)
               read (string,*)  val
               string = record(54:62)
               read (string,*)  var
               ekin = ekin + val
               vkin = vkin + var*var
            else if (record(2:12) .eq. 'Temperature') then
               string = record(21:36)
               read (string,*)  val
               string = record(54:62)
               read (string,*)  var
               temp = temp + val
               vtemp = vtemp + var*var
            else if (record(2:9) .eq. 'Pressure') then
               string = record(21:36)
               read (string,*)  val
               string = record(54:62)
               read (string,*)  var
               pres = pres + val
               vpres = vpres + var*var
            else if (record(2:8) .eq. 'Density') then
               string = record(21:36)
               read (string,*)  val
               string = record(54:62)
               read (string,*)  var
               dens = dens + val
               vdens = vdens + var*var
            end if
         end if
      end do
   80 continue
      block = dble(nblock)
      time = time * block
c
c     convert sums to average and standard deviation values
c
      if (nblock .ne. 0) then
         etot = etot / block
         vtot = sqrt(vtot/block)
         epot = epot / block
         vpot = sqrt(vpot/block)
         ekin = ekin / block
         vkin = sqrt(vkin/block)
         temp = temp / block
         vtemp = sqrt(vtemp/block)
         pres = pres / block
         vpres = sqrt(vpres/block)
         dens = dens / block
         vdens = sqrt(vdens/block)
      end if
c
c     print the averages and overall standard deviations
c
      if (nblock .ne. 0) then
         if (doinfo) then
            write (iout,90)
   90       format ()
         end if
         write (iout,100)  nblock
  100    format (' Total MD Blocks',8x,i12,' Blocks')
         if (doinfo) then
            write (iout,110)  nstep
  110       format (' Steps per Block',8x,i12,' Steps')
            if (time .ge. 1000000.0d0) then
               write (iout,120)  time/1000000.0d0
  120          format (' Simulation Time',8x,f12.4,' Microseconds',/)
            else if (time .ge. 1000.0d0) then
               write (iout,130)  time/1000.0d0
  130          format (' Simulation Time',8x,f12.4,' Nanoseconds',/)
            else
               write (iout,140)  time
  140          format (' Simulation Time',8x,f12.2,' Picoseconds',/)
            end if
         end if
         write (iout,150)  etot,vtot
  150    format (' Total Energy',7x,f16.4,' Kcal/mole   (+/-',
     &              f9.4,')')
         write (iout,160)  epot,vpot
  160    format (' Potential Energy',3x,f16.4,' Kcal/mole   (+/-',
     &              f9.4,')')
         write (iout,170)  ekin,vkin
  170    format (' Kinetic Energy',5x,f16.4,' Kcal/mole   (+/-',
     &              f9.4,')')
         write (iout,180)  temp,vtemp
  180    format (' Temperature',8x,f16.2,' Kelvin      (+/-',f9.2,')')
         write (iout,190)  pres,vpres
  190    format (' Pressure',11x,f16.2,' Atmosphere  (+/-',f9.2,')')
         write (iout,200)  dens,vdens
  200    format (' Density',12x,f16.4,' Grams/cc    (+/-',f9.4,')')
      else
         write (iout,210)
  210    format (/,' MDAVG  --  Input File Contains No Dynamics',
     &              ' Log Information')
      end if
      close (unit=ilog)
      end
