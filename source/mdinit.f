c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mdinit  --  initialize a dynamics trajectory  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mdinit" initializes the velocities and accelerations
c     for a molecular dynamics trajectory, including restarts
c
c
      subroutine mdinit
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'files.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'mdstuf.i'
      include 'molcul.i'
      include 'moldyn.i'
      include 'rgddyn.i'
      include 'rigid.i'
      include 'shake.i'
      include 'stodyn.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,idyn
      integer size,next
      integer lext,freeunit
      real*8 e,maxwell,speed
      real*8 vec(3)
      real*8, allocatable :: derivs(:,:)
      logical exist
      character*7 ext
      character*20 keyword
      character*120 dynfile
      character*120 record
      character*120 string
c
c
c     set default parameters for the dynamics trajectory
c
      integrate = 'BEEMAN'
      nfree = 0
      irest = 1
      velsave = .false.
      frcsave = .false.
      uindsave = .false.
      friction = 91.0d0
      use_sdarea = .false.
      iprint = 100
c
c     set default values for temperature and pressure control
c
      thermostat = 'BUSSI'
      tautemp = 0.2d0
      collide = 0.1d0
      do i = 1, maxnose
         xnh(i) = 0.0d0
         vnh(i) = 0.0d0
         qnh(i) = 0.0d0
         gnh(i) = 0.0d0
      end do
      barostat = 'BERENDSEN'
      anisotrop = .false.
      taupres = 2.0d0
      compress = 0.000046d0
      eta = 0.0d0
      voltrial = 20
      volmove = 100.0d0
      volscale = 'ATOMIC'
c
c     check for keywords containing any altered parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:10) .eq. 'INTEGRATE ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         else if (keyword(1:16) .eq. 'DEGREES-FREEDOM ') then
            read (string,*,err=10,end=10)  nfree
         else if (keyword(1:15) .eq. 'REMOVE-INERTIA ') then
            read (string,*,err=10,end=10)  irest
         else if (keyword(1:14) .eq. 'SAVE-VELOCITY ') then
            velsave = .true.
         else if (keyword(1:11) .eq. 'SAVE-FORCE ') then
            frcsave = .true.
         else if (keyword(1:13) .eq. 'SAVE-INDUCED ') then
            uindsave = .true.
         else if (keyword(1:9) .eq. 'FRICTION ') then
            read (string,*,err=10,end=10)  friction
         else if (keyword(1:17) .eq. 'FRICTION-SCALING ') then
            use_sdarea = .true.
         else if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (string,*,err=10,end=10)  iprint
         else if (keyword(1:11) .eq. 'THERMOSTAT ') then
            call getword (record,thermostat,next)
            call upcase (thermostat)
         else if (keyword(1:16) .eq. 'TAU-TEMPERATURE ') then
            read (string,*,err=10,end=10)  tautemp
         else if (keyword(1:10) .eq. 'COLLISION ') then
            read (string,*,err=10,end=10)  collide
         else if (keyword(1:10) .eq. 'NOSE-MASS ') then
            read (string,*,err=10,end=10)  (qnh(j),j=1,maxnose)
         else if (keyword(1:9) .eq. 'BAROSTAT ') then
            call getword (record,barostat,next)
            call upcase (barostat)
         else if (keyword(1:15) .eq. 'ANISO-PRESSURE ') then
            anisotrop = .true.
         else if (keyword(1:13) .eq. 'TAU-PRESSURE ') then
            read (string,*,err=10,end=10)  taupres
         else if (keyword(1:9) .eq. 'COMPRESS ') then
            read (string,*,err=10,end=10)  compress
         else if (keyword(1:13) .eq. 'VOLUME-TRIAL ') then
            read (string,*,err=10,end=10)  voltrial
         else if (keyword(1:12) .eq. 'VOLUME-MOVE ') then
            read (string,*,err=10,end=10)  volmove
         else if (keyword(1:13) .eq. 'VOLUME-SCALE ') then
            call getword (record,volscale,next)
            call upcase (volscale)
         end if
   10    continue
      end do
c
c     make sure all atoms or groups have a nonzero mass
c
      if (integrate .eq. 'RIGIDBODY') then
         do i = 1, ngrp
            if (grpmass(i) .le. 0.0d0) then
               grpmass(i) = 1.0d0
               if (igrp(1,i) .le. igrp(2,i)) then
                  totmass = totmass + 1.0d0
                  write (iout,20)  i
   20             format (/,' MDINIT  --  Warning, Mass of Group',i6,
     &                       ' Set to 1.0 for Dynamics')
               end if
            end if
         end do
      else
         do i = 1, n
            if (use(i) .and. mass(i).le.0.0d0) then
               mass(i) = 1.0d0
               totmass = totmass + 1.0d0
               write (iout,30)  i
   30          format (/,' MDINIT  --  Warning, Mass of Atom',i6,
     &                    ' Set to 1.0 for Dynamics')
            end if
         end do
      end if
c
c     enforce use of velocity Verlet with Andersen thermostat
c
      if (thermostat .eq. 'ANDERSEN') then
         if (integrate .eq. 'BEEMAN')  integrate = 'VERLET'
      end if
c
c     set masses for the thermostats in a Nose-Hoover chain
c
      if (thermostat .eq. 'NOSE-HOOVER') then
         if (qnh(1) .eq. 0.0d0)  qnh(1) = 0.1d0
         do j = 2, maxnose
            if (qnh(j) .eq. 0.0d0)  qnh(j) = qnh(j-1)
         end do
      end if
c
c     set the number of degrees of freedom for the system
c
      if (nfree .eq. 0) then
         if (integrate .eq. 'RIGIDBODY') then
            call grpline
            nfree = 6 * ngrp
            do i = 1, ngrp
               size = igrp(2,i) - igrp(1,i) + 1
               if (size .eq. 1)  nfree = nfree - 3
               if (linear(i))  nfree = nfree - 1
            end do
         else
            nfree = 3 * nuse
         end if
         if (use_rattle) then
            nfree = nfree - nrat
            do i = 1, nratx
               nfree = nfree - kratx(i)
            end do
         end if
         if (isothermal .and. integrate.ne.'STOCHASTIC'
     &       .and. thermostat.ne.'ANDERSEN') then
            if (use_bounds) then
               nfree = nfree - 3
            else
               nfree = nfree - 6
            end if
         end if
      end if
c
c     check for a nonzero number of degrees of freedom
c
      if (nfree .lt. 0)  nfree = 0
      if (debug) then
         write (iout,40)  nfree
   40    format (/,' Number of Degrees of Freedom for Dynamics :',i10)
      end if
      if (nfree .eq. 0) then
         write (iout,50)
   50    format (/,' MDINIT  --  No Degrees of Freedom for Dynamics')
         call fatal
      end if
c
c     try to restart using prior velocities and accelerations
c
      dynfile = filename(1:leng)//'.dyn'
      call version (dynfile,'old')
      inquire (file=dynfile,exist=exist)
      if (exist) then
         idyn = freeunit ()
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
         call readdyn (idyn)
         close (unit=idyn)
         call lattice
c
c     set translational velocities for rigid body dynamics
c
      else if (integrate .eq. 'RIGIDBODY') then
         do i = 1, ngrp
            speed = maxwell (grpmass(i),kelvin)
            call ranvec (vec)
            do j = 1, 3
               vcm(j,i) = speed * vec(j)
               wcm(j,i) = 0.0d0
               lm(j,i) = 0.0d0
            end do
         end do
         if (nuse .eq. n)  call mdrest
c
c     set velocities and accelerations for Cartesian dynamics
c
      else
         allocate (derivs(3,n))
         call gradient (e,derivs)
         do i = 1, n
            if (use(i)) then
               speed = maxwell (mass(i),kelvin)
               call ranvec (vec)
               do j = 1, 3
                  v(j,i) = speed * vec(j)
                  a(j,i) = -convert * derivs(j,i) / mass(i)
                  aold(j,i) = a(j,i)
               end do
            else
               do j = 1, 3
                  v(j,i) = 0.0d0
                  a(j,i) = 0.0d0
                  aold(j,i) = 0.0d0
               end do
            end if
         end do
         deallocate (derivs)
         if (nuse .eq. n)  call mdrest
      end if
c
c     check for any prior dynamics coordinate sets
c
      i = 0
      exist = .true.
      do while (exist)
         i = i + 1
         lext = 3
         call numeral (i,ext,lext)
         dynfile = filename(1:leng)//'.'//ext(1:lext)
         inquire (file=dynfile,exist=exist)
         if (.not.exist .and. i.lt.100) then
            lext = 2
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
         if (.not.exist .and. i.lt.10) then
            lext = 1
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
      end do
      nprior = i - 1
      return
      end
