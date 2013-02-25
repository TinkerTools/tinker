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
      include 'freeze.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'mdstuf.i'
      include 'molcul.i'
      include 'moldyn.i'
      include 'mpole.i'
      include 'rgddyn.i'
      include 'rigid.i'
      include 'stodyn.i'
      include 'units.i'
      include 'uprior.i'
      include 'usage.i'
      integer i,j,k,idyn
      integer size,next
      integer lext,freeunit
      real*8 e,ekt,qterm
      real*8 maxwell,speed
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
      bmnmix = 8
      nfree = 0
      irest = 1
      velsave = .false.
      frcsave = .false.
      uindsave = .false.
      friction = 91.0d0
      use_sdarea = .false.
      use_pred = .false.
      iprint = 100
c
c     set default values for temperature and pressure control
c
      thermostat = 'BUSSI'
      tautemp = 0.2d0
      collide = 0.1d0
      do i = 1, maxnose
         vnh(i) = 0.0d0
         qnh(i) = 0.0d0
         gnh(i) = 0.0d0
      end do
      barostat = 'BERENDSEN'
      anisotrop = .false.
      taupres = 2.0d0
      compress = 0.000046d0
      vbar = 0.0d0
      qbar = 0.0d0
      gbar = 0.0d0
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
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         else if (keyword(1:14) .eq. 'BEEMAN-MIXING ') then
            read (string,*,err=10,end=10)  bmnmix
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
         else if (keyword(1:14) .eq. 'POLAR-PREDICT ') then
            use_pred = .true.
         else if (keyword(1:11) .eq. 'THERMOSTAT ') then
            call getword (record,thermostat,next)
            call upcase (thermostat)
         else if (keyword(1:16) .eq. 'TAU-TEMPERATURE ') then
            read (string,*,err=10,end=10)  tautemp
         else if (keyword(1:10) .eq. 'COLLISION ') then
            read (string,*,err=10,end=10)  collide
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
         else if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (string,*,err=10,end=10)  iprint
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
c     perform dynamic allocation of some pointer arrays
c
      if (use_pred) then
         if (associated(udalt))  deallocate (udalt)
         if (associated(upalt))  deallocate (upalt)
         allocate (udalt(maxualt,3,n))
         allocate (upalt(maxualt,3,n))
c
c     set Gear predictor binomial coefficients
c
c        bpred(1) = 6.0d0
c        bpred(2) = -15.0d0
c        bpred(3) = 20.0d0
c        bpred(4) = -15.0d0
c        bpred(5) = 6.0d0
c        bpred(6) = -1.0d0
c
c     set always stable predictor-corrector (ASPC) coefficients
c
         bpred(1) = 22.0d0 / 7.0d0
         bpred(2) = -55.0d0 / 14.0d0
         bpred(3) = 55.0d0 / 21.0d0
         bpred(4) = -22.0d0 / 21.0d0
         bpred(5) = 5.0d0 / 21.0d0
         bpred(6) = -1.0d0 / 42.0d0
c
c     set Challacombe time-reversible coefficients
c
c        bpred(1) = 30.0d0 / 13.0d0
c        bpred(2) = -3.0d0 / 13.0d0
c        bpred(3) = -28.0d0 / 13.0d0
c        bpred(4) = -3.0d0 / 13.0d0
c        bpred(5) = 30.0d0 / 13.0d0
c        bpred(6) = -1.0d0
c
c    initialize prior values of induced dipole moments
c
         nualt = 0
         do i = 1, npole
            do j = 1, 3
               do k = 1, maxualt
                  udalt(k,j,i) = 0.0d0
                  upalt(k,j,i) = 0.0d0
               end do
            end do
         end do
      end if
c
c     enforce use of velocity Verlet with Andersen thermostat
c
      if (thermostat .eq. 'ANDERSEN') then
         if (integrate .eq. 'BEEMAN')  integrate = 'VERLET'
      end if
c
c     enforce use of Bussi thermostat/barostat with integrator
c
      if (integrate .eq. 'BUSSI') then
         thermostat = 'BUSSI'
         barostat = 'BUSSI'
      else if (thermostat.eq.'BUSSI' .and. barostat.eq.'BUSSI') then
         integrate = 'BUSSI'
      end if
c
c     enforce use of Nose-Hoover thermostat/barostat with integrator
c
      if (integrate .eq. 'NOSE-HOOVER') then
         thermostat = 'NOSE-HOOVER'
         barostat = 'NOSE-HOOVER'
      else if (thermostat.eq.'NOSE-HOOVER' .and.
     &         barostat.eq.'NOSE-HOOVER') then
         integrate = 'NOSE-HOOVER'
      end if
c
c     check for use of Monte Carlo barostat with constraints
c
      if (barostat.eq.'MONTECARLO' .and. volscale.eq.'ATOMIC') then
         if (use_rattle) then
            write (iout,40)
   40       format (/,' MDINIT  --  Atom-based Monte Carlo',
     &                 ' Barostat Incompatible with RATTLE')
            call fatal
         end if
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
         if (isothermal .and. thermostat.ne.'ANDERSEN'.and.
     &       integrate.ne.'STOCHASTIC' .and. integrate.ne.'GHMC') then
            if (use_bounds) then
               nfree = nfree - 3
            else
               nfree = nfree - 6
            end if
         end if
         if (barostat .eq. 'BUSSI')  nfree = nfree + 1
      end if
c
c     check for a nonzero number of degrees of freedom
c
      if (nfree .lt. 0)  nfree = 0
      if (debug) then
         write (iout,50)  nfree
   50    format (/,' Number of Degrees of Freedom for Dynamics :',i10)
      end if
      if (nfree .eq. 0) then
         write (iout,60)
   60    format (/,' MDINIT  --  No Degrees of Freedom for Dynamics')
         call fatal
      end if
c
c     set masses for Nose-Hoover thermostat and barostat
c
      if (thermostat .eq. 'NOSE-HOOVER') then
         ekt = gasconst * kelvin
         qterm = ekt * tautemp * tautemp
         do j = 1, maxnose
            if (qnh(j) .eq. 0.0d0)  qnh(j) = qterm
         end do
         qnh(1) = dble(nfree) * qnh(1)
      end if
      if (barostat .eq. 'NOSE-HOOVER') then
         ekt = gasconst * kelvin
         qterm = ekt * taupres * taupres
         qbar = dble(nfree+1) * qterm
      end if
c
c     decide whether to remove center of mass motion
c
      dorest = .true.
      if (irest .eq. 0)  dorest = .false.
      if (nuse. ne. n)  dorest = .false.
      if (integrate .eq. 'STOCHASTIC')  dorest = .false.
      if (integrate .eq. 'GHMC')  dorest = .false.
      if (isothermal .and. thermostat.eq.'ANDERSEN')  dorest = .false.
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
         if (nuse .eq. n)  call mdrest (0)
c
c     set velocities and fast/slow accelerations for RESPA method
c
      else if (integrate .eq. 'RESPA') then
         allocate (derivs(3,n))
         call gradslow (e,derivs)
         do i = 1, n
            if (use(i)) then
               speed = maxwell (mass(i),kelvin)
               call ranvec (vec)
               do j = 1, 3
                  v(j,i) = speed * vec(j)
                  a(j,i) = -convert * derivs(j,i) / mass(i)
               end do
            else
               do j = 1, 3
                  v(j,i) = 0.0d0
                  a(j,i) = 0.0d0
                  aalt(j,i) = 0.0d0
               end do
            end if
         end do
         call gradfast (e,derivs)
         do i = 1, n
            if (use(i)) then
               do j = 1, 3
                  aalt(j,i) = -convert * derivs(j,i) / mass(i)
               end do
            end if
         end do
         deallocate (derivs)
         if (nuse .eq. n)  call mdrest (0)
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
                  aalt(j,i) = a(j,i)
               end do
            else
               do j = 1, 3
                  v(j,i) = 0.0d0
                  a(j,i) = 0.0d0
                  aalt(j,i) = 0.0d0
               end do
            end if
         end do
         deallocate (derivs)
         if (nuse .eq. n)  call mdrest (0)
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
