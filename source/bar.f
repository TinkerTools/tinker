c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2013  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program bar  --  thermodynamic values from simulation data  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "bar" computes the free energy, enthalpy and entropy difference
c     between two states via Zwanzig free energy perturbation (FEP)
c     and Bennett acceptance ratio (BAR) methods
c
c     note in mode 1 the code takes as input trajectory archives A
c     and B, originally generated for states 0 and 1, respectively;
c     it then finds the total potential energy for all frames of both
c     trajectories under control of key files for both states 0 and 1;
c     in mode 2, the FEP and BAR algorithms are used to compute the
c     free energy for state 0 --> state 1, and similar estimators
c     provide the enthalpy and entropy for state 0 --> state 1, with
c     a PV term at a pressure asked for if the BAR file has volumes,
c     which must then be given for every frame of both trajectories,
c     and both trajectories must be at the same temperature; in
c     mode 3, a trial lambda energy file written by dynamics over a
c     schedule of lambda windows is read, and the FEP and BAR free
c     energies are found for each pair of adjacent windows and summed
c     over the whole schedule, taken in order of increasing lambda
c
c     modifications for NPT simulations by Chengwen Liu, University
c     of Texas at Austin, October 2015; enthalpy and entropy methods
c     by Aaron Gordon, Washington University, December 2016
c
c     literature references:
c
c     C. H. Bennett, "Efficient Estimation of Free Energy Differences
c     from Monte Carlo Data", Journal of Computational Physics, 22,
c     245-268 (1976)
c
c     M. A. Wyczalkowski, A. Vitalis and R. V. Pappu, "New Estimators
c     for Calculating Solvation Entropy and Enthalpy and Comparative
c     Assessments of Their Accuracy and Precision, Journal of Physical
c     Chemistry, 114, 8166-8180 (2010)  [entropy and enthalpy]
c
c
      program bar
      use iounit
      implicit none
      integer mode
      logical exist,query
      character*240 string
c
c
c     find thermodynamic perturbation procedure to perform
c
      call initial
      mode = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         if (mode.ge.1 .and. mode.le.3)  query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' The Tinker Thermodynamic Perturbation Utility',
     &              ' Can :',
     &           //,4x,'(1) Create BAR File with Perturbed Potential',
     &              ' Energies',
     &           /,4x,'(2) Compute Thermodynamic Values from Tinker',
     &              ' BAR File',
     &           /,4x,'(3) Compute Free Energy from Tinker FEP',
     &              ' File')
         do while (mode.lt.1 .or. mode.gt.3)
            mode = 0
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,40,err=50,end=50)  mode
   40       format (i10)
   50       continue
         end do
      end if
c
c     create Tinker BAR file or compute thermodynamic values
c
      if (mode .eq. 1)  call makebar
      if (mode .eq. 2)  call barcalc
      if (mode .eq. 3)  call fepcalc
      call final
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine makebar  --  store trajectory potential energy  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine makebar
      use boxes
      use files
      use inform
      use iounit
      use keys
      use output
      use titles
      implicit none
      integer i,next
      integer ibar,iarc
      integer ilog,imod
      integer lenga,lengb
      integer ltitlea,ltitleb
      integer nkey0,nkey1
      integer nfrma,nfrmb
      integer nloga,nlogb
      integer maxframe
      integer freeunit
      integer trimtext
      real*8 energy
      real*8 tempa,tempb
      real*8, allocatable :: ua0(:)
      real*8, allocatable :: ua1(:)
      real*8, allocatable :: ub0(:)
      real*8, allocatable :: ub1(:)
      real*8, allocatable :: vola(:)
      real*8, allocatable :: volb(:)
      logical exist,first
      logical recompute
      logical use_log
      character*1 answer
      character*240 filea
      character*240 fileb
      character*240 titlea
      character*240 titleb
      character*240 barfile
      character*240 arcfile
      character*240 logfile
      character*240 record
      character*240 string
      character*240, allocatable :: keys0(:)
      character*240, allocatable :: keys1(:)
c
c
c     get trajectory A and setup mechanics calculation
c
      call initial
      call getarc (iarc)
      close (unit=iarc)
      call mechanic
c
c     store the filename for trajectory A
c
      filea = filename
      lenga = leng
      titlea = title
      ltitlea = ltitle
c
c     perform dynamic allocation of some local arrays
c
      allocate (keys0(nkey))
c
c     store the keyword values for state 0
c
      nkey0 = nkey
      do i = 1, nkey0
         keys0(i) = keyline(i)
      end do
c
c     find the original temperature value for trajectory A
c
      tempa = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  tempa
   10 continue
      do while (tempa .lt. 0.0d0)
         write (iout,20)
   20    format (/,' Enter Trajectory A Temperature in Degrees',
     &              ' K [298] :  ',$)
         read (input,30,err=40)  tempa
   30    format (f20.0)
         if (tempa .le. 0.0d0)  tempa = 298.0d0
   40    continue
      end do
c
c     get trajectory B and setup mechanics calculation
c
      call getarc (iarc)
      close (unit=iarc)
      call mechanic
      silent = .true.
c
c     store the filename for trajectory B
c
      fileb = filename
      lengb = leng
      titleb = title
      ltitleb = ltitle
c
c     perform dynamic allocation of some local arrays
c
      allocate (keys1(nkey))
c
c     store the keyword values for state 1
c
      nkey1 = nkey
      do i = 1, nkey1
         keys1(i) = keyline(i)
      end do
c
c     make sure the keyword array can hold either state, since
c     "getkey" sized it to trajectory B alone and the saved
c     keyword sets are restored into it below
c
      if (allocated(keyline))  deallocate (keyline)
      allocate (keyline(max(nkey0,nkey1)))
c
c     find the original temperature value for trajectory B
c
      tempb = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  tempb
   50 continue
      do while (tempb .lt. 0.0d0)
         write (iout,60)
   60    format (/,' Enter Trajectory B Temperature in Degrees',
     &              ' K [298] :  ',$)
         read (input,70,err=80)  tempb
   70    format (f20.0)
         if (tempb .le. 0.0d0)  tempb = 298.0d0
   80    continue
      end do
c
c     decide whether to use energies from trajectory log files
c
      recompute = .true.
      answer = ' '
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  answer
   90 continue
      if (answer .eq. ' ') then
         answer = 'N'
         write (iout,100)
  100    format (/,' Obtain Energies from Trajectory Logs if',
     &              ' Available [N] :  ',$)
         read (input,110)  record
  110    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'Y')  recompute = .false.
c
c     perform dynamic allocation of some local arrays
c
      maxframe = 1000000
      allocate (ua0(maxframe))
      allocate (ua1(maxframe))
      allocate (ub0(maxframe))
      allocate (ub1(maxframe))
      allocate (vola(maxframe))
      allocate (volb(maxframe))
c
c     reopen the file corresponding to trajectory A
c
      iarc = freeunit ()
      arcfile = filea
      call suffix (arcfile,'arc','old')
      if (archive) then
         open (unit=iarc,file=arcfile,status ='old')
      else if (binary) then
         open (unit=iarc,file=arcfile,form='unformatted',status ='old')
      end if
c
c     read any log energies of trajectory A in state 0, where
c     a log without energies is not used
c
      use_log = .false.
      nloga = 0
      if (.not. recompute) then
         logfile = filea(1:lenga)//'.log'
         call version (logfile,'old')
         inquire (file=logfile,exist=exist)
         if (exist) then
            ilog = freeunit ()
            open (unit=ilog,file=logfile,status='old')
            do while (nloga .lt. maxframe)
               read (ilog,130,err=150,end=150)  record
  130          format (a240)
               if (record(1:18) .eq. ' Current Potential') then
                  string = record(19:240)
                  read (string,*,err=140,end=140)  ua0(nloga+1)
                  nloga = nloga + 1
  140             continue
               end if
            end do
  150       continue
            close (unit=ilog)
            use_log = (nloga .gt. 0)
         end if
      end if
c
c     reset trajectory A using the parameters for state 0
c
      if (.not. use_log) then
         rewind (unit=iarc)
         first = .true.
         call readcart (iarc,first)
         nkey = nkey0
         do i = 1, nkey
            keyline(i) = keys0(i)
         end do
         call mechanic
      end if
c
c     find potential energies for trajectory A in state 0
c
      write (iout,120)
  120 format (/,' Initial Processing for Trajectory A :',/)
      if (use_log) then
         do i = 1, nloga
            if (mod(i,100).eq.0 .or. i.eq.nloga) then
               write (iout,160)  i
               flush (iout)
            end if
         end do
      else
         i = 0
         do while (.not. abort)
            i = i + 1
            call cutoffs
            ua0(i) = energy ()
            call readcart (iarc,first)
            imod = mod(i,100)
            if ((.not.abort.and.imod.eq.0) .or.
     &          (abort.and.imod.ne.0)) then
               write (iout,160)  i
  160          format (7x,'Completed',i8,' Coordinate Frames')
               flush (iout)
            end if
            if (i .ge. maxframe)  abort = .true.
         end do
      end if
c
c     reset trajectory A using the parameters for state 1
c
      rewind (unit=iarc)
      first = .true.
      call readcart (iarc,first)
      nkey = nkey1
      do i = 1, nkey
         keyline(i) = keys1(i)
      end do
      call mechanic
c
c     find potential energies for trajectory A in state 1
c
      if (verbose) then
         write (iout,170)
  170    format (/,' Potential Energy Values for Trajectory A :',
     &           //,7x,'Frame',9x,'State 0',9x,'State 1',12x,'Delta',/)
      end if
      i = 0
      do while (.not. abort)
         i = i + 1
         call cutoffs
         ua1(i) = energy ()
         vola(i) = volbox
         if (verbose .and. (.not.use_log .or. i.le.nloga)) then
            write (iout,180)  i,ua0(i),ua1(i),ua1(i)-ua0(i)
  180       format (i11,2x,3f16.4)
         end if
         call readcart (iarc,first)
         if (i .ge. maxframe)  abort = .true.
      end do
      nfrma = i
      close (unit=iarc)
      if (nfrma .ge. maxframe) then
         write (iout,183)  maxframe,'A'
  183    format (/,' BAR  --  Reading Stopped at the Limit of',i8,
     &              ' Frames for Trajectory ',a)
      end if
c
c     stop if the log energies do not match the archive frames
c
      if (use_log .and. nloga.ne.nfrma) then
         write (iout,185)  nloga,nfrma
  185    format (/,' BAR  --  Log File has',i8,' Energies but Archive',
     &              ' has',i8,' Frames;',
     &           /,' Answer N to Recompute the Energies')
         call fatal
      end if
c
c     save potential energies and volumes for trajectory A
c
      ibar = freeunit ()
      barfile = filea(1:lenga)//'.bar'
      call version (barfile,'new')
      open (unit=ibar,file=barfile,status ='new')
      if (ltitlea .eq. 0) then
         write (ibar,190)  nfrma,tempa
  190    format (i8,f10.2)
      else
         write (ibar,200)  nfrma,tempa,titlea(1:ltitlea)
  200    format (i8,f10.2,2x,a)
      end if
      do i = 1, nfrma
         if (vola(i) .eq. 0.0d0) then
            write (ibar,210)  i,ua0(i),ua1(i)
  210       format (i8,2x,2f18.4)
         else
            write (ibar,220)  i,ua0(i),ua1(i),vola(i)
  220       format (i8,2x,3f18.4)
         end if
      end do
      flush (ibar)
c
c     reopen the file corresponding to trajectory B
c
      iarc = freeunit ()
      arcfile = fileb
      call suffix (arcfile,'arc','old')
      if (archive) then
         open (unit=iarc,file=arcfile,status ='old')
      else if (binary) then
         open (unit=iarc,file=arcfile,form='unformatted',status ='old')
      end if
c
c     read any log energies of trajectory B in state 1, where
c     a log without energies is not used
c
      use_log = .false.
      nlogb = 0
      if (.not. recompute) then
         logfile = fileb(1:lengb)//'.log'
         call version (logfile,'old')
         inquire (file=logfile,exist=exist)
         if (exist) then
            ilog = freeunit ()
            open (unit=ilog,file=logfile,status='old')
            do while (nlogb .lt. maxframe)
               read (ilog,240,err=260,end=260)  record
  240          format (a240)
               if (record(1:18) .eq. ' Current Potential') then
                  string = record(19:240)
                  read (string,*,err=250,end=250)  ub1(nlogb+1)
                  nlogb = nlogb + 1
  250             continue
               end if
            end do
  260       continue
            close (unit=ilog)
            use_log = (nlogb .gt. 0)
         end if
      end if
c
c     reset trajectory B using the parameters for state 1
c
      rewind (unit=iarc)
      first = .true.
      call readcart (iarc,first)
      nkey = nkey1
      do i = 1, nkey
         keyline(i) = keys1(i)
      end do
      call mechanic
c
c     find potential energies for trajectory B in state 1
c
      write (iout,230)
  230 format (/,' Initial Processing for Trajectory B :',/)
      if (use_log) then
         do i = 1, nlogb
            if (mod(i,100).eq.0 .or. i.eq.nlogb) then
               write (iout,270)  i
               flush (iout)
            end if
         end do
      else
         i = 0
         do while (.not. abort)
            i = i + 1
            call cutoffs
            ub1(i) = energy ()
            call readcart (iarc,first)
            imod = mod(i,100)
            if ((.not.abort.and.imod.eq.0) .or.
     &          (abort.and.imod.ne.0)) then
               write (iout,270)  i
  270          format (7x,'Completed',i8,' Coordinate Frames')
               flush (iout)
            end if
            if (i .ge. maxframe)  abort = .true.
         end do
      end if
c
c     reset trajectory B using the parameters for state 0
c
      rewind (unit=iarc)
      first = .true.
      call readcart (iarc,first)
      nkey = nkey0
      do i = 1, nkey
         keyline(i) = keys0(i)
      end do
      call mechanic
c
c     find potential energies for trajectory B in state 0
c
      if (verbose) then
         write (iout,280)
  280    format (/,' Potential Energy Values for Trajectory B :',
     &           //,7x,'Frame',9x,'State 0',9x,'State 1',12x,'Delta',/)
      end if
      i = 0
      do while (.not. abort)
         i = i + 1
         call cutoffs
         ub0(i) = energy ()
         volb(i) = volbox
         if (verbose .and. (.not.use_log .or. i.le.nlogb)) then
            write (iout,290)  i,ub0(i),ub1(i),ub0(i)-ub1(i)
  290       format (i11,2x,3f16.4)
         end if
         call readcart (iarc,first)
         if (i .ge. maxframe)  abort = .true.
      end do
      nfrmb = i
      close (unit=iarc)
      if (nfrmb .ge. maxframe)  write (iout,183)  maxframe,'B'
c
c     stop if the log energies do not match the archive frames,
c     and remove the BAR file begun for trajectory A
c
      if (use_log .and. nlogb.ne.nfrmb) then
         close (unit=ibar,status='delete')
         write (iout,185)  nlogb,nfrmb
         call fatal
      end if
c
c     save potential energies and volumes for trajectory B
c
      if (ltitleb .eq. 0) then
         write (ibar,300)  nfrmb,tempb
  300    format (i8,f10.2)
      else
         write (ibar,310)  nfrmb,tempb,titleb(1:ltitleb)
  310    format (i8,f10.2,2x,a)
      end if
      do i = 1, nfrmb
         if (volb(i) .eq. 0.0d0) then
            write (ibar,320)  i,ub0(i),ub1(i)
  320       format (i8,2x,2f18.4)
         else
            write (ibar,330)  i,ub0(i),ub1(i),volb(i)
  330       format (i8,2x,3f18.4)
         end if
      end do
      close (unit=ibar)
      write (iout,340)  barfile(1:trimtext(barfile))
  340 format (/,' Potential Energy Values Written To :  ',a)
c
c     perform deallocation of some local arrays
c
      deallocate (keys0)
      deallocate (keys1)
c
c     perform deallocation of some local arrays
c
      deallocate (ua0)
      deallocate (ua1)
      deallocate (ub0)
      deallocate (ub1)
      deallocate (vola)
      deallocate (volb)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine barcalc  --  get thermodynamics via FEP & BAR  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine barcalc
      use files
      use inform
      use iounit
      use titles
      implicit none
      integer i,j,k,ibar
      integer next
      integer ltitlea,ltitleb
      integer nfrma,nfrmb
      integer nfrm,nbst
      integer starta,startb
      integer stopa,stopb
      integer stepa,stepb
      integer maxa,maxb,istat
      integer nvola,nvolb
      integer barnbst
      integer freeunit
      integer trimtext
      integer nboot(2)
      real*8 tempa,tempb
      real*8 patm
      real*8 dgfep(4),dgbar(4)
      real*8 sdu(2),uave(4)
      real*8 hbar(3)
      real*8, allocatable :: ua0(:)
      real*8, allocatable :: ua1(:)
      real*8, allocatable :: ub0(:)
      real*8, allocatable :: ub1(:)
      real*8, allocatable :: vola(:)
      real*8, allocatable :: volb(:)
      logical exist,query
      logical usevol
      character*240 record
      character*240 string
      character*240 titlea
      character*240 titleb
      character*240 barfile
c
c
c     ask the user for file with potential energies and volumes
c
      call nextarg (barfile,exist)
      if (exist) then
         call basefile (barfile)
         call suffix (barfile,'bar','old')
         inquire (file=barfile,exist=exist)
      end if
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter Potential Energy BAR File Name :  ',$)
         read (input,20)  barfile
   20    format (a240)
         call basefile (barfile)
         call suffix (barfile,'bar','old')
         inquire (file=barfile,exist=exist)
      end do
c
c     find the number of frames listed for each trajectory
c
      maxa = 0
      maxb = 0
      ibar = freeunit ()
      open (unit=ibar,file=barfile,status='old')
      rewind (unit=ibar)
      read (ibar,*,err=25,end=25)  maxa
      do i = 1, maxa
         read (ibar,'()',err=25,end=25)
      end do
      read (ibar,*,err=25,end=25)  maxb
   25 continue
      close (unit=ibar)
      maxa = max(0,maxa)
      maxb = max(0,maxb)
c
c     perform dynamic allocation of some local arrays
c
      allocate (ua0(max(1,maxa)))
      allocate (ua1(max(1,maxa)))
      allocate (vola(max(1,maxa)))
      allocate (ub0(max(1,maxb)))
      allocate (ub1(max(1,maxb)))
      allocate (volb(max(1,maxb)))
c
c     set beginning and ending frame for trajectory A
c
      starta = 0
      stopa = 0
      stepa = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=30,end=30)  starta
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  stopa
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  stepa
   30 continue
      if (query) then
         write (iout,40)
   40    format (/,' First & Last Frame and Step Increment',
     &              ' for Trajectory A :  ',$)
         read (input,50)  record
   50    format (a120)
         read (record,*,err=60,end=60)  starta,stopa,stepa
   60    continue
      end if
      if (starta .lt. 1)  starta = 1
      if (stopa.lt.1 .or. stopa.gt.maxa)  stopa = maxa
      if (stepa .lt. 1)  stepa = 1
c
c     set beginning and ending frame for trajectory B
c
      startb = 0
      stopb = 0
      stepb = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=70,end=70)  startb
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=70,end=70)  stopb
      call nextarg (string,exist)
      if (exist)  read (string,*,err=70,end=70)  stepb
   70 continue
      if (query) then
         write (iout,80)
   80    format (/,' First & Last Frame and Step Increment',
     &              ' for Trajectory B :  ',$)
         read (input,90)  record
   90    format (a120)
         read (record,*,err=100,end=100)  startb,stopb,stepb
  100    continue
      end if
      if (startb .lt. 1)  startb = 1
      if (stopb.lt.1 .or. stopb.gt.maxb)  stopb = maxb
      if (stepb .lt. 1)  stepb = 1
c
c     read potential energies and volumes for trajectory A
c
      ibar = freeunit ()
      open (unit=ibar,file=barfile,status='old')
      rewind (unit=ibar)
      nfrma = 0
      j = 0
      tempa = 0.0d0
      next = 1
      read (ibar,110,err=250,end=250)  record
  110 format (a240)
      call gettext (record,string,next)
      read (string,*,err=250,end=250)  nfrma
      call gettext (record,string,next)
      read (string,*,err=115,end=115)  tempa
  115 continue
      titlea = record(next:trimtext(record))
      call trimhead (titlea)
      ltitlea = trimtext(titlea)
      do i = 1, starta-1
         read (ibar,120,err=160,end=160)
  120    format ()
      end do
      j = 0
      stopa = min(stopa,nfrma)
      if (stopa .ge. starta) then
         stopa = (stopa-starta)/stepa + 1
      else
         stopa = 0
      end if
      do i = 1, stopa
         read (ibar,130,err=160,end=160)  record
  130    format (a240)
         j = j + 1
         ua0(j) = 0.0d0
         ua1(j) = 0.0d0
         vola(j) = 0.0d0
         read (record,*,err=135,end=135)  k,ua0(j),ua1(j)
         read (record,*,err=135,end=140)  k,ua0(j),ua1(j),vola(j)
         goto 140
  135    continue
         write (iout,137)  record(1:trimtext(record))
  137    format (/,' BAR  --  Unreadable Line in BAR File :',//,' ',a)
         call fatal
  140    continue
         do k = 1, stepa-1
            read (ibar,150,err=160,end=160)
  150       format ()
         end do
      end do
  160 continue
      nfrma = j
c
c     reset the file position to the beginning of trajectory B
c
      rewind (unit=ibar)
      next = 1
      read (ibar,170,err=190,end=190)  record
  170 format (a240)
      call gettext (record,string,next)
      read (string,*,err=190,end=190)  k
      do i = 1, k
         read (ibar,180,err=190,end=190)
  180    format ()
      end do
  190 continue
c
c     read potential energies and volumes for trajectory B
c
      nfrmb = 0
      j = 0
      tempb = 0.0d0
      next = 1
      read (ibar,200,err=250,end=250)  record
  200 format (a240)
      call gettext (record,string,next)
      read (string,*,err=250,end=250)  nfrmb
      call gettext (record,string,next)
      read (string,*,err=205,end=205)  tempb
  205 continue
      titleb = record(next:trimtext(record))
      call trimhead (titleb)
      ltitleb = trimtext(titleb)
      do i = 1, startb-1
         read (ibar,210,err=250,end=250)
  210    format ()
      end do
      j = 0
      stopb = min(stopb,nfrmb)
      if (stopb .ge. startb) then
         stopb = (stopb-startb)/stepb + 1
      else
         stopb = 0
      end if
      do i = 1, stopb
         read (ibar,220,err=250,end=250)  record
  220    format (a240)
         j = j + 1
         ub0(j) = 0.0d0
         ub1(j) = 0.0d0
         volb(j) = 0.0d0
         read (record,*,err=225,end=225)  k,ub0(j),ub1(j)
         read (record,*,err=225,end=230)  k,ub0(j),ub1(j),volb(j)
         goto 230
  225    continue
         write (iout,137)  record(1:trimtext(record))
         call fatal
  230    continue
         do k = 1, stepb-1
            read (ibar,240,err=250,end=250)
  240       format ()
         end do
      end do
  250 continue
      nfrmb = j
      close (unit=ibar)
c
c     provide info about trajectories and number of frames
c
      write (iout,260)
  260 format (/,' Simulation Trajectory A and Thermodynamic',
     &           ' State 0 :  ',/)
      if (ltitlea .ne. 0) then
         write (iout,270)  titlea(1:ltitlea)
  270    format (' ',a)
      end if
      write (iout,280)  nfrma,tempa
  280 format (' Number of Frames :',4x,i8,/,' Temperature :',7x,f10.2)
      write (iout,290)
  290 format (/,' Simulation Trajectory B and Thermodynamic',
     &           ' State 1 :  ',/)
      if (ltitleb .ne. 0) then
         write (iout,300)  titleb(1:ltitleb)
  300    format (' ',a)
      end if
      write (iout,310)  nfrmb,tempb
  310 format (' Number of Frames :',4x,i8,/,' Temperature :',7x,f10.2)
c
c     stop if either trajectory has no frames selected
c
      if (nfrma .eq. 0) then
         write (iout,320)  'A'
  320    format (/,' BAR  --  No Frames Selected from Trajectory ',a)
         call fatal
      end if
      if (nfrmb .eq. 0) then
         write (iout,320)  'B'
         call fatal
      end if
c
c     stop if either trajectory has no temperature
c
      if (.not. (tempa .gt. 0.0d0)) then
         write (iout,325)  'A'
  325    format (/,' BAR  --  No Temperature Found for Trajectory ',a)
         call fatal
      end if
      if (.not. (tempb .gt. 0.0d0)) then
         write (iout,325)  'B'
         call fatal
      end if
c
c     stop if the trajectories differ in temperature
c
      if (abs(tempa-tempb) .gt. 0.01d0) then
         write (iout,330)  tempa,tempb
  330    format (/,' BAR  --  Trajectory Temperatures of',f8.2,
     &              ' and',f8.2,' K Differ;',
     &           /,' FEP and BAR Require a Single Temperature')
         call fatal
      end if
c
c     get the pressure for the PV terms if volumes are present
c
      nvola = 0
      do i = 1, nfrma
         if (vola(i) .ne. 0.0d0)  nvola = nvola + 1
      end do
      nvolb = 0
      do i = 1, nfrmb
         if (volb(i) .ne. 0.0d0)  nvolb = nvolb + 1
      end do
      usevol = (nvola+nvolb .gt. 0)
c
c     the PV terms need a volume for every frame of both trajectories,
c     or for none, since a missing volume would enter them as zero
c
      if (usevol) then
         if (nvola .lt. nfrma) then
            write (iout,335)  nfrma-nvola,nfrma,'A'
  335       format (/,' BAR  --  Volumes Missing for',i8,' of',i8,
     &                 ' Frames of Trajectory ',a,';',
     &              /,' PV Terms Need a Volume for Every Frame')
            call fatal
         end if
         if (nvolb .lt. nfrmb) then
            write (iout,335)  nfrmb-nvolb,nfrmb,'B'
            call fatal
         end if
      end if
      patm = 1.0d0
      if (usevol) then
         patm = -1.0d0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=340,end=340)  patm
  340    continue
         do while (patm .lt. 0.0d0)
            write (iout,350)
  350       format (/,' Enter Pressure in Atmospheres [1.0] :  ',$)
            read (input,360,err=370,end=380)  patm
  360       format (f20.0)
            if (patm .le. 0.0d0)  patm = 1.0d0
  370       continue
         end do
  380    continue
         if (patm .le. 0.0d0)  patm = 1.0d0
         write (iout,390)  patm
  390    format (/,' Pressure (Atm) :',4x,f10.2)
      end if
c
c     set the number of bootstrap trials to be generated
c
      nfrm = max(nfrma,nfrmb)
      nbst = barnbst (nfrm)
c
c     find the free energy, enthalpy and entropy via FEP and BAR
c
      call barpair (nfrma,nfrmb,ua0,ua1,ub0,ub1,vola,volb,tempa,
     &              tempb,patm,nbst,.true.,dgfep,dgbar,nboot,sdu,
     &              uave,hbar,istat)
      if (istat .ge. 1)  call fatal
c
c     perform deallocation of some local arrays
c
      deallocate (ua0)
      deallocate (ua1)
      deallocate (ub0)
      deallocate (ub1)
      deallocate (vola)
      deallocate (volb)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fepcalc  --  free energy from Tinker FEP file  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fepcalc" reads the trial lambda energies of a Tinker FEP file
c     written by dynamics over a schedule of lambda windows, finds the
c     free energy change between each pair of adjacent windows via the
c     FEP and BAR methods, and sums them over the whole schedule
c
c     note the windows are taken from the lambda values of the data
c     records, so a window without samples of its own is still found
c     from the trial energies of its neighbors; the windows are put
c     in order of increasing lambda for either direction of the
c     simulation, so the values are for lambda 0 --> 1, and all of
c     the samples of each window are used; the native and trial
c     energies of a sample share the same periodic box, whose volume
c     is in a fifth column of newer files, and the volumes are used
c     for the PV enthalpy terms only if the file header also gives
c     the pressure of an isobaric simulation; the errors come from
c     bootstrap trials that resample the samples of every window at
c     once, so the errors of the totals include the correlation of
c     pairs that share a window
c
c
      subroutine fepcalc
      use iounit
      use units
      implicit none
      integer i,j,k,m,n
      integer ifep,nrow
      integer nsamp,nwin,nwhead
      integer npair,norphan
      integer nvol,nzero,izero
      integer nbst,nbstin,nsmax
      integer na,nb,iw,istat
      integer kt,maxiter
      integer ntota,ntotb
      integer npboot,nwork
      integer ntf,ntc,nlist
      integer barnbst
      integer freeunit
      integer trimtext
      integer nboot(2)
      integer, allocatable :: rstep(:)
      integer, allocatable :: sfirst(:)
      integer, allocatable :: slast(:)
      integer, allocatable :: swin(:)
      integer, allocatable :: wsel(:)
      integer, allocatable :: wseen(:)
      integer, allocatable :: woff(:)
      integer, allocatable :: wsamp(:)
      integer, allocatable :: draw(:)
      integer, allocatable :: mapa(:)
      integer, allocatable :: mapb(:)
      integer, allocatable :: ia(:)
      integer, allocatable :: ib(:)
      integer, allocatable :: pstat(:)
      integer, allocatable :: pna(:)
      integer, allocatable :: pnb(:)
      integer, allocatable :: paoff(:)
      integer, allocatable :: pboff(:)
      integer, allocatable :: pskip(:)
      integer, allocatable :: pfail(:)
      integer, allocatable :: pconv(:)
      integer, allocatable :: ptry(:)
      integer, allocatable :: plist(:)
      real*8 temp,tol,lmda
      real*8 patm,pres,rt
      real*8 u0,u1,v0,v1
      real*8 lfirst,lsecond
      real*8 cumul,random
      real*8 barsd,hdir,tds
      real*8 cfore,cback
      real*8 hfore,hback
      real*8 ctrial,htrial
      real*8 sum(4),val(5)
      real*8 tsum(5),tsum2(5)
      real*8 dgfep(4),dgbar(4)
      real*8 sdu(2),uave(4)
      real*8 hbar(3)
      real*8, allocatable :: rsim(:)
      real*8, allocatable :: rtrl(:)
      real*8, allocatable :: rene(:)
      real*8, allocatable :: rvol(:)
      real*8, allocatable :: wlmda(:)
      real*8, allocatable :: pa0(:)
      real*8, allocatable :: pa1(:)
      real*8, allocatable :: pva(:)
      real*8, allocatable :: pb0(:)
      real*8, allocatable :: pb1(:)
      real*8, allocatable :: pvb(:)
      real*8, allocatable :: wa(:)
      real*8, allocatable :: hwa(:)
      real*8, allocatable :: ea(:)
      real*8, allocatable :: wb(:)
      real*8, allocatable :: hwb(:)
      real*8, allocatable :: eb(:)
      real*8, allocatable :: pxa(:)
      real*8, allocatable :: pxb(:)
      real*8, allocatable :: pcref(:)
      real*8, allocatable :: pfep(:,:)
      real*8, allocatable :: pbar(:,:)
      real*8, allocatable :: psdu(:,:)
      real*8, allocatable :: phs(:,:)
      real*8, allocatable :: psum(:,:)
      logical exist,header
      logical found,totals
      logical hasvol,usevol
      logical allf,allc,tconv
      character*240 record
      character*240 string
      character*240 fepfile
      character*80 line
      external random
c
c
c     ask the user for the file with the trial lambda energies
c
      call nextarg (fepfile,exist)
      if (exist) then
         call basefile (fepfile)
         call suffix (fepfile,'fep','old')
         inquire (file=fepfile,exist=exist)
      end if
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter Trial Lambda Energy FEP File Name :  ',$)
         read (input,20)  fepfile
   20    format (a240)
         call basefile (fepfile)
         call suffix (fepfile,'fep','old')
         inquire (file=fepfile,exist=exist)
      end do
c
c     get any requested number of bootstrap trials for each pair
c
      nbstin = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=70,end=70)  nbstin
   70 continue
      if (nbstin .lt. 0)  nbstin = 0
      if (nbstin .eq. 1)  nbstin = 2
c
c     read the header values and count the data records
c
      ifep = freeunit ()
      open (unit=ifep,file=fepfile,status='old')
      rewind (unit=ifep)
      temp = 0.0d0
      patm = 0.0d0
      nwhead = 0
      nrow = 0
      hasvol = .false.
      header = .true.
   80 continue
      read (ifep,90,err=100,end=100)  record
   90 format (a240)
      if (header) then
         k = index(record,'Number of Lambda Windows')
         if (k .ne. 0)  read (record(k+24:240),*,err=80,end=80)  nwhead
         k = index(record,'Simulation Temperature')
         if (k .ne. 0)  read (record(k+22:240),*,err=80,end=80)  temp
         k = index(record,'Simulation Pressure')
         if (k .ne. 0)  read (record(k+19:240),*,err=80,end=80)  patm
         if (index(record,'Sim Lambda') .ne. 0) then
            header = .false.
            hasvol = (index(record,'Volume') .ne. 0)
         end if
      else if (record .ne. ' ') then
         nrow = nrow + 1
      end if
      goto 80
  100 continue
      if (header) then
         write (iout,110)
  110    format (/,' FEPCALC  --  No Trial Lambda Energies Found',
     &              ' in the FEP File')
         call fatal
      end if
      if (temp .le. 0.0d0) then
         write (iout,120)
  120    format (/,' FEPCALC  --  Simulation Temperature not Found',
     &              ' in the FEP File')
         call fatal
      end if
      usevol = (hasvol .and. patm.gt.0.0d0)
      if (.not. usevol)  patm = 1.0d0
      pres = 0.0d0
      if (usevol)  pres = patm / prescon
c
c     perform dynamic allocation of some local arrays
c
      n = max(1,nrow)
      allocate (rstep(n))
      allocate (rsim(n))
      allocate (rtrl(n))
      allocate (rene(n))
      allocate (rvol(n))
      allocate (sfirst(n))
      allocate (slast(n))
      allocate (swin(n))
      allocate (wlmda(2*n))
c
c     read the data records, ignoring a final record that was cut
c     short; a whole record holds 77 columns as written by "savefep",
c     or 59 columns in older files without the volume column
c
      rewind (unit=ifep)
  130 continue
      read (ifep,140,err=180,end=180)  record
  140 format (a240)
      if (index(record,'Sim Lambda') .eq. 0)  goto 130
      m = 0
      n = 0
  150 continue
      read (ifep,160,err=180,end=180)  record
  160 format (a240)
      if (record .eq. ' ')  goto 150
      m = m + 1
      rvol(m) = 0.0d0
      if (hasvol) then
         if (m.eq.nrow .and. trimtext(record).lt.77)  goto 180
         read (record,*,err=170,end=170)  rstep(m),rsim(m),rtrl(m),
     &                                    rene(m),rvol(m)
      else
         if (m.eq.nrow .and. trimtext(record).lt.59)  goto 180
         read (record,*,err=170,end=170)  rstep(m),rsim(m),rtrl(m),
     &                                    rene(m)
      end if
      n = m
      goto 150
  170 continue
      if (m .lt. nrow) then
         write (iout,175)  record(1:trimtext(record))
  175    format (/,' FEPCALC  --  Unreadable Record in the FEP File :',
     &           //,' ',a)
         call fatal
      end if
  180 continue
      close (unit=ifep)
      nrow = n
c
c     a file with a pressure but only zero box volumes gets no PV
c     terms, while volumes missing from only some samples are an error
c
      tol = 1.0d-6
      if (usevol) then
         nvol = 0
         nzero = 0
         izero = 0
         do i = 1, nrow
            if (abs(rtrl(i)-rsim(i)) .lt. tol) then
               if (rvol(i) .ne. 0.0d0) then
                  nvol = nvol + 1
               else
                  nzero = nzero + 1
                  if (izero .eq. 0)  izero = i
               end if
            end if
         end do
         if (nvol .eq. 0) then
            usevol = .false.
            patm = 1.0d0
            pres = 0.0d0
            write (iout,182)
  182       format (/,' FEPCALC  --  File gives a Pressure but No Box',
     &                 ' Volumes;',
     &              /,' PV Terms are not Included')
         else if (nzero .gt. 0) then
            write (iout,184)  rstep(izero)
  184       format (/,' FEPCALC  --  Box Volume Missing at Step',i11,
     &              /,' in a File with Volumes for Other Samples')
            call fatal
         end if
      end if
c
c     group the records into samples, each with its native record
c     first, followed by those at the trial lambda values
c
      nsamp = 0
      norphan = 0
      do i = 1, nrow
         if (abs(rtrl(i)-rsim(i)) .lt. tol) then
            nsamp = nsamp + 1
            sfirst(nsamp) = i
            slast(nsamp) = i
         else if (nsamp .eq. 0) then
            norphan = norphan + 1
         else
            k = sfirst(nsamp)
            if (rstep(i).eq.rstep(k) .and. slast(nsamp).eq.i-1
     &             .and. abs(rsim(i)-rsim(k)).lt.tol) then
               slast(nsamp) = i
            else
               norphan = norphan + 1
            end if
         end if
      end do
c
c     at least two of the windows must have samples of their own
c
      found = .false.
      if (nsamp .gt. 0) then
         lfirst = rsim(sfirst(1))
         do i = 2, nsamp
            lsecond = rsim(sfirst(i))
            if (abs(lsecond-lfirst) .ge. tol) then
               found = .true.
               goto 190
            end if
         end do
      end if
  190 continue
      if (.not. found) then
         write (iout,200)
  200    format (/,' FEPCALC  --  Fewer than Two Lambda Windows',
     &              ' have Samples')
         return
      end if
c
c     find the windows from the simulation and trial lambda values
c
      nwin = 0
      do i = 1, nrow
         do j = 1, 2
            if (j .eq. 1)  lmda = rsim(i)
            if (j .eq. 2)  lmda = rtrl(i)
            found = .false.
            do k = 1, nwin
               if (abs(wlmda(k)-lmda) .lt. tol)  found = .true.
            end do
            if (.not. found) then
               nwin = nwin + 1
               wlmda(nwin) = lmda
            end if
         end do
      end do
c
c     sort the windows into order of increasing lambda
c
      do i = 2, nwin
         lmda = wlmda(i)
         j = i - 1
         do while (j .ge. 1)
            if (wlmda(j) .gt. lmda) then
               wlmda(j+1) = wlmda(j)
               j = j - 1
            else
               goto 210
            end if
         end do
  210    continue
         wlmda(j+1) = lmda
      end do
c
c     assign each sample to its window and count those of each
c
      allocate (wsel(nwin))
      allocate (wseen(nwin))
      do i = 1, nwin
         wsel(i) = 0
      end do
      do i = 1, nsamp
         lmda = rsim(sfirst(i))
         do k = 1, nwin
            if (abs(wlmda(k)-lmda) .lt. tol)  swin(i) = k
         end do
         wsel(swin(i)) = wsel(swin(i)) + 1
      end do
      nsmax = 0
      do i = 1, nwin
         nsmax = max(nsmax,wsel(i))
      end do
c
c     list the samples of each window in file order
c
      npair = nwin - 1
      allocate (woff(nwin+1))
      allocate (wsamp(max(1,nsamp)))
      woff(1) = 0
      do i = 1, nwin
         woff(i+1) = woff(i) + wsel(i)
         wseen(i) = 0
      end do
      do i = 1, nsamp
         iw = swin(i)
         wseen(iw) = wseen(iw) + 1
         wsamp(woff(iw)+wseen(iw)) = i
      end do
c
c     perform dynamic allocation of some local arrays
c
      n = max(1,woff(nwin+1))
      allocate (draw(n))
      allocate (mapa(n))
      allocate (mapb(n))
      allocate (pa0(n))
      allocate (pa1(n))
      allocate (pva(n))
      allocate (pb0(n))
      allocate (pb1(n))
      allocate (pvb(n))
      allocate (wa(n))
      allocate (hwa(n))
      allocate (ea(n))
      allocate (wb(n))
      allocate (hwb(n))
      allocate (eb(n))
      n = max(1,nsmax)
      allocate (ia(n))
      allocate (ib(n))
      allocate (pstat(npair))
      allocate (pna(npair))
      allocate (pnb(npair))
      allocate (paoff(npair))
      allocate (pboff(npair))
      allocate (pskip(npair))
      allocate (pfail(npair))
      allocate (pconv(npair))
      allocate (ptry(npair))
      allocate (plist(npair))
      allocate (pfep(4,npair))
      allocate (pbar(4,npair))
      allocate (psdu(2,npair))
      allocate (phs(3,npair))
      allocate (pxa(npair))
      allocate (pxb(npair))
      allocate (pcref(npair))
      allocate (psum(8,npair))
c
c     gather the frames of each pair of adjacent windows, where a
c     sample of window "iw" gives a state 0 frame if it has the energy
c     at window "iw+1", and a sample of window "iw+1" gives a state 1
c     frame if it has the energy at window "iw"; then find the free
c     energy values from all of the frames of each pair
c
      rt = gasconst * temp
      npboot = 0
      nwork = 0
      ntota = 0
      ntotb = 0
      do iw = 1, npair
         paoff(iw) = ntota
         pboff(iw) = ntotb
         na = 0
         nb = 0
         pskip(iw) = 0
         do m = 1, wsel(iw)
            i = wsamp(woff(iw)+m)
            mapa(woff(iw)+m) = 0
            j = 0
            do k = sfirst(i)+1, slast(i)
               if (j.eq.0 .and. abs(rtrl(k)-wlmda(iw+1)).lt.tol)  j = k
            end do
            if (j .ne. 0) then
               na = na + 1
               pa0(ntota+na) = rene(sfirst(i))
               pa1(ntota+na) = rene(j)
               pva(ntota+na) = 0.0d0
               if (usevol)  pva(ntota+na) = rvol(sfirst(i))
               mapa(woff(iw)+m) = na
            else
               pskip(iw) = pskip(iw) + 1
            end if
         end do
         do m = 1, wsel(iw+1)
            i = wsamp(woff(iw+1)+m)
            mapb(woff(iw+1)+m) = 0
            j = 0
            do k = sfirst(i)+1, slast(i)
               if (j.eq.0 .and. abs(rtrl(k)-wlmda(iw)).lt.tol)  j = k
            end do
            if (j .ne. 0) then
               nb = nb + 1
               pb1(ntotb+nb) = rene(sfirst(i))
               pb0(ntotb+nb) = rene(j)
               pvb(ntotb+nb) = 0.0d0
               if (usevol)  pvb(ntotb+nb) = rvol(sfirst(i))
               mapb(woff(iw+1)+m) = nb
            else
               pskip(iw) = pskip(iw) + 1
            end if
         end do
         ntota = ntota + na
         ntotb = ntotb + nb
         pna(iw) = na
         pnb(iw) = nb
         if (na.eq.0 .or. nb.eq.0) then
            pstat(iw) = -1
         else
            k = paoff(iw) + 1
            j = pboff(iw) + 1
            call barpair (na,nb,pa0(k),pa1(k),pb0(j),pb1(j),pva(k),
     &                    pvb(j),temp,temp,patm,0,.false.,dgfep,dgbar,
     &                    nboot,sdu,uave,hbar,istat)
            pstat(iw) = istat
            pfep(1,iw) = dgfep(1)
            pfep(3,iw) = dgfep(3)
            pbar(1,iw) = dgbar(1)
            pbar(2,iw) = dgbar(2)
            psdu(1,iw) = sdu(1)
            psdu(2,iw) = sdu(2)
            phs(1,iw) = hbar(1)
            phs(3,iw) = hbar(3)
            pcref(iw) = 0.0d0
            if (istat .eq. 0)  pcref(iw) = dgbar(1)
            call bootprep (na,nb,pa0(k),pa1(k),pb0(j),pb1(j),pva(k),
     &                     pvb(j),rt,rt,pres,pcref(iw),pxa(iw),pxb(iw),
     &                     wa(k),hwa(k),wb(j),hwb(j),ea(k),eb(j))
            if (istat .eq. 0) then
               npboot = npboot + 1
               nwork = nwork + na + nb
            end if
         end if
      end do
c
c     set the number of bootstrap trials, shared by all of the pairs,
c     from the frames of the pairs used in each trial, with no trials
c     if no pair has a BAR free energy
c
      nbst = nbstin
      if (nbst .eq. 0)  nbst = barnbst (nwork)
      if (npboot .eq. 0)  nbst = 0
      maxiter = 100
c
c     apply the estimators of every pair to bootstrap trials that
c     resample the samples of each window once, so that the totals
c     over the pairs include the correlation between pairs sharing a
c     window; each BAR trial starts from the BAR value of its pair
c
      do iw = 1, npair
         do k = 1, 8
            psum(k,iw) = 0.0d0
         end do
         pconv(iw) = 0
         ptry(iw) = 0
      end do
      do k = 1, 5
         tsum(k) = 0.0d0
         tsum2(k) = 0.0d0
      end do
      ntf = 0
      ntc = 0
      do kt = 1, nbst
         do iw = 1, nwin
            do m = 1, wsel(iw)
               draw(woff(iw)+m) = int(dble(wsel(iw))*random()) + 1
            end do
         end do
         allf = .true.
         allc = .true.
         do k = 1, 5
            val(k) = 0.0d0
         end do
         do iw = 1, npair
            if (pstat(iw) .ne. 0)  goto 230
            na = 0
            do m = 1, wsel(iw)
               j = mapa(woff(iw)+draw(woff(iw)+m))
               if (j .ne. 0) then
                  na = na + 1
                  ia(na) = j
               end if
            end do
            nb = 0
            do m = 1, wsel(iw+1)
               j = mapb(woff(iw+1)+draw(woff(iw+1)+m))
               if (j .ne. 0) then
                  nb = nb + 1
                  ib(nb) = j
               end if
            end do
            if (na.eq.0 .or. nb.eq.0) then
               allf = .false.
               allc = .false.
               goto 230
            end if
            k = paoff(iw) + 1
            j = pboff(iw) + 1
            call bootfep (na,ia,pa0(k),pa1(k),pva(k),wa(k),hwa(k),
     &                    pxa(iw),rt,pres,cfore,hfore,u0,v0)
            call bootfep (nb,ib,pb1(j),pb0(j),pvb(j),wb(j),hwb(j),
     &                    pxb(iw),rt,pres,cback,hback,u1,v1)
            cback = -cback
            ptry(iw) = ptry(iw) + 1
            psum(1,iw) = psum(1,iw) + cfore
            psum(2,iw) = psum(2,iw) + cfore*cfore
            psum(3,iw) = psum(3,iw) + cback
            psum(4,iw) = psum(4,iw) + cback*cback
            val(1) = val(1) + cfore
            val(2) = val(2) + cback
            call bootiter (na,ia,nb,ib,ea(k),eb(j),pcref(iw),rt,rt,
     &                     pbar(1,iw),maxiter,ctrial,tconv)
            if (tconv) then
               call boothbar (na,ia,nb,ib,pa0(k),pa1(k),pb0(j),pb1(j),
     &                        pva(k),pvb(j),ea(k),eb(j),pcref(iw),rt,
     &                        rt,pres,ctrial,htrial)
               pconv(iw) = pconv(iw) + 1
               psum(5,iw) = psum(5,iw) + ctrial
               psum(6,iw) = psum(6,iw) + ctrial*ctrial
               psum(7,iw) = psum(7,iw) + htrial
               psum(8,iw) = psum(8,iw) + htrial*htrial
               val(3) = val(3) + ctrial
               val(4) = val(4) + htrial
            else
               allc = .false.
            end if
  230       continue
         end do
c
c     the direct enthalpy uses the same samples of the end windows
c
         if (wsel(1).gt.0 .and. wsel(nwin).gt.0) then
            u0 = 0.0d0
            do m = 1, wsel(1)
               i = wsamp(woff(1)+draw(woff(1)+m))
               u0 = u0 + rene(sfirst(i)) + pres*rvol(sfirst(i))
            end do
            u1 = 0.0d0
            do m = 1, wsel(nwin)
               i = wsamp(woff(nwin)+draw(woff(nwin)+m))
               u1 = u1 + rene(sfirst(i)) + pres*rvol(sfirst(i))
            end do
            val(5) = u1/dble(wsel(nwin)) - u0/dble(wsel(1))
            tsum(5) = tsum(5) + val(5)
            tsum2(5) = tsum2(5) + val(5)*val(5)
         end if
         if (allf) then
            ntf = ntf + 1
            do k = 1, 2
               tsum(k) = tsum(k) + val(k)
               tsum2(k) = tsum2(k) + val(k)*val(k)
            end do
         end if
         if (allc) then
            ntc = ntc + 1
            do k = 3, 4
               tsum(k) = tsum(k) + val(k)
               tsum2(k) = tsum2(k) + val(k)*val(k)
            end do
         end if
      end do
c
c     find the bootstrap errors of each pair
c
      do iw = 1, npair
         pfail(iw) = 0
         if (pstat(iw) .eq. 0) then
            pfail(iw) = nbst - pconv(iw)
            pfep(2,iw) = barsd (ptry(iw),psum(1,iw),psum(2,iw))
            pfep(4,iw) = barsd (ptry(iw),psum(3,iw),psum(4,iw))
            pbar(3,iw) = 0.0d0
            if (pconv(iw) .ge. 1) then
               pbar(3,iw) = psum(5,iw) / dble(pconv(iw))
            end if
            pbar(4,iw) = barsd (pconv(iw),psum(5,iw),psum(6,iw))
            phs(2,iw) = barsd (pconv(iw),psum(7,iw),psum(8,iw))
         end if
      end do
c
c     print information about the file and the samples used
c
      write (iout,240)
  240 format (/,' Free Energy from Tinker FEP File :',/)
      write (iout,250)  fepfile(1:trimtext(fepfile))
  250 format (' FEP File Name :',12x,a)
      write (iout,260)  nwin
  260 format (' Lambda Windows Found :',7x,i8)
      if (nwhead.gt.0 .and. nwhead.ne.nwin) then
         write (iout,270)  nwhead
  270    format (' Lambda Windows in Header :',3x,i8,
     &              '  (Mismatch)')
      end if
      write (iout,290)  temp
  290 format (' Simulation Temperature :',4x,f9.2,' K')
      if (usevol) then
         write (iout,300)  patm
  300    format (' Simulation Pressure :',7x,f9.2,' Atm')
      end if
      write (iout,340)  nbst
  340 format (' Bootstrap Trials :',11x,i8)
      if (norphan .gt. 0) then
         write (iout,350)  norphan
  350    format (' Records without Native Energy :',i6,
     &              '  (Ignored)')
      end if
c
c     print the free energy for each pair of adjacent windows
c
      write (iout,360)
  360 format (/,' Free Energy by Window Pair :',
     &        //,'  Pair  Lambda A  Lambda B    BAR Iter     +/-',
     &           '    BAR Boot     +/-    Cumul dG',/)
      cumul = 0.0d0
      found = .true.
      do iw = 1, npair
         line = ' '
         write (line(1:26),370)  iw,wlmda(iw),wlmda(iw+1)
  370    format (i6,2f10.4)
         if (pstat(iw) .lt. 0) then
            k = iw
            if (wsel(iw) .ne. 0)  k = iw + 1
            if (wsel(k) .eq. 0) then
               write (line(27:80),380)  k
  380          format (3x,'Not Computed, Window',i4,' has no Samples')
            else
               write (line(27:80),390)
  390          format (3x,'Not Computed, No Neighbor Energies')
            end if
            found = .false.
         else if (pstat(iw) .eq. 1) then
            write (line(27:80),400)
  400       format (3x,'BAR Iteration not Converged')
            found = .false.
         else if (pstat(iw) .eq. 2) then
            write (line(27:80),410)
  410       format (3x,'BAR Failed, No Overlap')
            found = .false.
         else
            write (line(27:46),420)  pbar(1,iw),pbar(2,iw)
  420       format (f12.4,f8.4)
            if (pconv(iw) .ge. 2) then
               write (line(47:66),430)  pbar(3,iw),pbar(4,iw)
  430          format (f12.4,f8.4)
            else
               write (line(47:66),440)
  440          format (9x,'n/a')
            end if
            if (found) then
               cumul = cumul + pbar(1,iw)
               write (line(67:78),450)  cumul
  450          format (f12.4)
            end if
         end if
         write (iout,460)  line(1:trimtext(line))
  460    format (a)
      end do
c
c     note any samples skipped or bootstrap trials that failed
c
      n = 0
      do iw = 1, npair
         if (pskip(iw).gt.0 .or. pfail(iw).gt.0)  n = n + 1
      end do
      if (n .gt. 0)  write (iout,470)
  470 format ()
      do iw = 1, npair
         if (pskip(iw) .gt. 0) then
            write (iout,480)  iw,pskip(iw)
  480       format (' Pair',i4,' :',i8,' Samples Skipped, Missing',
     &                 ' a Neighbor Energy')
         end if
         if (pfail(iw) .gt. 0) then
            write (iout,490)  iw,pfail(iw),nbst
  490       format (' Pair',i4,' :',i8,' of',i8,' Bootstrap Trials',
     &                 ' did not Converge')
         end if
      end do
c
c     print sampling diagnostics for each pair of windows
c
      write (iout,500)
  500 format (/,' Sampling Diagnostics by Window Pair :',
     &        //,'  Pair    nA    nB    FEP Fwd    FEP Bwd    Fwd-Bwd',
     &           '  SD dU Fwd  SD dU Bwd',/)
      do iw = 1, npair
         line = ' '
         write (line(1:18),510)  iw,pna(iw),pnb(iw)
  510    format (3i6)
         if (pstat(iw) .lt. 0) then
            write (line(19:80),520)
  520       format (3x,'Not Computed')
         else if (pstat(iw) .eq. 1) then
            write (line(19:80),530)
  530       format (3x,'BAR Iteration not Converged')
         else if (pstat(iw) .eq. 2) then
            write (line(19:80),540)
  540       format (3x,'BAR Failed, No Overlap')
         else
            write (line(19:73),550)  pfep(1,iw),pfep(3,iw),
     &                               pfep(1,iw)-pfep(3,iw),
     &                               psdu(1,iw),psdu(2,iw)
  550       format (5f11.4)
         end if
         write (iout,560)  line(1:trimtext(line))
  560    format (a)
      end do
c
c     totals need every pair and the expected number of windows
c
      totals = .true.
      if (nwhead.gt.0 .and. nwhead.ne.nwin)  totals = .false.
      nlist = 0
      do iw = 1, npair
         if (pstat(iw) .ne. 0) then
            totals = .false.
            nlist = nlist + 1
            plist(nlist) = iw
         end if
      end do
      if (.not. totals) then
         if (nlist .gt. 0) then
            write (iout,570)
  570       format (/,' Totals not Computed; Missing or Failed Pairs :')
            write (iout,580)  (plist(k),k=1,nlist)
  580       format (5x,15i5)
         else
            write (iout,590)
  590       format (/,' Totals not Computed; Lambda Windows Differ',
     &                 ' from the File Header')
         end if
      end if
c
c     sum the free energies over pairs, with errors from the joint
c     bootstrap trials over all of the windows
c
      if (totals) then
         do k = 1, 4
            sum(k) = 0.0d0
         end do
         do iw = 1, npair
            sum(1) = sum(1) + pfep(1,iw)
            sum(2) = sum(2) + pfep(3,iw)
            sum(3) = sum(3) + pbar(1,iw)
            sum(4) = sum(4) + phs(1,iw)
         end do
         write (iout,600)
  600    format (/,' Free Energy Totals :',/)
         write (iout,610)  sum(1),barsd(ntf,tsum(1),tsum2(1))
  610    format (' Free Energy via Forward FEP',9x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         write (iout,620)  sum(2),barsd(ntf,tsum(2),tsum2(2))
  620    format (' Free Energy via Backward FEP',8x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         if (ntc .ge. 2) then
            write (iout,630)  sum(3),barsd(ntc,tsum(3),tsum2(3))
  630       format (' Free Energy via BAR Iteration',7x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
            write (iout,640)  tsum(3)/dble(ntc),
     &                        barsd(ntc,tsum(3),tsum2(3))
  640       format (' Free Energy via BAR Bootstrap',7x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         else
            write (iout,650)  sum(3)
  650       format (' Free Energy via BAR Iteration',7x,f12.4,
     &                 ' Kcal/mol, Error not Available')
            write (iout,660)
  660       format (' Free Energy via BAR Bootstrap',13x,
     &                 'not Available')
         end if
c
c     find the enthalpy directly from the end window energies,
c     and sum the BAR enthalpy and entropy over the pairs
c
         u0 = 0.0d0
         do m = 1, wsel(1)
            i = wsamp(woff(1)+m)
            u0 = u0 + rene(sfirst(i)) + pres*rvol(sfirst(i))
         end do
         u1 = 0.0d0
         do m = 1, wsel(nwin)
            i = wsamp(woff(nwin)+m)
            u1 = u1 + rene(sfirst(i)) + pres*rvol(sfirst(i))
         end do
         hdir = u1/dble(wsel(nwin)) - u0/dble(wsel(1))
         write (iout,670)
  670    format (/,' Enthalpy and Entropy Totals :',/)
         write (iout,680)  hdir,barsd(nbst,tsum(5),tsum2(5))
  680    format (' Enthalpy via Direct Estimate',8x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         if (ntc .ge. 2) then
            write (iout,690)  sum(4),barsd(ntc,tsum(4),tsum2(4))
  690       format (' Enthalpy via BAR Estimate',11x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         else
            write (iout,700)  sum(4)
  700       format (' Enthalpy via BAR Estimate',11x,f12.4,
     &                 ' Kcal/mol, Error not Available')
         end if
         tds = 0.0d0
         do iw = 1, npair
            tds = tds + phs(3,iw)
         end do
         write (iout,710)  tds
  710    format (' BAR Estimate -T*dS Value',12x,f12.4,' Kcal/mol')
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (rstep)
      deallocate (rsim)
      deallocate (rtrl)
      deallocate (rene)
      deallocate (rvol)
      deallocate (sfirst)
      deallocate (slast)
      deallocate (swin)
      deallocate (wlmda)
      deallocate (wsel)
      deallocate (wseen)
      deallocate (woff)
      deallocate (wsamp)
      deallocate (draw)
      deallocate (mapa)
      deallocate (mapb)
      deallocate (pa0)
      deallocate (pa1)
      deallocate (pva)
      deallocate (pb0)
      deallocate (pb1)
      deallocate (pvb)
      deallocate (wa)
      deallocate (hwa)
      deallocate (ea)
      deallocate (wb)
      deallocate (hwb)
      deallocate (eb)
      deallocate (ia)
      deallocate (ib)
      deallocate (pstat)
      deallocate (pna)
      deallocate (pnb)
      deallocate (paoff)
      deallocate (pboff)
      deallocate (pskip)
      deallocate (pfail)
      deallocate (pconv)
      deallocate (ptry)
      deallocate (plist)
      deallocate (pfep)
      deallocate (pbar)
      deallocate (psdu)
      deallocate (phs)
      deallocate (pxa)
      deallocate (pxb)
      deallocate (pcref)
      deallocate (psum)
      return
      end
