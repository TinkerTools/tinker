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
c     and both trajectories must be at the same temperature; a trial
c     lambda energy file written by dynamics over a schedule of lambda
c     windows is analyzed by the "analyzelmda" program instead
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
         if (mode.ge.1 .and. mode.le.2)  query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' The Tinker Thermodynamic Perturbation Utility',
     &              ' Can :',
     &           //,4x,'(1) Create BAR File with Perturbed Potential',
     &              ' Energies',
     &           /,4x,'(2) Compute Thermodynamic Values from Tinker',
     &              ' BAR File')
         do while (mode.lt.1 .or. mode.gt.2)
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
