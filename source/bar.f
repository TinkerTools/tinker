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
c     and both trajectories must be at the same temperature; in
c     mode 3, a trial lambda energy file written by dynamics over a
c     schedule of lambda windows is read, and the FEP and BAR free
c     energies are found for each pair of adjacent windows and summed
c     over the whole schedule
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
c     check for log with energies of trajectory A in state 0
c
      use_log = .false.
      if (.not. recompute) then
         logfile = filea(1:lenga)
         call suffix (logfile,'log','old')
         inquire (file=logfile,exist=exist)
         if  (exist) then
            use_log = .true.
            ilog = freeunit ()
            open (unit=ilog,file=logfile,status='old')
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
      i = 0
      do while (.not. abort)
         i = i + 1
         if (use_log) then
            abort = .true.
            do while (abort)
               read (ilog,130,err=150,end=150)  record
  130          format (a240)
               if (record(1:18) .eq. ' Current Potential') then
                  string = record(19:240)
                  read (string,*,err=140,end=140)  ua0(i)
                  abort = .false.
               end if
  140          continue
            end do
  150       continue
            if (abort)  i = i - 1
         else
            call cutoffs
            ua0(i) = energy ()
            call readcart (iarc,first)
         end if
         imod = mod(i,100)
         if ((.not.abort.and.imod.eq.0) .or. (abort.and.imod.ne.0)) then
            write (iout,160)  i
  160       format (7x,'Completed',i8,' Coordinate Frames')
            flush (iout)
         end if
         if (i .ge. maxframe)  abort = .true.
      end do
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
         if (verbose) then
            write (iout,180)  i,ua0(i),ua1(i),ua1(i)-ua0(i)
  180       format (i11,2x,3f16.4)
         end if
         call readcart (iarc,first)
         if (i .ge. maxframe)  abort = .true.
      end do
      nfrma = i
      close (unit=iarc)
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
c     check for log with energies of trajectory B in state 1
c
      use_log = .false.
      if (.not. recompute) then
         logfile = fileb(1:lengb)
         call suffix (logfile,'log','old')
         inquire (file=logfile,exist=exist)
         if (exist) then
            use_log = .true.
            ilog = freeunit ()
            open (unit=ilog,file=logfile,status='old')
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
      i = 0
      do while (.not. abort)
         i = i + 1
         if (use_log) then
            abort = .true.
            do while (abort)
               read (ilog,240,err=260,end=260)  record
  240          format (a240)
               if (record(1:18) .eq. ' Current Potential') then
                  string = record(19:240)
                  read (string,*,err=250,end=250)  ub1(i)
                  abort = .false.
               end if
  250          continue
            end do
  260       continue
            if (abort)  i = i - 1
         else
            call cutoffs
            ub1(i) = energy ()
            call readcart (iarc,first)
         end if
         imod = mod(i,100)
         if ((.not.abort.and.imod.eq.0) .or. (abort.and.imod.ne.0)) then
            write (iout,270)  i
  270       format (7x,'Completed',i8,' Coordinate Frames')
            flush (iout)
         end if
         if (i .ge. maxframe)  abort = .true.
      end do
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
         if (verbose) then
            write (iout,290)  i,ub0(i),ub1(i),ub0(i)-ub1(i)
  290       format (i11,2x,3f16.4)
         end if
         call readcart (iarc,first)
         if (i .ge. maxframe)  abort = .true.
      end do
      nfrmb = i
      close (unit=iarc)
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
      integer size,next
      integer ltitlea,ltitleb
      integer nfrma,nfrmb
      integer nfrm,nbst
      integer starta,startb
      integer stopa,stopb
      integer stepa,stepb
      integer maxa,maxb,istat
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
      tempa = 0.0d0
      next = 1
      read (ibar,110,err=250,end=250)  record
  110 format (a240)
      size = trimtext (record)
      call gettext (record,string,next)
      read (string,*,err=250,end=250)  nfrma
      call gettext (record,string,next)
      read (string,*,err=250,end=250)  tempa
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
         read (record,*,err=140,end=140)  k,ua0(j),ua1(j),vola(j)
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
      size = trimtext (record)
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
      tempb = 0.0d0
      next = 1
      read (ibar,200,err=250,end=250)  record
  200 format (a240)
      size = trimtext (record)
      call gettext (record,string,next)
      read (string,*,err=250,end=250)  nfrmb
      call gettext (record,string,next)
      read (string,*,err=250,end=250)  tempb
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
         read (record,*,err=230,end=230)  k,ub0(j),ub1(j),volb(j)
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
      usevol = .false.
      do i = 1, nfrma
         if (vola(i) .ne. 0.0d0)  usevol = .true.
      end do
      do i = 1, nfrmb
         if (volb(i) .ne. 0.0d0)  usevol = .true.
      end do
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
      nbst = min(100000,nint(1.0d8/dble(nfrm)))
c
c     find the free energy, enthalpy and entropy via FEP and BAR
c
      call barpair (nfrma,nfrmb,ua0,ua1,ub0,ub1,vola,volb,tempa,
     &              tempb,patm,nbst,.true.,dgfep,dgbar,nboot,sdu,
     &              uave,hbar,istat)
      if (istat .eq. 1)  call fatal
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
c     from the trial energies of its neighbors; the native and trial
c     energies of a sample share the same periodic box, whose volume
c     is in a fifth column of newer files, and the volumes are used
c     for the PV enthalpy terms only if the file header also gives
c     the pressure of an isobaric simulation
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
      integer first,last,step
      integer lastw,iord
      integer nbst,nbstin
      integer nbmin,nbmax
      integer nsmin,nsmax
      integer na,nb,iw,istat
      integer nlist,nmiss
      integer freeunit
      integer trimtext
      integer nboot(2)
      integer, allocatable :: rstep(:)
      integer, allocatable :: sfirst(:)
      integer, allocatable :: slast(:)
      integer, allocatable :: swin(:)
      integer, allocatable :: wcount(:)
      integer, allocatable :: wsel(:)
      integer, allocatable :: wseen(:)
      integer, allocatable :: pstat(:)
      integer, allocatable :: pna(:)
      integer, allocatable :: pnb(:)
      integer, allocatable :: pskip(:)
      integer, allocatable :: pfail(:)
      integer, allocatable :: pconv(:)
      integer, allocatable :: pnbst(:)
      integer, allocatable :: plist(:)
      real*8 temp,tol,lmda
      real*8 patm,pres
      real*8 v0,v1
      real*8 lfirst,lsecond
      real*8 cumul,random
      real*8 bst,ratio
      real*8 u0,u1,hdir
      real*8 hsum,hsum2
      real*8 sum(4),sum2(4)
      real*8 dgfep(4),dgbar(4)
      real*8 sdu(2),uave(4)
      real*8 hbar(3)
      real*8, allocatable :: rsim(:)
      real*8, allocatable :: rtrl(:)
      real*8, allocatable :: rene(:)
      real*8, allocatable :: rvol(:)
      real*8, allocatable :: wlmda(:)
      real*8, allocatable :: ua0(:)
      real*8, allocatable :: ua1(:)
      real*8, allocatable :: ub0(:)
      real*8, allocatable :: ub1(:)
      real*8, allocatable :: vola(:)
      real*8, allocatable :: volb(:)
      real*8, allocatable :: pfep(:,:)
      real*8, allocatable :: pbar(:,:)
      real*8, allocatable :: psdu(:,:)
      real*8, allocatable :: phs(:,:)
      logical exist,query
      logical header,found
      logical down,totals
      logical hasvol,usevol
      logical, allocatable :: ssel(:)
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
c     set the samples to be used within each lambda window
c
      first = 0
      last = 0
      step = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=30,end=30)  first
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  last
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  step
   30 continue
      if (query) then
         write (iout,40)
   40    format (/,' First & Last Sample and Step Increment',
     &              ' within Each Window :  ',$)
         read (input,50)  record
   50    format (a240)
         read (record,*,err=60,end=60)  first,last,step
   60    continue
      end if
      if (first .le. 0)  first = 1
      if (last .lt. 0)  last = 0
      if (step .le. 0)  step = 1
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
      allocate (ssel(n))
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
c     the direction of the schedule is set by the first two windows
c     with samples, since the windows must be strictly monotonic
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
      down = (lsecond .lt. lfirst)
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
c     sort the windows into the order of the schedule
c
      do i = 2, nwin
         lmda = wlmda(i)
         j = i - 1
         do while (j .ge. 1)
            if ((.not.down .and. wlmda(j).gt.lmda) .or.
     &          (down .and. wlmda(j).lt.lmda)) then
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
c     assign each sample to its window and apply the selection
c
      allocate (wcount(nwin))
      allocate (wsel(nwin))
      allocate (wseen(nwin))
      do i = 1, nwin
         wcount(i) = 0
         wsel(i) = 0
         wseen(i) = 0
      end do
      do i = 1, nsamp
         lmda = rsim(sfirst(i))
         do k = 1, nwin
            if (abs(wlmda(k)-lmda) .lt. tol)  swin(i) = k
         end do
         wcount(swin(i)) = wcount(swin(i)) + 1
      end do
      do i = 1, nsamp
         iw = swin(i)
         wseen(iw) = wseen(iw) + 1
         iord = wseen(iw)
         lastw = wcount(iw)
         if (last .gt. 0)  lastw = min(last,wcount(iw))
         ssel(i) = (iord.ge.first .and. iord.le.lastw .and.
     &                 mod(iord-first,step).eq.0)
         if (ssel(i))  wsel(iw) = wsel(iw) + 1
      end do
      nsmin = wsel(1)
      nsmax = wsel(1)
      do i = 2, nwin
         nsmin = min(nsmin,wsel(i))
         nsmax = max(nsmax,wsel(i))
      end do
c
c     perform dynamic allocation of some local arrays
c
      npair = nwin - 1
      n = max(1,nsmax)
      allocate (ua0(n))
      allocate (ua1(n))
      allocate (ub0(n))
      allocate (ub1(n))
      allocate (vola(n))
      allocate (volb(n))
      allocate (pstat(npair))
      allocate (pna(npair))
      allocate (pnb(npair))
      allocate (pskip(npair))
      allocate (pfail(npair))
      allocate (pconv(npair))
      allocate (pnbst(npair))
      allocate (plist(npair))
      allocate (pfep(4,npair))
      allocate (pbar(4,npair))
      allocate (psdu(2,npair))
      allocate (phs(3,npair))
      do i = 1, n
         vola(i) = 0.0d0
         volb(i) = 0.0d0
      end do
c
c     find the free energy between each pair of adjacent windows
c
      nbmin = 0
      nbmax = 0
      do iw = 1, npair
         na = 0
         nb = 0
         pskip(iw) = 0
         do i = 1, nsamp
            if (ssel(i) .and. swin(i).eq.iw) then
               found = .false.
               do k = sfirst(i)+1, slast(i)
                  if (abs(rtrl(k)-wlmda(iw+1)) .lt. tol) then
                     if (.not. found) then
                        na = na + 1
                        ua0(na) = rene(sfirst(i))
                        ua1(na) = rene(k)
                        if (usevol)  vola(na) = rvol(sfirst(i))
                     end if
                     found = .true.
                  end if
               end do
               if (.not. found)  pskip(iw) = pskip(iw) + 1
            else if (ssel(i) .and. swin(i).eq.iw+1) then
               found = .false.
               do k = sfirst(i)+1, slast(i)
                  if (abs(rtrl(k)-wlmda(iw)) .lt. tol) then
                     if (.not. found) then
                        nb = nb + 1
                        ub1(nb) = rene(sfirst(i))
                        ub0(nb) = rene(k)
                        if (usevol)  volb(nb) = rvol(sfirst(i))
                     end if
                     found = .true.
                  end if
               end do
               if (.not. found)  pskip(iw) = pskip(iw) + 1
            end if
         end do
         pna(iw) = na
         pnb(iw) = nb
         pfail(iw) = 0
         pconv(iw) = 0
         pnbst(iw) = 0
         if (na.eq.0 .or. nb.eq.0) then
            pstat(iw) = -1
         else
            nbst = nbstin
            if (nbst .eq. 0) then
               nbst = min(10000,nint(1.0d8/dble(max(na,nb))))
               nbst = max(2,nbst)
            end if
            call barpair (na,nb,ua0,ua1,ub0,ub1,vola,volb,temp,temp,
     &                    patm,nbst,.false.,dgfep,dgbar,nboot,sdu,
     &                    uave,hbar,istat)
            pstat(iw) = istat
            pnbst(iw) = nbst
            pconv(iw) = nboot(1)
            pfail(iw) = nboot(2)
            do k = 1, 4
               pfep(k,iw) = dgfep(k)
               pbar(k,iw) = dgbar(k)
            end do
            do k = 1, 2
               psdu(k,iw) = sdu(k)
            end do
            do k = 1, 3
               phs(k,iw) = hbar(k)
            end do
            if (nbmin .eq. 0)  nbmin = nbst
            nbmin = min(nbmin,nbst)
            nbmax = max(nbmax,nbst)
         end if
      end do
c
c     print information about the file and the samples used
c
      write (iout,220)
  220 format (/,' Free Energy from Tinker FEP File :',/)
      write (iout,230)  fepfile(1:trimtext(fepfile))
  230 format (' FEP File Name :',12x,a)
      write (iout,240)  nwin
  240 format (' Lambda Windows Found :',7x,i8)
      if (nwhead.gt.0 .and. nwhead.ne.nwin) then
         write (iout,250)  nwhead
  250    format (' Lambda Windows in Header :',3x,i8,
     &              '  (Mismatch)')
      end if
      write (iout,260)  wlmda(1),wlmda(nwin)
  260 format (' Lambda Schedule :',11x,f9.4,' -->',f9.4)
      write (iout,270)  temp
  270 format (' Simulation Temperature :',4x,f9.2,' K')
      if (usevol) then
         write (iout,275)  patm
  275    format (' Simulation Pressure :',7x,f9.2,' Atm')
      end if
      write (iout,280)  nsmin,nsmax
  280 format (' Samples per Window :',8x,i8,' to',i8)
      if (last .eq. 0) then
         write (iout,290)  first,step
  290    format (' Sample Selection :',10x,'First',i6,'  Last   All',
     &              '  Step',i6)
      else
         write (iout,300)  first,last,step
  300    format (' Sample Selection :',10x,'First',i6,'  Last',i6,
     &              '  Step',i6)
      end if
      if (nbmin .eq. nbmax) then
         write (iout,310)  nbmax
  310    format (' Bootstrap Trials per Pair :',1x,i8)
      else
         write (iout,320)  nbmin,nbmax
  320    format (' Bootstrap Trials per Pair :',1x,i8,' to',i8)
      end if
      if (norphan .gt. 0) then
         write (iout,330)  norphan
  330    format (' Records without a Native Energy :',i6,
     &              '  (Ignored)')
      end if
      write (iout,340)
  340 format (/,' Bootstrap errors assume independent samples,',
     &           ' which may not hold at the',
     &        /,' sampling interval used')
c
c     print the free energy for each pair of adjacent windows
c
      write (iout,350)
  350 format (/,' Free Energy by Window Pair :',
     &        //,'  Pair  Lambda A  Lambda B    BAR Iter     +/-',
     &           '    BAR Boot     +/-    Cumul dG',/)
      cumul = 0.0d0
      found = .true.
      do iw = 1, npair
         line = ' '
         write (line(1:26),360)  iw,wlmda(iw),wlmda(iw+1)
  360    format (i6,2f10.4)
         if (pstat(iw) .lt. 0) then
            k = iw
            if (wsel(iw) .ne. 0)  k = iw + 1
            if (wsel(k) .eq. 0) then
               write (line(27:80),370)  k
  370          format (3x,'Not Computed, Window',i4,' has no Samples')
            else
               write (line(27:80),380)
  380          format (3x,'Not Computed, No Neighbor Energies')
            end if
            found = .false.
         else if (pstat(iw) .eq. 1) then
            write (line(27:80),390)
  390       format (3x,'BAR Iteration not Converged')
            found = .false.
         else
            write (line(27:46),400)  pbar(1,iw),pbar(2,iw)
  400       format (f12.4,f8.4)
            if (pconv(iw) .ge. 2) then
               write (line(47:66),410)  pbar(3,iw),pbar(4,iw)
  410          format (f12.4,f8.4)
            else
               write (line(47:66),420)
  420          format (9x,'n/a')
            end if
            if (found) then
               cumul = cumul + pbar(1,iw)
               write (line(67:78),430)  cumul
  430          format (f12.4)
            end if
         end if
         write (iout,440)  line(1:trimtext(line))
  440    format (a)
      end do
c
c     note any samples skipped or bootstrap trials that failed
c
      n = 0
      do iw = 1, npair
         if (pskip(iw).gt.0 .or. pfail(iw).gt.0)  n = n + 1
      end do
      if (n .gt. 0)  write (iout,450)
  450 format ()
      do iw = 1, npair
         if (pskip(iw) .gt. 0) then
            write (iout,460)  iw,pskip(iw)
  460       format (' Pair',i4,' :',i8,' Samples Skipped, Missing',
     &                 ' a Neighbor Energy')
         end if
         if (pfail(iw) .gt. 0) then
            write (iout,470)  iw,pfail(iw),pnbst(iw)
  470       format (' Pair',i4,' :',i8,' of',i8,' Bootstrap Trials',
     &                 ' did not Converge')
         end if
      end do
c
c     print sampling diagnostics for each pair of windows
c
      write (iout,480)
  480 format (/,' Sampling Diagnostics by Window Pair :',
     &        //,'  Pair    nA    nB    FEP Fwd    FEP Bwd    Fwd-Bwd',
     &           '  SD dU Fwd  SD dU Bwd',/)
      do iw = 1, npair
         line = ' '
         write (line(1:18),490)  iw,pna(iw),pnb(iw)
  490    format (3i6)
         if (pstat(iw) .lt. 0) then
            write (line(19:80),500)
  500       format (3x,'Not Computed')
         else if (pstat(iw) .eq. 1) then
            write (line(19:80),510)
  510       format (3x,'BAR Iteration not Converged')
         else
            write (line(19:73),520)  pfep(1,iw),pfep(3,iw),
     &                               pfep(1,iw)-pfep(3,iw),
     &                               psdu(1,iw),psdu(2,iw)
  520       format (5f11.4)
         end if
         write (iout,530)  line(1:trimtext(line))
  530    format (a)
      end do
c
c     totals need every pair and the expected number of windows
c
      totals = .true.
      if (nwhead.gt.0 .and. nwhead.ne.nwin)  totals = .false.
      nlist = 0
      do iw = 1, npair
         if (pstat(iw).lt.0 .or. pstat(iw).eq.1) then
            totals = .false.
            nlist = nlist + 1
            plist(nlist) = iw
         end if
      end do
      if (.not. totals) then
         write (iout,540)
  540    format (/,' Totals not Computed; Missing or Failed Pairs :')
         if (nlist .gt. 0) then
            write (iout,550)  (plist(k),k=1,nlist)
  550       format (5x,15i5)
         end if
         if (nwhead.gt.0 .and. nwhead.ne.nwin) then
            write (iout,560)  nwhead,nwin
  560       format (' File Header gives',i6,' Lambda Windows, but',
     &                 i6,' are Found')
         end if
      end if
c
c     sum the free energies over pairs, with errors in quadrature
c
      if (totals) then
         nmiss = 0
         do k = 1, 4
            sum(k) = 0.0d0
            sum2(k) = 0.0d0
         end do
         do iw = 1, npair
            sum(1) = sum(1) + pfep(1,iw)
            sum2(1) = sum2(1) + pfep(2,iw)**2
            sum(2) = sum(2) + pfep(3,iw)
            sum2(2) = sum2(2) + pfep(4,iw)**2
            sum(3) = sum(3) + pbar(1,iw)
            sum2(3) = sum2(3) + pbar(2,iw)**2
            sum(4) = sum(4) + pbar(3,iw)
            sum2(4) = sum2(4) + pbar(4,iw)**2
            if (pconv(iw) .lt. 2)  nmiss = nmiss + 1
         end do
         write (iout,570)
  570    format (/,' Free Energy Totals :',/)
         write (iout,580)  sum(1),sqrt(sum2(1))
  580    format (' Free Energy via Forward FEP',9x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         write (iout,590)  sum(2),sqrt(sum2(2))
  590    format (' Free Energy via Backward FEP',8x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         write (iout,600)  sum(3),sqrt(sum2(3))
  600    format (' Free Energy via BAR Iteration',7x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         if (nmiss .eq. 0) then
            write (iout,610)  sum(4),sqrt(sum2(4))
  610       format (' Free Energy via BAR Bootstrap',7x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         else
            write (iout,620)
  620       format (' Free Energy via BAR Bootstrap',13x,
     &                 'not Available')
         end if
      end if
c
c     find the enthalpy directly from the end window energies
c
      if (totals) then
         na = 0
         nb = 0
         do i = 1, nsamp
            if (ssel(i) .and. swin(i).eq.1) then
               na = na + 1
               ua0(na) = rene(sfirst(i))
               vola(na) = rvol(sfirst(i))
            else if (ssel(i) .and. swin(i).eq.nwin) then
               nb = nb + 1
               ub1(nb) = rene(sfirst(i))
               volb(nb) = rvol(sfirst(i))
            end if
         end do
c
c     the direct enthalpy uses all of the end window samples, with
c     its error from bootstrap resamples of energy and volume
c
         u0 = 0.0d0
         v0 = 0.0d0
         do i = 1, na
            u0 = u0 + ua0(i)
            v0 = v0 + vola(i)
         end do
         u1 = 0.0d0
         v1 = 0.0d0
         do i = 1, nb
            u1 = u1 + ub1(i)
            v1 = v1 + volb(i)
         end do
         hdir = (u1+pres*v1)/dble(nb) - (u0+pres*v0)/dble(na)
         nbst = nbstin
         if (nbst .eq. 0) then
            nbst = min(10000,nint(1.0d8/dble(max(na,nb))))
            nbst = max(2,nbst)
         end if
         bst = dble(nbst)
         ratio = bst / (bst-1.0d0)
         hsum = 0.0d0
         hsum2 = 0.0d0
         do k = 1, nbst
            u0 = 0.0d0
            do i = 1, na
               j = int(dble(na)*random()) + 1
               u0 = u0 + ua0(j) + pres*vola(j)
            end do
            u1 = 0.0d0
            do i = 1, nb
               j = int(dble(nb)*random()) + 1
               u1 = u1 + ub1(j) + pres*volb(j)
            end do
            u1 = u1/dble(nb) - u0/dble(na)
            hsum = hsum + u1
            hsum2 = hsum2 + u1*u1
         end do
         u1 = hsum / bst
         write (iout,630)
  630    format (/,' Enthalpy and Entropy Totals :',/)
         write (iout,640)  hdir,
     &                     sqrt(max(0.0d0,ratio*(hsum2/bst-u1*u1)))
  640    format (' Enthalpy via Direct Estimate',8x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
c
c     sum the BAR enthalpy and entropy over the pairs, where the
c     enthalpy error needs two converged bootstrap trials per pair
c
         n = 0
         do k = 1, 4
            sum(k) = 0.0d0
            sum2(k) = 0.0d0
         end do
         do iw = 1, npair
            if (pconv(iw) .lt. 2)  n = n + 1
            sum(1) = sum(1) + phs(1,iw)
            sum2(1) = sum2(1) + phs(2,iw)**2
            sum(2) = sum(2) + phs(3,iw)
         end do
         if (n .eq. 0) then
            write (iout,650)  sum(1),sqrt(sum2(1))
  650       format (' Enthalpy via BAR Estimate',11x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         else
            write (iout,660)  sum(1)
  660       format (' Enthalpy via BAR Estimate',11x,f12.4,
     &                 ' Kcal/mol, Error not Available')
         end if
         write (iout,670)  sum(2)
  670    format (' BAR Estimate -T*dS Value',12x,f12.4,' Kcal/mol')
         write (iout,680)
  680    format (/,' Enthalpy and Entropy Estimates are much Noisier',
     &              ' than the Free Energy')
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
      deallocate (ssel)
      deallocate (wlmda)
      deallocate (wcount)
      deallocate (wsel)
      deallocate (wseen)
      deallocate (ua0)
      deallocate (ua1)
      deallocate (ub0)
      deallocate (ub1)
      deallocate (vola)
      deallocate (volb)
      deallocate (pstat)
      deallocate (pna)
      deallocate (pnb)
      deallocate (pskip)
      deallocate (pfail)
      deallocate (pconv)
      deallocate (pnbst)
      deallocate (plist)
      deallocate (pfep)
      deallocate (pbar)
      deallocate (psdu)
      deallocate (phs)
      return
      end
