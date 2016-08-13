c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2013  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program bar  --  free energy differences via FEP and BAR  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "bar" computes the free energy difference between two states
c     via the Zwanzig free energy perturbation (FEP) and Bennett
c     acceptance ratio (BAR) methods
c
c     note current version takes as input a trajectory archive and
c     key file for state A, and similar for state B; then finds the
c     total potential energy for all frames of each trajectory under
c     control of both key files; finally the FEP and BAR algorithms
c     are used to compute the free energy for state A --> state B
c
c     modifications to handle NPT simulation data by Chengwen Liu,
c     University of Texas at Austin during October 2015
c
c     literature references:
c
c     C. H. Bennett, "Efficient Estimation of Free Energy Differences
c     from Monte Carlo Data", Journal of Computational Physics, 22,
c     245-268 (1976)
c
c     K. B. Daly, J. B. Benziger, P. G. Debenedetti and
c     A. Z. Panagiotopoulos, "Massively Parallel Chemical Potential
c     Calculation on Graphics Processing Units", Computer Physics
c     Communications, 183, 2054-2062 (2012)  [modification for NPT]
c
c
      program bar
      use sizes
      use atoms
      use boxes
      use energi
      use files
      use inform
      use iounit
      use keys
      use titles
      use units
      implicit none
      integer i,j,k
      integer n1,n2
      integer ixyz,ibar,nbst
      integer iter,maxiter
      integer leng1,leng2
      integer ltitle1,ltitle2
      integer nkey1,nkey2
      integer nfrm1,nfrm2
      integer start1,start2
      integer stop1,stop2
      integer step1,step2
      integer maxframe
      integer freeunit
      integer trimtext
      integer, allocatable :: bst1(:)
      integer, allocatable :: bst2(:)
      real*8 energy
      real*8 rt,rt1,rt2
      real*8 delta,eps
      real*8 frm1,frm2
      real*8 temp1,temp2
      real*8 cold,cnew
      real*8 top,top2
      real*8 bot,bot2
      real*8 fterm,rfrm
      real*8 sum,sum2,vave
      real*8 mean,stdev
      real*8 random,ratio
      real*8 forward,backward
      real*8, allocatable :: ua0(:)
      real*8, allocatable :: ua1(:)
      real*8, allocatable :: ub0(:)
      real*8, allocatable :: ub1(:)
      real*8, allocatable :: vola(:)
      real*8, allocatable :: volb(:)
      logical exist,query,done
      character*240 record
      character*240 string
      character*240 fname1
      character*240 fname2
      character*240 title1
      character*240 title2
      character*240 xyzfile
      character*240 barfile
      character*240, allocatable :: keys1(:)
      character*240, allocatable :: keys2(:)
c
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
c     get trajectory A archive and setup mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     set beginning and ending frame and the step increment
c
      start1 = 0
      stop1 = 0
      step1 = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  start1
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  stop1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  step1
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Numbers of First & Last File and Step',
     &              ' Increment :  ',$)
         read (input,30)  record
   30    format (a240)
         read (record,*,err=40,end=40)  start1,stop1,step1
   40    continue
      end if
      if (start1 .eq. 0)  start1 = 1
      if (stop1 .eq. 0)  stop1 = maxframe
      if (step1 .eq. 0)  step1 = 1
c
c     find temperature at which trajectory A was originally run
c
      temp1 = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  temp1
   50 continue
      do while (temp1 .lt. 0.0d0)
         write (iout,60)
   60    format (/,' Enter the Trajectory Temperature in Degrees',
     &              ' K [298] :  ',$)
         read (input,70,err=80)  temp1
   70    format (f20.0)
         if (temp1 .le. 0.0d0)  temp1 = 298.0d0
   80    continue
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (keys1(nkey))
c
c     store the filename and the keyword file for state A
c
      n1 = n
      fname1 = filename
      leng1 = leng
      title1 = title
      ltitle1 = ltitle
      nkey1 = nkey
      do i = 1, nkey1
         keys1(i) = keyline(i)
      end do
c
c     get trajectory B archive and setup mechanics calculation
c
      call getxyz
      call mechanic
      silent = .true.
c
c     set beginning and ending frame and the step increment
c
      start2 = 0
      stop2 = 0
      step2 = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=90,end=90)  start2
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  stop2
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  step2
   90 continue
      if (query) then
         write (iout,100)
  100    format (/,' Numbers of First & Last File and Step',
     &              ' Increment :  ',$)
         read (input,110)  record
  110    format (a240)
         read (record,*,err=120,end=120)  start2,stop2,step2
  120    continue
      end if
      if (start2 .eq. 0)  start2 = 1
      if (stop2 .eq. 0)  stop2 = maxframe
      if (step2 .eq. 0)  step2 = 1
c
c     find temperature at which trajectory B was originally run
c
      temp2 = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=130,end=130)  temp2
  130 continue
      do while (temp2 .lt. 0.0d0)
         write (iout,140)
  140    format (/,' Enter the Trajectory Temperature in Degrees',
     &              ' K [298] :  ',$)
         read (input,150,err=160)  temp2
  150    format (f20.0)
         if (temp2 .le. 0.0d0)  temp2 = 298.0d0
  160    continue
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (keys2(nkey))
c
c     store the filename and the keyword file for state B
c
      n2 = n
      fname2 = filename
      leng2 = leng
      title2 = title
      ltitle2 = ltitle
      nkey2 = nkey
      do i = 1, nkey2
         keys2(i) = keyline(i)
      end do
c
c     reopen trajectory A and process the initial structure
c
      ixyz = freeunit ()
      xyzfile = fname1
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      do i = 1, start1-1
         call readxyz (ixyz)
      end do
      call readxyz (ixyz)
      nkey = nkey1
      do i = 1, nkey
         keyline(i) = keys1(i)
      end do
      if (abort) then
         write (iout,170)
  170    format (/,' BAR  --  No Coordinate Frames Available',
     &              ' from First Input File')
         call fatal
      end if
      call mechanic
c
c     find potential energies for trajectory A in state A
c
      write (iout,180)
  180 format (/,' Initial Processing for Trajectory A :',/)
      j = 0
      k = start1 - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         call cutoffs
         ua0(j) = energy ()
         vola(j) = volbox
         do i = 1, step1-1
            call readxyz (ixyz)
         end do
         call readxyz (ixyz)
         k = k + step1 - 1
         if (k .ge. stop1)  abort = .true.
         if (mod(j,100).eq.0 .or. abort) then
            write (iout,190)  j
  190       format (7x,'Completed',i8,' Coordinate Frames')
            flush (iout)
         end if
      end do
c
c     reset trajectory A and process the initial structure
c
      rewind (unit=ixyz)
      do i = 1, start1-1
         call readxyz (ixyz)
      end do
      call readxyz (ixyz)
      nkey = nkey2
      do i = 1, nkey
         keyline(i) = keys2(i)
      end do
      call mechanic
c
c     find potential energies for trajectory A in state B
c
      if (verbose) then
         write (iout,200)
  200    format (/,' Potential Energy Values for Trajectory A :',
     &           //,7x,'Frame',9x,'State A',9x,'State B',12x,'Delta',/)
      end if
      j = 0
      k = start1 - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         call cutoffs
         ua1(j) = energy ()
         if (verbose) then
            write (iout,210)  k,ua0(j),ua1(j),ua1(j)-ua0(j)
  210       format (i11,2x,3f16.4)
         end if
         do i = 1, step1-1
            call readxyz (ixyz)
         end do
         call readxyz (ixyz)
         k = k + step1 - 1
         if (k .ge. stop1)  abort = .true.
      end do
      nfrm1 = j
      frm1 = dble(nfrm1)
      close (unit=ixyz)
c
c     reopen trajectory B and process the initial structure
c
      ixyz = freeunit ()
      xyzfile = fname2
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      do i = 1, start2-1
         call readxyz (ixyz)
      end do
      call readxyz (ixyz)
      nkey = nkey1
      do i = 1, nkey
         keyline(i) = keys1(i)
      end do
      if (abort) then
         write (iout,220)
  220    format (/,' BAR  --  No Coordinate Frames Available',
     &              ' from Second Input File')
         call fatal
      end if
      call mechanic
c
c     find potential energies for trajectory B in state A
c
      write (iout,230)
  230 format (/,' Initial Processing for Trajectory B :',/)
      j = 0
      k = start2 - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         call cutoffs
         ub0(j) = energy ()
         volb(j) = volbox
         do i = 1, step2-1
            call readxyz (ixyz)
         end do
         call readxyz (ixyz)
         k = k + step2 - 1
         if (k .ge. stop2)  abort = .true.
         if (mod(j,100).eq.0 .or. abort) then
            write (iout,240)  j
  240       format (7x,'Completed',i8,' Coordinate Frames')
            flush (iout)
         end if
      end do
c
c     reset trajectory B and process the initial structure
c
      rewind (unit=ixyz)
      do i = 1, start2-1
         call readxyz (ixyz)
      end do
      call readxyz (ixyz)
      nkey = nkey2
      do i = 1, nkey
         keyline(i) = keys2(i)
      end do
      call mechanic
c
c     find potential energies for trajectory B in state B
c
      if (verbose) then
         write (iout,250)
  250    format (/,' Potential Energy Values for Trajectory B :',
     &           //,7x,'Frame',9x,'State A',9x,'State B',12x,'Delta',/)
      end if
      j = 0
      k = start2 - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         call cutoffs
         ub1(j) = energy ()
         if (verbose) then
            write (iout,260)  k,ub0(j),ub1(j),ub0(j)-ub1(j)
  260       format (i11,2x,3f16.4)
         end if
         do i = 1, step2-1
            call readxyz (ixyz)
         end do
         call readxyz (ixyz)
         k = k + step2 - 1
         if (k .ge. stop2)  abort = .true.
      end do
      nfrm2 = j
      frm2 = dble(nfrm2)
      close (unit=ixyz)
c
c     perform deallocation of some local arrays
c
      deallocate (keys1)
      deallocate (keys2)
c
c     set the frame ratio, temperature and Boltzmann factor
c
      rfrm = frm1 / frm2
      rt1 = gasconst * temp1
      rt2 = gasconst * temp2
      rt = 0.5d0 * (rt1+rt2)
c
c     compute the volume corrections for trajectories A and B
c
      vave = 0.0d0
      do i = 1, nfrm1
         vave = vave + vola(i)
      end do
      vave = vave / frm1
      if (vave .ne. 0.0d0) then
         do i = 1, nfrm1
            if (vola(i) .ne. 0.0d0)  vola(i) = -rt1 * log(vola(i)/vave)
         end do
      end if
      vave = 0.0d0
      do i = 1, nfrm2
         vave = vave + volb(i)
      end do
      vave = vave / frm2
      if (vave .ne. 0.0d0) then
         do i = 1, nfrm2
            if (volb(i) .ne. 0.0d0)  volb(i) = -rt2 * log(volb(i)/vave)
         end do
      end if
c
c     compute the free energy difference via Zwanzig equation
c
      write (iout,270)
  270 format (/,' Estimation of Free Energy Difference',
     &           ' via FEP Method :',/)
      sum = 0.0d0
      do i = 1, nfrm1
         sum = sum + exp((ua0(i)-ua1(i)+vola(i))/rt1)
      end do
      forward = -rt1 * log(sum/frm1)
      sum = 0.0d0
      do i = 1, nfrm2
         sum = sum + exp((ub1(i)-ub0(i)+volb(i))/rt2)
      end do
      backward = -rt2 * log(sum/frm2)
      write (iout,280)  forward
  280 format (' FEP Forward Free Energy',4x,f12.4,' Kcal/mol')
      write (iout,290)  backward
  290 format (' FEP Backward Free Energy',3x,f12.4,' Kcal/mol')
c
c     compute the initial free energy via the BAR equation
c
      write (iout,300)
  300 format (/,' Estimation of Free Energy Difference',
     &           ' via BAR Method :',/)
      maxiter = 100
      eps = 0.0001d0
      done = .false.
      iter = 0
      cold = 0.0d0
      top = 0.0d0
      top2 = 0.0d0
      do i = 1, nfrm2
         fterm = 1.0d0 / (1.0d0+exp((ub0(i)-ub1(i)+volb(i)+cold)/rt2))
         top = top + fterm
         top2 = top2 + fterm*fterm
      end do
      bot = 0.0d0
      bot2 = 0.0d0
      do i = 1, nfrm1
         fterm = 1.0d0 / (1.0d0+exp((ua1(i)-ua0(i)+vola(i)-cold)/rt1))
         bot = bot + fterm
         bot2 = bot2 + fterm*fterm
      end do
      cnew = rt*log(rfrm*top/bot) + cold
      stdev = sqrt((bot2-bot*bot/frm1)/(bot*bot)
     &                + (top2-top*top/frm2)/(top*top))
      delta = abs(cnew-cold)
      write (iout,310)  iter,cnew
  310 format (' BAR Iteration',i4,10x,f12.4,' Kcal/mol')
      if (delta .lt. eps) then
         done = .true.
         write (iout,320)  cnew,stdev
  320    format (' BAR Free Energy Estimate',3x,f12.4,
     &              ' +/-',f8.4,' Kcal/mol')
      end if
c
c     iterate the BAR equation to converge the free energy
c
      do while (.not. done)
         iter = iter + 1
         cold = cnew
         top = 0.0d0
         top2 = 0.0d0
         do i = 1, nfrm2
            fterm = 1.0d0 / (1.0d0+exp((ub0(i)-ub1(i)+volb(i)
     &                                     +cold)/rt2))
            top = top + fterm
            top2 = top2 + fterm*fterm
         end do
         bot = 0.0d0
         bot2 = 0.0d0
         do i = 1, nfrm1
            fterm = 1.0d0 / (1.0d0+exp((ua1(i)-ua0(i)+vola(i)
     &                                     -cold)/rt1))
            bot = bot + fterm
            bot2 = bot2 + fterm*fterm
         end do
         cnew = rt*log(rfrm*top/bot) + cold
         stdev = sqrt((bot2-bot*bot/frm1)/(bot*bot)
     &                   + (top2-top*top/frm2)/(top*top))
         delta = abs(cnew-cold)
         write (iout,330)  iter,cnew
  330    format (' BAR Iteration',i4,10x,f12.4,' Kcal/mol')
         if (delta .lt. eps) then
            done = .true.
            write (iout,340)  cnew,stdev
  340       format (/,' BAR Free Energy Estimate',3x,f12.4,
     &                 ' +/-',f8.4,' Kcal/mol')
         end if
         if (iter.ge.maxiter .and. .not.done) then
            done = .true.
            write (iout,350)  maxiter
  350       format (/,' BAR Free Energy Estimate not Converged',
     &                 ' after',i4,' Iterations')
            call fatal
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      nbst = 100000
      allocate (bst1(nbst))
      allocate (bst2(nbst))
c
c     use bootstrap analysis to estimate statistical error
c
      write (iout,360)
  360 format (/,' Bootstrap Error Analysis for BAR',
     &           ' Free Energy :')
      if (debug) then
         write (iout,370)
  370    format ()
      end if
      sum = 0.0d0
      sum2 = 0.0d0
      do k = 1, nbst
         done = .false.
         iter = 0
         cold = 0.0d0
         top = 0.0d0
         do i = 1, nfrm2
            bst2(i) = int(frm2*random()) + 1
            j = bst2(i)
            top = top + 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)+volb(i)
     &                                       +cold)/rt2))
         end do
         bot = 0.0d0
         do i = 1, nfrm1
            bst1(i) = int(frm1*random()) + 1
            j = bst1(i)
            bot = bot + 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)+vola(i)
     &                                       -cold)/rt1))
         end do
         cnew = rt*log(rfrm*top/bot) + cold
         delta = abs(cnew-cold)
         do while (.not. done)
            iter = iter + 1
            cold = cnew
            top = 0.0d0
            do i = 1, nfrm2
               j = bst2(i)
               top = top + 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)+volb(i)
     &                                          +cold)/rt2))
            end do
            bot = 0.0d0
            do i = 1, nfrm1
               j = bst1(i)
               bot = bot + 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)+vola(i)
     &                                          -cold)/rt1))
            end do
            cnew = rt*log(rfrm*top/bot) + cold
            delta = abs(cnew-cold)
            if (delta .lt. eps) then
               done = .true.
               sum = sum + cnew
               sum2 = sum2 + cnew*cnew
               if (debug) then
                  write (iout,380)  k,cnew,iter
  380             format (' Bootstrap Estimate',i7,2x,f12.4,
     &                       ' Kcal/mol at',i4,' Resamples')
               end if
            end if
         end do
      end do
      mean = sum / dble(nbst)
      ratio = dble(nbst/(nbst-1))
      stdev = sqrt(ratio*(sum2/dble(nbst)-mean*mean))
      write (iout,390)  mean,stdev
  390 format (/,' BAR Bootstrap Free Energy',2x,f12.4,
     &           ' +/-',f8.4,' Kcal/mol')
c
c     save the energies and volume corrections to a file
c
      ibar = freeunit ()
      filename = fname1
      barfile = filename(1:leng1)//'.bar'
      call version (barfile,'new')
      open (unit=ibar,file=barfile,status ='new')
      write (ibar,400)  nfrm1,title1(1:ltitle1)
  400 format (i6,2x,a)
      do i = 1, nfrm1
         write (ibar,410)  i,ua0(i),ua1(i),vola(i)
  410    format (i6,2x,3f16.4)
      end do
      write (ibar,420)  nfrm2,title2(1:ltitle2)
  420 format (i6,2x,a)
      do i = 1, nfrm2
         write (ibar,430)  i,ub0(i),ub1(i),volb(i)
  430    format (i6,2x,3f16.4)
      end do
      close (unit=ibar)
c
c     perform deallocation of some local arrays
c
      deallocate (bst1)
      deallocate (bst2)
      deallocate (ua0)
      deallocate (ua1)
      deallocate (ub0)
      deallocate (ub1)
      deallocate (vola)
      deallocate (volb)
c
c     perform any final tasks before program exit
c
      call final
      end
