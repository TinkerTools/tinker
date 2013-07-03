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
c     current version takes as input the trajectory archives and key
c     files for state A, and similar for state B; then computes the
c     total potential energy for all frames of each trajectory under
c     control of both key files; finally the FEP and BAR algorithms
c     are used to compute the free energy for state A --> state B
c
c     literature reference:
c
c     C. H. Bennett, "Efficient Estimation of Free Energy Differences
c     from Monte Carlo Data", Journal of Computational Physics, 22,
c     245-268 (1976)
c
c
      program bar
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'units.i'
      integer i,j,k
      integer ixyz,nbst
      integer iter,maxiter
      integer nkey1,nfrm1
      integer nkey2,nfrm2
      integer start1,stop1,step1
      integer start2,stop2,step2
      integer maxframe,freeunit
      integer, allocatable :: bst1(:)
      integer, allocatable :: bst2(:)
      real*8 energy
      real*8 temp,rt
      real*8 delta,eps
      real*8 cold,cnew
      real*8 top,bot,rfrm
      real*8 sum,sum2
      real*8 mean,stdev
      real*8 random,ratio
      real*8 forward,backward
      real*8, allocatable :: ua0(:)
      real*8, allocatable :: ua1(:)
      real*8, allocatable :: ub0(:)
      real*8, allocatable :: ub1(:)
      logical exist,query,done
      character*120 record
      character*120 string
      character*120 xyzfile
      character*120 fname1
      character*120 fname2
      character*120, allocatable :: keys1(:)
      character*120, allocatable :: keys2(:)
c
c
c     perform dynamic allocation of some local arrays
c
      maxframe = 20000
      allocate (ua0(maxframe))
      allocate (ua1(maxframe))
      allocate (ub0(maxframe))
      allocate (ub1(maxframe))
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
   30    format (a120)
         read (record,*,err=40,end=40)  start1,stop1,step1
   40    continue
      end if
      if (start1 .eq. 0)  start1 = 1
      if (stop1 .eq. 0)  stop1 = maxframe
      if (step1 .eq. 0)  step1 = 1
c
c     perform dynamic allocation of some local arrays
c
      allocate (keys1(nkey))
c
c     store the filename and the keyword file for state A
c
      fname1 = filename
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
         read (string,*,err=50,end=50)  start2
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  stop2
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  step2
   50 continue
      if (query) then
         write (iout,60)
   60    format (/,' Numbers of First & Last File and Step',
     &              ' Increment :  ',$)
         read (input,70)  record
   70    format (a120)
         read (record,*,err=80,end=80)  start2,stop2,step2
   80    continue
      end if
      if (start2 .eq. 0)  start2 = 1
      if (stop2 .eq. 0)  stop2 = maxframe
      if (step2 .eq. 0)  step2 = 1
c
c     perform dynamic allocation of some local arrays
c
      allocate (keys2(nkey))
c
c     store the filename and the keyword file for state B
c
      fname2 = filename
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
      do i = 1, start1
         call readxyz (ixyz)
      end do
      nkey = nkey1
      do i = 1, nkey
         keyline(i) = keys1(i)
      end do
      call mechanic
c
c     find potential energies for trajectory A in state A
c
      write (iout,90)
   90 format (/,' Initial Processing for Trajectory A :',/)
      j = 0
      k = start1 - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         ua0(j) = energy ()
         do i = 1, step1
            call readxyz (ixyz)
         end do
         k = k + step1 - 1
         if (k .ge. stop1)  abort = .true.
         if (mod(j,100).eq.0 .or. abort) then
            write (iout,100)  j
  100       format (7x,'Completed',i7,' Coordinate Frames')
         end if
      end do
      nfrm1 = j
c
c     reset trajectory A and process the initial structure
c
      rewind (unit=ixyz)
      do i = 1, start1
         call readxyz (ixyz)
      end do
      nkey = nkey2
      do i = 1, nkey
         keyline(i) = keys2(i)
      end do
      call mechanic
c
c     find potential energies for trajectory A in state B
c
      write (iout,110)
  110 format (/,' Potential Energy Values for Trajectory A :',
     &        //,7x,'Frame',9x,'State A',9x,'State B',12x,'Delta',/)
      j = 0
      k = start1 - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         ua1(j) = energy ()
         write (iout,120)  k,ua0(j),ua1(j),ua1(j)-ua0(j)
  120    format (i11,2x,3f16.4)
         do i = 1, step1
            call readxyz (ixyz)
         end do
         k = k + step1 - 1
         if (k .ge. stop1)  abort = .true.
      end do
      nfrm1 = j
      close (unit=ixyz)
c
c     reopen trajectory B and process the initial structure
c
      ixyz = freeunit ()
      xyzfile = fname2
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      do i = 1, start2
         call readxyz (ixyz)
      end do
      nkey = nkey1
      do i = 1, nkey
         keyline(i) = keys1(i)
      end do
      call mechanic
c
c     find potential energies for trajectory B in state A
c
      write (iout,130)
  130 format (/,' Initial Processing for Trajectory B :',/)
      j = 0
      k = start2 - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         ub0(j) = energy ()
         do i = 1, step2
            call readxyz (ixyz)
         end do
         k = k + step2 - 1
         if (k .ge. stop2)  abort = .true.
         if (mod(j,100).eq.0 .or. abort) then
            write (iout,140)  j
  140       format (7x,'Completed',i7,' Coordinate Frames')
         end if
      end do
      nfrm2 = j
c
c     reset trajectory B and process the initial structure
c
      rewind (unit=ixyz)
      do i = 1, start2
         call readxyz (ixyz)
      end do
      nkey = nkey2
      do i = 1, nkey
         keyline(i) = keys2(i)
      end do
      call mechanic
c
c     find potential energies for trajectory B in state B
c
      write (iout,150)
  150 format (/,' Potential Energy Values for Trajectory B :',
     &        //,7x,'Frame',9x,'State A',9x,'State B',12x,'Delta',/)
      j = 0
      k = start2 - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         ub1(j) = energy ()
         write (iout,160)  k,ub0(j),ub1(j),ub0(j)-ub1(j)
  160    format (i11,2x,3f16.4)
         do i = 1, step2
            call readxyz (ixyz)
         end do
         k = k + step2 - 1
         if (k .ge. stop2)  abort = .true.
      end do
      nfrm2 = j
      close (unit=ixyz)
c
c     perform deallocation of some local arrays
c
      deallocate (keys1)
      deallocate (keys2)
c
c     set the frame ratio, temperature and Boltzmann factor
c
      rfrm = dble(nfrm1) / dble(nfrm2)
      temp = 298.0d0
      rt = gasconst * temp
c
c     compute the free energy difference via Zwanzig equation
c
      write (iout,170)
  170 format (/,' Estimation of Free Energy Difference',
     &           ' via FEP Method :',/)
      sum = 0.0d0
      do i = 1, nfrm1
         sum = sum + exp((ua0(i)-ua1(i))/rt)
      end do
      forward = -rt * log(sum/dble(nfrm1))
      sum = 0.0d0
      do i = 1, nfrm2
         sum = sum + exp((ub1(i)-ub0(i))/rt)
      end do
      backward = -rt * log(sum/dble(nfrm2))
      write (iout,180)  forward
  180 format (' FEP Forward Free Energy',4x,f12.4,' Kcal/mol')
      write (iout,190)  backward
  190 format (' FEP Backward Free Energy',3x,f12.4,' Kcal/mol')
c
c     compute the initial free energy via the BAR equation
c
      write (iout,200)
  200 format (/,' Estimation of Free Energy Difference',
     &           ' via BAR Method :',/)
      maxiter = 100
      eps = 0.0001d0
      done = .false.
      iter = 0
      cold = 0.0d0
      top = 0.0d0
      do i = 1, nfrm2
         top = top + 1.0d0/(1.0d0+exp((ub0(i)-ub1(i)+cold)/rt))
      end do
      bot = 0.0d0
      do i = 1, nfrm1
         bot = bot + 1.0d0/(1.0d0+exp((ua1(i)-ua0(i)-cold)/rt))
      end do
      cnew = rt*log(rfrm*top/bot) + cold
      delta = abs(cnew-cold)
      if (delta .lt. eps) then
         done = .true.
         write (iout,210)  cnew,iter
  210    format (' BAR Free Energy Estimate',3x,f12.4,
     &              ' Kcal/mol at',i4,' Iterations')
      else
         write (iout,220)  iter,cnew
  220    format (' BAR Iteration',i4,10x,f12.4,' Kcal/mol')
      end if
c
c     iterate the BAR equation to converge the free energy
c
      do while (.not. done)
         iter = iter + 1
         cold = cnew
         top = 0.0d0
         do i = 1, nfrm2
            top = top + 1.0d0/(1.0d0+exp((ub0(i)-ub1(i)+cold)/rt))
         end do
         bot = 0.0d0
         do i = 1, nfrm1
            bot = bot + 1.0d0/(1.0d0+exp((ua1(i)-ua0(i)-cold)/rt))
         end do
         cnew = rt*log(rfrm*top/bot) + cold
         delta = abs(cnew-cold)
         if (delta .lt. eps) then
            done = .true.
            write (iout,230)  cnew,iter
  230       format (/,' BAR Free Energy Estimate',3x,f12.4,
     &                 ' Kcal/mol at',i4,' Iterations')
         else
            write (iout,240)  iter,cnew
  240       format (' BAR Iteration',i4,10x,f12.4,' Kcal/mol')
         end if
         if (iter.ge.maxiter .and. .not.done) then
            done = .true.
            write (iout,250)  maxiter
  250       format (/,' BAR Free Energy Estimate not Converged',
     &                 ' after',i4,' Iterations')
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      nbst = 100000
      allocate (bst1(nbst))
      allocate (bst2(nbst))
c
c     use bootstrapping analysis to estimate statistical error
c
      write (iout,260)
  260 format (/,' Bootstrapping Error Analysis of BAR',
     &           ' Free Energy :',/)
      sum = 0.0d0
      sum2 = 0.0d0
      do k = 1, nbst
         done = .false.
         iter = 0
         cold = 0.0d0
         top = 0.0d0
         do i = 1, nfrm2
            bst2(i) = int(dble(nfrm2)*random()) + 1
            j = bst2(i)
            top = top + 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)+cold)/rt))
         end do
         bot = 0.0d0
         do i = 1, nfrm1
            bst1(i) = int(dble(nfrm1)*random()) + 1
            j = bst1(i)
            bot = bot + 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)-cold)/rt))
         end do
         cnew = rt*log(rfrm*top/bot) + cold
         delta = abs(cnew-cold)
         do while (.not. done)
            iter = iter + 1
            cold = cnew
            top = 0.0d0
            do i = 1, nfrm2
               j = bst2(i)
               top = top + 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)+cold)/rt))
            end do
            bot = 0.0d0
            do i = 1, nfrm1
               j = bst1(i)
               bot = bot + 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)-cold)/rt))
            end do
            cnew = rt*log(rfrm*top/bot) + cold
            delta = abs(cnew-cold)
            if (delta .lt. eps) then
               done = .true.
               sum = sum + cnew
               sum2 = sum2 + cnew*cnew
               if (debug) then
                  write (iout,270)  k,cnew,iter
  270             format (' Bootstrap Estimate',i7,2x,f12.4,
     &                       ' Kcal/mol at',i4,' Resamples')
               end if
            end if
         end do
      end do
      mean = sum / dble(nbst)
      ratio = dble(nbst/(nbst-1))
      stdev = sqrt(ratio*(sum2/dble(nbst)-mean*mean))
      if (verbose) then
         write (iout,280)
  280    format ()
      end if
      write (iout,290)  mean,stdev
  290 format (' BAR Bootstrap Free Energy',2x,f12.4,
     &           ' +/-',f8.4,' Kcal/mol')
c
c     perform deallocation of some local arrays
c
      deallocate (ua0)
      deallocate (ua1)
      deallocate (ub0)
      deallocate (ub1)
      deallocate (bst1)
      deallocate (bst2)
c
c     perform any final tasks before program exit
c
      call final
      end
