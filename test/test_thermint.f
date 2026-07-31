c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_thermint  --  TI window bookkeeping  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_thermint" checks the lambda window schedule and the
c     block averaging of the lambda derivative used by the
c     thermodynamic integration method
c
c     these are pure module logic tests with no molecular system;
c     "etidyn" samples the global "dedl", so the accumulation is
c     driven by assigning dedl directly instead of by evaluating
c     any energy
c
c
      subroutine test_thermint
      implicit none
      logical skiptest
      character*(*) tname
      parameter (tname='test_thermint')
c
c
      if (skiptest(tname,'thermint'))  return
      call initial
      call test_thermint_avgstd
      call test_thermint_schedule
      call test_thermint_settisched
      call test_thermint_data
      call test_thermint_inittidyn
      call test_thermint_etidyn
      call test_thermint_partialblock
      call test_thermint_trailing
      call clearti
      call final
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_thermint_avgstd  --  block mean and std  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_thermint_avgstd" checks the average and population
c     standard deviation kernel that reduces each block of saved
c     dU/dlambda values
c
c
      subroutine test_thermint_avgstd
      implicit none
      integer i
      real*8 avg,std
      real*8 sd10
      real*8 eps
      real*8 v10(10)
      real*8 vc(4)
      real*8 v1(1)
      real*8 vbig(10)
c
c
c     population standard deviation of ten consecutive integers
c
      eps = 1.0d-12
      sd10 = sqrt(8.25d0)
      do i = 1, 10
         v10(i) = dble(i)
      end do
      call avgstd (v10,1,10,avg,std)
      call assert_real (avg,5.5d0,eps,'avgstd ten sample average')
      call assert_real (std,sd10,eps,'avgstd ten sample deviation')
c
c     the count must be honored, so only the first five participate
c
      call avgstd (v10,1,5,avg,std)
      call assert_real (avg,3.0d0,eps,'avgstd partial count average')
      call assert_real (std,sqrt(2.0d0),eps,
     &                  'avgstd partial count deviation')
c
c     a constant list must give exactly zero, not a roundoff residue
c
      do i = 1, 4
         vc(i) = 7.0d0
      end do
      call avgstd (vc,1,4,avg,std)
      call assert_real (avg,7.0d0,eps,'avgstd constant list average')
      call assert_real (std,0.0d0,0.0d0,
     &                  'avgstd constant list deviation')
c
c     a single sample has no spread
c
      v1(1) = 42.0d0
      call avgstd (v1,1,1,avg,std)
      call assert_real (avg,42.0d0,eps,'avgstd single sample average')
      call assert_real (std,0.0d0,0.0d0,
     &                  'avgstd single sample deviation')
c
c     an empty range returns through the early exit leaving zeros
c
      avg = -1.0d0
      std = -1.0d0
      call avgstd (v10,1,0,avg,std)
      call assert_real (avg,0.0d0,0.0d0,'avgstd empty range average')
      call assert_real (std,0.0d0,0.0d0,'avgstd empty range deviation')
c
c     a large common offset must not swamp the variance; a naive
c     sum of squares would lose roughly seven digits here
c
      do i = 1, 10
         vbig(i) = 1.0d8 + dble(i)
      end do
      call avgstd (vbig,1,10,avg,std)
      call assert_real (avg,1.0d8+5.5d0,1.0d-6,
     &                  'avgstd large offset average')
      call assert_real (std,sd10,1.0d-9,
     &                  'avgstd large offset deviation')
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_thermint_schedule  --  window schedule  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_thermint_schedule" checks that the main lambda walks
c     the schedule table one window at a time and stops moving once
c     the final window has been passed
c
c
      subroutine test_thermint_schedule
      use thrmint
      implicit none
      integer k
      real*8 eps
      real*8 lref5(4)
      character*40 label
c
c
c     twenty one windows step lambda by one twentieth each time
c
      eps = 1.0d-12
      call resetti (21,10,100,50)
      do k = 1, 20
         call tischedule
         write (label,10)  k
   10    format ('tischedule 21 bin index ',i0)
         call assert_int (tibin,k+1,label)
         write (label,20)  k
   20    format ('tischedule 21 bin lambda ',i0)
         call assert_real (tilmda,1.0d0-dble(k)/20.0d0,eps,label)
      end do
c
c     the final window must sit exactly on the endpoint
c
      call assert_real (tilmda,0.0d0,0.0d0,
     &                  'tischedule 21 bin endpoint')
c
c     one call past the end advances the index but leaves lambda
c     where the last window left it
c
      call tischedule
      call assert_real (tilmda,0.0d0,0.0d0,
     &                  'tischedule 21 bin past end lambda')
      call assert_int (tibin,22,'tischedule 21 bin past end index')
c
c     five windows give lambda of 1.00, 0.75, 0.50, 0.25 and 0.00
c
      call resetti (5,10,40,20)
      lref5(1) = 0.75d0
      lref5(2) = 0.50d0
      lref5(3) = 0.25d0
      lref5(4) = 0.00d0
      do k = 1, 4
         call tischedule
         write (label,30)  k
   30    format ('tischedule 5 bin step ',i0)
         call assert_real (tilmda,lref5(k),eps,label)
      end do
      call assert_real (tilmda,0.0d0,0.0d0,
     &                  'tischedule 5 bin endpoint')
c
c     two windows sample only the two endpoints
c
      call resetti (2,10,40,20)
      call assert_real (tilmda,1.0d0,0.0d0,'tischedule 2 bin start')
      call tischedule
      call assert_real (tilmda,0.0d0,0.0d0,'tischedule 2 bin end')
      call assert_int (tibin,2,'tischedule 2 bin index')
c
c     an ascending schedule must not be clamped back toward zero
c
      call resetti (4,10,40,20)
      tilmdalist(1) = 0.0d0
      tilmdalist(2) = 0.1d0
      tilmdalist(3) = 0.4d0
      tilmdalist(4) = 1.0d0
      tilmda = tilmdalist(1)
      do k = 2, 4
         call tischedule
         write (label,40)  k
   40    format ('tischedule ascending step ',i0)
         call assert_real (tilmda,tilmdalist(k),eps,label)
      end do
      call tischedule
      call assert_real (tilmda,1.0d0,0.0d0,
     &                  'tischedule ascending holds at one')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_thermint_settisched  --  table setup  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_thermint_settisched" checks that the schedule table is
c     generated from TI-NBIN when no windows are given explicitly,
c     and is taken verbatim from the TI-WINDOW values when they are
c
c     the rejection of out of range and non-monotonic schedules is
c     not covered here, since those paths end in "fatal" and would
c     stop the test binary
c
c
      subroutine test_thermint_settisched
      use thrmint
      implicit none
      integer k
      real*8 eps
      character*40 label
c
c
c     with no explicit windows the table comes from TI-NBIN and must
c     reproduce the evenly spaced schedule exactly
c
      eps = 1.0d-12
      tinbin = 21
      call settisched (0,.false.)
      call assert_int (tinbin,21,'settisched 21 bin count')
      call assert_int (size(tilmdalist),21,'settisched 21 bin size')
      do k = 1, 21
         write (label,10)  k
   10    format ('settisched 21 bin value ',i0)
         call assert_real (tilmdalist(k),1.0d0-dble(k-1)/20.0d0,
     &                     eps,label)
      end do
      call assert_real (tilmdalist(21),0.0d0,0.0d0,
     &                  'settisched 21 bin endpoint')
      call assert_int (tibin,1,'settisched 21 bin start index')
      call assert_real (tilmda,1.0d0,0.0d0,'settisched 21 bin start')
c
c     five and two window schedules hit their endpoints exactly
c
      tinbin = 5
      call settisched (0,.false.)
      call assert_real (tilmdalist(1),1.00d0,eps,'settisched 5 bin 1')
      call assert_real (tilmdalist(2),0.75d0,eps,'settisched 5 bin 2')
      call assert_real (tilmdalist(3),0.50d0,eps,'settisched 5 bin 3')
      call assert_real (tilmdalist(4),0.25d0,eps,'settisched 5 bin 4')
      call assert_real (tilmdalist(5),0.00d0,0.0d0,
     &                  'settisched 5 bin 5')
      tinbin = 2
      call settisched (0,.false.)
      call assert_real (tilmdalist(1),1.0d0,0.0d0,'settisched 2 bin 1')
      call assert_real (tilmdalist(2),0.0d0,0.0d0,'settisched 2 bin 2')
c
c     an explicit descending schedule sets the window count itself
c     and is compacted down from the parse buffer
c
      if (allocated(tilmdalist))  deallocate (tilmdalist)
      allocate (tilmdalist(40))
      tilmdalist(1) = 1.0d0
      tilmdalist(2) = 0.9d0
      tilmdalist(3) = 0.2d0
      tilmdalist(4) = 0.0d0
      tinbin = 0
      call settisched (4,.false.)
      call assert_int (tinbin,4,'settisched explicit count')
      call assert_int (size(tilmdalist),4,'settisched explicit size')
      call assert_real (tilmdalist(1),1.0d0,eps,'settisched down 1')
      call assert_real (tilmdalist(2),0.9d0,eps,'settisched down 2')
      call assert_real (tilmdalist(3),0.2d0,eps,'settisched down 3')
      call assert_real (tilmdalist(4),0.0d0,eps,'settisched down 4')
      call assert_int (tibin,1,'settisched explicit start index')
      call assert_real (tilmda,1.0d0,eps,'settisched explicit start')
c
c     an ascending schedule is equally valid and starts at its own
c     first value rather than at one
c
      if (allocated(tilmdalist))  deallocate (tilmdalist)
      allocate (tilmdalist(40))
      tilmdalist(1) = 0.0d0
      tilmdalist(2) = 0.1d0
      tilmdalist(3) = 0.4d0
      tilmdalist(4) = 1.0d0
      tinbin = 0
      call settisched (4,.false.)
      call assert_int (tinbin,4,'settisched ascending count')
      call assert_real (tilmdalist(1),0.0d0,eps,'settisched up 1')
      call assert_real (tilmdalist(2),0.1d0,eps,'settisched up 2')
      call assert_real (tilmdalist(3),0.4d0,eps,'settisched up 3')
      call assert_real (tilmdalist(4),1.0d0,eps,'settisched up 4')
      call assert_real (tilmda,0.0d0,eps,'settisched ascending start')
c
c     the schedule need not touch either endpoint; any monotonic
c     run of values inside [0,1] is a valid set of windows
c
      if (allocated(tilmdalist))  deallocate (tilmdalist)
      allocate (tilmdalist(40))
      tilmdalist(1) = 0.75d0
      tilmdalist(2) = 0.70d0
      tilmdalist(3) = 0.20d0
      tinbin = 0
      call settisched (3,.false.)
      call assert_int (tinbin,3,'settisched interior count')
      call assert_real (tilmdalist(1),0.75d0,eps,
     &                  'settisched interior 1')
      call assert_real (tilmdalist(2),0.70d0,eps,
     &                  'settisched interior 2')
      call assert_real (tilmdalist(3),0.20d0,eps,
     &                  'settisched interior 3')
      call assert_real (tilmda,0.75d0,eps,
     &                  'settisched interior start')
c
c     a single window is legal and covers the whole trajectory,
c     which the old closed form schedule could not express
c
      if (allocated(tilmdalist))  deallocate (tilmdalist)
      allocate (tilmdalist(40))
      tilmdalist(1) = 0.5d0
      tinbin = 0
      call settisched (1,.false.)
      call assert_int (tinbin,1,'settisched single count')
      call assert_real (tilmdalist(1),0.5d0,eps,'settisched single')
      tinstepavg = 10
      tieqratio = 0.0d0
      call inittidyn (200)
      call assert_int (tiwindow,200,'settisched single window')
      call assert_int (tinequil,0,'settisched single equilibration')
      call assert_int (tinblock,20,'settisched single blocks')
      call assert_real (tilmda,0.5d0,eps,'settisched single lambda')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_thermint_data  --  accumulator setup  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_thermint_data" checks that the block average arrays are
c     sized to one row per lambda window and are fully cleared each
c     time the windows are laid out
c
c
      subroutine test_thermint_data
      use dlmda
      use thrmint
      implicit none
      integer i,j
      logical ok
c
c
c     seven windows of sixty steps hold two blocks of thirteen
c
      call resetti (7,13,60,30)
      call assert_int (size(tidedllist),13,'tidata block buffer size')
      call inittidyn (420)
      call assert_int (tiwindow,60,'tidata window length')
      call assert_int (tinequil,30,'tidata equilibration steps')
      call assert_int (tinblock,2,'tidata blocks per window')
      call assert_int (size(tilmdadedl,1),7,'tidata dedl window rows')
      call assert_int (size(tilmdadedl,2),2,'tidata dedl block cols')
      call assert_int (size(tilmdadedlstd,1),7,'tidata std window rows')
      call assert_int (size(tilmdadedlstd,2),2,'tidata std block cols')
      call assert_int (size(tinbcount),7,'tidata block count size')
      call assert_int (tibin,1,'tidata initial window index')
      call assert_real (tilmda,1.0d0,0.0d0,'tidata initial lambda')
c
c     every window starts with no recorded blocks
c
      ok = .true.
      do i = 1, 7
         if (tinbcount(i) .ne. 0)  ok = .false.
         do j = 1, 2
            if (tilmdadedl(i,j) .ne. 0.0d0)  ok = .false.
            if (tilmdadedlstd(i,j) .ne. 0.0d0)  ok = .false.
         end do
      end do
      call assert_logical (ok,.true.,'tidata rows start empty')
c
c     laying out the windows again clears whatever had accumulated
c
      tinbcount(3) = 1
      tilmdadedl(3,1) = 5.0d0
      tilmdadedlstd(3,1) = 6.0d0
      tibin = 4
      tilmda = 0.25d0
      call inittidyn (420)
      call assert_int (size(tilmdadedl,1),7,'tidata reinit rows')
      call assert_int (tinbcount(3),0,'tidata reinit block count')
      call assert_real (tilmdadedl(3,1),0.0d0,0.0d0,
     &                  'tidata reinit dedl')
      call assert_real (tilmdadedlstd(3,1),0.0d0,0.0d0,
     &                  'tidata reinit std')
      call assert_int (tibin,1,'tidata reinit window index')
      call assert_real (tilmda,1.0d0,0.0d0,'tidata reinit lambda')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_thermint_inittidyn  --  window layout  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_thermint_inittidyn" checks the window length and
c     equilibration count derived from the total number of steps,
c     including the integer truncation of both divisions
c
c
      subroutine test_thermint_inittidyn
      use thrmint
      implicit none
c
c
c     two hundred steps over five windows of forty
c
      call resetti (5,10,40,20)
      tinbin = 5
      tieqratio = 0.5d0
      call inittidyn (200)
      call assert_int (tiwindow,40,'inittidyn 200 step window')
      call assert_int (tinequil,20,'inittidyn 200 step equilibration')
      call assert_int (tinblock,2,'inittidyn 200 step blocks')
      call assert_int (tibin,1,'inittidyn 200 step window index')
      call assert_real (tilmda,1.0d0,0.0d0,'inittidyn 200 step lambda')
c
c     a quarter of each window discarded over twenty one windows
c
      tinbin = 21
      tieqratio = 0.25d0
      call inittidyn (2100)
      call assert_int (tiwindow,100,'inittidyn 2100 step window')
      call assert_int (tinequil,25,'inittidyn 2100 step equilibration')
      call assert_int (tinblock,7,'inittidyn 2100 step blocks')
c
c     integer truncation, 205/5 gives 41 steps and 41*0.5 gives 20
c
      tinbin = 5
      tieqratio = 0.5d0
      call inittidyn (205)
      call assert_int (tiwindow,41,'inittidyn 205 step window')
      call assert_int (tinequil,20,'inittidyn 205 step equilibration')
      call assert_int (tinblock,2,'inittidyn 205 step blocks')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_thermint_etidyn  --  block accumulation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_thermint_etidyn" checks that only production steps are
c     averaged into blocks, that each window gets its own blocks,
c     and that the schedule advances exactly at the window boundary
c
c
      subroutine test_thermint_etidyn
      use dlmda
      use thrmint
      implicit none
      integer b,w
      integer istep,tistep
      real*8 eps,sd10
      real*8 lref,dmax
      real*8 ref(2,5)
      real*8 lmdaseen(200)
      character*40 label
c
c
c     five windows of forty steps, twenty of them equilibration,
c     leaving two blocks of ten production steps per window
c
      eps = 1.0d-12
      sd10 = sqrt(8.25d0)
      ref(1,1) = 25.5d0
      ref(2,1) = 35.5d0
      ref(1,2) = 65.5d0
      ref(2,2) = 75.5d0
      ref(1,3) = 105.5d0
      ref(2,3) = 115.5d0
      ref(1,4) = 145.5d0
      ref(2,4) = 155.5d0
      ref(1,5) = 185.5d0
      ref(2,5) = 195.5d0
c
      call resetti (5,10,40,20)
      do istep = 1, 200
         dedl = dble(istep)
         lmdaseen(istep) = tilmda
         call etidyn (istep)
      end do
c
c     each window holds the mean of ten consecutive integers
c
      do w = 1, 5
         write (label,10)  w
   10    format ('etidyn window ',i0,' block count')
         call assert_int (tinbcount(w),2,label)
         do b = 1, 2
            write (label,20)  w,b
   20       format ('etidyn window ',i0,' block ',i0,' average')
            call assert_real (tilmdadedl(w,b),ref(b,w),eps,label)
            write (label,30)  w,b
   30       format ('etidyn window ',i0,' block ',i0,' deviation')
            call assert_real (tilmdadedlstd(w,b),sd10,eps,label)
         end do
      end do
c
c     the schedule must advance at the window boundary, not one
c     step early or one step late
c
      dmax = 0.0d0
      do istep = 1, 200
         lref = 1.0d0 - dble((istep-1)/40)/4.0d0
         dmax = max(dmax,abs(lmdaseen(istep)-lref))
      end do
      call assert_real (dmax,0.0d0,1.0d-12,
     &                  'etidyn lambda schedule over 200 steps')
      call assert_int (tibin,6,'etidyn final window index')
      call assert_real (tilmda,0.0d0,0.0d0,'etidyn final lambda')
c
c     the same run with every equilibration step poisoned; the
c     block averages must be untouched, which makes the discard
c     while equilibrating rule explicit
c
      call resetti (5,10,40,20)
      do istep = 1, 200
         tistep = mod(istep-1,tiwindow) + 1
         if (tistep .le. tinequil) then
            dedl = -1.0d9
         else
            dedl = dble(istep)
         end if
         call etidyn (istep)
      end do
      do w = 1, 5
         write (label,50)  w
   50    format ('etidyn poisoned window ',i0,' block count')
         call assert_int (tinbcount(w),2,label)
         do b = 1, 2
            write (label,60)  w,b
   60       format ('etidyn poisoned window ',i0,' block ',i0)
            call assert_real (tilmdadedl(w,b),ref(b,w),eps,label)
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_thermint_partialblock  --  block resets  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_thermint_partialblock" checks that samples stranded in
c     an incomplete block at the end of a window do not leak into
c     the block averages of the next window
c
c
      subroutine test_thermint_partialblock
      use dlmda
      use thrmint
      implicit none
      integer w
      integer istep
      real*8 eps
      character*40 label
c
c
c     twenty five production steps per window with blocks of ten
c     flushes two blocks and strands five samples
c
      eps = 1.0d-12
      call resetti (5,10,40,15)
      do istep = 1, 200
         dedl = dble(istep)
         call etidyn (istep)
      end do
      do w = 1, 5
         write (label,10)  w
   10    format ('partialblock window ',i0,' block count')
         call assert_int (tinbcount(w),2,label)
      end do
c
c     window one covers steps 16-25 and 26-35, and steps 36-40
c     are orphaned in the unflushed block
c
      call assert_real (tilmdadedl(1,1),20.5d0,eps,
     &                  'partialblock window 1 block 1')
      call assert_real (tilmdadedl(1,2),30.5d0,eps,
     &                  'partialblock window 1 block 2')
c
c     window two production starts at step 56, so a leak from the
c     previous window would drag this below 60.5
c
      call assert_real (tilmdadedl(2,1),60.5d0,eps,
     &                  'partialblock window 2 block 1')
      call assert_real (tilmdadedl(2,2),70.5d0,eps,
     &                  'partialblock window 2 block 2')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_thermint_trailing  --  extra step guard  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_thermint_trailing" checks that dynamics steps past the
c     end of the last lambda window are ignored rather than writing
c     beyond the block average arrays
c
c
      subroutine test_thermint_trailing
      use dlmda
      use thrmint
      implicit none
      integer w
      integer istep
      character*40 label
c
c
c     two hundred ten steps over five windows of forty, so the last
c     ten steps fall past the schedule
c
      call resetti (5,10,40,20)
      do istep = 1, 210
         dedl = dble(istep)
         call etidyn (istep)
      end do
      call assert_int (tibin,6,'trailing final window index')
      do w = 1, 5
         write (label,10)  w
   10    format ('trailing window ',i0,' block count')
         call assert_int (tinbcount(w),2,label)
      end do
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  subroutine resetti  --  reset TI test state  ##
c     ##                                               ##
c     ###################################################
c
c
c     "resetti" sizes the thermodynamic integration accumulators
c     directly and sets the scalar state to deterministic unit test
c     defaults, bypassing "inittidyn" so the accumulation tests do
c     not depend on it
c
c
      subroutine resetti (nbin,nstepavg,window,nequil)
      use dlmda
      use thrmint
      implicit none
      integer nbin,nstepavg,window,nequil
      integer i,j
c
c
c     put the sublambda maps on the power law branch, which is
c     plain arithmetic and needs no molecular system
c
      use_relstage = .false.
      elmdamap = 'EXP'
      plmdamap = 'EXP'
      vlmdamap = 'EXP'
      elmdaexp = 1
      plmdaexp = 1
      vlmdaexp = 1
c
c     set the window layout implied by the requested geometry
c
      use_ti = .true.
      tinbin = nbin
      tinstepavg = nstepavg
      tiwindow = window
      tinequil = nequil
      tieqratio = 0.0d0
      if (window .gt. 0)  tieqratio = dble(nequil) / dble(window)
      tinblock = 0
      if (nstepavg .gt. 0)  tinblock = (window-nequil) / nstepavg
      dedl = 0.0d0
c
c     clear any previous allocation and size the accumulators
c
      if (allocated(tinbcount))  deallocate (tinbcount)
      if (allocated(tidedllist))  deallocate (tidedllist)
      if (allocated(tilmdadedl))  deallocate (tilmdadedl)
      if (allocated(tilmdadedlstd))  deallocate (tilmdadedlstd)
      allocate (tinbcount(nbin))
      allocate (tidedllist(nstepavg))
      allocate (tilmdadedl(nbin,tinblock))
      allocate (tilmdadedlstd(nbin,tinblock))
      do i = 1, nstepavg
         tidedllist(i) = 0.0d0
      end do
c
c     build the default evenly spaced schedule, which also sets
c     "tibin" and "tilmda" to the first window
c
      call settisched (0,.false.)
      do i = 1, nbin
         tinbcount(i) = 0
         do j = 1, tinblock
            tilmdadedl(i,j) = 0.0d0
            tilmdadedlstd(i,j) = 0.0d0
         end do
      end do
      return
      end
c
c
c     ##################################################
c     ##                                              ##
c     ##  subroutine clearti  --  drop TI test state  ##
c     ##                                              ##
c     ##################################################
c
c
c     "clearti" turns thermodynamic integration back off and frees
c     the accumulators; the test binary is a single process, so a
c     case that left "use_ti" set would change how later tests set
c     up their lambda terms
c
c
      subroutine clearti
      use dlmda
      use thrmint
      implicit none
c
c
      use_ti = .false.
      tibin = 0
      tinblock = 0
      tilmda = 1.0d0
      if (allocated(tinbcount))  deallocate (tinbcount)
      if (allocated(tidedllist))  deallocate (tidedllist)
      if (allocated(tilmdadedl))  deallocate (tilmdadedl)
      if (allocated(tilmdadedlstd))  deallocate (tilmdadedlstd)
      if (allocated(tilmdalist))  deallocate (tilmdalist)
      return
      end
