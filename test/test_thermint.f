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
      call test_thermint_phase
      call test_thermint_setsched
      call test_thermint_fraction
      call test_thermint_data
      call test_thermint_inittidyn
      call test_thermint_uneven
      call test_thermint_etidyn
      call test_thermint_partialblock
      call test_thermint_trailing
      call test_thermint_save
      call test_thermint_stop
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
      use dlmda
      use mutant
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
         call nextlmdawin
         write (label,10)  k
   10    format ('nextlmdawin 21 bin index ',i0)
         call assert_int (lmdawin,k+1,label)
         write (label,20)  k
   20    format ('nextlmdawin 21 bin lambda ',i0)
         call assert_real (lambda,1.0d0-dble(k)/20.0d0,eps,label)
      end do
c
c     the final window must sit exactly on the endpoint
c
      call assert_real (lambda,0.0d0,0.0d0,
     &                  'nextlmdawin 21 bin endpoint')
c
c     one call past the end advances the index but leaves lambda
c     where the last window left it
c
      call nextlmdawin
      call assert_real (lambda,0.0d0,0.0d0,
     &                  'nextlmdawin 21 bin past end lambda')
      call assert_int (lmdawin,22,'nextlmdawin 21 bin past end index')
c
c     further calls change nothing once the schedule has run out
c
      call nextlmdawin
      call nextlmdawin
      call assert_real (lambda,0.0d0,0.0d0,
     &                  'nextlmdawin 21 bin repeat lambda')
      call assert_int (lmdawin,22,'nextlmdawin 21 bin repeat index')
c
c     five windows give lambda of 1.00, 0.75, 0.50, 0.25 and 0.00
c
      call resetti (5,10,40,20)
      lref5(1) = 0.75d0
      lref5(2) = 0.50d0
      lref5(3) = 0.25d0
      lref5(4) = 0.00d0
      do k = 1, 4
         call nextlmdawin
         write (label,30)  k
   30    format ('nextlmdawin 5 bin step ',i0)
         call assert_real (lambda,lref5(k),eps,label)
      end do
      call assert_real (lambda,0.0d0,0.0d0,
     &                  'nextlmdawin 5 bin endpoint')
c
c     two windows sample only the two endpoints
c
      call resetti (2,10,40,20)
      call assert_real (lambda,1.0d0,0.0d0,'nextlmdawin 2 bin start')
      call nextlmdawin
      call assert_real (lambda,0.0d0,0.0d0,'nextlmdawin 2 bin end')
      call assert_int (lmdawin,2,'nextlmdawin 2 bin index')
c
c     an ascending schedule must not be clamped back toward zero
c
      call resetti (4,10,40,20)
      lmdawinlist(1) = 0.0d0
      lmdawinlist(2) = 0.1d0
      lmdawinlist(3) = 0.4d0
      lmdawinlist(4) = 1.0d0
      lambda = lmdawinlist(1)
      do k = 2, 4
         call nextlmdawin
         write (label,40)  k
   40    format ('nextlmdawin ascending step ',i0)
         call assert_real (lambda,lmdawinlist(k),eps,label)
      end do
      call nextlmdawin
      call assert_real (lambda,1.0d0,0.0d0,
     &                  'nextlmdawin ascending holds at one')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_thermint_phase  --  window step phase  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_thermint_phase" checks the position that "lmdawinphase"
c     reports for a dynamics step inside the current lambda window,
c     that "lmdawinstep" advances only on the final step of a window,
c     and that all of them do nothing when no schedule is laid out
c
c
      subroutine test_thermint_phase
      use dlmda
      use mutant
      use thrmint
      implicit none
      integer iprod
c
c
c     three windows of forty steps, the first ten equilibration
c
      call resetti (3,10,40,10)
      call assert_int (lmdawinlen,40,'phase window length')
      call assert_int (lmdawineq,10,'phase equilibration steps')
c
c     equilibration steps are not production, and the count starts
c     at one on the first step after equilibration
c
      call lmdawinphase (1,iprod)
      call assert_int (iprod,0,'phase first step')
      call lmdawinphase (10,iprod)
      call assert_int (iprod,0,'phase last equilibration step')
      call lmdawinphase (11,iprod)
      call assert_int (iprod,1,'phase first production step')
      call lmdawinphase (40,iprod)
      call assert_int (iprod,30,'phase boundary step')
c
c     looking at a step does not move the schedule
c
      call assert_int (lmdawin,1,'phase inspection keeps window')
      call assert_real (lambda,1.0d0,0.0d0,
     &                  'phase inspection keeps lambda')
c
c     the window advances on its final step and on no other
c
      call lmdawinstep (39)
      call assert_int (lmdawin,1,'phase no advance before boundary')
      call lmdawinstep (40)
      call assert_int (lmdawin,2,'phase advance at boundary')
      call assert_real (lambda,0.5d0,0.0d0,'phase second lambda')
      call lmdawinstep (40)
      call assert_int (lmdawin,2,'phase no second advance')
c
c     positions in a later window are measured from its own start
c
      call lmdawinphase (50,iprod)
      call assert_int (iprod,0,'phase second window equilibration')
      call lmdawinphase (51,iprod)
      call assert_int (iprod,1,'phase second window production')
      call lmdawinphase (80,iprod)
      call assert_int (iprod,30,'phase second window boundary')
c
c     once the schedule has run out no step is in production
c
      call lmdawinstep (80)
      call lmdawinstep (120)
      call assert_int (lmdawin,4,'phase schedule exhausted')
      call lmdawinphase (121,iprod)
      call assert_int (iprod,0,'phase exhausted step')
      call lmdawinstep (121)
      call assert_int (lmdawin,4,'phase exhausted no advance')
      call assert_real (lambda,0.0d0,0.0d0,'phase exhausted lambda')
c
c     a schedule that was built but never laid out over a run, as
c     in a program without lambda windows, is left alone
c
      call resetti (3,10,40,10)
      deallocate (lmdawinend)
      tinbcount = 0
      dedl = 5.0d0
      call lmdawinphase (11,iprod)
      call assert_int (iprod,0,'phase no layout step')
      call etidyn (11)
      call assert_int (tinbcount,0,'phase no layout sample')
      call lmdawinstep (40)
      call assert_int (lmdawin,1,'phase no layout advance')
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_thermint_setsched  --  table setup  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_thermint_setsched" checks that the schedule table is
c     generated from LAMBDA-NWINDOW when no windows are given, and
c     is taken verbatim from the LAMBDA-WINDOW values when they are
c
c     the rejection of out of range and non-monotonic schedules is
c     not covered here, since those paths end in "fatal" and would
c     stop the test binary
c
c
      subroutine test_thermint_setsched
      use dlmda
      use mutant
      use thrmint
      implicit none
      integer k
      integer nblock
      real*8 eps
      character*40 label
c
c
c     with no explicit windows the table comes from LAMBDA-NWINDOW
c     and must reproduce the evenly spaced schedule exactly
c
      eps = 1.0d-12
      nlmdawin = 21
      call setlmdasched (0,.false.)
      call assert_int (nlmdawin,21,'setlmdasched 21 bin count')
      call assert_int (size(lmdawinlist),21,'setlmdasched 21 bin size')
      do k = 1, 21
         write (label,10)  k
   10    format ('setlmdasched 21 bin value ',i0)
         call assert_real (lmdawinlist(k),1.0d0-dble(k-1)/20.0d0,
     &                     eps,label)
      end do
      call assert_real (lmdawinlist(21),0.0d0,0.0d0,
     &                  'setlmdasched 21 bin endpoint')
      call assert_int (lmdawin,1,'setlmdasched 21 bin start index')
      call assert_real (lambda,1.0d0,0.0d0,'setlmdasched 21 bin start')
c
c     five and two window schedules hit their endpoints exactly
c
      nlmdawin = 5
      call setlmdasched (0,.false.)
      call assert_real (lmdawinlist(1),1.00d0,eps,
     &                  'setlmdasched 5 bin 1')
      call assert_real (lmdawinlist(2),0.75d0,eps,
     &                  'setlmdasched 5 bin 2')
      call assert_real (lmdawinlist(3),0.50d0,eps,
     &                  'setlmdasched 5 bin 3')
      call assert_real (lmdawinlist(4),0.25d0,eps,
     &                  'setlmdasched 5 bin 4')
      call assert_real (lmdawinlist(5),0.00d0,0.0d0,
     &                  'setlmdasched 5 bin 5')
      nlmdawin = 2
      call setlmdasched (0,.false.)
      call assert_real (lmdawinlist(1),1.0d0,0.0d0,
     &                  'setlmdasched 2 bin 1')
      call assert_real (lmdawinlist(2),0.0d0,0.0d0,
     &                  'setlmdasched 2 bin 2')
c
c     an explicit descending schedule sets the window count itself
c     and is compacted down from the parse buffer
c
      call tibufinit
      lmdawinlist(1) = 1.0d0
      lmdawinlist(2) = 0.9d0
      lmdawinlist(3) = 0.2d0
      lmdawinlist(4) = 0.0d0
      nlmdawin = 0
      call setlmdasched (4,.false.)
      call assert_int (nlmdawin,4,'setlmdasched explicit count')
      call assert_int (size(lmdawinlist),4,'setlmdasched explicit size')
      call assert_int (size(lmdawinfrac),4,
     &                 'setlmdasched explicit fracs')
      call assert_real (lmdawinlist(1),1.0d0,eps,'setlmdasched down 1')
      call assert_real (lmdawinlist(2),0.9d0,eps,'setlmdasched down 2')
      call assert_real (lmdawinlist(3),0.2d0,eps,'setlmdasched down 3')
      call assert_real (lmdawinlist(4),0.0d0,eps,'setlmdasched down 4')
      call assert_int (lmdawin,1,'setlmdasched explicit start index')
      call assert_real (lambda,1.0d0,eps,'setlmdasched explicit start')
c
c     with no fractions asked for, the run is split evenly
c
      do k = 1, 4
         write (label,20)  k
   20    format ('setlmdasched even share ',i0)
         call assert_real (lmdawinfrac(k),0.25d0,eps,label)
      end do
c
c     an ascending schedule is equally valid and starts at its own
c     first value rather than at one
c
      call tibufinit
      lmdawinlist(1) = 0.0d0
      lmdawinlist(2) = 0.1d0
      lmdawinlist(3) = 0.4d0
      lmdawinlist(4) = 1.0d0
      nlmdawin = 0
      call setlmdasched (4,.false.)
      call assert_int (nlmdawin,4,'setlmdasched ascending count')
      call assert_real (lmdawinlist(1),0.0d0,eps,'setlmdasched up 1')
      call assert_real (lmdawinlist(2),0.1d0,eps,'setlmdasched up 2')
      call assert_real (lmdawinlist(3),0.4d0,eps,'setlmdasched up 3')
      call assert_real (lmdawinlist(4),1.0d0,eps,'setlmdasched up 4')
      call assert_real (lambda,0.0d0,eps,'setlmdasched ascending start')
c
c     the schedule need not touch either endpoint; any monotonic
c     run of values inside [0,1] is a valid set of windows
c
      call tibufinit
      lmdawinlist(1) = 0.75d0
      lmdawinlist(2) = 0.70d0
      lmdawinlist(3) = 0.20d0
      nlmdawin = 0
      call setlmdasched (3,.false.)
      call assert_int (nlmdawin,3,'setlmdasched interior count')
      call assert_real (lmdawinlist(1),0.75d0,eps,
     &                  'setlmdasched interior 1')
      call assert_real (lmdawinlist(2),0.70d0,eps,
     &                  'setlmdasched interior 2')
      call assert_real (lmdawinlist(3),0.20d0,eps,
     &                  'setlmdasched interior 3')
      call assert_real (lambda,0.75d0,eps,
     &                  'setlmdasched interior start')
c
c     a single window is legal and covers the whole trajectory,
c     which the old closed form schedule could not express
c
      call tibufinit
      lmdawinlist(1) = 0.5d0
      nlmdawin = 0
      call setlmdasched (1,.false.)
      call assert_int (nlmdawin,1,'setlmdasched single count')
      call assert_real (lmdawinlist(1),0.5d0,eps,'setlmdasched single')
      call assert_real (lmdawinfrac(1),1.0d0,eps,
     &                  'setlmdasched single share')
      call tisetavg (10)
      lmdawinratio = 0.0d0
      call inittidyn (200)
      call assert_int (lmdawinlen,200,'setlmdasched single window')
      call assert_int (lmdawineq,0,'setlmdasched single equilibration')
      nblock = (lmdawinlen-lmdawineq) / tinstepavg
      call assert_int (nblock,20,'setlmdasched single blocks')
      call assert_int (tinbtot,20,'setlmdasched single capacity')
      call assert_real (lambda,0.5d0,eps,'setlmdasched single lambda')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine tibufinit  --  fresh TI keyword parse buffer  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "tibufinit" sizes the schedule parse buffers the way "mutate"
c     does before reading keywords, with every time fraction left
c     marked as not specified
c
c
      subroutine tibufinit
      use dlmda
      implicit none
      integer i
c
c
      if (allocated(lmdawinlist))  deallocate (lmdawinlist)
      if (allocated(lmdawinfrac))  deallocate (lmdawinfrac)
      allocate (lmdawinlist(40))
      allocate (lmdawinfrac(40))
      do i = 1, 40
         lmdawinlist(i) = 0.0d0
         lmdawinfrac(i) = -1.0d0
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine tisetavg  --  resize the TI sample buffer  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "tisetavg" changes the number of steps averaged into a block
c     and resizes the sample buffer to match; "mutate" sizes that
c     buffer once from the keyword value, so a test that moves
c     "tinstepavg" afterwards has to resize it as well
c
c
      subroutine tisetavg (nstepavg)
      use thrmint
      implicit none
      integer nstepavg
      integer i
c
c
      tinstepavg = nstepavg
      if (allocated(tidedllist))  deallocate (tidedllist)
      allocate (tidedllist(nstepavg))
      do i = 1, nstepavg
         tidedllist(i) = 0.0d0
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_thermint_fraction  --  window time share  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_thermint_fraction" checks how the time fraction given on
c     a "LAMBDA-WINDOW" line is resolved: windows asking for a share
c     keep it, windows that stay silent split whatever is left, and
c     the whole table is rescaled to span exactly one run
c
c
      subroutine test_thermint_fraction
      use dlmda
      implicit none
      integer k
      real*8 eps,fsum
c
c
c     two of four windows name a share, the other two split the rest
c
      eps = 1.0d-12
      call tibufinit
      lmdawinlist(1) = 1.0d0
      lmdawinlist(2) = 0.7d0
      lmdawinlist(3) = 0.3d0
      lmdawinlist(4) = 0.0d0
      lmdawinfrac(1) = 0.4d0
      lmdawinfrac(3) = 0.2d0
      nlmdawin = 0
      call setlmdasched (4,.false.)
      call assert_real (lmdawinfrac(1),0.4d0,eps,'tifrac given share 1')
      call assert_real (lmdawinfrac(2),0.2d0,eps,
     &                  'tifrac spread share 2')
      call assert_real (lmdawinfrac(3),0.2d0,eps,'tifrac given share 3')
      call assert_real (lmdawinfrac(4),0.2d0,eps,
     &                  'tifrac spread share 4')
c
c     the resolved shares always cover the whole run
c
      fsum = 0.0d0
      do k = 1, 4
         fsum = fsum + lmdawinfrac(k)
      end do
      call assert_real (fsum,1.0d0,eps,'tifrac shares total one')
c
c     a single unspecified window absorbs everything left over
c
      call tibufinit
      lmdawinlist(1) = 1.0d0
      lmdawinlist(2) = 0.5d0
      lmdawinlist(3) = 0.0d0
      lmdawinfrac(1) = 0.25d0
      lmdawinfrac(2) = 0.25d0
      nlmdawin = 0
      call setlmdasched (3,.false.)
      call assert_real (lmdawinfrac(3),0.5d0,eps,
     &                  'tifrac single leftover')
c
c     shares that do not total one are rescaled, so the same ratios
c     describe the same schedule however they were written down
c
      call tibufinit
      lmdawinlist(1) = 1.0d0
      lmdawinlist(2) = 0.5d0
      lmdawinlist(3) = 0.0d0
      lmdawinfrac(1) = 0.1d0
      lmdawinfrac(2) = 0.2d0
      lmdawinfrac(3) = 0.1d0
      nlmdawin = 0
      call setlmdasched (3,.false.)
      call assert_real (lmdawinfrac(1),0.25d0,eps,'tifrac rescaled 1')
      call assert_real (lmdawinfrac(2),0.50d0,eps,'tifrac rescaled 2')
      call assert_real (lmdawinfrac(3),0.25d0,eps,'tifrac rescaled 3')
c
c     the same ratios written to total one give the same schedule
c
      call tibufinit
      lmdawinlist(1) = 1.0d0
      lmdawinlist(2) = 0.5d0
      lmdawinlist(3) = 0.0d0
      lmdawinfrac(1) = 0.25d0
      lmdawinfrac(2) = 0.50d0
      lmdawinfrac(3) = 0.25d0
      nlmdawin = 0
      call setlmdasched (3,.false.)
      call assert_real (lmdawinfrac(1),0.25d0,eps,'tifrac direct 1')
      call assert_real (lmdawinfrac(2),0.50d0,eps,'tifrac direct 2')
      call assert_real (lmdawinfrac(3),0.25d0,eps,'tifrac direct 3')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_thermint_uneven  --  uneven schedule  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_thermint_uneven" runs a schedule whose windows have very
c     different lengths, and checks that the step boundaries, the
c     equilibration split and the recorded blocks all follow the
c     requested time fractions
c
c
      subroutine test_thermint_uneven
      use dlmda
      use thrmint
      implicit none
      integer istep
      integer nblock
      real*8 eps
c
c
c     half the run at lambda one, most of the rest at lambda a half,
c     and a brief two percent visit to lambda zero
c
      eps = 1.0d-12
      call tibufinit
      lmdawinlist(1) = 1.0d0
      lmdawinlist(2) = 0.5d0
      lmdawinlist(3) = 0.0d0
      lmdawinfrac(1) = 0.5d0
      lmdawinfrac(3) = 0.02d0
      nlmdawin = 0
      call setlmdasched (3,.false.)
      call assert_real (lmdawinfrac(2),0.48d0,eps,'uneven middle share')
c
c     one thousand steps split 500, 480 and 20 with half of each
c     window spent equilibrating
c
      call tisetavg (10)
      lmdawinratio = 0.5d0
      call inittidyn (1000)
      call assert_int (lmdawinend(1),500,'uneven first boundary')
      call assert_int (lmdawinend(2),980,'uneven second boundary')
      call assert_int (lmdawinend(3),1000,'uneven last boundary')
c
c     the first window keeps 250 production steps, giving 25 blocks,
c     the second 240 steps giving 24, and the short window keeps 10
c
      call assert_int (lmdawinlen,500,'uneven first window length')
      call assert_int (lmdawineq,250,'uneven first equilibration')
      nblock = (lmdawinlen-lmdawineq) / tinstepavg
      call assert_int (nblock,25,'uneven first blocks')
      call assert_int (tinbtot,50,'uneven total capacity')
c
c     walking the whole run fills the accumulators exactly
c
      do istep = 1, 1000
         dedl = dble(istep)
         call tidynstep (istep)
      end do
      call assert_int (lmdawin,4,'uneven schedule exhausted')
      call assert_int (tinbcount,50,'uneven blocks recorded')
c
c     the lambda tags show where the time actually went
c
      call assert_real (tilmdahist(1),1.0d0,eps,'uneven first lambda')
      call assert_real (tilmdahist(25),1.0d0,eps,
     &                  'uneven lambda one end')
      call assert_real (tilmdahist(26),0.5d0,eps,'uneven middle lambda')
      call assert_real (tilmdahist(49),0.5d0,eps,'uneven middle end')
      call assert_real (tilmdahist(50),0.0d0,eps,'uneven short lambda')
c
c     the first block averages steps 251 to 260 of the run
c
      call assert_real (tilmdadedl(1),255.5d0,eps,'uneven first block')
c
c     the short window equilibrates for its first ten steps and
c     records a single block from steps 991 to 1000
c
      call assert_real (tilmdadedl(50),995.5d0,eps,'uneven short block')
c
c     a window too short for one block records nothing at all, but
c     the run still visits its lambda
c
      call tibufinit
      lmdawinlist(1) = 1.0d0
      lmdawinlist(2) = 0.5d0
      lmdawinlist(3) = 0.0d0
      lmdawinfrac(2) = 0.01d0
      nlmdawin = 0
      call setlmdasched (3,.false.)
      call tisetavg (50)
      lmdawinratio = 0.5d0
      call inittidyn (1000)
      call assert_int (lmdawinend(2)-lmdawinend(1),10,
     &                 'uneven tiny window')
      call assert_int (tinbtot,8,'uneven tiny window capacity')
      do istep = 1, 1000
         dedl = dble(istep)
         call tidynstep (istep)
      end do
      call assert_int (tinbcount,8,'uneven tiny window records')
      call assert_real (tilmdahist(4),1.0d0,eps,'uneven tiny before')
      call assert_real (tilmdahist(5),0.0d0,eps,'uneven tiny after')
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
c     "test_thermint_data" checks that the recording arrays are sized
c     to the exact number of blocks the schedule can hold and are
c     fully cleared each time the windows are laid out
c
c
      subroutine test_thermint_data
      use dlmda
      use mutant
      use thrmint
      implicit none
      integer i
      integer nblock
      logical ok
c
c
c     seven windows of sixty steps hold two blocks of thirteen
c
      call resetti (7,13,60,30)
      call assert_int (size(tidedllist),13,'tidata block buffer size')
      call inittidyn (420)
      call assert_int (lmdawinlen,60,'tidata window length')
      call assert_int (lmdawineq,30,'tidata equilibration steps')
      nblock = (lmdawinlen-lmdawineq) / tinstepavg
      call assert_int (nblock,2,'tidata blocks in first window')
      call assert_int (tinbtot,14,'tidata total block capacity')
      call assert_int (size(tilmdadedl),14,'tidata dedl length')
      call assert_int (size(tilmdadedlstd),14,'tidata std length')
      call assert_int (size(tilmdahist),14,'tidata lambda hist length')
      call assert_int (size(lmdawinend),7,'tidata window end size')
      call assert_int (tinbcount,0,'tidata initial block count')
      call assert_int (tinbsave,0,'tidata initial blocks saved')
      call assert_int (lmdawin,1,'tidata initial window index')
      call assert_real (lambda,1.0d0,0.0d0,'tidata initial lambda')
c
c     the boundaries march evenly to the end of the trajectory
c
      ok = .true.
      do i = 1, 7
         if (lmdawinend(i) .ne. 60*i)  ok = .false.
      end do
      call assert_logical (ok,.true.,'tidata window boundaries')
      call assert_int (lmdawinend(7),420,'tidata last boundary')
c
c     nothing is recorded before the run starts
c
      ok = .true.
      do i = 1, 14
         if (tilmdahist(i) .ne. 0.0d0)  ok = .false.
         if (tilmdadedl(i) .ne. 0.0d0)  ok = .false.
         if (tilmdadedlstd(i) .ne. 0.0d0)  ok = .false.
      end do
      call assert_logical (ok,.true.,'tidata rows start empty')
c
c     laying out the windows again clears whatever had accumulated
c
      tinbcount = 3
      tinbsave = 2
      tilmdahist(3) = 4.0d0
      tilmdadedl(3) = 5.0d0
      tilmdadedlstd(3) = 6.0d0
      lmdawin = 4
      lambda = 0.25d0
      call inittidyn (420)
      call assert_int (size(tilmdadedl),14,'tidata reinit length')
      call assert_int (tinbcount,0,'tidata reinit block count')
      call assert_int (tinbsave,0,'tidata reinit blocks saved')
      call assert_real (tilmdahist(3),0.0d0,0.0d0,
     &                  'tidata reinit lambda hist')
      call assert_real (tilmdadedl(3),0.0d0,0.0d0,
     &                  'tidata reinit dedl')
      call assert_real (tilmdadedlstd(3),0.0d0,0.0d0,
     &                  'tidata reinit std')
      call assert_int (lmdawin,1,'tidata reinit window index')
      call assert_real (lambda,1.0d0,0.0d0,'tidata reinit lambda')
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
      use dlmda
      use mutant
      use thrmint
      implicit none
      integer nblock
c
c
c     two hundred steps over five windows of forty
c
      call resetti (5,10,40,20)
      call inittidyn (200)
      call assert_int (lmdawinlen,40,'inittidyn 200 step window')
      call assert_int (lmdawineq,20,'inittidyn 200 step equilibration')
      nblock = (lmdawinlen-lmdawineq) / tinstepavg
      call assert_int (nblock,2,'inittidyn 200 step blocks')
      call assert_int (tinbtot,10,'inittidyn 200 step capacity')
      call assert_int (lmdawinend(5),200,'inittidyn 200 step coverage')
      call assert_int (lmdawin,1,'inittidyn 200 step window index')
      call assert_real (lambda,1.0d0,0.0d0,'inittidyn 200 step lambda')
c
c     a quarter of each window discarded over twenty one windows
c
      call resetti (21,10,100,25)
      call inittidyn (2100)
      call assert_int (lmdawinlen,100,'inittidyn 2100 step window')
      call assert_int (lmdawineq,25,'inittidyn 2100 step equilibration')
      nblock = (lmdawinlen-lmdawineq) / tinstepavg
      call assert_int (nblock,7,'inittidyn 2100 step blocks')
      call assert_int (tinbtot,147,'inittidyn 2100 step capacity')
      call assert_int (lmdawinend(21),2100,
     &                 'inittidyn 2100 step coverage')
c
c     an odd step count is absorbed by the boundaries rather than
c     left as a trailing remainder, and 41*0.5 truncates to 20
c
      call resetti (5,10,40,20)
      call inittidyn (205)
      call assert_int (lmdawinlen,41,'inittidyn 205 step window')
      call assert_int (lmdawineq,20,'inittidyn 205 step equilibration')
      nblock = (lmdawinlen-lmdawineq) / tinstepavg
      call assert_int (nblock,2,'inittidyn 205 step blocks')
      call assert_int (lmdawinend(1),41,'inittidyn 205 step boundary')
      call assert_int (lmdawinend(5),205,'inittidyn 205 step coverage')
      call assert_int (tinbtot,10,'inittidyn 205 step capacity')
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
c     averaged into blocks, that the blocks are appended in the order
c     they were recorded along with the lambda that produced them,
c     and that the schedule advances exactly at the window boundary
c
c
      subroutine test_thermint_etidyn
      use dlmda
      use mutant
      use thrmint
      implicit none
      integer b,i,w
      integer istep,tistep
      real*8 eps,sd10
      real*8 lref,dmax
      real*8 ref(10)
      real*8 lmdaseen(200)
      character*40 label
c
c
c     five windows of forty steps, twenty of them equilibration,
c     leaving two blocks of ten production steps per window
c
      eps = 1.0d-12
      sd10 = sqrt(8.25d0)
      ref(1) = 25.5d0
      ref(2) = 35.5d0
      ref(3) = 65.5d0
      ref(4) = 75.5d0
      ref(5) = 105.5d0
      ref(6) = 115.5d0
      ref(7) = 145.5d0
      ref(8) = 155.5d0
      ref(9) = 185.5d0
      ref(10) = 195.5d0
c
      call resetti (5,10,40,20)
      do istep = 1, 200
         dedl = dble(istep)
         lmdaseen(istep) = lambda
         call tidynstep (istep)
      end do
c
c     the blocks land end to end, each holding the mean of ten
c     consecutive integers and tagged with its own lambda
c
      call assert_int (tinbcount,10,'etidyn total block count')
      do i = 1, 10
         w = (i-1)/2 + 1
         b = mod(i-1,2) + 1
         write (label,20)  w,b
   20    format ('etidyn window ',i0,' block ',i0,' average')
         call assert_real (tilmdadedl(i),ref(i),eps,label)
         write (label,30)  w,b
   30    format ('etidyn window ',i0,' block ',i0,' deviation')
         call assert_real (tilmdadedlstd(i),sd10,eps,label)
         write (label,40)  w,b
   40    format ('etidyn window ',i0,' block ',i0,' lambda')
         lref = 1.0d0 - dble(w-1)/4.0d0
         call assert_real (tilmdahist(i),lref,eps,label)
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
      call assert_int (lmdawin,6,'etidyn final window index')
      call assert_real (lambda,0.0d0,0.0d0,'etidyn final lambda')
c
c     the same run with every equilibration step poisoned; the
c     block averages must be untouched, which makes the discard
c     while equilibrating rule explicit
c
      call resetti (5,10,40,20)
      do istep = 1, 200
         tistep = mod(istep-1,lmdawinlen) + 1
         if (tistep .le. lmdawineq) then
            dedl = -1.0d9
         else
            dedl = dble(istep)
         end if
         call tidynstep (istep)
      end do
      call assert_int (tinbcount,10,'etidyn poisoned block count')
      do i = 1, 10
         write (label,60)  i
   60    format ('etidyn poisoned block ',i0)
         call assert_real (tilmdadedl(i),ref(i),eps,label)
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
      integer istep
      real*8 eps
c
c
c     twenty five production steps per window with blocks of ten
c     flushes two blocks and strands five samples
c
      eps = 1.0d-12
      call resetti (5,10,40,15)
      do istep = 1, 200
         dedl = dble(istep)
         call tidynstep (istep)
      end do
      call assert_int (tinbtot,10,'partialblock capacity')
      call assert_int (tinbcount,10,'partialblock block count')
c
c     window one covers steps 16-25 and 26-35, and steps 36-40
c     are orphaned in the unflushed block
c
      call assert_real (tilmdadedl(1),20.5d0,eps,
     &                  'partialblock window 1 block 1')
      call assert_real (tilmdadedl(2),30.5d0,eps,
     &                  'partialblock window 1 block 2')
c
c     window two production starts at step 56, so a leak from the
c     previous window would drag this below 60.5
c
      call assert_real (tilmdadedl(3),60.5d0,eps,
     &                  'partialblock window 2 block 1')
      call assert_real (tilmdadedl(4),70.5d0,eps,
     &                  'partialblock window 2 block 2')
c
c     the stranded samples never earn a lambda tag of their own
c
      call assert_real (tilmdahist(3),0.75d0,eps,
     &                  'partialblock window 2 lambda')
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
c     "test_thermint_trailing" checks that the windows laid out by
c     "inittidyn" span the whole trajectory, and that steps pushed
c     past the end of the schedule are ignored rather than writing
c     beyond the block average arrays
c
c
      subroutine test_thermint_trailing
      use dlmda
      use thrmint
      implicit none
      integer istep
c
c
c     an awkward step count is still covered to the last step, so
c     no dynamics is left running outside the schedule
c
      call resetti (5,10,40,20)
      call inittidyn (203)
      call assert_int (lmdawinend(5),203,'trailing schedule coverage')
      do istep = 1, 203
         dedl = dble(istep)
         call tidynstep (istep)
      end do
      call assert_int (lmdawin,6,'trailing schedule exhausted')
      call assert_int (tinbcount,tinbtot,'trailing capacity filled')
c
c     two hundred ten steps over five windows of forty, so the last
c     ten steps fall past the schedule and must be dropped
c
      call resetti (5,10,40,20)
      do istep = 1, 210
         dedl = dble(istep)
         call tidynstep (istep)
      end do
      call assert_int (lmdawin,6,'trailing final window index')
      call assert_int (tinbcount,10,'trailing block count')
      call assert_int (tinbcount,tinbtot,'trailing no overflow')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_thermint_save  --  incremental output  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_thermint_save" checks that the block averages reach the
c     external file as the windows complete, that each row is
c     written exactly once, and that a call with nothing new to
c     report leaves the file alone
c
c
      subroutine test_thermint_save
      use dlmda
      use files
      use thrmint
      implicit none
      integer istep
      integer nblock
      integer nrow
      integer oldleng
      real*8 lam
      logical exist
      character*240 oldname
c
c
c     write into the test root under a name of our own, keeping the
c     caller's file name to restore on the way out
c
      oldname = filename
      oldleng = leng
      filename = 'tisave_tmp'
      leng = 10
      call tiwipe ('tisave_tmp.ti')
c
c     five windows of forty steps, twenty of them equilibration,
c     leaving two blocks of ten production steps per window
c
      call resetti (5,10,40,20)
      call inittidyn (200)
      nblock = (lmdawinlen-lmdawineq) / tinstepavg
      call assert_int (nblock,2,'saveti blocks per window')
c
c     "inittidyn" must not touch the file system; only "prttihead"
c     creates the file, and it starts with a header alone
c
      inquire (file='tisave_tmp.ti',exist=exist)
      call assert_logical (exist,.false.,
     &                     'saveti inittidyn writes no file')
      call prttihead
      call ticount ('tisave_tmp.ti',nrow)
      call assert_int (nrow,0,'saveti header only at start')
c
c     the first two windows contribute two blocks each
c
      do istep = 1, 80
         dedl = dble(istep)
         call tidynstep (istep)
      end do
      call saveti
      call ticount ('tisave_tmp.ti',nrow)
      call assert_int (nrow,4,'saveti rows after two windows')
c
c     two more windows append four more rows and no others
c
      do istep = 81, 160
         dedl = dble(istep)
         call tidynstep (istep)
      end do
      call saveti
      call ticount ('tisave_tmp.ti',nrow)
      call assert_int (nrow,8,'saveti rows after four windows')
c
c     a call with nothing new must not repeat anything
c
      call saveti
      call ticount ('tisave_tmp.ti',nrow)
      call assert_int (nrow,8,'saveti no rows when nothing is new')
c
c     the final window closes out the schedule
c
      do istep = 161, 200
         dedl = dble(istep)
         call tidynstep (istep)
      end do
      call saveti
      call ticount ('tisave_tmp.ti',nrow)
      call assert_int (nrow,10,'saveti rows after all windows')
c
c     every recorded block must appear exactly once, so the row
c     count has to match the number of blocks collected
c
      call assert_int (nrow,tinbcount,
     &                 'saveti rows match blocks recorded')
      call assert_int (tinbsave,tinbcount,'saveti fully written')
c
c     the rows carry the lambda that produced them, in the order
c     the blocks were recorded
c
      call tirowlmda ('tisave_tmp.ti',1,lam)
      call assert_real (lam,1.0d0,1.0d-8,'saveti first row lambda')
      call tirowlmda ('tisave_tmp.ti',10,lam)
      call assert_real (lam,0.0d0,1.0d-8,'saveti last row lambda')
      call tirowlmda ('tisave_tmp.ti',5,lam)
      call assert_real (lam,0.5d0,1.0d-8,'saveti middle row lambda')
c
c     clean up and put the caller's file name back
c
      call tiwipe ('tisave_tmp.ti')
      filename = oldname
      leng = oldleng
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_thermint_stop  --  requested stop rows  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_thermint_stop" runs dynamics that is asked to stop at
c     its first trajectory save, and checks that a block average
c     completed on that final step still reaches the output file,
c     while a block that is not yet complete writes nothing
c
c
      subroutine test_thermint_stop
      implicit none
      integer ist
      integer nrow
      logical exist
c
c
c     one window with no equilibration and blocks of ten steps; the
c     stop comes at step ten, which also completes the first block
c
      call tistopprep ('thermint_stop','10')
      call pushdir ('file/thermint_stop')
      call run_prog ('dynamic','water2 40 0.1 0.001 2 298',
     &               'out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'tistop boundary status')
         inquire (file='water2.end',exist=exist)
         call assert_logical (exist,.false.,'tistop stop file used')
         call ticount ('water2.ti',nrow)
         call assert_int (nrow,1,'tistop block at stop step')
      end if
      call popdir
      call tnist_clean ('thermint_stop')
c
c     with blocks of fifteen steps the stop falls inside the first
c     block, which is discarded and not written as a partial row
c
      call tistopprep ('thermint_stop','15')
      call pushdir ('file/thermint_stop')
      call run_prog ('dynamic','water2 40 0.1 0.001 2 298',
     &               'out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'tistop partial status')
         call ticount ('water2.ti',nrow)
         call assert_int (nrow,0,'tistop partial block dropped')
      end if
      call popdir
      call tnist_clean ('thermint_stop')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine tistopprep  --  fixture for a stopped TI run  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "tistopprep" creates a scratch directory holding a small water
c     box set up for a single window thermodynamic integration with
c     the requested block size, along with the file that asks the
c     dynamics to stop at its first trajectory save
c
c
      subroutine tistopprep (work,nstepavg)
      implicit none
      integer iend
      integer freeunit
      character*(*) work,nstepavg
      character*512 cmd
c
c
c     copy the structure and its keyfile into the scratch directory
c
      call pushdir ('file/mutate')
      cmd = 'rm -rf ../'//trim(work)//' ; mkdir -p ../'//trim(work)//
     &      ' ; cp water2.xyz ../'//trim(work)//'/ ; cp '//
     &      '175_water_vsoft_n15_ti_l00.key ../'//trim(work)//
     &      '/water2.key'
      call execute_command_line (cmd)
      call popdir
c
c     add the window schedule and the request to stop the run
c
      call pushdir ('file/'//trim(work))
      call tnist_append ('water2.key','integrator verlet')
      call tnist_append ('water2.key','lambda-window 1.0')
      call tnist_append ('water2.key','lambda-eqratio 0.0')
      call tnist_append ('water2.key','ti-nstepavg '//nstepavg)
      iend = freeunit ()
      open (unit=iend,file='water2.end',status='new')
      close (unit=iend)
      call popdir
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine ticount  --  count TI file data records  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "ticount" returns the number of data records in a block average
c     file, counting only the nonblank records that follow the line
c     labeling the columns of the block averages
c
c
      subroutine ticount (fname,nrow)
      implicit none
      integer nrow
      integer iti
      integer freeunit
      logical header
      character*(*) fname
      character*240 record
c
c
      nrow = 0
      header = .false.
      iti = freeunit ()
      open (unit=iti,file=fname,status='old')
   10 continue
      read (iti,20,end=30)  record
   20 format (a240)
      if (.not. header) then
         if (index(record,'Index') .ne. 0)  header = .true.
      else if (record .ne. ' ') then
         nrow = nrow + 1
      end if
      goto 10
   30 continue
      close (unit=iti)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine tirowlmda  --  lambda column of a TI row  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "tirowlmda" returns the lambda recorded in the requested data
c     record of a block average file, counting only the nonblank
c     rows that follow the line labeling the columns
c
c
      subroutine tirowlmda (fname,irow,lam)
      implicit none
      integer irow
      integer nrow,idx
      integer iti
      integer freeunit
      logical header
      real*8 lam
      character*(*) fname
      character*240 record
c
c
      lam = -1.0d0
      nrow = 0
      header = .false.
      iti = freeunit ()
      open (unit=iti,file=fname,status='old')
   10 continue
      read (iti,20,end=30)  record
   20 format (a240)
      if (.not. header) then
         if (index(record,'Index') .ne. 0)  header = .true.
      else if (record .ne. ' ') then
         nrow = nrow + 1
         if (nrow .eq. irow) then
            read (record,*,err=30,end=30)  idx,lam
            goto 30
         end if
      end if
      goto 10
   30 continue
      close (unit=iti)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine tiwipe  --  remove a TI file if present  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "tiwipe" deletes the named block average file when it exists,
c     so the test always starts from a predictable file name
c
c
      subroutine tiwipe (fname)
      implicit none
      integer iti
      integer freeunit
      logical exist
      character*(*) fname
c
c
      inquire (file=fname,exist=exist)
      if (exist) then
         iti = freeunit ()
         open (unit=iti,file=fname,status='old')
         close (unit=iti,status='delete')
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine tidynstep  --  TI work after a dynamics step  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "tidynstep" does for one dynamics step what the integrator
c     and the "dynamic" loop do between them, sampling the current
c     lambda window and then advancing at a window boundary
c
c
      subroutine tidynstep (istep)
      implicit none
      integer istep
c
c
      call lmdawindyn (istep)
      call lmdawinstep (istep)
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
c     "resetti" lays out a schedule of equal length lambda windows
c     and sets the scalar state to deterministic unit test defaults,
c     bypassing "inittidyn" so the accumulation tests do not depend
c     on how the step budget is divided
c
c
      subroutine resetti (nbin,nstepavg,window,nequil)
      implicit none
      integer nbin,nstepavg,window,nequil
      integer i
      integer, allocatable :: steps(:)
c
c
c     hand every window the same number of steps
c
      allocate (steps(nbin))
      do i = 1, nbin
         steps(i) = window
      end do
      call resettisteps (nbin,nstepavg,steps,nequil)
      deallocate (steps)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine resettisteps  --  uneven TI test state  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "resettisteps" is "resetti" for a schedule whose windows have
c     different lengths, given as the step count of each window; the
c     equilibration ratio is taken from the first window so the test
c     geometry stays easy to state
c
c
      subroutine resettisteps (nbin,nstepavg,steps,nequil)
      use dlmda
      use mutant
      use thrmint
      implicit none
      integer nbin,nstepavg,nequil
      integer steps(*)
      integer i
c
c
c     put the sublambda maps on the power law branch, which is
c     plain arithmetic and needs no molecular system
c
      use_rel = .false.
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
      nlmdawin = nbin
      tinstepavg = nstepavg
      lmdawinratio = 0.0d0
      if (steps(1) .gt. 0)  lmdawinratio = dble(nequil) / dble(steps(1))
      dedl = 0.0d0
c
c     clear any previous allocation and size the sample buffer
c
      if (allocated(tidedllist))  deallocate (tidedllist)
      allocate (tidedllist(nstepavg))
      do i = 1, nstepavg
         tidedllist(i) = 0.0d0
      end do
c
c     build the default evenly spaced schedule, which also sets
c     "lmdawin" and "lambda" to the first window
c
      call setlmdasched (0,.false.)
c
c     turn the requested window lengths into step boundaries, then
c     size the first window and the accumulators the same way the
c     production code does
c
      if (allocated(lmdawinend))  deallocate (lmdawinend)
      allocate (lmdawinend(nbin))
      lmdawinend(1) = steps(1)
      do i = 2, nbin
         lmdawinend(i) = lmdawinend(i-1) + steps(i)
      end do
      call setlmdawin
      call settiblocks
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
      use mutant
      use thrmint
      implicit none
c
c
      use_ti = .false.
      lmdawin = 0
      tinbcount = 0
      tinbsave = 0
      tinbtot = 0
      lambda = 1.0d0
      if (allocated(lmdawinend))  deallocate (lmdawinend)
      if (allocated(tidedllist))  deallocate (tidedllist)
      if (allocated(tilmdahist))  deallocate (tilmdahist)
      if (allocated(tilmdadedl))  deallocate (tilmdadedl)
      if (allocated(tilmdadedlstd))  deallocate (tilmdadedlstd)
      if (allocated(lmdawinlist))  deallocate (lmdawinlist)
      if (allocated(lmdawinfrac))  deallocate (lmdawinfrac)
      return
      end
