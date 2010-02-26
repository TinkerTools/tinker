c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program timer  --  timer for Cartesian energy functions  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "timer" measures the CPU time required for file reading and
c     parameter assignment, potential energy computation, energy
c     and gradient computation, and Hessian matrix evaluation
c
c
      program timer
      implicit none
      include 'sizes.i'
      include 'cutoff.i'
      include 'hescut.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,ncalls,next
      integer hinit(3,maxatm)
      integer hstop(3,maxatm)
      integer hindex(maxhess)
      real*8 value,energy,elapsed
      real*8 derivs(3,maxatm)
      real*8 h(maxhess)
      real*8 hdiag(3,maxatm)
      logical exist,query
      logical dohessian
      character*1 answer
      character*120 record
      character*120 string
c
c
c     read in the molecular system to be timed
c
      call initial
      call getxyz
c
c     get the number of calculation cycles to perform
c
      ncalls = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  ncalls
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter Desired Number of Repetitions [1] :  ',$)
         read (input,30)  ncalls
   30    format (i10)
      end if
      if (ncalls .eq. 0)  ncalls = 1
c
c     decide whether to include timing of Hessian evaluations
c
      dohessian = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,40)
   40    format (/,' Include Timing for Hessian Evaluations [Y] :  ',$)
         read (input,50)  record
   50    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  dohessian = .false.
c
c     get the timing for setup of the calculation
c
      call setime
      call mechanic
      if (use_list)  call nblist
      call getime (elapsed)
      write (iout,60)  elapsed
   60 format (/,' Computation Set-up :',f15.3,' Seconds')
c
c     set a large Hessian cutoff and turn off extra printing
c
      hesscut = 100.0d0
      verbose = .false.
c
c     run the potential energy only timing experiment
c
      call setime
      do i = 1, ncalls
         value = energy ()
      end do
      call getime (elapsed)
      write (iout,70)  elapsed,ncalls
   70 format (/,' Potential Energy :  ',f15.3,' Seconds for',
     &           i6,' Evaluations')
c
c     run the energy and gradient timing experiment
c
      call setime
      do i = 1, ncalls
         call gradient (value,derivs)
      end do
      call getime (elapsed)
      write (iout,80)  elapsed,ncalls
   80 format (/,' Energy & Gradient : ',f15.3,' Seconds for',
     &           i6,' Evaluations')
c
c     run the Hessian matrix only timing experiment
c
      if (dohessian) then
         call setime
         do i = 1, ncalls
            call hessian (h,hinit,hstop,hindex,hdiag)
         end do
         call getime (elapsed)
         write (iout,90)  elapsed,ncalls
   90    format (/,' Hessian Matrix :    ',f15.3,' Seconds for',
     &              i6,' Evaluations')
      end if
c
c     repeat the potential energy only timing experiment
c
      call setime
      do i = 1, ncalls
         value = energy ()
      end do
      call getime (elapsed)
      write (iout,100)  elapsed,ncalls
  100 format (/,' Potential Energy :  ',f15.3,' Seconds for',
     &           i6,' Evaluations')
c
c     repeat the energy and gradient timing experiment
c
      call setime
      do i = 1, ncalls
         call gradient (value,derivs)
      end do
      call getime (elapsed)
      write (iout,110)  elapsed,ncalls
  110 format (/,' Energy & Gradient : ',f15.3,' Seconds for',
     &           i6,' Evaluations')
c
c     repeat the Hessian matrix only timing experiment
c
      if (dohessian) then
         call setime
         do i = 1, ncalls
            call hessian (h,hinit,hstop,hindex,hdiag)
         end do
         call getime (elapsed)
         write (iout,120)  elapsed,ncalls
  120    format (/,' Hessian Matrix :    ',f15.3,' Seconds for',
     &              i6,' Evaluations')
      end if
c
c     perform any final tasks before program exit
c
      call final
      end
