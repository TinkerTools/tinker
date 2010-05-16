c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine clock  --  find elapsed time for current job  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "clock" determines elapsed CPU time in seconds since the
c     start of the job
c
c     note only one of the various implementations below should
c     be activated by removing comment characters
c
c
      subroutine clock (seconds)
      implicit none
      real*8 seconds
c
c
c     code for Fortran standard "cpu_time" intrinsic function
c
c     real time
c     call cpu_time (time)
c     seconds = dble(time)
c
c     code for the Unix standard "etime" intrinsic function
c
c     real etime,times(2)
c     seconds = dble(etime(times))
c
c     code for IBM xlf compiler on AIX appends an underscore
c
c     real etime_,times(2)
c     seconds = dble(etime_(times))
c
c     code for Hewlett-Packard HP-UX only gives wall clock time;
c     for the HP systems the +E1 compiler option is required
c
      real secnds,start
      logical first
      save first,start
      data first  / .true. /
      if (first) then
         first = .false.
         start = secnds (0.0)
      end if
      seconds = secnds (start)
      return
      end
