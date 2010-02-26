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
c     use the Unix-style "etime" intrinsic function, this code
c     works for all compilers except those noted below
c
c     real etime,times(2)
c     seconds = dble(etime(times))
c
c     code for Windows PC's under Compaq Visual Fortran
c
c     real time
c     call cpu_time (time)
c     seconds = dble(time)
c
c     code for the IBM AIX xlf compiler appends an underscore
c
      real etime_,times(2)
      seconds = dble(etime_(times))
c
c     code for Hewlett-Packard HP-UX only gives wall clock time;
c     for the HP systems the +E1 compiler option is required
c
c     real secnds,start
c     logical initial
c     save initial,start
c     data initial  / .true. /
c     if (initial) then
c        initial = .false.
c        start = secnds (0.0)
c     end if
c     seconds = secnds (start)
      return
      end
