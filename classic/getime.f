c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine getime  --  get elapsed CPU time in seconds  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "getime" gets elapsed CPU time in seconds for an interval
c
c
      subroutine getime (elapsed)
      implicit none
      include 'chrono.i'
      real*8 elapsed
c
c
c     elapsed time for the interval is the current total CPU
c     time minus the total time at the start of the interval
c
      call clock (elapsed)
      elapsed = elapsed - cputim
      return
      end
