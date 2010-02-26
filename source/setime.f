c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine setime  --  initialize elapsed CPU time clock  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "setime" initializes the elapsed interval CPU timer
c
c
      subroutine setime
      implicit none
      include 'chrono.i'
c
c
c     initialize interval at elapsed CPU time for current job
c
      call clock (cputim)
      return
      end
