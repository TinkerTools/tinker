c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  chrono.i  --  timing statistics for the current program  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     cputim   elapsed cpu time in seconds since start of program
c
c
      real*8 cputim
      common /chrono/ cputim
