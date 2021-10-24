c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module inform  --  program I/O and flow control values  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxask    maximum number of queries for interactive input
c
c     digits    decimal places output for energy and coordinates
c     iprint    steps between status printing (0=no printing)
c     iwrite    steps between coordinate saves (0=no saves)
c     isend     steps between socket communication (0=no sockets)
c     debug     debug printing level (0=none to 9=detailed output)
c     verbose   logical flag to turn on extra information printing
c     silent    logical flag to turn off all information printing
c     holdup    logical flag to wait for carriage return on exit
c     abort     logical flag to stop execution at next chance
c
c
      module inform
      implicit none
      integer maxask
      parameter (maxask=5)
      integer digits,iprint
      integer iwrite,isend
      integer debug
      logical verbose,silent
      logical holdup,abort
      save
      end
