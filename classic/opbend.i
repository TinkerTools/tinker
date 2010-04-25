c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  opbend.i  --  out-of-plane bends in the current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     opbk      force constant values for out-of-plane bending
c     nopbend   total number of out-of-plane bends in the system
c     iopb      bond angle numbers used in out-of-plane bending
c
c
      integer nopbend,iopb
      real*8 opbk
      common /opbend/ opbk(maxang),nopbend,iopb(maxang)
