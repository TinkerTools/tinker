c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  opdist.i  --  out-of-plane distances in current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     opdk      force constant values for out-of-plane distance
c     nopdist   total number of out-of-plane distances in the system
c     iopb      numbers of the atoms in each out-of-plane distance
c
c
      integer nopdist,iopd
      real*8 opdk
      common /opdist/ opdk(maxatm),nopdist,iopd(4,maxatm)
