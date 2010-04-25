c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  rotate.i  --  molecule partitions for rotation of a bond  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nrot        total number of atoms moving when bond rotates
c     rot         atom numbers of atoms moving when bond rotates
c     use_short   logical flag governing use of shortest atom list
c
c
      integer nrot,rot
      logical use_short
      common /rotate/ nrot,rot(maxatm),use_short
