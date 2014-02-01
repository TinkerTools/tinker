c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  rxnfld.i  --  reaction field matrix elements and indices  ##
c     ##                                                            ##
c     ################################################################
c
c
c     b1    first reaction field matrix element array
c     b2    second reaction field matrix element array
c     ijk   indices into the reaction field element arrays
c
c
      integer ijk
      real*8 b1,b2
      common /rxnfld/ b1(40,13),b2(40,13),ijk(0:5,0:5,0:5)
