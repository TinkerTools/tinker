c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  kopdst.i  --  forcefield parameters for out-plane distance  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     maxnopd   maximum number of out-of-plane distance entries
c
c     opds      force constant parameters for out-of-plane distance
c     kopd      string of atom classes for out-of-plane distance
c
c
      integer maxnopd
      parameter (maxnopd=500)
      real*8 opds
      character*16 kopd
      common /kopdst/ opds(maxnopd),kopd(maxnopd)
