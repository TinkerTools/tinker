c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module kopdst  --  out-of-plane distance forcefield params  ##
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
      module kopdst
      implicit none
      integer maxnopd
      real*8, allocatable :: opds(:)
      character*16, allocatable :: kopd(:)
      save
      end
