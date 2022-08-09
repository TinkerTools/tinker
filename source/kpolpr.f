c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2022  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module kpolpr  --  special Thole forcefield parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxnpp   maximum number of special pair polarization entries
c
c     thlpr    Thole damping values for special polarization pairs
c     thdpr    Thole direct damping for special polarization pairs
c     kppr     string of atom classes for special polarization pairs
c
c
      module kpolpr
      implicit none
      integer maxnpp
      real*8, allocatable :: thlpr(:)
      real*8, allocatable :: thdpr(:)
      character*8, allocatable :: kppr(:)
      save
      end
