c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module chgpen  --  charge penetration in current structure  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     ncp       total number of charge penetration sites in system
c     pcore     number of core electrons at each multipole site
c     pval      number of valence electrons at each multipole site
c     palpha    charge penetration damping at each multipole site
c
c
      module chgpen
      implicit none
      integer ncp
      real*8, allocatable :: pcore(:)
      real*8, allocatable :: pval(:)
      real*8, allocatable :: palpha(:)
      save
      end
