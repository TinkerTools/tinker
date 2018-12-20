c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module kcpen  --  charge penetration forcefield parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     cpele     valence electron magnitude for each atom class
c     cpalp     alpha charge penetration parameter for each atom class
c
c
      module kcpen
      implicit none
      real*8, allocatable :: cpele(:)
      real*8, allocatable :: cpalp(:)
      save
      end
