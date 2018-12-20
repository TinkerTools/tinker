c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module kdsp  --  damped dispersion forcefield parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     dspsix   C6 dispersion coefficient for each atom class
c     dspdmp   alpha dispersion parameter for each atom class
c
c
      module kdsp
      implicit none
      real*8, allocatable :: dspsix(:)
      real*8, allocatable :: dspdmp(:)
      save
      end
