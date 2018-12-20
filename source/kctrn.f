c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module kctrn  --  charge transfer forcefield parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     ctchg     charge transfer magnitude for each atom class
c     ctdmp     alpha charge transfer parameter for each atom class
c
c
      module kctrn
      implicit none
      real*8, allocatable :: ctchg(:)
      real*8, allocatable :: ctdmp(:)
      save
      end
