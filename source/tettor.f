c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module tettor  --  tetratorsions in the current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     ntettor  total number of tetratorsions in the system
c     itettor  numbers of the atoms in each tetratorsion
c
c
      module tettor
      implicit none
      integer ntettor
      integer, allocatable :: itettor(:,:)
      save
      end
