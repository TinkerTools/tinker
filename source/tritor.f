c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module tritor  --  tritorsions in the current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     ntritor  total number of tritorsions in the system
c     itritor  numbers of the atoms in each tritorsion
c
c
      module tritor
      implicit none
      integer ntritor
      integer, allocatable :: itritor(:,:)
      save
      end
