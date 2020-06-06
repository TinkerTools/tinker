c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module atmlst  --  bond and angle local geometry indices  ##
c     ##                                                            ##
c     ################################################################
c
c
c     bndlist   numbers of the bonds involving each atom
c     anglist   numbers of the angles centered on each atom
c     balist    numbers of the bonds comprising each angle
c
c
      module atmlst
      implicit none
      integer, allocatable :: bndlist(:,:)
      integer, allocatable :: anglist(:,:)
      integer, allocatable :: balist(:,:)
      save
      end
