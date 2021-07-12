c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module charge  --  partial charges in current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     nion      total number of partial charges in system
c     iion      number of the atom site for each partial charge
c     jion      neighbor generation site for each partial charge
c     kion      cutoff switching site for each partial charge
c     pchg      current atomic partial charge values (e-)
c     pchg0     original partial charge values for charge flux
c
c
      module charge
      implicit none
      integer nion
      integer, allocatable :: iion(:)
      integer, allocatable :: jion(:)
      integer, allocatable :: kion(:)
      real*8, allocatable :: pchg(:)
      real*8, allocatable :: pchg0(:)
      save
      end
