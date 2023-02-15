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
c     nion      total number of partial charge sites in the system
c     iion      number of the atom for each partial charge site
c     jion      neighbor generation site to use for each atom
c     kion      cutoff switching site to use for each atom
c     pchg      current partial charge value for each atom (e-)
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
