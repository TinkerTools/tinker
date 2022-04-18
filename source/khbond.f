c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module khbond  --  H-bonding term forcefield parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxnhb   maximum number of hydrogen bonding pair entries
c
c     radhb    radius parameter for hydrogen bonding pairs
c     epshb    well depth parameter for hydrogen bonding pairs
c     khb      string of atom types for hydrogen bonding pairs
c
c
      module khbond
      implicit none
      integer maxnhb
      real*8, allocatable :: radhb(:)
      real*8, allocatable :: epshb(:)
      character*8, allocatable :: khb(:)
      save
      end
