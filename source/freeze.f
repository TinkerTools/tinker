c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module freeze  --  definition of holonomic constraints  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nrat         number of holonomic distance constraints to apply
c     nratx        number of atom group holonomic constraints to apply
c     nwat         number of rigid three- and four-site water molecules
c     nwat4        number of rigid planar four-site water molecules
c     iratx        group number for each group holonomic constraint
c     kratx        spatial constraint type (1=plane, 2=line, 3=point)
c     irat         atom numbers of atoms in each holonomic constraint
c     iwat         atom numbers of O and H atoms in each rigid water
c     iwat4        atom numbers involved in each rigid four-site water
c     rateps       convergence tolerance for holonomic constraints
c     krat         ideal distance value for each holonomic constraint
c     kwat         ideal distances for O-H and H-H in rigid water
c     kwat4        geometry scaling values for rigid four-site water
c     use_freeze   logical flag to set use of holonomic contraints
c     frzimage     flag to use minimum image for holonomic constraint
c
c
      module freeze
      implicit none
      integer nrat,nratx
      integer nwat,nwat4
      integer, allocatable :: iratx(:)
      integer, allocatable :: kratx(:)
      integer, allocatable :: irat(:,:)
      integer, allocatable :: iwat(:,:)
      integer, allocatable :: iwat4(:,:)
      real*8 rateps
      real*8, allocatable :: krat(:)
      real*8, allocatable :: kwat(:,:)
      real*8, allocatable :: kwat4(:,:)
      logical use_freeze
      logical, allocatable :: frzimage(:)
      save
      end
