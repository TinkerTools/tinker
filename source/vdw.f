c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module vdw  --  van der Waals terms in current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nvdw       total number van der Waals sites in the system
c     ivdw       number of the atom for each van der Waals site
c     jvdw       index into the vdw parameter matrix for each atom
c     mvdw       index into the vdw parameter matrix for each class
c     ired       attached atom from which reduction factor is applied
c     kred       value of reduction factor parameter for each atom
c     xred       reduced x-coordinate for each atom in the system
c     yred       reduced y-coordinate for each atom in the system
c     zred       reduced z-coordinate for each atom in the system
c     radvdw     van der Waals radius for each atom in the system
c     epsvdw     van der Waals well depth for each atom in the system
c     radmin     minimum energy distance for each atom class pair
c     epsilon    well depth parameter for each atom class pair
c     radmin4    minimum energy distance for 1-4 interaction pairs
c     epsilon4   well depth parameter for 1-4 interaction pairs
c     radhbnd    minimum energy distance for hydrogen bonding pairs
c     epshbnd    well depth parameter for hydrogen bonding pairs
c
c
      module vdw
      implicit none
      integer nvdw
      integer, allocatable :: ivdw(:)
      integer, allocatable :: jvdw(:)
      integer, allocatable :: mvdw(:)
      integer, allocatable :: ired(:)
      real*8, allocatable :: kred(:)
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: radvdw(:)
      real*8, allocatable :: epsvdw(:)
      real*8, allocatable :: radmin(:,:)
      real*8, allocatable :: epsilon(:,:)
      real*8, allocatable :: radmin4(:,:)
      real*8, allocatable :: epsilon4(:,:)
      real*8, allocatable :: radhbnd(:,:)
      real*8, allocatable :: epshbnd(:,:)
      save
      end
