c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module kmulti  --  atomic multipole forcefield parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnmp   maximum number of atomic multipole parameter entries
c
c     multip   atomic monopole, dipole and quadrupole values
c     mpaxis   type of local axis definition for atomic multipoles
c     kmp      string of atom types for atomic multipoles
c
c
      module kmulti
      implicit none
      integer maxnmp
      real*8, allocatable :: multip(:,:)
      character*8, allocatable :: mpaxis(:)
      character*16, allocatable :: kmp(:)
      save
      end
