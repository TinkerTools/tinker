c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2020  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module ksolut  --  solvation term forcefield parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     pbr      Poisson-Boltzmann radius value for each atom type
c     csr      ddCOSMO solvation radius value for each atom type
c     gkr      Generalized Kirkwood radius value for each atom type
c
c
      module ksolut
      implicit none
      real*8, allocatable :: pbr(:)
      real*8, allocatable :: csr(:)
      real*8, allocatable :: gkr(:)
      save
      end
