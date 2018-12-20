c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module krepl  --  Pauli repulsion forcefield parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     prsiz     Pauli repulsion size value for each atom class
c     prdmp     alpha Pauli repulsion parameter for each atom class
c     prele     number of valence electrons for each atom class
c
c
      module krepl
      implicit none
      real*8, allocatable :: prsiz(:)
      real*8, allocatable :: prdmp(:)
      real*8, allocatable :: prele(:)
      save
      end
