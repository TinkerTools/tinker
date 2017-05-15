c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module kpolr  --  polarizability forcefield parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     polr   dipole polarizability parameters for each atom type
c     athl   Thole polarizability damping value for each atom type
c     pgrp   connected types in polarization group of each atom type
c
c
      module kpolr
      implicit none
      integer, allocatable :: pgrp(:,:)
      real*8, allocatable :: polr(:,:)
      real*8, allocatable :: athl(:)
      save
      end
