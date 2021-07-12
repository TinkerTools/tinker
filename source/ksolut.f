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
c     solrad   implicit solvation radius value for each atom class
c
c
      module ksolut
      implicit none
      real*8, allocatable :: solrad(:)
      save
      end
