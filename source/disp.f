c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module disp  --  damped dispersion in current structure  ##        
c     ##                                                           ##
c     ###############################################################
c
c
c     ndisp     total number of dispersion sites in the system
c     idisp     number of the atom for each dispersion site
c     csixpr    pairwise sum of C6 dispersion coefficients
c     csix      C6 dispersion coefficient value at each site
c     adisp     alpha dispersion damping value at each site
c
c
      module disp
      implicit none
      integer ndisp
      integer, allocatable :: idisp(:)
      real*8 csixpr
      real*8, allocatable :: csix(:)
      real*8, allocatable :: adisp(:)
      save
      end
