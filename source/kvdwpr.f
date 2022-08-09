c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module kvdwpr  --  special pair vdw forcefield parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnvp   maximum number of special pair van der Waals entries
c
c     radpr    radius parameter for special van der Waals pairs
c     epspr    well depth parameter for special van der Waals pairs
c     kvpr     string of atom classes for special van der Waals pairs
c
c
      module kvdwpr
      implicit none
      integer maxnvp
      real*8, allocatable :: radpr(:)
      real*8, allocatable :: epspr(:)
      character*8, allocatable :: kvpr(:)
      save
      end
