c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module kurybr  --  Urey-Bradley term forcefield parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     maxnu   maximum number of Urey-Bradley parameter entries
c
c     ucon    force constant parameters for Urey-Bradley terms
c     dst13   ideal 1-3 distance parameters for Urey-Bradley terms
c     ku      string of atom classes for Urey-Bradley terms
c
c
      module kurybr
      implicit none
      integer maxnu
      real*8, allocatable :: ucon(:)
      real*8, allocatable :: dst13(:)
      character*12, allocatable :: ku(:)
      save
      end
