c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module ewald  --  Ewald summation parameters and options  ##
c     ##                                                            ##
c     ################################################################
c
c
c     aewald     current value of Ewald convergence coefficient
c     aeewald    Ewald convergence coefficient for electrostatics
c     adewald    Ewald convergence coefficient for dispersion
c     boundary   Ewald boundary condition; none, tinfoil or vacuum
c
c
      module ewald
      implicit none
      real*8 aewald
      real*8 aeewald
      real*8 adewald
      character*7 boundary
      save
      end
