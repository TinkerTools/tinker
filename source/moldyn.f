c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module moldyn  --  MD trajectory velocity & acceleration  ##
c     ##                                                            ##
c     ################################################################
c
c
c     v       current velocity of each atom along the x,y,z-axes
c     a       current acceleration of each atom along x,y,z-axes
c     aalt    alternate acceleration of each atom along x,y,z-axes
c
c
      module moldyn
      use sizes
      implicit none
      real*8 v(3,maxatm)
      real*8 a(3,maxatm)
      real*8 aalt(3,maxatm)
      save
      end
