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
c     aslow   RESPA secondary slow acceleration of each atom
c     afast   RESPA secondary fast acceleration of each atom
c
c
      module moldyn
      implicit none
      real*8, allocatable :: v(:,:)
      real*8, allocatable :: a(:,:)
      real*8, allocatable :: aalt(:,:)
      real*8, allocatable :: aslow(:,:)
      real*8, allocatable :: afast(:,:)
      save
      end
