c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2020 by Chengwen Liu & Jay W. Ponder  ##
c     ##                 All Rights Reserved                  ##
c     ##########################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module cflux  --  charge flux terms in current system  ##
c     ##                                                         ##
c     #############################################################
c
c
c     jb      charge flux over unit bond length
c     b0      equilibrium bond length for charge flux
c     ja      charge flux over unit angle value
c     a0      equilibrium angle value for charge flux
c     pcflx   change in charge on each atom due to charge flux
c
c
      module cflux
      use sizes
      implicit none
      real*8, allocatable :: jb(:)
      real*8, allocatable :: b0(:)
      real*8, allocatable :: theta0(:)
      real*8, allocatable :: bp0(:,:)
      real*8, allocatable :: jbp(:,:)
      real*8, allocatable :: jtheta(:,:)
      real*8, allocatable :: pcflx(:)
      save
      end
