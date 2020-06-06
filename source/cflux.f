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
c     bflx    bond stretching charge flux constant (electrons/Ang)
c     aflx    angle bending charge flux constant (electrons/radian)
c     abflx   asymmetric stretch charge flux constant (electrons/Ang)
c
c
      module cflux
      implicit none
      real*8, allocatable :: bflx(:)
      real*8, allocatable :: aflx(:,:)
      real*8, allocatable :: abflx(:,:)
      save
      end
