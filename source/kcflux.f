c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2020 by Chengwen Liu & Jay W. Ponder  ##
c     ##                 All Rights Reserved                  ##
c     ##########################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module kcflux -- charge flux term forcefield parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxncfb   maximum number of bond stretch charge flux entries
c     maxncfa   maximum number of angle bend charge flux entries
c
c     cflb      charge flux over stretching of a bond length
c     cfla      charge flux over bending of a bond angle
c     cflab     charge flux over asymmetric bond within an angle
c     kcfb      string of atom classes for bond stretch charge flux
c     kcfa      string of atom classes for angle bend charge flux
c
c
      module kcflux
      implicit none
      integer maxncfb
      integer maxncfa
      real*8, allocatable :: cflb(:)
      real*8, allocatable :: cfla(:,:)
      real*8, allocatable :: cflab(:,:)
      character*8, allocatable :: kcfb(:)
      character*12, allocatable :: kcfa(:)
      save
      end
