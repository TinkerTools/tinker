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
      use sizes
      implicit none
      integer maxncfb
      integer maxncfa
      parameter (maxncfb=2000)
      parameter (maxncfa=2000)
      real*8 cflb(maxncfb)
      real*8 cfla(2,maxncfa)
      real*8 cflab(2,maxncfa)
      character*8 kcfb(maxncfb)
      character*12 kcfa(maxncfa)
      save
      end
