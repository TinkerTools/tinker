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
c     maxnbcf   maximum number of charge flux bond entries
c     maxnacf   maximum number of charge flux angle entries
c
c     jb        charge flux over unit bond length
c     beq       equilibrium bond length in charge flux
c     ja        charge flux over unit angle
c     theta0l   equilibrium angle in charge flux
c
c
      module kcflux
      use sizes
      implicit none
      integer maxnbcf
      integer maxnacf
      parameter (maxnbcf=2000)
      parameter (maxnacf=2000)
      real*8 jbnd(maxnbcf)
      real*8 beq(maxnbcf)
      real*8 theta0l(maxnacf)
      real*8 bp0l(2,maxnacf)
      real*8 jbpl(2,maxnacf)
      real*8 jthetal(2,maxnacf)
      character*8 kcfb(maxnbcf)
      character*12 kcfa(maxnacf)
      save
      end
