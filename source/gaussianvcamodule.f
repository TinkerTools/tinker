c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2023 by MKJ Chung, MJ Schnieders & JW Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module gaussianvcamodule  --  description of gaussian  ##
c     ##                                                         ##
c     #############################################################
c
c
c     v   gaussian volume
c     a   gaussian exponent
c     c   gaussian center
c
c
      module gaussianvcamodule
      implicit none
      type :: GaussianVca
      real*8 v
      real*8 a
      real*8 c(3)
      end type GaussianVca
      end module gaussianvcamodule
