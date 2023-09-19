c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2023 by MKJ Chung, MJ Schnieders & JW Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module gaussvolconst  --  constants used for GaussVol  ##
c     ##                                                         ##
c     #############################################################
c
c
c     max_order    max order of computing gaussian overlaps
c     kfc          conversion factor from spheres to gaussians
c     min_gvol     minimum threshold to include overlap in recursion
c     rad_offset   radius offset for finite diff volume calculation
c     volmina      volume cutoff A for switching function
c     volminb      volume cutoff B for switching function
c
c
      module gaussvolconst
      implicit none
      integer max_order
      real*8 kfc
      real*8 min_gvol
      real*8 rad_offset
      real*8 volmina
      real*8 volminb
      parameter (max_order = 16)
      parameter (kfc = 2.2269859253d0)
      parameter (min_gvol = tiny(1.0))
      parameter (rad_offset = 0.005d0)
      parameter (volmina = 0.01d0)
      parameter (volminb = 0.1d0)
      end module gaussvolconst
