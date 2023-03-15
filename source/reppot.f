c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module reppot  --  repulsion functional form details  ##
c     ##                                                        ##
c     ############################################################
c
c
c     r2scale    scale factor for 1-2 repulsion energy interactions
c     r3scale    scale factor for 1-3 repulsion energy interactions
c     r4scale    scale factor for 1-4 repulsion energy interactions
c     r5scale    scale factor for 1-5 repulsion energy interactions
c     reptyp     type of repulsion potential energy function
c
c
      module reppot
      implicit none
      real*8 r2scale
      real*8 r3scale
      real*8 r4scale
      real*8 r5scale
      character*8 reptyp
      save
      end
