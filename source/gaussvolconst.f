c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2023 by M. Chung, M. Schnieders & Jay Ponder  ##
c     ##                    All Rights Reserved                       ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module gaussvolconst  --  constants used in GaussVol  ##
c     ##                                                        ##
c     ############################################################
c
c
c     MAX_ORDER
c     KFC
c     PFC
c     MIN_GVOL
c     ANG3
c     VOLMINA
c     VOLMINB
c
c
      module gaussvolconst
      implicit none
      integer MAX_ORDER
      real*8 KFC,PFC
      real*8 MIN_GVOL
      real*8 ANG3,VOLMINA,VOLMINB
      parameter (MAX_ORDER = 16)
      parameter (KFC = 2.2269859253d0)
      parameter (PFC = 2.5d0)
      parameter (MIN_GVOL = tiny(1.0d0))
      parameter (ANG3 = 1.0d0)
      parameter (VOLMINA = 0.01d0 * ANG3)
      parameter (VOLMINB = 0.1d0 * ANG3)
      save
      end
