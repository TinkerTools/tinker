c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  mplpot.i  --  specifics of atomic multipole functions  ##
c     ##                                                         ##
c     #############################################################
c
c
c     m2scale   scale factor for 1-2 multipole energy interactions
c     m3scale   scale factor for 1-3 multipole energy interactions
c     m4scale   scale factor for 1-4 multipole energy interactions
c     m5scale   scale factor for 1-5 multipole energy interactions
c
c
      real*8 m2scale,m3scale
      real*8 m4scale,m5scale
      common /mplpot/ m2scale,m3scale,m4scale,m5scale
