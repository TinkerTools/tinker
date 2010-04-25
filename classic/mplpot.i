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
c     m2scale   factor by which 1-2 multipole interactions are scaled
c     m3scale   factor by which 1-3 multipole interactions are scaled
c     m4scale   factor by which 1-4 multipole interactions are scaled
c     m5scale   factor by which 1-5 multipole interactions are scaled
c
c
      real*8 m2scale,m3scale
      real*8 m4scale,m5scale
      common /mplpot/ m2scale,m3scale,m4scale,m5scale
