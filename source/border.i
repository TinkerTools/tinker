c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  border.i  --  bond orders for a conjugated pisystem  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     pbpl   pi-bond orders for bonds in "planar" pisystem
c     pnpl   pi-bond orders for bonds in "nonplanar" pisystem
c
c
      real*8 pbpl,pnpl
      common /border/ pbpl(maxpib),pnpl(maxpib)
