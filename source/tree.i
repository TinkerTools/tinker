c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  tree.i  --  potential smoothing and search tree levels  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxpss   maximum number of potential smoothing levels
c
c     etree    energy reference value at the top of the tree
c     ilevel   smoothing deformation value at each tree level
c     nlevel   number of levels of potential smoothing used
c
c
      integer maxpss
      parameter (maxpss=500)
      integer nlevel
      real*8 etree,ilevel
      common /tree/ etree,ilevel(0:maxpss),nlevel
