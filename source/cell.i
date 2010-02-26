c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  cell.i  --  periodic boundaries using replicated cells  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     xcell    length of the a-axis of the complete replicated cell
c     ycell    length of the b-axis of the complete replicated cell
c     zcell    length of the c-axis of the complete replicated cell
c     xcell2   half the length of the a-axis of the replicated cell
c     ycell2   half the length of the b-axis of the replicated cell
c     zcell2   half the length of the c-axis of the replicated cell
c     ncell    total number of cell replicates for periodic boundaries
c     icell    offset along axes for each replicate periodic cell
c
c
      integer ncell,icell
      real*8 xcell,ycell,zcell
      real*8 xcell2,ycell2,zcell2
      common /cell/ xcell,ycell,zcell,xcell2,ycell2,zcell2,ncell,
     &              icell(3,maxcell)
