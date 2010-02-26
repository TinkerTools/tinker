c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  moldyn.i  --  velocity and acceleration on MD trajectory  ##
c     ##                                                            ##
c     ################################################################
c
c
c     v       current velocity of each atom along the x,y,z-axes
c     a       current acceleration of each atom along x,y,z-axes
c     aold    previous acceleration of each atom along x,y,z-axes
c
c
      real*8 v,a,aold
      common /moldyn/ v(3,maxatm),a(3,maxatm),aold(3,maxatm)
