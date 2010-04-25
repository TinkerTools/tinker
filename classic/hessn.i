c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  hessn.i  --  Cartesian Hessian elements for a single atom  ##
c     ##                                                             ##
c     #################################################################
c
c
c     hessx   Hessian elements for x-component of current atom
c     hessy   Hessian elements for y-component of current atom
c     hessz   Hessian elements for z-component of current atom
c
c
      real*8 hessx,hessy,hessz
      common /hessn/ hessx(3,maxatm),hessy(3,maxatm),hessz(3,maxatm)
