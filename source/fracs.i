c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  fracs.i  --  atom distances to molecular center of mass  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     xfrac   fractional coordinate along a-axis of center of mass
c     yfrac   fractional coordinate along b-axis of center of mass
c     zfrac   fractional coordinate along c-axis of center of mass
c
c
      real*8 xfrac,yfrac,zfrac
      common /fracs/ xfrac(maxatm),yfrac(maxatm),zfrac(maxatm)
