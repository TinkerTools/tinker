c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  stodyn.i  --  frictional coefficients for SD trajectory  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     friction    global frictional coefficient for exposed particle
c     fgamma      atomic frictional coefficients for each atom
c     use_sdarea  logical flag to use surface area friction scaling
c
c
      real*8 friction,fgamma
      logical use_sdarea
      common /stodyn/ friction,fgamma(maxatm),use_sdarea
