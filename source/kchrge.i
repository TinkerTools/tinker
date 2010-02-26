c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  kchrge.i  --  forcefield parameters for partial charges  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     chg   partial charge parameters for each atom type
c
c
      real*8 chg
      common /kchrge/ chg(maxtyp)
