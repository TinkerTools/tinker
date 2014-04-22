c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2014 by Chao Lv & Jay William Ponder  ##
c     ##                 All Rights Reserved                  ##
c     ##########################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  angtor.i  --  angle-torsions in the current structure  ##
c     ##                                                         ##
c     #############################################################
c
c
c     kant      1-, 2- and 3-fold angle-torsion force constants
c     nangtor   total number of angle-torsion interactions
c     iat       torsion and angle numbers used in angle-torsion
c
c
      integer nangtor,iat
      real*8 kant
      common /angtor/ kant(6,maxtors),nangtor,iat(3,maxtors)
