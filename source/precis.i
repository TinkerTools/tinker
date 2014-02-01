c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  precis.i  --  values of machine precision tolerances  ##
c     ##                                                        ##
c     ############################################################
c
c
c     tiny    the smallest positive floating point value
c     small   the smallest relative floating point spacing
c     huge    the largest relative floating point spacing
c
c
      real*8 tiny,small,huge
      common /precis/ tiny,small,huge
