c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  syntrn.i  --  definition of synchronous transit path  ##
c     ##                                                        ##
c     ############################################################
c
c
c     t       value of the path coordinate (0=reactant, 1=product)
c     pm      path coordinate for extra point in quadratic transit
c     xmin1   reactant coordinates as array of optimization variables
c     xmin2   product coordinates as array of optimization variables
c     xm      extra coordinate set for quadratic synchronous transit
c
c
      real*8 t,pm,xmin1,xmin2,xm
      common /syntrn/ t,pm,xmin1(maxvar),xmin2(maxvar),xm(maxvar)
