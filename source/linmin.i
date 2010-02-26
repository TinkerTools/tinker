c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  linmin.i  --  parameters for line search minimization  ##
c     ##                                                         ##
c     #############################################################
c
c
c     stpmin   minimum step length in current line search direction
c     stpmax   maximum step length in current line search direction
c     cappa    stringency of line search (0=tight < cappa < 1=loose)
c     slpmax   projected gradient above which stepsize is reduced
c     angmax   maximum angle between search direction and -gradient
c     intmax   maximum number of interpolations during line search
c
c
      integer intmax
      real*8 stpmin,stpmax
      real*8 cappa,slpmax
      real*8 angmax
      common /linmin/ stpmin,stpmax,cappa,slpmax,angmax,intmax
