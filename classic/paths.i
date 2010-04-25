c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  paths.i  --  parameters for Elber reaction path method  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     p0       reactant Cartesian coordinates as variables
c     p1       product Cartesian coordinates as variables
c     pmid     midpoint between the reactant and product
c     pvect    vector connecting the reactant and product
c     pstep    step per cycle along reactant-product vector
c     pzet     current projection on reactant-product vector
c     pnorm    length of the reactant-product vector
c     acoeff   transformation matrix 'A' from Elber paper
c     gc       gradients of the path constraints
c
c
      real*8 p0,p1,pmid
      real*8 pvect,pstep
      real*8 pzet,pnorm
      real*8 acoeff,gc
      common /paths/ p0(maxvar),p1(maxvar),pmid(maxvar),pvect(maxvar),
     &               pstep(maxvar),pzet(maxvar),pnorm,acoeff(7,7),
     &               gc(maxvar,7)
