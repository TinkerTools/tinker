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
c     pc0      reactant Cartesian coordinates as variables
c     pc1      product Cartesian coordinates as variables
c     pvect    vector connecting the reactant and product
c     pstep    step per cycle along reactant-product vector
c     pzet     current projection on reactant-product vector
c     gc       gradient of the path constraints
c     pnorm    length of the reactant-product vector
c     acoeff   transformation matrix 'A' from Elber algorithm
c
c
      real*8, pointer :: pc0(:)
      real*8, pointer :: pc1(:)
      real*8, pointer :: pvect(:)
      real*8, pointer :: pstep(:)
      real*8, pointer :: pzet(:)
      real*8, pointer :: gc(:,:)
      real*8 pnorm,acoeff
      common /paths/ pc0,pc1,pvect,pstep,pzet,gc,pnorm,acoeff(7,7)
