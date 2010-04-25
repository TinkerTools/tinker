c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  bndpot.i  --  specifics of bond stretch functional forms  ##
c     ##                                                            ##
c     ################################################################
c
c
c     cbnd      cubic coefficient in bond stretch potential
c     qbnd      quartic coefficient in bond stretch potential
c     bndunit   convert bond stretch energy to kcal/mole
c     bndtyp    type of bond stretch potential energy function
c
c
      real*8 cbnd,qbnd
      real*8 bndunit
      character*8 bndtyp
      common /bndpot/ cbnd,qbnd,bndunit,bndtyp
