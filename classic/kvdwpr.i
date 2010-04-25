c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  kvdwpr.i  --  forcefield parameters for special vdw terms  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnvp   maximum number of special van der Waals pair entries
c
c     radpr    radius parameter for special van der Waals pairs
c     epspr    well depth parameter for special van der Waals pairs
c     kvpr     string of atom classes for special van der Waals pairs
c
c
      integer maxnvp
      parameter (maxnvp=500)
      real*8 radpr,epspr
      character*8 kvpr
      common /kvdwpr/ radpr(maxnvp),epspr(maxnvp),kvpr(maxnvp)
