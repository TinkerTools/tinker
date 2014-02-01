c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  kopbnd.i  --  forcefield parameters for out-of-plane bend  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnopb   maximum number of out-of-plane bending entries
c
c     opbn      force constant parameters for out-of-plane bending
c     kopb      string of atom classes for out-of-plane bending
c
c
      integer maxnopb
      parameter (maxnopb=500)
      real*8 opbn
      character*16 kopb
      common /kopbnd/ opbn(maxnopb),kopb(maxnopb)
