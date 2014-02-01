c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  kurybr.i  --  forcefield parameters for Urey-Bradley terms  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     maxnu   maximum number of Urey-Bradley parameter entries
c
c     ucon    force constant parameters for Urey-Bradley terms
c     dst13   ideal 1-3 distance parameters for Urey-Bradley terms
c     ku      string of atom classes for Urey-Bradley terms
c
c
      integer maxnu
      parameter (maxnu=2000)
      real*8 ucon,dst13
      character*12 ku
      common /kurybr/ ucon(maxnu),dst13(maxnu),ku(maxnu)
