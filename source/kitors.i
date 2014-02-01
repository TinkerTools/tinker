c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  kitors.i  --  forcefield parameters for improper torsions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnti   maximum number of improper torsion parameter entries
c
c     ti1      torsional parameters for improper 1-fold rotation
c     ti2      torsional parameters for improper 2-fold rotation
c     ti3      torsional parameters for improper 3-fold rotation
c     kti      string of atom classes for improper torsional parameters
c
c
      integer maxnti
      parameter (maxnti=500)
      real*8 ti1,ti2,ti3
      character*16 kti
      common /kitors/ ti1(2,maxnti),ti2(2,maxnti),ti3(2,maxnti),
     &                kti(maxnti)
