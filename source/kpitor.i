c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  kpitor.i  --  forcefield parameters for pi-orbit torsions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnpt   maximum number of pi-orbital torsion parameter entries
c
c     ptcon    force constant parameters for pi-orbital torsions
c     kpt      string of atom classes for pi-orbital torsion terms
c
c
      integer maxnpt
      parameter (maxnpt=500)
      real*8 ptcon
      character*8 kpt
      common /kpitor/ ptcon(maxnpt),kpt(maxnpt)
