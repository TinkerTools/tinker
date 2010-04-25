c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  kmulti.i  --  forcefield parameters for atomic multipoles  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnmp   maximum number of atomic multipole parameter entries
c
c     multip   atomic monopole, dipole and quadrupole values
c     mpaxis   type of local axis definition for atomic multipoles
c     kmp      string of atom types for atomic multipoles
c
c
      integer maxnmp
      parameter (maxnmp=2000)
      real*8 multip
      character*8 mpaxis
      character*16 kmp
      common /kmulti/ multip(13,maxnmp),mpaxis(maxnmp),kmp(maxnmp)
