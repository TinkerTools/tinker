c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  pitors.i  --  pi-orbital torsions in the current structure  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     kpit      2-fold pi-orbital torsional force constants
c     npitors   total number of pi-orbital torsional interactions
c     ipit      numbers of the atoms in each pi-orbital torsion
c
c
      integer npitors,ipit
      real*8 kpit
      common /pitors/ kpit(maxtors),npitors,ipit(6,maxtors)
