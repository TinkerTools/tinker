c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  strtor.i  --  stretch-torsions in the current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     kst       1-, 2- and 3-fold stretch-torsion force constants
c     nstrtor   total number of stretch-torsion interactions
c     ist       torsion and bond numbers used in stretch-torsion
c
c
      integer nstrtor,ist
      real*8 kst
      common /strtor/ kst(9,maxtors),nstrtor,ist(4,maxtors)
