c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  tortor.i  --  torsion-torsions in the current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     ntortor   total number of torsion-torsion interactions
c     itt       atoms and parameter indices for torsion-torsion
c
c
      integer ntortor,itt
      common /tortor/ ntortor,itt(3,maxbitor)
