c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  bond.i  --  covalent bonds in the current structure  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     bk      bond stretch force constants (kcal/mole/Ang**2)
c     bl      ideal bond length values in Angstroms
c     nbond   total number of bond stretches in the system
c     ibnd    numbers of the atoms in each bond stretch
c
c
      integer nbond,ibnd
      real*8 bk,bl
      common /bond/ bk(maxbnd),bl(maxbnd),nbond,ibnd(2,maxbnd)
