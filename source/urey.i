c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  urey.i  --  Urey-Bradley interactions in the structure  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     uk      Urey-Bradley force constants (kcal/mole/Ang**2)
c     ul      ideal 1-3 distance values in Angstroms
c     nurey   total number of Urey-Bradley terms in the system
c     iury    numbers of the atoms in each Urey-Bradley interaction
c
c
      integer nurey,iury
      real*8 uk,ul
      common /urey/ uk(maxang),ul(maxang),nurey,iury(3,maxang)
