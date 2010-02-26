c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  bitor.i  --  bitorsions within the current structure  ##
c     ##                                                        ##
c     ############################################################
c
c
c     nbitor  total number of bitorsions in the system
c     ibitor  numbers of the atoms in each bitorsion
c
c
      integer nbitor,ibitor
      common /bitor/ nbitor,ibitor(5,maxbitor)
