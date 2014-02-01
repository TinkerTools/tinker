c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  ptable.i  --  atomic symbols for the chemical elements  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     elemnt   atomic symbol for each chemical element
c
c
      character*3 elemnt
      common /ptable/ elemnt(maxele)
