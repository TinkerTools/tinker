c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module ptable  --  atomic symbols for chemical elements  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     elemnt   atomic symbol for each chemical element
c
c
      module ptable
      use sizes
      implicit none
      character*3 elemnt(maxele)
      save
      end
