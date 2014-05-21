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
c     maxele   maximum number of elements from periodic table
c
c     elemnt   atomic symbol for each chemical element
c
c
      module ptable
      implicit none
      integer maxele
      parameter (maxele=112)
      character*3 elemnt(maxele)
      save
      end
