c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module params  --  force field parameter file contents  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nprm      number of nonblank lines in the parameter file
c     prmline   contents of each individual parameter file line
c
c
      module params
      use sizes
      implicit none
      integer nprm
      character*120 prmline(maxprm)
      save
      end
