c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module bound  --  periodic boundary condition controls  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     polycut       cutoff distance for infinite polymer nonbonds
c     polycut2      square of infinite polymer nonbond cutoff
c     use_bounds    flag to use periodic boundary conditions
c     use_replica   flag to use replicates for periodic system
c     use_polymer   flag to mark presence of infinite polymer
c     use_wrap      flag to wrap coordinates in output files
c
c
      module bound
      implicit none
      real*8 polycut
      real*8 polycut2
      logical use_bounds
      logical use_replica
      logical use_polymer
      logical use_wrap
      save
      end
