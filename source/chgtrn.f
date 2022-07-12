c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module chgtrn  --  charge transfer in current structure  ##        
c     ##                                                           ##
c     ###############################################################
c
c
c     nct       total number of dispersion sites in the system
c     chgct     charge for charge transfer at each multipole site
c     dmpct     charge transfer damping factor at each multipole site
c
c
      module chgtrn
      implicit none
      integer nct
      real*8, allocatable :: chgct(:)
      real*8, allocatable :: dmpct(:)
      save
      end
