c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung, Zhi Wang & Jay Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module expol  --  exch polarization for current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     nexpol    total number of exch polarization sites in the system
c     kpep      exchange polarization spring constant at each site
c     prepep    exchange polarization constant prefactor at each site
c     dmppep    exchange polarization damping alpha at each site
c     lpep      exchange polarization logical at each site
c     scrtyp    type of screening
c
c
      module expol
      implicit none
      integer nexpol
      real*8, allocatable :: kpep(:)
      real*8, allocatable :: prepep(:)
      real*8, allocatable :: dmppep(:)
      logical, allocatable :: lpep(:)
      character*3 scrtyp
      save
      end