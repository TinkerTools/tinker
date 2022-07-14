c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2022 by Moses Chung, Zhi Wang & Jay Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module kexpl  --  exch-polarization forcefield parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     pepk     exchange-polarization spring constant for atom classes
c     peppre   exchange-polarization prefactor for atom classes
c     pepdmp   exchange-polarization damping alpha for atom classes
c     pepl     exchange-polarization logical flag for atom classes
c
c
      module kexpl
      implicit none
      real*8, allocatable :: pepk(:)
      real*8, allocatable :: peppre(:)
      real*8, allocatable :: pepdmp(:)
      logical, allocatable :: pepl(:)
      save
      end
