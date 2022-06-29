c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung, Zhi Wang & Jay Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module kexpl  --  exch polarization forcefield parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     pepk       exch polarization spring constant for each atom type
c     peppre     exch polarization constant prefactor for each atom type
c     pepdmp     exch polarization damping alpha for each atom type
c     pepl       exch polarization logical for each atom type
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
