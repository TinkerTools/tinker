c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2022 by Moses Chung, Zhi Wang & Jay Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module expol  --  exch-polarization in current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nexpol       total number of exch polarization sites in system
c     kpep         exchange polarization spring constant at each site
c     prepep       exchange polarization prefactor at each site
c     dmppep       exchange polarization damping alpha at each site
c     polscale     scale matrix for use in exchange polarization
c     invpolscale  scale matrix inverse for exchange polarization
c     lpep         exchange polarization logical at each site
c
c
      module expol
      implicit none
      integer nexpol
      real*8, allocatable :: kpep(:)
      real*8, allocatable :: prepep(:)
      real*8, allocatable :: dmppep(:)
      real*8, allocatable :: polscale(:,:,:)
      real*8, allocatable :: invpolscale(:,:,:)
      logical, allocatable :: lpep(:)
      save
      end
