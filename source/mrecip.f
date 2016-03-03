c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module mrecip  --  reciprocal PME for permanent multipoles  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     cmp     Cartesian permenent multipoles as polytensor vector
c     fmp     fractional permanent multipoles as polytensor vector
c     fphi    permanent multipole potential and field from FFT
c
c
      module mrecip
      implicit none
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: fphi(:,:)
      save
      end
