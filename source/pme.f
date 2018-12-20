c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module pme  --  values for particle mesh Ewald summation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nfft1      current number of PME grid points along a-axis
c     nfft2      current number of PME grid points along b-axis
c     nfft3      current number of PME grid points along c-axis
c     nefft1     number of grid points along electrostatic a-axis
c     nefft2     number of grid points along electrostatic b-axis
c     nefft3     number of grid points along electrostatic c-axis
c     ndfft1     number of grid points along dispersion a-axis
c     ndfft2     number of grid points along dispersion b-axis
c     ndfft3     number of grid points along dispersion c-axis
c     bsorder    current order of the PME B-spline values
c     bseorder   order of the electrostatic PME B-spline values
c     bsdorder   order of the dispersion PME B-spline values
c     igrid      initial Ewald grid values for B-spline
c     bsmod1     B-spline moduli along the a-axis direction
c     bsmod2     B-spline moduli along the b-axis direction
c     bsmod3     B-spline moduli along the c-axis direction
c     bsbuild    B-spline derivative coefficient temporary storage
c     thetai1    B-spline coefficients along the a-axis
c     thetai2    B-spline coefficients along the b-axis
c     thetai3    B-spline coefficients along the c-axis
c     qgrid      values on the particle mesh Ewald grid
c     qfac       prefactors for the particle mesh Ewald grid
c
c
      module pme
      implicit none
      integer nfft1,nfft2,nfft3
      integer nefft1,nefft2,nefft3
      integer ndfft1,ndfft2,ndfft3
      integer bsorder,bseorder,bsdorder
      integer, allocatable :: igrid(:,:)
      real*8, allocatable :: bsmod1(:)
      real*8, allocatable :: bsmod2(:)
      real*8, allocatable :: bsmod3(:)
      real*8, allocatable :: bsbuild(:,:)
      real*8, allocatable :: thetai1(:,:,:)
      real*8, allocatable :: thetai2(:,:,:)
      real*8, allocatable :: thetai3(:,:,:)
      real*8, allocatable :: qgrid(:,:,:,:)
      real*8, allocatable :: qfac(:,:,:)
      save
      end
