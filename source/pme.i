c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  pme.i  --  values for particle mesh Ewald summation  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     maxorder   maximum order of the B-spline approximation
c
c     bsmod1     B-spline moduli along the a-axis direction
c     bsmod2     B-spline moduli along the b-axis direction
c     bsmod3     B-spline moduli along the c-axis direction
c     thetai1    B-spline coefficients along the a-axis
c     thetai2    B-spline coefficients along the b-axis
c     thetai3    B-spline coefficients along the c-axis
c     qgrid      values on the particle mesh Ewald charge grid
c     qfac       prefactors for particle mesh Ewald charge grid
c     nfft1      number of grid points along the a-axis direction
c     nfft2      number of grid points along the b-axis direction
c     nfft3      number of grid points along the c-axis direction
c     bsorder    order of the PME B-spline approximation
c     igrid      initial Ewald charge grid values for B-spline
c
c
      integer maxorder
      parameter (maxorder=10)
      integer nfft1,nfft2,nfft3
      integer bsorder,igrid
      real*8 bsmod1,bsmod2,bsmod3
      real*8, pointer :: thetai1(:,:,:)
      real*8, pointer :: thetai2(:,:,:)
      real*8, pointer :: thetai3(:,:,:)
      real*8, pointer :: qgrid(:,:,:,:)
      real*8, pointer :: qfac(:,:,:)
      common /pme/ bsmod1(maxfft),bsmod2(maxfft),bsmod3(maxfft),
     &             thetai1,thetai2,thetai3,qgrid,qfac,
     &             nfft1,nfft2,nfft3,bsorder,igrid(3,maxatm)
