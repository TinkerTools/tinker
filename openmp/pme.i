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
c     maxprime   maximum number of prime factors of FFT dimension
c     maxtable   maximum size of the FFT table array
c     maxgrid    maximum dimension of the PME charge grid array
c
c     bsmod1     B-spline moduli along the a-axis direction
c     bsmod2     B-spline moduli along the b-axis direction
c     bsmod3     B-spline moduli along the c-axis direction
c     table      intermediate array used by the FFT calculation
c     qgrid      values on the particle mesh Ewald charge grid
c     qfac       prefactors for particle mesh Ewald charge grid
c     thetai1    B-spline coefficients along the a-axis
c     thetai2    B-spline coefficients along the b-axis
c     thetai3    B-spline coefficients along the c-axis
c     nfft1      number of grid points along the a-axis direction
c     nfft2      number of grid points along the b-axis direction
c     nfft3      number of grid points along the c-axis direction
c     bsorder    order of the PME B-spline approximation
c     iprime     prime factorization of each FFT dimension
c     initgrid   initialization of grid values for B-spline
c
c
      integer maxorder,maxprime
      integer maxtable,maxgrid
      parameter (maxorder=10)
      parameter (maxprime=15)
      parameter (maxtable=4*maxfft)
      parameter (maxgrid=2*(maxfft/2)+1)
      integer nfft1,nfft2,nfft3
      integer bsorder,iprime
      integer initgrid
      real*8 bsmod1,bsmod2,bsmod3
      real*8 table,qfac
      real*8 thetai1,thetai2,thetai3
      real*8, pointer :: qgrid(:,:,:,:)
      common /pme/ bsmod1(maxfft),bsmod2(maxfft),
     &             bsmod3(maxfft),table(maxtable,3),
     &		   qfac(maxfft,maxfft,maxfft),
     &             thetai1(4,maxorder,maxatm),
     &		   thetai2(4,maxorder,maxatm),
     &             thetai3(4,maxorder,maxatm),
     &             nfft1,nfft2,nfft3,bsorder,
     &             iprime(maxprime,3),initgrid(3,maxatm),
     &             qgrid
