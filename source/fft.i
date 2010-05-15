c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  fft.i  --  values and options for Fast Fourier transform  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxtable   maximum size of the FFT table intermediate array
c     maxprime   maximum number of prime factors of FFT dimension
c
c     table      intermediate array used by the FFT routine (fftpack)
c     iprime     prime factorization of each FFT dimension (fftpack)
c     planf      pointer to forward transform data structure (fftw)
c     planb      pointer to backward transform data structure (fftw)
c     ffttyp     type of FFT package; currently FFTPACK or FFTW
c
c
      integer maxtable
      integer maxprime
      parameter (maxtable=4*maxfft)
      parameter (maxprime=15)
      integer iprime
      integer*8 planf
      integer*8 planb
      real*8 table
      character*7 ffttyp
      common /fft/ table(maxtable,3),iprime(maxprime,3),planf,planb,
     &             ffttyp
