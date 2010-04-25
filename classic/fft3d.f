c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine fftsetup  --  3-D FFT initialization  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "fftsetup" does the initialization for a 3-D FFT via
c     three separate 1-D initializations
c
c
      subroutine fftsetup
      implicit none
      include 'sizes.i'
      include 'pme.i'
c
c
c     perform initialization along X, Y and Z directions
c
      call cffti (nfft1,table(1,1),iprime(1,1))
      call cffti (nfft2,table(1,2),iprime(1,2))
      call cffti (nfft3,table(1,3),iprime(1,3))
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine fftfront  --  3-D FFT forward transform  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "fftfront" does a 3-D FFT forward transform via three
c     separate 1-D transformations
c
c
      subroutine fftfront
      implicit none
      include 'sizes.i'
      include 'pme.i'
      integer i,j,k
      real*8 work(2,maxfft)
c
c
c     perform forward transform along X, Y and Z directions
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               work(1,i) = qgrid(1,i,j,k)
               work(2,i) = qgrid(2,i,j,k)
            end do
            call cfftf (nfft1,work,table(1,1),iprime(1,1))
            do i = 1, nfft1
               qgrid(1,i,j,k) = work(1,i)
               qgrid(2,i,j,k) = work(2,i)
            end do
         end do
      end do
      do k = 1, nfft3
         do i = 1, nfft1
            do j = 1, nfft2
               work(1,j) = qgrid(1,i,j,k)
               work(2,j) = qgrid(2,i,j,k)
            end do
            call cfftf (nfft2,work,table(1,2),iprime(1,2))
            do j = 1, nfft2
               qgrid(1,i,j,k) = work(1,j)
               qgrid(2,i,j,k) = work(2,j)
            end do
         end do
      end do
      do i = 1, nfft1
         do j = 1, nfft2
            do k = 1, nfft3
               work(1,k) = qgrid(1,i,j,k)
               work(2,k) = qgrid(2,i,j,k)
            end do
            call cfftf (nfft3,work,table(1,3),iprime(1,3))
            do k = 1, nfft3
               qgrid(1,i,j,k) = work(1,k)
               qgrid(2,i,j,k) = work(2,k)
            end do
         end do
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine fftback  --  3-D FFT backward transform  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "fftback" does a 3-D FFT backward transform via three
c     separate 1-D transformations
c
c
      subroutine fftback
      implicit none
      include 'sizes.i'
      include 'pme.i'
      integer i,j,k
      real*8 work(2,maxfft)
c
c
c     perform backward transform along X, Y and Z directions
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               work(1,i) = qgrid(1,i,j,k)
               work(2,i) = qgrid(2,i,j,k)
            end do
            call cfftb (nfft1,work,table(1,1),iprime(1,1))
            do i = 1, nfft1
               qgrid(1,i,j,k) = work(1,i)
               qgrid(2,i,j,k) = work(2,i)
            end do
         end do
      end do
      do k = 1, nfft3
         do i = 1, nfft1
            do j = 1, nfft2
               work(1,j) = qgrid(1,i,j,k)
               work(2,j) = qgrid(2,i,j,k)
            end do
            call cfftb (nfft2,work,table(1,2),iprime(1,2))
            do j = 1, nfft2
               qgrid(1,i,j,k) = work(1,j)
               qgrid(2,i,j,k) = work(2,j)
            end do
         end do
      end do
      do i = 1, nfft1
         do j = 1, nfft2
            do k = 1, nfft3
               work(1,k) = qgrid(1,i,j,k)
               work(2,k) = qgrid(2,i,j,k)
            end do
            call cfftb (nfft3,work,table(1,3),iprime(1,3))
            do k = 1, nfft3
               qgrid(1,i,j,k) = work(1,k)
               qgrid(2,i,j,k) = work(2,k)
            end do
         end do
      end do
      return
      end
