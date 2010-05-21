c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine fftsetup  --  setup 3-D Fast Fourier transform  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "fftsetup" does initialization for a 3-D FFT via a single 3-D
c     transform (FFTW) or three separate 1-D transforms (FFTPACK)
c
c
      subroutine fftsetup
      implicit none
      include 'sizes.i'
      include 'fft.i'
      include 'pme.i'
!$    integer ifront,iback
!$    integer error,iguess
!$    integer nprocs
!$    integer omp_get_num_procs
c
c
c     initialization of Fast Fourier transform using FFTW
c
!$    if (ffttyp .eq. 'FFTW') then
!$       ifront = -1
!$       iback = 1
!$       error = 0
!$       iguess = 0
!$       nprocs = 1
!$       nprocs = omp_get_num_procs ()
!$       call dfftw_init_threads (error)
!$       call dfftw_plan_with_nthreads (nprocs)
!$       call dfftw_plan_dft_3d (planf,nfft1,nfft2,nfft3,qgrid,
!$   &                              qgrid,ifront,iguess)
!$       call dfftw_plan_dft_3d (planb,nfft1,nfft2,nfft3,qgrid,
!$   &                              qgrid,iback,iguess)
!$    else
c
c     initialization of Fast Fourier transform using FFTPACK
c
         call cffti (nfft1,table(1,1),iprime(1,1))
         call cffti (nfft2,table(1,2),iprime(1,2))
         call cffti (nfft3,table(1,3),iprime(1,3))
!$    end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine fftfront  --  forward Fast Fourier transform  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "fftfront" performs a 3-D FFT forward transform via a single
c     3-D transform or three separate 1-D transforms
c
c
      subroutine fftfront
      implicit none
      include 'sizes.i'
      include 'fft.i'
      include 'pme.i'
      integer i,j,k
      real*8 work(2,maxfft)
c
c
c     perform a single 3-D forward transform using FFTW
c
!$    if (ffttyp .eq. 'FFTW') then
!$       call dfftw_execute_dft (planf,qgrid,qgrid)
!$    else
c
c     perform three 1-D forward transforms using FFTPACK
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
!$    end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine fftback  --  backward Fast Fourier transform  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "fftback" performs a 3-D FFT backward transform via a single
c     3-D transform or three separate 1-D transforms
c
c
      subroutine fftback
      implicit none
      include 'sizes.i'
      include 'fft.i'
      include 'pme.i'
      integer i,j,k
      real*8 work(2,maxfft)
c
c
c     perform a single 3-D backward transform using FFTW
c
!$    if (ffttyp .eq. 'FFTW') then
!$       call dfftw_execute_dft (planb,qgrid,qgrid)
!$    else
c
c     perform three 1-D backward transforms using FFTPACK
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
!$    end if
      return
      end
