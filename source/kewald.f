c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kewald  --  setup for particle mesh Ewald sum  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kewald" assigns particle mesh Ewald parameters and options
c     for a periodic system
c
c
      subroutine kewald
      use atoms
      use bound
      use boxes
      use chunks
      use ewald
      use fft
      use inform
      use iounit
      use keys
      use limits
      use openmp
      use pme
      implicit none
      integer maxpower
      integer maxfft
      parameter (maxpower=63)
      parameter (maxfft=864)
      integer i,k,next
      integer nbig,minfft
      integer iefft1,idfft1
      integer iefft2,idfft2
      integer iefft3,idfft3
      integer multi(maxpower)
      real*8 delta,rmax
      real*8 edens,ddens
      character*20 keyword
      character*240 record
      character*240 string
c
c     PME grid size must be even with factors of only 2, 3 and 5
c
      data multi  /   2,   4,   6,   8,  10,  12,  16,  18,  20,
     &               24,  30,  32,  36,  40,  48,  50,  54,  60,
     &               64,  72,  80,  90,  96, 100, 108, 120, 128,
     &              144, 150, 160, 162, 180, 192, 200, 216, 240,
     &              250, 256, 270, 288, 300, 320, 324, 360, 384,
     &              400, 432, 450, 480, 486, 500, 512, 540, 576,
     &              600, 640, 648, 720, 750, 768, 800, 810, 864 /
c
c
c     return if Ewald summation is not being used
c
      if (.not.use_ewald .and. .not.use_dewald)  return
c
c     set default values for Ewald options and parameters
c
      ffttyp = 'FFTPACK'
      if (nthread .gt. 1)  ffttyp = 'FFTW'
      boundary = 'TINFOIL'
      bseorder = 5
      bsdorder = 4
      edens = 1.2d0
      ddens = 0.8d0
      aeewald = 0.4d0
      adewald = 0.4d0
      minfft = 16
c
c     estimate an optimal value for the Ewald coefficient
c
      if (use_ewald)  call ewaldcof (aeewald,ewaldcut)
      if (use_dewald)  call ewaldcof (adewald,dewaldcut)
c
c     set the system extent for nonperiodic Ewald summation
c
      if (.not. use_bounds) then
         call extent (rmax)
         xbox = 2.0d0 * (rmax+max(ewaldcut,dewaldcut))
         ybox = xbox
         zbox = xbox
         alpha = 90.0d0
         beta = 90.0d0
         gamma = 90.0d0
         orthogonal = .true.
         call lattice
         boundary = 'NONE'
         edens = 0.7d0
         ddens = 0.7d0
      end if
c
c     set defaults for electrostatic and dispersion grid sizes
c
      nefft1 = 0
      nefft2 = 0
      nefft3 = 0
      ndfft1 = 0
      ndfft2 = 0
      ndfft3 = 0
c
c     get default grid counts from periodic system dimensions
c
      delta = 1.0d-8
      iefft1 = int(xbox*edens-delta) + 1
      iefft2 = int(ybox*edens-delta) + 1
      iefft3 = int(zbox*edens-delta) + 1
      idfft1 = int(xbox*ddens-delta) + 1
      idfft2 = int(ybox*ddens-delta) + 1
      idfft3 = int(zbox*ddens-delta) + 1
c
c     search keywords for Ewald summation commands
c
      do i = 1, nkey
         record = keyline(i)
         next = 1
         call upcase (record)
         call gettext (record,keyword,next)
         string = record(next:240)
         if (keyword(1:12) .eq. 'FFT-PACKAGE ') then
            call getword (record,ffttyp,next)
         else if (keyword(1:12) .eq. 'EWALD-ALPHA ') then
            read (string,*,err=20,end=20)  aeewald
         else if (keyword(1:13) .eq. 'DEWALD-ALPHA ') then
            read (string,*,err=20,end=20)  adewald
         else if (keyword(1:15) .eq. 'EWALD-BOUNDARY ') then
            boundary = 'VACUUM'
         else if (keyword(1:9) .eq. 'PME-GRID ') then
            iefft1 = 0
            iefft2 = 0
            iefft3 = 0
            read (string,*,err=10,end=10)  iefft1,iefft2,iefft3
   10       continue
            if (iefft2 .eq. 0)  iefft2 = iefft1
            if (iefft3 .eq. 0)  iefft3 = iefft1
         else if (keyword(1:10) .eq. 'DPME-GRID ') then
            idfft1 = 0
            idfft2 = 0
            idfft3 = 0
            read (string,*,err=15,end=15)  idfft1,idfft2,idfft3
   15       continue
            if (idfft2 .eq. 0)  idfft2 = idfft1
            if (idfft3 .eq. 0)  idfft3 = idfft1
         else if (keyword(1:10) .eq. 'PME-ORDER ') then
            read (string,*,err=20,end=20)  bseorder
         else if (keyword(1:11) .eq. 'DPME-ORDER ') then
            read (string,*,err=20,end=20)  bsdorder
         end if
   20    continue
      end do
c
c     determine electrostatic grid size from allowed values
c
      if (use_ewald) then
         nefft1 = maxfft
         nefft2 = maxfft
         nefft3 = maxfft
         do i = maxpower, 1, -1
            k = multi(i)
            if (k .le. maxfft) then
               if (k .ge. iefft1)  nefft1 = k
               if (k .ge. iefft2)  nefft2 = k
               if (k .ge. iefft3)  nefft3 = k
            end if
         end do
         if (nefft1 .lt. minfft)  nefft1 = minfft
         if (nefft2 .lt. minfft)  nefft2 = minfft
         if (nefft3 .lt. minfft)  nefft3 = minfft
      end if
c
c     determine dispersion grid size from allowed values
c
      if (use_dewald) then
         ndfft1 = maxfft
         ndfft2 = maxfft
         ndfft3 = maxfft
         do i = maxpower, 1, -1
            k = multi(i)
            if (k .le. maxfft) then
               if (k .ge. idfft1)  ndfft1 = k
               if (k .ge. idfft2)  ndfft2 = k
               if (k .ge. idfft3)  ndfft3 = k
            end if
         end do
         if (ndfft1 .lt. minfft)  ndfft1 = minfft
         if (ndfft2 .lt. minfft)  ndfft2 = minfft
         if (ndfft3 .lt. minfft)  ndfft3 = minfft
      end if
c
c     increase electrostatic grid to match dispersion grid
c
      if (use_ewald .and. use_dewald) then
         if (nefft1.lt.ndfft1 .or. nefft2.lt.ndfft2
     &          .or. nefft3.lt.ndfft3) then
            nefft1 = max(nefft1,ndfft1)
            nefft2 = max(nefft2,ndfft2)
            nefft3 = max(nefft3,ndfft3)
            write (iout,30)
   30       format (/,' KEWALD  --  Setting Electrostatic PME',
     &                 ' Grid to Match Dispersion Grid')
         end if
      end if
c
c     check the particle mesh Ewald grid dimensions
c
      nbig = max(nefft1,nefft2,nefft3,ndfft1,ndfft2,ndfft3)
      if (nbig .gt. maxfft) then
         write (iout,40)
   40    format (/,' KEWALD  --  PME Grid Size Too Large;',
     &              ' Increase MAXFFT')
         call fatal
      end if
      if (use_ewald .and. (nefft1.lt.iefft1.or.
     &       nefft2.lt.iefft2.or.nefft3.lt.iefft3)) then
         write (iout,50)
   50    format (/,' KEWALD  --  Warning, Small Electrostatic',
     &              'PME Grid Size')
      end if
      if (use_dewald .and. (ndfft1.lt.idfft1.or.
     &       ndfft2.lt.idfft2.or.ndfft3.lt.idfft3)) then
         write (iout,60)
   60    format (/,' KEWALD  --  Warning, Small Dispersion',
     &              'PME Grid Size')
      end if
c
c     set maximum sizes for PME grid and B-spline order
c
      nfft1 = max(nefft1,ndfft1)
      nfft2 = max(nefft2,ndfft2)
      nfft3 = max(nefft3,ndfft3)
      bsorder = max(bseorder,bsdorder)
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(bsmod1))  deallocate (bsmod1)
      if (allocated(bsmod2))  deallocate (bsmod2)
      if (allocated(bsmod3))  deallocate (bsmod3)
      if (allocated(bsbuild))  deallocate (bsbuild)
      if (allocated(thetai1))  deallocate (thetai1)
      if (allocated(thetai2))  deallocate (thetai2)
      if (allocated(thetai3))  deallocate (thetai3)
      if (allocated(qgrid))  deallocate (qgrid)
      if (allocated(qfac))  deallocate (qfac)
      if (allocated(pmetable))  deallocate (pmetable)
      allocate (bsmod1(nfft1))
      allocate (bsmod2(nfft2))
      allocate (bsmod3(nfft3))
      allocate (bsbuild(bsorder,bsorder))
      allocate (thetai1(4,bsorder,n))
      allocate (thetai2(4,bsorder,n))
      allocate (thetai3(4,bsorder,n))
      allocate (qgrid(2,nfft1,nfft2,nfft3))
      allocate (qfac(nfft1,nfft2,nfft3))
      allocate (pmetable(n,6*nthread))
c
c     print a message listing some of the Ewald parameters
c
      if (verbose) then
         write (iout,70)
   70    format (/,' Particle Mesh Ewald Parameters :',
     &           //,5x,'Type',16x,'Ewald Alpha',4x,'Grid',
     &              ' Dimensions',4x,'Spline Order',/)
         if (use_ewald) then
            write (iout,80)  aeewald,nefft1,nefft2,nefft3,bseorder
   80       format (3x,'Electrostatics',9x,f8.4,5x,3i5,7x,i5)
         end if
         if (use_dewald) then
            write (iout,90)  adewald,ndfft1,ndfft2,ndfft3,bsdorder
   90       format (3x,'Dispersion',13x,f8.4,5x,3i5,7x,i5)
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ewaldcof  --  estimation of Ewald coefficient  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ewaldcof" finds an Ewald coefficient such that all terms
c     beyond the specified cutoff distance will have a value less
c     than a specified tolerance
c
c
      subroutine ewaldcof (alpha,cutoff)
      implicit none
      integer i,k
      real*8 alpha,cutoff,eps
      real*8 x,xlo,xhi,y
      real*8 ratio,erfc
      external erfc
c
c
c     set the tolerance value; use of 1.0d-8 instead of 1.0d-6
c     gives large coefficients that ensure gradient continuity
c
      eps = 1.0d-8
c
c     get approximate value from cutoff and tolerance
c
      ratio = eps + 1.0d0
      x = 0.5d0
      i = 0
      do while (ratio .ge. eps)
         i = i + 1
         x = 2.0d0 * x
         y = x * cutoff
         ratio = erfc(y) / cutoff
      end do
c
c     use a binary search to refine the coefficient
c
      k = i + 60
      xlo = 0.0d0
      xhi = x
      do i = 1, k
         x = (xlo+xhi) / 2.0d0
         y = x * cutoff
         ratio = erfc(y) / cutoff
         if (ratio .ge. eps) then
            xlo = x
         else
            xhi = x
         end if
      end do
      alpha = x
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine extent  --  find maximum interatomic distance  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "extent" finds the largest interatomic distance in a system
c
c
      subroutine extent (rmax)
      use atoms
      implicit none
      integer i,k
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 r2,rmax
c
c
c     search all atom pairs to find the largest distance
c
      rmax = 0.0d0
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do k = i+1, n
            xk = x(k)
            yk = y(k)
            zk = z(k)
            r2 = (xi-xk)**2 + (yi-yk)**2 + (zi-zk)**2
            rmax = max(r2,rmax)
         end do
      end do
      rmax = sqrt(rmax)
      return
      end
