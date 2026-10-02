c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_barostat  --  Monte Carlo barostat dipoles  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "test_barostat" checks that a rejected Monte Carlo volume move
c     leaves the induced dipoles and their predictor history as they
c     were, while an accepted move keeps the trial dipoles
c
c
      subroutine test_barostat
      use atoms
      use bath
      use boxes
      use dlmda
      use mdstuf
      use polar
      use uprior
      implicit none
      integer i,step
      integer nualt0
      real*8 energy,e0
      real*8 epot,temp
      real*8 xbox0,box0(7)
      real*8, allocatable :: x0(:)
      real*8, allocatable :: y0(:)
      real*8, allocatable :: z0(:)
      real*8, allocatable :: uind0(:,:)
      real*8, allocatable :: uinp0(:,:)
      real*8, allocatable :: udir0(:,:)
      real*8, allocatable :: udirp0(:,:)
      real*8, allocatable :: udalt0(:,:,:)
      real*8, allocatable :: upalt0(:,:,:)
      logical skiptest
      character*(*) tname
      parameter (tname='test_barostat')
c
c
      if (skiptest(tname,'amoeba'))  return
      call pushdir ('file/barostat')
      call loadfix ('water30','water30.key')
      call predict
      call assert_int (n,2700,tname//' atoms')
c
c     set the Monte Carlo barostat state that "mdinit" would read
c
      kelvin = 298.0d0
      temp = kelvin
      isothermal = .true.
      atmsph = 1.0d0
      voltrial = 1
      volmove = 10.0d0
      volscale = 'MOLECULAR'
      prestyp = 'ISOTROPIC'
      integrate = 'VERLET'
      use_ostdyn = .false.
      use_metadyn = .false.
      use_abfdyn = .false.
c
c     match the C++ fixture's 19 solves and zero-based perturbations;
c     ASPC stores 17 slots here, including a zero-coefficient slot,
c     and 16 in C++, so both histories are full before the trial
c
      do step = 0, 18
         do i = 1, n
            x(i) = x(i) + 0.01d0*dble(mod(i-1+step,3)-1)
         end do
         e0 = energy ()
      end do
      call assert_int (nualt,maxualt,tname//' history full')
      allocate (x0(n))
      allocate (y0(n))
      allocate (z0(n))
      allocate (uind0(3,n))
      allocate (uinp0(3,n))
      allocate (udir0(3,n))
      allocate (udirp0(3,n))
      allocate (udalt0(maxualt,3,n))
      allocate (upalt0(maxualt,3,n))
c
c     use the same forced-rejection energy as the C++ test
c
      box0 = (/ xbox,ybox,zbox,alpha,beta,gamma,volbox /)
      x0 = x(1:n)
      y0 = y(1:n)
      z0 = z(1:n)
      uind0 = uind
      uinp0 = uinp
      udir0 = udir
      udirp0 = udirp
      nualt0 = nualt
      udalt0 = udalt
      upalt0 = upalt
      epot = -1.0d30
      call pmonte (epot,temp)
      call assert_real (epot,-1.0d30,0.0d0,tname//' reject epot')
      call assert_real (xbox,box0(1),0.0d0,tname//' reject xbox')
      call assert_real (ybox,box0(2),0.0d0,tname//' reject ybox')
      call assert_real (zbox,box0(3),0.0d0,tname//' reject zbox')
      call assert_real (alpha,box0(4),0.0d0,tname//' reject alpha')
      call assert_real (beta,box0(5),0.0d0,tname//' reject beta')
      call assert_real (gamma,box0(6),0.0d0,tname//' reject gamma')
      call assert_real (volbox,box0(7),0.0d0,tname//' reject volume')
      call assert_array1 (x,x0,n,0.0d0,tname//' reject x')
      call assert_array1 (y,y0,n,0.0d0,tname//' reject y')
      call assert_array1 (z,z0,n,0.0d0,tname//' reject z')
      call assert_array2 (uind,uind0,3,n,0.0d0,tname//' reject uind')
      call assert_array2 (uinp,uinp0,3,n,0.0d0,tname//' reject uinp')
      call assert_array2 (udir,udir0,3,n,0.0d0,tname//' reject udir')
      call assert_array2 (udirp,udirp0,3,n,0.0d0,
     &                    tname//' reject udirp')
      call assert_int (nualt,nualt0,tname//' reject nualt')
      call assert_array1 (udalt,udalt0,3*n*maxualt,0.0d0,
     &                    tname//' reject udalt')
      call assert_array1 (upalt,upalt0,3*n*maxualt,0.0d0,
     &                    tname//' reject upalt')
c
c     an old energy far above the trial energy accepts the move,
c     which adds the trial dipoles to the predictor history
c
      xbox0 = xbox
      x0 = x(1:n)
      y0 = y(1:n)
      z0 = z(1:n)
      uind0 = uind
      uinp0 = uinp
      udir0 = udir
      udirp0 = udirp
      nualt0 = nualt
      udalt0 = udalt
      upalt0 = upalt
      epot = e0 + 200.0d0
      call pmonte (epot,temp)
      call assert_logical (xbox.ne.xbox0,.true.,tname//' accept box')
      call assert_int (nualt,maxualt,tname//' accept nualt')
      call assert_array1 (udalt(1,:,:),uind,3*n,0.0d0,
     &                    tname//' accept udalt newest')
      call assert_array1 (upalt(1,:,:),uinp,3*n,0.0d0,
     &                    tname//' accept upalt newest')
      call assert_array1 (udalt(2:maxualt,:,:),
     &                    udalt0(1:maxualt-1,:,:),3*n*(maxualt-1),
     &                    0.0d0,tname//' accept udalt shifted')
      call assert_array1 (upalt(2:maxualt,:,:),
     &                    upalt0(1:maxualt-1,:,:),3*n*(maxualt-1),
     &                    0.0d0,tname//' accept upalt shifted')
      deallocate (x0)
      deallocate (y0)
      deallocate (z0)
      deallocate (uind0)
      deallocate (uinp0)
      deallocate (udir0)
      deallocate (udirp0)
      deallocate (udalt0)
      deallocate (upalt0)
      call popdir
      call final
      return
      end
