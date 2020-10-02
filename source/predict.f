c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2020  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine predict  --  induced dipole prediction values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "predict" checks for use of methods for predicting induced
c     dipoles, extrapolation coefficients and IELSCF parameters
c
c
      subroutine predict
      use atoms
      use ielscf
      use keys
      use uprior
      implicit none
      integer i,j,k
      integer next
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults for use of induced dipole prediction
c
      use_pred = .false.
      use_ielscf = .false.
      polpred = '    '
      maxualt = 6
      nualt = 0
c
c     get keywords containing induced dipole prediction options
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:14) .eq. 'POLAR-PREDICT ') then
            call getword (record,polpred,next)
            call upcase (polpred)
            use_pred = .true.
            if (polpred .eq. '    ') then
               polpred = 'LSQR'
            else if (polpred .eq. 'IEL ') then
               use_pred = .false.
               use_ielscf = .true.
            end if
         else if (keyword(1:8) .eq. 'IEL-SCF ') then
            use_ielscf = .true.
         end if
      end do
c
c     set the 6th-order Gear predictor binomial coefficients
c
      if (polpred .eq. 'GEAR') then
         maxualt = 7
         gear(1) = 6.0d0
         gear(2) = -15.0d0
         gear(3) = 20.0d0
         gear(4) = -15.0d0
         gear(5) = 6.0d0
         gear(6) = -1.0d0
         gear(7) = 0.0d0
      end if
c
c     set 16-step always stable predictor-corrector (ASPC) coefficients
c
      if (polpred .eq. 'ASPC') then
         maxualt = 17
         aspc(1) = 62.0d0 / 17.0d0
         aspc(2) = -310.0d0 / 51.0d0
         aspc(3) = 2170.0d0 / 323.0d0
         aspc(4) = -2329.0d0 / 400.0d0
         aspc(5) = 1701.0d0 / 409.0d0
         aspc(6) = -806.0d0 / 323.0d0
         aspc(7) = 1024.0d0 / 809.0d0
         aspc(8) = -479.0d0 / 883.0d0
         aspc(9) = 257.0d0 / 1316.0d0
         aspc(10) = -434.0d0 / 7429.0d0
         aspc(11) = 191.0d0 / 13375.0d0
         aspc(12) = -62.0d0 / 22287.0d0
         aspc(13) = 3.0d0 / 7217.0d0
         aspc(14) = -3.0d0 / 67015.0d0
         aspc(15) = 2.0d0 / 646323.0d0
         aspc(16) = -1.0d0 / 9694845.0d0
         aspc(17) = 0.0d0
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(udalt))  deallocate (udalt)
      if (allocated(upalt))  deallocate (upalt)
      if (allocated(usalt))  deallocate (usalt)
      if (allocated(upsalt))  deallocate (upsalt)
      if (use_pred) then
         allocate (udalt(maxualt,3,n))
         allocate (upalt(maxualt,3,n))
         allocate (usalt(maxualt,3,n))
         allocate (upsalt(maxualt,3,n))
      end if
c
c     initialize prior values of induced dipole moments
c
      if (use_pred) then
        do i = 1, n
            do j = 1, 3
               do k = 1, maxualt
                  udalt(k,j,i) = 0.0d0
                  upalt(k,j,i) = 0.0d0
                  usalt(k,j,i) = 0.0d0
                  upsalt(k,j,i) = 0.0d0
               end do
            end do
         end do
      end if
c
c     initialize inertial extended Lagrangian method
c
      if (use_ielscf)  call auxinit
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine auxinit  --  initialize auxiliary dipole values  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "auxinit" initializes auxiliary variables and settings for
c     inertial extended Lagrangian induced dipole prediction
c
c     literature reference:
c
c     A. Albaugh, O. Demerdash, and T. Head-Gordon, "An Efficient and
c     Stable Hybrid Extended Lagrangian/Self-Consistent Field Scheme
c     for Solving Classical Mutual Induction", Journal of Chemical
c     Physics, 143, 174104 (2015)
c
c
      subroutine auxinit
      use atomid
      use atoms
      use ielscf
      use keys
      use polar
      implicit none
      integer i,j,next
      real*8 speed
      real*8 weight
      real*8 maxwell
      real*8 vec(3)
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults for auxiliary thermostat control variables
c
      nfree_aux = 3 * npolar
      kelvin_aux = 100000.0d0
      tautemp_aux = 0.1d0
c
c     check for keywords containing auxiliary thermostat values
c 
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:13) .eq. 'AUX-TAUTEMP ') then
            read (string,*,err=10,end=10)  tautemp_aux
         else if (keyword(1:9) .eq. 'AUX-TEMP ') then
            read (string,*,err=10,end=10)  kelvin_aux
         end if
   10    continue
      end do
c
c     perform dynamic allocation of some global arrays
c 
      allocate (uaux(3,n))
      allocate (vaux(3,n))
      allocate (aaux(3,n))
      allocate (upaux(3,n))
      allocate (vpaux(3,n))
      allocate (apaux(3,n))
c
c     set auxiliary dipole values equal to induced dipoles
c
      use_ielscf = .false.
      call induce
      use_ielscf = .true.
      do i = 1, n
         do j = 1, 3
            uaux(j,i) = uind(j,i)
            upaux(j,i) = uinp(j,i)
         end do
      end do
c
c     set velocities and accelerations for auxiliary dipoles
c
      do i = 1, n
         weight = 1.0d0
         speed = maxwell (weight,kelvin_aux)
         call ranvec (vec)
         do j = 1, 3
            vaux(j,i) = speed * vec(j)
            aaux(j,i) = 0.0d0
            vpaux(j,i) = vaux(j,i)
            apaux(j,i) = 0.0d0
         end do
      end do
      return
      end
