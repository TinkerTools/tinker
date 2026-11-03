c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_eabf  --  adaptive biasing force tests  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_eabf" checks the adaptive biasing force bias, the interval
c     sampling built on the ost history, and the abf history restart
c
c
      subroutine test_eabf
      use bath
      implicit none
      logical skiptest
      character*(*) tname
      parameter (tname='test_eabf')
c
c
      if (skiptest(tname,'abf'))  return
      call initial
      kelvin = 300.0d0
      call test_eabf_bias
      call test_eabf_fill
      call test_eabf_fillbias
      call test_eabf_fillfull
      call test_eabf_dyn
      call test_eabf_mean
      call test_eabf_gate
      call test_eabf_resize
      call test_eabf_restart
      call test_eabf_orphan
      call final
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eabf_bias  --  lambda-only bias test  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eabf_bias" checks that eabfbias removes the integrated
c     mean force from the energy, saves its slope for the lambda
c     particle, and leaves the virial untouched
c
c
      subroutine test_eabf_bias
      use dlmda
      use energi
      use mutant
      use virial
      implicit none
      integer i,j
c
c
c     seed a linear mean force f(lambda) = 2 + 4*lambda on the grid
c
      call resetabf (5,4)
      use_abf = .true.
      do i = 1, nlmda
         lmdafmean(i) = 2.0d0 + 4.0d0*dble(i-1)*wlmda
      end do
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = 1.0d0
         end do
      end do
      lambda = 0.3d0
      esum = 1.0d0
      call eabfbias
c
c     the free energy 2*lambda + 2*lambda**2 leaves the energy and
c     its slope 2 + 4*lambda is saved for the lambda particle
c
      call assert_real (esum,1.0d0-0.78d0,1.0d-12,
     &                  'eabfbias removes the free energy')
      call assert_real (lmdavbias,-0.78d0,1.0d-12,
     &                  'eabfbias saves the bias energy')
      call assert_real (lmdadfdl,3.2d0,1.0d-12,
     &                  'eabfbias saves the mean force')
      call assert_real (vir(2,1),1.0d0,0.0d0,
     &                  'eabfbias leaves the virial')
      use_abf = .false.
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_eabf_fill  --  unsampled bin fill rule  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_eabf_fill" checks the fill of unsampled mean force bins:
c     a copy of the nearest sampled bin past either end of the data,
c     linear interpolation inside a gap, zero with no data, sampled
c     bins kept exactly, and the sampled values left unchanged, then
c     the reach: a gap up to twice the reach is interpolated, a wider
c     gap or an end is copied within reach and zero beyond it
c
c
      subroutine test_eabf_fill
      implicit none
      real*8 fsamp(5),fsave(5)
      real*8 fbias(5),ref(5)
      logical sampled(5)
c
c
c     bins below the first sampled bin copy its mean force
c
      fsamp = (/ 0.0d0, 0.0d0, 0.0d0, -20.0d0, -40.0d0 /)
      sampled = (/ .false., .false., .false., .true., .true. /)
      fsave = fsamp
      call lmdafill (5,fsamp,sampled,5,fbias)
      ref = (/ -20.0d0, -20.0d0, -20.0d0, -20.0d0, -40.0d0 /)
      call assert_array1 (fbias,ref,5,0.0d0,
     &                    'lmdafill copies the first sampled bin')
      call assert_array1 (fsamp,fsave,5,0.0d0,
     &                    'lmdafill leaves the sampled mean force')
c
c     bins past the last sampled bin copy its mean force
c
      fsamp = (/ -10.0d0, -20.0d0, 0.0d0, 0.0d0, 0.0d0 /)
      sampled = (/ .true., .true., .false., .false., .false. /)
      call lmdafill (5,fsamp,sampled,5,fbias)
      ref = (/ -10.0d0, -20.0d0, -20.0d0, -20.0d0, -20.0d0 /)
      call assert_array1 (fbias,ref,5,0.0d0,
     &                    'lmdafill copies the last sampled bin')
c
c     a gap between sampled bins is interpolated linearly
c
      fsamp = (/ -10.0d0, 0.0d0, 0.0d0, -40.0d0, -50.0d0 /)
      sampled = (/ .true., .false., .false., .true., .true. /)
      call lmdafill (5,fsamp,sampled,5,fbias)
      ref = (/ -10.0d0, -20.0d0, -30.0d0, -40.0d0, -50.0d0 /)
      call assert_array1 (fbias,ref,5,1.0d-12,
     &                    'lmdafill interpolates a gap')
c
c     a bin halfway between two sampled bins takes their average
c
      fsamp = (/ -10.0d0, 0.0d0, -30.0d0, 0.0d0, 0.0d0 /)
      sampled = (/ .true., .false., .true., .false., .false. /)
      call lmdafill (3,fsamp,sampled,3,fbias)
      call assert_real (fbias(2),-20.0d0,0.0d0,
     &                  'lmdafill averages a halfway bin')
c
c     both ends and a gap are filled in the same profile
c
      fsamp = (/ 0.0d0, -10.0d0, 0.0d0, 0.0d0, -40.0d0 /)
      sampled = (/ .false., .true., .false., .false., .true. /)
      call lmdafill (5,fsamp,sampled,5,fbias)
      ref = (/ -10.0d0, -10.0d0, -20.0d0, -30.0d0, -40.0d0 /)
      call assert_array1 (fbias,ref,5,1.0d-12,
     &                    'lmdafill fills an end and a gap')
c
c     a sampled zero is data, and values in unsampled bins are unused
c
      fsamp = (/ 7.0d0, 9.0d0, 0.0d0, -20.0d0, -40.0d0 /)
      sampled = (/ .false., .false., .true., .true., .true. /)
      call lmdafill (5,fsamp,sampled,5,fbias)
      ref = (/ 0.0d0, 0.0d0, 0.0d0, -20.0d0, -40.0d0 /)
      call assert_array1 (fbias,ref,5,0.0d0,
     &                    'lmdafill keeps a sampled zero')
c
c     with no sampled bin the profile is zero
c
      fsamp = (/ 3.0d0, 1.0d0, 4.0d0, 1.0d0, 5.0d0 /)
      sampled = (/ .false., .false., .false., .false., .false. /)
      call lmdafill (5,fsamp,sampled,5,fbias)
      ref = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
      call assert_array1 (fbias,ref,5,0.0d0,
     &                    'lmdafill gives zero without data')
c
c     a fully sampled profile is reproduced exactly
c
      fsamp = (/ 1.1d0, -2.3d0, 3.7d0, 0.0d0, -5.9d0 /)
      sampled = (/ .true., .true., .true., .true., .true. /)
      call lmdafill (5,fsamp,sampled,5,fbias)
      call assert_array1 (fbias,fsamp,5,0.0d0,
     &                    'lmdafill keeps a fully sampled profile')
c
c     a single bin keeps its data, or is zero without it
c
      fsamp(1) = 2.5d0
      sampled(1) = .true.
      call lmdafill (1,fsamp,sampled,1,fbias)
      call assert_real (fbias(1),2.5d0,0.0d0,
     &                  'lmdafill keeps a single sampled bin')
      sampled(1) = .false.
      call lmdafill (1,fsamp,sampled,1,fbias)
      call assert_real (fbias(1),0.0d0,0.0d0,
     &                  'lmdafill zeroes a single unsampled bin')
c
c     a gap of twice the reach is interpolated across; unsampled
c     bins hold junk values here, which must never reach the fill
c
      fsamp = (/ -10.0d0, 7.0d0, -30.0d0, 9.0d0, 11.0d0 /)
      sampled = (/ .true., .false., .true., .false., .false. /)
      call lmdafill (3,fsamp,sampled,1,fbias)
      ref(1:3) = (/ -10.0d0, -20.0d0, -30.0d0 /)
      call assert_array1 (fbias,ref,3,1.0d-12,
     &                    'lmdafill interpolates a gap of 2 reach 1')
      fsamp = (/ -10.0d0, 7.0d0, 9.0d0, 11.0d0, -50.0d0 /)
      sampled = (/ .true., .false., .false., .false., .true. /)
      call lmdafill (5,fsamp,sampled,2,fbias)
      ref = (/ -10.0d0, -20.0d0, -30.0d0, -40.0d0, -50.0d0 /)
      call assert_array1 (fbias,ref,5,1.0d-12,
     &                    'lmdafill interpolates a gap of 4 reach 2')
c
c     a gap wider than twice the reach copies each side within reach
c
      fsamp = (/ -10.0d0, 7.0d0, 9.0d0, -40.0d0, 11.0d0 /)
      sampled = (/ .true., .false., .false., .true., .false. /)
      call lmdafill (4,fsamp,sampled,1,fbias)
      ref(1:4) = (/ -10.0d0, -10.0d0, -40.0d0, -40.0d0 /)
      call assert_array1 (fbias,ref,4,0.0d0,
     &                    'lmdafill copies both sides of a gap of 3')
      fsamp = (/ -10.0d0, 7.0d0, 9.0d0, 11.0d0, -50.0d0 /)
      sampled = (/ .true., .false., .false., .false., .true. /)
      call lmdafill (5,fsamp,sampled,1,fbias)
      ref = (/ -10.0d0, -10.0d0, 0.0d0, -50.0d0, -50.0d0 /)
      call assert_array1 (fbias,ref,5,0.0d0,
     &                    'lmdafill zeroes the middle of a gap of 4')
c
c     past either end of the data a bin within reach copies it
c
      fsamp = (/ 7.0d0, 9.0d0, 11.0d0, -40.0d0, 13.0d0 /)
      sampled = (/ .false., .false., .false., .true., .false. /)
      call lmdafill (5,fsamp,sampled,1,fbias)
      ref = (/ 0.0d0, 0.0d0, -40.0d0, -40.0d0, -40.0d0 /)
      call assert_array1 (fbias,ref,5,0.0d0,
     &                    'lmdafill copies an end within reach')
c
c     a reach of zero fills nothing
c
      fsamp = (/ 7.0d0, -10.0d0, 9.0d0, -40.0d0, 11.0d0 /)
      sampled = (/ .false., .true., .false., .true., .false. /)
      call lmdafill (5,fsamp,sampled,0,fbias)
      ref = (/ 0.0d0, -10.0d0, 0.0d0, -40.0d0, 0.0d0 /)
      call assert_array1 (fbias,ref,5,0.0d0,
     &                    'lmdafill reach zero fills nothing')
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_eabf_fillbias  --  filled bias vs estimate  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "test_eabf_fillbias" checks that the abf bias and the slope saved
c     for the lambda particle and for Monte Carlo trials follow the
c     filled mean force, while the free energy estimate and the bin
c     mean force keep the sampled data only; a fully sampled profile
c     gives the same bias with or without the fill
c
c     note "resetabf" frees the ost kernels, so the test also shows
c     the abf fill never reads an ost array
c
c
      subroutine test_eabf_fillbias
      use dlmda
      use energi
      use mutant
      implicit none
      integer i
      real*8 vbias,slope
      real*8 vfill,dfill
      real*8 efreetot
c
c
c     sample the two top bins only, at lambda 0.75 and 1, with a
c     reach that fills every bin
c
      call resetabf (5,8)
      lmdafillr = 1.0d0
      use_abf = .true.
      use_abfdyn = .true.
      nlmdahist = 2
      lmdalhist(1) = 0.75d0
      lmdafhist(1) = -20.0d0
      lmdalhist(2) = 1.0d0
      lmdafhist(2) = -40.0d0
      do i = 1, nlmdahist
         call addabfhist (i)
      end do
c
c     the filled bias holds a mean force of -20 below lambda 0.75
c
      use_lmdafill = .true.
      lambda = 1.0d0
      esum = 0.0d0
      call eabfbias
      call assert_real (lmdavbias,22.5d0,1.0d-12,
     &                  'eabfbias fill bias energy at lambda 1')
      call assert_real (esum,22.5d0,1.0d-12,
     &                  'eabfbias fill energy at lambda 1')
      call assert_real (lmdadfdl,-40.0d0,0.0d0,
     &                  'eabfbias fill slope at lambda 1')
      lambda = 0.3d0
      esum = 0.0d0
      call eabfbias
      call assert_real (lmdavbias,6.0d0,1.0d-12,
     &                  'eabfbias fill bias energy in a gap')
      call assert_real (lmdadfdl,-20.0d0,0.0d0,
     &                  'eabfbias fill slope in a gap')
c
c     a Monte Carlo trial sees the same filled bias and slope
c
      call elmdabias (vbias,slope)
      call assert_real (vbias,lmdavbias,0.0d0,
     &                  'elmdabias matches the filled abf bias')
      call assert_real (slope,lmdadfdl,0.0d0,
     &                  'elmdabias matches the filled abf slope')
c
c     the estimate keeps zero mean force in the unsampled bins
c
      call assert_real (efreetot(),-10.0d0,1.0d-12,
     &                  'efreetot uses sampled bins only')
      call assert_real (lmdafmean(1),0.0d0,0.0d0,
     &                  'fill leaves the empty bin mean force')
      call assert_real (lmdafmean(4),-20.0d0,0.0d0,
     &                  'fill leaves the sampled bin mean force')
c
c     a reach of one bin fills only the bin below lambda 0.75, so
c     the mean force is 0, 0, -20, -20, -40
c
      lmdafillr = 0.25d0
      lambda = 0.3d0
      call eabfbias
      call assert_real (lmdavbias,0.1d0,1.0d-12,
     &                  'eabfbias reach bias energy past the ramp')
      call assert_real (lmdadfdl,-4.0d0,1.0d-12,
     &                  'eabfbias reach slope past the ramp')
      lambda = 0.2d0
      call eabfbias
      call assert_real (lmdavbias,0.0d0,0.0d0,
     &                  'eabfbias reach bias energy beyond reach')
      call assert_real (lmdadfdl,0.0d0,0.0d0,
     &                  'eabfbias reach slope beyond reach')
      lmdafillr = 1.0d0
c
c     without the fill the bias is the sampled free energy
c
      use_lmdafill = .false.
      lambda = 0.3d0
      call eabfbias
      call assert_real (lmdavbias,0.0d0,0.0d0,
     &                  'eabfbias nofill bias energy in a gap')
      call assert_real (lmdadfdl,0.0d0,0.0d0,
     &                  'eabfbias nofill slope in a gap')
      lambda = 0.6d0
      call eabfbias
      call assert_real (lmdavbias,0.4d0,1.0d-12,
     &                  'eabfbias nofill bias energy near data')
      call assert_real (lmdadfdl,-8.0d0,1.0d-12,
     &                  'eabfbias nofill slope near data')
      call elmdabias (vbias,slope)
      call assert_real (vbias,lmdavbias,0.0d0,
     &                  'elmdabias matches the nofill abf bias')
      call assert_real (slope,lmdadfdl,0.0d0,
     &                  'elmdabias matches the nofill abf slope')
c
c     once every bin is sampled the fill changes nothing
c
      nlmdahist = 5
      lmdalhist(3) = 0.0d0
      lmdafhist(3) = 2.0d0
      lmdalhist(4) = 0.25d0
      lmdafhist(4) = 3.0d0
      lmdalhist(5) = 0.5d0
      lmdafhist(5) = -5.0d0
      do i = 3, nlmdahist
         call addabfhist (i)
      end do
      lambda = 0.6d0
      use_lmdafill = .true.
      call eabfbias
      vfill = lmdavbias
      dfill = lmdadfdl
      use_lmdafill = .false.
      call eabfbias
      call assert_real (vfill,lmdavbias,0.0d0,
     &                  'fully sampled fill bias energy unchanged')
      call assert_real (dfill,lmdadfdl,0.0d0,
     &                  'fully sampled fill slope unchanged')
      lmdafillr = 0.1d0
      use_abf = .false.
      use_abfdyn = .false.
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eabf_fillfull  --  unit reach fills all  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eabf_fillfull" checks on the default grid of 201 bins that
c     a fill reach of one lambda, or more, gives the unlimited fill,
c     where each end copies the nearest sampled bin and a gap of any
c     width is interpolated, and that a reach of 0.25 copies each side
c     of a wider gap and the low end only within reach; each case is
c     compared with a profile built here by hand
c
c
      subroutine test_eabf_fillfull
      use dlmda
      implicit none
      integer i,j,k
      integer ilo,ihi,nr
      real*8 flo,fhi
      real*8 lmda
      real*8 eb,db,er,dr
      real*8 reach(3)
      real*8, allocatable :: fref(:)
      character*8 rtag
c
c
c     sample only lambda 0.3 and 0.9, leaving 60 bins below the data,
c     a gap of 120 bins and 20 bins above the data
c
      call resetabf (201,1)
      use_abf = .true.
      use_lmdafill = .true.
      ilo = 61
      ihi = 181
      flo = -20.0d0
      fhi = -80.0d0
      lmdafmean(ilo) = flo
      lmdafwt(ilo) = 1.0d0
      lmdafmean(ihi) = fhi
      lmdafwt(ihi) = 1.0d0
      allocate (fref(nlmda))
c
c     a reach of one or more matches the unlimited fill; a reach of
c     0.25 is 50 bins, so the gap is wider than twice the reach
c
      reach = (/ 1.0d0, 1.5d0, 0.25d0 /)
      do j = 1, 3
         lmdafillr = reach(j)
         write (rtag,'(f4.2)')  reach(j)
c
c     build the expected profile by hand, where a reach of one or
c     more copies each end and interpolates the gap, and the reach
c     of 0.25 is "nr" bins, which is less than half of the gap
c
         nr = 50
         do k = 1, nlmda
            if (reach(j) .ge. 1.0d0) then
               if (k .le. ilo) then
                  fref(k) = flo
               else if (k .ge. ihi) then
                  fref(k) = fhi
               else
                  fref(k) = flo + dble(k-ilo)/dble(ihi-ilo)*(fhi-flo)
               end if
            else
               fref(k) = 0.0d0
               if (k .le. ilo) then
                  if (ilo-k .le. nr)  fref(k) = flo
               else if (k .lt. ihi) then
                  if (k-ilo .le. nr)  fref(k) = flo
                  if (ihi-k .le. nr)  fref(k) = fhi
               else
                  if (k-ihi .le. nr)  fref(k) = fhi
               end if
            end if
         end do
c
c     compare the bias energy and slope over the lambda range
c
         do i = 0, 80
            lmda = dble(i) * 0.0125d0
            call efreebias (lmda,eb,db)
            call efreeprof (lmda,fref,er,dr)
            call assert_real (eb,er,1.0d-10,
     &         'efreebias reach '//trim(rtag)//' full grid energy')
            call assert_real (db,dr,1.0d-10,
     &         'efreebias reach '//trim(rtag)//' full grid slope')
         end do
      end do
c
c     check a few bins of the 0.25 reach directly, at the edge of the
c     low end and of each side of the gap
c
      call efreebias (0.045d0,eb,db)
      call assert_real (db,0.0d0,1.0d-9,
     &                  'efreebias reach 0.25 beyond low end')
      call efreebias (0.05d0,eb,db)
      call assert_real (db,flo,1.0d-9,
     &                  'efreebias reach 0.25 low end edge')
      call efreebias (0.55d0,eb,db)
      call assert_real (db,flo,1.0d-9,
     &                  'efreebias reach 0.25 gap lower edge')
      call efreebias (0.6d0,eb,db)
      call assert_real (db,0.0d0,1.0d-9,
     &                  'efreebias reach 0.25 gap middle')
      call efreebias (0.65d0,eb,db)
      call assert_real (db,fhi,1.0d-9,
     &                  'efreebias reach 0.25 gap upper edge')
      lmdafillr = 0.1d0
      deallocate (fref)
      use_abf = .false.
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_eabf_dyn  --  interval sample tests  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_eabf_dyn" drives elmdadyn over four sample intervals and
c     checks the sample, bin accumulators and free energy of a settled
c     interval, that a noisy interval is kept while the convergence
c     gate is off, and that the gate rejects it once turned on
c
c
      subroutine test_eabf_dyn
      use dlmda
      use mutant
      implicit none
      integer istep
      real*8 egsave
c
c
c     a settled interval records its lambda and dU/dlambda average
c
      call resetabf (5,4)
      use_abf = .true.
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      lmdadt = 0.0d0
      lmdadfdl = 0.5d0
      lmdastep = 0
      do istep = 1, lmdaintv
         lambda = 0.5d0
         dedl = 3.0d0
         call elmdadyn
         if (istep .lt. lmdaintv) then
            call assert_int (nlmdahist,0,
     &                       'elmdadyn abf waits for the interval end')
         end if
      end do
      call assert_real (deffdl,2.5d0,1.0d-12,
     &                  'elmdadyn abf removes the mean force')
      call assert_int (nlmdahist,1,'elmdadyn abf records one interval')
      call assert_int (lmdaihist(1),lmdaintv,
     &                 'elmdadyn abf stamps the step')
      call assert_real (lmdalhist(1),0.5d0,0.0d0,
     &                  'elmdadyn abf sample lambda')
      call assert_real (lmdafhist(1),3.0d0,1.0d-12,
     &                  'elmdadyn abf sample dU/dlambda')
      call assert_real (lmdafwt(3),1.0d0,0.0d0,'elmdadyn abf bin count')
      call assert_real (lmdafmean(3),3.0d0,1.0d-12,
     &                  'elmdadyn abf mean force')
      call assert_real (lmdadeltag,0.75d0,1.0d-12,
     &                  'elmdadyn abf free energy estimate')
c
c     a noisy interval is still recorded while the gate is off
c
      do istep = 1, lmdaintv
         lambda = 0.5d0
         dedl = 1.0d0
         if (mod(istep,2) .eq. 0)  dedl = 11.0d0
         call elmdadyn
      end do
      call assert_int (nlmdahist,2,
     &                 'elmdadyn abf keeps a noisy interval')
      call assert_real (lmdafwt(3),2.0d0,0.0d0,
     &                  'elmdadyn abf noisy bin count')
      call assert_real (lmdafmean(3),4.5d0,1.0d-12,
     &                  'elmdadyn abf noisy mean force')
c
c     with the gate on a noisy interval leaves the history unchanged
c
      use_lmdacv = .true.
      lmdacvstd = 1.0d0
      lmdacvrat = 0.0d0
      egsave = lmdadeltag
      do istep = 1, lmdaintv
         lambda = 0.5d0
         dedl = 1.0d0
         if (mod(istep,2) .eq. 0)  dedl = 11.0d0
         call elmdadyn
      end do
      call assert_int (nlmdahist,2,
     &                 'elmdadyn abf gate rejects a noisy one')
      call assert_real (lmdafwt(3),2.0d0,0.0d0,
     &                  'elmdadyn abf gated bin count')
      call assert_real (lmdafmean(3),4.5d0,1.0d-12,
     &                  'elmdadyn abf gated mean force')
      call assert_real (lmdadeltag,egsave,1.0d-12,
     &                  'elmdadyn abf gated free energy')
c
c     with the gate on a settled interval is still recorded
c
      do istep = 1, lmdaintv
         lambda = 0.5d0
         dedl = 3.0d0
         call elmdadyn
      end do
      call assert_int (nlmdahist,3,
     &                 'elmdadyn abf gate keeps a settled one')
      call assert_real (lmdafwt(3),3.0d0,0.0d0,
     &                  'elmdadyn abf gate settled bin count')
      call assert_real (lmdafmean(3),4.0d0,1.0d-12,
     &                  'elmdadyn abf gate settled mean force')
      use_lmdacv = .false.
      use_abf = .false.
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_eabf_mean  --  bin mean force tests  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_eabf_mean" checks that each sample sets the bin bias force
c     to the running mean of its bin, and that a rebuild from the
c     history reproduces the accumulated bins
c
c
      subroutine test_eabf_mean
      use dlmda
      implicit none
      integer i
      real*8 ref
c
c
c     the bias force follows the running mean of the bin samples
c
      call resetabf (5,8)
      do i = 1, 6
         lmdalhist(i) = 0.25d0
         lmdafhist(i) = dble(i)
         nlmdahist = i
         call addabfhist (i)
         ref = 0.5d0 * dble(i+1)
         call assert_real (lmdafmean(2),ref,1.0d-12,
     &                     'addabfhist sets the running mean force')
      end do
c
c     a rebuild from the history reproduces the accumulated bins
c
      call buildabfkernel
      call assert_real (lmdafwt(2),6.0d0,0.0d0,
     &                  'buildabfkernel bin count')
      call assert_real (lmdafsum(2),21.0d0,1.0d-12,
     &                  'buildabfkernel bin sum')
      call assert_real (lmdafmean(2),3.5d0,1.0d-12,
     &                  'buildabfkernel bin mean force')
      call assert_real (lmdafmean(1),0.0d0,0.0d0,
     &                  'buildabfkernel empty bin')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eabf_gate  --  lambda freezing tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eabf_gate" drives elmdadyn over one sample interval and
c     checks that the lambda particle moves only during the leading
c     propagation phase and that the sample sits on the frozen lambda
c
c
      subroutine test_eabf_gate
      use dlmda
      use mutant
      implicit none
      integer istep
      real*8 lam(6)
      real*8 frozen
c
c
c     drive an interval with a deterministic frictionless lambda
c     particle, so that any lambda motion comes from the gate alone
c
      call resetabf (5,4)
      use_abf = .true.
      lmdaintv = 6
      lmdanpa = 2
      lmdanpb = 2
      lmdanpc = 2
      lmdadt = 0.1d0
      lmdamass = 1.0d0
      lmdafric = 0.0d0
      lmdatheta = 0.25d0 * 3.14159265358979323846d0
      lmdavtheta = 0.0d0
      lambda = 0.5d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      do istep = 1, lmdaintv
         dedl = 1.0d0
         call elmdadyn
         lam(istep) = lambda
      end do
c
c     the particle moves only while the interval is in phase a
c
      call assert_logical (lam(1).ne.lam(2),.true.,
     &                     'elmdadyn abf propagates during phase a')
      frozen = lam(2)
      do istep = 3, lmdaintv
         call assert_real (lam(istep),frozen,0.0d0,
     &                     'elmdadyn abf holds lambda after phase a')
      end do
c
c     the sample is recorded exactly on the frozen lambda
c
      call assert_int (nlmdahist,1,
     &                 'elmdadyn abf samples a frozen interval')
      call assert_real (lmdalhist(1),frozen,0.0d0,
     &                  'elmdadyn abf sample on the frozen lambda')
      call assert_real (lmdafhist(1),1.0d0,1.0d-12,
     &                  'elmdadyn abf sample dU/dlambda on the gate')
      use_abf = .false.
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_eabf_resize  --  sample history growth  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_eabf_resize" records more interval samples than the abf
c     history holds and checks that the history doubles while every
c     sample and the bin mean force are kept
c
c
      subroutine test_eabf_resize
      use dlmda
      use mutant
      implicit none
      integer istep
c
c
c     three samples overflow a history sized for two
c
      call resetabf (5,2)
      use_abf = .true.
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      lmdadt = 0.0d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      do istep = 1, 3*lmdaintv
         lambda = 0.25d0
         dedl = dble((istep-1)/lmdaintv + 1)
         call elmdadyn
      end do
      call assert_int (nlmdahist,3,'resizeabfhist sample count')
      call assert_int (sizelmdahist,4,'resizeabfhist doubles the size')
      call assert_int (lmdaihist(1),4,'resizeabfhist keeps step 1')
      call assert_int (lmdaihist(3),12,'resizeabfhist keeps step 3')
      call assert_real (lmdafhist(1),1.0d0,1.0d-12,
     &                  'resizeabfhist keeps sample 1')
      call assert_real (lmdafhist(2),2.0d0,1.0d-12,
     &                  'resizeabfhist keeps sample 2')
      call assert_real (lmdafhist(3),3.0d0,1.0d-12,
     &                  'resizeabfhist keeps sample 3')
      call assert_real (lmdalhist(3),0.25d0,0.0d0,
     &                  'resizeabfhist keeps the sample lambda')
      call assert_real (lmdafmean(2),2.0d0,1.0d-12,
     &                  'resizeabfhist bin mean force')
      use_abf = .false.
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_eabf_restart  --  history restart test  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_eabf_restart" writes an abf history, checks that repeated
c     saves append only new samples, then continues the history in a
c     fresh run and checks the restored bins, lambda particle and
c     interval boundary, with the run settings winning over the file
c
c
      subroutine test_eabf_restart
      use bath
      use dlmda
      use files
      use mutant
      use ost
      implicit none
      integer i,ihis
      integer iline
      integer istep
      integer leng0
      integer size0,size1,size2,size3
      integer nbyte0,nbyte1
      integer freeunit
      integer trimtext
      real*8 kelvin0
      real*8 eostref
      real*8 fref(5)
      real*8 sref(5)
      real*8 pref(5)
      logical exist
      logical hastitle
      logical haslmda
      logical haslabel
      character*240 filename0
      character*240 abffile
      character*240 record
      character*8000 text0,text1
c
c
c     start an abf history in a temporary file
c
      filename0 = filename
      leng0 = leng
      kelvin0 = kelvin
      kelvin = 300.0d0
      filename = 'tinkertest-saveabf'
      leng = len_trim(filename)
      abffile = filename(1:leng)//'.abf'
      inquire (file=abffile,exist=exist)
      if (exist) then
         ihis = freeunit ()
         open (unit=ihis,file=abffile,status='old')
         close (unit=ihis,status='delete')
      end if
      call resetabf (5,4)
      use_ost = .false.
      use_abf = .true.
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      lmdadt = 0.0d0
      lmdamass = 2.0d0
      lmdafric = 0.5d0
      lmdatheta = 0.3d0
      lmdavtheta = 0.1d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      call initabffile
      inquire (file=abffile,exist=exist)
      call assert_logical (exist,.true.,'initabffile creates the file')
c
c     the history title marks the file as an abf history
c
      ihis = freeunit ()
      open (unit=ihis,file=abffile,status='old')
      read (ihis,10)  record
   10 format (a240)
      hastitle = (index(record,abftitle(1:trimtext(abftitle))) .gt. 0)
      do iline = 2, 4
         read (ihis,10)  record
      end do
      haslmda = (index(record,'Lambda State') .gt. 0)
      do iline = 5, 6
         read (ihis,10)  record
      end do
      close (unit=ihis)
      haslabel = (index(record,abflabel(1:trimtext(abflabel))) .gt. 0)
      call assert_logical (hastitle,.true.,
     &                     'initabffile writes the abf title')
      call assert_logical (haslmda,.true.,
     &                     'initabffile writes no ost grid state')
      call assert_logical (haslabel,.true.,
     &                     'initabffile writes the sample label')
c
c     save before any interval completes; the header-only file is
c     shorter than one header record buffer
c
      inquire (file=abffile,size=size0)
      call assert_logical (size0.lt.240,.true.,
     &                     'initabffile header is shorter than 240')
      call saveabf
      inquire (file=abffile,size=size1)
      call assert_int (size1,size0,'saveabf updates a header-only file')
c
c     record two intervals, then save twice with no new samples
c
      do istep = 1, 2*lmdaintv
         lambda = 0.25d0
         dedl = 1.0d0
         if (istep .gt. lmdaintv) then
            lambda = 0.75d0
            dedl = 5.0d0
         end if
         call elmdadyn
      end do
      call saveabf
      inquire (file=abffile,size=size1)
      call saveabf
      inquire (file=abffile,size=size2)
      call assert_int (size2,size1,'saveabf does not duplicate')
c
c     a third interval and part of a fourth, then save again
c
      do istep = 1, lmdaintv+2
         lambda = 0.75d0
         dedl = 7.0d0
         call elmdadyn
      end do
      call saveabf
      inquire (file=abffile,size=size3)
      call assert_logical (size3.gt.size2,.true.,
     &                     'saveabf appends new samples')
      call assert_int (nlmdahist,3,'elmdadyn abf restart sample count')
      call assert_int (lmdastep,14,'elmdadyn abf restart partial step')
      do i = 1, 5
         fref(i) = lmdafmean(i)
         sref(i) = lmdafsum(i)
         pref(i) = lmdafwt(i)
      end do
      eostref = lmdadeltag
c
c     a fresh run with matching grid and interval continues the file
c
      call resetabf (5,4)
      use_abf = .true.
      lmdaintv = 4
      kelvin = 310.0d0
      lmdamass = 9.0d0
      lmdafric = 0.2d0
      lmdadt = 0.002d0
      lmdatheta = 0.0d0
      lmdavtheta = 0.0d0
      call eostfilebytes (abffile,text0,nbyte0)
      call initabffile
      call eostfilebytes (abffile,text1,nbyte1)
      call assert_logical (histtail,.false.,
     &                     'rdabfhist finds no rows past the count')
      call assert_logical (nbyte1.eq.nbyte0 .and. text1.eq.text0,
     &                     .true.,'initabffile leaves a clean file')
      call assert_logical (lmdasavefile.eq.abffile,.true.,
     &                     'initabffile continues the same file')
      call assert_int (nlmdahist,3,'initabffile restores the samples')
      call assert_int (nlmdasave,3,'initabffile restores the save')
      call assert_int (lmdastep,12,'initabffile resumes on an interval')
      call assert_array1 (lmdafmean,fref,5,1.0d-12,
     &                    'initabffile restores the mean force')
      call assert_array1 (lmdafsum,sref,5,1.0d-12,
     &                    'initabffile restores the bin sums')
      call assert_array1 (lmdafwt,pref,5,0.0d0,
     &                    'initabffile restores the bin counts')
      call assert_real (lmdadeltag,eostref,1.0d-12,
     &                  'initabffile restores the free energy')
      call assert_real (lmdatheta,0.3d0,1.0d-14,
     &                  'initabffile restores the lambda theta')
      call assert_real (lmdavtheta,0.1d0,1.0d-14,
     &                  'initabffile restores the theta velocity')
      call assert_real (kelvin,310.0d0,0.0d0,
     &                  'initabffile keeps the run temperature')
      call assert_real (lmdamass,9.0d0,0.0d0,
     &                  'initabffile keeps the run lambda mass')
      call assert_real (lmdadt,0.002d0,0.0d0,
     &                  'initabffile keeps the run lambda step')
      call assert_logical (allocated(osthist),.false.,
     &                     'rdabf leaves the ost history unallocated')
      call assert_logical (allocated(gkernel),.false.,
     &                     'rdabf leaves the g kernel unallocated')
c
c     remove the temporary file and restore the global state
c
      ihis = freeunit ()
      open (unit=ihis,file=abffile,status='old')
      close (unit=ihis,status='delete')
      filename = filename0
      leng = leng0
      kelvin = kelvin0
      use_abf = .false.
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eabf_orphan  --  stray row restart test  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eabf_orphan" appends a sample row past the header count,
c     as left by a run stopped between appending samples and updating
c     the header, and checks that continuing the history drops it
c
c
      subroutine test_eabf_orphan
      use bath
      use dlmda
      use files
      use mutant
      implicit none
      integer ihis
      integer istep
      integer leng0
      integer freeunit
      real*8 kelvin0
      logical exist
      character*240 filename0
      character*240 abffile
c
c
c     start an abf history with two intervals in a temporary file
c
      filename0 = filename
      leng0 = leng
      kelvin0 = kelvin
      kelvin = 300.0d0
      filename = 'tinkertest-orphanabf'
      leng = len_trim(filename)
      abffile = filename(1:leng)//'.abf'
      inquire (file=abffile,exist=exist)
      if (exist) then
         ihis = freeunit ()
         open (unit=ihis,file=abffile,status='old')
         close (unit=ihis,status='delete')
      end if
      call resetabf (5,4)
      use_ost = .false.
      use_abf = .true.
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      lmdadt = 0.0d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      call initabffile
      do istep = 1, 2*lmdaintv
         lambda = 0.25d0
         dedl = 1.0d0
         call elmdadyn
      end do
      call saveabf
c
c     append a stray row past the header count
c
      ihis = freeunit ()
      open (unit=ihis,file=abffile,status='old',position='append')
      write (ihis,10)  999,0.9d0,1.0d6
   10 format (i12,2d26.16)
      close (unit=ihis)
c
c     continuing the history drops the stray row before new samples
c
      call resetabf (5,4)
      use_abf = .true.
      lmdaintv = 4
      kelvin = 300.0d0
      call initabffile
      call assert_logical (histtail,.true.,
     &                     'rdabfhist finds the stray row')
      call assert_int (nlmdahist,2,'initabffile skips the stray row')
      do istep = 1, lmdaintv
         lambda = 0.75d0
         dedl = 5.0d0
         call elmdadyn
      end do
      call saveabf
c
c     a fresh read holds the kept and new samples and no stray bin
c
      call resetabf (5,4)
      use_abf = .true.
      lmdasavefile = abffile
      call rdabf
      call assert_int (nlmdahist,3,'rdabf reads the kept samples')
      call assert_int (lmdaihist(3),12,'rdabf reads the new step')
      call assert_real (lmdalhist(3),0.75d0,0.0d0,
     &                  'rdabf reads the new sample lambda')
      call assert_real (lmdafhist(3),5.0d0,1.0d-12,
     &                  'rdabf reads the new sample average')
      call assert_real (lmdafwt(5),0.0d0,0.0d0,
     &                  'rdabf builds no bin from the stray row')
c
c     remove the temporary file and restore the global state
c
      ihis = freeunit ()
      open (unit=ihis,file=abffile,status='old')
      close (unit=ihis,status='delete')
      filename = filename0
      leng = leng0
      kelvin = kelvin0
      use_abf = .false.
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine resetabf  --  reset ABF test state  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "resetabf" sets up the ost test state and then frees every array
c     that the adaptive biasing force method must never touch, so any
c     stray access from an abf routine fails the test
c
c
      subroutine resetabf (nl,nhist)
      use ost
      implicit none
      integer nl,nhist
c
c
c     start from the full deterministic ost test state
c
      call resetost (nl,5,nhist)
c
c     free the arrays used only by ost and metadynamics
c
      deallocate (osthist)
      deallocate (ostnext)
      deallocate (osthead)
      deallocate (osthhist)
      deallocate (ostwlhist)
      deallocate (ostwfhist)
      deallocate (gfkernel)
      deallocate (gkernel)
      deallocate (glfkernel)
      deallocate (glkernel)
      deallocate (vkernelmax)
      return
      end
