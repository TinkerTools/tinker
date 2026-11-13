c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_eost  --  OST kernel and history tests  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_eost" checks OST helper routines against
c     deterministic inputs and expected outputs
c
c
      subroutine test_eost
      use bath
      use units
      implicit none
      logical skiptest,skipgroup
      character*(*) tname
      parameter (tname='test_eost')
c
c
      if (skipgroup(tname,'ost'))  return
      call initial
      kelvin = 300.0d0
      if (.not. skiptest(tname//'_index','ost',2))
     &   call test_eost_index
      if (.not. skiptest(tname//'_resize','ost',2))
     &   call test_eost_resize
      if (.not. skiptest(tname//'_buildindex','ost',2))
     &   call test_eost_buildindex
      if (.not. skiptest(tname//'_ensure','ost',2))
     &   call test_eost_ensure
      if (.not. skiptest(tname//'_gkernels','ost',2))
     &   call test_eost_gkernels
      if (.not. skiptest(tname//'_fkernel','ost',2))
     &   call test_eost_fkernel
      if (.not. skiptest(tname//'_kernelbuilds','ost',2))
     &   call test_eost_kernelbuilds
      if (.not. skiptest(tname//'_eginterpolate','ost',2))
     &   call test_eost_eginterpolate
      if (.not. skiptest(tname//'_avgstd','ost',2))
     &   call test_eost_avgstd
      if (.not. skiptest(tname//'_depcriteria','ost',2))
     &   call test_eost_depcriteria
      if (.not. skiptest(tname//'_efkernel','ost',2))
     &   call test_eost_efkernel
      if (.not. skiptest(tname//'_vkernelmax','ost',2))
     &   call test_eost_vkernelmax
      if (.not. skiptest(tname//'_fill','ost',2))
     &   call test_eost_fill
      if (.not. skiptest(tname//'_tempering','ost',2))
     &   call test_eost_tempering
      if (.not. skiptest(tname//'_ostdyn','ost',2))
     &   call test_eost_ostdyn
      if (.not. skiptest(tname//'_ostlocal','ost',2))
     &   call test_eost_ostlocal
      if (.not. skiptest(tname//'_ostphase','ost',2))
     &   call test_eost_ostphase
      if (.not. skiptest(tname//'_ostgate','ost',2))
     &   call test_eost_ostgate
      if (.not. skiptest(tname//'_save','ost',2))
     &   call test_eost_save
      if (.not. skiptest(tname//'_restart','ost',2))
     &   call test_eost_restart
      if (.not. skiptest(tname//'_meta','ost',2))
     &   call test_eost_meta
      if (.not. skiptest(tname//'_metaimage','ost',2))
     &   call test_eost_metaimage
      if (.not. skiptest(tname//'_metadyn','ost',2))
     &   call test_eost_metadyn
      if (.not. skiptest(tname//'_metatemper','ost',2))
     &   call test_eost_metatemper
      if (.not. skiptest(tname//'_metacut','ost',2))
     &   call test_eost_metacut
      if (.not. skiptest(tname//'_metarestart','ost',2))
     &   call test_eost_metarestart
      if (.not. skiptest(tname//'_temperkeys','ost',2))
     &   call test_eost_temperkeys
      if (.not. skiptest(tname//'_restartrun','ost',4))
     &   call test_eost_restartrun
      call final
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_eost_save  --  restart append test  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_eost_save" checks that repeated saves append only new
c     histories while updating the fixed-size restart header in place
c
c
      subroutine test_eost_save
      use bath
      use dlmda
      use files
      use ost
      implicit none
      integer i,ihis
      integer ios
      integer leng0
      integer nline
      integer nsave0
      integer size1,size2,size3
      integer freeunit
      real*8 brutevkmax
      real*8 kelvin0
      logical exist
      character*240 filename0
      character*240 ostfile
      character*240 record
c
c     create a small deterministic history and a temporary restart
c
      filename0 = filename
      leng0 = leng
      kelvin0 = kelvin
      kelvin = 321.0d0
      filename = 'tinkertest-saveost'
      leng = len_trim(filename)
      ostfile = filename(1:leng)//'.ost'
      inquire (file=ostfile,exist=exist)
      if (exist) then
         ihis = freeunit ()
         open (unit=ihis,file=ostfile,status='old')
         close (unit=ihis,status='delete')
      end if
      call resetost (3,5,3)
      use_ost = .true.
      nlmdasave = 0
      lmdastep = 10
      lmdatheta = 0.0d0
      lmdavtheta = 0.0d0
      lmdamass = 1.0d0
      lmdafric = 1.0d0
      lmdadt = 0.001d0
      call sethist (1,0.1d0,1.0d0,2.0d0,0.01d0,1.0d0)
      call sethist (2,0.2d0,2.0d0,3.0d0,0.01d0,1.0d0)
      call sethist (3,0.3d0,3.0d0,4.0d0,0.01d0,1.0d0)
c
c     open the history file the way "dynamic" does before sampling
c
      call initostfile
c
c     append two histories, then save again without adding history
c
      nlmdahist = 1
      call saveost
      inquire (file=ostfile,size=size1)
      nlmdahist = 3
      lmdastep = 30
      call saveost
      inquire (file=ostfile,size=size2)
      lmdastep = 31
      call saveost
      inquire (file=ostfile,size=size3)
      call assert_logical (size2.gt.size1,.true.,
     &                     'saveost appends new history')
      call assert_int (size3,size2,'saveost does not duplicate history')
c
c     the per-gaussian heights survive a round trip, and the derived
c     running maxima are rebuilt from the restored history
c
      kelvin = 111.0d0
      call rdost
      call assert_real (kelvin,321.0d0,1.0d-12,
     &                  'rdost restores the temperature')
      call assert_int (nlmdahist,3,'rdost restores the history count')
      call assert_real (osthhist(1),2.0d0,1.0d-12,
     &                  'rdost restores height 1')
      call assert_real (osthhist(2),3.0d0,1.0d-12,
     &                  'rdost restores height 2')
      call assert_real (osthhist(3),4.0d0,1.0d-12,
     &                  'rdost restores height 3')
      do i = 1, nlmda
         call assert_real (vkernelmax(i),brutevkmax(i),1.0d-12,
     &                     'rdost rebuilds the running maximum')
      end do
c
c     check the refreshed header and the sequential history records
c
      ihis = freeunit ()
      open (unit=ihis,file=ostfile,status='old')
      read (ihis,'(a)')  record
      read (ihis,'(a)')  record
      read (ihis,'(a)')  record
      read (record,*)  i,i,i,i,i,nsave0,i
      do i = 1, 5
         read (ihis,'(a)')  record
      end do
      nline = 0
      do
         read (ihis,'(a)',iostat=ios)  record
         if (ios .ne. 0)  exit
         nline = nline + 1
      end do
      close (unit=ihis,status='delete')
      call assert_int (nsave0,3,'saveost updates history count')
      call assert_int (nline,3,'saveost history record count')
      filename = filename0
      leng = leng0
      kelvin = kelvin0
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eost_restart  --  ost continuation test  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eost_restart" checks that a fresh run continues an
c     existing ost history file in place, restoring the gaussians,
c     the kernels, the flambda grid and the lambda particle, and
c     resuming at an interval boundary
c
c
      subroutine test_eost_restart
      use bath
      use dlmda
      use files
      use mutant
      use ost
      implicit none
      integer i,j,ihis
      integer ios
      integer leng0
      integer nline,nsave
      integer nbyte0,nbyte1
      integer freeunit
      real*8 kelvin0
      real*8 deltag0
      real*8 gmax,gdiff
      real*8 efreetot
      real*8 fmean0(5)
      real*8 gkern0(5,5)
      logical exist
      logical stray
      character*240 filename0
      character*240 ostfile
      character*240 record
      character*8000 text0,text1
c
c
c     start an ost history in a temporary file
c
      filename0 = filename
      leng0 = leng
      kelvin0 = kelvin
      kelvin = 300.0d0
      filename = 'tinkertest-ostrestart'
      leng = len_trim(filename)
      ostfile = filename(1:leng)//'.ost'
      inquire (file=ostfile,exist=exist)
      if (exist) then
         ihis = freeunit ()
         open (unit=ihis,file=ostfile,status='old')
         close (unit=ihis,status='delete')
      end if
      call resetost (5,5,8)
      use_ost = .true.
      fastkernel = .true.
      oststdev = 4.0d0
      lmdamass = 2.0d0
      lmdafric = 0.5d0
      lmdadt = 0.002d0
      call sethist (1,0.1d0,1.0d0,2.0d0,0.1d0,1.0d0)
      call sethist (2,0.5d0,2.0d0,3.0d0,0.1d0,1.0d0)
      call sethist (3,0.8d0,-1.0d0,4.0d0,0.1d0,1.0d0)
      do i = 1, 3
         lmdaihist(i) = 10 * i
      end do
      nlmdahist = 3
      call buildostindex
      call buildkernels
      lmdadeltag = efreetot ()
      deltag0 = lmdadeltag
      do i = 1, 5
         fmean0(i) = lmdafmean(i)
         do j = 1, 5
            gkern0(i,j) = gkernel(i,j)
         end do
      end do
c
c     open the file, then save the lambda particle partway into
c     the next interval
c
      lmdastep = 30
      call initostfile
      lmdastep = 37
      lambda = 0.4d0
      lmdatheta = 0.6d0
      lmdavtheta = 0.2d0
      call saveost
c
c     a fresh run with the same settings and another flambda grid
c     continues the file
c
      call resetost (5,9,8)
      use_ost = .true.
      fastkernel = .true.
      oststdev = 4.0d0
      lmdamass = 2.0d0
      lmdafric = 0.5d0
      lmdadt = 0.002d0
      lambda = 1.0d0
      lmdatheta = 0.0d0
      lmdavtheta = 0.0d0
      lmdastep = 0
      call eostfilebytes (ostfile,text0,nbyte0)
      call initostfile
      call eostfilebytes (ostfile,text1,nbyte1)
      call assert_logical (histrewrite,.false.,
     &                     'rdosthist finds no rows past the count')
      call assert_logical (nbyte1.eq.nbyte0 .and. text1.eq.text0,
     &                     .true.,'initostfile leaves a clean file')
      call assert_logical (lmdasavefile.eq.ostfile,.true.,
     &                     'initostfile continues the same file')
      call assert_int (nlmdahist,3,'initostfile restores the gaussians')
      call assert_int (nlmdasave,3,'initostfile restores the save')
      call assert_int (lmdastep,30,'initostfile resumes on an interval')
      call assert_int (nflmda,5,'initostfile restores the flambda bins')
      call assert_int (fli0,3,'initostfile restores the flambda zero')
      call assert_real (lambda,0.4d0,1.0d-14,
     &                  'initostfile restores the lambda')
      call assert_real (lmdatheta,0.6d0,1.0d-14,
     &                  'initostfile restores the lambda theta')
      call assert_real (lmdavtheta,0.2d0,1.0d-14,
     &                  'initostfile restores the theta velocity')
      call assert_array1 (lmdafmean,fmean0,5,1.0d-12,
     &                    'initostfile restores the mean force')
      gmax = 0.0d0
      gdiff = 0.0d0
      do i = 1, 5
         do j = 1, 5
            gmax = max(gmax,abs(gkern0(i,j)))
            gdiff = max(gdiff,abs(gkernel(i,j)-gkern0(i,j)))
         end do
      end do
      call assert_logical (gmax.gt.0.0d0,.true.,
     &                     'initostfile test g kernel is nonzero')
      call assert_real (gdiff,0.0d0,1.0d-12,
     &                  'initostfile restores the g kernel')
      call assert_real (lmdadeltag,deltag0,1.0d-12,
     &                  'initostfile restores the free energy')
c
c     a new gaussian is appended to the same file
c
      call sethist (4,0.3d0,0.0d0,1.0d0,0.1d0,1.0d0)
      lmdaihist(4) = 40
      nlmdahist = 4
      call saveost
      inquire (file=filename(1:leng)//'.ost_2',exist=exist)
      call assert_logical (exist,.false.,
     &                     'initostfile claims no new version')
      ihis = freeunit ()
      open (unit=ihis,file=ostfile,status='old')
      do i = 1, 3
         read (ihis,'(a)')  record
      end do
      read (record,*)  i,i,i,i,i,nsave,i
      do i = 1, 5
         read (ihis,'(a)')  record
      end do
      nline = 0
      do
         read (ihis,'(a)',iostat=ios)  record
         if (ios .ne. 0)  exit
         nline = nline + 1
      end do
      close (unit=ihis)
      call assert_int (nsave,4,'continued ost header history count')
      call assert_int (nline,4,'continued ost history record count')
c
c     a row past the header count, as left by a run stopped between
c     appending gaussians and updating the header, is dropped when
c     the history is continued
c
      ihis = freeunit ()
      open (unit=ihis,file=ostfile,status='old',position='append')
      write (ihis,'(a)')  '          50       1  stray row'
      close (unit=ihis)
      call resetost (5,9,8)
      use_ost = .true.
      fastkernel = .true.
      oststdev = 4.0d0
      lmdamass = 2.0d0
      lmdafric = 0.5d0
      lmdadt = 0.002d0
      call initostfile
      call assert_logical (histrewrite,.true.,
     &                     'rdosthist finds the stray row')
      call assert_int (nlmdahist,4,'initostfile skips the stray row')
      call sethist (5,0.6d0,0.0d0,1.0d0,0.1d0,1.0d0)
      lmdaihist(5) = 50
      nlmdahist = 5
      call saveost
      ihis = freeunit ()
      open (unit=ihis,file=ostfile,status='old')
      nline = 0
      stray = .false.
      do
         read (ihis,'(a)',iostat=ios)  record
         if (ios .ne. 0)  exit
         if (index(record,'stray') .gt. 0)  stray = .true.
         nline = nline + 1
      end do
      close (unit=ihis,status='delete')
      call assert_logical (stray,.false.,
     &                     'continued ost history drops the stray row')
      call assert_int (nline-8,5,'ost history keeps the other rows')
c
c     restore the global state
c
      filename = filename0
      leng = leng0
      kelvin = kelvin0
      use_ost = .false.
      fastkernel = .false.
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eost_metarestart  --  meta continuation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eost_metarestart" checks that a fresh run continues an
c     existing metadynamics history file in place, restoring the
c     gaussians, the bias grid and the lambda particle, and resuming
c     at an interval boundary
c
c
      subroutine test_eost_metarestart
      use bath
      use dlmda
      use files
      use mutant
      use ost
      implicit none
      integer i,ihis
      integer ios
      integer leng0
      integer nline,nsave
      integer nbyte0,nbyte1
      integer freeunit
      real*8 kelvin0
      real*8 deltag0
      real*8 metadeltag
      real*8 vgrid0(5)
      real*8 dvgrid0(5)
      real*8 lam(3)
      logical exist
      character*240 filename0
      character*240 metafile
      character*240 record
      character*8000 text0,text1
      data lam / 0.1d0,0.5d0,0.8d0 /
c
c
c     start a metadynamics history in a temporary file
c
      filename0 = filename
      leng0 = leng
      kelvin0 = kelvin
      kelvin = 300.0d0
      filename = 'tinkertest-metarestart'
      leng = len_trim(filename)
      metafile = filename(1:leng)//'.meta'
      inquire (file=metafile,exist=exist)
      if (exist) then
         ihis = freeunit ()
         open (unit=ihis,file=metafile,status='old')
         close (unit=ihis,status='delete')
      end if
      call resetost (5,5,1)
      use_meta = .true.
      call resetmeta (8)
      oststdev = 4.0d0
      lmdamass = 2.0d0
      lmdafric = 0.5d0
      lmdadt = 0.002d0
      do i = 1, 3
         metaihist(i) = 10 * i
         metalhist(i) = lam(i)
         metahhist(i) = dble(i+1)
         metawhist(i) = 0.1d0
      end do
      nmetahist = 3
      maxwlhist = 0.1d0
      call buildmetaindex
      do i = 1, nmetahist
         call addmetagrid (i)
      end do
      lmdadeltag = metadeltag ()
      deltag0 = lmdadeltag
      do i = 1, 5
         vgrid0(i) = vmetagrid(i)
         dvgrid0(i) = dvmetagrid(i)
      end do
c
c     open the file, then save the lambda particle partway into
c     the next interval
c
      lmdastep = 30
      call initmetafile
      lmdastep = 37
      lambda = 0.4d0
      lmdatheta = 0.6d0
      lmdavtheta = 0.2d0
      call savemeta
c
c     a fresh run with the same settings continues the file
c
      call resetost (5,5,1)
      use_meta = .true.
      call resetmeta (2)
      oststdev = 4.0d0
      lmdamass = 2.0d0
      lmdafric = 0.5d0
      lmdadt = 0.002d0
      lambda = 1.0d0
      lmdatheta = 0.0d0
      lmdavtheta = 0.0d0
      lmdastep = 0
      call eostfilebytes (metafile,text0,nbyte0)
      call initmetafile
      call eostfilebytes (metafile,text1,nbyte1)
      call assert_logical (histrewrite,.false.,
     &                     'rdmetahist finds no rows past the count')
      call assert_logical (nbyte1.eq.nbyte0 .and. text1.eq.text0,
     &                     .true.,'initmetafile leaves a clean file')
      call assert_logical (metasavefile.eq.metafile,.true.,
     &                     'initmetafile continues the same file')
      call assert_int (nmetahist,3,
     &                 'initmetafile restores the gaussians')
      call assert_int (nmethistsave,3,'initmetafile restores the save')
      call assert_int (lmdastep,30,
     &                 'initmetafile resumes on an interval')
      call assert_real (lambda,0.4d0,1.0d-14,
     &                  'initmetafile restores the lambda')
      call assert_real (lmdatheta,0.6d0,1.0d-14,
     &                  'initmetafile restores the lambda theta')
      call assert_real (lmdavtheta,0.2d0,1.0d-14,
     &                  'initmetafile restores the theta velocity')
      call assert_array1 (vmetagrid,vgrid0,5,1.0d-12,
     &                    'initmetafile restores the bias grid')
      call assert_array1 (dvmetagrid,dvgrid0,5,1.0d-12,
     &                    'initmetafile restores the slope grid')
      call assert_real (lmdadeltag,deltag0,1.0d-12,
     &                  'initmetafile restores the free energy')
c
c     a new gaussian is appended to the same file
c
      metaihist(4) = 40
      metalhist(4) = 0.3d0
      metahhist(4) = 1.0d0
      metawhist(4) = 0.1d0
      nmetahist = 4
      call savemeta
      inquire (file=filename(1:leng)//'.meta_2',exist=exist)
      call assert_logical (exist,.false.,
     &                     'initmetafile claims no new version')
      ihis = freeunit ()
      open (unit=ihis,file=metafile,status='old')
      do i = 1, 3
         read (ihis,'(a)')  record
      end do
      read (record,*)  i,i,nsave,i
      do i = 1, 5
         read (ihis,'(a)')  record
      end do
      nline = 0
      do
         read (ihis,'(a)',iostat=ios)  record
         if (ios .ne. 0)  exit
         nline = nline + 1
      end do
      close (unit=ihis,status='delete')
      call assert_int (nsave,4,'continued meta header history count')
      call assert_int (nline,4,'continued meta history record count')
c
c     restore the global state
c
      filename = filename0
      leng = leng0
      kelvin = kelvin0
      use_meta = .false.
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_eost_restartrun  --  dynamic restarts  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_eost_restartrun" runs "dynamic" twice for each adaptive
c     lambda bias method, the second run continuing the history and
c     trajectory of the first, and checks the runs that are refused
c
c
      subroutine test_eost_restartrun
      implicit none
c
c
c     continue each method from the files of a previous run
c
      call eostrstrun ('OST','ost')
      call eostrstrun ('META','meta')
      call eostrstrun ('ABF','abf')
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine eostrstrun  --  restart one dynamics run  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "eostrstrun" makes two 20 step dynamics runs of one lambda bias
c     method in a scratch directory, the second continuing the first;
c     it then checks that a run whose settings differ from the stored
c     history and a run without a thermostat are refused and leave the
c     history untouched, and that a run that does not end on a save
c     step writes its restart file at its last step, as the history
c
c
      subroutine eostrstrun (mode,ext)
      implicit none
      integer ist,istep
      integer nval,nrow,nsaved
      real*8 tlmda
      logical exist
      character*(*) mode,ext
      character*40 tag
      character*240 histfile
      character*512 cmd
c
c
c     copy the fixture, pointing its key at the parameter file
c
      tag = 'eost restart '//trim(mode)
      histfile = 'water2.'//ext
      call pushdir ('file/mutate')
      cmd = 'rm -rf ../eost_rstrun ; mkdir -p ../eost_rstrun ; '//
     &      'cp water2.xyz ../eost_rstrun/ ; '//
     &      'cp vsoft_default_ost.key ../eost_rstrun/water2.key'
      call execute_command_line (cmd)
      call popdir
      call pushdir ('file/eost_rstrun')
      call tnist_append ('water2.key','lambda-mode '//mode)
      call tnist_append ('water2.key','lambda-interval 10')
      call tnist_append ('water2.key','randomseed 12345')
c
c     a first run starts the history, and a second one continues it
c
      call run_prog ('dynamic','water2 20 1.0 0.01 2 300','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,trim(tag)//' first status')
         cmd = 'grep -q "Bias History Written To" out.txt'
         call execute_command_line (cmd,exitstat=ist)
         call assert_int (ist,0,trim(tag)//' first starts the history')
         call run_prog ('dynamic','water2 20 1.0 0.01 2 300',
     &                  'out.txt',ist)
         call assert_int (ist,0,trim(tag)//' second status')
         cmd = 'grep -q "Bias History Continued From" out.txt'
         call execute_command_line (cmd,exitstat=ist)
         call assert_int (ist,0,trim(tag)//' second continues')
         cmd = 'grep -q "Restarting Molecular Dynamics" out.txt'
         call execute_command_line (cmd,exitstat=ist)
         call assert_int (ist,0,trim(tag)//' second reads the dyn')
         inquire (file=trim(histfile)//'_2',exist=exist)
         call assert_logical (exist,.false.,
     &                        trim(tag)//' no new version')
c
c     a stored setting that differs stops the run, for abf the
c     interval and for the others the temperature
c
         if (mode .eq. 'ABF') then
            cmd = 'cp water2.key water2.bak ; '//
     &            'echo "lambda-interval 5" >> water2.key'
            call execute_command_line (cmd)
            call run_prog ('dynamic','water2 20 1.0 0.01 2 300',
     &                     'out.txt',ist)
            call execute_command_line ('mv water2.bak water2.key')
            cmd = 'grep -q "must match the ABF" out.txt'
         else
            call run_prog ('dynamic','water2 20 1.0 0.01 2 310',
     &                     'out.txt',ist)
            cmd = 'grep -q "Setting Differs from the Bias" out.txt'
         end if
         call execute_command_line (cmd,exitstat=ist)
         call assert_int (ist,0,trim(tag)//' setting refused')
c
c     an ost or metadynamics history stores the lambda temperature,
c     which follows the bath without the keyword, and a run with a
c     different one is refused
c
         if (mode .ne. 'ABF') then
            call eostlmdastate (histfile,nval,tlmda)
            call assert_int (nval,8,trim(tag)//' stores lambda temp')
            call assert_real (tlmda,300.0d0,0.0d0,
     &                        trim(tag)//' lambda temp follows bath')
            cmd = 'cp water2.key water2.bak ; '//
     &            'echo "lambda-temperature 600" >> water2.key'
            call execute_command_line (cmd)
            call run_prog ('dynamic','water2 20 1.0 0.01 2 300',
     &                     'out.txt',ist)
            call execute_command_line ('mv water2.bak water2.key')
            cmd = 'grep -q "LAMBDA-TEMPERATURE" out.txt'
            call execute_command_line (cmd,exitstat=ist)
            call assert_int (ist,0,trim(tag)//' lambda temp refused')
         end if
c
c     metadynamics refuses Monte Carlo lambda moves away from the
c     bath temperature
c
         if (mode .eq. 'META') then
            cmd = 'cp water2.key water2.bak ; '//
     &            'echo "lambda-propagate mc" >> water2.key ; '//
     &            'echo "lambda-temperature 600" >> water2.key'
            call execute_command_line (cmd)
            call run_prog ('dynamic','water2 20 1.0 0.01 2 300',
     &                     'out.txt',ist)
            call execute_command_line ('mv water2.bak water2.key')
            cmd = 'grep -q "Monte Carlo Lambda Moves" out.txt'
            call execute_command_line (cmd,exitstat=ist)
            call assert_int (ist,0,trim(tag)//' mc lambda temp refused')
         end if
c
c     a run without a thermostat is refused
c
         call run_prog ('dynamic','water2 20 1.0 0.01 1',
     &                  'out.txt',ist)
         cmd = 'grep -q "a Constant Temperature Ensemble" out.txt'
         call execute_command_line (cmd,exitstat=ist)
         call assert_int (ist,0,trim(tag)//' nve refused')
c
c     the refused runs leave the continued history at step 40
c
         call eosthiststep (histfile,istep)
         call assert_int (istep,40,trim(tag)//' history step')
c
c     a run shorter than the save interval saves nothing on the
c     way, so the restart file is written at its last step, along
c     with the history
c
         call execute_command_line ('rm -f water2.dyn')
         call run_prog ('dynamic','water2 5 1.0 0.01 2 300',
     &                  'out.txt',ist)
         call assert_int (ist,0,trim(tag)//' short status')
         inquire (file='water2.dyn',exist=exist)
         call assert_logical (exist,.true.,
     &                        trim(tag)//' short writes the dyn')
         call eosthiststep (histfile,istep)
         call assert_int (istep,45,trim(tag)//' short history step')
c
c     the next run continues both, the history from its last full
c     interval
c
         call run_prog ('dynamic','water2 20 1.0 0.01 2 300',
     &                  'out.txt',ist)
         call assert_int (ist,0,trim(tag)//' third status')
         cmd = 'grep -q "Restarting Molecular Dynamics" out.txt'
         call execute_command_line (cmd,exitstat=ist)
         call assert_int (ist,0,trim(tag)//' third reads the dyn')
         call eosthiststep (histfile,istep)
         call assert_int (istep,60,trim(tag)//' third history step')
c
c     a history written before the lambda temperature was stored
c     continues at the bath temperature, and its header is widened
c     by rewriting the whole file rather than over its first rows
c
         if (mode .ne. 'ABF') then
c
c     a lambda state record cut short, as by a save interrupted
c     while it rewrote the header, is refused and left in place
c
            cmd = 'cp '//trim(histfile)//' '//trim(histfile)//'.bak'
            call execute_command_line (cmd)
            call eostoldhead (histfile,5)
            cmd = 'cp '//trim(histfile)//' '//trim(histfile)//'.cut'
            call execute_command_line (cmd)
            call run_prog ('dynamic','water2 20 1.0 0.01 2 300',
     &                     'out.txt',ist)
            cmd = 'grep -q "Error while Reading Lambda Bias" out.txt'
            call execute_command_line (cmd,exitstat=ist)
            call assert_int (ist,0,trim(tag)//' short header refused')
            cmd = 'cmp -s '//trim(histfile)//' '//trim(histfile)//'.cut'
            call execute_command_line (cmd,exitstat=ist)
            call assert_int (ist,0,trim(tag)//' short header kept')
            cmd = 'mv '//trim(histfile)//'.bak '//trim(histfile)//
     &            ' ; rm -f '//trim(histfile)//'.cut'
            call execute_command_line (cmd)
c
c     with all seven older values the history continues
c
            call eostoldhead (histfile,7)
            call eostlmdastate (histfile,nval,tlmda)
            call assert_int (nval,7,trim(tag)//' old header made')
            call run_prog ('dynamic','water2 20 1.0 0.01 2 300',
     &                     'out.txt',ist)
            call assert_int (ist,0,trim(tag)//' old header status')
            cmd = 'grep -q "Bias History Continued From" out.txt'
            call execute_command_line (cmd,exitstat=ist)
            call assert_int (ist,0,trim(tag)//' old header continues')
            call eostlmdastate (histfile,nval,tlmda)
            call assert_int (nval,8,trim(tag)//' old header widened')
            call assert_real (tlmda,300.0d0,0.0d0,
     &                        trim(tag)//' old header at the bath')
            call eosthistrows (mode,histfile,nrow,nsaved)
            call assert_int (nrow,nsaved,trim(tag)//' old header rows')
            call run_prog ('dynamic','water2 20 1.0 0.01 2 300',
     &                     'out.txt',ist)
            call assert_int (ist,0,trim(tag)//' widened file reads')
            call eosthiststep (histfile,istep)
            call assert_int (istep,100,trim(tag)//' widened file step')
c
c     a new history run with its own lambda temperature stores it
c
            call execute_command_line ('rm -f water2.dyn '//
     &                                 trim(histfile)//'*')
            call tnist_append ('water2.key','lambda-temperature 900')
            call run_prog ('dynamic','water2 20 1.0 0.01 2 300',
     &                     'out.txt',ist)
            call assert_int (ist,0,trim(tag)//' hot lambda status')
            call eostlmdastate (histfile,nval,tlmda)
            call assert_real (tlmda,900.0d0,0.0d0,
     &                        trim(tag)//' hot lambda stored')
            call eosthiststep (histfile,istep)
            call assert_int (istep,20,trim(tag)//' hot lambda step')
         end if
      end if
      call popdir
      call tnist_clean ('eost_rstrun')
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine eosthiststep  --  history step count  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "eosthiststep" returns the dynamics step count stored in the
c     header of a lambda bias history file, the first value of its
c     integer state record
c
c
      subroutine eosthiststep (histfile,istep)
      implicit none
      integer ihis,istep
      integer freeunit
      character*240 record
      character*(*) histfile
c
c
c     read the third header record, which holds the integer state
c
      ihis = freeunit ()
      open (unit=ihis,file=histfile,status='old')
      read (ihis,'(a)')  record
      read (ihis,'(a)')  record
      read (ihis,'(a)')  record
      close (unit=ihis)
      read (record,*)  istep
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eostlmdastate  --  stored lambda temperature  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eostlmdastate" returns the number of values in the lambda
c     state record of an ost or metadynamics history header, and the
c     lambda particle temperature it stores when it holds eight
c
c
      subroutine eostlmdastate (histfile,nval,tlmda)
      implicit none
      integer i,ihis,nval
      integer freeunit
      real*8 tlmda
      real*8 state(8)
      character*240 record
      character*(*) histfile
c
c
c     the lambda state values follow the lambda state label
c
      ihis = freeunit ()
      open (unit=ihis,file=histfile,status='old')
      do i = 1, 6
         read (ihis,'(a)')  record
      end do
      read (ihis,'(a)')  record
      close (unit=ihis)
      nval = 7
      tlmda = 0.0d0
      read (record,*,err=10,end=10)  state
      nval = 8
      tlmda = state(8)
   10 continue
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eostoldhead  --  drop the stored lambda temp  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eostoldhead" rewrites an ost or metadynamics history keeping
c     only the first "nkeep" values of its lambda state record, which
c     with seven gives the record written before the lambda particle
c     temperature was stored, and with fewer a truncated record; every
c     other line is kept
c
c
      subroutine eostoldhead (histfile,nkeep)
      implicit none
      integer i,ihis,nline
      integer nkeep
      integer freeunit
      character*240 lines(5000)
      character*(*) histfile
c
c
c     read every line, then cut the lambda state values past the
c     ones kept
c
      ihis = freeunit ()
      open (unit=ihis,file=histfile,status='old')
      nline = 0
      do while (nline .lt. 5000)
         read (ihis,'(a)',end=10)  lines(nline+1)
         nline = nline + 1
      end do
   10 continue
      close (unit=ihis)
      call assert_logical (nline.lt.5000,.true.,
     &                     'eostoldhead buffer holds the file')
      lines(7)(nkeep*26+1:240) = ' '
      open (unit=ihis,file=histfile,status='replace')
      do i = 1, nline
         write (ihis,'(a)')  trim(lines(i))
      end do
      close (unit=ihis)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eosthistrows  --  saved and present row count  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eosthistrows" returns the number of rows after the eight
c     header records of an ost or metadynamics history, and the
c     count of saved rows its integer state record claims
c
c
      subroutine eosthistrows (mode,histfile,nrow,nsaved)
      implicit none
      integer i,ihis
      integer nrow,nsaved
      integer freeunit
      integer ival(7)
      character*240 record
      character*(*) mode
      character*(*) histfile
c
c
c     the saved count is the third integer of a metadynamics header
c     and the sixth of an ost header
c
      ihis = freeunit ()
      open (unit=ihis,file=histfile,status='old')
      do i = 1, 3
         read (ihis,'(a)')  record
      end do
      if (mode .eq. 'META') then
         read (record,*)  (ival(i),i=1,4)
         nsaved = ival(3)
      else
         read (record,*)  (ival(i),i=1,7)
         nsaved = ival(6)
      end if
      do i = 4, 8
         read (ihis,'(a)')  record
      end do
      nrow = 0
      do while (.true.)
         read (ihis,'(a)',end=10)  record
         nrow = nrow + 1
      end do
   10 continue
      close (unit=ihis)
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine eostfilebytes  --  whole file content  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "eostfilebytes" reads the whole content of a small history file
c     into a text buffer, so tests can check a file is left unchanged
c
c
      subroutine eostfilebytes (histfile,text,nbyte)
      implicit none
      integer ihis,nbyte
      integer freeunit
      character*(*) histfile
      character*(*) text
c
c
c     read the file as a single unformatted stream
c
      inquire (file=histfile,size=nbyte)
      text = ' '
      ihis = freeunit ()
      open (unit=ihis,file=histfile,status='old',access='stream',
     &      form='unformatted')
      read (ihis)  text(1:min(nbyte,len(text)))
      close (unit=ihis)
      call assert_logical (nbyte.le.len(text),.true.,
     &                     'eostfilebytes buffer holds the file')
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_eost_index  --  packed index tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_eost_index" checks conversion between lambda
c     grid coordinates and packed histogram indices
c
c
      subroutine test_eost_index
      implicit none
      integer i,j,k
c
c
c     test ij_to_k and k_to_ij
c
      call ij_to_k (3,4,7,k)
      call assert_int (k,24,
     &                 'ij_to_k input i=3 j=4 nrow=7')
      call k_to_ij (24,7,i,j)
      call assert_int (i,3,
     &                 'k_to_ij input k=24 nrow=7 i')
      call assert_int (j,4,
     &                 'k_to_ij input k=24 nrow=7 j')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_eost_resize  --  history resize tests  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_eost_resize" checks OST history resizing preserves
c     old entries and initializes new storage
c
c
      subroutine test_eost_resize
      use dlmda
      use ost
      implicit none
      integer i
c
c
c     set small history arrays and resize them
c
      call resetost (5,5,2)
      nlmdahist = 2
      do i = 1, 2
         osthist(i) = 10 + i
         ostnext(i) = i - 1
         lmdalhist(i) = 0.25d0 * dble(i)
         lmdafhist(i) = -3.0d0 + 2.0d0*dble(i)
         osthhist(i) = 1.0d0 + dble(i)
         ostwlhist(i) = 0.25d0
         ostwfhist(i) = 1.0d0
      end do
      call resizeosthist
c
c     check old data are preserved and new slots are initialized
c
      call assert_int (sizelmdahist,4,
     &                 'resizeosthist size input size=2')
      do i = 1, 4
         if (i .le. 2) then
      call assert_int (osthist(i),10+i,
     &                 'resizeosthist preserve osthist')
      call assert_int (ostnext(i),i-1,
     &                 'resizeosthist preserve ostnext')
      call assert_real (lmdalhist(i),0.25d0*dble(i),1.0d-12,
     &                  'resizeosthist preserve lmdalhist')
      call assert_real (lmdafhist(i),-3.0d0+2.0d0*dble(i),1.0d-12,
     &                  'resizeosthist preserve lmdafhist')
      call assert_real (osthhist(i),1.0d0+dble(i),1.0d-12,
     &                  'resizeosthist preserve osthhist')
      call assert_real (ostwlhist(i),0.25d0,1.0d-12,
     &                  'resizeosthist preserve ostwlhist')
      call assert_real (ostwfhist(i),1.0d0,1.0d-12,
     &                  'resizeosthist preserve ostwfhist')
         else
      call assert_int (osthist(i),0,
     &                 'resizeosthist init osthist')
      call assert_int (ostnext(i),0,
     &                 'resizeosthist init ostnext')
      call assert_real (lmdalhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init lmdalhist')
      call assert_real (lmdafhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init lmdafhist')
      call assert_real (osthhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init osthhist')
      call assert_real (ostwlhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init ostwlhist')
      call assert_real (ostwfhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init ostwfhist')
         end if
      end do
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine test_eost_buildindex  --  index tests  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "test_eost_buildindex" checks linked-list construction
c     for gaussian history bins
c
c
      subroutine test_eost_buildindex
      use dlmda
      use ost
      implicit none
      integer k,ilmda,iflmda
c
c
c     put two gaussians in the same bin and one in another bin
c
      call resetost (5,5,3)
      nlmdahist = 3
      call sethist (1,0.50d0,0.0d0,1.0d0,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,2.0d0,wlmda,wflmda)
      call sethist (3,0.75d0,1.0d0,3.0d0,wlmda,wflmda)
      call buildostindex
c
c     expected linked list head is newest entry in the shared bin
c
      ilmda = 3
      iflmda = 3
      call ij_to_k (ilmda,iflmda,nlmda,k)
      call assert_int (osthist(1),k,
     &                 'buildostindex packed shared bin')
      call assert_int (osthead(ilmda,iflmda),2,
     &                 'buildostindex head shared bin')
      call assert_int (ostnext(2),1,
     &                 'buildostindex next shared bin')
      call assert_int (ostnext(1),0,
     &                 'buildostindex tail shared bin')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_eost_ensure  --  flambda resize tests  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_eost_ensure" checks flambda grid expansion,
c     kernel preservation and index rebuilding
c
c
      subroutine test_eost_ensure
      use dlmda
      use ost
      implicit none
      integer k
      integer nold
      integer oldfli0
      integer ilmda,iflmda
c
c
c     high-side expansion preserves fli0 and old gkernel values
c
      call resetost (3,5,1)
      gkernel(2,3) = 7.0d0
      call ensureflambda (2000.0d0)
      call assert_int (nflmda,2105,
     &                 'ensureflambda high nflmda')
      call assert_int (fli0,3,
     &                 'ensureflambda high fli0')
      call assert_real (gkernel(2,3),7.0d0,1.0d-12,
     &                  'ensureflambda high preserve gkernel')
c
c     low-side expansion shifts fli0 and old gkernel values
c
      call resetost (3,5,1)
      gkernel(2,3) = 7.0d0
      call ensureflambda (-2000.0d0)
      call assert_int (nflmda,2105,
     &                 'ensureflambda low nflmda')
      call assert_int (fli0,2103,
     &                 'ensureflambda low fli0')
      call assert_real (gkernel(2,2103),7.0d0,1.0d-12,
     &                  'ensureflambda low shifted gkernel')
c
c     low-side expansion also rebuilds osthead, ostnext and osthist
c
      call resetost (3,5,2)
      nlmdahist = 2
      call sethist (1,0.50d0,0.0d0,1.0d0,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,2.0d0,wlmda,wflmda)
      call buildostindex
      call ensureflambda (-2000.0d0)
c
c     after low-side expansion, fli0 shifts from 3 to 2103, so
c     flambda=0.0 now belongs to iflmda=2103
c
      ilmda = 2
      iflmda = 2103
      call ij_to_k (ilmda,iflmda,nlmda,k)
      call assert_int (osthist(1),k,
     &                 'ensureflambda rebuild osthist(1)')
      call assert_int (osthist(2),k,
     &                 'ensureflambda rebuild osthist(2)')
      call assert_int (osthead(ilmda,iflmda),2,
     &                 'ensureflambda rebuild osthead')
      call assert_int (ostnext(2),1,
     &                 'ensureflambda rebuild ostnext head')
      call assert_int (ostnext(1),0,
     &                 'ensureflambda rebuild ostnext tail')
c
c     the resize copies four separate kernels, so seed all of them
c     with values that identify both the array and the bin, then
c     require every entry of every array to land where it belongs
c
      call resetost (3,5,1)
      nold = nflmda
      oldfli0 = fli0
      call seedkernels
      call ensureflambda (2000.0d0)
      call checkkernels ('ensureflambda high',nold,fli0-oldfli0)
      call resetost (3,5,1)
      nold = nflmda
      oldfli0 = fli0
      call seedkernels
      call ensureflambda (-2000.0d0)
      call checkkernels ('ensureflambda low',nold,fli0-oldfli0)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_eost_gkernels  --  g kernel tests  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "test_eost_gkernels" checks gaussian kernel updates,
c     mirroring, interpolation and lookup behavior
c
c
      subroutine test_eost_gkernels
      use dlmda
      use math
      use mutant
      use ost
      implicit none
      real*8 egbias,dgdl,dgdfl
      real*8 expected
      real*8 height
      real*8 targetl,targetf
c
c
c     choose height so the normalized gaussian prefactor is one
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
c
c     addkernelhist g from a single interior gaussian
c
      call addkernelhist (1,.false.)
      call assert_real (gkernel(3,3),1.0d0,1.0d-12,
     &                  'addkernelhist g center value')
      expected = exp(-1.0d0)
      call assert_real (gkernel(4,4),expected,1.0d-12,
     &                  'addkernelhist g off-center value')
c
c     addkernelhist g includes left-boundary mirror image
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,0.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      call addkernelhist (1,.false.)
      call assert_real (gkernel(1,3),2.0d0,1.0d-12,
     &                  'addkernelhist g left mirror center')
      expected = 2.0d0 * exp(-0.5d0)
      call assert_real (gkernel(2,3),expected,1.0d-12,
     &                  'addkernelhist g left mirror neighbor')
c
c     addkernelhist g includes right-boundary mirror image
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,1.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      call addkernelhist (1,.false.)
      call assert_real (gkernel(5,3),2.0d0,1.0d-12,
     &                  'addkernelhist g right mirror center')
      expected = 2.0d0 * exp(-0.5d0)
      call assert_real (gkernel(4,3),expected,1.0d-12,
     &                  'addkernelhist g right mirror neighbor')
c
c     rebuild the original interior gaussian for later checks
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
c
c     buildgkernel rebuilds the same full grid
c
      call buildgkernel
      call assert_real (gkernel(3,3),1.0d0,1.0d-12,
     &                  'buildgkernel center value')
      expected = exp(-1.0d0)
      call assert_real (gkernel(4,4),expected,1.0d-12,
     &                  'buildgkernel off-center value')
c
c     updategkernel adds only the newest gaussian to current grid
c
      gkernel = 0.0d0
      call updategkernel
      call assert_real (gkernel(3,3),1.0d0,1.0d-12,
     &                  'updategkernel center value')
c
c     adding a taller gaussian updates only from the new history entry
c
      nlmdahist = 2
      call sethist (2,0.75d0,1.0d0,2.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
      expected = 2.0d0 + exp(-1.0d0)
      call assert_real (gkernel(4,4),expected,1.0d-12,
     &                  'updategkernel added new center')
      expected = 1.0d0 + 2.0d0*exp(-1.0d0)
      call assert_real (gkernel(3,3),expected,1.0d-12,
     &                  'updategkernel keeps old plus new overlap')
c
c     egkernel evaluates the continuous gaussian and derivatives
c
      lambda = 0.75d0
      dedl = 1.0d0
      call egkernel (egbias,dgdl,dgdfl)
      expected = 2.0d0 + exp(-1.0d0)
      call assert_real (egbias,expected,1.0d-12,
     &                  'egkernel bias lambda=.75 flambda=1')
      call assert_real (dgdl,-4.0d0*exp(-1.0d0),1.0d-12,
     &                  'egkernel dgdl lambda=.75 flambda=1')
      call assert_real (dgdfl,-exp(-1.0d0),1.0d-12,
     &                  'egkernel dgdfl lambda=.75 flambda=1')
c
c     two gaussians in the same bin are both followed by ostnext
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 2
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      lambda = 0.50d0
      dedl = 0.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call assert_real (egbias,2.0d0,1.0d-12,
     &                  'egkernel same-bin linked bias')
      call assert_real (dgdl,0.0d0,1.0d-12,
     &                  'egkernel same-bin linked dgdl')
      call assert_real (dgdfl,0.0d0,1.0d-12,
     &                  'egkernel same-bin linked dgdfl')
c
c     multiple bins are found through osthead lookup, including
c     two gaussians in one bin and one gaussian in another bin
c
      call resetost (5,5,4)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 3
      call sethist (1,0.50d0,0.0d0,1.0d0*height,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,2.0d0*height,wlmda,wflmda)
      call sethist (3,0.75d0,1.0d0,3.0d0*height,wlmda,wflmda)
      call buildostindex
      lambda = 0.75d0
      dedl = 1.0d0
      call egkernel (egbias,dgdl,dgdfl)
      expected = 3.0d0 + 3.0d0*exp(-1.0d0)
      call assert_real (egbias,expected,1.0d-12,
     &                  'egkernel multi-bin mixed bias')
      call assert_real (dgdl,-12.0d0*exp(-1.0d0),1.0d-12,
     &                  'egkernel multi-bin mixed dgdl')
      call assert_real (dgdfl,-3.0d0*exp(-1.0d0),1.0d-12,
     &                  'egkernel multi-bin mixed dgdfl')
c
c     left endpoint includes both real and mirror gaussian images
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,0.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      lambda = 0.0d0
      dedl = 0.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call assert_real (egbias,2.0d0,1.0d-12,
     &                  'egkernel left mirror bias')
      call assert_real (dgdl,0.0d0,1.0d-12,
     &                  'egkernel left mirror dgdl')
c
c     right endpoint includes both real and mirror gaussian images
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,1.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      lambda = 1.0d0
      dedl = 0.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call assert_real (egbias,2.0d0,1.0d-12,
     &                  'egkernel right mirror bias')
      call assert_real (dgdl,0.0d0,1.0d-12,
     &                  'egkernel right mirror dgdl')
c
c     wide histogram widths are used for both grid and continuous bias
c
      call resetost (41,81,3)
      wlhist = 0.05d0
      wfhist = 10.0d0
      maxwlhist = wlhist
      maxwfhist = wfhist
      height = 2.0d0 * pi * wlhist * wfhist
      targetl = 0.525d0
      targetf = 5.0d0
      nlmdahist = 1
      call sethist (1,0.50d0,0.0d0,height,wlhist,wfhist)
      call buildostindex
      call buildgkernel
      expected = exp(-0.25d0)
      call assert_real (gkernel(22,46),expected,1.0d-12,
     &                  'buildgkernel wide hist width')
      lambda = targetl
      dedl = targetf
      call egkernel (egbias,dgdl,dgdfl)
      call assert_real (egbias,expected,1.0d-12,
     &                  'egkernel wide hist width bias')
      call assert_real (dgdl,-10.0d0*expected,1.0d-12,
     &                  'egkernel wide hist width dgdl')
      call assert_real (dgdfl,-0.05d0*expected,1.0d-12,
     &                  'egkernel wide hist width dgdfl')
c
c     if dU/dlambda is outside the current flambda grid,
c     egkernel returns zero
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      lambda = 0.50d0
      dedl = 100.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call assert_real (egbias,0.0d0,1.0d-12,
     &                  'egkernel outside flambda bias')
      call assert_real (dgdl,0.0d0,1.0d-12,
     &                  'egkernel outside flambda dgdl')
      call assert_real (dgdfl,0.0d0,1.0d-12,
     &                  'egkernel outside flambda dgdfl')
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine test_eost_fkernel  --  f kernel tests  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "test_eost_fkernel" checks mean-force kernel
c     construction from gaussian history data
c
c
      subroutine test_eost_fkernel
      use bath
      use dlmda
      use math
      use ost
      use units
      implicit none
      real*8 rt,rth
      real*8 expected
      real*8 hotmean
      real*8 height
      real*8 wratio
      real*8 w1,w2,w3,w4,w5
c
c     two nonzero gkernel weights at flambda=-1 and +1
c
      call resetost (5,5,1)
      rt = gasconst * kelvin
      gkernel(3,2) = log(2.0d0) * rt
      gkernel(3,4) = log(4.0d0) * rt
      vkernelmax(3) = log(4.0d0) * rt
      call buildfkernel
      expected = 1.0d0 / 3.0d0
      call assert_real (lmdafmean(3),expected,1.0d-12,
     &                  'buildfkernel weighted mean force')
      call assert_real (lmdafmean(1),0.0d0,1.0d-12,
     &                  'buildfkernel empty row')
c
c     build gkernel incrementally from multiple gaussians using
c     updategkernel, then construct the f kernel
c
      call resetost (5,5,4)
      rt = gasconst * kelvin
      height = 2.0d0 * pi * wlmda * wflmda
c
      nlmdahist = 1
      call sethist (1,0.25d0,-1.0d0,1.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
c
      nlmdahist = 2
      call sethist (2,0.50d0, 0.0d0,2.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
c
      nlmdahist = 3
      call sethist (3,0.75d0, 1.0d0,3.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
c
      call buildfkernel
c
c     expected weighted mean at lambda=0.5
c
      w1 = exp(gkernel(3,1)/rt)
      w2 = exp(gkernel(3,2)/rt)
      w3 = exp(gkernel(3,3)/rt)
      w4 = exp(gkernel(3,4)/rt)
      w5 = exp(gkernel(3,5)/rt)
      
      expected = (-2.d0*w1 - w2 + w4 + 2.d0*w5)
     &     /(w1+w2+w3+w4+w5)
      call assert_real (lmdafmean(3),expected,1.0d-12,
     &                  'buildfkernel from updated gkernel')
c
c     the weights undo the bias on the atoms at the bath temperature,
c     so a lambda particle three times hotter leaves the mean force
c     unchanged, in the separate and in the fused kernel builds, while
c     weighting at that temperature would give a different one
c
      rth = 3.0d0 * rt
      w1 = exp(gkernel(3,1)/rth)
      w2 = exp(gkernel(3,2)/rth)
      w3 = exp(gkernel(3,3)/rth)
      w4 = exp(gkernel(3,4)/rth)
      w5 = exp(gkernel(3,5)/rth)
      hotmean = (-2.0d0*w1 - w2 + w4 + 2.0d0*w5)
     &             / (w1+w2+w3+w4+w5)
      call assert_logical (abs(hotmean-expected).gt.1.0d-3,.true.,
     &                     'buildfkernel hot weights would differ')
      lmdakelvin = 3.0d0 * kelvin
      call buildfkernel
      call assert_real (lmdafmean(3),expected,1.0d-12,
     &                  'buildfkernel mean force at a hot lambda')
      call buildkernels
      call assert_real (lmdafmean(3),expected,1.0d-12,
     &                  'buildkernels mean force at a hot lambda')
      lmdakelvin = 0.0d0
c
c     a deeply filled bin overflows exp(g/kT) unless the largest bias
c     in the row is factored out; unshifted this returns a NaN
c
      call resetost (5,5,1)
      rt = gasconst * kelvin
      gkernel(3,2) = 500.0d0
      gkernel(3,4) = 499.0d0
      vkernelmax(3) = 500.0d0
      call buildfkernel
      wratio = exp((499.0d0-500.0d0)/rt)
      expected = (-1.0d0 + wratio) / (1.0d0 + wratio)
      call assert_real (lmdafmean(3),expected,1.0d-12,
     &                  'buildfkernel large bias no overflow')
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_eost_kernelbuilds  --  build tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_eost_kernelbuilds" checks full kernel rebuilds
c     and incremental updates against reference arrays
c
c
      subroutine test_eost_kernelbuilds
      use bath
      use dlmda
      use math
      use ost
      use units
      implicit none
      integer i,j
      integer ihist
      integer nhist
      real*8 rt
      real*8 flmda
      real*8 partfunc
      real*8 fsum
      real*8 weight
      real*8 vmaxref
      real*8 height
      real*8 fmanual(9)
      real*8 fref(9)
      real*8 fsumref(9)
      real*8 gref(9,9)
      real*8 pfref(9)
c
c
c     build a mixed history with overlapping gaussians, endpoint
c     mirror images, and multiple gaussian widths
c
      call resetost (9,9,8)
      rt = gasconst * kelvin
      nhist = 6
      nlmdahist = nhist
      height = 2.0d0 * pi * wlmda * wflmda
      call sethist (1,0.00d0, 0.0d0,0.7d0*height,wlmda,wflmda)
      call sethist (2,0.25d0,-1.0d0,1.1d0*height,wlmda,wflmda)
      call sethist (3,0.50d0, 0.0d0,1.6d0*height,wlmda,wflmda)
      call sethist (4,0.50d0, 0.0d0,0.4d0*height,wlmda,wflmda)
      call sethist (5,0.75d0, 1.0d0,2.3d0*height,wlmda,wflmda)
      call sethist (6,1.00d0, 0.0d0,0.9d0*height,wlmda,wflmda)
      ostwlhist(5) = 2.0d0 * wlmda
      ostwfhist(5) = 2.0d0 * wflmda
      call buildostindex
c
c     buildgkernel must clear stale values before rebuilding
c
      do i = 1, nlmda
         do j = 1, nflmda
            gkernel(i,j) = -123.0d0
         end do
      end do
      call buildgkernel
      call buildfkernel
      do i = 1, nlmda
         fref(i) = lmdafmean(i)
         do j = 1, nflmda
            gref(i,j) = gkernel(i,j)
         end do
      end do
      call assert_real (gkernel(1,1),0.0d0,1.0d-12,
     &                  'buildgkernel clears outside support')
c
c     compute the buildfkernel result independently from gkernel
c
      do i = 1, nlmda
         partfunc = 0.0d0
         fsum = 0.0d0
c
c     the accumulators factor out the largest bias in the row, so the
c     reference has to be built with the same vkernelmax shift
c
         vmaxref = 0.0d0
         do j = 1, nflmda
            if (gref(i,j) .ne. 0.0d0)
     &         vmaxref = max(vmaxref,gref(i,j))
         end do
         do j = 1, nflmda
            if (gref(i,j) .ne. 0.0d0) then
               flmda = dble(j-fli0) * wflmda
               weight = exp((gref(i,j)-vmaxref)/rt)
               fsum = fsum + flmda*weight
               partfunc = partfunc + weight
            end if
         end do
         if (partfunc .eq. 0.0d0) then
            fmanual(i) = 0.0d0
         else
            fmanual(i) = fsum / partfunc
         end if
         fsumref(i) = fsum
         pfref(i) = partfunc
      end do
      call assert_array1 (lmdafmean,fmanual,nlmda,1.0d-12,
     &                    'buildfkernel all-row weighted mean')
      call assert_array1 (lmdafsum,fsumref,nlmda,1.0d-12,
     &                    'buildfkernel fsum accumulator')
      call assert_array1 (lmdafwt,pfref,nlmda,1.0d-12,
     &                    'buildfkernel partfunc accumulator')
c
c     buildkernels must produce the same gkernel, lmdafmean and
c     free energy accumulators as the old full rebuild path
c
      do i = 1, nlmda
         lmdafmean(i) = -123.0d0
         lmdafsum(i) = -123.0d0
         lmdafwt(i) = -123.0d0
         do j = 1, nflmda
            gkernel(i,j) = -123.0d0
         end do
      end do
      call buildkernels
      call assert_array2 (gkernel,gref,nlmda,nflmda,1.0d-12,
     &                    'buildkernels gkernel reference')
      call assert_array1 (lmdafmean,fref,nlmda,1.0d-12,
     &                    'buildkernels lmdafmean reference')
      call assert_array1 (lmdafsum,fsumref,nlmda,1.0d-12,
     &                    'buildkernels fsum reference')
      call assert_array1 (lmdafwt,pfref,nlmda,1.0d-12,
     &                    'buildkernels partfunc reference')
c
c     incrementally updating kernels one history at a time must match
c     a full buildkernels rebuild
c
      call resetost (9,9,8)
      height = 2.0d0 * pi * wlmda * wflmda
      do ihist = 1, nhist
         nlmdahist = ihist
         if (ihist .eq. 1) then
            call sethist (ihist,0.00d0,0.0d0,0.7d0*height,wlmda,wflmda)
         else if (ihist .eq. 2) then
            call sethist (ihist,0.25d0,-1.0d0,1.1d0*height,wlmda,wflmda)
         else if (ihist .eq. 3) then
            call sethist (ihist,0.50d0,0.0d0,1.6d0*height,wlmda,wflmda)
         else if (ihist .eq. 4) then
            call sethist (ihist,0.50d0,0.0d0,0.4d0*height,wlmda,wflmda)
         else if (ihist .eq. 5) then
            call sethist (ihist,0.75d0,1.0d0,2.3d0*height,wlmda,wflmda)
            ostwlhist(ihist) = 2.0d0 * wlmda
            ostwfhist(ihist) = 2.0d0 * wflmda
         else
            call sethist (ihist,1.00d0,0.0d0,0.9d0*height,wlmda,wflmda)
         end if
         call buildostindex
         call updatekernels
      end do
      call assert_array2 (gkernel,gref,nlmda,nflmda,1.0d-12,
     &                    'updatekernels gkernel reference')
      call assert_array1 (lmdafmean,fref,nlmda,1.0d-12,
     &                    'updatekernels lmdafmean reference')
      call assert_array1 (lmdafsum,fsumref,nlmda,1.0d-12,
     &                    'updatekernels fsum reference')
      call assert_array1 (lmdafwt,pfref,nlmda,1.0d-12,
     &                    'updatekernels partfunc reference')
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_eost_avgstd  --  OST average tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_eost_avgstd" checks interval averages and
c     standard deviations for sampled lambda values
c
c
      subroutine test_eost_avgstd
      use dlmda
      implicit none
      integer i
      real*8 stdref
c
c
c     average all saved interval samples after equilibration prefix
c
      call resetost (5,5,1)
      lmdaintv = 6
      lmdanpa = 1
      lmdanpb = 1
      lmdanpc = 4
      do i = 1, lmdaintv
         lmdallist(i) = dble(i)
         lmdaflist(i) = 2.0d0*dble(i)
      end do
      call avgstd (lmdallist,lmdanpa+lmdanpb+1,lmdanpc,
     &             lmdaavg,lmdastd)
      call avgstd (lmdaflist,lmdanpa+lmdanpb+1,lmdanpc,
     &             dedlavg,dedlstd)
      stdref = sqrt(1.25d0)
      call assert_real (lmdaavg,4.5d0,1.0d-12,
     &                  'avgstd configurable lambda average')
      call assert_real (dedlavg,9.0d0,1.0d-12,
     &                  'avgstd configurable dE/dl average')
      call assert_real (lmdastd,stdref,1.0d-12,
     &                  'avgstd configurable lambda std')
      call assert_real (dedlstd,2.0d0*stdref,1.0d-12,
     &                  'avgstd configurable dE/dl std')
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_eost_eginterpolate  --  g interp tests  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_eost_eginterpolate" checks interpolated gaussian
c     bias values and derivatives
c
c
      subroutine test_eost_eginterpolate
      use dlmda
      use math
      use mutant
      use ost
      implicit none
      real*8 egbias0,dgdl0,dgdfl0
      real*8 egbias1,dgdl1,dgdfl1
      real*8 height
      real*8 sigl,sigf
c
c
c     grid-point interpolation should reproduce the analytic gaussian
c     sum and derivatives from egkernel exactly to roundoff
c
      call resetost (9,9,4)
      sigl = 2.0d0 * wlmda
      sigf = 2.0d0 * wflmda
      height = 2.0d0 * pi * sigl * sigf
      oststdev = 4.0d0
      nlmdahist = 3
      call sethist (1,0.25d0,-1.0d0,1.1d0*height,sigl,sigf)
      call sethist (2,0.50d0, 0.0d0,1.6d0*height,sigl,sigf)
      call sethist (3,0.75d0, 1.0d0,2.3d0*height,sigl,sigf)
      call buildostindex
      call buildkernels
      lambda = 0.50d0
      dedl = 0.0d0
      ostinterpol = .false.
      call egkernel (egbias0,dgdl0,dgdfl0)
      call egkernelinterpolate (egbias1,dgdl1,dgdfl1)
      call assert_real (egbias1,egbias0,1.0d-12,
     &                  'egkernel interpolate grid bias')
      call assert_real (dgdl1,dgdl0,1.0d-12,
     &                  'egkernel interpolate grid dgdl')
      call assert_real (dgdfl1,dgdfl0,1.0d-12,
     &                  'egkernel interpolate grid dgdfl')
c
c     off-grid interpolation is approximate; wide gaussians should
c     keep the value and derivative errors small
c
      call resetost (17,17,4)
      sigl = 4.0d0 * wlmda
      sigf = 4.0d0 * wflmda
      height = 2.0d0 * pi * sigl * sigf
      oststdev = 4.0d0
      nlmdahist = 3
      call sethist (1,0.25d0,-2.0d0,0.8d0*height,sigl,sigf)
      call sethist (2,0.50d0, 0.0d0,1.2d0*height,sigl,sigf)
      call sethist (3,0.75d0, 2.0d0,1.6d0*height,sigl,sigf)
      call buildostindex
      call buildkernels
      lambda = 0.53125d0
      dedl = 0.25d0
      ostinterpol = .false.
      call egkernel (egbias0,dgdl0,dgdfl0)
      call egkernelinterpolate (egbias1,dgdl1,dgdfl1)
      call assert_real (egbias1,egbias0,1.0d-3,
     &                  'egkernel interpolate offgrid bias')
      call assert_real (dgdl1,dgdl0,2.0d-2,
     &                  'egkernel interpolate offgrid dgdl')
      call assert_real (dgdfl1,dgdfl0,2.0d-2,
     &                  'egkernel interpolate offgrid dgdfl')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eost_efkernel  --  free-energy tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eost_efkernel" checks DeltaG interpolation and
c     total free-energy integration from f kernels
c
c
      subroutine test_eost_efkernel
      use dlmda
      use mutant
      implicit none
      integer i
      real*8 eostlmda,dfdl
      real*8 expected
      real*8 efreetot
c
c
c     use lmdafmean(lambda)=lambda and integrate to lambda=0.375
c
      call resetost (5,5,1)
      do i = 1, nlmda
         lmdafmean(i) = dble(i-1) * wlmda
      end do
      lambda = 0.375d0
      call efreelmda (lambda,eostlmda,dfdl)
      expected = 0.5d0 * lambda * lambda
      call assert_real (eostlmda,expected,1.0d-12,
     &                  'efreelmda DeltaG lambda=.375')
      call assert_real (dfdl,lambda,1.0d-12,
     &                  'efreelmda dDeltaG/dlambda')
c
c     use lmdafmean(lambda)=1+lambda so the endpoint mean forces are
c     nonzero and distinct from each other
c
      call resetost (5,5,1)
      do i = 1, nlmda
         lmdafmean(i) = 1.0d0 + dble(i-1)*wlmda
      end do
c
c     at and below lambda = 0 the free energy is zero and the
c     derivative comes from the first lambda bin
c
      lambda = 0.0d0
      call efreelmda (lambda,eostlmda,dfdl)
      call assert_real (eostlmda,0.0d0,1.0d-12,
     &                  'efreelmda DeltaG lambda=0')
      call assert_real (dfdl,1.0d0,1.0d-12,
     &                  'efreelmda dDeltaG/dlambda lambda=0')
      lambda = -0.25d0
      call efreelmda (lambda,eostlmda,dfdl)
      call assert_real (eostlmda,0.0d0,1.0d-12,
     &                  'efreelmda DeltaG lambda below 0')
      call assert_real (dfdl,1.0d0,1.0d-12,
     &                  'efreelmda dDeltaG/dlambda below 0')
c
c     at lambda = 1 the last interval is integrated in full, and
c     beyond lambda = 1 the loop falls through to the same result
c
      lambda = 1.0d0
      call efreelmda (lambda,eostlmda,dfdl)
      call assert_real (eostlmda,1.5d0,1.0d-12,
     &                  'efreelmda DeltaG lambda=1')
      call assert_real (dfdl,2.0d0,1.0d-12,
     &                  'efreelmda dDeltaG/dlambda lambda=1')
      lambda = 1.25d0
      call efreelmda (lambda,eostlmda,dfdl)
      call assert_real (eostlmda,1.5d0,1.0d-12,
     &                  'efreelmda DeltaG lambda above 1')
      call assert_real (dfdl,2.0d0,1.0d-12,
     &                  'efreelmda dDeltaG/dlambda above 1')
c
c     the trapezoid rule is exact for a linear mean force, so the
c     total must match the efreelmda value at lambda = 1
c
      call assert_real (efreetot(),1.5d0,1.0d-12,
     &                  'efreetot linear lmdafmean')
c
c     a constant mean force integrates to itself over unit lambda
c
      do i = 1, nlmda
         lmdafmean(i) = 2.5d0
      end do
      call assert_real (efreetot(),2.5d0,1.0d-12,
     &                  'efreetot constant lmdafmean')
c
c     an all-zero mean force gives no free energy change
c
      do i = 1, nlmda
         lmdafmean(i) = 0.0d0
      end do
      call assert_real (efreetot(),0.0d0,1.0d-12,
     &                  'efreetot zero lmdafmean')
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_eost_meta  --  metadynamics tests  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "test_eost_meta" checks metadynamics bias evaluation,
c     free-energy difference and history resizing
c
c
      subroutine test_eost_meta
      use math
      use ost
      implicit none
      real*8 vbias,dvdl
      real*8 pref
      real*8 expected
      real*8 metadeltag
c
c
c     one normalized 1D gaussian centered at lambda=0.5
c
      call resetost (5,5,1)
      call resetmeta (2)
      nmetahist = 1
      metalhist(1) = 0.5d0
      metahhist(1) = 2.0d0
      metawhist(1) = 0.25d0
      pref = metahhist(1) / (metawhist(1)*sqrt(2.0d0*pi))
c
c     a cutoff of many widths keeps every image term
c
      oststdev = 100.0d0
      maxwlhist = metawhist(1)
      call buildmetaindex
c
c     a gaussian centered midway is symmetric in its two images
c
      call emetabias (0.5d0,vbias,dvdl)
      expected = pref * (1.0d0 + 2.0d0*exp(-8.0d0))
      call assert_real (vbias,expected,1.0d-12,
     &                  'emetabias center value')
      call assert_real (dvdl,0.0d0,1.0d-12,
     &                  'emetabias center derivative')
      call emetabias (0.75d0,vbias,dvdl)
      expected = pref * (exp(-0.5d0)+exp(-12.5d0)+exp(-4.5d0))
      call assert_real (vbias,expected,1.0d-12,
     &                  'emetabias off-center value')
      expected = -16.0d0 * pref
     &              * (0.25d0*exp(-0.5d0) + 1.25d0*exp(-12.5d0)
     &                 - 0.75d0*exp(-4.5d0))
      call assert_real (dvdl,expected,1.0d-12,
     &                  'emetabias off-center derivative')
c
c     symmetric gaussian has zero endpoint free energy difference,
c     read from the end points of the bias grid
c
      call addmetagrid (1)
      call assert_real (vmetagrid(1),pref*(2.0d0*exp(-2.0d0)
     &                  +exp(-18.0d0)),1.0d-12,
     &                  'metadeltag grid end point bias')
      call assert_real (metadeltag(),0.0d0,1.0d-12,
     &                  'metadeltag symmetric gaussian')
c
c     a four width cutoff drops the image five widths away
c
      oststdev = 4.0d0
      call emetabias (0.75d0,vbias,dvdl)
      expected = pref * (exp(-0.5d0)+exp(-4.5d0))
      call assert_real (vbias,expected,1.0d-12,
     &                  'emetabias cutoff value')
      expected = -16.0d0 * pref
     &              * (0.25d0*exp(-0.5d0) - 0.75d0*exp(-4.5d0))
      call assert_real (dvdl,expected,1.0d-12,
     &                  'emetabias cutoff derivative')
c
c     resizing preserves old metadynamics history and zeros new slots
c
      nmetahist = 2
      metalhist(2) = 0.25d0
      metahhist(2) = 3.0d0
      metawhist(2) = 0.125d0
      call resizemeta
      call assert_int (sizemetahist,4,
     &                 'resizemeta size input size=2')
      call assert_real (metalhist(2),0.25d0,1.0d-12,
     &                  'resizemeta preserve metalhist(2)')
      call assert_real (metalhist(3),0.0d0,1.0d-12,
     &                  'resizemeta init metalhist(3)')
      call assert_int (size(metanext),4,
     &                 'resizemeta resizes the bin links')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_eost_depcriteria  --  deposit gate tests  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_eost_depcriteria" checks the deposition gate against a
c     tolerance built from an absolute and a relative part
c
c
      subroutine test_eost_depcriteria
      use dlmda
      implicit none
      logical depcriteria
c
c
c     with the gate off every interval is accepted
c
      use_lmdacv = .false.
      lmdacvstd = 0.0d0
      lmdacvrat = 0.0d0
      call assert_logical (depcriteria(0.0d0,1.0d6),.true.,
     &                     'depcriteria accepts with the gate off')
c
c     the tolerance is lmdacvstd plus lmdacvrat times the average
c
      use_lmdacv = .true.
      lmdacvstd = 10.0d0
      lmdacvrat = 0.2d0
      call assert_logical (depcriteria(0.0d0,9.9d0),.true.,
     &                     'depcriteria inside absolute tolerance')
      call assert_logical (depcriteria(0.0d0,10.0d0),.false.,
     &                     'depcriteria at absolute tolerance')
      call assert_logical (depcriteria(50.0d0,19.9d0),.true.,
     &                     'depcriteria inside relative tolerance')
      call assert_logical (depcriteria(50.0d0,20.0d0),.false.,
     &                     'depcriteria at relative tolerance')
      call assert_logical (depcriteria(-50.0d0,19.9d0),.true.,
     &                     'depcriteria relative tolerance is signless')
      call assert_logical (depcriteria(-50.0d0,20.0d0),.false.,
     &                     'depcriteria negative average at tolerance')
c
c     a vanishing tolerance rejects every interval
c
      lmdacvstd = 0.0d0
      lmdacvrat = 0.0d0
      call assert_logical (depcriteria(1.0d0,0.0d0),.false.,
     &                     'depcriteria zero tolerance')
      use_lmdacv = .false.
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eost_vkernelmax  --  bias level test  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eost_vkernelmax" checks the running maximum of the g
c     kernel along flambda and the global path bias level derived
c     from it, over every path that fills the kernel
c
c
      subroutine test_eost_vkernelmax
      use bath
      use dlmda
      use ost
      implicit none
      integer i
      real*8 vmin
      real*8 brutevkmax
      real*8 ostvminimax
      real*8, allocatable :: gsave(:,:)
c
c
c     two saved sources spread over the whole kernel
c
      kelvin = 300.0d0
      call resetost (5,5,4)
      nlmdahist = 2
      call sethist (1,0.25d0,0.0d0,1.0d0,0.25d0,1.0d0)
      call sethist (2,0.75d0,1.0d0,2.0d0,0.25d0,1.0d0)
      call buildostindex
      call buildkernels
      vmin = brutevkmax (1)
      do i = 1, nlmda
         call assert_real (vkernelmax(i),brutevkmax(i),1.0d-12,
     &                     'buildkernels running maximum')
         vmin = min(vmin,brutevkmax(i))
      end do
      call assert_real (ostvminimax(),vmin,1.0d-12,
     &                  'ostvminimax over lambda')
      call assert_logical (vmin.gt.0.0d0,.true.,
     &                     'ostvminimax reaches every lambda bin')
c
c     the incremental update path stays exact
c
      nlmdahist = 3
      call sethist (3,0.5d0,-1.0d0,1.5d0,0.25d0,1.0d0)
      call buildostindex
      call updatekernels
      do i = 1, nlmda
         call assert_real (vkernelmax(i),brutevkmax(i),1.0d-12,
     &                     'updatekernels running maximum')
      end do
c
c     the g kernel only path reproduces both kernel and maximum
c
      allocate (gsave(nlmda,nflmda))
      gsave = gkernel
      call buildgkernel
      call assert_array2 (gkernel,gsave,nlmda,nflmda,1.0d-12,
     &                    'buildgkernel reproduces gkernel')
      do i = 1, nlmda
         call assert_real (vkernelmax(i),brutevkmax(i),1.0d-12,
     &                     'buildgkernel running maximum')
      end do
      deallocate (gsave)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_eost_fill  --  ost filled bias f term  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_eost_fill" checks that the f term of the ost bias follows
c     the mean force with rows lacking g kernel support filled, while
c     the free energy estimate keeps the sampled rows only, and that
c     the rows a gaussian reaches count as sampled, with the same row
c     weights from both kernel build paths
c
c
      subroutine test_eost_fill
      use atoms
      use dlmda
      use energi
      use mutant
      use ost
      implicit none
      integer i
      integer nsave
      real*8 efreetot
      real*8 wtref(5)
      real*8 sumref(5)
      logical lmdasampled
c
c
c     with no atoms "eostbias" adds no force, so the unallocated
c     force arrays are never read
c
      nsave = n
      n = 0
c
c     give data and a mean force to the two top rows only, with no
c     gaussian at the current point so the g term is zero
c
      call resetost (5,5,4)
      use_ost = .true.
      use_metadyn = .false.
      lmdafwt(4) = 1.0d0
      lmdafwt(5) = 1.0d0
      lmdafmean(4) = -20.0d0
      lmdafmean(5) = -40.0d0
      dedl = 0.0d0
c
c     the f term holds a mean force of -20 below lambda 0.75 with a
c     reach that fills every row
c
      use_lmdafill = .true.
      lmdafillr = 1.0d0
      lambda = 1.0d0
      esum = 0.0d0
      call eostbias
      call assert_real (lmdavbias,22.5d0,1.0d-12,
     &                  'eostbias fill bias energy at lambda 1')
      call assert_real (esum,22.5d0,1.0d-12,
     &                  'eostbias fill energy at lambda 1')
      call assert_real (lmdadfdl,-40.0d0,0.0d0,
     &                  'eostbias fill slope at lambda 1')
      lambda = 0.3d0
      call eostbias
      call assert_real (lmdavbias,6.0d0,1.0d-12,
     &                  'eostbias fill bias energy in a gap')
      call assert_real (lmdadfdl,-20.0d0,0.0d0,
     &                  'eostbias fill slope in a gap')
c
c     the estimate keeps zero mean force in the unsupported rows
c
      call assert_real (efreetot(),-10.0d0,1.0d-12,
     &                  'ost efreetot uses sampled rows only')
      call assert_real (lmdafmean(1),0.0d0,0.0d0,
     &                  'ost fill leaves the empty row mean force')
c
c     a reach of one row fills only the row below lambda 0.75
c
      lmdafillr = 0.25d0
      call eostbias
      call assert_real (lmdavbias,0.1d0,1.0d-12,
     &                  'eostbias reach bias energy past the ramp')
      call assert_real (lmdadfdl,-4.0d0,1.0d-12,
     &                  'eostbias reach slope past the ramp')
c
c     without the fill the f term is the sampled free energy
c
      use_lmdafill = .false.
      call eostbias
      call assert_real (lmdavbias,0.0d0,0.0d0,
     &                  'eostbias nofill bias energy in a gap')
      call assert_real (lmdadfdl,0.0d0,0.0d0,
     &                  'eostbias nofill slope in a gap')
c
c     a gaussian at lambda 1 reaches only the two top rows
c
      call resetost (5,5,4)
      use_ost = .true.
      nlmdahist = 1
      call sethist (1,1.0d0,0.0d0,1.0d0,0.3d0,1.0d0)
      call buildostindex
      call buildkernels
      do i = 1, nlmda
         call assert_logical (lmdasampled(i),i.ge.4,
     &                        'buildkernels rows with support')
      end do
      call assert_real (lmdafmean(3),0.0d0,0.0d0,
     &                  'buildkernels empty row mean force')
c
c     the g kernel and f kernel path sets the same row weights
c
      do i = 1, nlmda
         wtref(i) = lmdafwt(i)
         lmdafwt(i) = -123.0d0
      end do
      call buildgkernel
      call buildfkernel
      call assert_array1 (lmdafwt,wtref,nlmda,1.0d-12,
     &                    'buildfkernel row weights')
      do i = 1, nlmda
         call assert_logical (lmdasampled(i),i.ge.4,
     &                        'buildfkernel rows with support')
      end do
c
c     a gaussian of zero height, as a fully tempered deposit gives,
c     leaves every row without data
c
      call resetost (5,5,4)
      use_ost = .true.
      nlmdahist = 1
      call sethist (1,1.0d0,0.0d0,0.0d0,0.3d0,1.0d0)
      call buildostindex
      call buildkernels
      do i = 1, nlmda
         call assert_logical (lmdasampled(i),.false.,
     &                        'zero height gaussian rows without data')
         call assert_real (lmdafsum(i),0.0d0,0.0d0,
     &                     'zero height gaussian row force sum')
         call assert_real (lmdafmean(i),0.0d0,0.0d0,
     &                     'zero height gaussian row mean force')
      end do
c
c     beside a gaussian of positive height, one of zero height that
c     reaches other cells of the same rows adds no weight, so both
c     kernel build paths agree on the rows
c
      nlmdahist = 2
      call sethist (1,1.0d0,0.0d0,1.0d0,0.3d0,1.0d0)
      call sethist (2,1.0d0,1.0d0,0.0d0,0.3d0,1.0d0)
      call buildostindex
      call buildkernels
      do i = 1, nlmda
         wtref(i) = lmdafwt(i)
         sumref(i) = lmdafsum(i)
      end do
      call buildgkernel
      call buildfkernel
      call assert_array1 (lmdafwt,wtref,nlmda,1.0d-12,
     &                    'zero height gaussian row weights')
      call assert_array1 (lmdafsum,sumref,nlmda,1.0d-12,
     &                    'zero height gaussian row force sums')
      do i = 1, nlmda
         call assert_logical (lmdasampled(i),i.ge.4,
     &                        'zero height gaussian rows with data')
      end do
      use_ost = .false.
      n = nsave
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function brutevkmax  --  direct kernel column maximum  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "brutevkmax" scans the g kernel directly for the largest
c     value at one lambda bin, ignoring the running maximum
c
c
      function brutevkmax (ilmda)
      use ost
      implicit none
      integer ilmda,j
      real*8 brutevkmax
c
c
      brutevkmax = 0.0d0
      do j = 1, nflmda
         brutevkmax = max(brutevkmax,gkernel(ilmda,j))
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_eost_tempering  --  height decay test  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_eost_tempering" checks the deposited gaussian height with
c     the global and local tempering factors disabled, below their
c     thresholds, active alone and active together, and that the
c     height stays positive and never exceeds the full height
c
c
      subroutine test_eost_tempering
      use bath
      use ost
      use units
      implicit none
      integer i,j,k
      real*8 rt
      real*8 h,prev
      real*8 h1,h2
      real*8 hg,hl
      real*8 vg,dl
      real*8 temperedheight
      real*8 vstar(4)
      real*8 delta(5)
      real*8 vgrid(6)
      real*8 dgrid(6)
      real*8 gpair(2,3)
      real*8 hgrid(6,6)
      data vstar / 1.5d0,2.0d0,3.0d0,5.0d0 /
      data delta / 0.0d0,1.0d0,1.5d0,3.0d0,5.0d0 /
      data vgrid / 0.0d0,0.5d0,1.0d0,2.0d0,5.0d0,20.0d0 /
      data dgrid / 0.0d0,0.5d0,1.0d0,2.0d0,5.0d0,20.0d0 /
      data gpair / 1.0d0,1.0d0,2.0d0,0.5d0,0.5d0,2.0d0 /
c
c
c     an untempered run deposits at the full height
c
      kelvin = 300.0d0
      call resetost (5,5,1)
      hbias = 1.0d-5
      rt = gasconst * kelvin
      use_ostgtemp = .false.
      use_ostltemp = .false.
      ostgthresh = 1.0d0
      ostgtempgamma = 1.0d0
      ostlthresh = 1.0d0
      ostltempgamma = 1.0d0
      call assert_real (temperedheight(0.0d0,0.0d0),hbias,1.0d-18,
     &                  'temperedheight disabled empty path')
      call assert_real (temperedheight(50.0d0,50.0d0),hbias,1.0d-18,
     &                  'temperedheight disabled filled path')
      call assert_real (temperedheight(50.0d0,500.0d0),hbias,1.0d-18,
     &                  'temperedheight disabled uneven path')
c
c     at or below both thresholds the height is still untempered
c
      use_ostgtemp = .true.
      use_ostltemp = .true.
      call assert_real (temperedheight(0.0d0,0.0d0),hbias,1.0d-18,
     &                  'temperedheight below both thresholds')
      call assert_real (temperedheight(1.0d0,1.0d0),hbias,1.0d-18,
     &                  'temperedheight at global threshold')
      call assert_real (temperedheight(1.0d0,2.0d0),hbias,1.0d-18,
     &                  'temperedheight at both thresholds')
      call assert_real (temperedheight(0.5d0,1.5d0),hbias,1.0d-18,
     &                  'temperedheight at local threshold')
c
c     the global factor alone decays with the path bias level and
c     ignores how far the deposit bin is ahead of that level
c
      use_ostltemp = .false.
      prev = hbias
      do i = 1, 4
         h = temperedheight(vstar(i),vstar(i))
         call assert_real (h,hbias*exp(-(vstar(i)-1.0d0)/rt),1.0d-18,
     &                     'temperedheight global above threshold')
         call assert_real (temperedheight(vstar(i),vstar(i)+50.0d0),
     &                     h,1.0d-18,'temperedheight global ignores '//
     &                     'the local excess')
         call assert_logical (h.lt.prev,.true.,
     &                        'temperedheight global decays')
         prev = h
      end do
c
c     a larger global tempering factor decays more slowly
c
      ostgtempgamma = 1.0d0
      h1 = temperedheight(3.0d0,3.0d0)
      ostgtempgamma = 2.0d0
      h2 = temperedheight(3.0d0,3.0d0)
      call assert_logical (h2.gt.h1,.true.,
     &                     'temperedheight larger gamma decays slower')
      call assert_real (h2,hbias*exp(-2.0d0/(2.0d0*rt)),1.0d-18,
     &                  'temperedheight gamma scaling')
c
c     the local factor alone depends only on the excess of the
c     deposit bin over the path bias level
c
      use_ostgtemp = .false.
      use_ostltemp = .true.
      ostgtempgamma = 1.0d0
      prev = hbias
      do i = 1, 5
         h = temperedheight(50.0d0,50.0d0+delta(i))
         hl = hbias * exp(-max(0.0d0,delta(i)-1.0d0)/rt)
         call assert_real (h,hl,1.0d-18,
     &                     'temperedheight local excess')
         call assert_real (temperedheight(0.0d0,delta(i)),h,1.0d-18,
     &                     'temperedheight local shift invariance')
         call assert_logical (h.le.prev,.true.,
     &                        'temperedheight local decays')
         prev = h
      end do
      call assert_real (temperedheight(50.0d0,50.0d0),hbias,1.0d-18,
     &                  'temperedheight least filled bin untempered')
c
c     both factors multiply, so with equal settings the path bias
c     level cancels and only the deposit bin bias level remains
c
      use_ostgtemp = .true.
      use_ostltemp = .true.
      ostgthresh = 1.0d0
      ostlthresh = 1.0d0
      ostgtempgamma = 2.0d0
      ostltempgamma = 2.0d0
      h = temperedheight(3.0d0,6.0d0)
      hg = exp(-2.0d0/(2.0d0*rt))
      hl = exp(-2.0d0/(2.0d0*rt))
      call assert_real (h,hbias*hg*hl,1.0d-18,
     &                  'temperedheight factors multiply')
      call assert_real (h,hbias*exp((2.0d0-6.0d0)/(2.0d0*rt)),1.0d-18,
     &                  'temperedheight equal settings')
      ostltempgamma = 0.5d0
      hg = exp(-2.0d0/(2.0d0*rt))
      hl = exp(-2.0d0/(0.5d0*rt))
      call assert_real (temperedheight(3.0d0,6.0d0),hbias*hg*hl,
     &                  1.0d-18,'temperedheight unequal gammas')
c
c     over a sweep of levels the height stays positive, never exceeds
c     the full height, and never grows with either level
c
      ostgthresh = 1.0d0
      ostlthresh = 1.0d0
      do k = 1, 3
         ostgtempgamma = gpair(1,k)
         ostltempgamma = gpair(2,k)
         do i = 1, 6
            do j = 1, 6
               vg = vgrid(i)
               dl = dgrid(j)
               h = temperedheight(vg,vg+dl)
               hgrid(i,j) = h
               call assert_logical (h.gt.0.0d0 .and. h.le.hbias,.true.,
     &                              'temperedheight bounded height')
               if (i .gt. 1) then
                  call assert_logical (h.le.hgrid(i-1,j),.true.,
     &                        'temperedheight nonincreasing in path')
               end if
               if (j .gt. 1) then
                  call assert_logical (h.le.hgrid(i,j-1),.true.,
     &                        'temperedheight nonincreasing in bin')
               end if
            end do
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eost_metaimage  --  meta image tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eost_metaimage" checks that the metadynamics bias sums
c     the same three reflected lambda images as the g kernel
c
c
      subroutine test_eost_metaimage
      use math
      use ost
      implicit none
      integer m,t
      real*8 vbias,dvdl
      real*8 vref,dref
      real*8 delta,bias
      real*8 pref,sig2
      real*8 single,far
      real*8 src(3)
      real*8 lam(4)
      data lam / 0.0d0,0.1d0,0.5d0,1.0d0 /
c
c
c     one gaussian near the lambda = 0 wall
c
      call resetost (5,5,1)
      call resetmeta (2)
      nmetahist = 1
      metalhist(1) = 0.1d0
      metahhist(1) = 2.0d0
      metawhist(1) = 0.25d0
      oststdev = 100.0d0
      maxwlhist = metawhist(1)
      call buildmetaindex
      pref = metahhist(1) / (metawhist(1)*sqrt(2.0d0*pi))
      sig2 = metawhist(1) * metawhist(1)
      src(1) = 0.1d0
      src(2) = -0.1d0
      src(3) = 1.9d0
c
c     the bias and derivative match a direct three image sum
c
      do t = 1, 4
         vref = 0.0d0
         dref = 0.0d0
         do m = 1, 3
            delta = lam(t) - src(m)
            bias = pref * exp(-0.5d0*delta*delta/sig2)
            vref = vref + bias
            dref = dref - delta*bias/sig2
         end do
         call emetabias (lam(t),vbias,dvdl)
         call assert_real (vbias,vref,1.0d-12,
     &                     'emetabias image sum value')
         call assert_real (dvdl,dref,1.0d-12,
     &                     'emetabias image sum derivative')
      end do
c
c     at the wall the nearby image doubles the bias
c
      single = pref * exp(-0.5d0*0.01d0/sig2)
      far = pref * exp(-0.5d0*1.9d0*1.9d0/sig2)
      call emetabias (0.0d0,vbias,dvdl)
      call assert_real (vbias,2.0d0*single+far,1.0d-12,
     &                  'emetabias doubles at the lambda wall')
      call assert_logical (vbias.gt.single,.true.,
     &                     'emetabias wall raises the bias')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eost_metadyn  --  meta deposit tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eost_metadyn" drives two deposit intervals and checks
c     the gaussian that elmdadyn stores, and that metadynamics ignores
c     the convergence gate and moves lambda on every step
c
c
      subroutine test_eost_metadyn
      use dlmda
      use mutant
      use ost
      implicit none
      integer istep
      real*8 avgref
      real*8 lam(4)
      real*8 lmv(4)
      data lam / 0.1d0,0.2d0,0.4d0,0.6d0 /
c
c
c     a frozen lambda particle keeps the sampled values controlled
c
      call resetost (5,5,1)
      use_meta = .true.
      call resetmeta (2)
      lmdaintv = 4
      lmdanpa = 1
      lmdanpb = 1
      lmdanpc = 2
      hbias = 2.0d0
      wlmda = 0.25d0
      wlhist = 0.3d0
      maxwlhist = wlhist
      dedl = 0.0d0
      lmdadt = 0.0d0
      lmdastep = 0
      do istep = 1, lmdaintv
         lambda = lam(istep)
         call elmdadyn
         if (istep .lt. lmdaintv) then
            call assert_int (nmetahist,0,
     &                       'elmdadyn meta waits for the interval end')
         end if
      end do
c
c     only the samples after the equilibration prefix are averaged
c
      avgref = (lam(3)+lam(4)) / dble(lmdanpc)
      call assert_int (nmetahist,1,
     &                 'elmdadyn meta deposits one gaussian')
      call assert_real (metalhist(1),avgref,1.0d-12,
     &                  'elmdadyn meta gaussian center')
      call assert_real (metahhist(1),hbias,1.0d-12,
     &                  'elmdadyn meta gaussian height')
      call assert_real (metawhist(1),wlhist,1.0d-12,
     &                  'elmdadyn meta gaussian width')
      call assert_int (metaihist(1),lmdaintv,
     &                 'elmdadyn meta gaussian step stamp')
c
c     meta ignores the convergence gate and moves a deterministic
c     frictionless lambda particle on every step of the interval
c
      use_lmdacv = .true.
      lmdacvstd = 1.0d0
      lmdacvrat = 0.0d0
      lmdadt = 0.1d0
      lmdamass = 1.0d0
      lmdafric = 0.0d0
      lmdatheta = 0.25d0 * 3.14159265358979323846d0
      lmdavtheta = 0.0d0
      do istep = 1, lmdaintv
         dedl = 1.0d0
         if (mod(istep,2) .eq. 0)  dedl = 11.0d0
         call elmdadyn
         lmv(istep) = lambda
      end do
      call assert_int (nmetahist,2,
     &                 'elmdadyn meta deposits a noisy interval')
      do istep = 2, lmdaintv
         call assert_logical (lmv(istep).ne.lmv(istep-1),.true.,
     &                        'elmdadyn meta moves lambda each step')
      end do
      use_lmdacv = .false.
      use_meta = .false.
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_eost_metatemper  --  meta decay tests  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_eost_metatemper" drives several deposit intervals with
c     tempering enabled and checks the deposited heights and widths
c     and the accumulated metadynamics grid; a local tempering request
c     must not change the heights, which follow the global level only
c
c
      subroutine test_eost_metatemper
      use bath
      use dlmda
      use mutant
      use ost
      implicit none
      integer il,istep,k
      integer ndep
      parameter (ndep=5)
      real*8 lmda
      real*8 vd,dd,vi,di
      real*8 refvstar
      real*8 temperedheight
c
c
c     deposit repeatedly at a fixed lambda with tempering on
c
      kelvin = 300.0d0
      call resetost (5,5,1)
      use_meta = .true.
      call resetmeta (8)
      lmdaintv = 4
      lmdanpa = 1
      lmdanpb = 1
      lmdanpc = 2
      hbias = 2.0d0
      dedl = 0.0d0
      lmdadt = 0.0d0
      use_ostgtemp = .true.
      ostgthresh = 0.5d0
      ostgtempgamma = 1.0d0
      wlhist = 0.3d0
      maxwlhist = wlhist
      oststdev = 100.0d0
c
c     a low local threshold would lower the heights if it were used
c
      use_ostltemp = .true.
      ostlthresh = 0.01d0
      ostltempgamma = 1.0d0
      lmdastep = 0
      do istep = 1, ndep*lmdaintv
         lambda = 0.5d0
         call elmdadyn
      end do
      call assert_int (nmetahist,ndep,'elmdadyn meta deposit count')
c
c     the first deposit sees an empty bias, so it is untempered
c
      call assert_real (metahhist(1),hbias,1.0d-12,
     &                  'elmdadyn meta first height untempered')
      call assert_logical (refvstar(1).gt.ostgthresh,.true.,
     &                     'elmdadyn meta crosses the threshold')
c
c     every later height follows the pre-deposit bias level
c
      do k = 2, ndep
         call assert_real (metahhist(k),
     &                     temperedheight(refvstar(k-1),refvstar(k-1)),
     &                     1.0d-12,'elmdadyn meta tempered height')
         call assert_logical (metahhist(k).lt.metahhist(k-1),.true.,
     &                        'elmdadyn meta heights decay')
         call assert_real (metawhist(k),wlhist,1.0d-12,
     &                     'elmdadyn meta width from wlhist')
      end do
      use_ostltemp = .false.
c
c     the grid matches the direct sum at each bin center, and the
c     interpolation reproduces it exactly at those nodes
c
      do il = 1, nlmda
         lmda = dble(il-1) * wlmda
         ostinterpol = .false.
         call emetabias (lmda,vd,dd)
         call assert_real (vmetagrid(il),vd,1.0d-12,
     &                     'addmetagrid matches the direct sum')
         call assert_real (dvmetagrid(il),dd,1.0d-12,
     &                     'addmetagrid matches the direct slope')
         ostinterpol = .true.
         call emetabias (lmda,vi,di)
         call assert_real (vi,vd,1.0d-12,
     &                     'emetabiasinterpolate at a node')
         call assert_real (di,dd,1.0d-12,
     &                     'emetabiasinterpolate slope at a node')
      end do
      ostinterpol = .false.
      use_meta = .false.
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function refvstar  --  direct metadynamics bias level  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "refvstar" sums the first "upto" saved metadynamics gaussians
c     directly and returns the smallest bias over the lambda bins
c
c
      function refvstar (upto)
      use dlmda
      use math
      use ost
      implicit none
      integer upto
      integer il,j,m
      real*8 refvstar
      real*8 lambda,delta
      real*8 sig,sig2,pref
      real*8 v
      real*8 src(3)
c
c
      refvstar = 0.0d0
      do il = 1, nlmda
         lambda = dble(il-1) * wlmda
         v = 0.0d0
         do j = 1, upto
            sig = metawhist(j)
            sig2 = sig * sig
            pref = metahhist(j) / (sig*sqrt(2.0d0*pi))
            src(1) = metalhist(j)
            src(2) = -metalhist(j)
            src(3) = 2.0d0 - metalhist(j)
            do m = 1, 3
               delta = lambda - src(m)
               v = v + pref*exp(-0.5d0*delta*delta/sig2)
            end do
         end do
         if (il .eq. 1) then
            refvstar = v
         else
            refvstar = min(refvstar,v)
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_eost_metacut  --  metadynamics bin lookup  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_eost_metacut" spreads gaussians of several widths over the
c     lambda range, including the end bins, and checks that the bias
c     found through the bin lookup, the grid built bin by bin and the
c     bias of a partly replayed history all match a direct sum over
c     every gaussian and image within the cutoff
c
c
      subroutine test_eost_metacut
      use dlmda
      use ost
      implicit none
      integer i,t
      integer il,ng
      real*8 lmda
      real*8 vbias,dvdl
      real*8 vref,dref
      real*8 dvmax,ddmax
      real*8 cen(6),hgt(6),wid(6)
      data cen / 0.0d0,0.02d0,0.33d0,0.5d0,0.97d0,1.0d0 /
      data hgt / 1.0d0,2.0d0,0.5d0,3.0d0,1.5d0,0.7d0 /
      data wid / 0.05d0,0.1d0,0.2d0,0.08d0,0.15d0,0.03d0 /
c
c
c     six gaussians on a grid of 21 lambda bins with a cutoff of three
c     widths, which reaches across several bins for the widest
c
      call resetost (21,5,1)
      ng = 6
      call resetmeta (ng)
      nmetahist = ng
      maxwlhist = 0.0d0
      do i = 1, ng
         metalhist(i) = cen(i)
         metahhist(i) = hgt(i)
         metawhist(i) = wid(i)
         metaihist(i) = i
         maxwlhist = max(maxwlhist,wid(i))
      end do
      oststdev = 3.0d0
      ostinterpol = .false.
      call buildmetaindex
c
c     the lookup sum matches the direct sum across lambda
c
      dvmax = 0.0d0
      ddmax = 0.0d0
      do t = 0, 40
         lmda = 0.025d0 * dble(t)
         call emetabias (lmda,vbias,dvdl)
         call brutemeta (lmda,ng,vref,dref)
         dvmax = max(dvmax,abs(vbias-vref))
         ddmax = max(ddmax,abs(dvdl-dref))
      end do
      call assert_real (dvmax,0.0d0,1.0d-11,
     &                  'emetabias lookup matches direct value')
      call assert_real (ddmax,0.0d0,1.0d-10,
     &                  'emetabias lookup matches direct slope')
c
c     the grid built within the cutoff matches the direct sum
c
      do il = 1, nlmda
         vmetagrid(il) = 0.0d0
         dvmetagrid(il) = 0.0d0
      end do
      do i = 1, ng
         call addmetagrid (i)
      end do
      dvmax = 0.0d0
      ddmax = 0.0d0
      do il = 1, nlmda
         lmda = dble(il-1) * wlmda
         call brutemeta (lmda,ng,vref,dref)
         dvmax = max(dvmax,abs(vmetagrid(il)-vref))
         ddmax = max(ddmax,abs(dvmetagrid(il)-dref))
      end do
      call assert_real (dvmax,0.0d0,1.0d-11,
     &                  'addmetagrid cutoff matches direct value')
      call assert_real (ddmax,0.0d0,1.0d-10,
     &                  'addmetagrid cutoff matches direct slope')
c
c     a partly replayed history counts only its first gaussians
c
      nmetahist = 3
      dvmax = 0.0d0
      do t = 0, 40
         lmda = 0.025d0 * dble(t)
         call emetabias (lmda,vbias,dvdl)
         call brutemeta (lmda,3,vref,dref)
         dvmax = max(dvmax,abs(vbias-vref))
      end do
      call assert_real (dvmax,0.0d0,1.0d-11,
     &                  'emetabias partial history value')
      nmetahist = ng
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine brutemeta  --  direct metadynamics cutoff  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "brutemeta" sums the first "upto" saved metadynamics gaussians
c     and their two images directly at one lambda, keeping only the
c     terms within the cutoff of "oststdev" widths
c
c
      subroutine brutemeta (lmda,upto,v,d)
      use math
      use ost
      implicit none
      integer upto
      integer j,m
      real*8 lmda,v,d
      real*8 delta,bias
      real*8 sig,sig2,pref
      real*8 src(3)
c
c
      v = 0.0d0
      d = 0.0d0
      do j = 1, upto
         sig = metawhist(j)
         sig2 = sig * sig
         pref = metahhist(j) / (sig*sqrt(2.0d0*pi))
         src(1) = metalhist(j)
         src(2) = -metalhist(j)
         src(3) = 2.0d0 - metalhist(j)
         do m = 1, 3
            delta = lmda - src(m)
            if (abs(delta) .le. oststdev*sig) then
               bias = pref * exp(-0.5d0*delta*delta/sig2)
               v = v + bias
               d = d - delta*bias/sig2
            end if
         end do
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_eost_ostdyn  --  ost deposit tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_eost_ostdyn" drives elmdadyn over three deposit intervals
c     and checks that a settled interval deposits, an unsettled one is
c     rejected, and only the averaging phase samples are averaged
c
c
      subroutine test_eost_ostdyn
      use bath
      use dlmda
      use mutant
      use ost
      implicit none
      integer istep
      real*8 eostsave
c
c
c     a settled interval deposits one gaussian at the interval end
c
      kelvin = 300.0d0
      call resetost (5,5,4)
      use_ost = .true.
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      use_lmdacv = .true.
      lmdacvstd = 1.0d0
      lmdacvrat = 0.0d0
      hbias = 1.0d0
      lmdadt = 0.0d0
      fastkernel = .true.
      d2edl2 = 0.0d0
      ostbdgdl = 0.0d0
      ostbdgdfl = 0.0d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      do istep = 1, lmdaintv
         lambda = 0.5d0
         dedl = 1.0d0
         call elmdadyn
         if (istep .lt. lmdaintv) then
            call assert_int (nlmdahist,0,
     &                       'elmdadyn ost waits for the interval end')
         end if
      end do
      call assert_int (nlmdahist,1,
     &                 'elmdadyn ost deposits a settled interval')
      call assert_int (lmdaihist(1),lmdaintv,
     &                 'elmdadyn ost stamps the step')
      call assert_real (lmdalhist(1),0.5d0,1.0d-12,
     &                  'elmdadyn ost gaussian lambda center')
      call assert_real (lmdafhist(1),1.0d0,1.0d-12,
     &                  'elmdadyn ost gaussian flambda center')
      call assert_real (osthhist(1),hbias,1.0d-12,
     &                  'elmdadyn ost untempered gaussian height')
      call assert_real (ostwlhist(1),wlhist,1.0d-12,
     &                  'elmdadyn ost gaussian lambda width')
c
c     an unsettled interval is rejected and changes nothing
c
      eostsave = lmdadeltag
      do istep = 1, lmdaintv
         lambda = 0.5d0
         dedl = 1.0d0
         if (mod(istep,2) .eq. 0)  dedl = 11.0d0
         call elmdadyn
      end do
      call assert_real (dedlavg,6.0d0,1.0d-12,
     &                  'elmdadyn ost unsettled interval average')
      call assert_real (dedlstd,5.0d0,1.0d-12,
     &                  'elmdadyn ost unsettled interval deviation')
      call assert_int (nlmdahist,1,
     &                 'elmdadyn ost rejects an unsettled interval')
      call assert_real (lmdadeltag,eostsave,1.0d-12,
     &                  'elmdadyn ost rejection leaves the free energy')
c
c     only the samples after the equilibration prefix are averaged
c
      lmdanpa = 1
      lmdanpb = 1
      lmdanpc = 2
      do istep = 1, lmdaintv
         lambda = 0.1d0
         dedl = 50.0d0
         if (istep .gt. lmdanpa+lmdanpb) then
            lambda = 0.5d0
            dedl = 2.0d0
         end if
         call elmdadyn
      end do
      call assert_int (nlmdahist,2,
     &                 'elmdadyn ost deposits the average phase')
      call assert_real (lmdalhist(2),0.5d0,1.0d-12,
     &                  'elmdadyn ost skips the prefix lambda')
      call assert_real (lmdafhist(2),2.0d0,1.0d-12,
     &                  'elmdadyn ost skips the prefix dU/dlambda')
      use_lmdacv = .false.
      use_ost = .false.
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eost_ostlocal  --  local tempering test  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eost_ostlocal" deposits into an unevenly filled kernel
c     with both tempering factors on, and checks that elmdadyn takes
c     the height from the pre-deposit bias levels, that the least
c     filled lambda bin deposits at the global height alone, and
c     that growing the flambda grid keeps the bin bias levels
c
c
      subroutine test_eost_ostlocal
      use bath
      use dlmda
      use mutant
      use ost
      use units
      implicit none
      integer i,istep
      integer imax,imin
      real*8 rt,gmin,gl
      real*8 hglobal
      real*8 brutevkmax
      real*8 ostvminimax
      real*8 temperedheight
c
c
c     fill the kernel much higher near lambda of zero
c
      kelvin = 300.0d0
      rt = gasconst * kelvin
      call resetost (5,5,8)
      use_ost = .true.
      oststdev = 4.0d0
      nlmdahist = 2
      call sethist (1,0.0d0,0.0d0,5.0d0,0.25d0,1.0d0)
      call sethist (2,0.75d0,0.0d0,1.0d0,0.25d0,1.0d0)
      call buildostindex
      call buildkernels
c
c     settle each deposit interval with both tempering factors on
c
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      use_lmdacv = .true.
      lmdacvstd = 1.0d0
      lmdacvrat = 0.0d0
      hbias = 1.0d0
      lmdadt = 0.0d0
      fastkernel = .true.
      d2edl2 = 0.0d0
      ostbdgdl = 0.0d0
      ostbdgdfl = 0.0d0
      lmdadfdl = 0.0d0
      use_ostgtemp = .true.
      use_ostltemp = .true.
      ostgthresh = 0.1d0
      ostlthresh = 0.1d0
      ostgtempgamma = 1.0d0
      ostltempgamma = 1.0d0
c
c     a deposit in the most filled bin is tempered by both factors
c
      imax = 1
      do i = 2, nlmda
         if (brutevkmax(i) .gt. brutevkmax(imax))  imax = i
      end do
      gmin = ostvminimax ()
      gl = vkernelmax(imax)
      call assert_logical (gmin.gt.ostgthresh,.true.,
     &                     'elmdadyn ost global factor active')
      call assert_logical (gl-gmin.gt.ostlthresh,.true.,
     &                     'elmdadyn ost local factor active')
      lmdastep = 0
      do istep = 1, lmdaintv
         lambda = dble(imax-1) * wlmda
         dedl = 0.0d0
         call elmdadyn
      end do
      hglobal = hbias * exp(-(gmin-ostgthresh)/rt)
      call assert_int (nlmdahist,3,
     &                 'elmdadyn ost deposits in the full bin')
      call assert_real (osthhist(3),temperedheight(gmin,gl),1.0d-12,
     &                  'elmdadyn ost height from pre-deposit levels')
      call assert_logical (osthhist(3).lt.hglobal,.true.,
     &                     'elmdadyn ost full bin below global height')
c
c     a deposit in the least filled bin sees no local excess
c
      imin = 1
      do i = 2, nlmda
         if (vkernelmax(i) .lt. vkernelmax(imin))  imin = i
      end do
      gmin = ostvminimax ()
      do istep = 1, lmdaintv
         lambda = dble(imin-1) * wlmda
         dedl = 0.0d0
         call elmdadyn
      end do
      hglobal = hbias * exp(-max(0.0d0,gmin-ostgthresh)/rt)
      call assert_int (nlmdahist,4,
     &                 'elmdadyn ost deposits in the least bin')
      call assert_real (osthhist(4),hglobal,1.0d-12,
     &                  'elmdadyn ost least filled bin global height')
c
c     growing the flambda grid keeps the running bin maxima
c
      call ensureflambda (500.0d0)
      do i = 1, nlmda
         call assert_real (vkernelmax(i),brutevkmax(i),1.0d-12,
     &                     'ensureflambda keeps bin bias levels')
      end do
      use_lmdacv = .false.
      use_ost = .false.
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eost_temperkeys  --  tempering keywords  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eost_temperkeys" checks that the tempering keywords set
c     their flags, thresholds and tempering factors, fill a missing
c     value with its default, replace a zero by the default and a
c     negative value by its magnitude
c
c
      subroutine test_eost_temperkeys
      use dlmda
      use keys
      use mutant
      use ost
      implicit none
c
c
c     both keywords set their flags and both values
c
      use_ost = .false.
      use_meta = .false.
      lambda = 0.5d0
      if (allocated(keyline))  deallocate (keyline)
      allocate (keyline(2))
      nkey = 2
      keyline(1) = 'OST-TEMPER-GLOBAL 2.0 3.0'
      keyline(2) = 'ost-temper-local 0.5 0.25'
      call mutate_ost
      call assert_logical (use_ostgtemp,.true.,
     &                     'OST-TEMPER-GLOBAL sets its flag')
      call assert_logical (use_ostltemp,.true.,
     &                     'OST-TEMPER-LOCAL sets its flag')
      call assert_real (ostgthresh,2.0d0,1.0d-12,
     &                  'OST-TEMPER-GLOBAL threshold')
      call assert_real (ostgtempgamma,3.0d0,1.0d-12,
     &                  'OST-TEMPER-GLOBAL tempering factor')
      call assert_real (ostlthresh,0.5d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL threshold')
      call assert_real (ostltempgamma,0.25d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL tempering factor')
c
c     a missing tempering factor keeps its default
c
      nkey = 1
      keyline(1) = 'OST-TEMPER-LOCAL 0.5'
      call mutate_ost
      call assert_logical (use_ostgtemp,.false.,
     &                     'OST-TEMPER-LOCAL leaves global off')
      call assert_logical (use_ostltemp,.true.,
     &                     'OST-TEMPER-LOCAL one value sets its flag')
      call assert_real (ostlthresh,0.5d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL one value threshold')
      call assert_real (ostltempgamma,1.0d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL default tempering factor')
c
c     zero values take the defaults
c
      keyline(1) = 'OST-TEMPER-GLOBAL 0 0'
      call mutate_ost
      call assert_real (ostgthresh,1.0d0,1.0d-12,
     &                  'OST-TEMPER-GLOBAL zero threshold')
      call assert_real (ostgtempgamma,1.0d0,1.0d-12,
     &                  'OST-TEMPER-GLOBAL zero tempering factor')
c
c     negative values take their magnitude
c
      keyline(1) = 'OST-TEMPER-LOCAL -2.0 -0.5'
      call mutate_ost
      call assert_real (ostlthresh,2.0d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL negative threshold')
      call assert_real (ostltempgamma,0.5d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL negative tempering factor')
c
c     no keywords and the removed keyword leave tempering off
c
      nkey = 0
      call mutate_ost
      call assert_logical (use_ostgtemp .or. use_ostltemp,.false.,
     &                     'tempering off without keywords')
      nkey = 1
      keyline(1) = 'OST-TEMPER'
      call mutate_ost
      call assert_logical (use_ostgtemp,.false.,
     &                     'OST-TEMPER no longer enables tempering')
c
c     restore an empty keyword list and a clean OST state
c
      nkey = 0
      deallocate (keyline)
      call resetost (5,5,1)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eost_ostphase  --  phase split tests     ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eost_ostphase" checks that setlmdaphase divides the deposit
c     interval into propagation, equilibration and averaging phases,
c     and that the clamps keep a propagation step and enough samples
c     to average without ever losing a sample from the interval
c
c
      subroutine test_eost_ostphase
      use dlmda
      implicit none
c
c
c     the requested ratios divide the interval by truncation
c
      lmdaintv = 10
      lmdaparatio = 0.3d0
      lmdapbratio = 0.3d0
      call setlmdaphase
      call assert_int (lmdanpa,3,'setlmdaphase propagation phase')
      call assert_int (lmdanpb,3,'setlmdaphase equilibration phase')
      call assert_int (lmdanpc,4,'setlmdaphase averaging phase')
      call assert_int (lmdanpa+lmdanpb+lmdanpc,lmdaintv,
     &                 'setlmdaphase spans the whole interval')
c
c     a zero propagation ratio still keeps one propagation step
c
      lmdaintv = 10
      lmdaparatio = 0.0d0
      lmdapbratio = 0.3d0
      call setlmdaphase
      call assert_int (lmdanpa,1,'setlmdaphase keeps one propagation')
      call assert_int (lmdanpa+lmdanpb+lmdanpc,lmdaintv,
     &                 'setlmdaphase spans a clamped interval')
c
c     a crowded interval gives back samples to the averaging phase,
c     taking them from the equilibration phase first
c
      lmdaintv = 10
      lmdaparatio = 0.4d0
      lmdapbratio = 0.5d0
      call setlmdaphase
      call assert_int (lmdanpc,2,'setlmdaphase restores the average')
      call assert_int (lmdanpa,4,'setlmdaphase spares the propagation')
      call assert_int (lmdanpb,4,'setlmdaphase trims the equilibration')
      call assert_int (lmdanpa+lmdanpb+lmdanpc,lmdaintv,
     &                 'setlmdaphase spans a crowded interval')
c
c     the shortest usable interval still holds all three phases
c
      lmdaintv = 3
      lmdaparatio = 0.4d0
      lmdapbratio = 0.4d0
      call setlmdaphase
      call assert_int (lmdanpa,1,'setlmdaphase minimum propagation')
      call assert_int (lmdanpc,2,'setlmdaphase minimum average')
      call assert_int (lmdanpa+lmdanpb+lmdanpc,lmdaintv,
     &                 'setlmdaphase spans the minimum interval')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eost_ostgate  --  lambda freezing tests  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eost_ostgate" drives elmdadyn over one deposit interval and
c     checks that the lambda particle moves only during the leading
c     propagation phase, that lambda is then held exactly fixed, and
c     that the deposited gaussian sits on that frozen lambda
c
c
      subroutine test_eost_ostgate
      use bath
      use dlmda
      use mutant
      use ost
      implicit none
      integer istep
      real*8 lam(6)
      real*8 frozen
c
c
c     drive an interval with a deterministic frictionless lambda
c     particle, so that any lambda motion comes from the gate alone
c
      kelvin = 300.0d0
      call resetost (5,5,4)
      use_ost = .true.
      lmdaintv = 6
      lmdanpa = 2
      lmdanpb = 2
      lmdanpc = 2
      use_lmdacv = .true.
      lmdacvstd = 1.0d0
      lmdacvrat = 0.0d0
      hbias = 1.0d0
      lmdadt = 0.1d0
      lmdamass = 1.0d0
      lmdafric = 0.0d0
      lmdatheta = 0.25d0 * 3.14159265358979323846d0
      lmdavtheta = 0.0d0
      lambda = 0.5d0
      fastkernel = .true.
      d2edl2 = 0.0d0
      ostbdgdl = 0.0d0
      ostbdgdfl = 0.0d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      do istep = 1, lmdaintv
         dedl = 1.0d0
         call elmdadyn
         lam(istep) = lambda
      end do
c
c     the particle moves while the interval is in its first phase
c
      call assert_logical (lam(1).ne.lam(2),.true.,
     &                     'elmdadyn ost propagates during phase a')
c
c     lambda is then bit identical for the rest of the interval
c
      frozen = lam(2)
      do istep = 3, lmdaintv
         call assert_real (lam(istep),frozen,0.0d0,
     &                     'elmdadyn ost holds lambda after phase a')
      end do
c
c     the averaged lambda is the frozen value, not a smear, so the
c     gaussian is deposited exactly on it
c
      call assert_int (nlmdahist,1,
     &                 'elmdadyn ost deposits a frozen interval')
      call assert_real (lmdaavg,frozen,0.0d0,
     &                  'elmdadyn ost averages the frozen lambda')
      call assert_real (lmdalhist(1),frozen,0.0d0,
     &                  'elmdadyn ost gaussian on the frozen lambda')
      call assert_real (lmdafhist(1),1.0d0,1.0d-12,
     &                  'elmdadyn ost centers the gaussian on flambda')
      use_lmdacv = .false.
      use_ost = .false.
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine resetost  --  reset OST test state  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "resetost" allocates OST arrays and initializes scalar
c     state to deterministic unit-test defaults
c
c
      subroutine resetost (nl,nf,nhist)
      use dlmda
      use mutant
      use ost
      implicit none
      integer nl,nf,nhist
      integer i,j
c
c
c     clear any previous allocation
c
      if (allocated(osthhist))  deallocate (osthhist)
      if (allocated(osthist))  deallocate (osthist)
      if (allocated(lmdaihist))  deallocate (lmdaihist)
      if (allocated(osthead))  deallocate (osthead)
      if (allocated(ostnext))  deallocate (ostnext)
      if (allocated(lmdallist))  deallocate (lmdallist)
      if (allocated(lmdaflist))  deallocate (lmdaflist)
      if (allocated(lmdalhist))  deallocate (lmdalhist)
      if (allocated(lmdafhist))  deallocate (lmdafhist)
      if (allocated(ostwlhist))  deallocate (ostwlhist)
      if (allocated(ostwfhist))  deallocate (ostwfhist)
      if (allocated(lmdafmean))  deallocate (lmdafmean)
      if (allocated(lmdafsum))  deallocate (lmdafsum)
      if (allocated(gfkernel))  deallocate (gfkernel)
      if (allocated(gkernel))  deallocate (gkernel)
      if (allocated(glfkernel))  deallocate (glfkernel)
      if (allocated(glkernel))  deallocate (glkernel)
      if (allocated(lmdafwt))  deallocate (lmdafwt)
      if (allocated(vkernelmax))  deallocate (vkernelmax)
      if (allocated(metalhist))  deallocate (metalhist)
      if (allocated(metahhist))  deallocate (metahhist)
      if (allocated(metawhist))  deallocate (metawhist)
c
c     set scalar state
c
      nlmda = nl
      nflmda = nf
      fli0 = (nflmda + 1) / 2
      wlmda = 1.0d0 / dble(nlmda-1)
      wflmda = 1.0d0
      wlmda2 = 0.5d0 * wlmda
      wflmda2 = 0.5d0 * wflmda
      wlhist = 0.005d0
      wfhist = 1.0d0
      maxwlhist = wlhist
      maxwfhist = wfhist
      nlmdahist = 0
      nlmdasave = 0
      sizelmdahist = nhist
      lmdathmap = 'SIN'
      lmdathalpha = 0.999999999d0
      lmdaintv = 10
      lmdanpa = 3
      lmdanpb = 3
      lmdanpc = 4
      lambda = 0.0d0
      lmdaavg = 0.0d0
      lmdastd = 0.0d0
      dedl = 0.0d0
      dedlavg = 0.0d0
      dedlstd = 0.0d0
      deffdl = 0.0d0
      use_lmdacv = .false.
      use_lmdafill = .false.
      lmdafillr = 0.1d0
      lmdacvstd = 0.0d0
      lmdacvrat = 0.0d0
      use_ost = .false.
      use_meta = .false.
      use_abf = .false.
      use_ostgtemp = .false.
      use_ostltemp = .false.
      ostgthresh = 0.0d0
      ostgtempgamma = 1.0d0
      ostlthresh = 0.0d0
      ostltempgamma = 1.0d0
      plmdamap = 'QNT'
      elmdamap = 'QNT'
      vlmdamap = 'QNT'
      plmdaexp = 1
      elmdaexp = 1
      vlmdaexp = 1
      plmdainvn = 1
      elmdainvn = 1
      vlmdainvn = 1
      plmdainveps = 0.0d0
      elmdainveps = 0.0d0
      vlmdainveps = 0.0d0
      lmdaparatio = 0.3d0
      lmdapbratio = 0.3d0
      hbias = 0.0d0
      lmdadeltag = 0.0d0
      oststdev = 1.0d0
      ostinterpol = .false.
      fastkernel = .false.
c
c     allocate arrays
c
      allocate (osthist(sizelmdahist))
      allocate (lmdaihist(sizelmdahist))
      allocate (osthead(nlmda,nflmda))
      allocate (ostnext(sizelmdahist))
      allocate (lmdallist(lmdaintv))
      allocate (lmdaflist(lmdaintv))
      allocate (lmdalhist(sizelmdahist))
      allocate (lmdafhist(sizelmdahist))
      allocate (osthhist(sizelmdahist))
      allocate (ostwlhist(sizelmdahist))
      allocate (ostwfhist(sizelmdahist))
      allocate (lmdafmean(nlmda))
      allocate (lmdafsum(nlmda))
      allocate (gfkernel(nlmda,nflmda))
      allocate (gkernel(nlmda,nflmda))
      allocate (glfkernel(nlmda,nflmda))
      allocate (glkernel(nlmda,nflmda))
      allocate (lmdafwt(nlmda))
      allocate (vkernelmax(nlmda))
c
c     initialize arrays
c
      do i = 1, sizelmdahist
         osthist(i) = 0
         lmdaihist(i) = 0
         ostnext(i) = 0
         lmdalhist(i) = 0.0d0
         lmdafhist(i) = 0.0d0
         osthhist(i) = 0.0d0
         ostwlhist(i) = 0.0d0
         ostwfhist(i) = 0.0d0
      end do
      do i = 1, lmdaintv
         lmdallist(i) = 0.0d0
         lmdaflist(i) = 0.0d0
      end do
      do i = 1, nlmda
         lmdafmean(i) = 0.0d0
         lmdafsum(i) = 0.0d0
         lmdafwt(i) = 0.0d0
         vkernelmax(i) = 0.0d0
         do j = 1, nflmda
            osthead(i,j) = 0
            gfkernel(i,j) = 0.0d0
            gkernel(i,j) = 0.0d0
            glfkernel(i,j) = 0.0d0
            glkernel(i,j) = 0.0d0
         end do
      end do
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine resetmeta  --  reset meta history  ##
c     ##                                                ##
c     ####################################################
c
c
c     "resetmeta" allocates metadynamics history arrays and the bin
c     lookup, and clears their contents for a test case
c
c
      subroutine resetmeta (nhist)
      use dlmda
      use ost
      implicit none
      integer nhist
      integer i
c
c
c     clear and allocate metadynamics history arrays
c
      if (allocated(metaihist))  deallocate (metaihist)
      if (allocated(metalhist))  deallocate (metalhist)
      if (allocated(metahhist))  deallocate (metahhist)
      if (allocated(metawhist))  deallocate (metawhist)
      if (allocated(metanext))  deallocate (metanext)
      if (allocated(metahead))  deallocate (metahead)
      if (allocated(vmetagrid))  deallocate (vmetagrid)
      if (allocated(dvmetagrid))  deallocate (dvmetagrid)
      sizemetahist = nhist
      nmetahist = 0
      allocate (metaihist(sizemetahist))
      allocate (metalhist(sizemetahist))
      allocate (metahhist(sizemetahist))
      allocate (metawhist(sizemetahist))
      allocate (metanext(sizemetahist))
      allocate (metahead(nlmda))
      allocate (vmetagrid(nlmda))
      allocate (dvmetagrid(nlmda))
      do i = 1, sizemetahist
         metaihist(i) = 0
         metalhist(i) = 0.0d0
         metahhist(i) = 0.0d0
         metawhist(i) = 0.0d0
         metanext(i) = 0
      end do
      do i = 1, nlmda
         metahead(i) = 0
         vmetagrid(i) = 0.0d0
         dvmetagrid(i) = 0.0d0
      end do
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine seedkernels  --  seed kernel tags  ##
c     ##                                                ##
c     ####################################################
c
c
c     "seedkernels" fills each flambda-dependent kernel
c     with values that encode its array and bin indices
c
c
      subroutine seedkernels
      use dlmda
      use ost
      implicit none
      integer i,j
      real*8 kerneltag
c
c
      do i = 1, nlmda
         do j = 1, nflmda
            gkernel(i,j) = kerneltag (1,i,j)
            gfkernel(i,j) = kerneltag (2,i,j)
            glkernel(i,j) = kerneltag (3,i,j)
            glfkernel(i,j) = kerneltag (4,i,j)
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  function kerneltag  --  encode a kernel bin  ##
c     ##                                               ##
c     ###################################################
c
c
c     "kerneltag" returns a unique numeric tag for a
c     kernel array and lambda/flambda bin pair
c
c
      function kerneltag (ikern,i,j)
      implicit none
      integer ikern,i,j
      real*8 kerneltag
c
c
      kerneltag = 1000.0d0*dble(ikern) + 10.0d0*dble(i) + dble(j)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine checkkernels  --  check resized kernels  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "checkkernels" compares resized kernels with reference
c     arrays whose seeded block has shifted by offset bins
c
c
      subroutine checkkernels (label,nold,offset)
      use dlmda
      use ost
      implicit none
      integer i,j
      integer nold,offset
      real*8 kerneltag
      real*8, allocatable :: gref(:,:)
      real*8, allocatable :: gfref(:,:)
      real*8, allocatable :: glref(:,:)
      real*8, allocatable :: glfref(:,:)
      character*(*) label
c
c
c     bins outside the copied block must be zero after the resize
c
      allocate (gref(nlmda,nflmda))
      allocate (gfref(nlmda,nflmda))
      allocate (glref(nlmda,nflmda))
      allocate (glfref(nlmda,nflmda))
      do i = 1, nlmda
         do j = 1, nflmda
            gref(i,j) = 0.0d0
            gfref(i,j) = 0.0d0
            glref(i,j) = 0.0d0
            glfref(i,j) = 0.0d0
         end do
      end do
c
c     the seeded block keeps its values but shifts by offset bins
c
      do i = 1, nlmda
         do j = 1, nold
            gref(i,j+offset) = kerneltag (1,i,j)
            gfref(i,j+offset) = kerneltag (2,i,j)
            glref(i,j+offset) = kerneltag (3,i,j)
            glfref(i,j+offset) = kerneltag (4,i,j)
         end do
      end do
      call assert_array2 (gkernel,gref,nlmda,nflmda,1.0d-12,
     &                    label//' gkernel')
      call assert_array2 (gfkernel,gfref,nlmda,nflmda,1.0d-12,
     &                    label//' gfkernel')
      call assert_array2 (glkernel,glref,nlmda,nflmda,1.0d-12,
     &                    label//' glkernel')
      call assert_array2 (glfkernel,glfref,nlmda,nflmda,1.0d-12,
     &                    label//' glfkernel')
      deallocate (gref)
      deallocate (gfref)
      deallocate (glref)
      deallocate (glfref)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine sethist  --  set one history entry  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "sethist" stores one gaussian history entry and its
c     packed lookup bin in the OST arrays
c
c
      subroutine sethist (ihist,lambda,flmda,height,sigl,sigf)
      use dlmda
      use ost
      implicit none
      integer ihist
      integer ilmda,iflmda,k
      integer lmdabin,flambdabin
      real*8 lambda,flmda
      real*8 height,sigl,sigf
c
c
c     save real gaussian center, parameters and packed lookup bin
c
      ilmda = lmdabin(lambda)
      iflmda = flambdabin(flmda)
      call ij_to_k (ilmda,iflmda,nlmda,k)
      osthist(ihist) = k
      lmdalhist(ihist) = lambda
      lmdafhist(ihist) = flmda
      osthhist(ihist) = height
      ostwlhist(ihist) = sigl
      ostwfhist(ihist) = sigf
      maxwlhist = max(maxwlhist,sigl)
      maxwfhist = max(maxwfhist,sigf)
      ostnext(ihist) = 0
      return
      end
