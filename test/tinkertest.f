c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  program tinkertest  --  Tinker unit test dispatcher  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "tinkertest" is the single Tinker unit test binary; it links
c     against the Tinker library and dispatches to the individual test
c     routines that live in the category files (amoeba.f, ...), in the
c     spirit of the tinker-gpu test executable
c
c     each test sets up its own structure with "loadfix", drives the
c     energy routines at levels 0, 1 and 3, and checks the results
c     against a stored reference with the assert primitives
c
c     usage:   tinkertest [filter] [-l level] [-v]
c
c     the binary is run from the test directory; the optional filter
c     selects tests whose name or tag contains the given string, e.g.
c     "tinkertest amoeba" runs only the AMOEBA-tagged cases
c
c     the level runs every case up to that cost level, from 1 (base)
c     to "maxlevel" (everything), as listed in the assert module; the default
c     is 1 with no filter and 5 with a filter, so that a filtered run
c     keeps all the cases of its tag
c
c
      program tinkertest
      use assert
      implicit none
      integer i,nargs
      integer level,ios
      character*240 root,filter,arg
      character*240 levstr
      logical detail,getlev
c
c
c     record the test root, optional name/tag filter, test level and
c     detail flag; "-v" (or "-d") turns on the per-check detail lines
c     and "-l N" (or "-lN") sets the test level
c
      call getcwd (root)
      filter = ' '
      levstr = ' '
      detail = .false.
      getlev = .false.
      nargs = command_argument_count ()
      do i = 1, nargs
         call get_command_argument (i,arg)
         if (getlev) then
            levstr = arg
            getlev = .false.
         else if (arg.eq.'-v' .or. arg.eq.'-d') then
            detail = .true.
         else if (arg .eq. '-l') then
            getlev = .true.
         else if (arg(1:2) .eq. '-l') then
            levstr = arg(3:)
         else if (len_trim(filter) .eq. 0) then
            filter = arg
         end if
      end do
c
c     set the default test level or check the requested one
c
      if (getlev)  levstr = '(missing)'
      if (len_trim(levstr) .eq. 0) then
         level = 1
         if (len_trim(filter) .ne. 0)  level = maxlevel
      else
         ios = 1
         level = 0
         if (verify(trim(levstr),'0123456789').eq.0 .and.
     &       len_trim(levstr).le.9) then
            read (levstr,'(i9)',iostat=ios)  level
         end if
         if (ios.ne.0 .or. level.lt.1 .or. level.gt.maxlevel) then
            write (*,10)  maxlevel,trim(levstr)
   10       format (/,' TINKERTEST  --  Test Level Must be 1 to',i2,
     &                 ', not "',a,'"')
            call exit (1)
         end if
      end if
      write (*,20)  level
   20 format (/,' Test Level',i3)
      call ta_init (root,filter,detail,level)
c
c     run unit tests
c
      call test_amoeba
      call test_angle
      call test_angtor
      call test_aplusgas
      call test_aplusliquid
      call test_barostat
      call test_barostat_nonprism
      call test_bond
      call test_bounds
      call test_box
      call test_chglj
      call test_chgtrn
      call test_disp
      call test_dodecahedron
      call test_eabf
      call test_eost
      call test_eostmap
      call test_emhippo
      call test_ephippo
      call test_expol
      call test_extfield
      call test_fep
      call test_geom
      call test_improp
      call test_imptor
      call test_kinetic
      call test_lmdamc
      call test_localframe
      call test_localframe2
      call test_localframe3
      call test_mpolestate
      call test_mutate
      call test_nacl
      call test_nblist
      call test_opbend
      call test_partialcharge
      call test_pitors
      call test_readxyz
      call test_polpair
      call test_rattle
      call test_rephippo
      call test_strbnd
      call test_strtor
      call test_thermint
      call test_torsion
      call test_tortor
      call test_truncated8
      call test_urey
      call test_vdw14
      call test_water4
c
c     run program tests
c
      call test_analyzelmda
      call test_bar
      call test_potential
      call test_pdbxyz
      call test_stochastic
      call test_testgrad
      call test_testlmda
      call test_testvir
      call test_tinkernist
      call test_verlet
      call test_xyzedit
      call test_xyzpdb
c
c     report the tally and set the process exit status
c
      call assert_summary ()
      end
