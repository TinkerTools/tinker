c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program tinkertest  --  Tinker unit test dispatcher         ##
c     ##                                                              ##
c     ##################################################################
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
c     usage:   tinkertest [filter]
c
c     the binary is run from the test directory; the optional filter
c     selects tests whose name or tag contains the given string, e.g.
c     "tinkertest amoeba" runs only the AMOEBA-tagged cases
c
c
      program tinkertest
      use assert
      implicit none
      integer i,nargs
      character*240 root,filter,arg
      logical detail
c
c
c     record the test root, optional name/tag filter and detail flag;
c     "-v" (or "-d") turns on the per-check detail lines
c
      call getcwd (root)
      filter = ' '
      detail = .false.
      nargs = command_argument_count ()
      do i = 1, nargs
         call get_command_argument (i,arg)
         if (arg.eq.'-v' .or. arg.eq.'-d') then
            detail = .true.
         else if (len_trim(filter) .eq. 0) then
            filter = arg
         end if
      end do
      call ta_init (root,filter,detail)
c
c     run unit tests
c
      call test_amoeba
      call test_angle
      call test_angtor
      call test_bond
      call test_geom
      call test_improp
      call test_imptor
      call test_opbend
      call test_pitors
      call test_strbnd
      call test_strtor
      call test_torsion
      call test_tortor
      call test_urey
c
c     run program tests
c
      call test_potential
      call test_pdbxyz
      call test_xyzpdb
c
c     report the tally and set the process exit status
c
      call assert_summary ()
      end
