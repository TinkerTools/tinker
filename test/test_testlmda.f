c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_testlmda  --  testlmda program tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     this file collects Tinker system tests that run the "testlmda"
c     program as a subprocess and compare its output against a stored
c     reference with a floating-point tolerance
c
c
      subroutine test_testlmda
      implicit none
c
c
      call test_testlmda_case ('01_water_adt_l05')
      call test_testlmda_case ('02_water_ast_l05')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_testlmda_case  --  one lambda derivative  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_testlmda_case" runs "testlmda" on the water mutation fixture
c     whose key file and reference share the given base name, computing
c     both the analytical and the numerical lambda derivatives, and
c     checks the full program output against the reference; the
c     tolerance is loose enough to absorb the finite difference noise in
c     the numerical derivatives
c
c
      subroutine test_testlmda_case (base)
      implicit none
      integer ist
      logical skiptest
      character*(*) base
      character*240 rpath
      character*512 args
c
c
      if (skiptest('test_testlmda_'//base,'testlmda'))  return
      call pushdir ('file/testlmda')
      call execute_command_line ('rm -f out.txt')
      args = '-k '//base//'.key water2.xyz Y Y 1e-4'
      call run_prog ('testlmda',trim(args),'out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call refpath ('testlmda',base//'.txt',rpath)
      call assert_files (rpath,'out.txt',1.0d-3,'test_testlmda '//base)
      call execute_command_line ('rm -f out.txt')
      call popdir
      return
      end
