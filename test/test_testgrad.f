c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_testgrad  --  testgrad program tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     this file collects Tinker system tests that run the "testgrad"
c     program as a subprocess and compare its output against a stored
c     reference with a floating-point tolerance
c
c
      subroutine test_testgrad
      implicit none
c
c
      call test_testgrad_case ('01_water_ye_m10v10','water2')
      call test_testgrad_case ('02_water_ye_m05v05','water2')
      call test_testgrad_case ('03_water_ye_m00v00','water2')
      call test_testgrad_case ('04_water_ast_ye_l10','water2')
      call test_testgrad_case ('05_water_ast_ye_l05','water2')
      call test_testgrad_case ('06_water_ast_ye_l00','water2')
      call test_testgrad_case ('09_g3','g3')
      call test_testgrad_case ('10_g3_nlist','g3')
      call test_testgrad_case ('11_aplus_thole0','tetramer')
      call test_testgrad_case ('12_aplus_polpair','tetramer')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_testgrad_case  --  one testgrad case  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_testgrad_case" runs "testgrad" on the given coordinates with
c     the key file and reference that share the given base name, computing
c     both the analytical and the numerical gradient, and checks the full
c     program output against the reference; the tolerance is loose enough
c     to absorb the finite difference noise in the numerical gradient
c
c
      subroutine test_testgrad_case (base,xyz)
      implicit none
      integer ist
      logical skiptest
      character*(*) base,xyz
      character*240 rpath
      character*512 args
c
c
      if (skiptest('test_testgrad_'//base,'testgrad'))  return
      call pushdir ('file/testgrad')
      call execute_command_line ('rm -f out.txt')
c
c     the trailing answer declines the per-component breakdown that
c     "testgrad" only offers for systems above one hundred atoms
c
      args = '-k '//base//'.key '//xyz//' Y Y 0.00001 N'
      call run_prog ('testgrad',trim(args),'out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call refpath ('testgrad',base//'.txt',rpath)
      call assert_files (rpath,'out.txt',1.0d-3,'test_testgrad '//base)
      call execute_command_line ('rm -f out.txt')
      call popdir
      return
      end
