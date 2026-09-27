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
      call test_testlmda_case ('01_water_adt_l05','1e-4')
      call test_testlmda_case ('02_water_ast_l05','1e-4')
      call test_testlmda_case ('03_water_adt_l06exp','1e-4')
      call test_testlmda_case ('04_water_ast_l06exp','1e-4')
      call test_testlmda_case ('05_water_ast_nodl_l05','1e-4')
      call test_testlmda_case ('06_water_ast_vonly_l05','1e-4')
      call test_testlmda_case ('07_water_ast_ye_l10','1e-4')
      call test_testlmda_case ('08_water_ast_ye_l05','1e-4')
      call test_testlmda_case ('09_water_ast_ye_l00','1e-4')
      call test_testlmda_case ('10_water_ast_vcorr_l05','1e-4')
      call test_testlmda_case ('11_water_ast_vcorr_annih_l05','1e-4')
      call test_testlmda_case ('12_water_ast_vcorr_l06exp','1e-4')
      call test_testlmda_case ('14_water_rels_vdwm_vcorr_annih_l05',
     &                         '1e-4')
      call test_testlmda_case ('16_water_rels_vdwm_lights_l05','1e-4')
      call test_testlmda_case ('17_water_rels_vdwm_vcorr_l050','1e-5')
      call test_testlmda_case ('18_water_rels_lig1_l085','1e-4')
      call test_testlmda_case ('19_water_rels_lig2_l015','1e-4')
      call test_testlmda_case ('20_water_rels_lig1_ne_l085','1e-4')
      call test_testlmda_case ('21_water_rels_lig1_nlist_exf_l085',
     &                         '1e-4')
      call test_testlmda_case ('22_water_rels_lig1_st_l085','1e-4')
      call test_testlmda_case ('23_water_rels_lig2_st_l015','1e-4')
      call test_testlmda_case ('24_water_rels_lig1_st_ne_l085','1e-4')
      call test_testlmda_case ('25_water_rels_lig1_st_polonly_l085',
     &                         '1e-4')
      call test_testlmda_case ('26_water_rels_lig1_dt_polonly_l078',
     &                         '1e-4')
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
c     both the analytical and the numerical lambda derivatives with the
c     given finite difference step, and checks the full program output
c     against the reference; the tolerance is loose enough to absorb the
c     finite difference noise in the numerical derivatives, and a steep
c     lambda map takes a smaller step to keep its truncation error low
c
c
      subroutine test_testlmda_case (base,step)
      implicit none
      integer ist
      logical skiptest
      character*(*) base,step
      character*240 rpath
      character*512 args
c
c
      if (skiptest('test_testlmda_'//base,'testlmda'))  return
      call pushdir ('file/testlmda')
      call execute_command_line ('rm -f out.txt')
      args = '-k '//base//'.key water2.xyz Y Y '//step
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
