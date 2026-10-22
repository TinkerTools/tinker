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
c     the ionwat cases decouple the chloride of a +2 environment under
c     Ewald, so the uniform background correction carries a lambda
c     scaled net charge into the dV/dL tensor from "empole4"; case 26
c     adds dual topology polarization alongside it
c
c     cases 28 to 40 cover multipoles with charge penetration and
c     charge flux: HIPPO water, an AMOEBA+ alcohol tetramer and gas
c     dimer, including a ligand atom whose charge flux bond reaches
c     into the environment, and a HIPPO organic mixture without flux
c
c     cases 41 to 54 add polarization to the same systems, taking the
c     single topology scalar derivative in direct, Ewald, neighbor list
c     and zero lambda forms, and the dual topology derivatives under
c     orthogonal space tempering
c
c     case 55 weights HIPPO water by atom group with charge flux, so
c     the flux potential and its lambda derivative carry the group
c     weight; cases 56 and 57 put the electrostatic and polarization
c     lambdas at different values, so the flux scale switches between
c     the two states; cases 58 and 59 use the multipole list without
c     Ewald; cases 60 to 63 are relative legs with charge penetration
c     and charge flux, for whole AMOEBA+ alcohols, single atoms whose
c     flux bonds reach into the environment, and HIPPO water
c
c     cases 64 and 65 drive HIPPO water polarization alone, with no
c     electrostatic lambda map, in single and dual topology
c
c     cases 73 to 80 cover charge transfer, which follows the
c     polarization lambda and is counted in dEP/dL: HIPPO water charge
c     transfer alone in single and dual topology, with a squared dual
c     topology weight, the multipole list and a zero lambda endpoint,
c     then with multipoles and polarization at split lambdas, under
c     Ewald, and for an AMOEBA+ alcohol tetramer; charge transfer
c     scales each pair by the product of its two group scales in
c     single topology and by the polarization weight in dual topology
c
c     cases 81 to 83 add replicate cells, group weights and a lambda
c     range that leaves the ligand pairs off
c
c     cases 84 to 90 put charge transfer on the staged relative legs
c     of an AMOEBA+ alcohol tetramer whose two ligands are in contact,
c     alone on the ligand 1 and ligand 2 legs in single and dual
c     topology and on the van der Waals leg, where it carries no
c     lambda derivative, then with multipoles, polarization and
c     charge flux under Ewald for the tetramer and for HIPPO water
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
      call test_testlmda_case ('07_water_ast_vcorr_l05','1e-4')
      call test_testlmda_case ('08_water_ast_vcorr_annih_l05','1e-4')
      call test_testlmda_case ('09_water_ast_vcorr_l06exp','1e-4')
      call test_testlmda_case ('10_water_rels_vdwm_vcorr_annih_l05',
     &                         '1e-4')
      call test_testlmda_case ('11_water_rels_vdwm_vcorr_l050','1e-5')
      call test_testlmda_case ('12_water_rels_lig1_l085','1e-4')
      call test_testlmda_case ('13_water_rels_lig2_l015','1e-4')
      call test_testlmda_case ('14_water_rels_lig1_ne_l085','1e-4')
      call test_testlmda_case ('15_water_rels_lig1_nlist_exf_l085',
     &                         '1e-4')
      call test_testlmda_case ('16_water_rels_lig1_st_l085','1e-4')
      call test_testlmda_case ('17_water_rels_lig2_st_l015','1e-4')
      call test_testlmda_case ('18_water_rels_lig1_st_ne_l085','1e-4')
      call test_testlmda_case ('19_water_rels_lig1_st_polonly_l085',
     &                         '1e-4')
      call test_testlmda_case ('20_water_rels_lig1_dt_polonly_l078',
     &                         '1e-4')
      call test_testlmda_case ('21_water_rels_lig1_st_l088','1e-4')
      call test_testlmda_case ('22_water_rels_lig1_dt_l088','1e-4')
      call test_testlmda_case ('23_water_rels_3stg_l085','1e-4')
      call test_testlmda_xyz ('24_ionwat_ewald_l05','ionwat','1e-4')
      call test_testlmda_xyz ('25_ionwat_ewald_nlist_l05','ionwat',
     &                        '1e-4')
      call test_testlmda_xyz ('26_ionwat_pol_ewald_l05','ionwat','1e-4')
      call test_testlmda_case ('27_water_rels_2stg_l075','1e-4')
      call test_testlmda_xyz ('28_hipwat_flux_l04','hipwat','1e-4')
      call test_testlmda_xyz ('29_hipwat_flux_ewald_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('30_hipwat_flux_ewald_nlist_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('31_hipwat_flux_ewald_lin_l00','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('32_aptetra_flux_l04','aptetra','1e-4')
      call test_testlmda_xyz ('33_aptetra_flux_ewald_l04','aptetra',
     &                        '1e-4')
      call test_testlmda_xyz ('34_aptetra_flux_ewald_nlist_l04',
     &                        'aptetra','1e-4')
      call test_testlmda_xyz ('35_aptetra_flux_ewald_lin_l00','aptetra',
     &                        '1e-4')
      call test_testlmda_xyz ('36_apgas_flux_l04','apgas','1e-4')
      call test_testlmda_xyz ('37_apgas_cross_l04','apgas','1e-4')
      call test_testlmda_xyz ('38_apgas_cross_lin_l00','apgas','1e-4')
      call test_testlmda_xyz ('39_apgas_cross_ewald_l04','apgas','1e-4')
      call test_testlmda_xyz ('40_hipmix_ewald_l04','hipmix','1e-4')
      call test_testlmda_xyz ('41_hipwat_pol_st_l04','hipwat','1e-4')
      call test_testlmda_xyz ('42_hipwat_pol_ewald_st_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('43_hipwat_pol_ewald_nlist_st_l04',
     &                        'hipwat','1e-4')
      call test_testlmda_xyz ('44_hipwat_pol_ewald_lin_st_l00','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('45_hipwat_pol_ewald_dt_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('46_aptetra_pol_st_l04','aptetra','1e-4')
      call test_testlmda_xyz ('47_aptetra_pol_ewald_st_l04','aptetra',
     &                        '1e-4')
      call test_testlmda_xyz ('48_aptetra_pol_ewald_nlist_st_l04',
     &                        'aptetra','1e-4')
      call test_testlmda_xyz ('49_aptetra_pol_ewald_lin_st_l00',
     &                        'aptetra','1e-4')
      call test_testlmda_xyz ('50_aptetra_pol_ewald_dt_l04','aptetra',
     &                        '1e-4')
      call test_testlmda_xyz ('51_apgas_pol_st_l04','apgas','1e-4')
      call test_testlmda_xyz ('52_apgas_pol_ewald_dt_l04','apgas',
     &                        '1e-4')
      call test_testlmda_xyz ('53_hipmix_pol_ewald_st_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('54_hipmix_pol_ewald_dt_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('55_hipwat_group_ost_l04','hipwat','1e-4')
      call test_testlmda_xyz ('56_hipwat_pol_split_dt_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('57_hipwat_pol_split_st_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('58_hipwat_flux_mlist_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('59_aptetra_pol_mlist_dt_l04','aptetra',
     &                        '1e-4')
      call test_testlmda_xyz ('60_aptetra_rel_lig1_l085','aptetra',
     &                        '1e-4')
      call test_testlmda_xyz ('61_aptetra_rel_lig2_l015','aptetra',
     &                        '1e-4')
      call test_testlmda_xyz ('62_aptetra_rel_cross_l085','aptetra',
     &                        '1e-4')
      call test_testlmda_xyz ('63_hipwat_rel_lig1_ewald_l085','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('64_hipwat_pol_only_st_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('65_hipwat_pol_only_dt_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('66_expnaclus_na_dt_l04','expnaclus',
     &                        '1e-4')
      call test_testlmda_xyz ('67_expnaclus_wat_dt_l04','expnaclus',
     &                        '1e-4')
      call test_testlmda_xyz ('68_expbox_ewald_dt_l04','expbox','1e-4')
      call test_testlmda_xyz ('69_expnaclus_wat_st_l04','expnaclus',
     &                        '1e-4')
      call test_testlmda_xyz ('70_expbox_ewald_st_l04','expbox','1e-4')
      call test_testlmda_xyz ('71_expnaclus_na_lin_st_l00','expnaclus',
     &                        '1e-4')
      call test_testlmda_xyz ('72_expbox_ewald_nlist_st_l04','expbox',
     &                        '1e-4')
      call test_testlmda_xyz ('73_hipwat_ct_only_st_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('74_hipwat_ct_only_dt_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('75_hipwat_ct_only_dtexp_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('76_hipwat_ct_mlist_dt_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('77_hipwat_ct_lin_dt_l00','hipwat','1e-4')
      call test_testlmda_xyz ('78_hipwat_pol_ct_split_st_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('79_hipwat_pol_ct_ewald_dt_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('80_aptetra_ct_st_l04','aptetra','1e-4')
      call test_testlmda_xyz ('81_hipwat_ct_replica_ost_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('82_hipwat_ct_group_ost_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('83_hipwat_ct_off_st_l04','hipwat',
     &                        '1e-4')
      call test_testlmda_xyz ('84_aptetra_ct_rel_lig1_st_l085',
     &                        'aptetra','1e-4')
      call test_testlmda_xyz ('85_aptetra_ct_rel_lig2_st_l015',
     &                        'aptetra','1e-4')
      call test_testlmda_xyz ('86_aptetra_ct_rel_lig1_dt_l080',
     &                        'aptetra','1e-4')
      call test_testlmda_xyz ('87_aptetra_ct_rel_lig2_dt_l020',
     &                        'aptetra','1e-4')
      call test_testlmda_xyz ('88_aptetra_ct_rel_vdwm_l050','aptetra',
     &                        '1e-4')
      call test_testlmda_xyz ('89_aptetra_pol_ct_rel_lig1_ewald_l080',
     &                        'aptetra','1e-4')
      call test_testlmda_xyz ('90_hipwat_pol_ct_rel_lig1_ewald_l080',
     &                        'hipwat','1e-4')
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
c     via "test_testlmda_xyz"
c
c
      subroutine test_testlmda_case (base,step)
      implicit none
      character*(*) base,step
c
c
      call test_testlmda_xyz (base,'water2',step)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_testlmda_xyz  --  one lambda derivative  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_testlmda_xyz" runs "testlmda" on the given coordinates with
c     the key file and reference that share the given base name,
c     computing both the analytical and the numerical lambda derivatives
c     with the given finite difference step, and checks the full program
c     output against the reference; the tolerance is loose enough to
c     absorb the finite difference noise in the numerical derivatives,
c     and a steep lambda map takes a smaller step to keep its truncation
c     error low
c
c
      subroutine test_testlmda_xyz (base,xyz,step)
      implicit none
      integer ist
      logical skiptest
      character*(*) base,xyz,step
      character*240 rpath
      character*512 args
c
c
      if (skiptest('test_testlmda_'//base,'testlmda'))  return
      call pushdir ('file/testlmda')
      call execute_command_line ('rm -f out.txt')
      args = '-k '//base//'.key '//xyz//'.xyz Y Y '//step
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
