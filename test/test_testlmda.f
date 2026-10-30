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
c     cases 91 to 100 cover HIPPO Pauli repulsion, which follows the
c     van der Waals lambda and is counted in dEV/dL, for a pentane and
c     acetamide mixture: the default soft core by double loop and by
c     neighbor list, annihilation, an early taper, the shift shapes
c     (1-lambda)^2 and a non-integer one from REPEL-SOFTCORE, both
c     lambda endpoints, and first derivatives alone with an exponent
c     below two
c
c     cases 101 to 110 cover HIPPO damped dispersion, whose soft core
c     also follows the van der Waals lambda and is counted in dEV/dL,
c     for the same mixture: the default soft core, then a switch radius
c     that reaches the closest mutated pairs by neighbor list, under
c     annihilation and with an early taper, Ewald summation by double
c     loop and by neighbor list with annihilation, non-integer soft
c     core exponents, both lambda endpoints under Ewald, and first
c     derivatives alone with a lambda exponent of one
c
c     cases 111 and 112 cover Ewald dispersion in replicate cells, for
c     two DMSO molecules in a small box, and under atom group weights
c
c     cases 113 to 116 add the long range dispersion correction to a
c     soft core with a 2.5 Angstrom switch radius, alone and under
c     annihilation, with a lambda exponent of one, and with an exponent
c     of two at zero lambda
c
c     case 117 drives repulsion through a cubic van der Waals lambda
c     map, so both map derivatives enter the chain rule, which every
c     linear map case leaves at one and zero
c
c     cases 118 and 119 give the buffered 14-7 soft core a softening
c     of alpha*lambda^m*(1-lambda)^p from VDW-SOFTCORE, with m of 1.5
c     and p of 3, then a non-integer m of 0.5 and p of one
c
c     case 120 drives Ewald dispersion through the same cubic map, so
c     the soft core powers of the real space, reciprocal space and self
c     terms all carry both map derivatives
c
c     case 121 drives case 113 through the cubic map, so the long range
c     dispersion correction carries both map derivatives as well
c
c     cases 122 to 125 take the other van der Waals map shapes: 122 and
c     123 drive the repulsion of case 117 through a quintic taper over
c     0.2 to 0.9 and an inverse power map with n of 3, 124 drives the
c     dispersion correction of case 113 through the same taper, and 125
c     the Ewald dispersion of case 120 through an inverse power map with
c     n of 2, so the curvature of each map enters with its own sign
c
c     cases 126 to 133 put HIPPO repulsion on the staged relative legs
c     of a water box, where the second ligand couples as one minus the
c     van der Waals lambda and the two ligands never interact: the van
c     der Waals leg by double loop and by multipole list, two water
c     ligands under annihilation, the quintic map, the three legs with
c     lambda on each one, and replicate cells in a small box; two
c     ligands at the same site are left to "test_mutate", since the
c     coordinate check reports them in thread order
c
c     cases 134 to 145 put HIPPO damped dispersion on the same relative
c     legs, where reciprocal space carries the environment and the two
c     ligands as separate densities: without Ewald by double loop and
c     by neighbor list, two water ligands under annihilation, and the
c     long range correction; then under Ewald by double loop and by
c     neighbor list, under annihilation, the quintic map, the three
c     legs with lambda on each one, where the ligand legs merge or drop
c     a ligand density, and replicate cells in a small box
c
c     case 147 decouples a water whose hal tail is handed to an
c     undamped Ewald dispersion sum, which follows the hal soft core
c     lambda power in place of the long range correction
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
      call test_testlmda_xyz ('91_hipmix_rep_l04','hipmix','1e-4')
      call test_testlmda_xyz ('92_hipmix_rep_nlist_l04','hipmix','1e-4')
      call test_testlmda_xyz ('93_hipmix_rep_annih_l04','hipmix','1e-4')
      call test_testlmda_xyz ('94_hipmix_rep_taper_l04','hipmix','1e-4')
      call test_testlmda_xyz ('95_hipmix_rep_sc_m0p2_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('96_hipmix_rep_sc_real_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('97_hipmix_rep_l00','hipmix','1e-4')
      call test_testlmda_xyz ('98_hipmix_rep_l10','hipmix','1e-4')
      call test_testlmda_xyz ('99_hipmix_rep_sc_m0p2_l10','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('100_hipmix_rep_d1_n12_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('101_hipmix_dsp_l04','hipmix','1e-4')
      call test_testlmda_xyz ('102_hipmix_dsp_nlist_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('103_hipmix_dsp_annih_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('104_hipmix_dsp_taper_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('105_hipmix_dsp_ewald_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('106_hipmix_dsp_ewald_nlist_annih_l04',
     &                        'hipmix','1e-4')
      call test_testlmda_xyz ('107_hipmix_dsp_sc_real_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('108_hipmix_dsp_ewald_n2_l00','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('109_hipmix_dsp_ewald_l10','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('110_hipmix_dsp_d1_n1_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('111_dmso2_dsp_ewald_replica_l04',
     &                        '../testgrad/dmso2','1e-4')
      call test_testlmda_xyz ('112_hipmix_dsp_ewald_group_l04',
     &                        'hipmix','1e-4')
      call test_testlmda_xyz ('113_hipmix_dsp_corr_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('114_hipmix_dsp_corr_annih_l04',
     &                        'hipmix','1e-4')
      call test_testlmda_xyz ('115_hipmix_dsp_corr_n1_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('116_hipmix_dsp_corr_n2_l00','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('117_hipmix_rep_exp3_l08','hipmix',
     &                        '1e-4')
      call test_testlmda_case ('118_water_vsoft_m15_p3_l05','1e-4')
      call test_testlmda_case ('119_water_vsoft_m05_p1_l03','1e-4')
      call test_testlmda_xyz ('120_hipmix_dsp_ewald_exp3_l08','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('121_hipmix_dsp_corr_exp3_l06','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('122_hipmix_rep_qnt_l045','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('123_hipmix_rep_inv_l03','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('124_hipmix_dsp_corr_qnt_l06','hipmix',
     &                        '1e-4')
      call test_testlmda_xyz ('125_hipmix_dsp_ewald_inv_l04','hipmix',
     &                        '1e-4')
      call test_testlmda_case ('126_water_rep_rels_vdwm_l05','1e-4')
      call test_testlmda_case ('127_water_rep_rels_vdwm_nlist_l05',
     &                         '1e-4')
      call test_testlmda_case ('128_water_rep_rels_vdwm_annih_l05',
     &                         '1e-4')
      call test_testlmda_case ('129_water_rep_rels_vdwm_qnt_l035',
     &                         '1e-4')
      call test_testlmda_case ('130_water_rep_rels3_l050','1e-4')
      call test_testlmda_case ('131_water_rep_rels3_l085','1e-4')
      call test_testlmda_case ('132_water_rep_rels3_l015','1e-4')
      call test_testlmda_case ('133_water_rep_rels_vdwm_replica_l05',
     &                         '1e-4')
      call test_testlmda_case ('134_water_dsp_rels_vdwm_l05','1e-4')
      call test_testlmda_case ('135_water_dsp_rels_vdwm_nlist_l05',
     &                         '1e-4')
      call test_testlmda_case ('136_water_dsp_rels_vdwm_annih_l05',
     &                         '1e-4')
      call test_testlmda_case ('137_water_dsp_rels_vdwm_corr_l03',
     &                         '1e-4')
      call test_testlmda_case ('138_water_dsp_ewald_rels_vdwm_l04',
     &                         '1e-4')
      call test_testlmda_case
     &   ('139_water_dsp_ewald_rels_vdwm_nlist_l04','1e-4')
      call test_testlmda_case
     &   ('140_water_dsp_ewald_rels_vdwm_annih_l04','1e-4')
      call test_testlmda_case ('141_water_dsp_ewald_rels_vdwm_qnt_l035',
     &                         '1e-4')
      call test_testlmda_case ('142_water_dsp_ewald_rels3_l050','1e-4')
      call test_testlmda_case ('143_water_dsp_ewald_rels3_l085','1e-4')
      call test_testlmda_case ('144_water_dsp_ewald_rels3_l015','1e-4')
      call test_testlmda_case
     &   ('145_water_dsp_ewald_rels_vdwm_replica_l04','1e-4')
      call test_testlmda_case
     &   ('146_water_dsp_ewald_annih_all_l00','1e-4')
      call test_testlmda_case ('147_water_dhal_l05','1e-4')
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
