c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################
c     ##                                              ##
c     ##  subroutine test_mutate  --  mutation tests  ##
c     ##                                              ##
c     ##################################################
c
c
c     "test_mutate" checks AMOEBA mutation energy, gradient, virial
c     and named energy component regressions
c
c
      subroutine test_mutate
      implicit none
c
c
      call test_mutate_mv
      call test_mutate_mp
      call test_mutate_ast
      call test_mutate_adt
      call test_mutate_rdt
      call test_mutate_qnt
      call test_mutate_exp
      call test_mutate_inv
      call test_mutate_exf
      call test_mutate_emplar
      call test_mutate_rels
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_mv  --  mutation lambda-scan cases  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_mv" runs the nine water mutation fixtures; cases
c     001-006 isolate electrostatics at three ele-lambda values with
c     Ewald on and off, while cases 007-009 isolate van der Waals at
c     three vdw-lambda values
c
c
      subroutine test_mutate_mv
      implicit none
c
c
      call test_mutate_fixed ('001_water_ye_m10.key',
     &   '001_water_ye_m10.txt','001_water_ye_m10',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('002_water_ne_m10.key',
     &   '002_water_ne_m10.txt','002_water_ne_m10',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('003_water_ye_m05.key',
     &   '003_water_ye_m05.txt','003_water_ye_m05',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('004_water_ne_m05.key',
     &   '004_water_ne_m05.txt','004_water_ne_m05',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('005_water_ye_m00.key',
     &   '005_water_ye_m00.txt','005_water_ye_m00',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('006_water_ne_m00.key',
     &   '006_water_ne_m00.txt','006_water_ne_m00',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('007_water_v10.key',
     &   '007_water_v10.txt','007_water_v10',
     &   .false., .false., .true.,  .true.)
      call test_mutate_fixed ('008_water_v05.key',
     &   '008_water_v05.txt','008_water_v05',
     &   .false., .false., .true.,  .true.)
      call test_mutate_fixed ('009_water_v00.key',
     &   '009_water_v00.txt','009_water_v00',
     &   .false., .false., .true.,  .true.)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_mp  --  electrostatic lambda cases  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_mp" runs the twenty water mutation fixtures that
c     scan the electrostatic lambda values; cases 010-015 keep only the
c     multipole term at three ele-lambda values with Ewald on and off,
c     cases 016-021 keep only the polarization term at three pol-lambda
c     values, and cases 022-029 leave both terms active while varying
c     ele-lambda and pol-lambda together
c
c
      subroutine test_mutate_mp
      implicit none
c
c
      call test_mutate_fixed ('010_water_ye_m10.key',
     &   '010_water_ye_m10.txt','010_water_ye_m10',
     &   .true.,  .false., .false., .true.)
      call test_mutate_fixed ('011_water_ne_m10.key',
     &   '011_water_ne_m10.txt','011_water_ne_m10',
     &   .true.,  .false., .false., .false.)
      call test_mutate_fixed ('012_water_ye_m05.key',
     &   '012_water_ye_m05.txt','012_water_ye_m05',
     &   .true.,  .false., .false., .true.)
      call test_mutate_fixed ('013_water_ne_m05.key',
     &   '013_water_ne_m05.txt','013_water_ne_m05',
     &   .true.,  .false., .false., .false.)
      call test_mutate_fixed ('014_water_ye_m00.key',
     &   '014_water_ye_m00.txt','014_water_ye_m00',
     &   .true.,  .false., .false., .true.)
      call test_mutate_fixed ('015_water_ne_m00.key',
     &   '015_water_ne_m00.txt','015_water_ne_m00',
     &   .true.,  .false., .false., .false.)
      call test_mutate_fixed ('016_water_ye_p10.key',
     &   '016_water_ye_p10.txt','016_water_ye_p10',
     &   .false., .true.,  .false., .true.)
      call test_mutate_fixed ('017_water_ne_p10.key',
     &   '017_water_ne_p10.txt','017_water_ne_p10',
     &   .false., .true.,  .false., .false.)
      call test_mutate_fixed ('018_water_ye_p05.key',
     &   '018_water_ye_p05.txt','018_water_ye_p05',
     &   .false., .true.,  .false., .true.)
      call test_mutate_fixed ('019_water_ne_p05.key',
     &   '019_water_ne_p05.txt','019_water_ne_p05',
     &   .false., .true.,  .false., .false.)
      call test_mutate_fixed ('020_water_ye_p00.key',
     &   '020_water_ye_p00.txt','020_water_ye_p00',
     &   .false., .true.,  .false., .true.)
      call test_mutate_fixed ('021_water_ne_p00.key',
     &   '021_water_ne_p00.txt','021_water_ne_p00',
     &   .false., .true.,  .false., .false.)
      call test_mutate_fixed ('022_water_ye_m10p05.key',
     &   '022_water_ye_m10p05.txt','022_water_ye_m10p05',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('023_water_ne_m10p05.key',
     &   '023_water_ne_m10p05.txt','023_water_ne_m10p05',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('024_water_ye_m05p10.key',
     &   '024_water_ye_m05p10.txt','024_water_ye_m05p10',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('025_water_ne_m05p10.key',
     &   '025_water_ne_m05p10.txt','025_water_ne_m05p10',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('026_water_ye_m05p00.key',
     &   '026_water_ye_m05p00.txt','026_water_ye_m05p00',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('027_water_ne_m05p00.key',
     &   '027_water_ne_m05p00.txt','027_water_ne_m05p00',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('028_water_ye_m00p05.key',
     &   '028_water_ye_m00p05.txt','028_water_ye_m00p05',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('029_water_ne_m00p05.key',
     &   '029_water_ne_m00p05.txt','029_water_ne_m00p05',
     &   .true.,  .true.,  .false., .false.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_mutate_ast  --  single topology lambda  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_mutate_ast" runs the eleven absolute single topology water
c     fixtures 030-040, each carrying the "lambda-deriv" keyword, and
c     drives them through "test_mutate_calc" with the level 4 lambda
c     derivative checks enabled; cases 030-035 keep only the multipole
c     term at three ele-lambda values with Ewald on and off, cases
c     036-038 keep only the van der Waals term at three vdw-lambda
c     values, and cases 039-040 leave the multipole and polarization
c     terms active with Ewald on and off; the no-Ewald cases cannot use
c     a pairwise neighbor list
c
c
      subroutine test_mutate_ast
      implicit none
c
c
      call test_mutate_calc ('water2','030_water_ast_ye_m10.key',
     &   '030_water_ast_ye_m10.txt','030_water_ast_ye_m10',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','031_water_ast_ne_m10.key',
     &   '031_water_ast_ne_m10.txt','031_water_ast_ne_m10',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','032_water_ast_ye_m05.key',
     &   '032_water_ast_ye_m05.txt','032_water_ast_ye_m05',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','033_water_ast_ne_m05.key',
     &   '033_water_ast_ne_m05.txt','033_water_ast_ne_m05',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','034_water_ast_ye_m00.key',
     &   '034_water_ast_ye_m00.txt','034_water_ast_ye_m00',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','035_water_ast_ne_m00.key',
     &   '035_water_ast_ne_m00.txt','035_water_ast_ne_m00',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','036_water_ast_v10.key',
     &   '036_water_ast_v10.txt','036_water_ast_v10',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','037_water_ast_v05.key',
     &   '037_water_ast_v05.txt','037_water_ast_v05',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','038_water_ast_v00.key',
     &   '038_water_ast_v00.txt','038_water_ast_v00',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','039_water_ast_ye_m05p10.key',
     &   '039_water_ast_ye_m05p10.txt','039_water_ast_ye_m05p10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','040_water_ast_ne_m05p10.key',
     &   '040_water_ast_ne_m05p10.txt','040_water_ast_ne_m05p10',
     &   .true.,  .true.,  .false., .false., .true.)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_mutate_adt  --  abs dual topo lambda  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_mutate_adt" runs the seventeen absolute dual topology water
c     fixtures 041-057, each carrying the "lambda-deriv" keyword and a
c     dual topology keyword, and drives them through "test_mutate_calc"
c     with the level 4 lambda derivative checks enabled; cases 041-046
c     keep only the multipole term at three ele-lambda values with Ewald
c     on and off, cases 047-052 keep only the polarization term at three
c     pol-lambda values, cases 053-055 keep only the van der Waals term
c     at three vdw-lambda values, and cases 056-057 leave the multipole
c     and polarization terms active with Ewald on and off; the no-Ewald
c     cases cannot use a pairwise neighbor list
c
c
      subroutine test_mutate_adt
      implicit none
c
c
      call test_mutate_calc ('water2','041_water_adt_ye_m10.key',
     &   '041_water_adt_ye_m10.txt','041_water_adt_ye_m10',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','042_water_adt_ne_m10.key',
     &   '042_water_adt_ne_m10.txt','042_water_adt_ne_m10',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','043_water_adt_ye_m05.key',
     &   '043_water_adt_ye_m05.txt','043_water_adt_ye_m05',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','044_water_adt_ne_m05.key',
     &   '044_water_adt_ne_m05.txt','044_water_adt_ne_m05',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','045_water_adt_ye_m00.key',
     &   '045_water_adt_ye_m00.txt','045_water_adt_ye_m00',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','046_water_adt_ne_m00.key',
     &   '046_water_adt_ne_m00.txt','046_water_adt_ne_m00',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','047_water_adt_ye_p10.key',
     &   '047_water_adt_ye_p10.txt','047_water_adt_ye_p10',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','048_water_adt_ne_p10.key',
     &   '048_water_adt_ne_p10.txt','048_water_adt_ne_p10',
     &   .false., .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','049_water_adt_ye_p05.key',
     &   '049_water_adt_ye_p05.txt','049_water_adt_ye_p05',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','050_water_adt_ne_p05.key',
     &   '050_water_adt_ne_p05.txt','050_water_adt_ne_p05',
     &   .false., .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','051_water_adt_ye_p00.key',
     &   '051_water_adt_ye_p00.txt','051_water_adt_ye_p00',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','052_water_adt_ne_p00.key',
     &   '052_water_adt_ne_p00.txt','052_water_adt_ne_p00',
     &   .false., .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','053_water_adt_v10.key',
     &   '053_water_adt_v10.txt','053_water_adt_v10',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','054_water_adt_v05.key',
     &   '054_water_adt_v05.txt','054_water_adt_v05',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','055_water_adt_v00.key',
     &   '055_water_adt_v00.txt','055_water_adt_v00',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','056_water_adt_ye_m05p10.key',
     &   '056_water_adt_ye_m05p10.txt','056_water_adt_ye_m05p10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','057_water_adt_ne_m05p10.key',
     &   '057_water_adt_ne_m05p10.txt','057_water_adt_ne_m05p10',
     &   .true.,  .true.,  .false., .false., .true.)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_mutate_rdt  --  rel dual topo lambda  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_mutate_rdt" runs the seventeen relative dual topology water
c     fixtures 058-074, each carrying the "lambda-deriv" keyword, a dual
c     topology keyword and a pair of "ligand1" and "ligand2" groups, and
c     drives them through "test_mutate_calc" with the level 4 lambda
c     derivative checks enabled; cases 058-063 keep only the multipole
c     term at three ele-lambda values with Ewald on and off, cases
c     064-069 keep only the polarization term at three pol-lambda
c     values, cases 070-072 keep only the van der Waals term at three
c     vdw-lambda values, and cases 073-074 leave the multipole and
c     polarization terms active with Ewald on and off; the no-Ewald
c     cases cannot use a pairwise neighbor list
c
c
      subroutine test_mutate_rdt
      implicit none
c
c
      call test_mutate_calc ('water2','058_water_rdt_ye_m10.key',
     &   '058_water_rdt_ye_m10.txt','058_water_rdt_ye_m10',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','059_water_rdt_ne_m10.key',
     &   '059_water_rdt_ne_m10.txt','059_water_rdt_ne_m10',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','060_water_rdt_ye_m05.key',
     &   '060_water_rdt_ye_m05.txt','060_water_rdt_ye_m05',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','061_water_rdt_ne_m05.key',
     &   '061_water_rdt_ne_m05.txt','061_water_rdt_ne_m05',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','062_water_rdt_ye_m00.key',
     &   '062_water_rdt_ye_m00.txt','062_water_rdt_ye_m00',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','063_water_rdt_ne_m00.key',
     &   '063_water_rdt_ne_m00.txt','063_water_rdt_ne_m00',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','064_water_rdt_ye_p10.key',
     &   '064_water_rdt_ye_p10.txt','064_water_rdt_ye_p10',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','065_water_rdt_ne_p10.key',
     &   '065_water_rdt_ne_p10.txt','065_water_rdt_ne_p10',
     &   .false., .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','066_water_rdt_ye_p05.key',
     &   '066_water_rdt_ye_p05.txt','066_water_rdt_ye_p05',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','067_water_rdt_ne_p05.key',
     &   '067_water_rdt_ne_p05.txt','067_water_rdt_ne_p05',
     &   .false., .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','068_water_rdt_ye_p00.key',
     &   '068_water_rdt_ye_p00.txt','068_water_rdt_ye_p00',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','069_water_rdt_ne_p00.key',
     &   '069_water_rdt_ne_p00.txt','069_water_rdt_ne_p00',
     &   .false., .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','070_water_rdt_v10.key',
     &   '070_water_rdt_v10.txt','070_water_rdt_v10',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','071_water_rdt_v05.key',
     &   '071_water_rdt_v05.txt','071_water_rdt_v05',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','072_water_rdt_v00.key',
     &   '072_water_rdt_v00.txt','072_water_rdt_v00',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','073_water_rdt_ye_m05p10.key',
     &   '073_water_rdt_ye_m05p10.txt','073_water_rdt_ye_m05p10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','074_water_rdt_ne_m05p10.key',
     &   '074_water_rdt_ne_m05p10.txt','074_water_rdt_ne_m05p10',
     &   .true.,  .true.,  .false., .false., .true.)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_mutate_qnt  --  quintic lambda map  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_mutate_qnt" runs the nine water fixtures 075-083 that map
c     the main lambda to the electrostatics, polarization and van der
c     Waals sub-lambdas with the quintic "qnt" scheme; each carries the
c     "ost" keyword and drives "test_mutate_calc" with the level 4
c     lambda derivative checks enabled, verifying the multipole,
c     polarization and van der Waals components together; cases 075-077
c     use single topology, 078-080 absolute dual topology and 081-083
c     relative dual topology, each at ost-lambda values 1.0, 0.5 and
c     0.0; all fixtures use Ewald and support a pairwise neighbor list
c
c
      subroutine test_mutate_qnt
      implicit none
c
c
      call test_mutate_calc ('water2','075_water_qnt_ast_l10.key',
     &   '075_water_qnt_ast_l10.txt','075_water_qnt_ast_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','076_water_qnt_ast_l05.key',
     &   '076_water_qnt_ast_l05.txt','076_water_qnt_ast_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','077_water_qnt_ast_l00.key',
     &   '077_water_qnt_ast_l00.txt','077_water_qnt_ast_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','078_water_qnt_adt_l10.key',
     &   '078_water_qnt_adt_l10.txt','078_water_qnt_adt_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','079_water_qnt_adt_l05.key',
     &   '079_water_qnt_adt_l05.txt','079_water_qnt_adt_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','080_water_qnt_adt_l00.key',
     &   '080_water_qnt_adt_l00.txt','080_water_qnt_adt_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','081_water_qnt_rdt_l10.key',
     &   '081_water_qnt_rdt_l10.txt','081_water_qnt_rdt_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','082_water_qnt_rdt_l05.key',
     &   '082_water_qnt_rdt_l05.txt','082_water_qnt_rdt_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','083_water_qnt_rdt_l00.key',
     &   '083_water_qnt_rdt_l00.txt','083_water_qnt_rdt_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      return
      end
c
c
c     ######################################################
c     ##                                                  ##
c     ##  subroutine test_mutate_exp  --  exp lambda map  ##
c     ##                                                  ##
c     ######################################################
c
c
c     "test_mutate_exp" runs the nine water fixtures 084-092 that map
c     the main lambda to the electrostatics, polarization and van der
c     Waals sub-lambdas with the exponential "exp" scheme; each carries
c     the "ost" keyword and drives "test_mutate_calc" with the level 4
c     lambda derivative checks enabled, verifying the multipole,
c     polarization and van der Waals components together; cases 084-086
c     use single topology, 087-089 absolute dual topology and 090-092
c     relative dual topology, each at ost-lambda values 1.0, 0.5 and
c     0.0; all fixtures use Ewald and support a pairwise neighbor list
c
c
      subroutine test_mutate_exp
      implicit none
c
c
      call test_mutate_calc ('water2','084_water_exp_ast_l10.key',
     &   '084_water_exp_ast_l10.txt','084_water_exp_ast_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','085_water_exp_ast_l05.key',
     &   '085_water_exp_ast_l05.txt','085_water_exp_ast_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','086_water_exp_ast_l00.key',
     &   '086_water_exp_ast_l00.txt','086_water_exp_ast_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','087_water_exp_adt_l10.key',
     &   '087_water_exp_adt_l10.txt','087_water_exp_adt_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','088_water_exp_adt_l05.key',
     &   '088_water_exp_adt_l05.txt','088_water_exp_adt_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','089_water_exp_adt_l00.key',
     &   '089_water_exp_adt_l00.txt','089_water_exp_adt_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','090_water_exp_rdt_l10.key',
     &   '090_water_exp_rdt_l10.txt','090_water_exp_rdt_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','091_water_exp_rdt_l05.key',
     &   '091_water_exp_rdt_l05.txt','091_water_exp_rdt_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','092_water_exp_rdt_l00.key',
     &   '092_water_exp_rdt_l00.txt','092_water_exp_rdt_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_mutate_inv  --  inverse lambda map  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_mutate_inv" runs the nine water fixtures 093-101 that map
c     the main lambda to the electrostatics, polarization and van der
c     Waals sub-lambdas with the inverse "inv" scheme; each carries the
c     "ost" keyword and drives "test_mutate_calc" with the level 4
c     lambda derivative checks enabled, verifying the multipole,
c     polarization and van der Waals components together; cases 093-095
c     use single topology, 096-098 absolute dual topology and 099-101
c     relative dual topology, each at ost-lambda values 1.0, 0.5 and
c     0.0; all fixtures use Ewald and support a pairwise neighbor list
c
c
      subroutine test_mutate_inv
      implicit none
c
c
      call test_mutate_calc ('water2','093_water_inv_ast_l10.key',
     &   '093_water_inv_ast_l10.txt','093_water_inv_ast_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','094_water_inv_ast_l05.key',
     &   '094_water_inv_ast_l05.txt','094_water_inv_ast_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','095_water_inv_ast_l00.key',
     &   '095_water_inv_ast_l00.txt','095_water_inv_ast_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','096_water_inv_adt_l10.key',
     &   '096_water_inv_adt_l10.txt','096_water_inv_adt_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','097_water_inv_adt_l05.key',
     &   '097_water_inv_adt_l05.txt','097_water_inv_adt_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','098_water_inv_adt_l00.key',
     &   '098_water_inv_adt_l00.txt','098_water_inv_adt_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','099_water_inv_rdt_l10.key',
     &   '099_water_inv_rdt_l10.txt','099_water_inv_rdt_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','100_water_inv_rdt_l05.key',
     &   '100_water_inv_rdt_l05.txt','100_water_inv_rdt_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','101_water_inv_rdt_l00.key',
     &   '101_water_inv_rdt_l00.txt','101_water_inv_rdt_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_mutate_exf  --  applied external field  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_mutate_exf" runs the seventeen water fixtures 102-118 that
c     apply an external electric field to a mutated system, each
c     carrying the "lambda-deriv" keyword so the level 4 lambda
c     derivative checks are enabled; cases 102-110 keep only the
c     multipole term, with 102-104 single topology, 105-107 absolute
c     dual topology and 108-110 relative dual topology, each at
c     ele-lambda values 1.0, 0.5 and 0.0; cases 111-116 keep only the
c     polarization term, with 111-113 absolute and 114-116 relative
c     dual topology at pol-lambda values 1.0, 0.5 and 0.0, since the
c     induced dipoles respond to the applied field through the direct
c     field; cases 117-118 leave the multipole and polarization terms
c     active together, under absolute and relative dual topology, and
c     use unequal interpolation exponents of three and four so the two
c     terms cannot mask an error in each other; the remaining dual
c     topology cases use an exponent of two, since with the default
c     exponent of one the dual topology weighting of the external
c     field term cannot be distinguished from a linear scaling; all
c     fixtures use Ewald and support a pairwise neighbor list
c
c
      subroutine test_mutate_exf
      implicit none
c
c
      call test_mutate_calc ('water2','102_water_exf_ast_m10.key',
     &   '102_water_exf_ast_m10.txt','102_water_exf_ast_m10',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','103_water_exf_ast_m05.key',
     &   '103_water_exf_ast_m05.txt','103_water_exf_ast_m05',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','104_water_exf_ast_m00.key',
     &   '104_water_exf_ast_m00.txt','104_water_exf_ast_m00',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','105_water_exf_adt_m10.key',
     &   '105_water_exf_adt_m10.txt','105_water_exf_adt_m10',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','106_water_exf_adt_m05.key',
     &   '106_water_exf_adt_m05.txt','106_water_exf_adt_m05',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','107_water_exf_adt_m00.key',
     &   '107_water_exf_adt_m00.txt','107_water_exf_adt_m00',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','108_water_exf_rdt_m10.key',
     &   '108_water_exf_rdt_m10.txt','108_water_exf_rdt_m10',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','109_water_exf_rdt_m05.key',
     &   '109_water_exf_rdt_m05.txt','109_water_exf_rdt_m05',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','110_water_exf_rdt_m00.key',
     &   '110_water_exf_rdt_m00.txt','110_water_exf_rdt_m00',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','111_water_exf_adt_p10.key',
     &   '111_water_exf_adt_p10.txt','111_water_exf_adt_p10',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','112_water_exf_adt_p05.key',
     &   '112_water_exf_adt_p05.txt','112_water_exf_adt_p05',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','113_water_exf_adt_p00.key',
     &   '113_water_exf_adt_p00.txt','113_water_exf_adt_p00',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','114_water_exf_rdt_p10.key',
     &   '114_water_exf_rdt_p10.txt','114_water_exf_rdt_p10',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','115_water_exf_rdt_p05.key',
     &   '115_water_exf_rdt_p05.txt','115_water_exf_rdt_p05',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','116_water_exf_rdt_p00.key',
     &   '116_water_exf_rdt_p00.txt','116_water_exf_rdt_p00',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','117_water_exf_adt_m05p10.key',
     &   '117_water_exf_adt_m05p10.txt','117_water_exf_adt_m05p10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','118_water_exf_rdt_m05p05.key',
     &   '118_water_exf_rdt_m05p05.txt','118_water_exf_rdt_m05p05',
     &   .true.,  .true.,  .false., .true.,  .true.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_mutate_emplar  --  mpole plus polar dt  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_mutate_emplar" runs the twelve water fixtures 119-130 that
c     keep the multipole and polarization terms active together at
c     matched ele-lambda and pol-lambda values, each carrying the
c     "lambda-deriv" keyword and a dual topology keyword, and drives
c     them through "test_mutate_calc" with the level 4 lambda derivative
c     checks enabled; cases 119-124 use absolute dual topology and cases
c     125-130 relative dual topology, each at ele-/pol-lambda values
c     1.0, 0.5, and 0.0 with Ewald on and off; both terms use dual
c     topology interpolation exponents of three so neither term can
c     mask an error in the other; the no-Ewald cases cannot use a
c     pairwise neighbor list
c
c
      subroutine test_mutate_emplar
      implicit none
c
c
      call test_mutate_calc ('water2','119_water_adt_ye_m10p10v05.key',
     &   '119_water_adt_ye_m10p10v05.txt','119_water_adt_ye_m10p10v05',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','120_water_adt_ne_m10p10v05.key',
     &   '120_water_adt_ne_m10p10v05.txt','120_water_adt_ne_m10p10v05',
     &   .true.,  .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','121_water_adt_ye_m05p05v00.key',
     &   '121_water_adt_ye_m05p05v00.txt','121_water_adt_ye_m05p05v00',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','122_water_adt_ne_m05p05v00.key',
     &   '122_water_adt_ne_m05p05v00.txt','122_water_adt_ne_m05p05v00',
     &   .true.,  .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','123_water_adt_ye_m00p00v10.key',
     &   '123_water_adt_ye_m00p00v10.txt','123_water_adt_ye_m00p00v10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','124_water_adt_ne_m00p00v10.key',
     &   '124_water_adt_ne_m00p00v10.txt','124_water_adt_ne_m00p00v10',
     &   .true.,  .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','125_water_rdt_ye_m10p10v05.key',
     &   '125_water_rdt_ye_m10p10v05.txt','125_water_rdt_ye_m10p10v05',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','126_water_rdt_ne_m10p10v05.key',
     &   '126_water_rdt_ne_m10p10v05.txt','126_water_rdt_ne_m10p10v05',
     &   .true.,  .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','127_water_rdt_ye_m05p05v00.key',
     &   '127_water_rdt_ye_m05p05v00.txt','127_water_rdt_ye_m05p05v00',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','128_water_rdt_ne_m05p05v00.key',
     &   '128_water_rdt_ne_m05p05v00.txt','128_water_rdt_ne_m05p05v00',
     &   .true.,  .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','129_water_rdt_ye_m00p00v10.key',
     &   '129_water_rdt_ye_m00p00v10.txt','129_water_rdt_ye_m00p00v10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','130_water_rdt_ne_m00p00v10.key',
     &   '130_water_rdt_ne_m00p00v10.txt','130_water_rdt_ne_m00p00v10',
     &   .true.,  .true.,  .false., .false., .true.)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_mutate_rels  --  staged rel dual topo  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_mutate_rels" runs the seven water fixtures 135-141 that
c     drive the staged relative free energy schedule with a main lambda,
c     so the two ligands are discharged and recharged one at a time while
c     van der Waals morphs between them in the middle window; the cases
c     sit one per regime of the schedule, at main lambda values 1.0 and
c     0.0 where a leg is flat and only the coupled endpoint is built,
c     0.85 and 0.15 inside the two mixing legs, 0.7 and 0.3 at the leg
c     boundaries where the van der Waals morph window opens and closes,
c     and 0.5 in the middle of that window where both ligands are
c     decoupled and every electrostatic lambda derivative vanishes; all
c     seven carry the "lambda-deriv" keyword and run the level 4 lambda
c     derivative checks, and the reference values agree with the tinker9
c     CUDA implementation of the same schedule
c
c
      subroutine test_mutate_rels
      implicit none
c
c
      call test_mutate_calc ('water2','135_water_rels_ye_l100.key',
     &   '135_water_rels_ye_l100.txt','135_water_rels_ye_l100',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','136_water_rels_ye_l085.key',
     &   '136_water_rels_ye_l085.txt','136_water_rels_ye_l085',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','137_water_rels_ye_l070.key',
     &   '137_water_rels_ye_l070.txt','137_water_rels_ye_l070',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','138_water_rels_ye_l050.key',
     &   '138_water_rels_ye_l050.txt','138_water_rels_ye_l050',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','139water_rels_ye_l030.key',
     &   '139water_rels_ye_l030.txt','139water_rels_ye_l030',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','140_water_rels_ye_l015.key',
     &   '140_water_rels_ye_l015.txt','140_water_rels_ye_l015',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','141_water_rels_ye_l000.key',
     &   '141_water_rels_ye_l000.txt','141_water_rels_ye_l000',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_mutate_fixed  --  one mutation case  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_mutate_fixed" runs a single mutation fixture with and
c     without the neighbor-list keyword, checking the level 0/1/3
c     energy, gradient, virial and named component regressions; it is
c     a thin wrapper around "test_mutate_calc" with the level 4 lambda
c     derivative checks disabled, and backs the "test_mutate_mv" and
c     "test_mutate_mp" case lists
c
c
      subroutine test_mutate_fixed
     &   (key,ref,cname,checkm,checkp,checkv,canlist)
      implicit none
      logical checkm,checkp,checkv,canlist
      character*(*) key,ref,cname
c
c
      call test_mutate_calc ('water',key,ref,cname,checkm,checkp,
     &                       checkv,canlist,.false.)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_mutate_calc  --  one mutation case  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_mutate_calc" runs a single mutation fixture with and
c     without the neighbor-list keyword; for each neighbor-list
c     variant the force field is built once, then the checks are
c     repeated twice before teardown, giving four passes per fixture;
c     the "checkm", "checkp" and "checkv" flags select which level 3
c     energy components are verified (Atomic Multipoles, Polarization
c     and Van der Waals); the "canlist" flag records whether the
c     neighbor-list variant is compatible with this fixture (false for
c     the no-Ewald cases); when "dolmda" is set, the level 4 checks
c     also verify the analytical lambda derivatives, second lambda
c     derivatives, per-atom lambda gradient and dV/dL tensor against
c     the reference read by "load_lmdaref", reusing the values already
c     produced by the level 1 gradient call
c
c
      subroutine test_mutate_calc
     &   (base,key,ref,cname,checkm,checkp,checkv,canlist,dolmda)
      use action
      use atoms
      use dlmda
      use energi
      use virial
      implicit none
      integer nat,natl,irun,ilist,nlist
      real*8 energy,e,ref_e,ref_ei
      real*8 eps_e,eps_g,eps_v,eps_l,refv(3,3)
      real*8 ref_dedl(4),ref_d2edl2(4),ref_dvdl(3,3)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: refg(:,:)
      real*8, allocatable :: ref_lg(:,:)
      logical skiptest,checkm,checkp,checkv,canlist,uselist,dolmda
      character*(*) base,key,ref,cname
      character*240 rpath,pre
      character*8 rtag
c
c
      if (skiptest('test_'//trim(cname),'mutate'))  return
c
c     run each fixture without the neighbor-list keyword, and with it
c     when the fixture supports a pairwise neighbor list
c
      nlist = 2
      if (.not. canlist)  nlist = 1
      do ilist = 1, nlist
         uselist = (ilist .eq. 2)
c
c     set up the force field once for this neighbor-list variant
c
         call pushdir ('file/mutate')
         if (uselist) then
            call loadfix_keyadd (base,key,'neighbor-list')
            pre = 'test_'//trim(cname)//' list'
         else
            call loadfix (base,key)
            pre = 'test_'//trim(cname)//' nolist'
         end if
         allocate (derivs(3,n))
         allocate (refg(3,n))
         call refpath ('mutate',ref,rpath)
         call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
         eps_e = 1.0d-4
         eps_g = 1.0d-4
         eps_v = 1.0d-3
c
c     read the reference lambda derivatives for the level 4 checks
c
         if (dolmda) then
            allocate (ref_lg(3,n))
            call load_lmdaref (rpath,n,ref_dedl,ref_d2edl2,ref_lg,
     &                         ref_dvdl,natl)
            eps_l = 1.0d-4
         end if
c
c     repeat the checks twice against the built system
c
         do irun = 1, 2
            if (irun .eq. 1) then
               rtag = ' run1'
            else
               rtag = ' run2'
            end if
c
c     level 0  --  total potential energy
c
            e = energy ()
            call assert_real (esum,ref_e,eps_e,
     &                        trim(pre)//' energy (v0)'//trim(rtag))
c
c     level 1  --  total energy, Cartesian gradient and virial
c
            call gradient (e,derivs)
            call assert_real (esum,ref_e,eps_e,
     &                        trim(pre)//' grad-e (v1)'//trim(rtag))
            call assert_grad (derivs,refg,n,eps_g,
     &                        trim(pre)//' grad (v1)'//trim(rtag))
            call assert_grad (vir,refv,3,eps_v,
     &                        trim(pre)//' virial (v1)'//trim(rtag))
c
c     level 4  --  lambda derivatives from the level 1 gradient call
c
            if (dolmda) then
               call assert_real (dedl,ref_dedl(1),eps_l,
     &                        trim(pre)//' dE/dL (v4)'//trim(rtag))
               call assert_real (devdl,ref_dedl(2),eps_l,
     &                        trim(pre)//' dEV/dL (v4)'//trim(rtag))
               call assert_real (demdl,ref_dedl(3),eps_l,
     &                        trim(pre)//' dEM/dL (v4)'//trim(rtag))
               call assert_real (depdl,ref_dedl(4),eps_l,
     &                        trim(pre)//' dEP/dL (v4)'//trim(rtag))
               call assert_real (d2edl2,ref_d2edl2(1),eps_l,
     &                        trim(pre)//' d2E/dL2 (v4)'//trim(rtag))
               call assert_real (d2evdl2,ref_d2edl2(2),eps_l,
     &                        trim(pre)//' d2EV/dL2 (v4)'//trim(rtag))
               call assert_real (d2emdl2,ref_d2edl2(3),eps_l,
     &                        trim(pre)//' d2EM/dL2 (v4)'//trim(rtag))
               call assert_real (d2epdl2,ref_d2edl2(4),eps_l,
     &                        trim(pre)//' d2EP/dL2 (v4)'//trim(rtag))
               call assert_grad (dfsumdl,ref_lg,n,eps_g,
     &                        trim(pre)//' lgrad (v4)'//trim(rtag))
               call assert_grad (dvirdl,ref_dvdl,3,eps_v,
     &                        trim(pre)//' dV/dL (v4)'//trim(rtag))
            end if
c
c     level 3  --  total and named AMOEBA energy components
c
            call analysis (e)
            call assert_real (esum,ref_e,eps_e,
     &                        trim(pre)//' analysis (v3)'//trim(rtag))
            if (checkm) then
               call check_engcnt (rpath,'Atomic Multipoles',em,nem,
     &                       eps_e,trim(pre)//' mpole (v3)'//trim(rtag))
            end if
            if (checkp) then
               call check_engcnt (rpath,'Polarization',ep,nep,
     &                       eps_e,trim(pre)//' polar (v3)'//trim(rtag))
            end if
            if (checkv) then
               call check_engcnt (rpath,'Van der Waals',ev,nev,
     &                       eps_e,trim(pre)//' vdw (v3)'//trim(rtag))
            end if
         end do
c
c     clean up this neighbor-list variant
c
         deallocate (derivs)
         deallocate (refg)
         if (dolmda)  deallocate (ref_lg)
         call popdir
         call final
      end do
      return
      end
