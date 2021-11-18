#pragma once

#include "macro.hh"
#include "sizes.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(merck, nligne);
extern int TINKER_MOD(merck, bt_1)[2][500];
extern int TINKER_MOD(merck, eqclass)[5][500];
extern int TINKER_MOD(merck, crd)[100];
extern int TINKER_MOD(merck, val)[100];
extern int TINKER_MOD(merck, pilp)[100];
extern int TINKER_MOD(merck, mltb)[100];
extern int TINKER_MOD(merck, arom)[100];
extern int TINKER_MOD(merck, lin)[100];
extern int TINKER_MOD(merck, sbmb)[100];
extern int TINKER_MOD(merck, mmffarom)[6][TINKER_MOD__maxtyp];
extern int TINKER_MOD(merck, mmffaromc)[6][TINKER_MOD__maxtyp];
extern int TINKER_MOD(merck, mmffaroma)[6][TINKER_MOD__maxtyp];
extern double TINKER_MOD(merck, rad0)[100];
extern double TINKER_MOD(merck, paulel)[100];
extern double TINKER_MOD(merck, r0ref)[100][100];
extern double TINKER_MOD(merck, kbref)[100][100];
extern double TINKER_MOD(merck, mmff_kb)[100][100];
extern double TINKER_MOD(merck, mmff_kb1)[100][100];
extern double TINKER_MOD(merck, mmff_b0)[100][100];
extern double TINKER_MOD(merck, mmff_b1)[100][100];
extern double* TINKER_MOD(merck, mmff_ka);
extern double* TINKER_MOD(merck, mmff_ka1);
extern double* TINKER_MOD(merck, mmff_ka2);
extern double* TINKER_MOD(merck, mmff_ka3);
extern double* TINKER_MOD(merck, mmff_ka4);
extern double* TINKER_MOD(merck, mmff_ka5);
extern double* TINKER_MOD(merck, mmff_ka6);
extern double* TINKER_MOD(merck, mmff_ka7);
extern double* TINKER_MOD(merck, mmff_ka8);
extern double* TINKER_MOD(merck, mmff_ang0);
extern double* TINKER_MOD(merck, mmff_ang1);
extern double* TINKER_MOD(merck, mmff_ang2);
extern double* TINKER_MOD(merck, mmff_ang3);
extern double* TINKER_MOD(merck, mmff_ang4);
extern double* TINKER_MOD(merck, mmff_ang5);
extern double* TINKER_MOD(merck, mmff_ang6);
extern double* TINKER_MOD(merck, mmff_ang7);
extern double* TINKER_MOD(merck, mmff_ang8);
extern double* TINKER_MOD(merck, stbn_abc);
extern double* TINKER_MOD(merck, stbn_cba);
extern double* TINKER_MOD(merck, stbn_abc1);
extern double* TINKER_MOD(merck, stbn_cba1);
extern double* TINKER_MOD(merck, stbn_abc2);
extern double* TINKER_MOD(merck, stbn_cba2);
extern double* TINKER_MOD(merck, stbn_abc3);
extern double* TINKER_MOD(merck, stbn_cba3);
extern double* TINKER_MOD(merck, stbn_abc4);
extern double* TINKER_MOD(merck, stbn_cba4);
extern double* TINKER_MOD(merck, stbn_abc5);
extern double* TINKER_MOD(merck, stbn_cba5);
extern double* TINKER_MOD(merck, stbn_abc6);
extern double* TINKER_MOD(merck, stbn_cba6);
extern double* TINKER_MOD(merck, stbn_abc7);
extern double* TINKER_MOD(merck, stbn_cba7);
extern double* TINKER_MOD(merck, stbn_abc8);
extern double* TINKER_MOD(merck, stbn_cba8);
extern double* TINKER_MOD(merck, stbn_abc9);
extern double* TINKER_MOD(merck, stbn_cba9);
extern double* TINKER_MOD(merck, stbn_abc10);
extern double* TINKER_MOD(merck, stbn_cba10);
extern double* TINKER_MOD(merck, stbn_abc11);
extern double* TINKER_MOD(merck, stbn_cba11);
extern double TINKER_MOD(merck, defstbn_abc)[5][5][5];
extern double TINKER_MOD(merck, defstbn_cba)[5][5][5];
extern double TINKER_MOD(merck, t1_1)[2001][2];
extern double TINKER_MOD(merck, t2_1)[2001][2];
extern double TINKER_MOD(merck, t3_1)[2001][2];
extern double TINKER_MOD(merck, t1_2)[2001][2];
extern double TINKER_MOD(merck, t2_2)[2001][2];
extern double TINKER_MOD(merck, t3_2)[2001][2];
extern char TINKER_MOD(merck, kt_1)[2001][16];
extern char TINKER_MOD(merck, kt_2)[2001][16];
extern double TINKER_MOD(merck, g)[TINKER_MOD__maxclass];
extern double TINKER_MOD(merck, alph)[TINKER_MOD__maxclass];
extern double TINKER_MOD(merck, nn)[TINKER_MOD__maxclass];
extern char TINKER_MOD(merck, da)[TINKER_MOD__maxclass][1];
extern double TINKER_MOD(merck, bci)[100][100];
extern double TINKER_MOD(merck, bci_1)[100][100];
extern double TINKER_MOD(merck, pbci)[TINKER_MOD__maxclass];
extern double TINKER_MOD(merck, fcadj)[TINKER_MOD__maxclass];
#ifdef __cplusplus
}
#endif
