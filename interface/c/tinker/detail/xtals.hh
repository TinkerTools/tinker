#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxlsq 1000
#define TINKER_MOD__maxrsd 1000
extern int TINKER_MOD(xtals, nxtal);
extern int TINKER_MOD(xtals, nvary);
extern int TINKER_MOD(xtals, ivary)[TINKER_MOD__maxlsq];
extern int TINKER_MOD(xtals, iresid)[TINKER_MOD__maxrsd];
extern int TINKER_MOD(xtals, vary)[TINKER_MOD__maxlsq][2];
extern double TINKER_MOD(xtals, e0_lattice);
extern char TINKER_MOD(xtals, varxtl)[TINKER_MOD__maxlsq][16];
extern char TINKER_MOD(xtals, rsdxtl)[TINKER_MOD__maxrsd][16];
#ifdef __cplusplus
}
#endif
