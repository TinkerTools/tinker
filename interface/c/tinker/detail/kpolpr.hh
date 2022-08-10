#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(kpolpr, maxnpp);
extern double* TINKER_MOD(kpolpr, thlpr);
extern double* TINKER_MOD(kpolpr, thdpr);
extern char (*TINKER_MOD(kpolpr, kppr))[8];
#ifdef __cplusplus
}
#endif
