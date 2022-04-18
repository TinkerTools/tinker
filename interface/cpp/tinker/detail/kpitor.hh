#pragma once

#include "macro.hh"

namespace tinker { namespace kpitor {
extern int& maxnpt;
extern double*& ptcon;
extern char (*&kpt)[8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kpitor, maxnpt);
extern "C" double* TINKER_MOD(kpitor, ptcon);
extern "C" char (*TINKER_MOD(kpitor, kpt))[8];

int& maxnpt = TINKER_MOD(kpitor, maxnpt);
double*& ptcon = TINKER_MOD(kpitor, ptcon);
char (*&kpt)[8] = TINKER_MOD(kpitor, kpt);
#endif
} }
