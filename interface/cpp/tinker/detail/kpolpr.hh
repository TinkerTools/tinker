#pragma once

#include "macro.hh"

namespace tinker { namespace kpolpr {
extern int& maxnpp;
extern double*& thlpr;
extern double*& thdpr;
extern char (*&kppr)[8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kpolpr, maxnpp);
extern "C" double* TINKER_MOD(kpolpr, thlpr);
extern "C" double* TINKER_MOD(kpolpr, thdpr);
extern "C" char (*TINKER_MOD(kpolpr, kppr))[8];

int& maxnpp = TINKER_MOD(kpolpr, maxnpp);
double*& thlpr = TINKER_MOD(kpolpr, thlpr);
double*& thdpr = TINKER_MOD(kpolpr, thdpr);
char (*&kppr)[8] = TINKER_MOD(kpolpr, kppr);
#endif
} }
