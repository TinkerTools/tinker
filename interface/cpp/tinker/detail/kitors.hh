#pragma once

#include "macro.hh"

namespace tinker { namespace kitors {
extern int& maxnti;
extern double*& ti1;
extern double*& ti2;
extern double*& ti3;
extern char (*&kti)[16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kitors, maxnti);
extern "C" double* TINKER_MOD(kitors, ti1);
extern "C" double* TINKER_MOD(kitors, ti2);
extern "C" double* TINKER_MOD(kitors, ti3);
extern "C" char (*TINKER_MOD(kitors, kti))[16];

int& maxnti = TINKER_MOD(kitors, maxnti);
double*& ti1 = TINKER_MOD(kitors, ti1);
double*& ti2 = TINKER_MOD(kitors, ti2);
double*& ti3 = TINKER_MOD(kitors, ti3);
char (*&kti)[16] = TINKER_MOD(kitors, kti);
#endif
} }
