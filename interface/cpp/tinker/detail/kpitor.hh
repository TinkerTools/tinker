#pragma once

#include "macro.hh"

namespace tinker { namespace kpitor {
const int maxnpt = 500;
extern double (&ptcon)[maxnpt];
extern char (&kpt)[maxnpt][8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(kpitor, ptcon)[maxnpt];
extern "C" char TINKER_MOD(kpitor, kpt)[maxnpt][8];

double (&ptcon)[maxnpt] = TINKER_MOD(kpitor, ptcon);
char (&kpt)[maxnpt][8] = TINKER_MOD(kpitor, kpt);
#endif
} }
