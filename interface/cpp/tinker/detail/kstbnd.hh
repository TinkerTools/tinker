#pragma once

#include "macro.hh"

namespace tinker { namespace kstbnd {
const int maxnsb = 2000;
extern double (&stbn)[maxnsb][2];
extern char (&ksb)[maxnsb][12];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(kstbnd, stbn)[maxnsb][2];
extern "C" char TINKER_MOD(kstbnd, ksb)[maxnsb][12];

double (&stbn)[maxnsb][2] = TINKER_MOD(kstbnd, stbn);
char (&ksb)[maxnsb][12] = TINKER_MOD(kstbnd, ksb);
#endif
} }
