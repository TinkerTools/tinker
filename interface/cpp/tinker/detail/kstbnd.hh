#pragma once

#include "macro.hh"

namespace tinker { namespace kstbnd {
extern int& maxnsb;
extern double*& stbn;
extern char (*&ksb)[12];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kstbnd, maxnsb);
extern "C" double* TINKER_MOD(kstbnd, stbn);
extern "C" char (*TINKER_MOD(kstbnd, ksb))[12];

int& maxnsb = TINKER_MOD(kstbnd, maxnsb);
double*& stbn = TINKER_MOD(kstbnd, stbn);
char (*&ksb)[12] = TINKER_MOD(kstbnd, ksb);
#endif
} }
