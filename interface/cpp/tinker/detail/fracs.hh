#pragma once

#include "macro.hh"

namespace tinker { namespace fracs {
extern double*& xfrac;
extern double*& yfrac;
extern double*& zfrac;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(fracs, xfrac);
extern "C" double* TINKER_MOD(fracs, yfrac);
extern "C" double* TINKER_MOD(fracs, zfrac);

double*& xfrac = TINKER_MOD(fracs, xfrac);
double*& yfrac = TINKER_MOD(fracs, yfrac);
double*& zfrac = TINKER_MOD(fracs, zfrac);
#endif
} }
