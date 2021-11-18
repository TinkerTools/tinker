#pragma once

#include "macro.hh"

namespace tinker { namespace pistuf {
extern double*& bkpi;
extern double*& blpi;
extern double*& kslope;
extern double*& lslope;
extern double*& torsp2;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(pistuf, bkpi);
extern "C" double* TINKER_MOD(pistuf, blpi);
extern "C" double* TINKER_MOD(pistuf, kslope);
extern "C" double* TINKER_MOD(pistuf, lslope);
extern "C" double* TINKER_MOD(pistuf, torsp2);

double*& bkpi = TINKER_MOD(pistuf, bkpi);
double*& blpi = TINKER_MOD(pistuf, blpi);
double*& kslope = TINKER_MOD(pistuf, kslope);
double*& lslope = TINKER_MOD(pistuf, lslope);
double*& torsp2 = TINKER_MOD(pistuf, torsp2);
#endif
} }
