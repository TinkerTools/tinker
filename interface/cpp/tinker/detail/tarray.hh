#pragma once

#include "macro.hh"

namespace tinker { namespace tarray {
extern int& ntpair;
extern int*& tindex;
extern double*& tdipdip;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(tarray, ntpair);
extern "C" int* TINKER_MOD(tarray, tindex);
extern "C" double* TINKER_MOD(tarray, tdipdip);

int& ntpair = TINKER_MOD(tarray, ntpair);
int*& tindex = TINKER_MOD(tarray, tindex);
double*& tdipdip = TINKER_MOD(tarray, tdipdip);
#endif
} }
