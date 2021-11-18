#pragma once

#include "macro.hh"

namespace tinker { namespace charge {
extern int& nion;
extern int*& iion;
extern int*& jion;
extern int*& kion;
extern double*& pchg;
extern double*& pchg0;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(charge, nion);
extern "C" int* TINKER_MOD(charge, iion);
extern "C" int* TINKER_MOD(charge, jion);
extern "C" int* TINKER_MOD(charge, kion);
extern "C" double* TINKER_MOD(charge, pchg);
extern "C" double* TINKER_MOD(charge, pchg0);

int& nion = TINKER_MOD(charge, nion);
int*& iion = TINKER_MOD(charge, iion);
int*& jion = TINKER_MOD(charge, jion);
int*& kion = TINKER_MOD(charge, kion);
double*& pchg = TINKER_MOD(charge, pchg);
double*& pchg0 = TINKER_MOD(charge, pchg0);
#endif
} }
