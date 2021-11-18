#pragma once

#include "macro.hh"

namespace tinker { namespace chgpen {
extern int& ncp;
extern double*& pcore;
extern double*& pval;
extern double*& pval0;
extern double*& palpha;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(chgpen, ncp);
extern "C" double* TINKER_MOD(chgpen, pcore);
extern "C" double* TINKER_MOD(chgpen, pval);
extern "C" double* TINKER_MOD(chgpen, pval0);
extern "C" double* TINKER_MOD(chgpen, palpha);

int& ncp = TINKER_MOD(chgpen, ncp);
double*& pcore = TINKER_MOD(chgpen, pcore);
double*& pval = TINKER_MOD(chgpen, pval);
double*& pval0 = TINKER_MOD(chgpen, pval0);
double*& palpha = TINKER_MOD(chgpen, palpha);
#endif
} }
