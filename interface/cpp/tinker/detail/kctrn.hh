#pragma once

#include "macro.hh"

namespace tinker { namespace kctrn {
extern double*& ctchg;
extern double*& ctdmp;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(kctrn, ctchg);
extern "C" double* TINKER_MOD(kctrn, ctdmp);

double*& ctchg = TINKER_MOD(kctrn, ctchg);
double*& ctdmp = TINKER_MOD(kctrn, ctdmp);
#endif
} }
