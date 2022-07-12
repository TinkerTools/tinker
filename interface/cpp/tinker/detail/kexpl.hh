#pragma once

#include "macro.hh"

namespace tinker { namespace kexpl {
extern double*& pepk;
extern double*& peppre;
extern double*& pepdmp;
extern int*& pepl;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(kexpl, pepk);
extern "C" double* TINKER_MOD(kexpl, peppre);
extern "C" double* TINKER_MOD(kexpl, pepdmp);
extern "C" int* TINKER_MOD(kexpl, pepl);

double*& pepk = TINKER_MOD(kexpl, pepk);
double*& peppre = TINKER_MOD(kexpl, peppre);
double*& pepdmp = TINKER_MOD(kexpl, pepdmp);
int*& pepl = TINKER_MOD(kexpl, pepl);
#endif
} }
