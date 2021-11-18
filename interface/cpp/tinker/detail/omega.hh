#pragma once

#include "macro.hh"

namespace tinker { namespace omega {
extern int& nomega;
extern int*& iomega;
extern int*& zline;
extern double*& dihed;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(omega, nomega);
extern "C" int* TINKER_MOD(omega, iomega);
extern "C" int* TINKER_MOD(omega, zline);
extern "C" double* TINKER_MOD(omega, dihed);

int& nomega = TINKER_MOD(omega, nomega);
int*& iomega = TINKER_MOD(omega, iomega);
int*& zline = TINKER_MOD(omega, zline);
double*& dihed = TINKER_MOD(omega, dihed);
#endif
} }
