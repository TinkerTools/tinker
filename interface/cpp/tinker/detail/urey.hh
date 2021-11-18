#pragma once

#include "macro.hh"

namespace tinker { namespace urey {
extern int& nurey;
extern int*& iury;
extern double*& uk;
extern double*& ul;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(urey, nurey);
extern "C" int* TINKER_MOD(urey, iury);
extern "C" double* TINKER_MOD(urey, uk);
extern "C" double* TINKER_MOD(urey, ul);

int& nurey = TINKER_MOD(urey, nurey);
int*& iury = TINKER_MOD(urey, iury);
double*& uk = TINKER_MOD(urey, uk);
double*& ul = TINKER_MOD(urey, ul);
#endif
} }
