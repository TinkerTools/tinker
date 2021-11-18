#pragma once

#include "macro.hh"

namespace tinker { namespace pitors {
extern int& npitors;
extern int*& ipit;
extern double*& kpit;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(pitors, npitors);
extern "C" int* TINKER_MOD(pitors, ipit);
extern "C" double* TINKER_MOD(pitors, kpit);

int& npitors = TINKER_MOD(pitors, npitors);
int*& ipit = TINKER_MOD(pitors, ipit);
double*& kpit = TINKER_MOD(pitors, kpit);
#endif
} }
