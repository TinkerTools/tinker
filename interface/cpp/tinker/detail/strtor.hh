#pragma once

#include "macro.hh"

namespace tinker { namespace strtor {
extern int& nstrtor;
extern int*& ist;
extern double*& kst;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(strtor, nstrtor);
extern "C" int* TINKER_MOD(strtor, ist);
extern "C" double* TINKER_MOD(strtor, kst);

int& nstrtor = TINKER_MOD(strtor, nstrtor);
int*& ist = TINKER_MOD(strtor, ist);
double*& kst = TINKER_MOD(strtor, kst);
#endif
} }
