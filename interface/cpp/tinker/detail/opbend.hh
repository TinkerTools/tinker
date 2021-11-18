#pragma once

#include "macro.hh"

namespace tinker { namespace opbend {
extern int& nopbend;
extern int*& iopb;
extern double*& opbk;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(opbend, nopbend);
extern "C" int* TINKER_MOD(opbend, iopb);
extern "C" double* TINKER_MOD(opbend, opbk);

int& nopbend = TINKER_MOD(opbend, nopbend);
int*& iopb = TINKER_MOD(opbend, iopb);
double*& opbk = TINKER_MOD(opbend, opbk);
#endif
} }
