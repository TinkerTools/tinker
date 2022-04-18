#pragma once

#include "macro.hh"

namespace tinker { namespace kiprop {
extern int& maxndi;
extern double*& dcon;
extern double*& tdi;
extern char (*&kdi)[16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kiprop, maxndi);
extern "C" double* TINKER_MOD(kiprop, dcon);
extern "C" double* TINKER_MOD(kiprop, tdi);
extern "C" char (*TINKER_MOD(kiprop, kdi))[16];

int& maxndi = TINKER_MOD(kiprop, maxndi);
double*& dcon = TINKER_MOD(kiprop, dcon);
double*& tdi = TINKER_MOD(kiprop, tdi);
char (*&kdi)[16] = TINKER_MOD(kiprop, kdi);
#endif
} }
