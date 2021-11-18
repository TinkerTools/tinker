#pragma once

#include "macro.hh"

namespace tinker { namespace kiprop {
const int maxndi = 500;
extern double (&dcon)[maxndi];
extern double (&tdi)[maxndi];
extern char (&kdi)[maxndi][16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(kiprop, dcon)[maxndi];
extern "C" double TINKER_MOD(kiprop, tdi)[maxndi];
extern "C" char TINKER_MOD(kiprop, kdi)[maxndi][16];

double (&dcon)[maxndi] = TINKER_MOD(kiprop, dcon);
double (&tdi)[maxndi] = TINKER_MOD(kiprop, tdi);
char (&kdi)[maxndi][16] = TINKER_MOD(kiprop, kdi);
#endif
} }
