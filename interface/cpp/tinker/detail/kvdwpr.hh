#pragma once

#include "macro.hh"

namespace tinker { namespace kvdwpr {
extern int& maxnvp;
extern double*& radpr;
extern double*& epspr;
extern char (*&kvpr)[8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kvdwpr, maxnvp);
extern "C" double* TINKER_MOD(kvdwpr, radpr);
extern "C" double* TINKER_MOD(kvdwpr, epspr);
extern "C" char (*TINKER_MOD(kvdwpr, kvpr))[8];

int& maxnvp = TINKER_MOD(kvdwpr, maxnvp);
double*& radpr = TINKER_MOD(kvdwpr, radpr);
double*& epspr = TINKER_MOD(kvdwpr, epspr);
char (*&kvpr)[8] = TINKER_MOD(kvdwpr, kvpr);
#endif
} }
