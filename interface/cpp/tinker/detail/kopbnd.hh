#pragma once

#include "macro.hh"

namespace tinker { namespace kopbnd {
extern int& maxnopb;
extern double*& opbn;
extern char (*&kopb)[16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kopbnd, maxnopb);
extern "C" double* TINKER_MOD(kopbnd, opbn);
extern "C" char (*TINKER_MOD(kopbnd, kopb))[16];

int& maxnopb = TINKER_MOD(kopbnd, maxnopb);
double*& opbn = TINKER_MOD(kopbnd, opbn);
char (*&kopb)[16] = TINKER_MOD(kopbnd, kopb);
#endif
} }
