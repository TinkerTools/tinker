#pragma once

#include "macro.hh"

namespace tinker { namespace kopbnd {
const int maxnopb = 500;
extern double (&opbn)[maxnopb];
extern char (&kopb)[maxnopb][16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(kopbnd, opbn)[maxnopb];
extern "C" char TINKER_MOD(kopbnd, kopb)[maxnopb][16];

double (&opbn)[maxnopb] = TINKER_MOD(kopbnd, opbn);
char (&kopb)[maxnopb][16] = TINKER_MOD(kopbnd, kopb);
#endif
} }
