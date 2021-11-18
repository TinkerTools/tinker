#pragma once

#include "macro.hh"

namespace tinker { namespace kopdst {
const int maxnopd = 500;
extern double (&opds)[maxnopd];
extern char (&kopd)[maxnopd][16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(kopdst, opds)[maxnopd];
extern "C" char TINKER_MOD(kopdst, kopd)[maxnopd][16];

double (&opds)[maxnopd] = TINKER_MOD(kopdst, opds);
char (&kopd)[maxnopd][16] = TINKER_MOD(kopdst, kopd);
#endif
} }
