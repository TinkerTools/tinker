#pragma once

#include "macro.hh"

namespace tinker { namespace kopdst {
extern int& maxnopd;
extern double*& opds;
extern char (*&kopd)[16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kopdst, maxnopd);
extern "C" double* TINKER_MOD(kopdst, opds);
extern "C" char (*TINKER_MOD(kopdst, kopd))[16];

int& maxnopd = TINKER_MOD(kopdst, maxnopd);
double*& opds = TINKER_MOD(kopdst, opds);
char (*&kopd)[16] = TINKER_MOD(kopdst, kopd);
#endif
} }
