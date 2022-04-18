#pragma once

#include "macro.hh"

namespace tinker { namespace kcflux {
extern int& maxncfb;
extern int& maxncfa;
extern double*& cflb;
extern double*& cfla;
extern double*& cflab;
extern char (*&kcfb)[8];
extern char (*&kcfa)[12];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kcflux, maxncfb);
extern "C" int TINKER_MOD(kcflux, maxncfa);
extern "C" double* TINKER_MOD(kcflux, cflb);
extern "C" double* TINKER_MOD(kcflux, cfla);
extern "C" double* TINKER_MOD(kcflux, cflab);
extern "C" char (*TINKER_MOD(kcflux, kcfb))[8];
extern "C" char (*TINKER_MOD(kcflux, kcfa))[12];

int& maxncfb = TINKER_MOD(kcflux, maxncfb);
int& maxncfa = TINKER_MOD(kcflux, maxncfa);
double*& cflb = TINKER_MOD(kcflux, cflb);
double*& cfla = TINKER_MOD(kcflux, cfla);
double*& cflab = TINKER_MOD(kcflux, cflab);
char (*&kcfb)[8] = TINKER_MOD(kcflux, kcfb);
char (*&kcfa)[12] = TINKER_MOD(kcflux, kcfa);
#endif
} }
