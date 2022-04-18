#pragma once

#include "macro.hh"

namespace tinker { namespace kmulti {
extern int& maxnmp;
extern double*& multip;
extern char (*&mpaxis)[8];
extern char (*&kmp)[16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kmulti, maxnmp);
extern "C" double* TINKER_MOD(kmulti, multip);
extern "C" char (*TINKER_MOD(kmulti, mpaxis))[8];
extern "C" char (*TINKER_MOD(kmulti, kmp))[16];

int& maxnmp = TINKER_MOD(kmulti, maxnmp);
double*& multip = TINKER_MOD(kmulti, multip);
char (*&mpaxis)[8] = TINKER_MOD(kmulti, mpaxis);
char (*&kmp)[16] = TINKER_MOD(kmulti, kmp);
#endif
} }
