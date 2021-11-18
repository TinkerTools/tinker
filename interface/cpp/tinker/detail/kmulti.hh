#pragma once

#include "macro.hh"

namespace tinker { namespace kmulti {
const int maxnmp = 2000;
extern double (&multip)[maxnmp][13];
extern char (&mpaxis)[maxnmp][8];
extern char (&kmp)[maxnmp][16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(kmulti, multip)[maxnmp][13];
extern "C" char TINKER_MOD(kmulti, mpaxis)[maxnmp][8];
extern "C" char TINKER_MOD(kmulti, kmp)[maxnmp][16];

double (&multip)[maxnmp][13] = TINKER_MOD(kmulti, multip);
char (&mpaxis)[maxnmp][8] = TINKER_MOD(kmulti, mpaxis);
char (&kmp)[maxnmp][16] = TINKER_MOD(kmulti, kmp);
#endif
} }
