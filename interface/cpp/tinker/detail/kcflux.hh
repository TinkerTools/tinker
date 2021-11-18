#pragma once

#include "macro.hh"

namespace tinker { namespace kcflux {
const int maxncfb = 2000;
const int maxncfa = 2000;
extern double (&cflb)[maxncfb];
extern double (&cfla)[maxncfa][2];
extern double (&cflab)[maxncfa][2];
extern char (&kcfb)[maxncfb][8];
extern char (&kcfa)[maxncfa][12];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(kcflux, cflb)[maxncfb];
extern "C" double TINKER_MOD(kcflux, cfla)[maxncfa][2];
extern "C" double TINKER_MOD(kcflux, cflab)[maxncfa][2];
extern "C" char TINKER_MOD(kcflux, kcfb)[maxncfb][8];
extern "C" char TINKER_MOD(kcflux, kcfa)[maxncfa][12];

double (&cflb)[maxncfb] = TINKER_MOD(kcflux, cflb);
double (&cfla)[maxncfa][2] = TINKER_MOD(kcflux, cfla);
double (&cflab)[maxncfa][2] = TINKER_MOD(kcflux, cflab);
char (&kcfb)[maxncfb][8] = TINKER_MOD(kcflux, kcfb);
char (&kcfa)[maxncfa][12] = TINKER_MOD(kcflux, kcfa);
#endif
} }
