#pragma once

#include "macro.hh"

namespace tinker { namespace korbs {
const int maxnpi = 500;
const int maxnpi5 = 200;
const int maxnpi4 = 200;
extern double (&sslope)[maxnpi];
extern double (&sslope5)[maxnpi5];
extern double (&sslope4)[maxnpi4];
extern double (&tslope)[maxnpi];
extern double (&tslope5)[maxnpi5];
extern double (&tslope4)[maxnpi4];
extern double*& electron;
extern double*& ionize;
extern double*& repulse;
extern char (&kpi)[maxnpi][8];
extern char (&kpi5)[maxnpi5][8];
extern char (&kpi4)[maxnpi4][8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(korbs, sslope)[maxnpi];
extern "C" double TINKER_MOD(korbs, sslope5)[maxnpi5];
extern "C" double TINKER_MOD(korbs, sslope4)[maxnpi4];
extern "C" double TINKER_MOD(korbs, tslope)[maxnpi];
extern "C" double TINKER_MOD(korbs, tslope5)[maxnpi5];
extern "C" double TINKER_MOD(korbs, tslope4)[maxnpi4];
extern "C" double* TINKER_MOD(korbs, electron);
extern "C" double* TINKER_MOD(korbs, ionize);
extern "C" double* TINKER_MOD(korbs, repulse);
extern "C" char TINKER_MOD(korbs, kpi)[maxnpi][8];
extern "C" char TINKER_MOD(korbs, kpi5)[maxnpi5][8];
extern "C" char TINKER_MOD(korbs, kpi4)[maxnpi4][8];

double (&sslope)[maxnpi] = TINKER_MOD(korbs, sslope);
double (&sslope5)[maxnpi5] = TINKER_MOD(korbs, sslope5);
double (&sslope4)[maxnpi4] = TINKER_MOD(korbs, sslope4);
double (&tslope)[maxnpi] = TINKER_MOD(korbs, tslope);
double (&tslope5)[maxnpi5] = TINKER_MOD(korbs, tslope5);
double (&tslope4)[maxnpi4] = TINKER_MOD(korbs, tslope4);
double*& electron = TINKER_MOD(korbs, electron);
double*& ionize = TINKER_MOD(korbs, ionize);
double*& repulse = TINKER_MOD(korbs, repulse);
char (&kpi)[maxnpi][8] = TINKER_MOD(korbs, kpi);
char (&kpi5)[maxnpi5][8] = TINKER_MOD(korbs, kpi5);
char (&kpi4)[maxnpi4][8] = TINKER_MOD(korbs, kpi4);
#endif
} }
