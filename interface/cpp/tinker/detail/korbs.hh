#pragma once

#include "macro.hh"

namespace tinker { namespace korbs {
extern int& maxnpi;
extern int& maxnpi5;
extern int& maxnpi4;
extern double*& electron;
extern double*& ionize;
extern double*& repulse;
extern double*& sslope;
extern double*& sslope5;
extern double*& sslope4;
extern double*& tslope;
extern double*& tslope5;
extern double*& tslope4;
extern char (*&kpi)[8];
extern char (*&kpi5)[8];
extern char (*&kpi4)[8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(korbs, maxnpi);
extern "C" int TINKER_MOD(korbs, maxnpi5);
extern "C" int TINKER_MOD(korbs, maxnpi4);
extern "C" double* TINKER_MOD(korbs, electron);
extern "C" double* TINKER_MOD(korbs, ionize);
extern "C" double* TINKER_MOD(korbs, repulse);
extern "C" double* TINKER_MOD(korbs, sslope);
extern "C" double* TINKER_MOD(korbs, sslope5);
extern "C" double* TINKER_MOD(korbs, sslope4);
extern "C" double* TINKER_MOD(korbs, tslope);
extern "C" double* TINKER_MOD(korbs, tslope5);
extern "C" double* TINKER_MOD(korbs, tslope4);
extern "C" char (*TINKER_MOD(korbs, kpi))[8];
extern "C" char (*TINKER_MOD(korbs, kpi5))[8];
extern "C" char (*TINKER_MOD(korbs, kpi4))[8];

int& maxnpi = TINKER_MOD(korbs, maxnpi);
int& maxnpi5 = TINKER_MOD(korbs, maxnpi5);
int& maxnpi4 = TINKER_MOD(korbs, maxnpi4);
double*& electron = TINKER_MOD(korbs, electron);
double*& ionize = TINKER_MOD(korbs, ionize);
double*& repulse = TINKER_MOD(korbs, repulse);
double*& sslope = TINKER_MOD(korbs, sslope);
double*& sslope5 = TINKER_MOD(korbs, sslope5);
double*& sslope4 = TINKER_MOD(korbs, sslope4);
double*& tslope = TINKER_MOD(korbs, tslope);
double*& tslope5 = TINKER_MOD(korbs, tslope5);
double*& tslope4 = TINKER_MOD(korbs, tslope4);
char (*&kpi)[8] = TINKER_MOD(korbs, kpi);
char (*&kpi5)[8] = TINKER_MOD(korbs, kpi5);
char (*&kpi4)[8] = TINKER_MOD(korbs, kpi4);
#endif
} }
