#pragma once

#include "macro.hh"

namespace tinker { namespace ktrtor {
extern int& maxntt;
const int maxtgrd = 30;
const int maxtgrd2 = maxtgrd*maxtgrd;
extern int*& tnx;
extern int*& tny;
extern double*& ttx;
extern double*& tty;
extern double*& tbf;
extern double*& tbx;
extern double*& tby;
extern double*& tbxy;
extern char (*&ktt)[20];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(ktrtor, maxntt);
extern "C" int* TINKER_MOD(ktrtor, tnx);
extern "C" int* TINKER_MOD(ktrtor, tny);
extern "C" double* TINKER_MOD(ktrtor, ttx);
extern "C" double* TINKER_MOD(ktrtor, tty);
extern "C" double* TINKER_MOD(ktrtor, tbf);
extern "C" double* TINKER_MOD(ktrtor, tbx);
extern "C" double* TINKER_MOD(ktrtor, tby);
extern "C" double* TINKER_MOD(ktrtor, tbxy);
extern "C" char (*TINKER_MOD(ktrtor, ktt))[20];

int& maxntt = TINKER_MOD(ktrtor, maxntt);
int*& tnx = TINKER_MOD(ktrtor, tnx);
int*& tny = TINKER_MOD(ktrtor, tny);
double*& ttx = TINKER_MOD(ktrtor, ttx);
double*& tty = TINKER_MOD(ktrtor, tty);
double*& tbf = TINKER_MOD(ktrtor, tbf);
double*& tbx = TINKER_MOD(ktrtor, tbx);
double*& tby = TINKER_MOD(ktrtor, tby);
double*& tbxy = TINKER_MOD(ktrtor, tbxy);
char (*&ktt)[20] = TINKER_MOD(ktrtor, ktt);
#endif
} }
