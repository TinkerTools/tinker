#pragma once

#include "macro.hh"

namespace tinker { namespace ktrtor {
const int maxntt = 100;
const int maxtgrd = 30;
const int maxtgrd2 = maxtgrd*maxtgrd;
extern int (&tnx)[maxntt];
extern int (&tny)[maxntt];
extern double (&ttx)[maxntt][maxtgrd];
extern double (&tty)[maxntt][maxtgrd];
extern double (&tbf)[maxntt][maxtgrd2];
extern double (&tbx)[maxntt][maxtgrd2];
extern double (&tby)[maxntt][maxtgrd2];
extern double (&tbxy)[maxntt][maxtgrd2];
extern char (&ktt)[maxntt][20];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(ktrtor, tnx)[maxntt];
extern "C" int TINKER_MOD(ktrtor, tny)[maxntt];
extern "C" double TINKER_MOD(ktrtor, ttx)[maxntt][maxtgrd];
extern "C" double TINKER_MOD(ktrtor, tty)[maxntt][maxtgrd];
extern "C" double TINKER_MOD(ktrtor, tbf)[maxntt][maxtgrd2];
extern "C" double TINKER_MOD(ktrtor, tbx)[maxntt][maxtgrd2];
extern "C" double TINKER_MOD(ktrtor, tby)[maxntt][maxtgrd2];
extern "C" double TINKER_MOD(ktrtor, tbxy)[maxntt][maxtgrd2];
extern "C" char TINKER_MOD(ktrtor, ktt)[maxntt][20];

int (&tnx)[maxntt] = TINKER_MOD(ktrtor, tnx);
int (&tny)[maxntt] = TINKER_MOD(ktrtor, tny);
double (&ttx)[maxntt][maxtgrd] = TINKER_MOD(ktrtor, ttx);
double (&tty)[maxntt][maxtgrd] = TINKER_MOD(ktrtor, tty);
double (&tbf)[maxntt][maxtgrd2] = TINKER_MOD(ktrtor, tbf);
double (&tbx)[maxntt][maxtgrd2] = TINKER_MOD(ktrtor, tbx);
double (&tby)[maxntt][maxtgrd2] = TINKER_MOD(ktrtor, tby);
double (&tbxy)[maxntt][maxtgrd2] = TINKER_MOD(ktrtor, tbxy);
char (&ktt)[maxntt][20] = TINKER_MOD(ktrtor, ktt);
#endif
} }
