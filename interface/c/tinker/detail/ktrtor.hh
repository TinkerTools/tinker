#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(ktrtor, maxntt);
#define TINKER_MOD__maxtgrd 30
#define TINKER_MOD__maxtgrd2 (30*30)
extern int* TINKER_MOD(ktrtor, tnx);
extern int* TINKER_MOD(ktrtor, tny);
extern double* TINKER_MOD(ktrtor, ttx);
extern double* TINKER_MOD(ktrtor, tty);
extern double* TINKER_MOD(ktrtor, tbf);
extern double* TINKER_MOD(ktrtor, tbx);
extern double* TINKER_MOD(ktrtor, tby);
extern double* TINKER_MOD(ktrtor, tbxy);
extern char (*TINKER_MOD(ktrtor, ktt))[20];
#ifdef __cplusplus
}
#endif
