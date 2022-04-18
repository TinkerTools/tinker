#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(ktrtor, maxntt);
extern int TINKER_MOD(ktrtor, maxtgrd);
extern int TINKER_MOD(ktrtor, maxtgrd2);
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
