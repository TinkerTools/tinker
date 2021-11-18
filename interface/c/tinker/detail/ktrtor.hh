#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxntt 100
#define TINKER_MOD__maxtgrd 30
#define TINKER_MOD__maxtgrd2 maxtgrd*maxtgrd
extern int TINKER_MOD(ktrtor, tnx)[TINKER_MOD__maxntt];
extern int TINKER_MOD(ktrtor, tny)[TINKER_MOD__maxntt];
extern double TINKER_MOD(ktrtor, ttx)[TINKER_MOD__maxntt][TINKER_MOD__maxtgrd];
extern double TINKER_MOD(ktrtor, tty)[TINKER_MOD__maxntt][TINKER_MOD__maxtgrd];
extern double TINKER_MOD(ktrtor, tbf)[TINKER_MOD__maxntt][TINKER_MOD__maxtgrd2];
extern double TINKER_MOD(ktrtor, tbx)[TINKER_MOD__maxntt][TINKER_MOD__maxtgrd2];
extern double TINKER_MOD(ktrtor, tby)[TINKER_MOD__maxntt][TINKER_MOD__maxtgrd2];
extern double TINKER_MOD(ktrtor, tbxy)[TINKER_MOD__maxntt][TINKER_MOD__maxtgrd2];
extern char TINKER_MOD(ktrtor, ktt)[TINKER_MOD__maxntt][20];
#ifdef __cplusplus
}
#endif
