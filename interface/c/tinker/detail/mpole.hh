#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxpole 13
extern int TINKER_MOD(mpole, npole);
extern int* TINKER_MOD(mpole, ipole);
extern int* TINKER_MOD(mpole, polsiz);
extern int* TINKER_MOD(mpole, pollist);
extern int* TINKER_MOD(mpole, zaxis);
extern int* TINKER_MOD(mpole, xaxis);
extern int* TINKER_MOD(mpole, yaxis);
extern double* TINKER_MOD(mpole, pole);
extern double* TINKER_MOD(mpole, rpole);
extern double* TINKER_MOD(mpole, mono0);
extern char (*TINKER_MOD(mpole, polaxe))[8];
#ifdef __cplusplus
}
#endif
