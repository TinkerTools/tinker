#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(kdipol, maxnd);
extern int TINKER_MOD(kdipol, maxnd5);
extern int TINKER_MOD(kdipol, maxnd4);
extern int TINKER_MOD(kdipol, maxnd3);
extern double* TINKER_MOD(kdipol, dpl);
extern double* TINKER_MOD(kdipol, dpl5);
extern double* TINKER_MOD(kdipol, dpl4);
extern double* TINKER_MOD(kdipol, dpl3);
extern double* TINKER_MOD(kdipol, pos);
extern double* TINKER_MOD(kdipol, pos5);
extern double* TINKER_MOD(kdipol, pos4);
extern double* TINKER_MOD(kdipol, pos3);
extern char (*TINKER_MOD(kdipol, kd))[8];
extern char (*TINKER_MOD(kdipol, kd5))[8];
extern char (*TINKER_MOD(kdipol, kd4))[8];
extern char (*TINKER_MOD(kdipol, kd3))[8];
#ifdef __cplusplus
}
#endif
