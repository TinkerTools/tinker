#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(kmulti, maxnmp);
extern double* TINKER_MOD(kmulti, multip);
extern char (*TINKER_MOD(kmulti, mpaxis))[8];
extern char (*TINKER_MOD(kmulti, kmp))[16];
#ifdef __cplusplus
}
#endif
