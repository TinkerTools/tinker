#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(kvdwpr, maxnvp);
extern double* TINKER_MOD(kvdwpr, radpr);
extern double* TINKER_MOD(kvdwpr, epspr);
extern char (*TINKER_MOD(kvdwpr, kvpr))[8];
#ifdef __cplusplus
}
#endif
