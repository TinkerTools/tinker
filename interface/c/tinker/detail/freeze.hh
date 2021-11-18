#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(freeze, nrat);
extern int TINKER_MOD(freeze, nratx);
extern int* TINKER_MOD(freeze, iratx);
extern int* TINKER_MOD(freeze, kratx);
extern int* TINKER_MOD(freeze, irat);
extern double TINKER_MOD(freeze, rateps);
extern double* TINKER_MOD(freeze, krat);
extern int TINKER_MOD(freeze, use_rattle);
extern int* TINKER_MOD(freeze, ratimage);
#ifdef __cplusplus
}
#endif
