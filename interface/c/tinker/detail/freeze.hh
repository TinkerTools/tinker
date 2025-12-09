#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(freeze, nrat);
extern int TINKER_MOD(freeze, nratx);
extern int TINKER_MOD(freeze, nwat);
extern int TINKER_MOD(freeze, nwat4);
extern int* TINKER_MOD(freeze, iratx);
extern int* TINKER_MOD(freeze, kratx);
extern int* TINKER_MOD(freeze, irat);
extern int* TINKER_MOD(freeze, iwat);
extern int* TINKER_MOD(freeze, iwat4);
extern double TINKER_MOD(freeze, rateps);
extern double* TINKER_MOD(freeze, krat);
extern double* TINKER_MOD(freeze, kwat);
extern double* TINKER_MOD(freeze, kwat4);
extern int TINKER_MOD(freeze, use_freeze);
extern int* TINKER_MOD(freeze, frzimage);
#ifdef __cplusplus
}
#endif
