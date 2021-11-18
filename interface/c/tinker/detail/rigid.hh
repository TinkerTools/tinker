#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double* TINKER_MOD(rigid, xrb);
extern double* TINKER_MOD(rigid, yrb);
extern double* TINKER_MOD(rigid, zrb);
extern double* TINKER_MOD(rigid, rbc);
extern int TINKER_MOD(rigid, use_rigid);
#ifdef __cplusplus
}
#endif
