#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxprm 25000
extern int TINKER_MOD(params, nprm);
extern char TINKER_MOD(params, prmline)[TINKER_MOD__maxprm][240];
#ifdef __cplusplus
}
#endif
