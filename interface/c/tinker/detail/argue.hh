#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxarg 20
extern int TINKER_MOD(argue, narg);
extern int TINKER_MOD(argue, listarg)[TINKER_MOD__maxarg+1];
extern char TINKER_MOD(argue, arg)[TINKER_MOD__maxarg+1][240];
#ifdef __cplusplus
}
#endif
