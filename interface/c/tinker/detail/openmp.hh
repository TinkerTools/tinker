#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(openmp, nproc);
extern int TINKER_MOD(openmp, nthread);
#ifdef __cplusplus
}
#endif
