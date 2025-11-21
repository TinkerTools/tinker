#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(openmp, nproc);
extern int TINKER_MOD(openmp, nthread);
extern int TINKER_MOD(openmp, nnest);
#ifdef __cplusplus
}
#endif
