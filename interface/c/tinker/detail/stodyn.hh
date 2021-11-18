#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(stodyn, friction);
extern double* TINKER_MOD(stodyn, fgamma);
extern int TINKER_MOD(stodyn, use_sdarea);
#ifdef __cplusplus
}
#endif
