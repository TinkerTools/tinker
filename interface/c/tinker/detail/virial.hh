#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(virial, vir)[3][3];
extern int TINKER_MOD(virial, use_virial);
#ifdef __cplusplus
}
#endif
