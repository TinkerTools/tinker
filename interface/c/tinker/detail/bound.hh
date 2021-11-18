#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(bound, polycut);
extern double TINKER_MOD(bound, polycut2);
extern int TINKER_MOD(bound, use_bounds);
extern int TINKER_MOD(bound, use_replica);
extern int TINKER_MOD(bound, use_polymer);
#ifdef __cplusplus
}
#endif
