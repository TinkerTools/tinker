#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(ewald, aewald);
extern double TINKER_MOD(ewald, aeewald);
extern double TINKER_MOD(ewald, apewald);
extern double TINKER_MOD(ewald, adewald);
extern char TINKER_MOD(ewald, boundary)[7];
#ifdef __cplusplus
}
#endif
