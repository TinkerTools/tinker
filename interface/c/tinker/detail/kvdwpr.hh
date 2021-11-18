#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxnvp 500
extern double TINKER_MOD(kvdwpr, radpr)[TINKER_MOD__maxnvp];
extern double TINKER_MOD(kvdwpr, epspr)[TINKER_MOD__maxnvp];
extern char TINKER_MOD(kvdwpr, kvpr)[TINKER_MOD__maxnvp][8];
#ifdef __cplusplus
}
#endif
