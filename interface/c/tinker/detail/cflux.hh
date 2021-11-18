#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(cflux, nbflx);
extern int TINKER_MOD(cflux, naflx);
extern double* TINKER_MOD(cflux, bflx);
extern double* TINKER_MOD(cflux, aflx);
extern double* TINKER_MOD(cflux, abflx);
#ifdef __cplusplus
}
#endif
