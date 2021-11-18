#pragma once

#include "macro.hh"
#include "sizes.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(atoms, n);
extern int TINKER_MOD(atoms, type)[TINKER_MOD__maxatm];
extern double TINKER_MOD(atoms, x)[TINKER_MOD__maxatm];
extern double TINKER_MOD(atoms, y)[TINKER_MOD__maxatm];
extern double TINKER_MOD(atoms, z)[TINKER_MOD__maxatm];
#ifdef __cplusplus
}
#endif
