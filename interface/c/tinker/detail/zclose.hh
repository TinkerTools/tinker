#pragma once

#include "macro.hh"
#include "sizes.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(zclose, nadd);
extern int TINKER_MOD(zclose, ndel);
extern int TINKER_MOD(zclose, iadd)[TINKER_MOD__maxatm][2];
extern int TINKER_MOD(zclose, idel)[TINKER_MOD__maxatm][2];
#ifdef __cplusplus
}
#endif
