#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double* TINKER_MOD(moldyn, v);
extern double* TINKER_MOD(moldyn, a);
extern double* TINKER_MOD(moldyn, aalt);
extern double* TINKER_MOD(moldyn, aslow);
extern double* TINKER_MOD(moldyn, afast);
#ifdef __cplusplus
}
#endif
