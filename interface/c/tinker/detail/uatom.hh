#pragma once

#include "macro.hh"
#include "sizes.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(uatom, nunique);
extern int TINKER_MOD(uatom, utype)[TINKER_MOD__maxtyp];
extern int TINKER_MOD(uatom, utypeinv)[TINKER_MOD__maxtyp];
extern double TINKER_MOD(uatom, utv1)[TINKER_MOD__maxtyp][3];
extern double TINKER_MOD(uatom, utv2)[TINKER_MOD__maxtyp][3];
#ifdef __cplusplus
}
#endif
