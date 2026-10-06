#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(fep, fepintv);
extern int TINKER_MOD(fep, nfep);
extern int TINKER_MOD(fep, nfepsave);
extern int TINKER_MOD(fep, nfeptot);
extern int* TINKER_MOD(fep, fepstep);
extern double* TINKER_MOD(fep, fepene);
extern double* TINKER_MOD(fep, feplmda);
extern double* TINKER_MOD(fep, feptrial);
#ifdef __cplusplus
}
#endif
