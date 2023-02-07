#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(repel, nrep);
extern int* TINKER_MOD(repel, irep);
extern double* TINKER_MOD(repel, sizpr);
extern double* TINKER_MOD(repel, dmppr);
extern double* TINKER_MOD(repel, elepr);
extern double* TINKER_MOD(repel, repole);
extern double* TINKER_MOD(repel, rrepole);
#ifdef __cplusplus
}
#endif
