#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(piorbs, norbit);
extern int TINKER_MOD(piorbs, nconj);
extern int TINKER_MOD(piorbs, reorbit);
extern int TINKER_MOD(piorbs, nbpi);
extern int TINKER_MOD(piorbs, ntpi);
extern int* TINKER_MOD(piorbs, iorbit);
extern int* TINKER_MOD(piorbs, iconj);
extern int* TINKER_MOD(piorbs, kconj);
extern int* TINKER_MOD(piorbs, piperp);
extern int* TINKER_MOD(piorbs, ibpi);
extern int* TINKER_MOD(piorbs, itpi);
extern double* TINKER_MOD(piorbs, pbpl);
extern double* TINKER_MOD(piorbs, pnpl);
extern int* TINKER_MOD(piorbs, listpi);
#ifdef __cplusplus
}
#endif
