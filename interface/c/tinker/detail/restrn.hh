#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(restrn, npfix);
extern int TINKER_MOD(restrn, ndfix);
extern int TINKER_MOD(restrn, nafix);
extern int TINKER_MOD(restrn, ntfix);
extern int TINKER_MOD(restrn, ngfix);
extern int TINKER_MOD(restrn, nchir);
extern int* TINKER_MOD(restrn, ipfix);
extern int* TINKER_MOD(restrn, kpfix);
extern int* TINKER_MOD(restrn, idfix);
extern int* TINKER_MOD(restrn, iafix);
extern int* TINKER_MOD(restrn, itfix);
extern int* TINKER_MOD(restrn, igfix);
extern int* TINKER_MOD(restrn, ichir);
extern double TINKER_MOD(restrn, depth);
extern double TINKER_MOD(restrn, width);
extern double TINKER_MOD(restrn, rwall);
extern double* TINKER_MOD(restrn, xpfix);
extern double* TINKER_MOD(restrn, ypfix);
extern double* TINKER_MOD(restrn, zpfix);
extern double* TINKER_MOD(restrn, pfix);
extern double* TINKER_MOD(restrn, dfix);
extern double* TINKER_MOD(restrn, afix);
extern double* TINKER_MOD(restrn, tfix);
extern double* TINKER_MOD(restrn, gfix);
extern double* TINKER_MOD(restrn, chir);
extern int TINKER_MOD(restrn, use_basin);
extern int TINKER_MOD(restrn, use_wall);
#ifdef __cplusplus
}
#endif
