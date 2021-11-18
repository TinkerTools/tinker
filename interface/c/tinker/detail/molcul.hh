#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(molcul, nmol);
extern int* TINKER_MOD(molcul, imol);
extern int* TINKER_MOD(molcul, kmol);
extern int* TINKER_MOD(molcul, molcule);
extern double TINKER_MOD(molcul, totmass);
extern double* TINKER_MOD(molcul, molmass);
#ifdef __cplusplus
}
#endif
