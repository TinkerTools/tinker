#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(omega, nomega);
extern int* TINKER_MOD(omega, iomega);
extern int* TINKER_MOD(omega, zline);
extern double* TINKER_MOD(omega, dihed);
#ifdef __cplusplus
}
#endif
