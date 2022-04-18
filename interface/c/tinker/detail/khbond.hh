#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(khbond, maxnhb);
extern double* TINKER_MOD(khbond, radhb);
extern double* TINKER_MOD(khbond, epshb);
extern char (*TINKER_MOD(khbond, khb))[8];
#ifdef __cplusplus
}
#endif
