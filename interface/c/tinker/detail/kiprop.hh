#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(kiprop, maxndi);
extern double* TINKER_MOD(kiprop, dcon);
extern double* TINKER_MOD(kiprop, tdi);
extern char (*TINKER_MOD(kiprop, kdi))[16];
#ifdef __cplusplus
}
#endif
