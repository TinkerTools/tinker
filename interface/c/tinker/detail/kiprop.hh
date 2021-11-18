#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxndi 500
extern double TINKER_MOD(kiprop, dcon)[TINKER_MOD__maxndi];
extern double TINKER_MOD(kiprop, tdi)[TINKER_MOD__maxndi];
extern char TINKER_MOD(kiprop, kdi)[TINKER_MOD__maxndi][16];
#ifdef __cplusplus
}
#endif
