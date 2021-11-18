#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxnhb 500
extern double TINKER_MOD(khbond, radhb)[TINKER_MOD__maxnhb];
extern double TINKER_MOD(khbond, epshb)[TINKER_MOD__maxnhb];
extern char TINKER_MOD(khbond, khb)[TINKER_MOD__maxnhb][8];
#ifdef __cplusplus
}
#endif
