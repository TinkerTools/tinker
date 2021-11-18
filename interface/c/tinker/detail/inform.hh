#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxask 5
extern int TINKER_MOD(inform, digits);
extern int TINKER_MOD(inform, iprint);
extern int TINKER_MOD(inform, iwrite);
extern int TINKER_MOD(inform, isend);
extern int TINKER_MOD(inform, verbose);
extern int TINKER_MOD(inform, debug);
extern int TINKER_MOD(inform, silent);
extern int TINKER_MOD(inform, holdup);
extern int TINKER_MOD(inform, abort);
#ifdef __cplusplus
}
#endif
