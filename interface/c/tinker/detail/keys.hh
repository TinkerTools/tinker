#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxkey 25000
extern int TINKER_MOD(keys, nkey);
extern char TINKER_MOD(keys, keyline)[TINKER_MOD__maxkey][240];
#ifdef __cplusplus
}
#endif
