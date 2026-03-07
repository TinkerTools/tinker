#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(keys, nkey);
extern char (*TINKER_MOD(keys, keyline))[240];
#ifdef __cplusplus
}
#endif
