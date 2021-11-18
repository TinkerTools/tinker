#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(mdstuf, nfree);
extern int TINKER_MOD(mdstuf, irest);
extern int TINKER_MOD(mdstuf, bmnmix);
extern int TINKER_MOD(mdstuf, nrespa);
extern double TINKER_MOD(mdstuf, arespa);
extern int TINKER_MOD(mdstuf, dorest);
extern char TINKER_MOD(mdstuf, integrate)[11];
#ifdef __cplusplus
}
#endif
