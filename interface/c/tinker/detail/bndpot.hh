#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(bndpot, cbnd);
extern double TINKER_MOD(bndpot, qbnd);
extern double TINKER_MOD(bndpot, bndunit);
extern char TINKER_MOD(bndpot, bndtyp)[8];
#ifdef __cplusplus
}
#endif
