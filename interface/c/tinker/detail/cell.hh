#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(cell, ncell);
extern int* TINKER_MOD(cell, icell);
extern double TINKER_MOD(cell, xcell);
extern double TINKER_MOD(cell, ycell);
extern double TINKER_MOD(cell, zcell);
extern double TINKER_MOD(cell, xcell2);
extern double TINKER_MOD(cell, ycell2);
extern double TINKER_MOD(cell, zcell2);
#ifdef __cplusplus
}
#endif
