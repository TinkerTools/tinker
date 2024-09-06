#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(alfmol, alfthread);
extern double TINKER_MOD(alfmol, delcxeps);
extern int TINKER_MOD(alfmol, alfhydro);
extern int TINKER_MOD(alfmol, alfsosgmp);
extern char TINKER_MOD(alfmol, alfmethod)[6];
extern char TINKER_MOD(alfmol, alfsort)[6];
#ifdef __cplusplus
}
#endif
