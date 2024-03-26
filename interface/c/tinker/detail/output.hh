#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(output, nonly);
extern int* TINKER_MOD(output, ionly);
extern int* TINKER_MOD(output, ionlyinv);
extern int TINKER_MOD(output, archive);
extern int TINKER_MOD(output, binary);
extern int TINKER_MOD(output, noversion);
extern int TINKER_MOD(output, overwrite);
extern int TINKER_MOD(output, coordsave);
extern int TINKER_MOD(output, dynsave);
extern int TINKER_MOD(output, cyclesave);
extern int TINKER_MOD(output, onlysave);
extern int TINKER_MOD(output, arcsave);
extern int TINKER_MOD(output, dcdsave);
extern int TINKER_MOD(output, velsave);
extern int TINKER_MOD(output, frcsave);
extern int TINKER_MOD(output, uindsave);
extern int TINKER_MOD(output, ustcsave);
extern int TINKER_MOD(output, usyssave);
extern int* TINKER_MOD(output, usysuse);
extern char TINKER_MOD(output, coordtype)[9];
#ifdef __cplusplus
}
#endif
