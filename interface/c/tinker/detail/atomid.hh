#pragma once

#include "macro.hh"
#include "sizes.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(atomid, tag)[TINKER_MOD__maxatm];
extern int TINKER_MOD(atomid, class)[TINKER_MOD__maxatm];
extern int TINKER_MOD(atomid, atomic)[TINKER_MOD__maxatm];
extern int TINKER_MOD(atomid, valence)[TINKER_MOD__maxatm];
extern double TINKER_MOD(atomid, mass)[TINKER_MOD__maxatm];
extern char TINKER_MOD(atomid, name)[TINKER_MOD__maxatm][3];
extern char TINKER_MOD(atomid, story)[TINKER_MOD__maxatm][24];
#ifdef __cplusplus
}
#endif
