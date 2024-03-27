#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(extfld, exfreq);
extern double TINKER_MOD(extfld, exfld)[3];
extern double TINKER_MOD(extfld, texfld)[3];
extern int TINKER_MOD(extfld, use_exfld);
extern int TINKER_MOD(extfld, use_exfreq);
#ifdef __cplusplus
}
#endif
