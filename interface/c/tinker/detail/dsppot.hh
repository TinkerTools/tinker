#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(dsppot, dsp2scale);
extern double TINKER_MOD(dsppot, dsp3scale);
extern double TINKER_MOD(dsppot, dsp4scale);
extern double TINKER_MOD(dsppot, dsp5scale);
extern double TINKER_MOD(dsppot, dhcut);
extern double TINKER_MOD(dsppot, dhoff);
extern int TINKER_MOD(dsppot, use_dcorr);
extern int TINKER_MOD(dsppot, use_dhal);
extern int TINKER_MOD(dsppot, dhdamp);
#ifdef __cplusplus
}
#endif
