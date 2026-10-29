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
extern double TINKER_MOD(dsppot, dhc0);
extern double TINKER_MOD(dsppot, dhc1);
extern double TINKER_MOD(dsppot, dhc2);
extern double TINKER_MOD(dsppot, dhc3);
extern double TINKER_MOD(dsppot, dhc4);
extern double TINKER_MOD(dsppot, dhc5);
extern int TINKER_MOD(dsppot, use_dcorr);
extern int TINKER_MOD(dsppot, use_dhal);
extern int TINKER_MOD(dsppot, dhdamp);
#ifdef __cplusplus
}
#endif
