#pragma once

#include "macro.hh"

namespace tinker { namespace dsppot {
extern double& dsp2scale;
extern double& dsp3scale;
extern double& dsp4scale;
extern double& dsp5scale;
extern int& use_dcorr;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(dsppot, dsp2scale);
extern "C" double TINKER_MOD(dsppot, dsp3scale);
extern "C" double TINKER_MOD(dsppot, dsp4scale);
extern "C" double TINKER_MOD(dsppot, dsp5scale);
extern "C" int TINKER_MOD(dsppot, use_dcorr);

double& dsp2scale = TINKER_MOD(dsppot, dsp2scale);
double& dsp3scale = TINKER_MOD(dsppot, dsp3scale);
double& dsp4scale = TINKER_MOD(dsppot, dsp4scale);
double& dsp5scale = TINKER_MOD(dsppot, dsp5scale);
int& use_dcorr = TINKER_MOD(dsppot, use_dcorr);
#endif
} }
