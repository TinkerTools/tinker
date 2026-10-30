#pragma once

#include "macro.hh"

namespace tinker { namespace dsppot {
extern double& dsp2scale;
extern double& dsp3scale;
extern double& dsp4scale;
extern double& dsp5scale;
extern double& dhcut;
extern double& dhoff;
extern int& use_dcorr;
extern int& use_dhal;
extern int& dhdamp;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(dsppot, dsp2scale);
extern "C" double TINKER_MOD(dsppot, dsp3scale);
extern "C" double TINKER_MOD(dsppot, dsp4scale);
extern "C" double TINKER_MOD(dsppot, dsp5scale);
extern "C" double TINKER_MOD(dsppot, dhcut);
extern "C" double TINKER_MOD(dsppot, dhoff);
extern "C" int TINKER_MOD(dsppot, use_dcorr);
extern "C" int TINKER_MOD(dsppot, use_dhal);
extern "C" int TINKER_MOD(dsppot, dhdamp);

double& dsp2scale = TINKER_MOD(dsppot, dsp2scale);
double& dsp3scale = TINKER_MOD(dsppot, dsp3scale);
double& dsp4scale = TINKER_MOD(dsppot, dsp4scale);
double& dsp5scale = TINKER_MOD(dsppot, dsp5scale);
double& dhcut = TINKER_MOD(dsppot, dhcut);
double& dhoff = TINKER_MOD(dsppot, dhoff);
int& use_dcorr = TINKER_MOD(dsppot, use_dcorr);
int& use_dhal = TINKER_MOD(dsppot, use_dhal);
int& dhdamp = TINKER_MOD(dsppot, dhdamp);
#endif
} }
