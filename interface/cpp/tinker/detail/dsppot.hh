#pragma once

#include "macro.hh"

namespace tinker { namespace dsppot {
extern double& dsp2scale;
extern double& dsp3scale;
extern double& dsp4scale;
extern double& dsp5scale;
extern double& dhcut;
extern double& dhoff;
extern double& dhc0;
extern double& dhc1;
extern double& dhc2;
extern double& dhc3;
extern double& dhc4;
extern double& dhc5;
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
extern "C" double TINKER_MOD(dsppot, dhc0);
extern "C" double TINKER_MOD(dsppot, dhc1);
extern "C" double TINKER_MOD(dsppot, dhc2);
extern "C" double TINKER_MOD(dsppot, dhc3);
extern "C" double TINKER_MOD(dsppot, dhc4);
extern "C" double TINKER_MOD(dsppot, dhc5);
extern "C" int TINKER_MOD(dsppot, use_dcorr);
extern "C" int TINKER_MOD(dsppot, use_dhal);
extern "C" int TINKER_MOD(dsppot, dhdamp);

double& dsp2scale = TINKER_MOD(dsppot, dsp2scale);
double& dsp3scale = TINKER_MOD(dsppot, dsp3scale);
double& dsp4scale = TINKER_MOD(dsppot, dsp4scale);
double& dsp5scale = TINKER_MOD(dsppot, dsp5scale);
double& dhcut = TINKER_MOD(dsppot, dhcut);
double& dhoff = TINKER_MOD(dsppot, dhoff);
double& dhc0 = TINKER_MOD(dsppot, dhc0);
double& dhc1 = TINKER_MOD(dsppot, dhc1);
double& dhc2 = TINKER_MOD(dsppot, dhc2);
double& dhc3 = TINKER_MOD(dsppot, dhc3);
double& dhc4 = TINKER_MOD(dsppot, dhc4);
double& dhc5 = TINKER_MOD(dsppot, dhc5);
int& use_dcorr = TINKER_MOD(dsppot, use_dcorr);
int& use_dhal = TINKER_MOD(dsppot, use_dhal);
int& dhdamp = TINKER_MOD(dsppot, dhdamp);
#endif
} }
