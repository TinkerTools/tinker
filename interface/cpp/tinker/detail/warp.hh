#pragma once

#include "macro.hh"

namespace tinker { namespace warp {
extern double& deform;
extern double& difft;
extern double& diffv;
extern double& diffc;
extern double*& m2;
extern int& use_smooth;
extern int& use_dem;
extern int& use_gda;
extern int& use_tophat;
extern int& use_stophat;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(warp, deform);
extern "C" double TINKER_MOD(warp, difft);
extern "C" double TINKER_MOD(warp, diffv);
extern "C" double TINKER_MOD(warp, diffc);
extern "C" double* TINKER_MOD(warp, m2);
extern "C" int TINKER_MOD(warp, use_smooth);
extern "C" int TINKER_MOD(warp, use_dem);
extern "C" int TINKER_MOD(warp, use_gda);
extern "C" int TINKER_MOD(warp, use_tophat);
extern "C" int TINKER_MOD(warp, use_stophat);

double& deform = TINKER_MOD(warp, deform);
double& difft = TINKER_MOD(warp, difft);
double& diffv = TINKER_MOD(warp, diffv);
double& diffc = TINKER_MOD(warp, diffc);
double*& m2 = TINKER_MOD(warp, m2);
int& use_smooth = TINKER_MOD(warp, use_smooth);
int& use_dem = TINKER_MOD(warp, use_dem);
int& use_gda = TINKER_MOD(warp, use_gda);
int& use_tophat = TINKER_MOD(warp, use_tophat);
int& use_stophat = TINKER_MOD(warp, use_stophat);
#endif
} }
