#pragma once

#include "macro.hh"

namespace tinker { namespace reppot {
extern double& r2scale;
extern double& r3scale;
extern double& r4scale;
extern double& r5scale;
extern int& reppolar;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(reppot, r2scale);
extern "C" double TINKER_MOD(reppot, r3scale);
extern "C" double TINKER_MOD(reppot, r4scale);
extern "C" double TINKER_MOD(reppot, r5scale);
extern "C" int TINKER_MOD(reppot, reppolar);

double& r2scale = TINKER_MOD(reppot, r2scale);
double& r3scale = TINKER_MOD(reppot, r3scale);
double& r4scale = TINKER_MOD(reppot, r4scale);
double& r5scale = TINKER_MOD(reppot, r5scale);
int& reppolar = TINKER_MOD(reppot, reppolar);
#endif
} }
