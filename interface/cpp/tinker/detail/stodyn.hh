#pragma once

#include "macro.hh"

namespace tinker { namespace stodyn {
extern double& friction;
extern double*& fgamma;
extern int& use_sdarea;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(stodyn, friction);
extern "C" double* TINKER_MOD(stodyn, fgamma);
extern "C" int TINKER_MOD(stodyn, use_sdarea);

double& friction = TINKER_MOD(stodyn, friction);
double*& fgamma = TINKER_MOD(stodyn, fgamma);
int& use_sdarea = TINKER_MOD(stodyn, use_sdarea);
#endif
} }
