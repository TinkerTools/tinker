#pragma once

#include "macro.hh"

namespace tinker { namespace virial {
extern double (&vir)[3][3];
extern int& use_virial;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(virial, vir)[3][3];
extern "C" int TINKER_MOD(virial, use_virial);

double (&vir)[3][3] = TINKER_MOD(virial, vir);
int& use_virial = TINKER_MOD(virial, use_virial);
#endif
} }
