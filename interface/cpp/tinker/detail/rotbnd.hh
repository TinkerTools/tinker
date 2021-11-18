#pragma once

#include "macro.hh"

namespace tinker { namespace rotbnd {
extern int& nrot;
extern int*& rot;
extern int& use_short;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(rotbnd, nrot);
extern "C" int* TINKER_MOD(rotbnd, rot);
extern "C" int TINKER_MOD(rotbnd, use_short);

int& nrot = TINKER_MOD(rotbnd, nrot);
int*& rot = TINKER_MOD(rotbnd, rot);
int& use_short = TINKER_MOD(rotbnd, use_short);
#endif
} }
