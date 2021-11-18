#pragma once

#include "macro.hh"

namespace tinker { namespace scales {
extern double*& scale;
extern int& set_scale;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(scales, scale);
extern "C" int TINKER_MOD(scales, set_scale);

double*& scale = TINKER_MOD(scales, scale);
int& set_scale = TINKER_MOD(scales, set_scale);
#endif
} }
