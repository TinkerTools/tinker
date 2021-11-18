#pragma once

#include "macro.hh"
#include "sizes.hh"

namespace tinker { namespace gkstuf {
using namespace sizes;

extern double& gkc;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(gkstuf, gkc);

double& gkc = TINKER_MOD(gkstuf, gkc);
#endif
} }
