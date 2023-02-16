#pragma once

#include "macro.hh"

namespace tinker { namespace extfld {
extern double (&exfld)[3];
extern int& use_exfld;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(extfld, exfld)[3];
extern "C" int TINKER_MOD(extfld, use_exfld);

double (&exfld)[3] = TINKER_MOD(extfld, exfld);
int& use_exfld = TINKER_MOD(extfld, use_exfld);
#endif
} }
