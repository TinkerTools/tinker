#pragma once

#include "macro.hh"

namespace tinker { namespace extfld {
extern double& exfreq;
extern double (&exfld)[3];
extern double (&texfld)[3];
extern int& use_exfld;
extern int& use_exfreq;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(extfld, exfreq);
extern "C" double TINKER_MOD(extfld, exfld)[3];
extern "C" double TINKER_MOD(extfld, texfld)[3];
extern "C" int TINKER_MOD(extfld, use_exfld);
extern "C" int TINKER_MOD(extfld, use_exfreq);

double& exfreq = TINKER_MOD(extfld, exfreq);
double (&exfld)[3] = TINKER_MOD(extfld, exfld);
double (&texfld)[3] = TINKER_MOD(extfld, texfld);
int& use_exfld = TINKER_MOD(extfld, use_exfld);
int& use_exfreq = TINKER_MOD(extfld, use_exfreq);
#endif
} }
