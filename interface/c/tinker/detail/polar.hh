#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(polar, npolar);
extern int* TINKER_MOD(polar, ipolar);
extern int* TINKER_MOD(polar, jpolar);
extern double* TINKER_MOD(polar, polarity);
extern double* TINKER_MOD(polar, thole);
extern double* TINKER_MOD(polar, tholed);
extern double* TINKER_MOD(polar, pdamp);
extern double* TINKER_MOD(polar, thlval);
extern double* TINKER_MOD(polar, thdval);
extern double* TINKER_MOD(polar, udir);
extern double* TINKER_MOD(polar, udirp);
extern double* TINKER_MOD(polar, udirs);
extern double* TINKER_MOD(polar, udirps);
extern double* TINKER_MOD(polar, uind);
extern double* TINKER_MOD(polar, uinp);
extern double* TINKER_MOD(polar, uinds);
extern double* TINKER_MOD(polar, uinps);
extern double* TINKER_MOD(polar, uexact);
extern int* TINKER_MOD(polar, douind);
#ifdef __cplusplus
}
#endif
