#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxopt 6
extern int TINKER_MOD(polopt, optorder);
extern int TINKER_MOD(polopt, optlevel);
extern double* TINKER_MOD(polopt, copt);
extern double* TINKER_MOD(polopt, copm);
extern double* TINKER_MOD(polopt, uopt);
extern double* TINKER_MOD(polopt, uoptp);
extern double* TINKER_MOD(polopt, uopts);
extern double* TINKER_MOD(polopt, uoptps);
extern double* TINKER_MOD(polopt, fopt);
extern double* TINKER_MOD(polopt, foptp);
#ifdef __cplusplus
}
#endif
