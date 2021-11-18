#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(ielscf, nfree_aux);
extern double TINKER_MOD(ielscf, tautemp_aux);
extern double TINKER_MOD(ielscf, kelvin_aux);
extern double* TINKER_MOD(ielscf, uaux);
extern double* TINKER_MOD(ielscf, upaux);
extern double* TINKER_MOD(ielscf, vaux);
extern double* TINKER_MOD(ielscf, vpaux);
extern double* TINKER_MOD(ielscf, aaux);
extern double* TINKER_MOD(ielscf, apaux);
extern int TINKER_MOD(ielscf, use_ielscf);
#ifdef __cplusplus
}
#endif
