#pragma once

#include "macro.hh"
#include "sizes.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(phipsi, chiral)[TINKER_MOD__maxres];
extern int TINKER_MOD(phipsi, disulf)[TINKER_MOD__maxres];
extern double TINKER_MOD(phipsi, phi)[TINKER_MOD__maxres];
extern double TINKER_MOD(phipsi, psi)[TINKER_MOD__maxres];
extern double TINKER_MOD(phipsi, omg)[TINKER_MOD__maxres];
extern double TINKER_MOD(phipsi, chi)[TINKER_MOD__maxres][4];
#ifdef __cplusplus
}
#endif
