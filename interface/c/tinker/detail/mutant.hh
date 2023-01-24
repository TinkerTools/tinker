#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(mutant, nmut);
extern int TINKER_MOD(mutant, vcouple);
extern int* TINKER_MOD(mutant, imut);
extern int* TINKER_MOD(mutant, type0);
extern int* TINKER_MOD(mutant, class0);
extern int* TINKER_MOD(mutant, type1);
extern int* TINKER_MOD(mutant, class1);
extern double TINKER_MOD(mutant, lambda);
extern double TINKER_MOD(mutant, vlambda);
extern double TINKER_MOD(mutant, elambda);
extern double TINKER_MOD(mutant, tlambda);
extern double TINKER_MOD(mutant, scexp);
extern double TINKER_MOD(mutant, scalpha);
extern int* TINKER_MOD(mutant, mut);
#ifdef __cplusplus
}
#endif
