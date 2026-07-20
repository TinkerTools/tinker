#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(mutant, nmut);
extern int TINKER_MOD(mutant, nmutb);
extern int TINKER_MOD(mutant, vcouple);
extern int* TINKER_MOD(mutant, imut);
extern int* TINKER_MOD(mutant, type0);
extern int* TINKER_MOD(mutant, class0);
extern int* TINKER_MOD(mutant, type1);
extern int* TINKER_MOD(mutant, class1);
extern int* TINKER_MOD(mutant, mutg);
extern double TINKER_MOD(mutant, lambda);
extern double TINKER_MOD(mutant, vlambda);
extern double TINKER_MOD(mutant, elambda);
extern double TINKER_MOD(mutant, tlambda);
extern double TINKER_MOD(mutant, scexp);
extern double TINKER_MOD(mutant, scalpha);
extern int TINKER_MOD(mutant, use_rel);
extern int TINKER_MOD(mutant, use_subsys);
extern int* TINKER_MOD(mutant, mut);
extern int* TINKER_MOD(mutant, subon);
#ifdef __cplusplus
}
#endif
