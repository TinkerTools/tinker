#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxpred 17
extern int TINKER_MOD(uprior, nualt);
extern int TINKER_MOD(uprior, maxualt);
extern double TINKER_MOD(uprior, gear)[TINKER_MOD__maxpred];
extern double TINKER_MOD(uprior, aspc)[TINKER_MOD__maxpred];
extern double TINKER_MOD(uprior, bpred)[TINKER_MOD__maxpred];
extern double TINKER_MOD(uprior, bpredp)[TINKER_MOD__maxpred];
extern double TINKER_MOD(uprior, bpreds)[TINKER_MOD__maxpred];
extern double TINKER_MOD(uprior, bpredps)[TINKER_MOD__maxpred];
extern double* TINKER_MOD(uprior, udalt);
extern double* TINKER_MOD(uprior, upalt);
extern double* TINKER_MOD(uprior, usalt);
extern double* TINKER_MOD(uprior, upsalt);
extern int TINKER_MOD(uprior, use_pred);
extern char TINKER_MOD(uprior, polpred)[4];
#ifdef __cplusplus
}
#endif
