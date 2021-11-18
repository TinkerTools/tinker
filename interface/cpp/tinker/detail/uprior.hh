#pragma once

#include "macro.hh"

namespace tinker { namespace uprior {
const int maxpred = 17;
extern int& nualt;
extern int& maxualt;
extern double (&gear)[maxpred];
extern double (&aspc)[maxpred];
extern double (&bpred)[maxpred];
extern double (&bpredp)[maxpred];
extern double (&bpreds)[maxpred];
extern double (&bpredps)[maxpred];
extern double*& udalt;
extern double*& upalt;
extern double*& usalt;
extern double*& upsalt;
extern int& use_pred;
extern char (&polpred)[4];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(uprior, nualt);
extern "C" int TINKER_MOD(uprior, maxualt);
extern "C" double TINKER_MOD(uprior, gear)[maxpred];
extern "C" double TINKER_MOD(uprior, aspc)[maxpred];
extern "C" double TINKER_MOD(uprior, bpred)[maxpred];
extern "C" double TINKER_MOD(uprior, bpredp)[maxpred];
extern "C" double TINKER_MOD(uprior, bpreds)[maxpred];
extern "C" double TINKER_MOD(uprior, bpredps)[maxpred];
extern "C" double* TINKER_MOD(uprior, udalt);
extern "C" double* TINKER_MOD(uprior, upalt);
extern "C" double* TINKER_MOD(uprior, usalt);
extern "C" double* TINKER_MOD(uprior, upsalt);
extern "C" int TINKER_MOD(uprior, use_pred);
extern "C" char TINKER_MOD(uprior, polpred)[4];

int& nualt = TINKER_MOD(uprior, nualt);
int& maxualt = TINKER_MOD(uprior, maxualt);
double (&gear)[maxpred] = TINKER_MOD(uprior, gear);
double (&aspc)[maxpred] = TINKER_MOD(uprior, aspc);
double (&bpred)[maxpred] = TINKER_MOD(uprior, bpred);
double (&bpredp)[maxpred] = TINKER_MOD(uprior, bpredp);
double (&bpreds)[maxpred] = TINKER_MOD(uprior, bpreds);
double (&bpredps)[maxpred] = TINKER_MOD(uprior, bpredps);
double*& udalt = TINKER_MOD(uprior, udalt);
double*& upalt = TINKER_MOD(uprior, upalt);
double*& usalt = TINKER_MOD(uprior, usalt);
double*& upsalt = TINKER_MOD(uprior, upsalt);
int& use_pred = TINKER_MOD(uprior, use_pred);
char (&polpred)[4] = TINKER_MOD(uprior, polpred);
#endif
} }
