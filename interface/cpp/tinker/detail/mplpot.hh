#pragma once

#include "macro.hh"

namespace tinker { namespace mplpot {
extern double& m2scale;
extern double& m3scale;
extern double& m4scale;
extern double& m5scale;
extern int& use_chgpen;
extern char (&pentyp)[7];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(mplpot, m2scale);
extern "C" double TINKER_MOD(mplpot, m3scale);
extern "C" double TINKER_MOD(mplpot, m4scale);
extern "C" double TINKER_MOD(mplpot, m5scale);
extern "C" int TINKER_MOD(mplpot, use_chgpen);
extern "C" char TINKER_MOD(mplpot, pentyp)[7];

double& m2scale = TINKER_MOD(mplpot, m2scale);
double& m3scale = TINKER_MOD(mplpot, m3scale);
double& m4scale = TINKER_MOD(mplpot, m4scale);
double& m5scale = TINKER_MOD(mplpot, m5scale);
int& use_chgpen = TINKER_MOD(mplpot, use_chgpen);
char (&pentyp)[7] = TINKER_MOD(mplpot, pentyp);
#endif
} }
