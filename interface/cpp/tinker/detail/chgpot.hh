#pragma once

#include "macro.hh"

namespace tinker { namespace chgpot {
extern double& electric;
extern double& dielec;
extern double& ebuffer;
extern double& c1scale;
extern double& c2scale;
extern double& c3scale;
extern double& c4scale;
extern double& c5scale;
extern int& neutnbr;
extern int& neutcut;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(chgpot, electric);
extern "C" double TINKER_MOD(chgpot, dielec);
extern "C" double TINKER_MOD(chgpot, ebuffer);
extern "C" double TINKER_MOD(chgpot, c1scale);
extern "C" double TINKER_MOD(chgpot, c2scale);
extern "C" double TINKER_MOD(chgpot, c3scale);
extern "C" double TINKER_MOD(chgpot, c4scale);
extern "C" double TINKER_MOD(chgpot, c5scale);
extern "C" int TINKER_MOD(chgpot, neutnbr);
extern "C" int TINKER_MOD(chgpot, neutcut);

double& electric = TINKER_MOD(chgpot, electric);
double& dielec = TINKER_MOD(chgpot, dielec);
double& ebuffer = TINKER_MOD(chgpot, ebuffer);
double& c1scale = TINKER_MOD(chgpot, c1scale);
double& c2scale = TINKER_MOD(chgpot, c2scale);
double& c3scale = TINKER_MOD(chgpot, c3scale);
double& c4scale = TINKER_MOD(chgpot, c4scale);
double& c5scale = TINKER_MOD(chgpot, c5scale);
int& neutnbr = TINKER_MOD(chgpot, neutnbr);
int& neutcut = TINKER_MOD(chgpot, neutcut);
#endif
} }
