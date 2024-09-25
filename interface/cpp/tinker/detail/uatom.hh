#pragma once

#include "macro.hh"
#include "sizes.hh"

namespace tinker { namespace uatom {
using namespace sizes;

extern int& nunique;
extern int (&utype)[maxtyp];
extern int (&utypeinv)[maxtyp];
extern double (&utv1)[maxtyp][3];
extern double (&utv2)[maxtyp][3];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(uatom, nunique);
extern "C" int TINKER_MOD(uatom, utype)[maxtyp];
extern "C" int TINKER_MOD(uatom, utypeinv)[maxtyp];
extern "C" double TINKER_MOD(uatom, utv1)[maxtyp][3];
extern "C" double TINKER_MOD(uatom, utv2)[maxtyp][3];

int& nunique = TINKER_MOD(uatom, nunique);
int (&utype)[maxtyp] = TINKER_MOD(uatom, utype);
int (&utypeinv)[maxtyp] = TINKER_MOD(uatom, utypeinv);
double (&utv1)[maxtyp][3] = TINKER_MOD(uatom, utv1);
double (&utv2)[maxtyp][3] = TINKER_MOD(uatom, utv2);
#endif
} }
