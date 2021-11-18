#pragma once

#include "macro.hh"
#include "sizes.hh"

namespace tinker { namespace atomid {
using namespace sizes;

extern int (&tag)[maxatm];
extern int (&class_)[maxatm];
extern int (&atomic)[maxatm];
extern int (&valence)[maxatm];
extern double (&mass)[maxatm];
extern char (&name)[maxatm][3];
extern char (&story)[maxatm][24];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(atomid, tag)[maxatm];
extern "C" int TINKER_MOD(atomid, class)[maxatm];
extern "C" int TINKER_MOD(atomid, atomic)[maxatm];
extern "C" int TINKER_MOD(atomid, valence)[maxatm];
extern "C" double TINKER_MOD(atomid, mass)[maxatm];
extern "C" char TINKER_MOD(atomid, name)[maxatm][3];
extern "C" char TINKER_MOD(atomid, story)[maxatm][24];

int (&tag)[maxatm] = TINKER_MOD(atomid, tag);
int (&class_)[maxatm] = TINKER_MOD(atomid, class);
int (&atomic)[maxatm] = TINKER_MOD(atomid, atomic);
int (&valence)[maxatm] = TINKER_MOD(atomid, valence);
double (&mass)[maxatm] = TINKER_MOD(atomid, mass);
char (&name)[maxatm][3] = TINKER_MOD(atomid, name);
char (&story)[maxatm][24] = TINKER_MOD(atomid, story);
#endif
} }
