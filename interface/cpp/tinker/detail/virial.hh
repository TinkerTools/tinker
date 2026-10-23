#pragma once

#include "macro.hh"

namespace tinker { namespace virial {
extern double (&vir)[3][3];
extern double (&epvir)[3][3];
extern double (&emvir)[3][3];
extern double (&evvir)[3][3];
extern double (&ectvir)[3][3];
extern double (&ervir)[3][3];
extern double (&edspvir)[3][3];
extern int& use_virial;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(virial, vir)[3][3];
extern "C" double TINKER_MOD(virial, epvir)[3][3];
extern "C" double TINKER_MOD(virial, emvir)[3][3];
extern "C" double TINKER_MOD(virial, evvir)[3][3];
extern "C" double TINKER_MOD(virial, ectvir)[3][3];
extern "C" double TINKER_MOD(virial, ervir)[3][3];
extern "C" double TINKER_MOD(virial, edspvir)[3][3];
extern "C" int TINKER_MOD(virial, use_virial);

double (&vir)[3][3] = TINKER_MOD(virial, vir);
double (&epvir)[3][3] = TINKER_MOD(virial, epvir);
double (&emvir)[3][3] = TINKER_MOD(virial, emvir);
double (&evvir)[3][3] = TINKER_MOD(virial, evvir);
double (&ectvir)[3][3] = TINKER_MOD(virial, ectvir);
double (&ervir)[3][3] = TINKER_MOD(virial, ervir);
double (&edspvir)[3][3] = TINKER_MOD(virial, edspvir);
int& use_virial = TINKER_MOD(virial, use_virial);
#endif
} }
