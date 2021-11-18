#pragma once

#include "macro.hh"

namespace tinker { namespace ielscf {
extern int& nfree_aux;
extern double& tautemp_aux;
extern double& kelvin_aux;
extern double*& uaux;
extern double*& upaux;
extern double*& vaux;
extern double*& vpaux;
extern double*& aaux;
extern double*& apaux;
extern int& use_ielscf;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(ielscf, nfree_aux);
extern "C" double TINKER_MOD(ielscf, tautemp_aux);
extern "C" double TINKER_MOD(ielscf, kelvin_aux);
extern "C" double* TINKER_MOD(ielscf, uaux);
extern "C" double* TINKER_MOD(ielscf, upaux);
extern "C" double* TINKER_MOD(ielscf, vaux);
extern "C" double* TINKER_MOD(ielscf, vpaux);
extern "C" double* TINKER_MOD(ielscf, aaux);
extern "C" double* TINKER_MOD(ielscf, apaux);
extern "C" int TINKER_MOD(ielscf, use_ielscf);

int& nfree_aux = TINKER_MOD(ielscf, nfree_aux);
double& tautemp_aux = TINKER_MOD(ielscf, tautemp_aux);
double& kelvin_aux = TINKER_MOD(ielscf, kelvin_aux);
double*& uaux = TINKER_MOD(ielscf, uaux);
double*& upaux = TINKER_MOD(ielscf, upaux);
double*& vaux = TINKER_MOD(ielscf, vaux);
double*& vpaux = TINKER_MOD(ielscf, vpaux);
double*& aaux = TINKER_MOD(ielscf, aaux);
double*& apaux = TINKER_MOD(ielscf, apaux);
int& use_ielscf = TINKER_MOD(ielscf, use_ielscf);
#endif
} }
