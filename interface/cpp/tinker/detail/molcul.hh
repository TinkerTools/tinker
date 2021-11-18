#pragma once

#include "macro.hh"

namespace tinker { namespace molcul {
extern int& nmol;
extern int*& imol;
extern int*& kmol;
extern int*& molcule;
extern double& totmass;
extern double*& molmass;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(molcul, nmol);
extern "C" int* TINKER_MOD(molcul, imol);
extern "C" int* TINKER_MOD(molcul, kmol);
extern "C" int* TINKER_MOD(molcul, molcule);
extern "C" double TINKER_MOD(molcul, totmass);
extern "C" double* TINKER_MOD(molcul, molmass);

int& nmol = TINKER_MOD(molcul, nmol);
int*& imol = TINKER_MOD(molcul, imol);
int*& kmol = TINKER_MOD(molcul, kmol);
int*& molcule = TINKER_MOD(molcul, molcule);
double& totmass = TINKER_MOD(molcul, totmass);
double*& molmass = TINKER_MOD(molcul, molmass);
#endif
} }
