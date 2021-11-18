#pragma once

#include "macro.hh"

namespace tinker { namespace mutant {
extern int& nmut;
extern int& vcouple;
extern int*& imut;
extern int*& type0;
extern int*& class0;
extern int*& type1;
extern int*& class1;
extern double& lambda;
extern double& tlambda;
extern double& vlambda;
extern double& elambda;
extern double& scexp;
extern double& scalpha;
extern int*& mut;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(mutant, nmut);
extern "C" int TINKER_MOD(mutant, vcouple);
extern "C" int* TINKER_MOD(mutant, imut);
extern "C" int* TINKER_MOD(mutant, type0);
extern "C" int* TINKER_MOD(mutant, class0);
extern "C" int* TINKER_MOD(mutant, type1);
extern "C" int* TINKER_MOD(mutant, class1);
extern "C" double TINKER_MOD(mutant, lambda);
extern "C" double TINKER_MOD(mutant, tlambda);
extern "C" double TINKER_MOD(mutant, vlambda);
extern "C" double TINKER_MOD(mutant, elambda);
extern "C" double TINKER_MOD(mutant, scexp);
extern "C" double TINKER_MOD(mutant, scalpha);
extern "C" int* TINKER_MOD(mutant, mut);

int& nmut = TINKER_MOD(mutant, nmut);
int& vcouple = TINKER_MOD(mutant, vcouple);
int*& imut = TINKER_MOD(mutant, imut);
int*& type0 = TINKER_MOD(mutant, type0);
int*& class0 = TINKER_MOD(mutant, class0);
int*& type1 = TINKER_MOD(mutant, type1);
int*& class1 = TINKER_MOD(mutant, class1);
double& lambda = TINKER_MOD(mutant, lambda);
double& tlambda = TINKER_MOD(mutant, tlambda);
double& vlambda = TINKER_MOD(mutant, vlambda);
double& elambda = TINKER_MOD(mutant, elambda);
double& scexp = TINKER_MOD(mutant, scexp);
double& scalpha = TINKER_MOD(mutant, scalpha);
int*& mut = TINKER_MOD(mutant, mut);
#endif
} }
