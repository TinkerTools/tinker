#pragma once

#include "macro.hh"

namespace tinker { namespace mutant {
extern int& nmut;
extern int& nmutb;
extern int& vcouple;
extern int*& imut;
extern int*& type0;
extern int*& class0;
extern int*& type1;
extern int*& class1;
extern int*& mutg;
extern double& lambda;
extern double& vlambda;
extern double& elambda;
extern double& plambda;
extern double& tlambda;
extern double& scexp;
extern double& scalpha;
extern int& mutfield;
extern int& setelambda;
extern int& setplambda;
extern int& setvlambda;
extern int& use_past;
extern int& use_rel;
extern int& use_subsys;
extern int*& mut;
extern int*& subon;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(mutant, nmut);
extern "C" int TINKER_MOD(mutant, nmutb);
extern "C" int TINKER_MOD(mutant, vcouple);
extern "C" int* TINKER_MOD(mutant, imut);
extern "C" int* TINKER_MOD(mutant, type0);
extern "C" int* TINKER_MOD(mutant, class0);
extern "C" int* TINKER_MOD(mutant, type1);
extern "C" int* TINKER_MOD(mutant, class1);
extern "C" int* TINKER_MOD(mutant, mutg);
extern "C" double TINKER_MOD(mutant, lambda);
extern "C" double TINKER_MOD(mutant, vlambda);
extern "C" double TINKER_MOD(mutant, elambda);
extern "C" double TINKER_MOD(mutant, plambda);
extern "C" double TINKER_MOD(mutant, tlambda);
extern "C" double TINKER_MOD(mutant, scexp);
extern "C" double TINKER_MOD(mutant, scalpha);
extern "C" int TINKER_MOD(mutant, mutfield);
extern "C" int TINKER_MOD(mutant, setelambda);
extern "C" int TINKER_MOD(mutant, setplambda);
extern "C" int TINKER_MOD(mutant, setvlambda);
extern "C" int TINKER_MOD(mutant, use_past);
extern "C" int TINKER_MOD(mutant, use_rel);
extern "C" int TINKER_MOD(mutant, use_subsys);
extern "C" int* TINKER_MOD(mutant, mut);
extern "C" int* TINKER_MOD(mutant, subon);

int& nmut = TINKER_MOD(mutant, nmut);
int& nmutb = TINKER_MOD(mutant, nmutb);
int& vcouple = TINKER_MOD(mutant, vcouple);
int*& imut = TINKER_MOD(mutant, imut);
int*& type0 = TINKER_MOD(mutant, type0);
int*& class0 = TINKER_MOD(mutant, class0);
int*& type1 = TINKER_MOD(mutant, type1);
int*& class1 = TINKER_MOD(mutant, class1);
int*& mutg = TINKER_MOD(mutant, mutg);
double& lambda = TINKER_MOD(mutant, lambda);
double& vlambda = TINKER_MOD(mutant, vlambda);
double& elambda = TINKER_MOD(mutant, elambda);
double& plambda = TINKER_MOD(mutant, plambda);
double& tlambda = TINKER_MOD(mutant, tlambda);
double& scexp = TINKER_MOD(mutant, scexp);
double& scalpha = TINKER_MOD(mutant, scalpha);
int& mutfield = TINKER_MOD(mutant, mutfield);
int& setelambda = TINKER_MOD(mutant, setelambda);
int& setplambda = TINKER_MOD(mutant, setplambda);
int& setvlambda = TINKER_MOD(mutant, setvlambda);
int& use_past = TINKER_MOD(mutant, use_past);
int& use_rel = TINKER_MOD(mutant, use_rel);
int& use_subsys = TINKER_MOD(mutant, use_subsys);
int*& mut = TINKER_MOD(mutant, mut);
int*& subon = TINKER_MOD(mutant, subon);
#endif
} }
