#pragma once

#include "macro.hh"

namespace tinker { namespace mutant {
extern int& nmut;
extern int& nmutb;
extern int& vcouple;
extern int*& class0;
extern int*& class1;
extern int*& imut;
extern int*& mutg;
extern int*& type0;
extern int*& type1;
extern double& elambda;
extern double& lambda;
extern double& plambda;
extern double& scalpham;
extern double& scalphav;
extern double& scexp;
extern double& tlambda;
extern double& vlambda;
extern int& mutfield;
extern int& setelambda;
extern int& setesoft;
extern int& setplambda;
extern int& setvlambda;
extern int& use_esoft;
extern int& use_past;
extern int& use_rel;
extern int& use_subsys;
extern int*& mut;
extern int*& subon;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(mutant, nmut);
extern "C" int TINKER_MOD(mutant, nmutb);
extern "C" int TINKER_MOD(mutant, vcouple);
extern "C" int* TINKER_MOD(mutant, class0);
extern "C" int* TINKER_MOD(mutant, class1);
extern "C" int* TINKER_MOD(mutant, imut);
extern "C" int* TINKER_MOD(mutant, mutg);
extern "C" int* TINKER_MOD(mutant, type0);
extern "C" int* TINKER_MOD(mutant, type1);
extern "C" double TINKER_MOD(mutant, elambda);
extern "C" double TINKER_MOD(mutant, lambda);
extern "C" double TINKER_MOD(mutant, plambda);
extern "C" double TINKER_MOD(mutant, scalpham);
extern "C" double TINKER_MOD(mutant, scalphav);
extern "C" double TINKER_MOD(mutant, scexp);
extern "C" double TINKER_MOD(mutant, tlambda);
extern "C" double TINKER_MOD(mutant, vlambda);
extern "C" int TINKER_MOD(mutant, mutfield);
extern "C" int TINKER_MOD(mutant, setelambda);
extern "C" int TINKER_MOD(mutant, setesoft);
extern "C" int TINKER_MOD(mutant, setplambda);
extern "C" int TINKER_MOD(mutant, setvlambda);
extern "C" int TINKER_MOD(mutant, use_esoft);
extern "C" int TINKER_MOD(mutant, use_past);
extern "C" int TINKER_MOD(mutant, use_rel);
extern "C" int TINKER_MOD(mutant, use_subsys);
extern "C" int* TINKER_MOD(mutant, mut);
extern "C" int* TINKER_MOD(mutant, subon);

int& nmut = TINKER_MOD(mutant, nmut);
int& nmutb = TINKER_MOD(mutant, nmutb);
int& vcouple = TINKER_MOD(mutant, vcouple);
int*& class0 = TINKER_MOD(mutant, class0);
int*& class1 = TINKER_MOD(mutant, class1);
int*& imut = TINKER_MOD(mutant, imut);
int*& mutg = TINKER_MOD(mutant, mutg);
int*& type0 = TINKER_MOD(mutant, type0);
int*& type1 = TINKER_MOD(mutant, type1);
double& elambda = TINKER_MOD(mutant, elambda);
double& lambda = TINKER_MOD(mutant, lambda);
double& plambda = TINKER_MOD(mutant, plambda);
double& scalpham = TINKER_MOD(mutant, scalpham);
double& scalphav = TINKER_MOD(mutant, scalphav);
double& scexp = TINKER_MOD(mutant, scexp);
double& tlambda = TINKER_MOD(mutant, tlambda);
double& vlambda = TINKER_MOD(mutant, vlambda);
int& mutfield = TINKER_MOD(mutant, mutfield);
int& setelambda = TINKER_MOD(mutant, setelambda);
int& setesoft = TINKER_MOD(mutant, setesoft);
int& setplambda = TINKER_MOD(mutant, setplambda);
int& setvlambda = TINKER_MOD(mutant, setvlambda);
int& use_esoft = TINKER_MOD(mutant, use_esoft);
int& use_past = TINKER_MOD(mutant, use_past);
int& use_rel = TINKER_MOD(mutant, use_rel);
int& use_subsys = TINKER_MOD(mutant, use_subsys);
int*& mut = TINKER_MOD(mutant, mut);
int*& subon = TINKER_MOD(mutant, subon);
#endif
} }
