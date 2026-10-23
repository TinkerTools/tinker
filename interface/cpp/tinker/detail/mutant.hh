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
extern double& scalphar;
extern double& scalphav;
extern double& scexp;
extern double& scexpm;
extern double& scexpp;
extern double& scexpr;
extern double& scexprm;
extern double& scexprp;
extern double& scrdsp;
extern double& scexpd;
extern double& scexpdq;
extern double& tlambda;
extern double& vlambda;
extern double*& fluxsc;
extern double*& pdltflx;
extern int& mutfield;
extern int& mutpsc;
extern int& setelambda;
extern int& setplambda;
extern int& setvlambda;
extern int& use_rel;

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
extern "C" double TINKER_MOD(mutant, scalphar);
extern "C" double TINKER_MOD(mutant, scalphav);
extern "C" double TINKER_MOD(mutant, scexp);
extern "C" double TINKER_MOD(mutant, scexpm);
extern "C" double TINKER_MOD(mutant, scexpp);
extern "C" double TINKER_MOD(mutant, scexpr);
extern "C" double TINKER_MOD(mutant, scexprm);
extern "C" double TINKER_MOD(mutant, scexprp);
extern "C" double TINKER_MOD(mutant, scrdsp);
extern "C" double TINKER_MOD(mutant, scexpd);
extern "C" double TINKER_MOD(mutant, scexpdq);
extern "C" double TINKER_MOD(mutant, tlambda);
extern "C" double TINKER_MOD(mutant, vlambda);
extern "C" double* TINKER_MOD(mutant, fluxsc);
extern "C" double* TINKER_MOD(mutant, pdltflx);
extern "C" int TINKER_MOD(mutant, mutfield);
extern "C" int TINKER_MOD(mutant, mutpsc);
extern "C" int TINKER_MOD(mutant, setelambda);
extern "C" int TINKER_MOD(mutant, setplambda);
extern "C" int TINKER_MOD(mutant, setvlambda);
extern "C" int TINKER_MOD(mutant, use_rel);

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
double& scalphar = TINKER_MOD(mutant, scalphar);
double& scalphav = TINKER_MOD(mutant, scalphav);
double& scexp = TINKER_MOD(mutant, scexp);
double& scexpm = TINKER_MOD(mutant, scexpm);
double& scexpp = TINKER_MOD(mutant, scexpp);
double& scexpr = TINKER_MOD(mutant, scexpr);
double& scexprm = TINKER_MOD(mutant, scexprm);
double& scexprp = TINKER_MOD(mutant, scexprp);
double& scrdsp = TINKER_MOD(mutant, scrdsp);
double& scexpd = TINKER_MOD(mutant, scexpd);
double& scexpdq = TINKER_MOD(mutant, scexpdq);
double& tlambda = TINKER_MOD(mutant, tlambda);
double& vlambda = TINKER_MOD(mutant, vlambda);
double*& fluxsc = TINKER_MOD(mutant, fluxsc);
double*& pdltflx = TINKER_MOD(mutant, pdltflx);
int& mutfield = TINKER_MOD(mutant, mutfield);
int& mutpsc = TINKER_MOD(mutant, mutpsc);
int& setelambda = TINKER_MOD(mutant, setelambda);
int& setplambda = TINKER_MOD(mutant, setplambda);
int& setvlambda = TINKER_MOD(mutant, setvlambda);
int& use_rel = TINKER_MOD(mutant, use_rel);
#endif
} }
