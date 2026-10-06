#pragma once

#include "macro.hh"

namespace tinker { namespace fep {
extern int& fepintv;
extern int& nfep;
extern int& nfepsave;
extern int& nfeptot;
extern int*& fepstep;
extern double*& fepene;
extern double*& feplmda;
extern double*& feptrial;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(fep, fepintv);
extern "C" int TINKER_MOD(fep, nfep);
extern "C" int TINKER_MOD(fep, nfepsave);
extern "C" int TINKER_MOD(fep, nfeptot);
extern "C" int* TINKER_MOD(fep, fepstep);
extern "C" double* TINKER_MOD(fep, fepene);
extern "C" double* TINKER_MOD(fep, feplmda);
extern "C" double* TINKER_MOD(fep, feptrial);

int& fepintv = TINKER_MOD(fep, fepintv);
int& nfep = TINKER_MOD(fep, nfep);
int& nfepsave = TINKER_MOD(fep, nfepsave);
int& nfeptot = TINKER_MOD(fep, nfeptot);
int*& fepstep = TINKER_MOD(fep, fepstep);
double*& fepene = TINKER_MOD(fep, fepene);
double*& feplmda = TINKER_MOD(fep, feplmda);
double*& feptrial = TINKER_MOD(fep, feptrial);
#endif
} }
