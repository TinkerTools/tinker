#pragma once

#include "macro.hh"

namespace tinker { namespace poltcg {
extern int& tcgorder;
extern int& tcgnab;
extern double& tcgpeek;
extern double*& uad;
extern double*& uap;
extern double*& ubd;
extern double*& ubp;
extern int& tcgguess;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(poltcg, tcgorder);
extern "C" int TINKER_MOD(poltcg, tcgnab);
extern "C" double TINKER_MOD(poltcg, tcgpeek);
extern "C" double* TINKER_MOD(poltcg, uad);
extern "C" double* TINKER_MOD(poltcg, uap);
extern "C" double* TINKER_MOD(poltcg, ubd);
extern "C" double* TINKER_MOD(poltcg, ubp);
extern "C" int TINKER_MOD(poltcg, tcgguess);

int& tcgorder = TINKER_MOD(poltcg, tcgorder);
int& tcgnab = TINKER_MOD(poltcg, tcgnab);
double& tcgpeek = TINKER_MOD(poltcg, tcgpeek);
double*& uad = TINKER_MOD(poltcg, uad);
double*& uap = TINKER_MOD(poltcg, uap);
double*& ubd = TINKER_MOD(poltcg, ubd);
double*& ubp = TINKER_MOD(poltcg, ubp);
int& tcgguess = TINKER_MOD(poltcg, tcgguess);
#endif
} }
