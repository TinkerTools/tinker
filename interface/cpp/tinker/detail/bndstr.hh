#pragma once

#include "macro.hh"

namespace tinker { namespace bndstr {
extern int& nbond;
extern int*& ibnd;
extern double*& bk;
extern double*& bl;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(bndstr, nbond);
extern "C" int* TINKER_MOD(bndstr, ibnd);
extern "C" double* TINKER_MOD(bndstr, bk);
extern "C" double* TINKER_MOD(bndstr, bl);

int& nbond = TINKER_MOD(bndstr, nbond);
int*& ibnd = TINKER_MOD(bndstr, ibnd);
double*& bk = TINKER_MOD(bndstr, bk);
double*& bl = TINKER_MOD(bndstr, bl);
#endif
} }
