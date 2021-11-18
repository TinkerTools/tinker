#pragma once

#include "macro.hh"

namespace tinker { namespace krepl {
extern double*& prsiz;
extern double*& prdmp;
extern double*& prele;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(krepl, prsiz);
extern "C" double* TINKER_MOD(krepl, prdmp);
extern "C" double* TINKER_MOD(krepl, prele);

double*& prsiz = TINKER_MOD(krepl, prsiz);
double*& prdmp = TINKER_MOD(krepl, prdmp);
double*& prele = TINKER_MOD(krepl, prele);
#endif
} }
