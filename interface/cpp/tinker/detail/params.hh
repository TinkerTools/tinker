#pragma once

#include "macro.hh"

namespace tinker { namespace params {
const int maxprm = 25000;
extern int& nprm;
extern char (&prmline)[maxprm][240];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(params, nprm);
extern "C" char TINKER_MOD(params, prmline)[maxprm][240];

int& nprm = TINKER_MOD(params, nprm);
char (&prmline)[maxprm][240] = TINKER_MOD(params, prmline);
#endif
} }
