#pragma once

#include "macro.hh"

namespace tinker { namespace params {
extern int& nprm;
extern char (*&prmline)[240];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(params, nprm);
extern "C" char (*TINKER_MOD(params, prmline))[240];

int& nprm = TINKER_MOD(params, nprm);
char (*&prmline)[240] = TINKER_MOD(params, prmline);
#endif
} }
