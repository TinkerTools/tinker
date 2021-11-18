#pragma once

#include "macro.hh"

namespace tinker { namespace usage {
extern int& nuse;
extern int*& iuse;
extern int*& use;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(usage, nuse);
extern "C" int* TINKER_MOD(usage, iuse);
extern "C" int* TINKER_MOD(usage, use);

int& nuse = TINKER_MOD(usage, nuse);
int*& iuse = TINKER_MOD(usage, iuse);
int*& use = TINKER_MOD(usage, use);
#endif
} }
