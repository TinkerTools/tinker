#pragma once

#include "macro.hh"

namespace tinker { namespace tortor {
extern int& ntortor;
extern int*& itt;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(tortor, ntortor);
extern "C" int* TINKER_MOD(tortor, itt);

int& ntortor = TINKER_MOD(tortor, ntortor);
int*& itt = TINKER_MOD(tortor, itt);
#endif
} }
