#pragma once

#include "macro.hh"

namespace tinker { namespace kanang {
extern double*& anan;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(kanang, anan);

double*& anan = TINKER_MOD(kanang, anan);
#endif
} }
