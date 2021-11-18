#pragma once

#include "macro.hh"

namespace tinker { namespace kcpen {
extern double*& cpele;
extern double*& cpalp;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(kcpen, cpele);
extern "C" double* TINKER_MOD(kcpen, cpalp);

double*& cpele = TINKER_MOD(kcpen, cpele);
double*& cpalp = TINKER_MOD(kcpen, cpalp);
#endif
} }
