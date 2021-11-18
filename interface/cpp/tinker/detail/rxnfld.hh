#pragma once

#include "macro.hh"

namespace tinker { namespace rxnfld {
extern int (&ijk)[6][6][6];
extern double (&b1)[13][40];
extern double (&b2)[13][40];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(rxnfld, ijk)[6][6][6];
extern "C" double TINKER_MOD(rxnfld, b1)[13][40];
extern "C" double TINKER_MOD(rxnfld, b2)[13][40];

int (&ijk)[6][6][6] = TINKER_MOD(rxnfld, ijk);
double (&b1)[13][40] = TINKER_MOD(rxnfld, b1);
double (&b2)[13][40] = TINKER_MOD(rxnfld, b2);
#endif
} }
