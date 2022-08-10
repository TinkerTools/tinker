#pragma once

#include "macro.hh"

namespace tinker { namespace kpolr {
extern int*& pgrp;
extern double*& polr;
extern double*& athl;
extern double*& dthl;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int* TINKER_MOD(kpolr, pgrp);
extern "C" double* TINKER_MOD(kpolr, polr);
extern "C" double* TINKER_MOD(kpolr, athl);
extern "C" double* TINKER_MOD(kpolr, dthl);

int*& pgrp = TINKER_MOD(kpolr, pgrp);
double*& polr = TINKER_MOD(kpolr, polr);
double*& athl = TINKER_MOD(kpolr, athl);
double*& dthl = TINKER_MOD(kpolr, dthl);
#endif
} }
