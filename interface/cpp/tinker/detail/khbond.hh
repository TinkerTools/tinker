#pragma once

#include "macro.hh"

namespace tinker { namespace khbond {
const int maxnhb = 500;
extern double (&radhb)[maxnhb];
extern double (&epshb)[maxnhb];
extern char (&khb)[maxnhb][8];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(khbond, radhb)[maxnhb];
extern "C" double TINKER_MOD(khbond, epshb)[maxnhb];
extern "C" char TINKER_MOD(khbond, khb)[maxnhb][8];

double (&radhb)[maxnhb] = TINKER_MOD(khbond, radhb);
double (&epshb)[maxnhb] = TINKER_MOD(khbond, epshb);
char (&khb)[maxnhb][8] = TINKER_MOD(khbond, khb);
#endif
} }
