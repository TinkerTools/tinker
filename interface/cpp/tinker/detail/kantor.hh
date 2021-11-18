#pragma once

#include "macro.hh"

namespace tinker { namespace kantor {
const int maxnat = 500;
extern double (&atcon)[maxnat][6];
extern char (&kat)[maxnat][16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(kantor, atcon)[maxnat][6];
extern "C" char TINKER_MOD(kantor, kat)[maxnat][16];

double (&atcon)[maxnat][6] = TINKER_MOD(kantor, atcon);
char (&kat)[maxnat][16] = TINKER_MOD(kantor, kat);
#endif
} }
