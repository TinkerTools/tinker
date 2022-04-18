#pragma once

#include "macro.hh"

namespace tinker { namespace kantor {
extern int& maxnat;
extern double*& atcon;
extern char (*&kat)[16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(kantor, maxnat);
extern "C" double* TINKER_MOD(kantor, atcon);
extern "C" char (*TINKER_MOD(kantor, kat))[16];

int& maxnat = TINKER_MOD(kantor, maxnat);
double*& atcon = TINKER_MOD(kantor, atcon);
char (*&kat)[16] = TINKER_MOD(kantor, kat);
#endif
} }
