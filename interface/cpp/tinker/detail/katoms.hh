#pragma once

#include "macro.hh"

namespace tinker { namespace katoms {
extern int*& atmcls;
extern int*& atmnum;
extern int*& ligand;
extern double*& weight;
extern char (*&symbol)[3];
extern char (*&describe)[24];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int* TINKER_MOD(katoms, atmcls);
extern "C" int* TINKER_MOD(katoms, atmnum);
extern "C" int* TINKER_MOD(katoms, ligand);
extern "C" double* TINKER_MOD(katoms, weight);
extern "C" char (*TINKER_MOD(katoms, symbol))[3];
extern "C" char (*TINKER_MOD(katoms, describe))[24];

int*& atmcls = TINKER_MOD(katoms, atmcls);
int*& atmnum = TINKER_MOD(katoms, atmnum);
int*& ligand = TINKER_MOD(katoms, ligand);
double*& weight = TINKER_MOD(katoms, weight);
char (*&symbol)[3] = TINKER_MOD(katoms, symbol);
char (*&describe)[24] = TINKER_MOD(katoms, describe);
#endif
} }
