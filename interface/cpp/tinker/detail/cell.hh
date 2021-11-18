#pragma once

#include "macro.hh"

namespace tinker { namespace cell {
extern int& ncell;
extern int*& icell;
extern double& xcell;
extern double& ycell;
extern double& zcell;
extern double& xcell2;
extern double& ycell2;
extern double& zcell2;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(cell, ncell);
extern "C" int* TINKER_MOD(cell, icell);
extern "C" double TINKER_MOD(cell, xcell);
extern "C" double TINKER_MOD(cell, ycell);
extern "C" double TINKER_MOD(cell, zcell);
extern "C" double TINKER_MOD(cell, xcell2);
extern "C" double TINKER_MOD(cell, ycell2);
extern "C" double TINKER_MOD(cell, zcell2);

int& ncell = TINKER_MOD(cell, ncell);
int*& icell = TINKER_MOD(cell, icell);
double& xcell = TINKER_MOD(cell, xcell);
double& ycell = TINKER_MOD(cell, ycell);
double& zcell = TINKER_MOD(cell, zcell);
double& xcell2 = TINKER_MOD(cell, xcell2);
double& ycell2 = TINKER_MOD(cell, ycell2);
double& zcell2 = TINKER_MOD(cell, zcell2);
#endif
} }
