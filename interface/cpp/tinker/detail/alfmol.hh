#pragma once

#include "macro.hh"

namespace tinker { namespace alfmol {
extern int& alfthread;
extern double& delcxeps;
extern int& alfhydro;
extern int& alfsosgmp;
extern char (&alfmethod)[6];
extern char (&alfsort)[6];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(alfmol, alfthread);
extern "C" double TINKER_MOD(alfmol, delcxeps);
extern "C" int TINKER_MOD(alfmol, alfhydro);
extern "C" int TINKER_MOD(alfmol, alfsosgmp);
extern "C" char TINKER_MOD(alfmol, alfmethod)[6];
extern "C" char TINKER_MOD(alfmol, alfsort)[6];

int& alfthread = TINKER_MOD(alfmol, alfthread);
double& delcxeps = TINKER_MOD(alfmol, delcxeps);
int& alfhydro = TINKER_MOD(alfmol, alfhydro);
int& alfsosgmp = TINKER_MOD(alfmol, alfsosgmp);
char (&alfmethod)[6] = TINKER_MOD(alfmol, alfmethod);
char (&alfsort)[6] = TINKER_MOD(alfmol, alfsort);
#endif
} }
