#pragma once

#include "macro.hh"

namespace tinker { namespace xtals {
const int maxlsq = 1000;
const int maxrsd = 1000;
extern int& nxtal;
extern int& nvary;
extern int (&ivary)[maxlsq];
extern int (&iresid)[maxrsd];
extern int (&vary)[maxlsq][2];
extern double& e0_lattice;
extern char (&vartyp)[maxlsq][16];
extern char (&rsdtyp)[maxrsd][16];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(xtals, nxtal);
extern "C" int TINKER_MOD(xtals, nvary);
extern "C" int TINKER_MOD(xtals, ivary)[maxlsq];
extern "C" int TINKER_MOD(xtals, iresid)[maxrsd];
extern "C" int TINKER_MOD(xtals, vary)[maxlsq][2];
extern "C" double TINKER_MOD(xtals, e0_lattice);
extern "C" char TINKER_MOD(xtals, vartyp)[maxlsq][16];
extern "C" char TINKER_MOD(xtals, rsdtyp)[maxrsd][16];

int& nxtal = TINKER_MOD(xtals, nxtal);
int& nvary = TINKER_MOD(xtals, nvary);
int (&ivary)[maxlsq] = TINKER_MOD(xtals, ivary);
int (&iresid)[maxrsd] = TINKER_MOD(xtals, iresid);
int (&vary)[maxlsq][2] = TINKER_MOD(xtals, vary);
double& e0_lattice = TINKER_MOD(xtals, e0_lattice);
char (&vartyp)[maxlsq][16] = TINKER_MOD(xtals, vartyp);
char (&rsdtyp)[maxrsd][16] = TINKER_MOD(xtals, rsdtyp);
#endif
} }
