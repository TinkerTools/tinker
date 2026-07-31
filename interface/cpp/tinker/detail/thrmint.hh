#pragma once

#include "macro.hh"

namespace tinker { namespace thrmint {
extern int& tibin;
extern int& tinbin;
extern int& tinblock;
extern int& tinequil;
extern int& tinstepavg;
extern int& tiwindow;
extern int*& tinbcount;
extern double& tieqratio;
extern double& tilmda;
extern double*& tidedllist;
extern double*& tilmdalist;
extern double*& tilmdadedl;
extern double*& tilmdadedlstd;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(thrmint, tibin);
extern "C" int TINKER_MOD(thrmint, tinbin);
extern "C" int TINKER_MOD(thrmint, tinblock);
extern "C" int TINKER_MOD(thrmint, tinequil);
extern "C" int TINKER_MOD(thrmint, tinstepavg);
extern "C" int TINKER_MOD(thrmint, tiwindow);
extern "C" int* TINKER_MOD(thrmint, tinbcount);
extern "C" double TINKER_MOD(thrmint, tieqratio);
extern "C" double TINKER_MOD(thrmint, tilmda);
extern "C" double* TINKER_MOD(thrmint, tidedllist);
extern "C" double* TINKER_MOD(thrmint, tilmdalist);
extern "C" double* TINKER_MOD(thrmint, tilmdadedl);
extern "C" double* TINKER_MOD(thrmint, tilmdadedlstd);

int& tibin = TINKER_MOD(thrmint, tibin);
int& tinbin = TINKER_MOD(thrmint, tinbin);
int& tinblock = TINKER_MOD(thrmint, tinblock);
int& tinequil = TINKER_MOD(thrmint, tinequil);
int& tinstepavg = TINKER_MOD(thrmint, tinstepavg);
int& tiwindow = TINKER_MOD(thrmint, tiwindow);
int*& tinbcount = TINKER_MOD(thrmint, tinbcount);
double& tieqratio = TINKER_MOD(thrmint, tieqratio);
double& tilmda = TINKER_MOD(thrmint, tilmda);
double*& tidedllist = TINKER_MOD(thrmint, tidedllist);
double*& tilmdalist = TINKER_MOD(thrmint, tilmdalist);
double*& tilmdadedl = TINKER_MOD(thrmint, tilmdadedl);
double*& tilmdadedlstd = TINKER_MOD(thrmint, tilmdadedlstd);
#endif
} }
