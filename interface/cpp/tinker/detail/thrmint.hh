#pragma once

#include "macro.hh"

namespace tinker { namespace thrmint {
extern int& tibin;
extern int& tinbin;
extern int& tinequil;
extern int& tinstepavg;
extern int& tiwindow;
extern double& tieqratio;
extern double& tilmda;
extern double*& tidedllist;
extern double*& tilmdadedl;
extern double*& tilmdadedlstd;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(thrmint, tibin);
extern "C" int TINKER_MOD(thrmint, tinbin);
extern "C" int TINKER_MOD(thrmint, tinequil);
extern "C" int TINKER_MOD(thrmint, tinstepavg);
extern "C" int TINKER_MOD(thrmint, tiwindow);
extern "C" double TINKER_MOD(thrmint, tieqratio);
extern "C" double TINKER_MOD(thrmint, tilmda);
extern "C" double* TINKER_MOD(thrmint, tidedllist);
extern "C" double* TINKER_MOD(thrmint, tilmdadedl);
extern "C" double* TINKER_MOD(thrmint, tilmdadedlstd);

int& tibin = TINKER_MOD(thrmint, tibin);
int& tinbin = TINKER_MOD(thrmint, tinbin);
int& tinequil = TINKER_MOD(thrmint, tinequil);
int& tinstepavg = TINKER_MOD(thrmint, tinstepavg);
int& tiwindow = TINKER_MOD(thrmint, tiwindow);
double& tieqratio = TINKER_MOD(thrmint, tieqratio);
double& tilmda = TINKER_MOD(thrmint, tilmda);
double*& tidedllist = TINKER_MOD(thrmint, tidedllist);
double*& tilmdadedl = TINKER_MOD(thrmint, tilmdadedl);
double*& tilmdadedlstd = TINKER_MOD(thrmint, tilmdadedlstd);
#endif
} }
