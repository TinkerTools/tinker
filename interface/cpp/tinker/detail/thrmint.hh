#pragma once

#include "macro.hh"

namespace tinker { namespace thrmint {
extern int& tibin;
extern int& tinbcount;
extern int& tinbin;
extern int& tinblock;
extern int& tinbsave;
extern int& tinbtot;
extern int& tinequil;
extern int& tinstepavg;
extern int& tiwindow;
extern int*& tiwinend;
extern double& tieqratio;
extern double& tilmda;
extern double*& tidedllist;
extern double*& tifraclist;
extern double*& tilmdadedl;
extern double*& tilmdadedlstd;
extern double*& tilmdahist;
extern double*& tilmdalist;
extern char (&tifile)[240];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(thrmint, tibin);
extern "C" int TINKER_MOD(thrmint, tinbcount);
extern "C" int TINKER_MOD(thrmint, tinbin);
extern "C" int TINKER_MOD(thrmint, tinblock);
extern "C" int TINKER_MOD(thrmint, tinbsave);
extern "C" int TINKER_MOD(thrmint, tinbtot);
extern "C" int TINKER_MOD(thrmint, tinequil);
extern "C" int TINKER_MOD(thrmint, tinstepavg);
extern "C" int TINKER_MOD(thrmint, tiwindow);
extern "C" int* TINKER_MOD(thrmint, tiwinend);
extern "C" double TINKER_MOD(thrmint, tieqratio);
extern "C" double TINKER_MOD(thrmint, tilmda);
extern "C" double* TINKER_MOD(thrmint, tidedllist);
extern "C" double* TINKER_MOD(thrmint, tifraclist);
extern "C" double* TINKER_MOD(thrmint, tilmdadedl);
extern "C" double* TINKER_MOD(thrmint, tilmdadedlstd);
extern "C" double* TINKER_MOD(thrmint, tilmdahist);
extern "C" double* TINKER_MOD(thrmint, tilmdalist);
extern "C" char TINKER_MOD(thrmint, tifile)[240];

int& tibin = TINKER_MOD(thrmint, tibin);
int& tinbcount = TINKER_MOD(thrmint, tinbcount);
int& tinbin = TINKER_MOD(thrmint, tinbin);
int& tinblock = TINKER_MOD(thrmint, tinblock);
int& tinbsave = TINKER_MOD(thrmint, tinbsave);
int& tinbtot = TINKER_MOD(thrmint, tinbtot);
int& tinequil = TINKER_MOD(thrmint, tinequil);
int& tinstepavg = TINKER_MOD(thrmint, tinstepavg);
int& tiwindow = TINKER_MOD(thrmint, tiwindow);
int*& tiwinend = TINKER_MOD(thrmint, tiwinend);
double& tieqratio = TINKER_MOD(thrmint, tieqratio);
double& tilmda = TINKER_MOD(thrmint, tilmda);
double*& tidedllist = TINKER_MOD(thrmint, tidedllist);
double*& tifraclist = TINKER_MOD(thrmint, tifraclist);
double*& tilmdadedl = TINKER_MOD(thrmint, tilmdadedl);
double*& tilmdadedlstd = TINKER_MOD(thrmint, tilmdadedlstd);
double*& tilmdahist = TINKER_MOD(thrmint, tilmdahist);
double*& tilmdalist = TINKER_MOD(thrmint, tilmdalist);
char (&tifile)[240] = TINKER_MOD(thrmint, tifile);
#endif
} }
