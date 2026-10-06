#pragma once

#include "macro.hh"

namespace tinker { namespace thrmint {
extern int& tinbcount;
extern int& tinbsave;
extern int& tinbtot;
extern int& tinstepavg;
extern double*& tidedllist;
extern double*& tilmdadedl;
extern double*& tilmdadedlstd;
extern double*& tilmdahist;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(thrmint, tinbcount);
extern "C" int TINKER_MOD(thrmint, tinbsave);
extern "C" int TINKER_MOD(thrmint, tinbtot);
extern "C" int TINKER_MOD(thrmint, tinstepavg);
extern "C" double* TINKER_MOD(thrmint, tidedllist);
extern "C" double* TINKER_MOD(thrmint, tilmdadedl);
extern "C" double* TINKER_MOD(thrmint, tilmdadedlstd);
extern "C" double* TINKER_MOD(thrmint, tilmdahist);

int& tinbcount = TINKER_MOD(thrmint, tinbcount);
int& tinbsave = TINKER_MOD(thrmint, tinbsave);
int& tinbtot = TINKER_MOD(thrmint, tinbtot);
int& tinstepavg = TINKER_MOD(thrmint, tinstepavg);
double*& tidedllist = TINKER_MOD(thrmint, tidedllist);
double*& tilmdadedl = TINKER_MOD(thrmint, tilmdadedl);
double*& tilmdadedlstd = TINKER_MOD(thrmint, tilmdadedlstd);
double*& tilmdahist = TINKER_MOD(thrmint, tilmdahist);
#endif
} }
