#pragma once

#include "macro.hh"

namespace tinker { namespace openmp {
extern int& nproc;
extern int& nthread;
extern int& nnest;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(openmp, nproc);
extern "C" int TINKER_MOD(openmp, nthread);
extern "C" int TINKER_MOD(openmp, nnest);

int& nproc = TINKER_MOD(openmp, nproc);
int& nthread = TINKER_MOD(openmp, nthread);
int& nnest = TINKER_MOD(openmp, nnest);
#endif
} }
