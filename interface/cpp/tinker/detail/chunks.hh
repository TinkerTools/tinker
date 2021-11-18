#pragma once

#include "macro.hh"

namespace tinker { namespace chunks {
extern int& nchunk;
extern int& nchk1;
extern int& nchk2;
extern int& nchk3;
extern int& ngrd1;
extern int& ngrd2;
extern int& ngrd3;
extern int& nlpts;
extern int& nrpts;
extern int& grdoff;
extern int*& pmetable;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(chunks, nchunk);
extern "C" int TINKER_MOD(chunks, nchk1);
extern "C" int TINKER_MOD(chunks, nchk2);
extern "C" int TINKER_MOD(chunks, nchk3);
extern "C" int TINKER_MOD(chunks, ngrd1);
extern "C" int TINKER_MOD(chunks, ngrd2);
extern "C" int TINKER_MOD(chunks, ngrd3);
extern "C" int TINKER_MOD(chunks, nlpts);
extern "C" int TINKER_MOD(chunks, nrpts);
extern "C" int TINKER_MOD(chunks, grdoff);
extern "C" int* TINKER_MOD(chunks, pmetable);

int& nchunk = TINKER_MOD(chunks, nchunk);
int& nchk1 = TINKER_MOD(chunks, nchk1);
int& nchk2 = TINKER_MOD(chunks, nchk2);
int& nchk3 = TINKER_MOD(chunks, nchk3);
int& ngrd1 = TINKER_MOD(chunks, ngrd1);
int& ngrd2 = TINKER_MOD(chunks, ngrd2);
int& ngrd3 = TINKER_MOD(chunks, ngrd3);
int& nlpts = TINKER_MOD(chunks, nlpts);
int& nrpts = TINKER_MOD(chunks, nrpts);
int& grdoff = TINKER_MOD(chunks, grdoff);
int*& pmetable = TINKER_MOD(chunks, pmetable);
#endif
} }
