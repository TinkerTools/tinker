#pragma once

#include "macro.hh"

namespace tinker { namespace group {
extern int& ngrp;
extern int*& kgrp;
extern int*& grplist;
extern int*& igrp;
extern double*& grpmass;
extern double*& wgrp;
extern int& use_group;
extern int& use_intra;
extern int& use_inter;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(group, ngrp);
extern "C" int* TINKER_MOD(group, kgrp);
extern "C" int* TINKER_MOD(group, grplist);
extern "C" int* TINKER_MOD(group, igrp);
extern "C" double* TINKER_MOD(group, grpmass);
extern "C" double* TINKER_MOD(group, wgrp);
extern "C" int TINKER_MOD(group, use_group);
extern "C" int TINKER_MOD(group, use_intra);
extern "C" int TINKER_MOD(group, use_inter);

int& ngrp = TINKER_MOD(group, ngrp);
int*& kgrp = TINKER_MOD(group, kgrp);
int*& grplist = TINKER_MOD(group, grplist);
int*& igrp = TINKER_MOD(group, igrp);
double*& grpmass = TINKER_MOD(group, grpmass);
double*& wgrp = TINKER_MOD(group, wgrp);
int& use_group = TINKER_MOD(group, use_group);
int& use_intra = TINKER_MOD(group, use_intra);
int& use_inter = TINKER_MOD(group, use_inter);
#endif
} }
