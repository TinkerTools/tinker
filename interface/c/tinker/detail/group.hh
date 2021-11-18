#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(group, ngrp);
extern int* TINKER_MOD(group, kgrp);
extern int* TINKER_MOD(group, grplist);
extern int* TINKER_MOD(group, igrp);
extern double* TINKER_MOD(group, grpmass);
extern double* TINKER_MOD(group, wgrp);
extern int TINKER_MOD(group, use_group);
extern int TINKER_MOD(group, use_intra);
extern int TINKER_MOD(group, use_inter);
#ifdef __cplusplus
}
#endif
