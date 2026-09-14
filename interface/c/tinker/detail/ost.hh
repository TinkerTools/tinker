#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(ost, fli0);
extern int TINKER_MOD(ost, nflmda);
extern int TINKER_MOD(ost, nmetahist);
extern int TINKER_MOD(ost, nmethistsave);
extern int TINKER_MOD(ost, ostcvbin);
extern int TINKER_MOD(ost, sizemetahist);
extern int* TINKER_MOD(ost, metaihist);
extern int* TINKER_MOD(ost, osthist);
extern int* TINKER_MOD(ost, ostnext);
extern int* TINKER_MOD(ost, osthead);
extern double TINKER_MOD(ost, hbias);
extern double TINKER_MOD(ost, maxwfhist);
extern double TINKER_MOD(ost, maxwlhist);
extern double TINKER_MOD(ost, ostbdgdfl);
extern double TINKER_MOD(ost, ostbdgdl);
extern double TINKER_MOD(ost, ostcvdif);
extern double TINKER_MOD(ost, ostcvrat);
extern double TINKER_MOD(ost, ostcvslp);
extern double TINKER_MOD(ost, ostcvstd);
extern double TINKER_MOD(ost, ostdedlslp);
extern double TINKER_MOD(ost, ostdgdl);
extern double TINKER_MOD(ost, ostgtempgamma);
extern double TINKER_MOD(ost, ostgthresh);
extern double TINKER_MOD(ost, ostlambdaslp);
extern double TINKER_MOD(ost, ostltempgamma);
extern double TINKER_MOD(ost, ostlthresh);
extern double TINKER_MOD(ost, oststdev);
extern double TINKER_MOD(ost, wfhist);
extern double TINKER_MOD(ost, wflmda);
extern double TINKER_MOD(ost, wflmda2);
extern double TINKER_MOD(ost, wlhist);
extern double* TINKER_MOD(ost, dvmetagrid);
extern double* TINKER_MOD(ost, metahhist);
extern double* TINKER_MOD(ost, metalhist);
extern double* TINKER_MOD(ost, metawhist);
extern double* TINKER_MOD(ost, ostdedlavgbin);
extern double* TINKER_MOD(ost, ostdedlslpbin);
extern double* TINKER_MOD(ost, ostdedlstdbin);
extern double* TINKER_MOD(ost, osthhist);
extern double* TINKER_MOD(ost, ostlmdaavgbin);
extern double* TINKER_MOD(ost, ostlmdaslpbin);
extern double* TINKER_MOD(ost, ostlmdastdbin);
extern double* TINKER_MOD(ost, ostwfhist);
extern double* TINKER_MOD(ost, ostwlhist);
extern double* TINKER_MOD(ost, vkernelmax);
extern double* TINKER_MOD(ost, vmetagrid);
extern double* TINKER_MOD(ost, gfkernel);
extern double* TINKER_MOD(ost, gkernel);
extern double* TINKER_MOD(ost, glfkernel);
extern double* TINKER_MOD(ost, glkernel);
extern int TINKER_MOD(ost, fastkernel);
extern int TINKER_MOD(ost, ostinterpol);
extern int TINKER_MOD(ost, use_ostgtemp);
extern int TINKER_MOD(ost, use_ostltemp);
extern char TINKER_MOD(ost, ostlabel)[40];
extern char TINKER_MOD(ost, osttitle)[40];
extern char TINKER_MOD(ost, metasavefile)[240];
#ifdef __cplusplus
}
#endif
