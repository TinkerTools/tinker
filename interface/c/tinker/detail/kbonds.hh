#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(kbonds, maxnb);
extern int TINKER_MOD(kbonds, maxnb5);
extern int TINKER_MOD(kbonds, maxnb4);
extern int TINKER_MOD(kbonds, maxnb3);
extern int TINKER_MOD(kbonds, maxnel);
extern double* TINKER_MOD(kbonds, bcon);
extern double* TINKER_MOD(kbonds, bcon5);
extern double* TINKER_MOD(kbonds, bcon4);
extern double* TINKER_MOD(kbonds, bcon3);
extern double* TINKER_MOD(kbonds, blen);
extern double* TINKER_MOD(kbonds, blen5);
extern double* TINKER_MOD(kbonds, blen4);
extern double* TINKER_MOD(kbonds, blen3);
extern double* TINKER_MOD(kbonds, dlen);
extern char (*TINKER_MOD(kbonds, kb))[8];
extern char (*TINKER_MOD(kbonds, kb5))[8];
extern char (*TINKER_MOD(kbonds, kb4))[8];
extern char (*TINKER_MOD(kbonds, kb3))[8];
extern char (*TINKER_MOD(kbonds, kel))[12];
#ifdef __cplusplus
}
#endif
