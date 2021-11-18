#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxnb 2000
#define TINKER_MOD__maxnb5 500
#define TINKER_MOD__maxnb4 500
#define TINKER_MOD__maxnb3 500
#define TINKER_MOD__maxnel 500
extern double TINKER_MOD(kbonds, bcon)[TINKER_MOD__maxnb];
extern double TINKER_MOD(kbonds, bcon5)[TINKER_MOD__maxnb5];
extern double TINKER_MOD(kbonds, bcon4)[TINKER_MOD__maxnb4];
extern double TINKER_MOD(kbonds, bcon3)[TINKER_MOD__maxnb3];
extern double TINKER_MOD(kbonds, blen)[TINKER_MOD__maxnb];
extern double TINKER_MOD(kbonds, blen5)[TINKER_MOD__maxnb5];
extern double TINKER_MOD(kbonds, blen4)[TINKER_MOD__maxnb4];
extern double TINKER_MOD(kbonds, blen3)[TINKER_MOD__maxnb3];
extern double TINKER_MOD(kbonds, dlen)[TINKER_MOD__maxnel];
extern char TINKER_MOD(kbonds, kb)[TINKER_MOD__maxnb][8];
extern char TINKER_MOD(kbonds, kb5)[TINKER_MOD__maxnb5][8];
extern char TINKER_MOD(kbonds, kb4)[TINKER_MOD__maxnb4][8];
extern char TINKER_MOD(kbonds, kb3)[TINKER_MOD__maxnb3][8];
extern char TINKER_MOD(kbonds, kel)[TINKER_MOD__maxnel][12];
#ifdef __cplusplus
}
#endif
