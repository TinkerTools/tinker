#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxele 112
extern double TINKER_MOD(ptable, atmass)[TINKER_MOD__maxele];
extern double TINKER_MOD(ptable, vdwrad)[TINKER_MOD__maxele];
extern double TINKER_MOD(ptable, covrad)[TINKER_MOD__maxele];
extern char TINKER_MOD(ptable, elemnt)[TINKER_MOD__maxele][3];
#ifdef __cplusplus
}
#endif
