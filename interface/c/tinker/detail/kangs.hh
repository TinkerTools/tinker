#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxna 2000
#define TINKER_MOD__maxna5 500
#define TINKER_MOD__maxna4 500
#define TINKER_MOD__maxna3 500
#define TINKER_MOD__maxnap 2000
#define TINKER_MOD__maxnaf 500
extern double TINKER_MOD(kangs, acon)[TINKER_MOD__maxna];
extern double TINKER_MOD(kangs, acon5)[TINKER_MOD__maxna5];
extern double TINKER_MOD(kangs, acon4)[TINKER_MOD__maxna4];
extern double TINKER_MOD(kangs, acon3)[TINKER_MOD__maxna3];
extern double TINKER_MOD(kangs, aconp)[TINKER_MOD__maxnap];
extern double TINKER_MOD(kangs, aconf)[TINKER_MOD__maxnaf];
extern double TINKER_MOD(kangs, ang)[TINKER_MOD__maxna][3];
extern double TINKER_MOD(kangs, ang5)[TINKER_MOD__maxna5][3];
extern double TINKER_MOD(kangs, ang4)[TINKER_MOD__maxna4][3];
extern double TINKER_MOD(kangs, ang3)[TINKER_MOD__maxna3][3];
extern double TINKER_MOD(kangs, angp)[TINKER_MOD__maxnap][2];
extern double TINKER_MOD(kangs, angf)[TINKER_MOD__maxnaf][2];
extern char TINKER_MOD(kangs, ka)[TINKER_MOD__maxna][12];
extern char TINKER_MOD(kangs, ka5)[TINKER_MOD__maxna5][12];
extern char TINKER_MOD(kangs, ka4)[TINKER_MOD__maxna4][12];
extern char TINKER_MOD(kangs, ka3)[TINKER_MOD__maxna3][12];
extern char TINKER_MOD(kangs, kap)[TINKER_MOD__maxnap][12];
extern char TINKER_MOD(kangs, kaf)[TINKER_MOD__maxnaf][12];
#ifdef __cplusplus
}
#endif
