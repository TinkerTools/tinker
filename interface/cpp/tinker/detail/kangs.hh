#pragma once

#include "macro.hh"

namespace tinker { namespace kangs {
const int maxna = 2000;
const int maxna5 = 500;
const int maxna4 = 500;
const int maxna3 = 500;
const int maxnap = 2000;
const int maxnaf = 500;
extern double (&acon)[maxna];
extern double (&acon5)[maxna5];
extern double (&acon4)[maxna4];
extern double (&acon3)[maxna3];
extern double (&aconp)[maxnap];
extern double (&aconf)[maxnaf];
extern double (&ang)[maxna][3];
extern double (&ang5)[maxna5][3];
extern double (&ang4)[maxna4][3];
extern double (&ang3)[maxna3][3];
extern double (&angp)[maxnap][2];
extern double (&angf)[maxnaf][2];
extern char (&ka)[maxna][12];
extern char (&ka5)[maxna5][12];
extern char (&ka4)[maxna4][12];
extern char (&ka3)[maxna3][12];
extern char (&kap)[maxnap][12];
extern char (&kaf)[maxnaf][12];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(kangs, acon)[maxna];
extern "C" double TINKER_MOD(kangs, acon5)[maxna5];
extern "C" double TINKER_MOD(kangs, acon4)[maxna4];
extern "C" double TINKER_MOD(kangs, acon3)[maxna3];
extern "C" double TINKER_MOD(kangs, aconp)[maxnap];
extern "C" double TINKER_MOD(kangs, aconf)[maxnaf];
extern "C" double TINKER_MOD(kangs, ang)[maxna][3];
extern "C" double TINKER_MOD(kangs, ang5)[maxna5][3];
extern "C" double TINKER_MOD(kangs, ang4)[maxna4][3];
extern "C" double TINKER_MOD(kangs, ang3)[maxna3][3];
extern "C" double TINKER_MOD(kangs, angp)[maxnap][2];
extern "C" double TINKER_MOD(kangs, angf)[maxnaf][2];
extern "C" char TINKER_MOD(kangs, ka)[maxna][12];
extern "C" char TINKER_MOD(kangs, ka5)[maxna5][12];
extern "C" char TINKER_MOD(kangs, ka4)[maxna4][12];
extern "C" char TINKER_MOD(kangs, ka3)[maxna3][12];
extern "C" char TINKER_MOD(kangs, kap)[maxnap][12];
extern "C" char TINKER_MOD(kangs, kaf)[maxnaf][12];

double (&acon)[maxna] = TINKER_MOD(kangs, acon);
double (&acon5)[maxna5] = TINKER_MOD(kangs, acon5);
double (&acon4)[maxna4] = TINKER_MOD(kangs, acon4);
double (&acon3)[maxna3] = TINKER_MOD(kangs, acon3);
double (&aconp)[maxnap] = TINKER_MOD(kangs, aconp);
double (&aconf)[maxnaf] = TINKER_MOD(kangs, aconf);
double (&ang)[maxna][3] = TINKER_MOD(kangs, ang);
double (&ang5)[maxna5][3] = TINKER_MOD(kangs, ang5);
double (&ang4)[maxna4][3] = TINKER_MOD(kangs, ang4);
double (&ang3)[maxna3][3] = TINKER_MOD(kangs, ang3);
double (&angp)[maxnap][2] = TINKER_MOD(kangs, angp);
double (&angf)[maxnaf][2] = TINKER_MOD(kangs, angf);
char (&ka)[maxna][12] = TINKER_MOD(kangs, ka);
char (&ka5)[maxna5][12] = TINKER_MOD(kangs, ka5);
char (&ka4)[maxna4][12] = TINKER_MOD(kangs, ka4);
char (&ka3)[maxna3][12] = TINKER_MOD(kangs, ka3);
char (&kap)[maxnap][12] = TINKER_MOD(kangs, kap);
char (&kaf)[maxnaf][12] = TINKER_MOD(kangs, kaf);
#endif
} }
