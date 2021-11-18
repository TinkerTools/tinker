#pragma once

#include "macro.hh"

namespace tinker { namespace kbonds {
const int maxnb = 2000;
const int maxnb5 = 500;
const int maxnb4 = 500;
const int maxnb3 = 500;
const int maxnel = 500;
extern double (&bcon)[maxnb];
extern double (&bcon5)[maxnb5];
extern double (&bcon4)[maxnb4];
extern double (&bcon3)[maxnb3];
extern double (&blen)[maxnb];
extern double (&blen5)[maxnb5];
extern double (&blen4)[maxnb4];
extern double (&blen3)[maxnb3];
extern double (&dlen)[maxnel];
extern char (&kb)[maxnb][8];
extern char (&kb5)[maxnb5][8];
extern char (&kb4)[maxnb4][8];
extern char (&kb3)[maxnb3][8];
extern char (&kel)[maxnel][12];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(kbonds, bcon)[maxnb];
extern "C" double TINKER_MOD(kbonds, bcon5)[maxnb5];
extern "C" double TINKER_MOD(kbonds, bcon4)[maxnb4];
extern "C" double TINKER_MOD(kbonds, bcon3)[maxnb3];
extern "C" double TINKER_MOD(kbonds, blen)[maxnb];
extern "C" double TINKER_MOD(kbonds, blen5)[maxnb5];
extern "C" double TINKER_MOD(kbonds, blen4)[maxnb4];
extern "C" double TINKER_MOD(kbonds, blen3)[maxnb3];
extern "C" double TINKER_MOD(kbonds, dlen)[maxnel];
extern "C" char TINKER_MOD(kbonds, kb)[maxnb][8];
extern "C" char TINKER_MOD(kbonds, kb5)[maxnb5][8];
extern "C" char TINKER_MOD(kbonds, kb4)[maxnb4][8];
extern "C" char TINKER_MOD(kbonds, kb3)[maxnb3][8];
extern "C" char TINKER_MOD(kbonds, kel)[maxnel][12];

double (&bcon)[maxnb] = TINKER_MOD(kbonds, bcon);
double (&bcon5)[maxnb5] = TINKER_MOD(kbonds, bcon5);
double (&bcon4)[maxnb4] = TINKER_MOD(kbonds, bcon4);
double (&bcon3)[maxnb3] = TINKER_MOD(kbonds, bcon3);
double (&blen)[maxnb] = TINKER_MOD(kbonds, blen);
double (&blen5)[maxnb5] = TINKER_MOD(kbonds, blen5);
double (&blen4)[maxnb4] = TINKER_MOD(kbonds, blen4);
double (&blen3)[maxnb3] = TINKER_MOD(kbonds, blen3);
double (&dlen)[maxnel] = TINKER_MOD(kbonds, dlen);
char (&kb)[maxnb][8] = TINKER_MOD(kbonds, kb);
char (&kb5)[maxnb5][8] = TINKER_MOD(kbonds, kb5);
char (&kb4)[maxnb4][8] = TINKER_MOD(kbonds, kb4);
char (&kb3)[maxnb3][8] = TINKER_MOD(kbonds, kb3);
char (&kel)[maxnel][12] = TINKER_MOD(kbonds, kel);
#endif
} }
