#pragma once

#include "macro.hh"

namespace tinker { namespace polpot {
extern int& politer;
extern double& poleps;
extern double& p2scale;
extern double& p3scale;
extern double& p4scale;
extern double& p5scale;
extern double& p2iscale;
extern double& p3iscale;
extern double& p4iscale;
extern double& p5iscale;
extern double& d1scale;
extern double& d2scale;
extern double& d3scale;
extern double& d4scale;
extern double& u1scale;
extern double& u2scale;
extern double& u3scale;
extern double& u4scale;
extern double& w2scale;
extern double& w3scale;
extern double& w4scale;
extern double& w5scale;
extern double& uaccel;
extern int& polprt;
extern int& dpequal;
extern int& use_thole;
extern int& use_tholed;
extern int& use_expol;
extern char (&scrtyp)[3];
extern char (&poltyp)[6];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(polpot, politer);
extern "C" double TINKER_MOD(polpot, poleps);
extern "C" double TINKER_MOD(polpot, p2scale);
extern "C" double TINKER_MOD(polpot, p3scale);
extern "C" double TINKER_MOD(polpot, p4scale);
extern "C" double TINKER_MOD(polpot, p5scale);
extern "C" double TINKER_MOD(polpot, p2iscale);
extern "C" double TINKER_MOD(polpot, p3iscale);
extern "C" double TINKER_MOD(polpot, p4iscale);
extern "C" double TINKER_MOD(polpot, p5iscale);
extern "C" double TINKER_MOD(polpot, d1scale);
extern "C" double TINKER_MOD(polpot, d2scale);
extern "C" double TINKER_MOD(polpot, d3scale);
extern "C" double TINKER_MOD(polpot, d4scale);
extern "C" double TINKER_MOD(polpot, u1scale);
extern "C" double TINKER_MOD(polpot, u2scale);
extern "C" double TINKER_MOD(polpot, u3scale);
extern "C" double TINKER_MOD(polpot, u4scale);
extern "C" double TINKER_MOD(polpot, w2scale);
extern "C" double TINKER_MOD(polpot, w3scale);
extern "C" double TINKER_MOD(polpot, w4scale);
extern "C" double TINKER_MOD(polpot, w5scale);
extern "C" double TINKER_MOD(polpot, uaccel);
extern "C" int TINKER_MOD(polpot, polprt);
extern "C" int TINKER_MOD(polpot, dpequal);
extern "C" int TINKER_MOD(polpot, use_thole);
extern "C" int TINKER_MOD(polpot, use_tholed);
extern "C" int TINKER_MOD(polpot, use_expol);
extern "C" char TINKER_MOD(polpot, scrtyp)[3];
extern "C" char TINKER_MOD(polpot, poltyp)[6];

int& politer = TINKER_MOD(polpot, politer);
double& poleps = TINKER_MOD(polpot, poleps);
double& p2scale = TINKER_MOD(polpot, p2scale);
double& p3scale = TINKER_MOD(polpot, p3scale);
double& p4scale = TINKER_MOD(polpot, p4scale);
double& p5scale = TINKER_MOD(polpot, p5scale);
double& p2iscale = TINKER_MOD(polpot, p2iscale);
double& p3iscale = TINKER_MOD(polpot, p3iscale);
double& p4iscale = TINKER_MOD(polpot, p4iscale);
double& p5iscale = TINKER_MOD(polpot, p5iscale);
double& d1scale = TINKER_MOD(polpot, d1scale);
double& d2scale = TINKER_MOD(polpot, d2scale);
double& d3scale = TINKER_MOD(polpot, d3scale);
double& d4scale = TINKER_MOD(polpot, d4scale);
double& u1scale = TINKER_MOD(polpot, u1scale);
double& u2scale = TINKER_MOD(polpot, u2scale);
double& u3scale = TINKER_MOD(polpot, u3scale);
double& u4scale = TINKER_MOD(polpot, u4scale);
double& w2scale = TINKER_MOD(polpot, w2scale);
double& w3scale = TINKER_MOD(polpot, w3scale);
double& w4scale = TINKER_MOD(polpot, w4scale);
double& w5scale = TINKER_MOD(polpot, w5scale);
double& uaccel = TINKER_MOD(polpot, uaccel);
int& polprt = TINKER_MOD(polpot, polprt);
int& dpequal = TINKER_MOD(polpot, dpequal);
int& use_thole = TINKER_MOD(polpot, use_thole);
int& use_tholed = TINKER_MOD(polpot, use_tholed);
int& use_expol = TINKER_MOD(polpot, use_expol);
char (&scrtyp)[3] = TINKER_MOD(polpot, scrtyp);
char (&poltyp)[6] = TINKER_MOD(polpot, poltyp);
#endif
} }
