#pragma once

#include "macro.hh"

namespace tinker { namespace solute {
const int maxneck = 45;
extern double& doffset;
extern double& onipr;
extern double& p1;
extern double& p2;
extern double& p3;
extern double& p4;
extern double& p5;
extern double& descoff;
extern double (&rneck)[maxneck];
extern double (&aneck)[maxneck][maxneck];
extern double (&bneck)[maxneck][maxneck];
extern double*& rsolv;
extern double*& rdescr;
extern double*& asolv;
extern double*& rborn;
extern double*& drb;
extern double*& drbp;
extern double*& drobc;
extern double*& gpol;
extern double*& shct;
extern double*& aobc;
extern double*& bobc;
extern double*& gobc;
extern double*& vsolv;
extern double*& wace;
extern double*& s2ace;
extern double*& uace;
extern double*& sneck;
extern double*& bornint;
extern int& useneck;
extern int& usetanh;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(solute, doffset);
extern "C" double TINKER_MOD(solute, onipr);
extern "C" double TINKER_MOD(solute, p1);
extern "C" double TINKER_MOD(solute, p2);
extern "C" double TINKER_MOD(solute, p3);
extern "C" double TINKER_MOD(solute, p4);
extern "C" double TINKER_MOD(solute, p5);
extern "C" double TINKER_MOD(solute, descoff);
extern "C" double TINKER_MOD(solute, rneck)[maxneck];
extern "C" double TINKER_MOD(solute, aneck)[maxneck][maxneck];
extern "C" double TINKER_MOD(solute, bneck)[maxneck][maxneck];
extern "C" double* TINKER_MOD(solute, rsolv);
extern "C" double* TINKER_MOD(solute, rdescr);
extern "C" double* TINKER_MOD(solute, asolv);
extern "C" double* TINKER_MOD(solute, rborn);
extern "C" double* TINKER_MOD(solute, drb);
extern "C" double* TINKER_MOD(solute, drbp);
extern "C" double* TINKER_MOD(solute, drobc);
extern "C" double* TINKER_MOD(solute, gpol);
extern "C" double* TINKER_MOD(solute, shct);
extern "C" double* TINKER_MOD(solute, aobc);
extern "C" double* TINKER_MOD(solute, bobc);
extern "C" double* TINKER_MOD(solute, gobc);
extern "C" double* TINKER_MOD(solute, vsolv);
extern "C" double* TINKER_MOD(solute, wace);
extern "C" double* TINKER_MOD(solute, s2ace);
extern "C" double* TINKER_MOD(solute, uace);
extern "C" double* TINKER_MOD(solute, sneck);
extern "C" double* TINKER_MOD(solute, bornint);
extern "C" int TINKER_MOD(solute, useneck);
extern "C" int TINKER_MOD(solute, usetanh);

double& doffset = TINKER_MOD(solute, doffset);
double& onipr = TINKER_MOD(solute, onipr);
double& p1 = TINKER_MOD(solute, p1);
double& p2 = TINKER_MOD(solute, p2);
double& p3 = TINKER_MOD(solute, p3);
double& p4 = TINKER_MOD(solute, p4);
double& p5 = TINKER_MOD(solute, p5);
double& descoff = TINKER_MOD(solute, descoff);
double (&rneck)[maxneck] = TINKER_MOD(solute, rneck);
double (&aneck)[maxneck][maxneck] = TINKER_MOD(solute, aneck);
double (&bneck)[maxneck][maxneck] = TINKER_MOD(solute, bneck);
double*& rsolv = TINKER_MOD(solute, rsolv);
double*& rdescr = TINKER_MOD(solute, rdescr);
double*& asolv = TINKER_MOD(solute, asolv);
double*& rborn = TINKER_MOD(solute, rborn);
double*& drb = TINKER_MOD(solute, drb);
double*& drbp = TINKER_MOD(solute, drbp);
double*& drobc = TINKER_MOD(solute, drobc);
double*& gpol = TINKER_MOD(solute, gpol);
double*& shct = TINKER_MOD(solute, shct);
double*& aobc = TINKER_MOD(solute, aobc);
double*& bobc = TINKER_MOD(solute, bobc);
double*& gobc = TINKER_MOD(solute, gobc);
double*& vsolv = TINKER_MOD(solute, vsolv);
double*& wace = TINKER_MOD(solute, wace);
double*& s2ace = TINKER_MOD(solute, s2ace);
double*& uace = TINKER_MOD(solute, uace);
double*& sneck = TINKER_MOD(solute, sneck);
double*& bornint = TINKER_MOD(solute, bornint);
int& useneck = TINKER_MOD(solute, useneck);
int& usetanh = TINKER_MOD(solute, usetanh);
#endif
} }
