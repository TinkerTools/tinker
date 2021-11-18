#pragma once

#include "macro.hh"

namespace tinker { namespace solute {
extern double& doffset;
extern double& onipr;
extern double& p1;
extern double& p2;
extern double& p3;
extern double& p4;
extern double& p5;
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

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(solute, doffset);
extern "C" double TINKER_MOD(solute, onipr);
extern "C" double TINKER_MOD(solute, p1);
extern "C" double TINKER_MOD(solute, p2);
extern "C" double TINKER_MOD(solute, p3);
extern "C" double TINKER_MOD(solute, p4);
extern "C" double TINKER_MOD(solute, p5);
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

double& doffset = TINKER_MOD(solute, doffset);
double& onipr = TINKER_MOD(solute, onipr);
double& p1 = TINKER_MOD(solute, p1);
double& p2 = TINKER_MOD(solute, p2);
double& p3 = TINKER_MOD(solute, p3);
double& p4 = TINKER_MOD(solute, p4);
double& p5 = TINKER_MOD(solute, p5);
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
#endif
} }
