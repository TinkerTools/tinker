#pragma once

#include "macro.hh"

namespace tinker { namespace shapes {
extern int& maxedge;
extern int& maxtetra;
extern int& npoint;
extern int& nvertex;
extern int& ntetra;
extern int& nnew;
extern int& nfree;
extern int& nkill;
extern int& nlinkfacet;
extern int*& newlist;
extern int*& freespace;
extern int*& killspace;
extern int*& vinfo;
extern int*& tedge;
extern int*& tinfo;
extern int*& tnindex;
extern int*& tetra;
extern int*& tneighbor;
extern int*& linkfacet;
extern int*& linkindex;
extern double& epsln2;
extern double& epsln3;
extern double& epsln4;
extern double& epsln5;
extern double*& crdball;
extern double*& radball;
extern double*& wghtball;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(shapes, maxedge);
extern "C" int TINKER_MOD(shapes, maxtetra);
extern "C" int TINKER_MOD(shapes, npoint);
extern "C" int TINKER_MOD(shapes, nvertex);
extern "C" int TINKER_MOD(shapes, ntetra);
extern "C" int TINKER_MOD(shapes, nnew);
extern "C" int TINKER_MOD(shapes, nfree);
extern "C" int TINKER_MOD(shapes, nkill);
extern "C" int TINKER_MOD(shapes, nlinkfacet);
extern "C" int* TINKER_MOD(shapes, newlist);
extern "C" int* TINKER_MOD(shapes, freespace);
extern "C" int* TINKER_MOD(shapes, killspace);
extern "C" int* TINKER_MOD(shapes, vinfo);
extern "C" int* TINKER_MOD(shapes, tedge);
extern "C" int* TINKER_MOD(shapes, tinfo);
extern "C" int* TINKER_MOD(shapes, tnindex);
extern "C" int* TINKER_MOD(shapes, tetra);
extern "C" int* TINKER_MOD(shapes, tneighbor);
extern "C" int* TINKER_MOD(shapes, linkfacet);
extern "C" int* TINKER_MOD(shapes, linkindex);
extern "C" double TINKER_MOD(shapes, epsln2);
extern "C" double TINKER_MOD(shapes, epsln3);
extern "C" double TINKER_MOD(shapes, epsln4);
extern "C" double TINKER_MOD(shapes, epsln5);
extern "C" double* TINKER_MOD(shapes, crdball);
extern "C" double* TINKER_MOD(shapes, radball);
extern "C" double* TINKER_MOD(shapes, wghtball);

int& maxedge = TINKER_MOD(shapes, maxedge);
int& maxtetra = TINKER_MOD(shapes, maxtetra);
int& npoint = TINKER_MOD(shapes, npoint);
int& nvertex = TINKER_MOD(shapes, nvertex);
int& ntetra = TINKER_MOD(shapes, ntetra);
int& nnew = TINKER_MOD(shapes, nnew);
int& nfree = TINKER_MOD(shapes, nfree);
int& nkill = TINKER_MOD(shapes, nkill);
int& nlinkfacet = TINKER_MOD(shapes, nlinkfacet);
int*& newlist = TINKER_MOD(shapes, newlist);
int*& freespace = TINKER_MOD(shapes, freespace);
int*& killspace = TINKER_MOD(shapes, killspace);
int*& vinfo = TINKER_MOD(shapes, vinfo);
int*& tedge = TINKER_MOD(shapes, tedge);
int*& tinfo = TINKER_MOD(shapes, tinfo);
int*& tnindex = TINKER_MOD(shapes, tnindex);
int*& tetra = TINKER_MOD(shapes, tetra);
int*& tneighbor = TINKER_MOD(shapes, tneighbor);
int*& linkfacet = TINKER_MOD(shapes, linkfacet);
int*& linkindex = TINKER_MOD(shapes, linkindex);
double& epsln2 = TINKER_MOD(shapes, epsln2);
double& epsln3 = TINKER_MOD(shapes, epsln3);
double& epsln4 = TINKER_MOD(shapes, epsln4);
double& epsln5 = TINKER_MOD(shapes, epsln5);
double*& crdball = TINKER_MOD(shapes, crdball);
double*& radball = TINKER_MOD(shapes, radball);
double*& wghtball = TINKER_MOD(shapes, wghtball);
#endif
} }
