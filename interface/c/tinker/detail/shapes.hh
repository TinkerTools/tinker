#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(shapes, maxedge);
extern int TINKER_MOD(shapes, maxtetra);
extern int TINKER_MOD(shapes, npoint);
extern int TINKER_MOD(shapes, nvertex);
extern int TINKER_MOD(shapes, ntetra);
extern int TINKER_MOD(shapes, nnew);
extern int TINKER_MOD(shapes, nfree);
extern int TINKER_MOD(shapes, nkill);
extern int TINKER_MOD(shapes, nlinkfacet);
extern int* TINKER_MOD(shapes, newlist);
extern int* TINKER_MOD(shapes, freespace);
extern int* TINKER_MOD(shapes, killspace);
extern int* TINKER_MOD(shapes, vinfo);
extern int* TINKER_MOD(shapes, tedge);
extern int* TINKER_MOD(shapes, tinfo);
extern int* TINKER_MOD(shapes, tnindex);
extern int* TINKER_MOD(shapes, tetra);
extern int* TINKER_MOD(shapes, tneighbor);
extern int* TINKER_MOD(shapes, linkfacet);
extern int* TINKER_MOD(shapes, linkindex);
extern double TINKER_MOD(shapes, epsln2);
extern double TINKER_MOD(shapes, epsln3);
extern double TINKER_MOD(shapes, epsln4);
extern double TINKER_MOD(shapes, epsln5);
extern double* TINKER_MOD(shapes, crdball);
extern double* TINKER_MOD(shapes, radball);
extern double* TINKER_MOD(shapes, wghtball);
#ifdef __cplusplus
}
#endif
