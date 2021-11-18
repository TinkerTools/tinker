#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(neigh, maxvlst);
extern int TINKER_MOD(neigh, maxelst);
extern int TINKER_MOD(neigh, maxulst);
extern int* TINKER_MOD(neigh, nvlst);
extern int* TINKER_MOD(neigh, vlst);
extern int* TINKER_MOD(neigh, nelst);
extern int* TINKER_MOD(neigh, elst);
extern int* TINKER_MOD(neigh, nulst);
extern int* TINKER_MOD(neigh, ulst);
extern double TINKER_MOD(neigh, lbuffer);
extern double TINKER_MOD(neigh, pbuffer);
extern double TINKER_MOD(neigh, lbuf2);
extern double TINKER_MOD(neigh, pbuf2);
extern double TINKER_MOD(neigh, vbuf2);
extern double TINKER_MOD(neigh, vbufx);
extern double TINKER_MOD(neigh, dbuf2);
extern double TINKER_MOD(neigh, dbufx);
extern double TINKER_MOD(neigh, cbuf2);
extern double TINKER_MOD(neigh, cbufx);
extern double TINKER_MOD(neigh, mbuf2);
extern double TINKER_MOD(neigh, mbufx);
extern double TINKER_MOD(neigh, ubuf2);
extern double TINKER_MOD(neigh, ubufx);
extern double* TINKER_MOD(neigh, xvold);
extern double* TINKER_MOD(neigh, yvold);
extern double* TINKER_MOD(neigh, zvold);
extern double* TINKER_MOD(neigh, xeold);
extern double* TINKER_MOD(neigh, yeold);
extern double* TINKER_MOD(neigh, zeold);
extern double* TINKER_MOD(neigh, xuold);
extern double* TINKER_MOD(neigh, yuold);
extern double* TINKER_MOD(neigh, zuold);
extern int TINKER_MOD(neigh, dovlst);
extern int TINKER_MOD(neigh, dodlst);
extern int TINKER_MOD(neigh, doclst);
extern int TINKER_MOD(neigh, domlst);
extern int TINKER_MOD(neigh, doulst);
#ifdef __cplusplus
}
#endif
