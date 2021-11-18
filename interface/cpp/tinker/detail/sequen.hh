#pragma once

#include "macro.hh"
#include "sizes.hh"

namespace tinker { namespace sequen {
using namespace sizes;

extern int& nseq;
extern int& nchain;
extern int (&ichain)[maxres][2];
extern int (&seqtyp)[maxres];
extern char (&chnnam)[maxres][1];
extern char (&seq)[maxres][3];
extern char (&chntyp)[maxres][7];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(sequen, nseq);
extern "C" int TINKER_MOD(sequen, nchain);
extern "C" int TINKER_MOD(sequen, ichain)[maxres][2];
extern "C" int TINKER_MOD(sequen, seqtyp)[maxres];
extern "C" char TINKER_MOD(sequen, chnnam)[maxres][1];
extern "C" char TINKER_MOD(sequen, seq)[maxres][3];
extern "C" char TINKER_MOD(sequen, chntyp)[maxres][7];

int& nseq = TINKER_MOD(sequen, nseq);
int& nchain = TINKER_MOD(sequen, nchain);
int (&ichain)[maxres][2] = TINKER_MOD(sequen, ichain);
int (&seqtyp)[maxres] = TINKER_MOD(sequen, seqtyp);
char (&chnnam)[maxres][1] = TINKER_MOD(sequen, chnnam);
char (&seq)[maxres][3] = TINKER_MOD(sequen, seq);
char (&chntyp)[maxres][7] = TINKER_MOD(sequen, chntyp);
#endif
} }
