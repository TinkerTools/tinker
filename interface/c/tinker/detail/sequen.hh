#pragma once

#include "macro.hh"
#include "sizes.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(sequen, nseq);
extern int TINKER_MOD(sequen, nchain);
extern int TINKER_MOD(sequen, ichain)[TINKER_MOD__maxres][2];
extern int TINKER_MOD(sequen, seqtyp)[TINKER_MOD__maxres];
extern char TINKER_MOD(sequen, chnnam)[TINKER_MOD__maxres][1];
extern char TINKER_MOD(sequen, seq)[TINKER_MOD__maxres][3];
extern char TINKER_MOD(sequen, chntyp)[TINKER_MOD__maxres][7];
#ifdef __cplusplus
}
#endif
