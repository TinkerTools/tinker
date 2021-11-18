#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxgauss 10
extern int TINKER_MOD(vdwpot, ngauss);
extern double TINKER_MOD(vdwpot, igauss)[TINKER_MOD__maxgauss][2];
extern double TINKER_MOD(vdwpot, abuck);
extern double TINKER_MOD(vdwpot, bbuck);
extern double TINKER_MOD(vdwpot, cbuck);
extern double TINKER_MOD(vdwpot, ghal);
extern double TINKER_MOD(vdwpot, dhal);
extern double TINKER_MOD(vdwpot, v2scale);
extern double TINKER_MOD(vdwpot, v3scale);
extern double TINKER_MOD(vdwpot, v4scale);
extern double TINKER_MOD(vdwpot, v5scale);
extern int TINKER_MOD(vdwpot, use_vcorr);
extern char TINKER_MOD(vdwpot, vdwindex)[5];
extern char TINKER_MOD(vdwpot, radtyp)[5];
extern char TINKER_MOD(vdwpot, radsiz)[8];
extern char TINKER_MOD(vdwpot, gausstyp)[8];
extern char TINKER_MOD(vdwpot, radrule)[10];
extern char TINKER_MOD(vdwpot, epsrule)[10];
extern char TINKER_MOD(vdwpot, vdwtyp)[13];
#ifdef __cplusplus
}
#endif
