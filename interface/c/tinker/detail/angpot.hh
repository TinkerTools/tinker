#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(angpot, angunit);
extern double TINKER_MOD(angpot, stbnunit);
extern double TINKER_MOD(angpot, aaunit);
extern double TINKER_MOD(angpot, opbunit);
extern double TINKER_MOD(angpot, opdunit);
extern double TINKER_MOD(angpot, cang);
extern double TINKER_MOD(angpot, qang);
extern double TINKER_MOD(angpot, pang);
extern double TINKER_MOD(angpot, sang);
extern double TINKER_MOD(angpot, copb);
extern double TINKER_MOD(angpot, qopb);
extern double TINKER_MOD(angpot, popb);
extern double TINKER_MOD(angpot, sopb);
extern double TINKER_MOD(angpot, copd);
extern double TINKER_MOD(angpot, qopd);
extern double TINKER_MOD(angpot, popd);
extern double TINKER_MOD(angpot, sopd);
extern char TINKER_MOD(angpot, opbtyp)[8];
extern char (*TINKER_MOD(angpot, angtyp))[8];
#ifdef __cplusplus
}
#endif
