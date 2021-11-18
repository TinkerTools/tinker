#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern int TINKER_MOD(pme, nfft1);
extern int TINKER_MOD(pme, nfft2);
extern int TINKER_MOD(pme, nfft3);
extern int TINKER_MOD(pme, nefft1);
extern int TINKER_MOD(pme, nefft2);
extern int TINKER_MOD(pme, nefft3);
extern int TINKER_MOD(pme, ndfft1);
extern int TINKER_MOD(pme, ndfft2);
extern int TINKER_MOD(pme, ndfft3);
extern int TINKER_MOD(pme, bsorder);
extern int TINKER_MOD(pme, bseorder);
extern int TINKER_MOD(pme, bsporder);
extern int TINKER_MOD(pme, bsdorder);
extern int* TINKER_MOD(pme, igrid);
extern double* TINKER_MOD(pme, bsmod1);
extern double* TINKER_MOD(pme, bsmod2);
extern double* TINKER_MOD(pme, bsmod3);
extern double* TINKER_MOD(pme, bsbuild);
extern double* TINKER_MOD(pme, thetai1);
extern double* TINKER_MOD(pme, thetai2);
extern double* TINKER_MOD(pme, thetai3);
extern double* TINKER_MOD(pme, qgrid);
extern double* TINKER_MOD(pme, qfac);
#ifdef __cplusplus
}
#endif
