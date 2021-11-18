#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__maxprime 15
extern int TINKER_MOD(fft, iprime)[3][TINKER_MOD__maxprime];
extern unsigned long long TINKER_MOD(fft, planf);
extern unsigned long long TINKER_MOD(fft, planb);
extern double* TINKER_MOD(fft, ffttable);
extern char TINKER_MOD(fft, ffttyp)[7];
#ifdef __cplusplus
}
#endif
