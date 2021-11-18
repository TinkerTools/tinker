#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
#define TINKER_MOD__rcarbon 1.7e0
#define TINKER_MOD__rwater 1.4e0
#define TINKER_MOD__acsurf 120.7628e0
#define TINKER_MOD__safact 0.3516e0
#define TINKER_MOD__tgrad 100.0e0
#define TINKER_MOD__toffset 6.0e0
#define TINKER_MOD__hpmfcut 11.0e0
#define TINKER_MOD__hd1 -0.7308004860404441194e0
#define TINKER_MOD__hd2 0.2001645051578760659e0
#define TINKER_MOD__hd3 -0.0905499953418473502e0
#define TINKER_MOD__hc1 3.8167879266271396155e0
#define TINKER_MOD__hc2 5.4669162286016419472e0
#define TINKER_MOD__hc3 7.1167694861385353278e0
#define TINKER_MOD__hw1 1.6858993102248638341e0
#define TINKER_MOD__hw2 1.3906405621629980285e0
#define TINKER_MOD__hw3 1.5741657341338335385e0
extern int TINKER_MOD(hpmf, npmf);
extern int* TINKER_MOD(hpmf, ipmf);
extern double* TINKER_MOD(hpmf, rpmf);
extern double* TINKER_MOD(hpmf, acsa);
#ifdef __cplusplus
}
#endif
