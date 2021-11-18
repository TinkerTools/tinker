#pragma once

#include "macro.hh"

namespace tinker { namespace bath {
const int maxnose = 4;
extern int& voltrial;
extern double& kelvin;
extern double& atmsph;
extern double& tautemp;
extern double& taupres;
extern double& compress;
extern double& collide;
extern double& eta;
extern double& volmove;
extern double& vbar;
extern double& qbar;
extern double& gbar;
extern double (&vnh)[maxnose];
extern double (&qnh)[maxnose];
extern double (&gnh)[maxnose];
extern int& isothermal;
extern int& isobaric;
extern int& anisotrop;
extern char (&volscale)[9];
extern char (&barostat)[11];
extern char (&thermostat)[11];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(bath, voltrial);
extern "C" double TINKER_MOD(bath, kelvin);
extern "C" double TINKER_MOD(bath, atmsph);
extern "C" double TINKER_MOD(bath, tautemp);
extern "C" double TINKER_MOD(bath, taupres);
extern "C" double TINKER_MOD(bath, compress);
extern "C" double TINKER_MOD(bath, collide);
extern "C" double TINKER_MOD(bath, eta);
extern "C" double TINKER_MOD(bath, volmove);
extern "C" double TINKER_MOD(bath, vbar);
extern "C" double TINKER_MOD(bath, qbar);
extern "C" double TINKER_MOD(bath, gbar);
extern "C" double TINKER_MOD(bath, vnh)[maxnose];
extern "C" double TINKER_MOD(bath, qnh)[maxnose];
extern "C" double TINKER_MOD(bath, gnh)[maxnose];
extern "C" int TINKER_MOD(bath, isothermal);
extern "C" int TINKER_MOD(bath, isobaric);
extern "C" int TINKER_MOD(bath, anisotrop);
extern "C" char TINKER_MOD(bath, volscale)[9];
extern "C" char TINKER_MOD(bath, barostat)[11];
extern "C" char TINKER_MOD(bath, thermostat)[11];

int& voltrial = TINKER_MOD(bath, voltrial);
double& kelvin = TINKER_MOD(bath, kelvin);
double& atmsph = TINKER_MOD(bath, atmsph);
double& tautemp = TINKER_MOD(bath, tautemp);
double& taupres = TINKER_MOD(bath, taupres);
double& compress = TINKER_MOD(bath, compress);
double& collide = TINKER_MOD(bath, collide);
double& eta = TINKER_MOD(bath, eta);
double& volmove = TINKER_MOD(bath, volmove);
double& vbar = TINKER_MOD(bath, vbar);
double& qbar = TINKER_MOD(bath, qbar);
double& gbar = TINKER_MOD(bath, gbar);
double (&vnh)[maxnose] = TINKER_MOD(bath, vnh);
double (&qnh)[maxnose] = TINKER_MOD(bath, qnh);
double (&gnh)[maxnose] = TINKER_MOD(bath, gnh);
int& isothermal = TINKER_MOD(bath, isothermal);
int& isobaric = TINKER_MOD(bath, isobaric);
int& anisotrop = TINKER_MOD(bath, anisotrop);
char (&volscale)[9] = TINKER_MOD(bath, volscale);
char (&barostat)[11] = TINKER_MOD(bath, barostat);
char (&thermostat)[11] = TINKER_MOD(bath, thermostat);
#endif
} }
