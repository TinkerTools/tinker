#pragma once

#include "macro.hh"

namespace tinker { namespace vdwpot {
const int maxgauss = 10;
extern int& ngauss;
extern double (&igauss)[maxgauss][2];
extern double& abuck;
extern double& bbuck;
extern double& cbuck;
extern double& ghal;
extern double& dhal;
extern double& v2scale;
extern double& v3scale;
extern double& v4scale;
extern double& v5scale;
extern int& use_vcorr;
extern char (&vdwindex)[5];
extern char (&radtyp)[5];
extern char (&radsiz)[8];
extern char (&gausstyp)[8];
extern char (&radrule)[10];
extern char (&epsrule)[10];
extern char (&vdwtyp)[13];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(vdwpot, ngauss);
extern "C" double TINKER_MOD(vdwpot, igauss)[maxgauss][2];
extern "C" double TINKER_MOD(vdwpot, abuck);
extern "C" double TINKER_MOD(vdwpot, bbuck);
extern "C" double TINKER_MOD(vdwpot, cbuck);
extern "C" double TINKER_MOD(vdwpot, ghal);
extern "C" double TINKER_MOD(vdwpot, dhal);
extern "C" double TINKER_MOD(vdwpot, v2scale);
extern "C" double TINKER_MOD(vdwpot, v3scale);
extern "C" double TINKER_MOD(vdwpot, v4scale);
extern "C" double TINKER_MOD(vdwpot, v5scale);
extern "C" int TINKER_MOD(vdwpot, use_vcorr);
extern "C" char TINKER_MOD(vdwpot, vdwindex)[5];
extern "C" char TINKER_MOD(vdwpot, radtyp)[5];
extern "C" char TINKER_MOD(vdwpot, radsiz)[8];
extern "C" char TINKER_MOD(vdwpot, gausstyp)[8];
extern "C" char TINKER_MOD(vdwpot, radrule)[10];
extern "C" char TINKER_MOD(vdwpot, epsrule)[10];
extern "C" char TINKER_MOD(vdwpot, vdwtyp)[13];

int& ngauss = TINKER_MOD(vdwpot, ngauss);
double (&igauss)[maxgauss][2] = TINKER_MOD(vdwpot, igauss);
double& abuck = TINKER_MOD(vdwpot, abuck);
double& bbuck = TINKER_MOD(vdwpot, bbuck);
double& cbuck = TINKER_MOD(vdwpot, cbuck);
double& ghal = TINKER_MOD(vdwpot, ghal);
double& dhal = TINKER_MOD(vdwpot, dhal);
double& v2scale = TINKER_MOD(vdwpot, v2scale);
double& v3scale = TINKER_MOD(vdwpot, v3scale);
double& v4scale = TINKER_MOD(vdwpot, v4scale);
double& v5scale = TINKER_MOD(vdwpot, v5scale);
int& use_vcorr = TINKER_MOD(vdwpot, use_vcorr);
char (&vdwindex)[5] = TINKER_MOD(vdwpot, vdwindex);
char (&radtyp)[5] = TINKER_MOD(vdwpot, radtyp);
char (&radsiz)[8] = TINKER_MOD(vdwpot, radsiz);
char (&gausstyp)[8] = TINKER_MOD(vdwpot, gausstyp);
char (&radrule)[10] = TINKER_MOD(vdwpot, radrule);
char (&epsrule)[10] = TINKER_MOD(vdwpot, epsrule);
char (&vdwtyp)[13] = TINKER_MOD(vdwpot, vdwtyp);
#endif
} }
