#pragma once

#include "macro.hh"

namespace tinker { namespace boxes {
extern double& xbox;
extern double& ybox;
extern double& zbox;
extern double& alpha;
extern double& beta;
extern double& gamma;
extern double& xbox2;
extern double& ybox2;
extern double& zbox2;
extern double& box34;
extern double& volbox;
extern double& alpha_sin;
extern double& alpha_cos;
extern double& beta_sin;
extern double& beta_cos;
extern double& gamma_sin;
extern double& gamma_cos;
extern double& beta_term;
extern double& gamma_term;
extern double (&lvec)[3][3];
extern double (&recip)[3][3];
extern int& orthogonal;
extern int& monoclinic;
extern int& triclinic;
extern int& octahedron;
extern int& dodecadron;
extern int& nonprism;
extern char (&spacegrp)[10];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(boxes, xbox);
extern "C" double TINKER_MOD(boxes, ybox);
extern "C" double TINKER_MOD(boxes, zbox);
extern "C" double TINKER_MOD(boxes, alpha);
extern "C" double TINKER_MOD(boxes, beta);
extern "C" double TINKER_MOD(boxes, gamma);
extern "C" double TINKER_MOD(boxes, xbox2);
extern "C" double TINKER_MOD(boxes, ybox2);
extern "C" double TINKER_MOD(boxes, zbox2);
extern "C" double TINKER_MOD(boxes, box34);
extern "C" double TINKER_MOD(boxes, volbox);
extern "C" double TINKER_MOD(boxes, alpha_sin);
extern "C" double TINKER_MOD(boxes, alpha_cos);
extern "C" double TINKER_MOD(boxes, beta_sin);
extern "C" double TINKER_MOD(boxes, beta_cos);
extern "C" double TINKER_MOD(boxes, gamma_sin);
extern "C" double TINKER_MOD(boxes, gamma_cos);
extern "C" double TINKER_MOD(boxes, beta_term);
extern "C" double TINKER_MOD(boxes, gamma_term);
extern "C" double TINKER_MOD(boxes, lvec)[3][3];
extern "C" double TINKER_MOD(boxes, recip)[3][3];
extern "C" int TINKER_MOD(boxes, orthogonal);
extern "C" int TINKER_MOD(boxes, monoclinic);
extern "C" int TINKER_MOD(boxes, triclinic);
extern "C" int TINKER_MOD(boxes, octahedron);
extern "C" int TINKER_MOD(boxes, dodecadron);
extern "C" int TINKER_MOD(boxes, nonprism);
extern "C" char TINKER_MOD(boxes, spacegrp)[10];

double& xbox = TINKER_MOD(boxes, xbox);
double& ybox = TINKER_MOD(boxes, ybox);
double& zbox = TINKER_MOD(boxes, zbox);
double& alpha = TINKER_MOD(boxes, alpha);
double& beta = TINKER_MOD(boxes, beta);
double& gamma = TINKER_MOD(boxes, gamma);
double& xbox2 = TINKER_MOD(boxes, xbox2);
double& ybox2 = TINKER_MOD(boxes, ybox2);
double& zbox2 = TINKER_MOD(boxes, zbox2);
double& box34 = TINKER_MOD(boxes, box34);
double& volbox = TINKER_MOD(boxes, volbox);
double& alpha_sin = TINKER_MOD(boxes, alpha_sin);
double& alpha_cos = TINKER_MOD(boxes, alpha_cos);
double& beta_sin = TINKER_MOD(boxes, beta_sin);
double& beta_cos = TINKER_MOD(boxes, beta_cos);
double& gamma_sin = TINKER_MOD(boxes, gamma_sin);
double& gamma_cos = TINKER_MOD(boxes, gamma_cos);
double& beta_term = TINKER_MOD(boxes, beta_term);
double& gamma_term = TINKER_MOD(boxes, gamma_term);
double (&lvec)[3][3] = TINKER_MOD(boxes, lvec);
double (&recip)[3][3] = TINKER_MOD(boxes, recip);
int& orthogonal = TINKER_MOD(boxes, orthogonal);
int& monoclinic = TINKER_MOD(boxes, monoclinic);
int& triclinic = TINKER_MOD(boxes, triclinic);
int& octahedron = TINKER_MOD(boxes, octahedron);
int& dodecadron = TINKER_MOD(boxes, dodecadron);
int& nonprism = TINKER_MOD(boxes, nonprism);
char (&spacegrp)[10] = TINKER_MOD(boxes, spacegrp);
#endif
} }
