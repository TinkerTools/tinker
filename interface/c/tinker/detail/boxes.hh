#pragma once

#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif
extern double TINKER_MOD(boxes, xbox);
extern double TINKER_MOD(boxes, ybox);
extern double TINKER_MOD(boxes, zbox);
extern double TINKER_MOD(boxes, alpha);
extern double TINKER_MOD(boxes, beta);
extern double TINKER_MOD(boxes, gamma);
extern double TINKER_MOD(boxes, xbox2);
extern double TINKER_MOD(boxes, ybox2);
extern double TINKER_MOD(boxes, zbox2);
extern double TINKER_MOD(boxes, box34);
extern double TINKER_MOD(boxes, volbox);
extern double TINKER_MOD(boxes, alpha_sin);
extern double TINKER_MOD(boxes, alpha_cos);
extern double TINKER_MOD(boxes, beta_sin);
extern double TINKER_MOD(boxes, beta_cos);
extern double TINKER_MOD(boxes, gamma_sin);
extern double TINKER_MOD(boxes, gamma_cos);
extern double TINKER_MOD(boxes, beta_term);
extern double TINKER_MOD(boxes, gamma_term);
extern double TINKER_MOD(boxes, lvec)[3][3];
extern double TINKER_MOD(boxes, recip)[3][3];
extern int TINKER_MOD(boxes, orthogonal);
extern int TINKER_MOD(boxes, monoclinic);
extern int TINKER_MOD(boxes, triclinic);
extern int TINKER_MOD(boxes, octahedron);
extern int TINKER_MOD(boxes, dodecadron);
extern int TINKER_MOD(boxes, nonprism);
extern int TINKER_MOD(boxes, nosymm);
extern char TINKER_MOD(boxes, spacegrp)[10];
#ifdef __cplusplus
}
#endif
