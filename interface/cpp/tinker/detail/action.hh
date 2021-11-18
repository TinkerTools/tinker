#pragma once

#include "macro.hh"

namespace tinker { namespace action {
extern int& neb;
extern int& nea;
extern int& neba;
extern int& neub;
extern int& neaa;
extern int& neopb;
extern int& neopd;
extern int& neid;
extern int& neit;
extern int& net;
extern int& nept;
extern int& nebt;
extern int& neat;
extern int& nett;
extern int& nev;
extern int& ner;
extern int& nedsp;
extern int& nec;
extern int& necd;
extern int& ned;
extern int& nem;
extern int& nep;
extern int& nect;
extern int& new_;
extern int& nerxf;
extern int& nes;
extern int& nelf;
extern int& neg;
extern int& nex;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(action, neb);
extern "C" int TINKER_MOD(action, nea);
extern "C" int TINKER_MOD(action, neba);
extern "C" int TINKER_MOD(action, neub);
extern "C" int TINKER_MOD(action, neaa);
extern "C" int TINKER_MOD(action, neopb);
extern "C" int TINKER_MOD(action, neopd);
extern "C" int TINKER_MOD(action, neid);
extern "C" int TINKER_MOD(action, neit);
extern "C" int TINKER_MOD(action, net);
extern "C" int TINKER_MOD(action, nept);
extern "C" int TINKER_MOD(action, nebt);
extern "C" int TINKER_MOD(action, neat);
extern "C" int TINKER_MOD(action, nett);
extern "C" int TINKER_MOD(action, nev);
extern "C" int TINKER_MOD(action, ner);
extern "C" int TINKER_MOD(action, nedsp);
extern "C" int TINKER_MOD(action, nec);
extern "C" int TINKER_MOD(action, necd);
extern "C" int TINKER_MOD(action, ned);
extern "C" int TINKER_MOD(action, nem);
extern "C" int TINKER_MOD(action, nep);
extern "C" int TINKER_MOD(action, nect);
extern "C" int TINKER_MOD(action, new);
extern "C" int TINKER_MOD(action, nerxf);
extern "C" int TINKER_MOD(action, nes);
extern "C" int TINKER_MOD(action, nelf);
extern "C" int TINKER_MOD(action, neg);
extern "C" int TINKER_MOD(action, nex);

int& neb = TINKER_MOD(action, neb);
int& nea = TINKER_MOD(action, nea);
int& neba = TINKER_MOD(action, neba);
int& neub = TINKER_MOD(action, neub);
int& neaa = TINKER_MOD(action, neaa);
int& neopb = TINKER_MOD(action, neopb);
int& neopd = TINKER_MOD(action, neopd);
int& neid = TINKER_MOD(action, neid);
int& neit = TINKER_MOD(action, neit);
int& net = TINKER_MOD(action, net);
int& nept = TINKER_MOD(action, nept);
int& nebt = TINKER_MOD(action, nebt);
int& neat = TINKER_MOD(action, neat);
int& nett = TINKER_MOD(action, nett);
int& nev = TINKER_MOD(action, nev);
int& ner = TINKER_MOD(action, ner);
int& nedsp = TINKER_MOD(action, nedsp);
int& nec = TINKER_MOD(action, nec);
int& necd = TINKER_MOD(action, necd);
int& ned = TINKER_MOD(action, ned);
int& nem = TINKER_MOD(action, nem);
int& nep = TINKER_MOD(action, nep);
int& nect = TINKER_MOD(action, nect);
int& new_ = TINKER_MOD(action, new);
int& nerxf = TINKER_MOD(action, nerxf);
int& nes = TINKER_MOD(action, nes);
int& nelf = TINKER_MOD(action, nelf);
int& neg = TINKER_MOD(action, neg);
int& nex = TINKER_MOD(action, nex);
#endif
} }
