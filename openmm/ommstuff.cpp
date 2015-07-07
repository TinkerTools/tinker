/*
 *
 *    ############################################################
 *    ##                  COPYRIGHT (C) 2015                    ##
 *    ##     by Mark Friedrichs, Lee-Ping Wang & Jay Ponder     ##
 *    ##                  All Rights Reserved                   ##
 *    ############################################################
 *
 *    ############################################################
 *    ##                                                        ##
 *    ##  ommstuff.cpp  --  TINKER interface to the OpenMM API  ##
 *    ##                                                        ##
 *    ############################################################
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

// to convert from .c to .cpp, many things must be enclosed in extern "C" {}

extern "C" { 
    void kinetic_(double*, double(*)[3][3]); 
    void mdstat_(int*, double*, double*, double*, double*, double*, double*);
    void mdsave_(int*, double*, double*, double*);
    void ewca1_(double*);
    void born_();
    void bounds_();
    void empole1_();
    void egk1_();
    void enp1_(double*, double*);
}

typedef struct OpenMMData_s OpenMMData;

typedef struct {
    char s20[20];
} char20;

/*
 *    ############################################################
 *                    OpenMM Wrapper include Files
 *    ############################################################
 */

#include "OpenMMCWrapper.h"
#include "AmoebaOpenMMCWrapper.h"
#include "OpenMM.h"
using namespace OpenMM;

struct OpenMMData_s {
    OpenMM_System* system;
    OpenMM_Context* context;
    OpenMM_Integrator* integrator;
};

/*
 *    ############################################################
 *                   Map the TINKER Data Structures
 *    ############################################################
 */

int maxval__;

static struct {
    double* ak;
    double* anat;
    double* afld;
    int nangle;
    int* iang;
} angbnd__;
 
static struct {
    double angunit;
    double stbnunit;
    double aaunit;
    double opbunit;
    double opdunit;
    double cang;
    double qang;
    double pang;
    double sang;
    double copb;
    double qopb;
    double popb;
    double sopb;
    double copd;
    double qopd;
    double popd;
    double sopd;
    char* angtyp;
    char opbtyp[9];
} angpot__;

static struct {
    double* mass;
    int* tag;
    int* classs;   // variable "class" not allowed in C++; add an extra "s"
    int* atomic;
    int* valence;
    char* name; 
    char* story;
} atomid__;

static struct {
    double* x;
    double* y;
    double* z;
    int n;
    int* type;
} atoms__;

static struct {
    double kelvin;
    double atmsph;
    double tautemp;
    double taupres;
    double compress;
    double collide;
    double* vnh;
    double* qnh;
    double* gnh;
    double volmove;
    int voltrial;
    int isothermal;
    int isobaric;
    int anisotrop;
    char thermostat[12];
    char barostat[11];
    char volscale[10];
} bath__;

static struct {
    int* ibitor;
    int nbitor;
} bitor__;

static struct {
    double cbnd;
    double qbnd;
    double bndunit;
    char bndtyp[9];
} bndpot__;

static struct {
    double* bk;
    double* bl;
    int nbond;
    int* ibnd;
} bndstr__;

static struct {
    double* xbox;
    double* ybox;
    double* zbox;
    double alpha,beta,gamma;
    double* xbox2;
    double* ybox2;
    double* zbox2;
    double box34;
    double lvec[3][3];
    double recip[3][3];
    double volbox;
    double beta_sin,beta_cos;
    double gamma_sin,gamma_cos;
    double beta_term,gamma_term;
    int orthogonal,monoclinic;
    int triclinic,octahedron;
    char spacegrp[11];
} boxes__;

static struct {
    double electric;
    double dielec,ebuffer;
    double c2scale,c3scale;
    double c4scale,c5scale;
    int neutnbr,neutcut;
} chgpot__;

static struct {
    int* n12;
    int* i12;
    int* n13;
    int* i13;
    int* n14;
    int* i14;
    int* n15;
    int* i15;
    int maxval;
} couple__;

static struct {
    double* desum;
    double* deb;
    double* dea;
    double* deba;
    double* deub;
    double* deaa;
    double* deopb;
    double* deopd;
    double* deid;
    double* deit;
    double* det;
    double* dept;
    double* debt;
    double* deat;
    double* dett;
    double* dev;
    double* dec;
    double* decd;
    double* ded;
    double* dem;
    double* dep;
    double* der;
    double* des;
    double* delf;
    double* deg;
    double* dex;
} deriv__;

static struct {
    double* esum;
    double* eb;
    double* ea;
    double* eba;
    double* eub;
    double* eaa;
    double* eopb;
    double* eopd;
    double* eid;
    double* eit;
    double* et;
    double* ept;
    double* ebt;
    double* eat;
    double* ett;
    double* ev;
    double* ec;
    double* ecd;
    double* ed;
    double* em;
    double* ep;
    double* er;
    double* es;
    double* elf;
    double* eg;
    double* ex;
} energi__;

static struct {
    double aewald;
    char boundary[8];
} ewald__;

static struct {
    double* krat; 
    int nrat;
    int nratx;
    int* irat;
    int* iratx;
    int* kratx;
    int* ratimage;
    int use_rattle;
} freeze__;

static struct {
    int digits,iprint;
    int iwrite,isend;
    int verbose,debug;
    int holdup,abort;
} inform__;

static struct {
    double* ttx;
    double* tty;
    double* tbf;
    double* tbx;
    double* tby;
    double* tbxy;
    int* tnx;
    int* tny;
    char20* ktt;
    int maxntt;
    int maxtgrd;
} ktrtor__;

static struct {
    double* rad;
    double* eps;
    double* rad4;
    double* eps4;
    double* reduct;
} kvdws__;

static struct {
    double vdwcut,chgcut;
    double dplcut,mpolecut;
    double vdwtaper,chgtaper;
    double dpltaper,mpoletaper;
    double ewaldcut;
    int use_ewald,use_lights;
    int use_list,use_vlist;
    int use_clist,use_mlist;
} limits__; 

static struct {
    int nfree;
    int irest;
    int velsave;
    int frcsave;
    int uindsave;
    char integrate[11];
} mdstuf__;

static struct {
    double* v;
    double* a;
    double* aalt;
} moldyn__;

static struct {
    double m2scale,m3scale;
    double m4scale,m5scale;
} mplpot__;

static struct {
    double* pole;
    double* rpole;
    int npole;
    int* ipole;
    int* polsiz;
    int* pollist;
    char* polaxe;
    int* zaxis;
    int* xaxis;
    int* yaxis;
    int maxpole;
} mpole__;

static struct {
    double solvprs,surften;
    double spcut,spoff;
    double stcut,stoff;
    double* rcav;
    double* rdisp;
    double* cdisp;
} nonpol__;

static struct {
    double* opbk; 
    int* iopb;
    int nopbend;
} opbend__;

static struct {
    double* kpit;
    int* ipit;
    int npitors;
} pitors__;

static struct {
    double* bsmod1;
    double* bsmod2;
    double* bsmod3;
    double*** thetai1;
    double*** thetai2;
    double*** thetai3;
    double**** qgrid;
    double*** qfac;
    int nfft1,nfft2,nfft3;
    int bsorder;
    int* igrid;
} pme__;

static struct {
    double* polarity;
    double* thole;
    double* pdamp;
    double* uind;
    double* uinp;
    double* uinds;
    double* uinps;
    int npolar;
} polar__;

static struct {
    int* np11;
    int* ip11;
    int* np12;
    int* ip12;
    int* np13;
    int* ip13;
    int* np14;
    int* ip14;
    int maxp11;
    int maxp12;
    int maxp13;
    int maxp14;
} polgrp__;

static struct {
    double poleps,p2scale;
    double p3scale,p4scale;
    double p41scale,p5scale;
    double d1scale,d2scale;
    double d3scale,d4scale;
    double u1scale,u2scale;
    double u3scale,u4scale;
    char poltyp[7];
} polpot__;

static struct {
    int use_bond,use_angle,use_strbnd;
    int use_urey,use_angang,use_opbend;
    int use_opdist,use_improp,use_imptor;
    int use_tors,use_pitors,use_strtor;
    int use_angtor,use_tortor,use_vdw;
    int use_charge,use_chgdpl,use_dipole;
    int use_mpole,use_polar,use_rxnfld;
    int use_solv,use_metal,use_geom;
    int use_extra,use_born,use_orbit;
} potent__;

static struct {
    double* rsolv;
    double* asolv;
    double* rborn;
    double* drb;
    double* drbp;
    double* drobc;
    double doffset;
    double p1,p2,p3,p4,p5;
    double* gpol;
    double* shct;
    double* aobc;
    double* bobc;
    double* gobc;
    double* vsolv;
    double* wace;
    double* s2ace;
    double* uace;
    char solvtyp[9];
    char borntyp[9];
} solute__;

static struct {
    double friction;
    double* fgamma;
    int use_sdarea;
} stodyn__;

static struct {
    double* sbk;
    int* isb;
    int nstrbnd;
} strbnd__;

static struct {
    double idihunit,itorunit,torsunit;
    double ptorunit,storunit;
    double atorunit,ttorunit;
} torpot__;

static struct {
    double* tors1;
    double* tors2;
    double* tors3;
    double* tors4;
    double* tors5;
    double* tors6;
    int ntors;
    int* itors;
} tors__;

static struct {
    int* itt;
    int ntortor;
} tortor__;

static struct {
    double* uk;
    double* ul;
    int* iury;
    int nurey;
} urey__;

static struct {
    double cury;
    double qury;
    double ureyunit;
} urypot__;

static struct {
    int* nuse;
    int* iuse;
    int* use;
} usage__;

static struct {
    double* radmin;
    double* epsilon;
    double* radmin4;
    double* epsilon4;
    double* radhbnd;
    double* epshbnd;
    double* kred;
    int* ired;
    int nvdw;
    int* ivdw;
    int* jvdw;
    int maxclass;
} vdw__;

static struct {
    double abuck,bbuck,cbuck;
    double ghal,dhal;
    double v2scale,v3scale;
    double v4scale,v5scale;
    double* igauss;
    int ngauss;
    int use_vcorr;
    char vdwindex[6];
    char vdwtyp[14];
    char radtyp[6];
    char radsiz[9];
    char radrule[11];
    char epsrule[11];
    char gausstyp[9];
} vdwpot__;

static void setNullTerminator (char* string, int maxLength, char* buffer) {

    int count;
    int ptr_i;
    char* ptr_c;
    ptr_c = string;
    ptr_i = (int) (*ptr_c);
    count = 0;

    // add contents of string to buffer until a non-character
    // is encountered or the end of the string is reached,
    // then add NULL to the end of the buffer 

    while (ptr_i > 33 && ptr_i < 126 && count < maxLength) {
        buffer[count++] = (*ptr_c);
        ptr_c++;
        ptr_i = (int) (*ptr_c);
    } 

    buffer[count] = '\0';

    return;
}

extern "C" {
void set_parameters_ (int* maxval) {

   maxval__ = *maxval;
}

void set_angbnd_data_ (double* ak, double* anat, double* afld,
                       int* nangle, int* iang) {

   angbnd__.ak = ak;
   angbnd__.anat = anat;
   angbnd__.afld = afld;
   angbnd__.nangle = *nangle;
   angbnd__.iang = iang;
}

void set_angpot_data_ (double* angunit, double* stbnunit,
                       double* aaunit, double* opbunit, double* opdunit,
                       double* cang, double* qang, double* pang, double* sang,
                       double* copb, double* qopb, double* popb, double* sopb,
                       double* copd, double* qopd, double* popd, double* sopd,
                       char* angtyp, char* opbtyp) {

   angpot__.angunit = *angunit;
   angpot__.stbnunit = *stbnunit;
   angpot__.aaunit = *aaunit;
   angpot__.opbunit = *opbunit;
   angpot__.opdunit = *opdunit;
   angpot__.cang = *cang;
   angpot__.qang = *qang;
   angpot__.pang = *pang;
   angpot__.sang = *sang;
   angpot__.copb = *copb;
   angpot__.qopb = *qopb;
   angpot__.popb = *popb;
   angpot__.sopb = *sopb;
   angpot__.copd = *copd;
   angpot__.qopd = *qopd;
   angpot__.popd = *popd;
   angpot__.sopd = *sopd;
   angpot__.angtyp = angtyp;
   setNullTerminator (opbtyp, 8, angpot__.opbtyp);
}

void set_atomid_data_ (double* mass, int* tag, int* classs, int* atomic,
                       int* valence, char* name, char* story) {

   atomid__.mass = mass;
   atomid__.tag = tag;
   atomid__.classs = classs;
   atomid__.atomic = atomic;
   atomid__.valence = valence;
   atomid__.name = name;
   atomid__.story = story;
}

void set_atoms_data_ (double* x, double* y, double* z, int* n, int* type) {

   int ii;
   atoms__.x = x;
   atoms__.y = y;
   atoms__.z = z;
   atoms__.n = *n;
   atoms__.type = type;
}

void set_bath_data_ (double* kelvin, double* atmsph, double* tautemp,
                     double* taupres, double* compress, double* collide,
                     double* vnh, double* qnh, double* gnh, double* volmove,
                     int* voltrial, int* isothermal, int* isobaric,
                     int* anisotrop, char* thermostat, char* barostat,
                     char* volscale) {

   bath__.kelvin = *kelvin;
   bath__.atmsph = *atmsph;
   bath__.tautemp = *tautemp;
   bath__.taupres = *taupres;
   bath__.compress = *compress;
   bath__.collide = *collide;
   bath__.vnh = vnh;
   bath__.qnh = qnh;
   bath__.gnh = gnh;
   bath__.volmove = *volmove;
   bath__.voltrial = *voltrial;
   bath__.isothermal = *isothermal;
   bath__.isobaric = *isobaric;
   bath__.anisotrop = *anisotrop;
   setNullTerminator (thermostat, 11, bath__.thermostat);
   setNullTerminator (barostat, 10, bath__.barostat);
   setNullTerminator (volscale, 9, bath__.volscale);
}

void set_bitor_data_ (int* ibitor, int* nbitor) {

   bitor__.ibitor = ibitor;
   bitor__.nbitor = *nbitor;
}

void set_bndpot_data_ (double* cbnd, double* qbnd, double* bndunit,
                       char* bndtyp) {

   bndpot__.cbnd = *cbnd;
   bndpot__.qbnd = *qbnd;
   bndpot__.bndunit = *bndunit;
   setNullTerminator (bndtyp, 8, bndpot__.bndtyp);
}

void set_bndstr_data_ (double* bk, double* bl, int* nbond, int* ibnd) {

   bndstr__.bk = bk;
   bndstr__.bl = bl;
   bndstr__.nbond = *nbond;
   bndstr__.ibnd = ibnd;
}

void set_boxes_data_ (double* xbox, double* ybox, double* zbox,
                      double* alpha, double* beta, double* gamma,
                      double* xbox2, double* ybox2, double* zbox2, 
                      double* box34, double* lvec, double* recip, 
                      double* volbox, double* beta_sin, double* beta_cos, 
                      double* gamma_sin, double* gamma_cos,
                      double* beta_term, double* gamma_term, 
                      int* orthogonal, int* monoclinic, int* triclinic,
                      int* octahedron, char* spacegrp) {

   int ii, jj;
   double* ptr_d;

   boxes__.xbox = xbox;
   boxes__.ybox = ybox;
   boxes__.zbox = zbox;
   boxes__.alpha = *alpha;
   boxes__.beta = *beta;
   boxes__.gamma = *gamma;
   boxes__.xbox2 = xbox2;
   boxes__.ybox2 = ybox2;
   boxes__.zbox2 = zbox2;
   boxes__.box34 = *box34;

   ptr_d = lvec;
   for (ii = 0; ii < 3; ii++) {
       for (jj = 0; jj < 3; jj++) {
           boxes__.lvec[ii][jj] = *ptr_d;
           ptr_d++;
       }
   }
   ptr_d = recip;
   for (ii = 0; ii < 3; ii++) {
       for (jj = 0; jj < 3; jj++) {
           boxes__.recip[ii][jj] = *ptr_d;
           ptr_d++;
       }
   }

   boxes__.volbox = *volbox;
   boxes__.beta_sin = *beta_sin;
   boxes__.beta_cos = *beta_cos;
   boxes__.beta_term = *beta_term;
   boxes__.gamma_sin = *gamma_sin;
   boxes__.gamma_cos = *gamma_cos;
   boxes__.gamma_term = *gamma_term;
   boxes__.orthogonal = *orthogonal;
   boxes__.monoclinic = *monoclinic;
   boxes__.triclinic = *triclinic;
   boxes__.octahedron = *octahedron;

   setNullTerminator (spacegrp, 10, boxes__.spacegrp);
}

void set_chgpot_data_( double* electric,double* dielec, double* ebuffer,
                       double* c2scale, double* c3scale, double* c4scale,
                       double* c5scale, int* neutnbr, int* neutcut ) {

   chgpot__.electric = *electric;
   chgpot__.dielec = *dielec;
   chgpot__.ebuffer = *ebuffer;
   chgpot__.c2scale = *c2scale;
   chgpot__.c3scale = *c3scale;
   chgpot__.c4scale = *c4scale;
   chgpot__.c5scale = *c5scale;
   chgpot__.neutnbr = *neutnbr;
   chgpot__.neutcut = *neutcut;
}

void set_couple_data_( int* n12, int* i12, int* n13, int* i13, int* n14,
                       int* i14, int* n15, int* i15, int* maxval ) { 

   couple__.n12 = n12;
   couple__.i12 = i12;
   couple__.n13 = n13;
   couple__.i13 = i13;
   couple__.n14 = n14;
   couple__.i14 = i14;
   couple__.n15 = n15;
   couple__.i15 = i15;
   couple__.maxval = *maxval;
}

void set_deriv_data_( double* desum, double* deb, double* dea,
                      double* deba, double* deub, double* deaa,
                      double* deopb, double* deopd, double* deid,
                      double* deit, double* det, double* dept,
                      double* debt, double* deat, double* dett,
                      double* dev, double* dec, double* decd,
                      double* ded, double* dem, double* dep,
                      double* der, double* des, double* delf,
                      double* deg, double* dex ) {

   deriv__.desum = desum;
   deriv__.deb = deb;
   deriv__.dea = dea;
   deriv__.deba = deba;
   deriv__.deub = deub;
   deriv__.deaa = deaa;
   deriv__.deopb = deopb;
   deriv__.deopd = deopd;
   deriv__.deid = deid;
   deriv__.deit = deit;
   deriv__.det = det;
   deriv__.dept = dept;
   deriv__.debt = debt;
   deriv__.deat = deat;
   deriv__.dett = dett;
   deriv__.dev = dev;
   deriv__.dec = dec;
   deriv__.decd = decd;
   deriv__.ded = ded;
   deriv__.dem = dem;
   deriv__.dep = dep;
   deriv__.der = der;
   deriv__.des = des;
   deriv__.delf = delf;
   deriv__.deg = deg;
   deriv__.dex = dex;
}

void set_energi_data_( double* esum, double* eb, double* ea,
                       double* eba, double* eub, double* eaa,
                       double* eopb, double* eopd, double* eid,
                       double* eit, double* et, double* ept,
                       double* ebt, double* eat, double* ett,
                       double* ev, double* ec, double* ecd,
                       double* ed, double* em, double* ep,
                       double* er, double* es, double* elf,
                       double* eg, double* ex ) {

   energi__.esum = esum;
   energi__.eb = eb;
   energi__.ea = ea;
   energi__.eba = eba;
   energi__.eub = eub;
   energi__.eaa = eaa;
   energi__.eopb = eopb;
   energi__.eopd = eopd;
   energi__.eid = eid;
   energi__.eit = eit;
   energi__.et = et;
   energi__.ept = ept;
   energi__.ebt = ebt;
   energi__.eat = eat;
   energi__.ett = ett;
   energi__.ev = ev;
   energi__.ec = ec;
   energi__.ecd = ecd;
   energi__.ed = ed;
   energi__.em = em;
   energi__.ep = ep;
   energi__.er = er;
   energi__.es = es;
   energi__.elf = elf;
   energi__.eg = eg;
   energi__.ex = ex;
}

void set_ewald_data_( double* aewald, char* boundary ) {

   ewald__.aewald = *aewald;
   setNullTerminator( boundary, 7, ewald__.boundary );
}

void set_freeze_data_( double* krat, int* nrat, int* nratx, 
                       int* irat, int* iratx, int* kratx,
                       int* ratimage, int* use_rattle ) {

   freeze__.krat = krat;
   freeze__.nrat = *nrat;
   freeze__.nratx = *nratx;
   freeze__.irat = irat;
   freeze__.iratx = iratx;
   freeze__.kratx = kratx;
   freeze__.ratimage = ratimage;
   freeze__.use_rattle = *use_rattle;
}

void set_inform_data_( int* digits, int* iprint, int* iwrite, int* isend,
                       int* verbose, int* debug, int* holdup, int* abort ) {

   inform__.digits = *digits;
   inform__.iprint = *iprint;
   inform__.iwrite = *iwrite;
   inform__.isend = *isend;
   inform__.verbose = *verbose;
   inform__.debug = *debug;
   inform__.holdup = *holdup;
   inform__.abort = *abort;
}

void set_ktrtor_data_( double* ttx, double* tty, double* tbf, double* tbx,
                       double* tby, double* tbxy, int* tnx, int* tny,
                       char20* ktt, int* maxntt, int* maxtgrd ) {

   ktrtor__.ttx = ttx;
   ktrtor__.tty = tty;
   ktrtor__.tbf = tbf;
   ktrtor__.tbx = tbx;
   ktrtor__.tby = tby;
   ktrtor__.tbxy = tbxy;
   ktrtor__.tnx = tnx;
   ktrtor__.tny = tny;
   ktrtor__.ktt = ktt;
   ktrtor__.maxntt = *maxntt;
   ktrtor__.maxtgrd = *maxtgrd;
}

void set_kvdws_data_( double* rad, double* eps, double* rad4,
                      double* eps4, double* reduct ) {

   kvdws__.rad = rad;
   kvdws__.eps = eps;
   kvdws__.rad4 = rad4;
   kvdws__.eps4 = eps4;
   kvdws__.reduct = reduct;
}

void set_limits_data_( double* vdwcut, double* chgcut, double* dplcut,
                       double* mpolecut, double* vdwtaper, double* chgtaper,
                       double* dpltaper, double* mpoletaper, double* ewaldcut,
                       int* use_ewald, int* use_lights, int* use_list,
                       int* use_vlist, int* use_clist, int* use_mlist ) {

   limits__.vdwcut = *vdwcut;
   limits__.chgcut = *chgcut;
   limits__.dplcut = *dplcut;
   limits__.mpolecut = *mpolecut;
   limits__.vdwtaper = *vdwtaper;
   limits__.chgtaper = *chgtaper;
   limits__.dpltaper = *dpltaper;
   limits__.mpoletaper = *mpoletaper;
   limits__.ewaldcut = *ewaldcut;
   limits__.use_ewald = *use_ewald;
   limits__.use_lights = *use_lights;
   limits__.use_list = *use_list;
   limits__.use_vlist = *use_vlist;
   limits__.use_clist = *use_clist;
   limits__.use_mlist = *use_mlist;
}

void set_mdstuf_data_( int* nfree, int* irest, int* velsave,
                       int* frcsave, int* uindsave, char* integrate ) {

   mdstuf__.nfree = *nfree;
   mdstuf__.irest = *irest;
   mdstuf__.velsave = *velsave;
   mdstuf__.frcsave = *frcsave;
   mdstuf__.uindsave = *uindsave;

   setNullTerminator( integrate, 10, mdstuf__.integrate );
}

void set_moldyn_data_( double* v, double* a, double* aalt ) {

   moldyn__.v = v;
   moldyn__.a = a;
   moldyn__.aalt = aalt;
}

void set_mplpot_data_( double* m2scale, double* m3scale,
                       double* m4scale, double* m5scale ) {

   mplpot__.m2scale = *m2scale;
   mplpot__.m3scale = *m3scale;
   mplpot__.m4scale = *m4scale;
   mplpot__.m5scale = *m5scale;
}

void set_mpole_data_( double* pole, double* rpole, int* npole, 
                      int* ipole, int* polsiz, int* pollist, char* polaxe,
                      int* zaxis, int* xaxis, int* yaxis, int* maxpole ) {

   mpole__.pole = pole;
   mpole__.rpole = rpole;
   mpole__.npole = *npole;
   mpole__.ipole = ipole;
   mpole__.polsiz = polsiz;
   mpole__.pollist = pollist;
   mpole__.polaxe = polaxe;
   mpole__.zaxis = zaxis;
   mpole__.xaxis = xaxis;
   mpole__.yaxis = yaxis;
   mpole__.maxpole = *maxpole;
}

void set_nonpol_data_( double* solvprs, double* surften, double* spcut,
                       double* spoff, double* stcut, double* stoff,
                       double* rcav, double* rdisp, double* cdisp ) {

   nonpol__.solvprs = *solvprs;
   nonpol__.surften = *surften;
   nonpol__.spcut = *spcut;
   nonpol__.spoff = *spoff;
   nonpol__.stcut = *stcut;
   nonpol__.stoff = *stoff;
   nonpol__.rcav = rcav;
   nonpol__.rdisp = rdisp;
   nonpol__.cdisp = cdisp;
}

void set_opbend_data_( double* opbk, int* iopb, int* nopbend ) {

   opbend__.opbk = opbk;
   opbend__.iopb = iopb;
   opbend__.nopbend = *nopbend;
}

void set_pitors_data_( double* kpit, int* ipit, int* npitors ) {

   pitors__.kpit = kpit;
   pitors__.ipit = ipit;
   pitors__.npitors = *npitors;
}

void set_pme_data_( double* bsmod1, double* bsmod2, double* bsmod3,
                    double*** thetai1, double*** thetai2, double*** thetai3,
                    double**** qgrid, double*** qfac, int* nfft1,
                    int* nfft2, int* nfft3, int* bsorder, int* igrid ) {

   pme__.bsmod1 = bsmod1;
   pme__.bsmod2 = bsmod2;
   pme__.bsmod3 = bsmod3;
   pme__.thetai1 = thetai1;
   pme__.thetai2 = thetai2;
   pme__.thetai3 = thetai3;
   pme__.qfac = qfac;
   pme__.qgrid = qgrid;
   pme__.nfft1 = *nfft1;
   pme__.nfft2 = *nfft2;
   pme__.nfft3 = *nfft3;
   pme__.bsorder = *bsorder;
   pme__.igrid = igrid;
}

void set_polar_data_( double* polarity, double* thole, double* pdamp,
                      double* uind, double* uinp, double* uinds,
                      double* uinps, int* npolar ) {

   polar__.polarity = polarity;
   polar__.thole = thole;
   polar__.pdamp = pdamp;
   polar__.uind = uind;
   polar__.uinp = uinp;
   polar__.uinds = uinds;
   polar__.uinps = uinps;
   polar__.npolar = *npolar;
}

void set_polgrp_data_( int* np11, int* ip11, int* np12, int* ip12, 
                       int* np13, int* ip13, int* np14, int* ip14,
                       int* maxp11, int* maxp12, int* maxp13, int* maxp14 ) { 

   polgrp__.np11 = np11;
   polgrp__.ip11 = ip11;
   polgrp__.np12 = np12;
   polgrp__.ip12 = ip12;
   polgrp__.np13 = np13;
   polgrp__.ip13 = ip13;
   polgrp__.np14 = np14;
   polgrp__.ip14 = ip14;
   polgrp__.maxp11 = *maxp11;
   polgrp__.maxp12 = *maxp12;
   polgrp__.maxp13 = *maxp13;
   polgrp__.maxp14 = *maxp14;
}

void set_polpot_data_( double* poleps, double* p2scale, double* p3scale,
                       double* p4scale, double* p41scale, double* p5scale,
                       double* d1scale, double* d2scale, double* d3scale,
                       double* d4scale, double* u1scale, double* u2scale,
                       double* u3scale, double* u4scale, char* poltyp ) {

   polpot__.poleps = *poleps;
   polpot__.p2scale = *p2scale;
   polpot__.p3scale = *p3scale;
   polpot__.p4scale = *p4scale;
   polpot__.p41scale = *p41scale;
   polpot__.p5scale = *p5scale;
   polpot__.d1scale = *d1scale;
   polpot__.d2scale = *d2scale;
   polpot__.d3scale = *d3scale;
   polpot__.d4scale = *d4scale;
   polpot__.u1scale = *u1scale;
   polpot__.u2scale = *u2scale;
   polpot__.u3scale = *u3scale;
   polpot__.u4scale = *u4scale;
   setNullTerminator( poltyp, 7, polpot__.poltyp );
}

void set_potent_data_( int* use_bond, int* use_angle, int* use_strbnd,
                       int* use_urey, int* use_angang, int* use_opbend,
                       int* use_opdist, int* use_improp, int* use_imptor,
                       int* use_tors, int* use_pitors, int* use_strtor,
                       int* use_angtor, int* use_tortor, int* use_vdw,
                       int* use_charge, int* use_chgdpl, int* use_dipole,
                       int* use_mpole, int* use_polar, int* use_rxnfld,
                       int* use_solv, int* use_metal, int* use_geom,
                       int* use_extra, int* use_born, int* use_orbit ) {

   potent__.use_bond = *use_bond;
   potent__.use_angle = *use_angle;
   potent__.use_urey = *use_urey;
   potent__.use_strbnd = *use_strbnd;
   potent__.use_angang = *use_angang;
   potent__.use_opbend = *use_opbend;
   potent__.use_opdist = *use_opdist;
   potent__.use_improp = *use_improp;
   potent__.use_imptor = *use_improp;
   potent__.use_tors = *use_tors;
   potent__.use_pitors = *use_pitors;
   potent__.use_strtor = *use_strtor;
   potent__.use_angtor = *use_angtor;
   potent__.use_tortor = *use_tortor;
   potent__.use_vdw = *use_vdw;
   potent__.use_charge = *use_charge;
   potent__.use_chgdpl = *use_chgdpl;
   potent__.use_dipole = *use_dipole;
   potent__.use_mpole = *use_mpole;
   potent__.use_polar = *use_polar;
   potent__.use_rxnfld = *use_rxnfld;
   potent__.use_solv = *use_solv;
   potent__.use_metal = *use_metal;
   potent__.use_geom = *use_geom;
   potent__.use_extra = *use_extra;
   potent__.use_born = *use_born;
   potent__.use_orbit = *use_orbit;
}

void set_solute_data_( double* rsolv, double* asolv, double* rborn,
                       double* drb, double* drbp, double* drobc,
                       double* doffset, double* p1, double* p2, double* p3,
                       double* p4, double* p5, double* gpol, double* shct,
                       double* aobc, double* bobc, double* gobc,
                       double* vsolv, double* wace, double* s2ace,
                       double* uace, char* solvtyp, char* borntyp ) {

   solute__.rsolv = rsolv;
   solute__.asolv = asolv;
   solute__.rborn = rborn;
   solute__.drb = drb;
   solute__.drbp = drbp;
   solute__.drobc = drobc;
   solute__.doffset = *doffset;
   solute__.p1 = *p1;
   solute__.p2 = *p2;
   solute__.p3 = *p3;
   solute__.p4 = *p4;
   solute__.p5 = *p5;
   solute__.gpol = gpol;
   solute__.shct = shct;
   solute__.aobc = aobc;
   solute__.bobc = bobc;
   solute__.gobc = gobc;
   solute__.vsolv = vsolv;
   solute__.wace = wace;
   solute__.s2ace = s2ace;
   solute__.uace = uace;
   setNullTerminator( solvtyp, 8, solute__.solvtyp );
   setNullTerminator( borntyp, 8, solute__.borntyp);
}

void set_stodyn_data_( double* friction, double* fgamma, int* use_sdarea ) {

   stodyn__.friction = *friction;
   stodyn__.fgamma = fgamma;
   stodyn__.use_sdarea = *use_sdarea;
}

void set_strbnd_data_( double* sbk, int* isb, int* nstrbnd ) {

   strbnd__.sbk = sbk;
   strbnd__.isb = isb;
   strbnd__.nstrbnd = *nstrbnd;
}

void set_torpot_data_ (double* idihunit, double* itorunit, double* torsunit,
                       double* ptorunit, double* storunit, double* atorunit,
                       double* ttorunit) {

   torpot__.idihunit = *idihunit;
   torpot__.itorunit = *itorunit;
   torpot__.torsunit = *torsunit;
   torpot__.ptorunit = *ptorunit;
   torpot__.storunit = *storunit;
   torpot__.atorunit = *atorunit;
   torpot__.ttorunit = *ttorunit;
}

void set_tors_data_ (double* tors1, double* tors2, double* tors3,
                     double* tors4, double* tors5, double* tors6,
                     int* ntors, int* itors) {

   tors__.tors1 = tors1;
   tors__.tors2 = tors2;
   tors__.tors3 = tors3;
   tors__.tors4 = tors4;
   tors__.tors5 = tors5;
   tors__.tors6 = tors6;
   tors__.ntors = *ntors;
   tors__.itors = itors;
}

void set_tortor_data_ (int* itt, int* ntortor) {

   tortor__.itt = itt;
   tortor__.ntortor = *ntortor;
}

void set_urey_data_ (double* uk, double* ul, int* iury, int* nurey) {

   urey__.uk = uk;
   urey__.ul = ul;
   urey__.iury = iury;
   urey__.nurey = *nurey;
}

void set_urypot_data_ (double* cury, double* qury, double* ureyunit) {

   urypot__.cury = *cury;
   urypot__.qury = *qury;
   urypot__.ureyunit = *ureyunit;
}

void set_usage_data_ (int* nuse, int* iuse, int* use) {

   usage__.nuse = nuse;
   usage__.iuse = iuse;
   usage__.use = use;
}

void set_vdw_data_ (double* radmin, double* epsilon, double* radmin4,
                    double* epsilon4, double* radhbnd, double* epshbnd,
                    double* kred, int* ired, int* nvdw, int* ivdw, int* jvdw) {

   vdw__.radmin = radmin;
   vdw__.epsilon = epsilon;
   vdw__.radmin4 = radmin4;
   vdw__.epsilon4 = epsilon4;
   vdw__.radhbnd = radhbnd;
   vdw__.epshbnd = epshbnd;
   vdw__.kred = kred;
   vdw__.ired = ired;
   vdw__.nvdw = *nvdw;
   vdw__.ivdw = ivdw;
   vdw__.jvdw = jvdw;
}

void set_vdwpot_data_ (double* abuck, double* bbuck, double* cbuck,
                       double* ghal, double* dhal, double* v2scale,
                       double* v3scale, double* v4scale, double* v5scale,
                       double* igauss, int* ngauss, int* use_vcorr,
                       char* vdwindex, char* vdwtyp, char* radtyp,
                       char* radsiz, char* radrule, char* epsrule,
                       char* gausstyp) {

   vdwpot__.abuck = *abuck;
   vdwpot__.bbuck = *bbuck;
   vdwpot__.cbuck = *cbuck;
   vdwpot__.ghal = *ghal;
   vdwpot__.dhal = *dhal;
   vdwpot__.v2scale = *v2scale;
   vdwpot__.v3scale = *v3scale;
   vdwpot__.v4scale = *v4scale;
   vdwpot__.v5scale = *v5scale;
   vdwpot__.igauss = igauss;
   vdwpot__.ngauss = *ngauss;
   vdwpot__.use_vcorr = *use_vcorr;
   setNullTerminator (vdwindex, 5, vdwpot__.vdwindex);
   setNullTerminator (vdwtyp, 13, vdwpot__.vdwtyp);
   setNullTerminator (radtyp, 5, vdwpot__.radtyp);
   setNullTerminator (radsiz, 8, vdwpot__.radsiz);
   setNullTerminator (radrule, 10, vdwpot__.radrule);
   setNullTerminator (epsrule, 10, vdwpot__.epsrule);
   setNullTerminator (gausstyp, 8, vdwpot__.gausstyp);
}
}

/*
 *    ############################################################
 *             Setup Masses, COM Removal and Constraints
 *    ############################################################
 */

static void setupSystemParticles (OpenMM_System* system, FILE* log) {

    int ii;
    for (ii = 0; ii < atoms__.n; ii++) {
        OpenMM_System_addParticle (system, atomid__.mass[ii]);
    }
}

static void setupCMMotionRemover (OpenMM_System* system, FILE* log) {

    int frequency = mdstuf__.irest > 0 ? mdstuf__.irest : 100;
    OpenMM_CMMotionRemover* cMMotionRemover;
    cMMotionRemover = OpenMM_CMMotionRemover_create (frequency);
    OpenMM_System_addForce (system, (OpenMM_Force*) cMMotionRemover);

}

struct ConstraintMap {

    int* constraintOffset;  // offset into constraint list
    int* constraintCount;   // number of constraints for atom i
    int* constraintList;    // list of constraints sorted by atom index
     
    // For constraint involving atom i and j with i < j,
    //     constraintList[offset+kk] = j
    // where offset=constraintOffset[i], and 0 <= kk < constraintCount[i]
    // Note: one constraintOffset or constraintCount could be eliminated
    // since constraintCount[i] = constraintOffset[i+1] - constraintOffset[i]
};

static void freeConstraintMap (struct ConstraintMap* map) {

    free (map->constraintOffset);
    free (map->constraintCount);
    free (map->constraintList);
}

static void mapConstraints (struct ConstraintMap* map, FILE* log) {

    int ii, jj;
    int p1, p2;
    int offset, count;

    int numberOfParticles = atoms__.n;

    int* constraintCount = (int*) malloc (sizeof(int)*numberOfParticles);
    int* constraintOffset = (int*) malloc (sizeof(int)*numberOfParticles);

    memset (constraintCount, 0, sizeof(int)*numberOfParticles);
    memset (constraintOffset, 0, sizeof(int)*numberOfParticles);

    // count number of constraints for each particle where that particle
    // has the smaller index of the two constrainted particles

    for (ii = 0; ii < freeze__.nrat; ii++) {
        p1 = *(freeze__.irat+2*ii) - 1;
        p2 = *(freeze__.irat + 2*ii +1) - 1;
        if (p1 > p2) {
            p1 = p2;
        }
        constraintCount[p1]++;
    }

    // set the offset value

    constraintOffset[0] = 0;
    for (ii = 1; ii < numberOfParticles; ii++){
        constraintOffset[ii] = constraintOffset[ii-1] + constraintCount[ii-1];
    }
        
    // allocate constraint list and load

    int* constraintList = (int*)  malloc (sizeof(int)*freeze__.nrat);
    memset (constraintCount, 0, sizeof(int)*numberOfParticles);
    for (ii = 0; ii < freeze__.nrat; ii++) {
        p1 = *(freeze__.irat+2*ii) - 1;
        p2 = *(freeze__.irat + 2*ii +1) - 1;
        if (p1 > p2) {
            int p3 = p2;
            p2 = p1;
            p1 = p3;
        }
        offset = constraintOffset[p1];
        count = constraintCount[p1];
        constraintCount[p1]++;
        constraintList[offset+count] = p2;
    }

    if (log && 0) {
        for (ii = 0; ii < numberOfParticles; ii++) {
            offset = constraintOffset[ii];
            count = constraintCount[ii];
            (void) fprintf (stderr, "%5d Offset=%5d count=%5d: ",
                            ii, offset, count);
            for (jj = 0; jj < count; jj++) {
                 (void) fprintf (stderr, "%5d ", constraintList[offset+jj] );
            }
            (void) fprintf (stderr, "\n ");
        }
    }

    map->constraintCount = constraintCount;
    map->constraintOffset = constraintOffset;
    map->constraintList = constraintList;
}

static int checkForConstraint (struct ConstraintMap* map,
                               int p1, int p2, FILE* log) {

    int ii, jj;
    int offset;
    int match = 0;

    if (p1 > p2) {
        int p3 = p2;
        p2 = p1;
        p1 = p3;
    }

    offset = map->constraintOffset[p1];
    for (jj = 0; jj < map->constraintCount[p1] && match == 0; jj++) {
        if (map->constraintList[offset+jj] == p2) {
            match = 1;
        }
    }
    return match;
}

/*
 *    ############################################################
 *                  Setup Individual Potential Terms
 *    ############################################################
 */

static void setupAmoebaBondForce (OpenMM_System* system,
                                  int removeConstrainedBonds, FILE* log) {

    int ii, jj;
    int match;
    int* bondPtr;
    double kParameterConversion;

    struct ConstraintMap map;
    if (removeConstrainedBonds) {
        mapConstraints (&map, log);
    }

    OpenMM_AmoebaBondForce* amoebaBondForce;
    amoebaBondForce = OpenMM_AmoebaBondForce_create ();

    kParameterConversion = OpenMM_KJPerKcal
                              / (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom);
    bondPtr = bndstr__.ibnd;
    for (ii = 0; ii < bndstr__.nbond; ii++) {
        match = removeConstrainedBonds ? checkForConstraint (&map, (*bondPtr)-1,
                                                 *(bondPtr+1)-1, log ) : 0;
        if (match == 0) {
            OpenMM_AmoebaBondForce_addBond (amoebaBondForce, (*bondPtr)-1,
                      *(bondPtr+1)-1, bndstr__.bl[ii]*OpenMM_NmPerAngstrom,
                      kParameterConversion*bndpot__.bndunit*bndstr__.bk[ii]); 
        }
        bondPtr += 2;
    }
    OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic (amoebaBondForce,
                      bndpot__.cbnd/OpenMM_NmPerAngstrom);
    OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic (amoebaBondForce,
                  bndpot__.qbnd/(OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom));

    // if (OpenMM_AmoebaBondForce_getNumBonds (amoebaBondForce) > 0 ) {
         OpenMM_System_addForce (system, (OpenMM_Force*) amoebaBondForce);
    // }

    if (removeConstrainedBonds) {
        freeConstraintMap (&map);
    }
}

static void setupAmoebaAngleForce (OpenMM_System* system,
                                   int removeConstrainedAngles, FILE* log) {

    int ii, jj;
    int* angleIndexPtr;
    int match;
    char* angleTypPtr;

    struct ConstraintMap map;
    if (removeConstrainedAngles) {
        mapConstraints (&map, log);
    }

    OpenMM_AmoebaAngleForce* amoebaAngleForce;
    amoebaAngleForce = OpenMM_AmoebaAngleForce_create ();
    OpenMM_System_addForce (system, (OpenMM_Force*) amoebaAngleForce);

    // TINKER includes both harmonic and in-plane angles in these
    // data structs; they are separate in the OpenMM implementation

    angleIndexPtr = angbnd__.iang;
    angleTypPtr = angpot__.angtyp;
    for (ii = 0; ii < angbnd__.nangle; ii++) {
        if (strncasecmp( "HARMONIC", angleTypPtr, 8) == 0 ) { 
            match = removeConstrainedAngles ? checkForConstraint (&map,
                       *(angleIndexPtr)-1, *(angleIndexPtr+2)-1, log) : 0;
            if (match == 0) {
                OpenMM_AmoebaAngleForce_addAngle (amoebaAngleForce,
                       *(angleIndexPtr)-1, (*(angleIndexPtr+1))-1,
                       (*(angleIndexPtr+2))-1, angbnd__.anat[ii],
                       OpenMM_KJPerKcal*angpot__.angunit*angbnd__.ak[ii]); 
            }
        }
        angleIndexPtr += 4;
        angleTypPtr += 8;
    }
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic (amoebaAngleForce,
                                                       angpot__.cang);
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic (amoebaAngleForce,
                                                         angpot__.qang);
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic (amoebaAngleForce,
                                                        angpot__.pang);
    OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic (amoebaAngleForce,
                                                        angpot__.sang);

    if (removeConstrainedAngles) {
        freeConstraintMap (&map);
    }
}

static void setupAmoebaInPlaneAngleForce (OpenMM_System* system,
                                  int removeConstrainedBonds, FILE* log) {

    int ii, jj;
    int* angleIndexPtr;
    char* angleTypPtr;
    int match;
    double kParameterConversion;

    struct ConstraintMap map;
    if (removeConstrainedBonds) {
        mapConstraints (&map, log);
    }

    OpenMM_AmoebaInPlaneAngleForce* amoebaInPlaneAngleForce;
    amoebaInPlaneAngleForce = OpenMM_AmoebaInPlaneAngleForce_create ();
    OpenMM_System_addForce (system, (OpenMM_Force*) amoebaInPlaneAngleForce);

    // TINKER has both harmonic and in-plane angles in these structures;
    // they are separate in the OpenMM implementation

    angleIndexPtr = angbnd__.iang;
    angleTypPtr = angpot__.angtyp;
    for (ii = 0; ii < angbnd__.nangle; ii++) {
        if (strncasecmp( "IN-PLANE", angleTypPtr, 8 ) == 0) { 
            match = removeConstrainedBonds ? checkForConstraint (&map,
                        *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1, log) : 0;
            if (match == 0) {
              OpenMM_AmoebaInPlaneAngleForce_addAngle (amoebaInPlaneAngleForce,
                        *(angleIndexPtr)-1, (*(angleIndexPtr+1))-1,
                        (*(angleIndexPtr+2))-1, (*(angleIndexPtr+3))-1,
                        angbnd__.anat[ii],
                        OpenMM_KJPerKcal*angpot__.angunit*angbnd__.ak[ii]); 
            }
        }
        angleIndexPtr += 4;
        angleTypPtr += 8;
    }
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleCubic 
                        (amoebaInPlaneAngleForce, angpot__.cang);
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleQuartic
                        (amoebaInPlaneAngleForce, angpot__.qang);
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAnglePentic
                        (amoebaInPlaneAngleForce, angpot__.pang);
    OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleSextic
                        (amoebaInPlaneAngleForce, angpot__.sang);

    if (removeConstrainedBonds) {
        freeConstraintMap (&map);
    }
}

static void setupAmoebaStretchBendForce (OpenMM_System* system,
                                 int removeConstrainedBonds, FILE* log) {

    int ii, jj, index, abIndex, cbIndex;
    int* angleIndexPtr;
    double bondLengthAB;
    double bondLengthCB;
    double* bondLengthPtr;
    int match;

    struct ConstraintMap map;
    if (removeConstrainedBonds) {
        mapConstraints (&map, log);
    }

    OpenMM_AmoebaStretchBendForce* amoebaStretchBendForce;
    amoebaStretchBendForce = OpenMM_AmoebaStretchBendForce_create ();
    OpenMM_System_addForce (system, (OpenMM_Force*) amoebaStretchBendForce);

    for (ii = 0; ii < strbnd__.nstrbnd; ii++) {

        index = *(strbnd__.isb + 3*ii) - 1;

        abIndex = *(strbnd__.isb + 3*ii + 1);
        cbIndex = *(strbnd__.isb + 3*ii + 2);

        if (abIndex != 0) {
            bondLengthPtr = bndstr__.bl + abIndex - 1;
            bondLengthAB = (*bondLengthPtr)*OpenMM_NmPerAngstrom;
        } else {
            bondLengthAB = -1.0;
        }

        if (cbIndex != 0) {
            bondLengthPtr = bndstr__.bl + cbIndex - 1;
            bondLengthCB = (*bondLengthPtr)*OpenMM_NmPerAngstrom;
        } else {
            bondLengthCB = -1.0;
        }

        angleIndexPtr = angbnd__.iang + 4*index;

        match = removeConstrainedBonds ? checkForConstraint (&map,
                        *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1, log) : 0;
        if (match == 0) {
             OpenMM_AmoebaStretchBendForce_addStretchBend
                        (amoebaStretchBendForce, *(angleIndexPtr)-1,
                        (*(angleIndexPtr+1))-1, (*(angleIndexPtr+2))-1,
                        bondLengthAB, bondLengthCB,
                        OpenMM_RadiansPerDegree*(*(angbnd__.anat +index)),
                        (OpenMM_KJPerKcal/
            OpenMM_NmPerAngstrom)*angpot__.stbnunit*(*(strbnd__.sbk+2*ii)),
                        (OpenMM_KJPerKcal/
            OpenMM_NmPerAngstrom)*angpot__.stbnunit*(*(strbnd__.sbk+2*ii+1)) );
        }
    }

    if (removeConstrainedBonds) {
        freeConstraintMap (&map);
    }
}

static void setupAmoebaUreyBradleyForce (OpenMM_System* system,
                                 int removeConstrainedBonds, FILE* log) {

    int ii, jj;
    int* angleIndexPtr;
    double kParameterConversion;
    int match;

    struct ConstraintMap map;
    if (removeConstrainedBonds) {
        mapConstraints (&map, log);
    }

    OpenMM_HarmonicBondForce* harmonicBondForce;
    harmonicBondForce = OpenMM_HarmonicBondForce_create ();
    OpenMM_System_addForce (system, (OpenMM_Force*) harmonicBondForce);

    angleIndexPtr = urey__.iury;
    kParameterConversion = urypot__.ureyunit*OpenMM_KJPerKcal/
                               (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom)*2;
    for (ii = 0; ii < urey__.nurey; ii++) {
        match = removeConstrainedBonds ? checkForConstraint (&map,
                        *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1, log ) : 0;
        if (match == 0) {
            OpenMM_HarmonicBondForce_addBond (harmonicBondForce,
                                *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1,
                                OpenMM_NmPerAngstrom*urey__.ul[ii],
                                kParameterConversion*urey__.uk[ii]); 
        }
        angleIndexPtr += 3;
    }
    // OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic (amoebaBondForce,
    //                                                  urypot__.cury);
    // OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic (amoebaBondForce,
    //                                                    urypot__.qury);

    if (removeConstrainedBonds) {
        freeConstraintMap (&map);
    }
}

static void setupAmoebaOutOfPlaneBendForce (OpenMM_System* system, FILE* log) {

    int ii, index;
    int* angleIndexPtr;

    OpenMM_AmoebaOutOfPlaneBendForce* amoebaOutOfPlaneBendForce;
    amoebaOutOfPlaneBendForce = OpenMM_AmoebaOutOfPlaneBendForce_create ();
    OpenMM_System_addForce (system, (OpenMM_Force*) amoebaOutOfPlaneBendForce);

    for (ii = 0; ii < opbend__.nopbend; ii++) {
        index = *(opbend__.iopb + ii) - 1;
        angleIndexPtr = angbnd__.iang + 4*index;
        OpenMM_AmoebaOutOfPlaneBendForce_addOutOfPlaneBend
                    (amoebaOutOfPlaneBendForce, *(angleIndexPtr)-1,
                    (*(angleIndexPtr+1))-1, (*(angleIndexPtr+2))-1,
                    (*(angleIndexPtr+3))-1,
                    OpenMM_KJPerKcal*angpot__.opbunit*(*(opbend__.opbk +ii)));
    }

    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendCubic
                    (amoebaOutOfPlaneBendForce, angpot__.cang);
    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendQuartic
                    (amoebaOutOfPlaneBendForce, angpot__.qang);
    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendPentic
                    (amoebaOutOfPlaneBendForce, angpot__.pang);
    OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendSextic
                    (amoebaOutOfPlaneBendForce, angpot__.sang);
}

static void setupAmoebaTorsionForce (OpenMM_System* system, FILE* log) {

    int ii;
    double torsunit;
    int* torsIndexPtr;
    double* torsPtr;

    OpenMM_DoubleArray* torsion1;
    OpenMM_DoubleArray* torsion2;
    OpenMM_DoubleArray* torsion3;

    OpenMM_PeriodicTorsionForce* amoebaTorsionForce;

    torsion1 = OpenMM_DoubleArray_create(2);
    torsion2 = OpenMM_DoubleArray_create(2);
    torsion3 = OpenMM_DoubleArray_create(2);

    amoebaTorsionForce = OpenMM_PeriodicTorsionForce_create ();
    OpenMM_System_addForce (system, (OpenMM_Force*) amoebaTorsionForce);

    torsunit = OpenMM_KJPerKcal*torpot__.torsunit;
    torsIndexPtr = tors__.itors;
    for (ii = 0; ii < tors__.ntors; ii++) {

        torsPtr = tors__.tors1 + ii*4;
        OpenMM_DoubleArray_set (torsion1, 0, torsunit*(*torsPtr));
        OpenMM_DoubleArray_set (torsion1, 1, acos((*(torsPtr+2))));

        torsPtr = tors__.tors2 + ii*4;
        OpenMM_DoubleArray_set (torsion2, 0, torsunit*(*torsPtr));
        OpenMM_DoubleArray_set (torsion2, 1, acos((*(torsPtr+2))));

        torsPtr = tors__.tors3 + ii*4;
        OpenMM_DoubleArray_set (torsion3, 0, torsunit*(*torsPtr));
        OpenMM_DoubleArray_set (torsion3, 1, acos((*(torsPtr+2))));

        OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                  (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                  (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 1,
                  OpenMM_DoubleArray_get (torsion1,1),
                  OpenMM_DoubleArray_get (torsion1,0)); 

        OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                  (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                  (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 2,
                  OpenMM_DoubleArray_get(torsion2,1),
                  OpenMM_DoubleArray_get(torsion2,0)); 

        OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                  (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                  (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 3,
                  OpenMM_DoubleArray_get(torsion3,1),
                  OpenMM_DoubleArray_get(torsion3,0)); 

        torsIndexPtr += 4;
    }

    OpenMM_DoubleArray_destroy (torsion1);
    OpenMM_DoubleArray_destroy (torsion2);
    OpenMM_DoubleArray_destroy (torsion3);
}

static void setupAmoebaPiTorsionForce( OpenMM_System* system, FILE* log ) {

    int ii;
    int* piTorsIndexPtr;

    OpenMM_AmoebaPiTorsionForce* amoebaPiTorsionForce;
    amoebaPiTorsionForce = OpenMM_AmoebaPiTorsionForce_create ();
    OpenMM_System_addForce (system, (OpenMM_Force*) amoebaPiTorsionForce);

    piTorsIndexPtr = pitors__.ipit;
    for (ii = 0; ii < pitors__.npitors; ii++) {
        OpenMM_AmoebaPiTorsionForce_addPiTorsion (amoebaPiTorsionForce,
                      (*piTorsIndexPtr) -1, (*(piTorsIndexPtr+1))-1,
                      (*(piTorsIndexPtr+2))-1, (*(piTorsIndexPtr+3))-1,
                      (*(piTorsIndexPtr+4))-1, (*(piTorsIndexPtr+5))-1,
                      OpenMM_KJPerKcal*torpot__.ptorunit*pitors__.kpit[ii]);
        piTorsIndexPtr += 6;
    }
}

static int getChiralIndex( int atomB, int atomC, int atomD ){
 
    int ii, j, m, k;
    int chiralAtom;

    // test for chirality at the central torsion-torsion site
   
    chiralAtom = -1;
    if (*(couple__.n12 + atomC) == 4) {
        j = 0;
        for (ii = 0; ii < 4; ii++) {
            m = *(couple__.i12 + atomC*maxval__ + ii) - 1;
            if (m != atomB && m != atomD) {
                if (j == 0) {
                  j = m;
               } else {
                  k = m;
               }
            }
        }
        if (atoms__.type[j] > atoms__.type[k])  chiralAtom = j;
        if (atoms__.type[k] > atoms__.type[j])  chiralAtom = k;
        if (atomid__.atomic[j] > atomid__.atomic[k])  chiralAtom = j;
        if (atomid__.atomic[k] > atomid__.atomic[j])  chiralAtom = k;
    }
    return chiralAtom;
}

static void setupAmoebaTorsionTorsionForce (OpenMM_System* system, FILE* log) {

    int ii, jj, kk, index, count;
    int ia, ib, ic, id, ie, ichiral;
    int gridIndex;
    int xIndex, yIndex;
    int* ibitorPtr;
    int* ittPtr;
    int maxntt, maxtgrd, maxtgrd2;
    int addIndex;
    int numberOfTorsionTorsions;
    OpenMM_DoubleArray* values;
    OpenMM_3D_DoubleArray* grid;

    OpenMM_AmoebaTorsionTorsionForce* amoebaTorsionTorsionForce;
    amoebaTorsionTorsionForce = OpenMM_AmoebaTorsionTorsionForce_create ();
    OpenMM_System_addForce (system, (OpenMM_Force*) amoebaTorsionTorsionForce);

    // atoms/grid in torsion-torsion

    ittPtr = tortor__.itt;
    for (ii = 0; ii < tortor__.ntortor; ii++) {
        index = *(ittPtr) - 1;
        gridIndex = *(ittPtr+1) - 1;
        ibitorPtr = bitor__.ibitor + index*5;
        if (*(ittPtr+2) == 1) {
            count = 0;
            ia = *(ibitorPtr + count++) - 1;
            ib = *(ibitorPtr + count++) - 1;
            ic = *(ibitorPtr + count++) - 1;
            id = *(ibitorPtr + count++) - 1;
            ie = *(ibitorPtr + count++) - 1;
        } else {
            count = 4;
            ia = *(ibitorPtr + count--) - 1;
            ib = *(ibitorPtr + count--) - 1;
            ic = *(ibitorPtr + count--) - 1;
            id = *(ibitorPtr + count--) - 1;
            ie = *(ibitorPtr + count--) - 1;
        }
        ichiral = getChiralIndex (ib, ic, id);
        OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion
                        (amoebaTorsionTorsionForce, ia, ib, ic, id, ie,
                         ichiral, gridIndex );
        ittPtr += 3;
    }

    numberOfTorsionTorsions  =
                     OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions
                         (amoebaTorsionTorsionForce);

    // grids

    maxntt = ktrtor__.maxntt;
    maxtgrd = ktrtor__.maxtgrd;
    maxtgrd2 = maxtgrd*maxtgrd;
    values = OpenMM_DoubleArray_create(6);

    if (numberOfTorsionTorsions) {

        int kk = 0;
        for (int i = 0; i < maxntt; i++) {
            char char21[21];
            for (int j = 0; j < 20; j++) {
                char21[j] = ktrtor__.ktt[i].s20[j];
            }
            char21[20] = '\n';
            if (char21[0] != ' ') {
                kk++;
            }
        }

        for (ii = 0; ii < kk; ii++) {
            grid = OpenMM_3D_DoubleArray_create (*(ktrtor__.tnx+ii),
                                                 *(ktrtor__.tny+ii), 6);
            xIndex = 0;
            yIndex = 0;
    
            for (jj = 0; jj < *(ktrtor__.tnx+ii)*(*(ktrtor__.tny+ii)); jj++) {
                addIndex  = 0;
                OpenMM_DoubleArray_set (values, addIndex++,
                    *(ktrtor__.ttx + maxtgrd*ii+xIndex));
                OpenMM_DoubleArray_set (values, addIndex++,
                    *(ktrtor__.tty + maxtgrd*ii+yIndex));
    
                OpenMM_DoubleArray_set (values, addIndex++,
                    OpenMM_KJPerKcal*(*(ktrtor__.tbf  + maxtgrd2*ii + jj)));
                OpenMM_DoubleArray_set (values, addIndex++,
                    OpenMM_KJPerKcal*(*(ktrtor__.tbx  + maxtgrd2*ii + jj)));
                OpenMM_DoubleArray_set (values, addIndex++,
                    OpenMM_KJPerKcal*(*(ktrtor__.tby  + maxtgrd2*ii + jj)));
                OpenMM_DoubleArray_set (values, addIndex++,
                    OpenMM_KJPerKcal*(*(ktrtor__.tbxy + maxtgrd2*ii + jj)));

                OpenMM_3D_DoubleArray_set (grid, yIndex, xIndex, values);

                xIndex++;
                if (xIndex == *(ktrtor__.tnx+ii)) {
                    xIndex = 0;
                    yIndex++;
                }
            }

            OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid
                    (amoebaTorsionTorsionForce, ii, grid);
            OpenMM_3D_DoubleArray_destroy (grid);
        }
    }

    OpenMM_DoubleArray_destroy (values);
}

static void setupAmoebaWcaDispersionForce (OpenMM_System* system, FILE* log) {

    int ii;
    double cdispTotal = 0.0;
    double epso = 0.1100;
    double epsh = 0.0135;
    double rmino = 1.7025;
    double rminh = 1.3275;
    double awater = 0.033428;
    double slevy = 1.0;
//  double dispoff = 0.26;
    double dispoff = 0.0;
    double shctd = 0.81;

    OpenMM_AmoebaWcaDispersionForce*  amoebaWcaDispersionForce;
    amoebaWcaDispersionForce = OpenMM_AmoebaWcaDispersionForce_create ();
    OpenMM_System_addForce (system, (OpenMM_Force*) amoebaWcaDispersionForce);

    for (ii = 0; ii < atoms__.n; ii++) {
        cdispTotal += nonpol__.cdisp[ii];
        OpenMM_AmoebaWcaDispersionForce_addParticle (amoebaWcaDispersionForce,
                OpenMM_NmPerAngstrom*kvdws__.rad[atomid__.classs[ii]-1],
                OpenMM_KJPerKcal*kvdws__.eps[atomid__.classs[ii]-1]);
    }

    OpenMM_AmoebaWcaDispersionForce_setEpso (amoebaWcaDispersionForce,
                                             epso*OpenMM_KJPerKcal);
    OpenMM_AmoebaWcaDispersionForce_setEpsh (amoebaWcaDispersionForce,
                                             epsh*OpenMM_KJPerKcal);
    OpenMM_AmoebaWcaDispersionForce_setRmino (amoebaWcaDispersionForce,
                                              rmino*OpenMM_NmPerAngstrom);
    OpenMM_AmoebaWcaDispersionForce_setRminh (amoebaWcaDispersionForce,
                                              rminh*OpenMM_NmPerAngstrom);
    OpenMM_AmoebaWcaDispersionForce_setDispoff (amoebaWcaDispersionForce,
                                                dispoff*OpenMM_NmPerAngstrom);
    OpenMM_AmoebaWcaDispersionForce_setAwater (amoebaWcaDispersionForce,
       awater/(OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom));
    OpenMM_AmoebaWcaDispersionForce_setSlevy (amoebaWcaDispersionForce, slevy);
    OpenMM_AmoebaWcaDispersionForce_setShctd (amoebaWcaDispersionForce, shctd);
}

static void computePeriodicBoxVectors(OpenMM_Vec3 &a, OpenMM_Vec3 &b, OpenMM_Vec3 &c) {
  /* a, b, and c are the output vectors, computed from the unit cell dimensions
   * in boxes__. This should not be called when no box is defined
   */
  const double DEG_TO_RAD = M_PI / 180.0;

  a.x = OpenMM_NmPerAngstrom*(*boxes__.xbox);
  a.y = a.z = 0.0;

  b.x = OpenMM_NmPerAngstrom*(*boxes__.ybox)*cos(DEG_TO_RAD*boxes__.gamma);
  b.y = OpenMM_NmPerAngstrom*(*boxes__.ybox)*sin(DEG_TO_RAD*boxes__.gamma);
  b.z = 0.0;

  c.x = (*boxes__.zbox)*cos(DEG_TO_RAD*boxes__.beta);
  c.y = *boxes__.zbox*(cos(DEG_TO_RAD*boxes__.alpha) -
                cos(DEG_TO_RAD*boxes__.beta)*cos(DEG_TO_RAD*boxes__.gamma)) /
                sin(DEG_TO_RAD*boxes__.gamma);
  c.z = sqrt((*boxes__.zbox)*(*boxes__.zbox)-c.x*c.x-c.y*c.y);

  c.x *= OpenMM_NmPerAngstrom;
  c.y *= OpenMM_NmPerAngstrom;
  c.z *= OpenMM_NmPerAngstrom;

  // Make numbers that should be zero exactly zero
  if (abs(a.x) < 1e-6) a.x = 0.0;
  if (abs(a.y) < 1e-6) a.y = 0.0;
  if (abs(a.z) < 1e-6) a.z = 0.0;
  if (abs(b.x) < 1e-6) b.x = 0.0;
  if (abs(b.y) < 1e-6) b.y = 0.0;
  if (abs(b.z) < 1e-6) b.z = 0.0;
  if (abs(c.x) < 1e-6) c.x = 0.0;
  if (abs(c.y) < 1e-6) c.y = 0.0;
  if (abs(c.z) < 1e-6) c.z = 0.0;

  // Reduce the box vectors if necessary
  c.x -= b.x*round(c.y/b.y) - a.x*round(c.x/a.x);;
  c.y -= b.y*round(c.y/b.y);

  b.x -= a.x*round(b.x/a.x);

#if 0
  // DEBUG
  fprintf(stderr, "Box vectors:\n");
  fprintf(stderr, "[ %10.4f  %10.4f  %10.4f ]\n", a.x, a.y, a.z);
  fprintf(stderr, "[ %10.4f  %10.4f  %10.4f ]\n", b.x, b.y, b.z);
  fprintf(stderr, "[ %10.4f  %10.4f  %10.4f ]\n", c.x, c.y, c.z);
#endif /* 0 */
}

static void setDefaultPeriodicBoxVectors (OpenMM_System* system, FILE* log) {

  OpenMM_Vec3 a;
  OpenMM_Vec3 b;
  OpenMM_Vec3 c;

  // Nothing to do if we have no box
  if (*boxes__.xbox == 0) return;

  computePeriodicBoxVectors(a, b, c);
  OpenMM_System_setDefaultPeriodicBoxVectors (system, &a, &b, &c);
}

static void printDefaultPeriodicBoxVectors (OpenMM_System* system, FILE* log) {

    OpenMM_Vec3 a;
    OpenMM_Vec3 b;
    OpenMM_Vec3 c;
   
    OpenMM_System_getDefaultPeriodicBoxVectors (system, &a, &b, &c);
    a.x = a.x / OpenMM_NmPerAngstrom;
    a.y = a.y / OpenMM_NmPerAngstrom;
    a.z = a.z / OpenMM_NmPerAngstrom;
    b.x = b.x / OpenMM_NmPerAngstrom;
    b.y = b.y / OpenMM_NmPerAngstrom;
    b.z = b.z / OpenMM_NmPerAngstrom;
    c.x = c.x / OpenMM_NmPerAngstrom;
    c.y = c.y / OpenMM_NmPerAngstrom;
    c.z = c.z / OpenMM_NmPerAngstrom;
    (void) fprintf (log, "\n Box Size: [%12.4f %12.4f %12.4f ]",
                    a.x, a.y, a.z );
    (void) fprintf (log, "\n  (Ang)    [%12.4f %12.4f %12.4f ]",
                    b.x, b.y, b.z );
    (void) fprintf (log, "\n           [%12.4f %12.4f %12.4f ]\n",
                    c.x, c.y, c.z );
}

static void setupAmoebaVdwForce (OpenMM_System* system, FILE* log) {

    char buffer[128];
    int ii,jj,i;
    OpenMM_Boolean useCorrection;
    OpenMM_IntArray* exclusions;

    OpenMM_AmoebaVdwForce* amoebaVdwForce;
    amoebaVdwForce = OpenMM_AmoebaVdwForce_create ();
    OpenMM_System_addForce (system, (OpenMM_Force*) amoebaVdwForce);

    for (ii = 0; ii < atoms__.n; ii++) {
        i = vdw__.jvdw[ii];
        OpenMM_AmoebaVdwForce_addParticle (amoebaVdwForce, vdw__.ired[ii]-1,
                             OpenMM_NmPerAngstrom*(kvdws__.rad[i-1]),
                             OpenMM_KJPerKcal*(kvdws__.eps[i-1]),
                             vdw__.kred[ii]);
    }

    setNullTerminator (vdwpot__.radrule, 10, buffer);
    OpenMM_AmoebaVdwForce_setSigmaCombiningRule (amoebaVdwForce, buffer);

    setNullTerminator (vdwpot__.epsrule, 10, buffer);
    OpenMM_AmoebaVdwForce_setEpsilonCombiningRule (amoebaVdwForce, buffer);

    OpenMM_AmoebaVdwForce_setCutoff (amoebaVdwForce,
                                     limits__.vdwcut*OpenMM_NmPerAngstrom);

    useCorrection = OpenMM_False;
    if (vdwpot__.use_vcorr)  useCorrection = OpenMM_True;
    OpenMM_AmoebaVdwForce_setUseDispersionCorrection (amoebaVdwForce,
                                                      useCorrection);

    if (*boxes__.xbox == 0)
        OpenMM_AmoebaVdwForce_setNonbondedMethod (amoebaVdwForce,
                                    OpenMM_AmoebaVdwForce_NoCutoff);
    else
        OpenMM_AmoebaVdwForce_setNonbondedMethod (amoebaVdwForce,
                                    OpenMM_AmoebaVdwForce_CutoffPeriodic);

    setDefaultPeriodicBoxVectors (system, log);

    exclusions = OpenMM_IntArray_create (0);
    for (ii = 0; ii < atoms__.n; ii++) {
        OpenMM_IntArray_append (exclusions, ii);
        if (fabs( vdwpot__.v2scale ) <= 0.0) {
            for (jj = 0; jj < *(couple__.n12 + ii); jj++) {
                OpenMM_IntArray_append (exclusions,
                        *(couple__.i12 + couple__.maxval*ii + jj)-1 );
            }
         }
         if (fabs(vdwpot__.v3scale) <= 0.0) {
            for (jj = 0; jj < *(couple__.n13 + ii); jj++) {
                OpenMM_IntArray_append (exclusions,
                        *(couple__.i13 + 3*couple__.maxval*ii + jj)-1 );
            }
        }
        OpenMM_AmoebaVdwForce_setParticleExclusions (amoebaVdwForce,
                                                     ii, exclusions);
        OpenMM_IntArray_resize (exclusions, 0);
    }
    OpenMM_IntArray_destroy (exclusions);
}

static void setupAndersenThermostat (OpenMM_System* system, FILE* log) {

    OpenMM_AndersenThermostat* andersenThermostat;
    
    andersenThermostat = OpenMM_AndersenThermostat_create (bath__.kelvin,
                                                           bath__.collide);
    OpenMM_System_addForce (system, (OpenMM_Force*) andersenThermostat);

    if (log) {
        (void) fprintf (log, "\n Andersen Thermostat:\n" );
        (void) fprintf (log, "\n Temperature          %15.4f K",
                        OpenMM_AndersenThermostat_getDefaultTemperature
                        (andersenThermostat));
        (void) fprintf (log, "\n Collision Frequency  %15.7e ps^(-1)",
                        OpenMM_AndersenThermostat_getDefaultCollisionFrequency
                        (andersenThermostat));
        (void) fprintf (log, "\n Random Number Seed     %d\n",
                        OpenMM_AndersenThermostat_getRandomNumberSeed
                        (andersenThermostat));
    }
}

static void setupMonteCarloBarostat (OpenMM_System* system, FILE* log) {

    OpenMM_MonteCarloBarostat* monteCarloBarostat;

    int frequency = 25;
    monteCarloBarostat = OpenMM_MonteCarloBarostat_create
                             (bath__.atmsph*1.01295, bath__.kelvin, frequency);
    OpenMM_System_addForce (system, (OpenMM_Force*) monteCarloBarostat);

    if (log) {
        (void) fprintf (log, "\n MonteCarlo Barostat:\n");
        (void) fprintf (log, "\n Temperature          %15.4f K",
                        OpenMM_MonteCarloBarostat_getTemperature
                        (monteCarloBarostat));
        (void) fprintf (log, "\n Pressure             %15.4f atm",
                        OpenMM_MonteCarloBarostat_getDefaultPressure
                        (monteCarloBarostat));
        (void) fprintf (log, "\n Frequency              %d",
                        OpenMM_MonteCarloBarostat_getFrequency
                        (monteCarloBarostat));
        (void) fprintf (log, "\n Random Number Seed     %d\n",
                        OpenMM_MonteCarloBarostat_getRandomNumberSeed
                        (monteCarloBarostat));
    }
}

static void loadCovalentArray (int numberToLoad, int* valuesToLoad,
                               OpenMM_IntArray* covaletMap) {

    int ii;
    OpenMM_IntArray_resize (covaletMap, numberToLoad);
    for (ii = 0; ii < numberToLoad; ii++) {
        OpenMM_IntArray_set (covaletMap, ii, *(valuesToLoad +ii) - 1);
    }
}

static void setupAmoebaMultipoleForce (OpenMM_System* system, FILE* log) {

    char buffer[128];
    char* axisPtr;
    int ii, jj, index;
    int invalidAxis, errorReport;
    double* polePtr;

    OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes axisType;
    OpenMM_DoubleArray* dipoles;
    OpenMM_DoubleArray* quadrupoles;

    OpenMM_IntArray* covalentMap;
    OpenMM_IntArray* gridDimensions;

    double dipoleConversion = OpenMM_NmPerAngstrom;
    double quadrupoleConversion = OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom;
    double polarityConversion = OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom
                                         *OpenMM_NmPerAngstrom;
    double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

    OpenMM_AmoebaMultipoleForce* amoebaMultipoleForce;
    amoebaMultipoleForce = OpenMM_AmoebaMultipoleForce_create();
    OpenMM_System_addForce(system, (OpenMM_Force*) amoebaMultipoleForce);

    dipoles = OpenMM_DoubleArray_create(3);
    quadrupoles = OpenMM_DoubleArray_create(9);
    errorReport = 0;
    for (ii = 0; ii < atoms__.n; ii++) {

        axisPtr = mpole__.polaxe + ii*8;
        if (strncasecmp (axisPtr, "Z-then-X", 8) == 0) {
            axisType = OpenMM_AmoebaMultipoleForce_ZThenX;
        } else if (strncasecmp(axisPtr, "Bisector", 8) == 0) { 
            axisType = OpenMM_AmoebaMultipoleForce_Bisector;
        } else if (strncasecmp(axisPtr, "Z-Bisect", 8) == 0) { 
            axisType = OpenMM_AmoebaMultipoleForce_ZBisect;
        } else if (strncasecmp(axisPtr, "3-Fold", 6) == 0) { 
            axisType = OpenMM_AmoebaMultipoleForce_ThreeFold;
        } else if (strncasecmp(axisPtr, "Z-Only", 6) == 0) { 
            axisType = OpenMM_AmoebaMultipoleForce_ZOnly;
        } else if (strncasecmp(axisPtr, "None", 4) == 0
                      || strncasecmp( axisPtr, "    ",   4 ) == 0 ) { 
            axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
        } else {
            errorReport++;
            setNullTerminator (axisPtr, 8, buffer);
            (void) fprintf (stderr,
      "setupAmoebaMultipoleForce: Axis Type=%s for Atom %7d Not Supported\n",
                   buffer, ii );
            if (errorReport > 20) {
                (void) fflush (stderr);
                exit (-1);
            }
        }

        polePtr = mpole__.pole + ii*mpole__.maxpole + 1;
        for (jj = 0; jj < 3; jj++) {
            OpenMM_DoubleArray_set (dipoles, jj,
                                      (*(polePtr))*dipoleConversion);
            polePtr++;
        }
        for (jj = 0; jj < 9; jj++) {
            OpenMM_DoubleArray_set (quadrupoles, jj,
                                      (*(polePtr))*quadrupoleConversion);
            polePtr++;
        }

        polePtr = mpole__.pole + ii*mpole__.maxpole;
        OpenMM_AmoebaMultipoleForce_addMultipole (amoebaMultipoleForce,
                                     *polePtr, dipoles, quadrupoles, axisType,
                                     *(mpole__.zaxis+ii)-1,
                                     *(mpole__.xaxis+ii)-1,
                                     *(mpole__.yaxis+ii)-1,
                                     polar__.thole[ii],
                                     polar__.pdamp[ii]*dampingFactorConversion,
                                     polar__.polarity[ii]*polarityConversion);
    }

    if (errorReport) {
        exit (-1);
    }

    if (limits__.use_ewald) {

        double ewaldTolerance = 1.0e-04;
        OpenMM_AmoebaMultipoleForce_setNonbondedMethod (amoebaMultipoleForce,
                                     OpenMM_AmoebaMultipoleForce_PME);
        OpenMM_AmoebaMultipoleForce_setCutoffDistance (amoebaMultipoleForce,
                                     limits__.ewaldcut*OpenMM_NmPerAngstrom);
        OpenMM_AmoebaMultipoleForce_setAEwald (amoebaMultipoleForce,
                                     ewald__.aewald/OpenMM_NmPerAngstrom);
        // OpenMM_AmoebaMultipoleForce_setPmeBSplineOrder (amoebaMultipoleForce,
        //                           pme__.bsorder);
    
        // grid dimensions
    
        gridDimensions = OpenMM_IntArray_create (3);
    
        OpenMM_IntArray_set (gridDimensions, 0, pme__.nfft1);
        OpenMM_IntArray_set (gridDimensions, 1, pme__.nfft2);
        OpenMM_IntArray_set (gridDimensions, 2, pme__.nfft3);
    
        OpenMM_AmoebaMultipoleForce_setPmeGridDimensions (amoebaMultipoleForce,
                                                          gridDimensions);
        OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance
                                     (amoebaMultipoleForce, ewaldTolerance);
        OpenMM_IntArray_destroy (gridDimensions);
        setDefaultPeriodicBoxVectors (system, log);
    
    } else {
        OpenMM_AmoebaMultipoleForce_setNonbondedMethod (amoebaMultipoleForce,
                                     OpenMM_AmoebaMultipoleForce_NoCutoff);
    }

    if (strncasecmp (polpot__.poltyp, "DIRECT", 6) == 0) { 
        OpenMM_AmoebaMultipoleForce_setPolarizationType (amoebaMultipoleForce,
                                     OpenMM_AmoebaMultipoleForce_Direct);
    } else {
        OpenMM_AmoebaMultipoleForce_setPolarizationType (amoebaMultipoleForce,
                                     OpenMM_AmoebaMultipoleForce_Mutual);
    }

    int PolType_out = OpenMM_AmoebaMultipoleForce_getPolarizationType
                                     (amoebaMultipoleForce);

    OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations
                                     (amoebaMultipoleForce, 500);
    OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon
                                     (amoebaMultipoleForce, polpot__.poleps);

    covalentMap = OpenMM_IntArray_create (0);
    for (ii = 0; ii < atoms__.n; ii++) {

        loadCovalentArray (*(couple__.n12 + ii),
                           (couple__.i12 + couple__.maxval*ii), covalentMap);  
        OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                           OpenMM_AmoebaMultipoleForce_Covalent12,
                           covalentMap);

        loadCovalentArray (*(couple__.n13 + ii),
                           (couple__.i13 + 3*couple__.maxval*ii), covalentMap);  
        OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                           OpenMM_AmoebaMultipoleForce_Covalent13,
                           covalentMap);

        loadCovalentArray (*(couple__.n14 + ii),
                           (couple__.i14 + 9*couple__.maxval*ii), covalentMap);  
        OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                           OpenMM_AmoebaMultipoleForce_Covalent14,
                           covalentMap);

        loadCovalentArray (*(couple__.n15 + ii),
                           (couple__.i15 + 27*couple__.maxval*ii), covalentMap);  
        OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                           OpenMM_AmoebaMultipoleForce_Covalent15,
                           covalentMap);

        loadCovalentArray (*(polgrp__.np11 + ii),
                           (polgrp__.ip11 + polgrp__.maxp11*ii), covalentMap);  
        OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                           OpenMM_AmoebaMultipoleForce_PolarizationCovalent11,
                           covalentMap);

        loadCovalentArray (*(polgrp__.np12 + ii),
                           (polgrp__.ip12 + polgrp__.maxp12*ii), covalentMap);  
        OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                           OpenMM_AmoebaMultipoleForce_PolarizationCovalent12,
                           covalentMap);

        loadCovalentArray (*(polgrp__.np13 + ii),
                           (polgrp__.ip13 + polgrp__.maxp13*ii), covalentMap);  
        OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                           OpenMM_AmoebaMultipoleForce_PolarizationCovalent13,
                           covalentMap);

        loadCovalentArray (*(polgrp__.np14 + ii),
                           (polgrp__.ip14 + polgrp__.maxp14*ii), covalentMap);  
        OpenMM_AmoebaMultipoleForce_setCovalentMap( amoebaMultipoleForce, ii,
                           OpenMM_AmoebaMultipoleForce_PolarizationCovalent14,
                           covalentMap);
    }

    OpenMM_DoubleArray_destroy (dipoles);
    OpenMM_DoubleArray_destroy (quadrupoles);

    OpenMM_IntArray_destroy (covalentMap);
}

static double getObcShct (int atmnum) {

    double shct = 0.80;
    if (atmnum == 1)  shct = 0.85;
    if (atmnum == 6)  shct = 0.72;
    if (atmnum == 7)  shct = 0.79;
    if (atmnum == 8)  shct = 0.85;
    if (atmnum == 9)  shct = 0.88;
    if (atmnum == 15)  shct = 0.86;
    if (atmnum == 16)  shct = 0.96;
    if (atmnum == 26)  shct = 0.88;
    return shct;
}

static void setupAmoebaGeneralizedKirkwoodForce (OpenMM_System* system,
                                         int includeCavityTerm, FILE* log) {

    int ii;
    int useGrycuk;
    int useObc;
    double* polePtr;
    double shct;
    char buffer[80];

    // check that Born radius type is recognized -- if not exit
    // after printing messaage; force Born type='Grycuk' for now

    setNullTerminator (solute__.borntyp, 8, buffer);
    useGrycuk = 1;
    if (strncasecmp (buffer, "GRYCUK", 6 ) != 0) {
        if (log) {
            (void) fprintf (log, "setupAmoebaGeneralizedKirkwoodForce: Born type=%s -- forcing Born type to 'Grycuk'.\n", buffer);
        } else {
            (void) fprintf (stderr, "setupAmoebaGeneralizedKirkwoodForce: Born type=%s -- forcing Born type to 'Grycuk'.\n", buffer);
        }
    }

    OpenMM_AmoebaGeneralizedKirkwoodForce* amoebaGeneralizedKirkwoodForce;
    amoebaGeneralizedKirkwoodForce =
        OpenMM_AmoebaGeneralizedKirkwoodForce_create ();
    OpenMM_System_addForce (system, (OpenMM_Force*)
                            amoebaGeneralizedKirkwoodForce);

    OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric
                           (amoebaGeneralizedKirkwoodForce, 78.3);
    OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric
                           (amoebaGeneralizedKirkwoodForce, chgpot__.dielec);
    OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm
                           (amoebaGeneralizedKirkwoodForce, includeCavityTerm);

    // use default values for following fields

    // OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce_setDielectricOffset
    //                     (amoebaGeneralizedGeneralizedKirkwoodForce,
    //                      solute__.doffset*OpenMM_NmPerAngstrom);
    // OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce_setProbeRadius
    //                     (OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce,
    //                      0.14);
    // OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce_setSurfaceAreaFactor
    //                     (amoebaGeneralizedGeneralizedKirkwoodForce,
    //                      surfaceAreaFactor);

    // set parameters for GeneralizedKirkwood

    for (ii = 0; ii < atoms__.n; ii++) {
        polePtr = mpole__.pole + ii*mpole__.maxpole;
        if (useGrycuk) {
            shct = solute__.shct[ii];
        } else {
            shct = getObcShct (atomid__.atomic[ii]);
        }
        OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle
                           (amoebaGeneralizedKirkwoodForce, *(polePtr),
                            OpenMM_NmPerAngstrom*solute__.rsolv[ii], shct);
    }
}

/*
 *    ############################################################
 *            Setup Positions, Velocities and Constraints
 *    ############################################################
 */

static void setupPositions (OpenMM_Vec3Array* initialPosInNm, FILE* log) {

    int ii;
    for (ii = 0; ii < atoms__.n; ii++) {
        OpenMM_Vec3 posInNm;
        posInNm.x = atoms__.x[ii]*OpenMM_NmPerAngstrom;
        posInNm.y = atoms__.y[ii]*OpenMM_NmPerAngstrom;
        posInNm.z = atoms__.z[ii]*OpenMM_NmPerAngstrom;
        OpenMM_Vec3Array_append (initialPosInNm, posInNm);
    }
}

static void setupVelocities (OpenMM_Vec3Array* initialVelInNm, FILE* log) {

    int ii;
    for (ii = 0; ii < atoms__.n; ii++) {
        OpenMM_Vec3 velInNm;
        int offset;
        offset = 3*ii;
        velInNm.x = moldyn__.v[offset]*OpenMM_NmPerAngstrom;
        velInNm.y = moldyn__.v[offset+1]*OpenMM_NmPerAngstrom;
        velInNm.z = moldyn__.v[offset+2]*OpenMM_NmPerAngstrom;
        OpenMM_Vec3Array_append (initialVelInNm, velInNm);
    }
}

static void setupConstraints (OpenMM_System* system, FILE* log) {

    int ii;
    for (ii = 0; ii < freeze__.nrat; ii++) {
        OpenMM_System_addConstraint (system, *(freeze__.irat+2*ii) -1,
                                     *(freeze__.irat + 2*ii +1)-1,
                            (*(freeze__.krat +ii))*OpenMM_NmPerAngstrom);
    }
}

/*
 *    ############################################################
 *            Platform for Calculation: Reference or CUDA
 *    ############################################################
 */

static OpenMM_Platform* getReferencePlatform (FILE* log) {
   
    OpenMM_Platform* platform = OpenMM_Platform_getPlatformByName ("Reference");
    if (platform == NULL) {
        if (log) {
            (void) fprintf (log, "Reference Platform Unavailable\n");
        }
        return platform;
    }
    return platform;
}

static OpenMM_Platform* getCUDAPlatform (FILE* log) {
   
    OpenMM_Platform* platform = OpenMM_Platform_getPlatformByName ("CUDA");
    if (platform == NULL) {
        if (log) {
            (void) fprintf (log, "\n CUDA Platform Unavailable\n");
        }
        return platform;
    }

    const char* deviceId = getenv ("CUDA_DEVICE");
    if (deviceId != NULL) {
        OpenMM_Platform_setPropertyDefaultValue (platform, "CUDADevice",
                                                 deviceId);
        if (log) {
            (void) fprintf (log, "\n Platform CUDA: Setting Device ID to %s from Env Variable CUDA_DEVICE\n", deviceId);
        }
    } else if (log) {
        (void) fprintf (log, "\n Platform CUDA: Using Default Value for CUDA Device ID\n");
    }

    OpenMM_Platform_setPropertyDefaultValue (platform, "CudaPrecision",
                                             "double" );

    return platform;
}

/*
 *    ############################################################
 *                 Initialize OpenMM Data Structures
 *    ############################################################
 *
 *    Following actions are performed here:
 *    (1) Load any available OpenMM plugins
 *    (2) Allocate a OpenMMData structure to hang on to OpenMM data
 *        structures in a manner which is opaque to the caller
 *    (3) Initialize the OpenMM::System with the AMOEBA force field
 *        objects 
 *    (4) Create an Integrator and a Context associating the
 *        Integrator with the System
 *    (5) Select the OpenMM platform to be used.
 *    (6) Return an opaque pointer to the OpenMMData struct
 */

extern "C" {
void openmm_init_ (void** ommHandle, double* dt) {

    int ii;
    int mdMode = 0;
    int isoThermal = 0;
    int removeConstrainedCovalentIxns = 0;
    char buffer[128];
    FILE* log = stderr;

    // Allocate space for opaque handle to hold OpenMM objects
    // such as system, integrator, context, etc.

    OpenMMData* omm = (OpenMMData*) malloc(sizeof(struct OpenMMData_s));

    // These are temporary OpenMM objects used and discarded here

    OpenMM_Vec3Array*       initialPosInNm;
    OpenMM_Vec3Array*       initialVelInNm;
    OpenMM_StringArray*     pluginList;
    OpenMM_Platform*        platform;

    // Load all OpenMM plugin libraries from their default location;
    // Call the plugin loading routine twice to fix an issue with OSX
    // where the first library in the alphabetical list gets skipped

    pluginList = OpenMM_Platform_loadPluginsFromDirectory
                     (OpenMM_Platform_getDefaultPluginsDirectory());
    pluginList = OpenMM_Platform_loadPluginsFromDirectory
                     (OpenMM_Platform_getDefaultPluginsDirectory());
    (void) fprintf (stderr, "\n Default OpenMM Plugin Directory: %s\n\n",
                        OpenMM_Platform_getDefaultPluginsDirectory());
    for (ii = 0; ii < OpenMM_StringArray_getSize(pluginList); ii++)
    {
        (void) fprintf (stderr, " Plugin Library: %s\n",
                            OpenMM_StringArray_get(pluginList, ii));
    }
    OpenMM_StringArray_destroy (pluginList);
    (void) fflush (NULL);

    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: the OpenMM
    // System takes ownership of the force objects; don't delete them yourself

    omm->system = OpenMM_System_create ();
    setupSystemParticles (omm->system, log);
    setupCMMotionRemover (omm->system, log);

    if (potent__.use_bond) {
        setupAmoebaBondForce (omm->system,
                              removeConstrainedCovalentIxns, log);
    }

    if (potent__.use_angle) {
        setupAmoebaAngleForce (omm->system,
                               removeConstrainedCovalentIxns, log);
        setupAmoebaInPlaneAngleForce (omm->system,
                                      removeConstrainedCovalentIxns, log);
    }

    if (potent__.use_strbnd) {
        setupAmoebaStretchBendForce (omm->system,
                                     removeConstrainedCovalentIxns, log);
    }

    if (potent__.use_urey) {
        setupAmoebaUreyBradleyForce (omm->system,
                                     removeConstrainedCovalentIxns, log);
    }

    if (potent__.use_opbend) {
        setupAmoebaOutOfPlaneBendForce (omm->system, log);
    }

    if (potent__.use_tors) {
        setupAmoebaTorsionForce (omm->system, log);
    }

    if (potent__.use_pitors) {
        setupAmoebaPiTorsionForce (omm->system, log);
    }

    if (potent__.use_tortor) {
        setupAmoebaTorsionTorsionForce (omm->system, log);
    }

    if (potent__.use_vdw) {
        setupAmoebaVdwForce (omm->system, log);
    }

    if (potent__.use_mpole){
        setupAmoebaMultipoleForce (omm->system, log);
        if (potent__.use_solv) {
	    setupAmoebaGeneralizedKirkwoodForce (omm->system, 1, log);
        }
    }

    if (potent__.use_solv) {
        setupAmoebaWcaDispersionForce (omm->system, log);
    }

    if (bath__.isobaric && bath__.atmsph > 0.0) {
        mdMode = 4;
        setupMonteCarloBarostat (omm->system, log);
    }

    setNullTerminator (bath__.thermostat, 11, buffer);
    if (strncasecmp( buffer, "ANDERSEN", 8 ) == 0) {
        isoThermal = 1;
        setupAndersenThermostat (omm->system, log);
    }

    // setup of constraints, positions and velocities

    setupConstraints (omm->system, log) ;

    initialPosInNm = OpenMM_Vec3Array_create (0);
    setupPositions (initialPosInNm, log);

    initialVelInNm = OpenMM_Vec3Array_create (0);
    setupVelocities (initialVelInNm, log);

    // Choose an Integrator, and a Context connecting the System with the
    // Integrator. Let the Context choose the best available Platform.
    // Initialize the configuration from default positions collected above.

    setNullTerminator (mdstuf__.integrate, 10, buffer);
    if (strncasecmp (buffer, "VERLET", 6) == 0) {

        omm->integrator = (OpenMM_Integrator*)OpenMM_VerletIntegrator_create
                               (*dt);
        if (mdMode == 4 && isoThermal == 0) {
            setupAndersenThermostat (omm->system, log);
        }

    } else if (strncasecmp (buffer, "STOCHASTIC", 10) == 0) {

        if (log) {
            (void) fprintf (log, "\n Stochastic Integrator:\n");
            (void) fprintf (log, "\n Temperature          %15.4f K",
                                bath__.kelvin );
            (void) fprintf (log, "\n Friction             %15.4f ps^(-1)",
                                stodyn__.friction );
            (void) fprintf (log, "\n TimeStep             %15.4f ps\n",
                                *dt );
            (void) fflush (log);
        }
        omm->integrator = (OpenMM_Integrator*)OpenMM_LangevinIntegrator_create
                               (bath__.kelvin, stodyn__.friction, *dt);

    } else {
        (void) fprintf (stderr, "\n Integrator %s is Not Supported\n", buffer);
        (void) fflush (stderr);
        exit (-1);
    }

    //platform = getReferencePlatform (log);
    platform = getCUDAPlatform (log);
    if (platform == NULL) {
        exit (-1);
    }

    omm->context = OpenMM_Context_create_2 (omm->system, omm->integrator,
                                            platform);
    //(void) fprintf (log, "\n OpenMMDataHandle:  %x\n", (void*)(omm));
    //(void) fprintf (log, "\n Integrator:  %x\n", (void*)(omm->integrator));

    OpenMM_Context_setPositions (omm->context, initialPosInNm);
    OpenMM_Context_setVelocities (omm->context, initialVelInNm);
    {
        int arraySz;
        int maxPrint;
        double x1, x2, x3;
        double v1, v2, v3;
        arraySz = OpenMM_Vec3Array_getSize (initialPosInNm);
        maxPrint = 5;
        (void) fprintf (log, "\n Initial Positions and Velocities:\n\n"); 
        for (ii = 0; ii < arraySz; ii++)
        {
            x1 = (*OpenMM_Vec3Array_get(initialPosInNm, ii)).x
                      / OpenMM_NmPerAngstrom;
            x2 = (*OpenMM_Vec3Array_get(initialPosInNm, ii)).y
                      / OpenMM_NmPerAngstrom;
            x3 = (*OpenMM_Vec3Array_get(initialPosInNm, ii)).z
                      / OpenMM_NmPerAngstrom;
            v1 = (*OpenMM_Vec3Array_get(initialVelInNm, ii)).x
                      / OpenMM_NmPerAngstrom;
            v2 = (*OpenMM_Vec3Array_get(initialVelInNm, ii)).y
                      / OpenMM_NmPerAngstrom;
            v3 = (*OpenMM_Vec3Array_get(initialVelInNm, ii)).z
                     / OpenMM_NmPerAngstrom;
            (void) fprintf (log, "%7d  POS    %16.7e %16.7e %16.7e\n",
                            ii+1, x1, x2, x3);
            (void) fprintf (log, "%7d  VEL    %16.7e %16.7e %16.7e\n",
                            ii+1, v1, v2, v3);
            if (ii == maxPrint-1 && ii < (arraySz-maxPrint-1))
                ii = arraySz - maxPrint - 1;
        }
    }

    *ommHandle = (void*) omm;
}

/*
 *    ############################################################
 *          Copy State from OpenMM to TINKER Data Structures
 *    ############################################################
 *
 *    @param omm                    handle with OpenMM data structures
 *    @param timeInPs               output simulation time in ps
 *    @param kineticEnergyInKcal    output kinetic energy in kcal
 *    @param potentialEnergyInKcal  output potential energy in kcal
 */

void openmm_copy_state_ (void** omm, double *timeInPs,
                         double* kineticEnergyInKcal,
                         double* potentialEnergyInKcal) {

    OpenMM_State* state;
    const OpenMM_Vec3Array* positionArray;
    const OpenMM_Vec3Array* velocityArray;
    const OpenMM_Vec3Array* forceArray;
    OpenMM_Vec3 aBox;
    OpenMM_Vec3 bBox;
    OpenMM_Vec3 cBox;
    int infoMask;
    int ii;
    int offset;
    int debug = 0;
    double amass;
    double positionConvert;
    double velocityConvert;
    double forceConvert;
    OpenMMData* openMMDataHandle;

    openMMDataHandle = (OpenMMData*) (*omm);

    infoMask = OpenMM_State_Positions;
    infoMask += OpenMM_State_Velocities;
    infoMask += OpenMM_State_Forces;    
    infoMask += OpenMM_State_Energy;

    // State object is created here and must be explicitly destroyed below

    state = OpenMM_Context_getState (openMMDataHandle->context, infoMask, 0);
    *timeInPs = OpenMM_State_getTime (state);

    OpenMM_State_getPeriodicBoxVectors (state, &aBox, &bBox, &cBox);

    *(boxes__.xbox) = aBox.x / OpenMM_NmPerAngstrom;
    *(boxes__.ybox) = bBox.y / OpenMM_NmPerAngstrom;
    *(boxes__.zbox) = cBox.z / OpenMM_NmPerAngstrom;

    *(boxes__.xbox2) = 0.5 * (*(boxes__.xbox));
    *(boxes__.ybox2) = 0.5 * (*(boxes__.ybox));
    *(boxes__.zbox2) = 0.5 * (*(boxes__.zbox));

    // Positions, velocities and forces are maintained as Vec3Arrays
    // inside the State. This gives us access, but don't destroy them
    // yourself as they will go away with the State

    positionConvert = 1.0 / OpenMM_NmPerAngstrom;
    velocityConvert = 1.0 / OpenMM_NmPerAngstrom;
    forceConvert = 10.0;

    positionArray = OpenMM_State_getPositions (state);
    velocityArray = OpenMM_State_getVelocities (state);
    forceArray = OpenMM_State_getForces (state);

    for (ii = 0; ii < atoms__.n; ii++) {
        atoms__.x[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).x
                             * positionConvert;
        atoms__.y[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).y
                             * positionConvert;
        atoms__.z[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).z
                             * positionConvert;
    }

    for (ii = 0; ii < atoms__.n; ii++) {
        offset = 3*ii;
        moldyn__.v[offset] = (*OpenMM_Vec3Array_get(velocityArray, ii)).x
                                  * velocityConvert;
        moldyn__.v[offset+1] = (*OpenMM_Vec3Array_get(velocityArray, ii)).y
                                    * velocityConvert;
        moldyn__.v[offset+2] = (*OpenMM_Vec3Array_get(velocityArray, ii)).z
                                    * velocityConvert;
    }

    for (ii = 0; ii < atoms__.n; ii++) {
        offset = 3*ii;
        amass = 1.0 / atomid__.mass[ii];
        moldyn__.a[offset] = (*OpenMM_Vec3Array_get(forceArray, ii)).x
                                    * amass * forceConvert;
        moldyn__.a[offset+1] = (*OpenMM_Vec3Array_get(forceArray, ii)).y
                                    * amass * forceConvert;
        moldyn__.a[offset+2] = (*OpenMM_Vec3Array_get(forceArray, ii)).z
                                    * amass * forceConvert;
    }

    for (ii = 0; ii < atoms__.n; ii++) {
        offset = 3*ii;
        moldyn__.aalt[offset] = 0.0;
        moldyn__.aalt[offset+1] = 0.0;
        moldyn__.aalt[offset+2] = 0.0;
    }

    if (debug ) {

        (void) fprintf (stderr, "State: E=%15.7e [%15.7e %15.7e] t=%15.7e ps\n",
                            (*kineticEnergyInKcal + *potentialEnergyInKcal),
                            *kineticEnergyInKcal, *potentialEnergyInKcal,
                            *timeInPs);

        for (ii = 0; ii < atoms__.n; ii++) {
            (void) fprintf (stderr, "%7d  POS    %16.7e %16.7e %16.7e\n", ii+1,
                (*OpenMM_Vec3Array_get(positionArray,ii)).x*positionConvert,
                (*OpenMM_Vec3Array_get(positionArray,ii)).y*positionConvert,
                (*OpenMM_Vec3Array_get(positionArray,ii)).z*positionConvert);
            (void) fprintf (stderr, "%7d  VEL    %16.7e %16.7e %16.7e\n", ii+1,
                (*OpenMM_Vec3Array_get(velocityArray,ii)).x*velocityConvert,
                (*OpenMM_Vec3Array_get(velocityArray,ii)).y*velocityConvert,
                (*OpenMM_Vec3Array_get(velocityArray,ii)).z*velocityConvert);
            (void) fprintf (stderr, "%7d  FRC    %16.7e %16.7e %16.7e\n", ii+1,
                (*OpenMM_Vec3Array_get(forceArray,ii)).x*forceConvert,
                (*OpenMM_Vec3Array_get(forceArray,ii)).y*forceConvert,
                (*OpenMM_Vec3Array_get(forceArray,ii)).z*forceConvert);
        }
        (void) fflush (stderr);
    }

    // convert energies from kJ/mol to kcal/mol

    *kineticEnergyInKcal = OpenMM_State_getKineticEnergy (state)
                               * OpenMM_KcalPerKJ;
    *potentialEnergyInKcal = OpenMM_State_getPotentialEnergy (state)
                                 * OpenMM_KcalPerKJ;
 
    OpenMM_State_destroy (state);
}

/*
 *    ############################################################
 *       Update TINKER Data Structures; Call mdstat and mdsave
 *    ############################################################
 *
 *    @param omm          handle containing OpenMM data structures
 *    @param dt           simulation time step in ps
 *    @param istep        current step in MD loop
 *    @param callMdStat   if nonzero, call TINKER mdstat routine
 *    @param callMdSave   if nonzero, call TINKER mdsave routine
 */

void openmm_update_ (void** omm, double* dt, int* istep,
                     int* callMdStat, int* callMdSave ) {

    OpenMM_State* state;
    const OpenMM_Vec3Array* positionArray;
    const OpenMM_Vec3Array* velocityArray;
    const OpenMM_Vec3Array* forceArray;
    OpenMM_Vec3 aBox;
    OpenMM_Vec3 bBox;
    OpenMM_Vec3 cBox;
    int infoMask;
    int ii;
    double amass;
    double positionConvert;
    double velocityConvert;
    double forceConvert;

    static const double gasconst = 1.98720415e-3;
    int debug = 0;

    double totalEnergy, potentialEnergy, kineticEnergy;

    OpenMMData* openMMDataHandle;

    openMMDataHandle = (OpenMMData*) (*omm);

    infoMask = OpenMM_State_Positions;
    infoMask += OpenMM_State_Velocities;
    infoMask += OpenMM_State_Forces;    
    infoMask += OpenMM_State_Energy;

    // state object is created here and must be destroyed below

    state = OpenMM_Context_getState (openMMDataHandle->context, infoMask, 0);
    OpenMM_State_getPeriodicBoxVectors (state, &aBox, &bBox, &cBox);

    *(boxes__.xbox) = aBox.x / OpenMM_NmPerAngstrom;
    *(boxes__.ybox) = bBox.y / OpenMM_NmPerAngstrom;
    *(boxes__.zbox) = cBox.z / OpenMM_NmPerAngstrom;

    *(boxes__.xbox2) = 0.5 * (*(boxes__.xbox));
    *(boxes__.ybox2) = 0.5 * (*(boxes__.ybox));
    *(boxes__.zbox2) = 0.5 * (*(boxes__.zbox));

    // fprintf (stderr, "openmm_update_ %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n",
    // *(boxes__.xbox), *(boxes__.ybox), *(boxes__.zbox),
    // aBox.x, bBox.y, cBox.z );

    // load positions/velocities and energies

    positionConvert = 1.0 / OpenMM_NmPerAngstrom;
    velocityConvert = 1.0 / OpenMM_NmPerAngstrom;
    forceConvert = 10.0;

    positionArray = OpenMM_State_getPositions (state);
    velocityArray = OpenMM_State_getVelocities (state);
    forceArray = OpenMM_State_getForces (state);

    for (ii = 0; ii < atoms__.n; ii++) {
        atoms__.x[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).x
                            * positionConvert;
        atoms__.y[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).y
                            * positionConvert;
        atoms__.z[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).z
                            * positionConvert;
    }

    for (ii = 0; ii < atoms__.n; ii++) {
        int offset = ii*3;
        moldyn__.v[offset] = (*OpenMM_Vec3Array_get(velocityArray, ii)).x
                                  * velocityConvert;
        moldyn__.v[offset+1] = (*OpenMM_Vec3Array_get(velocityArray, ii)).y
                                    * velocityConvert;
        moldyn__.v[offset+2] = (*OpenMM_Vec3Array_get(velocityArray, ii)).z
                                    * velocityConvert;
    }

    for (ii = 0; ii < atoms__.n; ii++) {
        int offset = ii*3;
        amass = 1.0 / atomid__.mass[ii];
        moldyn__.a[offset] = (*OpenMM_Vec3Array_get(forceArray, ii)).x
                                  * amass * forceConvert;
        moldyn__.a[offset+1] = (*OpenMM_Vec3Array_get(forceArray, ii)).y
                                    * amass * forceConvert;
        moldyn__.a[offset+2] = (*OpenMM_Vec3Array_get(forceArray, ii)).z
                                    * amass * forceConvert;
    }

    for (ii = 0; ii < atoms__.n; ii++) {
        int offset = ii*3;
        moldyn__.aalt[offset] = 0.0;
        moldyn__.aalt[offset+1] = 0.0;
        moldyn__.aalt[offset+2] = 0.0;
    }

    potentialEnergy = (OpenMM_State_getPotentialEnergy(state))
                          * OpenMM_KcalPerKJ;
    kineticEnergy = (OpenMM_State_getKineticEnergy(state))
                        * OpenMM_KcalPerKJ;
    totalEnergy = potentialEnergy + kineticEnergy;
 
    if (debug) {

        double eksum, temp, pres;
        double ekin[3][3];

        kinetic_ (&eksum, &ekin);
        temp = 2.0 * eksum / ((double)(mdstuf__.nfree)*gasconst);

        (void) fprintf (stderr, "State: E=%15.7e [%15.7e %15.7e] t=%15.7e ps T=%12.5e InitT=%12.3f\n",
                        totalEnergy, potentialEnergy, kineticEnergy,
                        (*dt)*(*istep), temp, bath__.kelvin );

        for (ii = 0; ii < atoms__.n; ii++) {
            (void) fprintf (stderr, "%7d  POS    %16.7e %16.7e %16.7e\n", ii+1,
                (*OpenMM_Vec3Array_get(positionArray,ii)).x*positionConvert,
                (*OpenMM_Vec3Array_get(positionArray,ii)).y*positionConvert,
                (*OpenMM_Vec3Array_get(positionArray,ii)).z*positionConvert);
            (void) fprintf (stderr, "%7d  VEL    %16.7e %16.7e %16.7e\n", ii+1,
                (*OpenMM_Vec3Array_get(velocityArray,ii)).x*velocityConvert,
                (*OpenMM_Vec3Array_get(velocityArray,ii)).y*velocityConvert,
                (*OpenMM_Vec3Array_get(velocityArray,ii)).z*velocityConvert);
            (void) fprintf (stderr, "%7d  FRC    %16.7e %16.7e %16.7e\n", ii+1,
                (*OpenMM_Vec3Array_get(forceArray,ii)).x*forceConvert,
                (*OpenMM_Vec3Array_get(forceArray,ii)).y*forceConvert,
                (*OpenMM_Vec3Array_get(forceArray,ii)).z*forceConvert);
        }
        (void) fflush (stderr);
    }

    // make calls to mdstat and/or mdsave if flags are set

    if (*callMdStat || *callMdSave) {
        double eksum, temp, pres;
        double ekin[3][3];
        kinetic_ (&eksum, &ekin);

        if (*callMdStat) {
            temp = 2.0 * eksum / ((double)(mdstuf__.nfree)*gasconst);
            pres = 0.0;
            mdstat_ (istep,dt,&totalEnergy,&potentialEnergy,&eksum,&temp,&pres);
        }
        if (*callMdSave) {
            bounds_ ();
            mdsave_ (istep,dt,&potentialEnergy,&eksum);
        }
    }
    OpenMM_State_destroy (state);
}

void openmm_get_energies_ (void** omm, double* ke, double* pe ) {

    OpenMM_State* state;
    const OpenMM_Vec3Array* positionArray;
    OpenMM_Vec3 aBox;
    OpenMM_Vec3 bBox;
    OpenMM_Vec3 cBox;
    int infoMask;
    int ii;

    static const double gasconst = 1.98720415e-3;
    int debug = 0;

    double totalEnergy, potentialEnergy, kineticEnergy;
    OpenMMData* openMMDataHandle;

    openMMDataHandle = (OpenMMData*) (*omm);

    infoMask = OpenMM_State_Positions;
    infoMask += OpenMM_State_Velocities;
    if (debug) {
        infoMask += OpenMM_State_Forces;    
    }
    infoMask += OpenMM_State_Energy;

    // State object is created here and must be explicitly destroyed below

    state = OpenMM_Context_getState (openMMDataHandle->context, infoMask, 0);
    OpenMM_State_getPeriodicBoxVectors (state, &aBox, &bBox, &cBox);

    // load positions, velocities and energies

    //positionArray = OpenMM_State_getPositions (state);
    //velocityArray = OpenMM_State_getVelocities (state);

    potentialEnergy = (OpenMM_State_getPotentialEnergy (state))
                          * OpenMM_KcalPerKJ;
    kineticEnergy = (OpenMM_State_getKineticEnergy (state)) * OpenMM_KcalPerKJ;
    totalEnergy = potentialEnergy + kineticEnergy;

    *ke = kineticEnergy;
    *pe = potentialEnergy;

    OpenMM_State_destroy (state);
}

void openmm_update_box_ (void** omm ) {

    OpenMM_Vec3 aBox;
    OpenMM_Vec3 bBox;
    OpenMM_Vec3 cBox;

    OpenMMData* openMMDataHandle;
    openMMDataHandle = (OpenMMData*) (*omm);

    computePeriodicBoxVectors(aBox, bBox, cBox);

    OpenMM_Context_setPeriodicBoxVectors (openMMDataHandle->context,
                                          &aBox, &bBox, &cBox);
}

void openmm_update_positions_ (void** omm ) {

    OpenMMData* openMMDataHandle;
    OpenMM_Vec3Array* initialPosInNm;
    FILE* log = stderr;

    openMMDataHandle = (OpenMMData*) (*omm);

    initialPosInNm = OpenMM_Vec3Array_create (0);
    setupPositions (initialPosInNm, log);
    OpenMM_Context_setPositions (openMMDataHandle->context, initialPosInNm);
    OpenMM_Vec3Array_destroy (initialPosInNm);
}

void openmm_take_steps_ (void** omm, int* numSteps) {

    OpenMMData* openMMDataHandle = (OpenMMData*) (*omm);
    OpenMM_Integrator_step (openMMDataHandle->integrator, *numSteps);
}

void openmm_cleanup_ (void** omm) {

    // clean up top-level heap allocated objects we are done with

    OpenMMData* openMMDataHandle = (OpenMMData*) (*omm);
    OpenMM_Context_destroy (openMMDataHandle->context);
    OpenMM_Integrator_destroy (openMMDataHandle->integrator);
    OpenMM_System_destroy (openMMDataHandle->system);
    free (openMMDataHandle);
}
}

static void zeroTinkerForce (double* tinkerForce) {

    int ii;
    for (ii = 0; ii < atoms__.n; ii++){
        *(tinkerForce + 3*ii + 0) = 0.0;
        *(tinkerForce + 3*ii + 1) = 0.0;
        *(tinkerForce + 3*ii + 2) = 0.0;
    }
}

static void zeroVec3Force (OpenMM_Vec3Array* tinkerForce) {

    int ii;
    for (ii = 0; ii < atoms__.n; ii++) {
        OpenMM_Vec3 force;
        force.x = force.y = force.z = 0.0;
        OpenMM_Vec3Array_set (tinkerForce, ii, force);
    }
}

static void setTinker1DArray (int size, double* tinkerArray, double value) {

    int ii;
    for (ii = 0; ii < size; ii++) {
        tinkerArray[ii] = value;
    }
}

static void loadTinkerForce (double* tinkerForce, int add,
                             OpenMM_Vec3Array* arrayToLoad) {

    int ii;
    if (add == 0) {
        for (ii = 0; ii < atoms__.n; ii++) {
            OpenMM_Vec3 force;
            force.x = *(tinkerForce + 3*ii);
            force.y = *(tinkerForce + 3*ii + 1);
            force.z = *(tinkerForce + 3*ii + 2);
            OpenMM_Vec3Array_set (arrayToLoad, ii, force);
        }
    } else {
        double factor = (double) add;
        for (ii = 0; ii < atoms__.n; ii++) {
            OpenMM_Vec3 force;
            const OpenMM_Vec3* currentForce = OpenMM_Vec3Array_get
                                                  (arrayToLoad, ii);
            force.x = currentForce->x + factor * (*(tinkerForce + 3*ii));
            force.y = currentForce->y + factor * (*(tinkerForce + 3*ii + 1));
            force.z = currentForce->z + factor * (*(tinkerForce + 3*ii + 2));
            OpenMM_Vec3Array_set (arrayToLoad, ii, force);
        }
    }
}

static int usingImplicitSolvent (void) {

    int implicitSolventActive = -1;
    char solvatationType[16];
    char bornType[16];

    setNullTerminator (solute__.solvtyp, 8, solvatationType);
    setNullTerminator (solute__.borntyp, 8, bornType);

    // return <0 if parameter/option combination is unsupported
    //         0 if explicit solvent (Ewald is in use)
    //         1 if GK implicit solvent via OBC
    //         2 if GK implicit solvent via Grycuk

    if( strncasecmp (solvatationType, "GK", 2) == 0 &&
            ((strncasecmp (bornType, "OBC", 3) == 0) ||
             (strncasecmp (bornType, "GRYCUK", 6) == 0))) {
        if (limits__.use_ewald) {
            implicitSolventActive = -2;
        } else {
            if ((strncasecmp( bornType, "OBC", 3) == 0)) {
                implicitSolventActive = 1;
            } else {
                implicitSolventActive = 2;
            }
        }
    } else if (limits__.use_ewald) {
        implicitSolventActive = 0;
    }
    return implicitSolventActive;
}

/*
 *    ############################################################
 *       Map Data Structures and Check for Supported Potentials
 *    ############################################################
 *
 *    (1) Map TINKER data structures with arrays to C data structs
 *    (2) Check for TINKER parameter settings consistent with OpenMM
 *        values and potential types
 *    (3) If invalid match or setting found, exit
 */

extern "C" {
void openmm_validate_ (int* testingActive) {

    char buffer[128];
    int ii;
    int implicitSolventActive;
    int invalidAxisType;
    int invalid = 0;
    int testing = *testingActive;
    const char* ixnString = "Interaction Not Currently Supported by OpenMM.\n";
    FILE* log = stderr;

    if (testing && log) {
        (void) fprintf (log, "Testing Mode is Active (=%d)\n", testing); 
    }

    // check that requested bond potential type is supported

    if (strncasecmp (bndpot__.bndtyp, "HARMONIC", 8) != 0) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Only HARMONIC bond potential supported: %s not available\n",
                            bndpot__.bndtyp); 
        }
    }
 
    // check that requested angle potential type is supported

    for (ii = 0; ii < angbnd__.nangle; ii++) {
        if (strncasecmp (angpot__.angtyp + 8*ii, "HARMONIC", 8) != 0 &&
            strncasecmp (angpot__.angtyp + 8*ii, "IN-PLANE", 8) != 0 ) { 
            invalid++;
            if (log) {
                setNullTerminator (angpot__.angtyp + 8*ii, 8, buffer);
                (void) fprintf (log, "Only HARMONIC or IN-PLANE angle potential supported: %s not available\n",
                                buffer);
            }
        }
    }

    // check that scaling parameters agree with hardcoded values

    if (fabs( 0.0 - vdwpot__.v2scale)  > 0.0 || 
        fabs( 0.0 - vdwpot__.v3scale)  > 0.0 || 
        fabs( 1.0 - vdwpot__.v4scale)  > 0.0 || 
        fabs( 1.0 - vdwpot__.v5scale)  > 0.0 ) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Vdw Scale Parameters [%8.1f %8.1f %8.1f %8.1f] are not [0.0, 0.0, 1.0, 1.0]\n",
                            vdwpot__.v2scale, vdwpot__.v3scale,
                            vdwpot__.v4scale, vdwpot__.v5scale);
        }
    }

    if (fabs( 0.0 - polpot__.p2scale) > 0.0 ||
        fabs( 0.0 - polpot__.p3scale) > 0.0 ||
        fabs( 1.0 - polpot__.p4scale) > 0.0 ||
        fabs( 0.5 - polpot__.p41scale) > 0.0 ||
        fabs( 1.0 - polpot__.p5scale) > 0.0) {
        invalid++;
        if (log) {
            (void) fprintf (log, "P-Scale Parameters [%8.1f %8.1f %8.1f %8.1f %8.1f] are not [0.0, 0.0, 1.0, 0.5, 1.0]\n",
                            polpot__.p2scale, polpot__.p3scale,
                            polpot__.p4scale, polpot__.p41scale,
                            polpot__.p5scale);
        }
    }

    if (fabs( 0.0 - polpot__.d1scale) > 0.0 ||
        fabs( 1.0 - polpot__.d2scale) > 0.0 ||
        fabs( 1.0 - polpot__.d3scale) > 0.0 ||
        fabs( 1.0 - polpot__.d4scale) > 0.0) {
        invalid++;
        if (log) {
            (void) fprintf (log, "D-Scale Parameters [%8.1f %8.1f %8.1f %8.1f] are not [0.0, 1.0, 1.0, 1.0]\n",
                            polpot__.d1scale, polpot__.d2scale,
                            polpot__.d3scale, polpot__.d4scale);
        }
    }

    if (fabs( 1.0 - polpot__.u1scale) > 0.0 ||
        fabs( 1.0 - polpot__.u2scale) > 0.0 ||
        fabs( 1.0 - polpot__.u3scale) > 0.0 ||
        fabs( 1.0 - polpot__.u4scale) > 0.0) {
        invalid++;
        if (log) {
            (void) fprintf (log, "U-Scale Parameters [%8.1f %8.1f %8.1f %8.1f] are not [1.0, 1.0, 1.0, 1.0]\n",
                            polpot__.u1scale, polpot__.u2scale,
                            polpot__.u3scale, polpot__.u4scale);
        }
    }

    if (fabs( 0.0 - mplpot__.m2scale) > 0.0 ||
        fabs( 0.0 - mplpot__.m3scale) > 0.0 ||
        fabs( 0.4 - mplpot__.m4scale) > 0.0 ||
        fabs( 0.8 - mplpot__.m5scale) > 0.0) {
        invalid++;
        if (log) {
            (void) fprintf (log, "M-Scale Parameters [%8.1f %8.1f %8.1f %8.1f] are not [0.0, 0.0, 0.4, 0.8]\n",
                            mplpot__.m2scale, mplpot__.m3scale,
                            mplpot__.m4scale, mplpot__.m5scale);
        }
    }

    // check that local coordinate frame axis type is supported

    if (potent__.use_mpole) {
        for (ii = 0; ii < atoms__.n; ii++) {
            invalidAxisType = 0;
            if (strncasecmp (mpole__.polaxe + 8*ii, "Z-then-X", 8) != 0 &&
                strncasecmp (mpole__.polaxe + 8*ii, "Bisector", 8) != 0 &&
                strncasecmp (mpole__.polaxe + 8*ii, "Z-Bisect", 8) != 0 &&
                strncasecmp (mpole__.polaxe + 8*ii, "3-Fold", 6) != 0 &&
                strncasecmp (mpole__.polaxe + 8*ii, "None", 4) != 0 &&
                strncasecmp (mpole__.polaxe + 8*ii, "    ", 4) != 0 &&
                strncasecmp (mpole__.polaxe + 8*ii, "Z-Only", 6) != 0) {
                invalid++;
                if (log) {
                    setNullTerminator (mpole__.polaxe + 8*ii, 8, buffer);
                    (void) fprintf (log, "Axis Type=%s for Atom %7d Not Supported\n",
                                    buffer, ii);
                    if (invalid > 20) {
                        (void) fflush (log);
                        exit (-1);
                    }
                }
            }
        }
    }

    // check for any particles that are "inactive" (ie, frozen)

    if (*(usage__.nuse) !=  atoms__.n) {
        if (log) {
            (void) fprintf (log, "Inactive Atoms Not Supported: nuse=%d, N=%d\n",
                            *(usage__.nuse), atoms__.n);
        }
        invalid++;
    }

    // check for torsional values beyond the supported 3-fold term

    if (potent__.use_tors) {
        double* torsPtr;
        for (ii = 0; ii < tors__.ntors; ii++) {
            torsPtr = tors__.tors4 + ii*4;
            if (*torsPtr != 0.0) {
               if (log) {
                    (void) fprintf (log, "4-Fold Parameter for Torsion %d is NonZero %15.7e\n",
                                    ii, *torsPtr);
                }
                invalid++;
            }
            torsPtr = tors__.tors5 + ii*4;
            if (*torsPtr != 0.0) {
                if (log) {
                    (void) fprintf (log, "5-Fold Parameter for Torsion %d is NonZero %15.7e\n",
                                    ii, *torsPtr);
                }
                invalid++;
            }
            torsPtr = tors__.tors6 + ii*4;
            if (*torsPtr != 0.0) {
                if (log) {
                    (void) fprintf (log, "6-Fold Parameter for Torsion %d is NonZero %15.7e\n",
                                    ii, *torsPtr);
                }
                invalid++;
            }
        }
    }

    // check for use of potential terms not implemented in OpenMM

    if (potent__.use_angang) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Angle-Angle %s", ixnString);
        }
    }
    if (potent__.use_opdist) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Out-of-Plane Distance %s", ixnString);
        }
    }
    if (potent__.use_improp) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Improper Dihedral %s", ixnString);
        }
    }
    if (potent__.use_imptor) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Improper Torsion %s", ixnString);
        }
    }
    if (potent__.use_strtor) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Stretch-Torsion %s", ixnString);
        }
    }
    if (potent__.use_angtor) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Angle-Torsion %s", ixnString);
        }
    }
    if (potent__.use_charge) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Charge-Charge %s", ixnString);
        }
    }
    if (potent__.use_chgdpl) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Atomwise Charge-Dipole %s", ixnString);
        }
    }
    if (potent__.use_dipole) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Dipole-Dipole %s", ixnString);
        }
    }
    if (potent__.use_rxnfld) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Reaction Field %s", ixnString);
        }
    }
    if (potent__.use_metal) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Ligand Field %s", ixnString);
        }
    }
    if (potent__.use_geom) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Geometric Restraint %s", ixnString);
        }
    }
/*
    if (potent__.use_extra) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Extra Potentials Not Supported\n");
        }
    }
*/
    if (potent__.use_orbit) {
        invalid++;
        if (log) {
            (void) fprintf (log, "SCF-Pi-MO Calculation Not Supported\n");
        }
    }
    implicitSolventActive = usingImplicitSolvent ();
    if (implicitSolventActive == -1 && !testing) {
        if (log) {
            (void) fprintf (log, "Either (Solvate= 'GK' && Born-Radius='GRYCUK') or the Ewald Parameter Must be Set (but not both)\n");
        }
        invalid++;
    }
    if (implicitSolventActive == -2 && !testing) {
        if (log) {
            (void) fprintf (log, "Unsupported Force Combination: Both solvate= 'GK' && born-radius='GRYCUK' and Ewald parameter cannot be set\n");
        }
        invalid++;
    }
    setNullTerminator (mdstuf__.integrate, 10, buffer);
    if (strncasecmp (buffer, "VERLET", 6) != 0 &&
        strncasecmp (buffer, "STOCHASTIC", 6) != 0 && !testing ) {
        invalid++;
        if (log) {
            (void) fprintf (log, "Integrator %s is Not Supported\n", buffer);
        }
    }
    if (log) {
        (void) fflush (log);
    }

    if (invalid && !testing) {
        exit (-1);
    }
}
}

/*
 *    ############################################################
 *           Compare TINKER and OpenMM Energies and Forces
 *    ############################################################
 */

extern "C" {
int openmm_test_ (void) {

    OpenMM_Vec3Array* initialPosInNm;
    OpenMM_StringArray* pluginList;
    OpenMM_Platform* platform;
    OpenMM_System* system;
    OpenMM_Context* context;
    OpenMM_Integrator* integrator;
    OpenMM_State* state;

    const OpenMM_Vec3Array* posArrayInNm;
    const OpenMM_Vec3Array* openMMForces;

    int infoMask;
    int ii, jj;
    int countActiveForces;
    string testName;
    double conversion, delta, dot;
    double tinkerNorm, openMMNorm;
    double openMMPotentialEnergy;
    OpenMM_Vec3Array* tinkerForce;
    double tinkerEnergy;
    double maxEDelta, maxFDelta;
    double minDot, avgDot;
    int maxFDeltaIndex, minDotIndex;
    int implicitSolventActive;

    int removeConstrainedCovalentIxns = 0;
    FILE* log = stderr;

    if (log) {
        (void) fprintf (log, "\n Testing TINKER vs OpenMM for AMOEBA:\n");
    }

    //if (log) {
    //    (void) fprintf (log, "\n Default OpenMM Plugin Directory: %s\n\n",
    //                    OpenMM_Platform_getDefaultPluginsDirectory());
    //    (void) fflush (log);
    //    pluginList = OpenMM_Platform_loadPluginsFromDirectory
    //                     (OpenMM_Platform_getDefaultPluginsDirectory());
    //    pluginList = OpenMM_Platform_loadPluginsFromDirectory
    //                     (OpenMM_Platform_getDefaultPluginsDirectory());
    //    for (ii = 0; ii < OpenMM_StringArray_getSize( pluginList ); ii++) {
    //        (void) fprintf (log, " Plugin Library: %s\n",
    //                        OpenMM_StringArray_get(pluginList, ii));
    //    }
    //    (void) fflush (log);
    //}
    //OpenMM_StringArray_destroy (pluginList);

    implicitSolventActive = usingImplicitSolvent ();

    tinkerForce = OpenMM_Vec3Array_create (atoms__.n);

    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: the OpenMM
    // System takes ownership of the force objects; don't delete them yourself.

    testName = "";
    system = OpenMM_System_create ();
    setupSystemParticles (system, log);

    countActiveForces = 0;
    if (potent__.use_bond)  countActiveForces++;
    if (potent__.use_angle)  countActiveForces++;
    if (potent__.use_strbnd)  countActiveForces++;
    if (potent__.use_urey)  countActiveForces++;
    if (potent__.use_opbend)  countActiveForces++;
    if (potent__.use_opdist)  countActiveForces++;
    if (potent__.use_improp)  countActiveForces++;
    if (potent__.use_imptor)  countActiveForces++;
    if (potent__.use_tors)  countActiveForces++;
    if (potent__.use_pitors)  countActiveForces++;
    if (potent__.use_strtor)  countActiveForces++;
    if (potent__.use_angtor)  countActiveForces++;
    if (potent__.use_tortor)  countActiveForces++;
    if (potent__.use_vdw)  countActiveForces++;
    if (potent__.use_charge)  countActiveForces++;
    if (potent__.use_chgdpl)  countActiveForces++;
    if (potent__.use_dipole)  countActiveForces++;
    if (potent__.use_mpole)  countActiveForces++;
    if (potent__.use_polar)  countActiveForces++;
    if (potent__.use_rxnfld)  countActiveForces++;
    if (potent__.use_metal)  countActiveForces++;
    if (potent__.use_geom)  countActiveForces++;
    if (potent__.use_solv)  countActiveForces++;

    if (log) {
        (void) fprintf (log, "\n Potential Terms Used in TINKER:\n" );
        (void) fprintf (log, "\n    Bond=    %d", abs(potent__.use_bond));
        (void) fprintf (log, "    Angle=   %d", abs(potent__.use_angle));
        (void) fprintf (log, "    StrBnd=  %d", abs(potent__.use_strbnd));
        (void) fprintf (log, "    Urey=    %d", abs(potent__.use_urey));
        (void) fprintf (log, "    AngAng=  %d", abs(potent__.use_angang));
        (void) fprintf (log, "\n    OPBend=  %d", abs(potent__.use_opbend));
        (void) fprintf (log, "    OPDist=  %d", abs(potent__.use_opdist));
        (void) fprintf (log, "    ImProp=  %d", abs(potent__.use_improp));
        (void) fprintf (log, "    ImpTor=  %d", abs(potent__.use_imptor));
        (void) fprintf (log, "    Tors=    %d", abs(potent__.use_tors));
        (void) fprintf (log, "\n    PiTors=  %d", abs(potent__.use_pitors));
        (void) fprintf (log, "    StrTor=  %d", abs(potent__.use_strtor));
        (void) fprintf (log, "    AngTor=  %d", abs(potent__.use_angtor));
        (void) fprintf (log, "    TorTor=  %d", abs(potent__.use_tortor));
        (void) fprintf (log, "    Vdw=     %d", abs(potent__.use_vdw));
        (void) fprintf (log, "\n    Charge=  %d", abs(potent__.use_charge));
        (void) fprintf (log, "    ChgDpl=  %d", abs(potent__.use_chgdpl));
        (void) fprintf (log, "    Dipole=  %d", abs(potent__.use_dipole));
        (void) fprintf (log, "    MPole=   %d", abs(potent__.use_mpole));
        (void) fprintf (log, "    Polar=   %d", abs(potent__.use_polar));
        (void) fprintf (log, "\n    RxnFld=  %d", abs(potent__.use_rxnfld));
        (void) fprintf (log, "    LigFld=  %d", abs(potent__.use_metal));
        (void) fprintf (log, "    Restrn=  %d", abs(potent__.use_geom));
        (void) fprintf (log, "    Solv=    %d", abs(potent__.use_solv));
        (void) fprintf (log, "    Extra=   %d\n", abs(potent__.use_extra));
    }

    if (countActiveForces > 1) {

        setupAmoebaBondForce (system, removeConstrainedCovalentIxns, log);
        setupAmoebaAngleForce (system, removeConstrainedCovalentIxns, log);
        setupAmoebaInPlaneAngleForce (system, removeConstrainedCovalentIxns,
                                      log);
        setupAmoebaTorsionForce (system, log);
        setupAmoebaPiTorsionForce (system, log);
        setupAmoebaStretchBendForce (system, removeConstrainedCovalentIxns,
                                     log);
        setupAmoebaOutOfPlaneBendForce (system, log);
        setupAmoebaTorsionTorsionForce (system, log);
        setupAmoebaVdwForce (system, log);
        setupAmoebaMultipoleForce (system, log);

        if (potent__.use_solv) {
            setupAmoebaWcaDispersionForce (system, log);
            setupAmoebaGeneralizedKirkwoodForce (system, 1, log);
        }

        if (potent__.use_urey) {
            setupAmoebaUreyBradleyForce (system,
                                         removeConstrainedCovalentIxns, log);
        } 

        loadTinkerForce (deriv__.desum, 0, tinkerForce);
        tinkerEnergy = *energi__.esum;
        testName = "AmoebaAllTest";

        if (log) {
            (void) fprintf (log,
                       "\n Potential Energy Components from TINKER:\n\n"
                       "    EB=  %15.7e   EA=  %15.7e   EBA= %15.7e\n"
                       "    EUB= %15.7e   EAA= %15.7e   EOPB=%15.7e\n"
                       "    EOPD=%15.7e   EID= %15.7e   EIT= %15.7e\n"
                       "    ET=  %15.7e   EPT= %15.7e   EBT= %15.7e\n"
                       "    EAT= %15.7e   ETT= %15.7e   EV=  %15.7e\n"
                       "    EC=  %15.7e   ECD= %15.7e   ED=  %15.7e\n"
                       "    EM=  %15.7e   EP=  %15.7e   ER=  %15.7e\n"
                       "    ES=  %15.7e   ELF= %15.7e   EG=  %15.7e\n"
                       "    EX=  %15.7e\n",
                       *energi__.eb,   *energi__.ea,  *energi__.eba,
                       *energi__.eub,  *energi__.eaa, *energi__.eopb,
                       *energi__.eopd, *energi__.eid, *energi__.eit,
                       *energi__.et,   *energi__.ept, *energi__.ebt,
                       *energi__.eat,  *energi__.ett, *energi__.ev,
                       *energi__.ec,   *energi__.ecd, *energi__.ed,
                       *energi__.em,   *energi__.ep,  *energi__.er,
                       *energi__.es,   *energi__.elf, *energi__.eg,
                       *energi__.ex );
            (void) fflush (log);
        }
    } else if (potent__.use_bond && !potent__.use_mpole) {

        setupAmoebaBondForce (system, removeConstrainedCovalentIxns, log);
        loadTinkerForce (deriv__.deb, 0, tinkerForce);
        tinkerEnergy = *energi__.eb;
        testName = "AmoebaHarmonicBondTest";

    } else if (potent__.use_angle) {
    
        // note TINKER angle = OpenMM (Angle + InPlaneAngle)

        setupAmoebaAngleForce (system, removeConstrainedCovalentIxns, log);
        setupAmoebaInPlaneAngleForce (system, removeConstrainedCovalentIxns,
                                      log);
        loadTinkerForce (deriv__.dea, 0, tinkerForce);
        tinkerEnergy = *energi__.ea;
        testName = "AmoebaHarmonicAngleTest";

    } else if (potent__.use_strbnd) {

        setupAmoebaStretchBendForce (system, removeConstrainedCovalentIxns,
                                     log);
        loadTinkerForce (deriv__.deba, 0, tinkerForce);
        tinkerEnergy = *energi__.eba;
        testName = "AmoebaStretchBendTest";

    } else if (potent__.use_urey) {

        setupAmoebaUreyBradleyForce (system, removeConstrainedCovalentIxns,
                                     log);
        loadTinkerForce (deriv__.deub, 0, tinkerForce);
        tinkerEnergy = *energi__.eub;
        testName = "AmoebaUreyBradleyForceTest";

    } else if (potent__.use_opbend) {

        setupAmoebaOutOfPlaneBendForce (system, log);
        loadTinkerForce (deriv__.deopb, 0, tinkerForce);
        tinkerEnergy = *energi__.eopb;
        testName = "AmoebaOutOfPlaneBendTest";

    } else if (potent__.use_tors) {

        setupAmoebaTorsionForce (system, log);
        loadTinkerForce (deriv__.det, 0, tinkerForce);
        tinkerEnergy = *energi__.et;
        testName = "AmoebaTorsionTest";

    } else if (potent__.use_pitors) {

        setupAmoebaPiTorsionForce (system, log);
        loadTinkerForce (deriv__.dept, 0, tinkerForce);
        tinkerEnergy = *energi__.ept;
        testName = "AmoebaPiTorsionTest";

    } else if (potent__.use_tortor) {

        setupAmoebaTorsionTorsionForce (system, log);
        loadTinkerForce (deriv__.dett, 0, tinkerForce);
        tinkerEnergy = *energi__.ett;
        testName = "AmoebaTorsionTorsionTest";

    } else if (potent__.use_vdw) {

        setupAmoebaVdwForce (system, log);
        loadTinkerForce (deriv__.dev, 0, tinkerForce);
        tinkerEnergy = *energi__.ev;
        if (limits__.use_vlist) {
            testName = "AmoebaVdwTest";
        } else {
            testName = "AmoebaVdwNoCutoffTest";
        }
        if (log) {
            (void) fprintf (log, "\n Use Vdw Neighbor List=%d Vdw Cutoff=%12.3f\n",
                            limits__.use_vlist, limits__.vdwcut);
        }

    } else if (potent__.use_mpole && !potent__.use_solv) {

        if (log) {
            (void) fprintf (log, "implicitSolventActive=%d\n",
                            implicitSolventActive);
        }
        if (implicitSolventActive == 0) {
            setupAmoebaBondForce (system, removeConstrainedCovalentIxns, log);
            loadTinkerForce (deriv__.deb, 0, tinkerForce);
            tinkerEnergy = *energi__.eb;
        } else {
            zeroVec3Force (tinkerForce);
            tinkerEnergy = 0.0;
        }
        setupAmoebaMultipoleForce (system, log);
        loadTinkerForce (deriv__.dem, 1, tinkerForce);
        loadTinkerForce (deriv__.dep, 1, tinkerForce);
        tinkerEnergy += *energi__.em + *energi__.ep;

        if (strncasecmp (polpot__.poltyp, "DIRECT", 6) == 0) { 
            testName = "AmoebaMultipoleDirectTest";
        } else {
            testName = "AmoebaMultipoleMutualTest";
        }

    } else if (potent__.use_solv && implicitSolventActive > 0 &&
               !potent__.use_mpole) {

        // to get only TINKER WCA, zero deriv__.des, then call ewca1 

        zeroTinkerForce (deriv__.des);
        ewca1_ (&tinkerEnergy);
        loadTinkerForce (deriv__.des, 0, tinkerForce);
        setupAmoebaWcaDispersionForce (system, log);
        testName = "AmoebaWcaDispersionTest";

    } else if (implicitSolventActive > 0 && potent__.use_mpole) {

        // generalized Kirkwood; OpenMM should have cavity term turned off

        double ecav, edisp;
        if (log) { 
            char buffer[128];
            setNullTerminator (solute__.borntyp, 8, buffer);
            (void) fprintf (log, "Born Radius type=%s ", buffer);
            setNullTerminator (solute__.solvtyp, 8, buffer);
            (void) fprintf (log, "Solvation Type=%s\n", buffer);
        }
        setupAmoebaMultipoleForce (system, log);
        setupAmoebaGeneralizedKirkwoodForce (system, 1, log);
        zeroTinkerForce (deriv__.des);
        setTinker1DArray (atoms__.n, solute__.drb, 0.0);
        setTinker1DArray (atoms__.n, solute__.drbp, 0.0);

        *energi__.es = 0.0;
        born_ ();
        empole1_ ();
        egk1_ ();

        loadTinkerForce (deriv__.des, 0, tinkerForce);
        loadTinkerForce (deriv__.dem, 1, tinkerForce);
        loadTinkerForce (deriv__.dep, 1, tinkerForce);
        tinkerEnergy = *energi__.es + *energi__.em + *energi__.ep;

        enp1_ (&ecav,&edisp);
        if (log) { 
            (void) fprintf (log, "Energies total=%15.7e Gk=%15.7e (cavity=%15.7e dispersion=%15.7e) Multipole=%15.7e\n",
                            tinkerEnergy, *energi__.es, ecav, edisp,
                            *energi__.em + *energi__.ep );

            //for (ii = 0; ii < atoms__.n; ii++) {
            //    fprintf (stderr, "Fs %5d [%15.7e %15.7e %15.7e] [%15.7e %15.7e%15.7e] [%15.7e %15.7e %15.7e]\n",
            //             ii, *(deriv__.des + 3*ii), *(deriv__.des + 3*ii+1),
            //             *(deriv__.des + 3*ii+2), *(deriv__.dem + 3*ii),
            //             *(deriv__.dem + 3*ii+1), *(deriv__.dem + 3*ii+2),
            //             *(deriv__.dep + 3*ii), *(deriv__.dep + 3*ii+1),
            //             *(deriv__.dep + 3*ii+2));
            //}

        }
        if( strncasecmp (polpot__.poltyp, "DIRECT", 6) == 0) {
            testName = "AmoebaKirkwoodDirectTest";
        } else {
            testName = "AmoebaKirkwoodMutualTest";
        }
    }

    if (log) {
        if (!testName.empty()) {
            (void) fprintf (log, "\n Test Option:  %s\n", testName.c_str());
            (void) fflush (NULL);
        } else {
            (void) fprintf (log, "\n Test Option not Recognized; Exiting\n");
            (void) fflush (log);
            exit (-1);
        }
    }

    initialPosInNm = OpenMM_Vec3Array_create (0);
    setupPositions (initialPosInNm, log);

    integrator = (OpenMM_Integrator*) OpenMM_VerletIntegrator_create (0.001);

    //platform = getReferencePlatform (log);
    platform = getCUDAPlatform (log);
    context = OpenMM_Context_create_2 (system, integrator, platform);

    OpenMM_Context_setPositions (context, initialPosInNm);

    infoMask = OpenMM_State_Positions;
    infoMask += OpenMM_State_Forces;
    infoMask += OpenMM_State_Energy;

    state = OpenMM_Context_getState (context, infoMask, 0);

    openMMForces = OpenMM_State_getForces (state);
    openMMPotentialEnergy = OpenMM_State_getPotentialEnergy(state)
                                / OpenMM_KJPerKcal;

    conversion = OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;

    //if (log) {
    //    (void) fprintf (log, "\n Unit Conversion (kJ/nm->kcal/Ang):%15.5f\n",
    //                    conversion );
    //    printDefaultPeriodicBoxVectors (system, log);
    //    (void) fflush (log);
    //}

    conversion = -1.0 / conversion;

    // find any differences between TINKER and OpenMM forces

    maxFDelta = 0.0;
    maxFDeltaIndex = -1;
    minDot = 1.0;
    avgDot = 0.0;
    minDotIndex = -1;

    if (log) {

        (void) fprintf (log, "\n TINKER vs OpenMM Energy Values:\n");
        (void) fprintf (log, "\n TINKER Potential Energy   %18.4f",
                        tinkerEnergy);
        (void) fprintf (log, "\n OpenMM Potential Energy   %18.4f\n",
                        openMMPotentialEnergy);
        maxEDelta = fabs(tinkerEnergy - openMMPotentialEnergy);

        maxFDelta = 0.0;
        for (ii = 0; ii < atoms__.n; ii++) {
            double relxNrm;
            double dot;
            OpenMM_Vec3 force;
            const OpenMM_Vec3* tinkerF;
            force.x = conversion*(*OpenMM_Vec3Array_get(openMMForces, ii)).x;
            force.y = conversion*(*OpenMM_Vec3Array_get(openMMForces, ii)).y;
            force.z = conversion*(*OpenMM_Vec3Array_get(openMMForces, ii)).z;
            tinkerF = OpenMM_Vec3Array_get(tinkerForce, ii);

            tinkerNorm = sqrt(tinkerF->x*tinkerF->x + tinkerF->y*tinkerF->y
                                  + tinkerF->z*tinkerF->z);
            openMMNorm = sqrt(force.x*force.x + force.y*force.y
                                  + force.z*force.z );

            delta = sqrt((tinkerF->x-force.x)*(tinkerF->x-force.x)
                       + (tinkerF->y-force.y)*(tinkerF->y-force.y)
                       + (tinkerF->z-force.z)*(tinkerF->z-force.z));
            dot = ((tinkerNorm > 0.0) && (openMMNorm > 0.0)) ?
                   (tinkerF->x*force.x + tinkerF->y*force.y
                        + tinkerF->z*force.z) / (tinkerNorm*openMMNorm) : 0.0;

            if (delta > maxFDelta) {
                maxFDelta = delta;
                maxFDeltaIndex = ii + 1;
            }
            if (dot < minDot) {
                minDot = dot;
                minDotIndex = ii + 1;
            }
            avgDot += dot;

            //(void) fprintf (stderr, "%6d   [%15.7e %15.7e %15.7e]",
            //                ii+1, tinkerF->x, tinkerF->y, tinkerF->z);
            //(void) fprintf (stderr, "\n         [%15.7e %15.7e %15.7e]\n",
            //                force.x, force.y, force.z); 
        }
    
        if (atoms__.n) {
            avgDot /= (double)(atoms__.n); 
        }

        (void) fprintf (log, "\n Summary of TINKER vs OpenMM Comparison:\n");
        (void) fprintf (log, "\n EnergyDelta                  %15.7e\n",
                        maxEDelta);
        (void) fprintf (log, " MaxForceDelta at   %7d   %15.7e\n",
                        maxFDeltaIndex, maxFDelta);
        (void) fprintf (log, " MinDotProduct at   %7d   %15.7e\n",
                        minDotIndex, minDot);
        (void) fprintf (log, " AvgDotProduct                %15.7e\n",
                        avgDot);
    }

    OpenMM_Vec3Array_destroy (tinkerForce);
    OpenMM_State_destroy (state);
    OpenMM_Context_destroy (context);
    OpenMM_Integrator_destroy (integrator);
    OpenMM_System_destroy (system);
}

void openmm_serialize_(void** omm) {
    OpenMMData* ommHandle = (OpenMMData*) (*omm);

    char *stuff;

    FILE *f = fopen("system.xml", "w");

    stuff = OpenMM_XmlSerializer_serializeSystem(ommHandle->system);

    fputs(stuff, f);
    fputs("\n", f);

    // We need to free the serialized content and close the file
    free(stuff);
    fclose(f);
}
}
