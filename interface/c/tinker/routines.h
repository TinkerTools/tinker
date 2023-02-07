#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#if defined(TINKER_GFORTRAN) && (__GNUC__ <= 7)
// https://gcc.gnu.org/onlinedocs/gfortran/Argument-passing-conventions.html
typedef int tinker_fchar_len_t;
#else
#include <string.h>
typedef size_t tinker_fchar_len_t;
#endif

typedef struct tinker_fchars_st { char* string; tinker_fchar_len_t capacity; } tinker_fchars;

void tinkerFortranRuntimeBegin(int, char**);
void tinkerFortranRuntimeEnd();

// active.f
void active_();
#define tinker_f_active active_

// alterchg.f
void alterchg_();
#define tinker_f_alterchg alterchg_
void bndchg_(double* pdelta);
#define tinker_f_bndchg bndchg_
void angchg_(double* pdelta);
#define tinker_f_angchg angchg_

// alterpol.f
void alterpol_();
#define tinker_f_alterpol alterpol_
void altpol0a_();
#define tinker_f_altpol0a altpol0a_
void altpol0b_();
#define tinker_f_altpol0b altpol0b_
void rotexpl_(double* r, double* xr, double* yr, double* zr, double* p33i, double* p33k, double* ks2i, double* ks2k);
#define tinker_f_rotexpl rotexpl_

// analysis.f
void analysis_(double* energy);
#define tinker_f_analysis analysis_

// angles.f
void angles_();
#define tinker_f_angles angles_

// attach.f
void attach_();
#define tinker_f_attach attach_

// baoab.f
void baoab_(int* istep, double* dt);
#define tinker_f_baoab baoab_
void oprep_(double* dt, double* vfric, double* vrand);
#define tinker_f_oprep oprep_

// basefile.f
void basefile_(char* string, tinker_fchar_len_t string_cap);
inline void tinker_f_basefile(tinker_fchars string) {
    return basefile_(string.string, string.capacity);
}

// beeman.f
void beeman_(int* istep, double* dt);
#define tinker_f_beeman beeman_

// bicubic.f
void bcuint_(double* y, double* y1, double* y2, double* y12, double* x1l, double* x1u, double* x2l, double* x2u, double* x1, double* x2, double* ansy);
#define tinker_f_bcuint bcuint_
void bcuint1_(double* y, double* y1, double* y2, double* y12, double* x1l, double* x1u, double* x2l, double* x2u, double* x1, double* x2, double* ansy, double* ansy1, double* ansy2);
#define tinker_f_bcuint1 bcuint1_
void bcuint2_(double* y, double* y1, double* y2, double* y12, double* x1l, double* x1u, double* x2l, double* x2u, double* x1, double* x2, double* ansy, double* ansy1, double* ansy2, double* ansy12, double* ansy11, double* ansy22);
#define tinker_f_bcuint2 bcuint2_
void bcucof_(double* y, double* y1, double* y2, double* y12, double* d1, double* d2, double* c);
#define tinker_f_bcucof bcucof_

// bitors.f
void bitors_();
#define tinker_f_bitors bitors_

// bonds.f
void bonds_();
#define tinker_f_bonds bonds_

// born.f
void born_();
#define tinker_f_born born_
void born1_();
#define tinker_f_born1 born1_

// bounds.f
void bounds_();
#define tinker_f_bounds bounds_

// bussi.f
void bussi_(int* istep, double* dt);
#define tinker_f_bussi bussi_

// calendar.f
void calendar_(int* year, int* month, int* day, int* hour, int* minute, int* second);
#define tinker_f_calendar calendar_

// center.f
void center_(int* n1, double* x1, double* y1, double* z1, int* n2, double* x2, double* y2, double* z2, double* xmid, double* ymid, double* zmid);
#define tinker_f_center center_

// chkpole.f
void chkpole_();
#define tinker_f_chkpole chkpole_

// chkring.f
void chkring_(int* iring, int* ia, int* ib, int* ic, int* id);
#define tinker_f_chkring chkring_

// chkxyz.f
void chkxyz_(int* clash);
#define tinker_f_chkxyz chkxyz_

// cholesky.f
void cholesky_(int* nvar, double* a, double* b);
#define tinker_f_cholesky cholesky_

// clock.f
void settime_();
#define tinker_f_settime settime_
void gettime_(double* wall, double* cpu);
#define tinker_f_gettime gettime_

// cluster.f
void cluster_();
#define tinker_f_cluster cluster_

// column.f
void column_(int* nvar, int* hinit, int* hstop, int* hindex, int* cinit, int* cstop, int* cindex, int* cvalue);
#define tinker_f_column column_

// command.f
void command_();
#define tinker_f_command command_

// connect.f
void connect_();
#define tinker_f_connect connect_

// connolly.f
void connolly_(double* volume, double* area, double* radius, double* probe, double* exclude);
#define tinker_f_connolly connolly_
void wiggle_();
#define tinker_f_wiggle wiggle_
void nearby_();
#define tinker_f_nearby nearby_
void torus_();
#define tinker_f_torus torus_
void place_();
#define tinker_f_place place_
void inedge_(int* ien, int* itt);
#define tinker_f_inedge inedge_
void compress_();
#define tinker_f_compress compress_
void saddles_();
#define tinker_f_saddles saddles_
void gettor_(int* ia, int* ja, int* ttok, double* torcen, double* torad, double* torax);
#define tinker_f_gettor gettor_
void getprb_(int* ia, int* ja, int* ka, int* prbok, int* tb, double* bijk, double* hijk, double* uijk);
#define tinker_f_getprb getprb_
void ipedge_(int* ieq, int* ia);
#define tinker_f_ipedge ipedge_
void contact_();
#define tinker_f_contact contact_
void vam_(double* volume, double* area);
#define tinker_f_vam vam_
double depth_(int* ip, double* alt);
#define tinker_f_depth depth_
void measpm_(int* ifn, double* prism);
#define tinker_f_measpm measpm_
void measfq_(int* ifq, double* areaq, double* volq);
#define tinker_f_measfq measfq_
void measfs_(int* ifs, double* areas, double* vols, double* areasp, double* volsp);
#define tinker_f_measfs measfs_
void measfn_(int* ifn, double* arean, double* voln);
#define tinker_f_measfn measfn_
void projct_(double* pnt, double* unvect, int* icy, int* ia, double* spv, int* nedge, int* fail);
#define tinker_f_projct projct_
int ptincy_(double* pnt, double* unvect, int* icy);
#define tinker_f_ptincy ptincy_
void equclc_(double* spv, int* nedge, double* equ);
#define tinker_f_equclc equclc_
double rotang_(double* equ, int* nedge, double* unvect);
#define tinker_f_rotang rotang_
void vcross_(double* x, double* y, double* z);
#define tinker_f_vcross vcross_
double dot_(double* x, double* y);
#define tinker_f_dot dot_
double anorm_(double* x);
#define tinker_f_anorm anorm_
void vnorm_(double* x, double* xn);
#define tinker_f_vnorm vnorm_
double dist2_(double* x, double* y);
#define tinker_f_dist2 dist2_
double triple_(double* x, double* y, double* z);
#define tinker_f_triple triple_
double vecang_(double* v1, double* v2, double* axis, double* hand);
#define tinker_f_vecang vecang_
void cirpln_(double* circen, double* cirrad, double* cirvec, double* plncen, double* plnvec, int* cinsp, int* cintp, double* xpnt1, double* xpnt2);
#define tinker_f_cirpln cirpln_
void gendot_(int* ndots, double* dots, double* radius, double* xcenter, double* ycenter, double* zcenter);
#define tinker_f_gendot gendot_
void cerror_(char* string, tinker_fchar_len_t string_cap);
inline void tinker_f_cerror(tinker_fchars string) {
    return cerror_(string.string, string.capacity);
}

// control.f
void control_();
#define tinker_f_control control_

// cspline.f
void cspline_(int* n, double* xn, double* fn, double* b, double* c, double* d, double* h, double* du, double* dm, double* rc, double* rs);
#define tinker_f_cspline cspline_
void cytsy_(int* n, double* dm, double* du, double* cr, double* rs, double* x, int* iflag);
#define tinker_f_cytsy cytsy_
void cytsyp_(int* n, double* dm, double* du, double* cr, int* iflag);
#define tinker_f_cytsyp cytsyp_
void cytsys_(int* n, double* dm, double* du, double* cr, double* rs, double* x);
#define tinker_f_cytsys cytsys_

// cutoffs.f
void cutoffs_();
#define tinker_f_cutoffs cutoffs_

// damping.f
void dampewald_(int* rorder, double* r, double* r2, double* scale, double* dmpe);
#define tinker_f_dampewald dampewald_
void dampthole_(int* i, int* k, int* rorder, double* r, double* dmpik);
#define tinker_f_dampthole dampthole_
void damptholed_(int* i, int* k, int* rorder, double* r, double* dmpik);
#define tinker_f_damptholed damptholed_
void damppole_(double* r, int* rorder, double* alphai, double* alphak, double* dmpi, double* dmpk, double* dmpik);
#define tinker_f_damppole damppole_
void dampdir_(double* r, double* alphai, double* alphak, double* dmpi, double* dmpk);
#define tinker_f_dampdir dampdir_
void dampmut_(double* r, double* alphai, double* alphak, double* dmpik);
#define tinker_f_dampmut dampmut_
void damppot_(double* r, double* alphak, double* dmpk);
#define tinker_f_damppot damppot_
void damprep_(double* r, double* r2, double* rr1, double* rr3, double* rr5, double* rr7, double* rr9, double* rr11, int* rorder, double* dmpi, double* dmpk, double* dmpik);
#define tinker_f_damprep damprep_
void dampexpl_(double* r, double* preik, double* alphai, double* alphak, double* s2, double* ds2);
#define tinker_f_dampexpl dampexpl_

// dcflux.f
void dcflux_(double* pot, double* dcfx, double* dcfy, double* dcfz);
#define tinker_f_dcflux dcflux_

// deflate.f
void deflate_(int* n, int* nv, double* a, double* ev, double* vec);
#define tinker_f_deflate deflate_

// delete.f
void delete_(int* iatom);
#define tinker_f_delete delete_

// dexpol.f
void dexpol_();
#define tinker_f_dexpol dexpol_
void dexpol1a_();
#define tinker_f_dexpol1a dexpol1a_
void dexpol1b_();
#define tinker_f_dexpol1b dexpol1b_
void rotdexpl_(double* r, double* xr, double* yr, double* zr, double* ai, double* ak);
#define tinker_f_rotdexpl rotdexpl_

// diagq.f
void diagq_(int* n, int* nv, double* dd, double* ev, double* vec);
#define tinker_f_diagq diagq_

// diffeq.f
void diffeq_(int* nvar, double* y, double* x1, double* x2, double* eps, double* h1, double* hmin, int* nok, int* nbad, void (*gvalue)(double*, double*, double*));
#define tinker_f_diffeq diffeq_
void bsstep_(int* nvar, double* x, double* dydx, double* y, double* htry, double* eps, double* yscal, double* hdid, double* hnext, void (*gvalue)(double*, double*, double*));
#define tinker_f_bsstep bsstep_
void mmid_(int* nstep, double* htot, int* nvar, double* xs, double* dydx, double* y, double* yout, void (*gvalue)(double*, double*, double*));
#define tinker_f_mmid mmid_
void pzextr_(int* iest, int* nvar, double* xest, double* yest, double* yz, double* dy);
#define tinker_f_pzextr pzextr_
void gdastat_(int* nstep, double* beta, double* xx, char* status, tinker_fchar_len_t status_cap);
inline void tinker_f_gdastat(int* nstep, double* beta, double* xx, tinker_fchars status) {
    return gdastat_(nstep, beta, xx, status.string, status.capacity);
}

// eangang.f
void eangang_();
#define tinker_f_eangang eangang_

// eangang1.f
void eangang1_();
#define tinker_f_eangang1 eangang1_

// eangang2.f
void eangang2_(int* i);
#define tinker_f_eangang2 eangang2_
void eangang2a_(int* i, double* de);
#define tinker_f_eangang2a eangang2a_

// eangang3.f
void eangang3_();
#define tinker_f_eangang3 eangang3_

// eangle.f
void eangle_();
#define tinker_f_eangle eangle_

// eangle1.f
void eangle1_();
#define tinker_f_eangle1 eangle1_

// eangle2.f
void eangle2_(int* i);
#define tinker_f_eangle2 eangle2_
void eangle2a_(int* iatom);
#define tinker_f_eangle2a eangle2a_
void eangle2b_(int* i, double* de);
#define tinker_f_eangle2b eangle2b_

// eangle3.f
void eangle3_();
#define tinker_f_eangle3 eangle3_

// eangtor.f
void eangtor_();
#define tinker_f_eangtor eangtor_

// eangtor1.f
void eangtor1_();
#define tinker_f_eangtor1 eangtor1_

// eangtor2.f
void eangtor2_(int* i);
#define tinker_f_eangtor2 eangtor2_

// eangtor3.f
void eangtor3_();
#define tinker_f_eangtor3 eangtor3_

// ebond.f
void ebond_();
#define tinker_f_ebond ebond_

// ebond1.f
void ebond1_();
#define tinker_f_ebond1 ebond1_

// ebond2.f
void ebond2_(int* i);
#define tinker_f_ebond2 ebond2_

// ebond3.f
void ebond3_();
#define tinker_f_ebond3 ebond3_

// ebuck.f
void ebuck_();
#define tinker_f_ebuck ebuck_
void ebuck0a_();
#define tinker_f_ebuck0a ebuck0a_
void ebuck0b_();
#define tinker_f_ebuck0b ebuck0b_
void ebuck0c_();
#define tinker_f_ebuck0c ebuck0c_
void ebuck0d_();
#define tinker_f_ebuck0d ebuck0d_

// ebuck1.f
void ebuck1_();
#define tinker_f_ebuck1 ebuck1_
void ebuck1a_();
#define tinker_f_ebuck1a ebuck1a_
void ebuck1b_();
#define tinker_f_ebuck1b ebuck1b_
void ebuck1c_();
#define tinker_f_ebuck1c ebuck1c_
void ebuck1d_();
#define tinker_f_ebuck1d ebuck1d_

// ebuck2.f
void ebuck2_(int* i);
#define tinker_f_ebuck2 ebuck2_
void ebuck2a_(int* iatom);
#define tinker_f_ebuck2a ebuck2a_
void ebuck2b_(int* i);
#define tinker_f_ebuck2b ebuck2b_

// ebuck3.f
void ebuck3_();
#define tinker_f_ebuck3 ebuck3_
void ebuck3a_();
#define tinker_f_ebuck3a ebuck3a_
void ebuck3b_();
#define tinker_f_ebuck3b ebuck3b_
void ebuck3c_();
#define tinker_f_ebuck3c ebuck3c_
void ebuck3d_();
#define tinker_f_ebuck3d ebuck3d_

// echarge.f
void echarge_();
#define tinker_f_echarge echarge_
void echarge0a_();
#define tinker_f_echarge0a echarge0a_
void echarge0b_();
#define tinker_f_echarge0b echarge0b_
void echarge0c_();
#define tinker_f_echarge0c echarge0c_
void echarge0d_();
#define tinker_f_echarge0d echarge0d_
void echarge0e_();
#define tinker_f_echarge0e echarge0e_
void echarge0f_();
#define tinker_f_echarge0f echarge0f_
void echarge0g_();
#define tinker_f_echarge0g echarge0g_
void ecrecip_();
#define tinker_f_ecrecip ecrecip_

// echarge1.f
void echarge1_();
#define tinker_f_echarge1 echarge1_
void echarge1a_();
#define tinker_f_echarge1a echarge1a_
void echarge1b_();
#define tinker_f_echarge1b echarge1b_
void echarge1c_();
#define tinker_f_echarge1c echarge1c_
void echarge1d_();
#define tinker_f_echarge1d echarge1d_
void echarge1e_();
#define tinker_f_echarge1e echarge1e_
void echarge1f_();
#define tinker_f_echarge1f echarge1f_
void echarge1g_();
#define tinker_f_echarge1g echarge1g_
void ecrecip1_();
#define tinker_f_ecrecip1 ecrecip1_

// echarge2.f
void echarge2_(int* i);
#define tinker_f_echarge2 echarge2_
void echarge2a_(int* i);
#define tinker_f_echarge2a echarge2a_
void echarge2b_(int* i);
#define tinker_f_echarge2b echarge2b_
void echarge2c_(int* i);
#define tinker_f_echarge2c echarge2c_
void echarge2r_();
#define tinker_f_echarge2r echarge2r_
void echarge2d_(int* i);
#define tinker_f_echarge2d echarge2d_
void echarge2e_(int* i);
#define tinker_f_echarge2e echarge2e_
void echarge2f_(int* i);
#define tinker_f_echarge2f echarge2f_

// echarge3.f
void echarge3_();
#define tinker_f_echarge3 echarge3_
void echarge3a_();
#define tinker_f_echarge3a echarge3a_
void echarge3b_();
#define tinker_f_echarge3b echarge3b_
void echarge3c_();
#define tinker_f_echarge3c echarge3c_
void echarge3d_();
#define tinker_f_echarge3d echarge3d_
void echarge3e_();
#define tinker_f_echarge3e echarge3e_
void echarge3f_();
#define tinker_f_echarge3f echarge3f_
void echarge3g_();
#define tinker_f_echarge3g echarge3g_
void ecrecip3_();
#define tinker_f_ecrecip3 ecrecip3_

// echgdpl.f
void echgdpl_();
#define tinker_f_echgdpl echgdpl_

// echgdpl1.f
void echgdpl1_();
#define tinker_f_echgdpl1 echgdpl1_

// echgdpl2.f
void echgdpl2_(int* i);
#define tinker_f_echgdpl2 echgdpl2_

// echgdpl3.f
void echgdpl3_();
#define tinker_f_echgdpl3 echgdpl3_

// echgtrn.f
void echgtrn_();
#define tinker_f_echgtrn echgtrn_
void echgtrn0a_();
#define tinker_f_echgtrn0a echgtrn0a_
void echgtrn0b_();
#define tinker_f_echgtrn0b echgtrn0b_
void echgtrn0c_();
#define tinker_f_echgtrn0c echgtrn0c_

// echgtrn1.f
void echgtrn1_();
#define tinker_f_echgtrn1 echgtrn1_
void echgtrn1a_();
#define tinker_f_echgtrn1a echgtrn1a_
void echgtrn1b_();
#define tinker_f_echgtrn1b echgtrn1b_

// echgtrn2.f
void echgtrn2_(int* iatom);
#define tinker_f_echgtrn2 echgtrn2_

// echgtrn3.f
void echgtrn3_();
#define tinker_f_echgtrn3 echgtrn3_
void echgtrn3a_();
#define tinker_f_echgtrn3a echgtrn3a_
void echgtrn3b_();
#define tinker_f_echgtrn3b echgtrn3b_
void echgtrn3c_();
#define tinker_f_echgtrn3c echgtrn3c_

// edipole.f
void edipole_();
#define tinker_f_edipole edipole_

// edipole1.f
void edipole1_();
#define tinker_f_edipole1 edipole1_

// edipole2.f
void edipole2_(int* i);
#define tinker_f_edipole2 edipole2_

// edipole3.f
void edipole3_();
#define tinker_f_edipole3 edipole3_

// edisp.f
void edisp_();
#define tinker_f_edisp edisp_
void edisp0a_();
#define tinker_f_edisp0a edisp0a_
void edisp0b_();
#define tinker_f_edisp0b edisp0b_
void edisp0c_();
#define tinker_f_edisp0c edisp0c_
void edreal0c_();
#define tinker_f_edreal0c edreal0c_
void edisp0d_();
#define tinker_f_edisp0d edisp0d_
void edreal0d_();
#define tinker_f_edreal0d edreal0d_
void edrecip_();
#define tinker_f_edrecip edrecip_

// edisp1.f
void edisp1_();
#define tinker_f_edisp1 edisp1_
void edisp1a_();
#define tinker_f_edisp1a edisp1a_
void edisp1b_();
#define tinker_f_edisp1b edisp1b_
void edisp1c_();
#define tinker_f_edisp1c edisp1c_
void edreal1c_();
#define tinker_f_edreal1c edreal1c_
void edisp1d_();
#define tinker_f_edisp1d edisp1d_
void edreal1d_();
#define tinker_f_edreal1d edreal1d_
void edrecip1_();
#define tinker_f_edrecip1 edrecip1_

// edisp2.f
void edisp2_(int* iatom);
#define tinker_f_edisp2 edisp2_

// edisp3.f
void edisp3_();
#define tinker_f_edisp3 edisp3_
void edisp3a_();
#define tinker_f_edisp3a edisp3a_
void edisp3b_();
#define tinker_f_edisp3b edisp3b_
void edisp3c_();
#define tinker_f_edisp3c edisp3c_
void edreal3c_();
#define tinker_f_edreal3c edreal3c_
void edisp3d_();
#define tinker_f_edisp3d edisp3d_
void edreal3d_();
#define tinker_f_edreal3d edreal3d_

// egauss.f
void egauss_();
#define tinker_f_egauss egauss_
void egauss0a_();
#define tinker_f_egauss0a egauss0a_
void egauss0b_();
#define tinker_f_egauss0b egauss0b_
void egauss0c_();
#define tinker_f_egauss0c egauss0c_
void egauss0d_();
#define tinker_f_egauss0d egauss0d_

// egauss1.f
void egauss1_();
#define tinker_f_egauss1 egauss1_
void egauss1a_();
#define tinker_f_egauss1a egauss1a_
void egauss1b_();
#define tinker_f_egauss1b egauss1b_
void egauss1c_();
#define tinker_f_egauss1c egauss1c_
void egauss1d_();
#define tinker_f_egauss1d egauss1d_

// egauss2.f
void egauss2_(int* i);
#define tinker_f_egauss2 egauss2_
void egauss2a_(int* iatom);
#define tinker_f_egauss2a egauss2a_
void egauss2b_(int* iatom);
#define tinker_f_egauss2b egauss2b_

// egauss3.f
void egauss3_();
#define tinker_f_egauss3 egauss3_
void egauss3a_();
#define tinker_f_egauss3a egauss3a_
void egauss3b_();
#define tinker_f_egauss3b egauss3b_
void egauss3c_();
#define tinker_f_egauss3c egauss3c_
void egauss3d_();
#define tinker_f_egauss3d egauss3d_

// egeom.f
void egeom_();
#define tinker_f_egeom egeom_

// egeom1.f
void egeom1_();
#define tinker_f_egeom1 egeom1_

// egeom2.f
void egeom2_(int* i);
#define tinker_f_egeom2 egeom2_

// egeom3.f
void egeom3_();
#define tinker_f_egeom3 egeom3_

// ehal.f
void ehal_();
#define tinker_f_ehal ehal_
void ehal0a_();
#define tinker_f_ehal0a ehal0a_
void ehal0b_();
#define tinker_f_ehal0b ehal0b_
void ehal0c_();
#define tinker_f_ehal0c ehal0c_

// ehal1.f
void ehal1_();
#define tinker_f_ehal1 ehal1_
void ehal1a_();
#define tinker_f_ehal1a ehal1a_
void ehal1b_();
#define tinker_f_ehal1b ehal1b_
void ehal1c_();
#define tinker_f_ehal1c ehal1c_

// ehal2.f
void ehal2_(int* iatom);
#define tinker_f_ehal2 ehal2_

// ehal3.f
void ehal3_();
#define tinker_f_ehal3 ehal3_
void ehal3a_();
#define tinker_f_ehal3a ehal3a_
void ehal3b_();
#define tinker_f_ehal3b ehal3b_
void ehal3c_();
#define tinker_f_ehal3c ehal3c_

// eimprop.f
void eimprop_();
#define tinker_f_eimprop eimprop_

// eimprop1.f
void eimprop1_();
#define tinker_f_eimprop1 eimprop1_

// eimprop2.f
void eimprop2_(int* i);
#define tinker_f_eimprop2 eimprop2_

// eimprop3.f
void eimprop3_();
#define tinker_f_eimprop3 eimprop3_

// eimptor.f
void eimptor_();
#define tinker_f_eimptor eimptor_

// eimptor1.f
void eimptor1_();
#define tinker_f_eimptor1 eimptor1_

// eimptor2.f
void eimptor2_(int* i);
#define tinker_f_eimptor2 eimptor2_

// eimptor3.f
void eimptor3_();
#define tinker_f_eimptor3 eimptor3_

// elj.f
void elj_();
#define tinker_f_elj elj_
void elj0a_();
#define tinker_f_elj0a elj0a_
void elj0b_();
#define tinker_f_elj0b elj0b_
void elj0c_();
#define tinker_f_elj0c elj0c_
void elj0d_();
#define tinker_f_elj0d elj0d_
void elj0e_();
#define tinker_f_elj0e elj0e_

// elj1.f
void elj1_();
#define tinker_f_elj1 elj1_
void elj1a_();
#define tinker_f_elj1a elj1a_
void elj1b_();
#define tinker_f_elj1b elj1b_
void elj1c_();
#define tinker_f_elj1c elj1c_
void elj1d_();
#define tinker_f_elj1d elj1d_
void elj1e_();
#define tinker_f_elj1e elj1e_

// elj2.f
void elj2_(int* i);
#define tinker_f_elj2 elj2_
void elj2a_(int* iatom);
#define tinker_f_elj2a elj2a_
void elj2b_(int* i);
#define tinker_f_elj2b elj2b_
void elj2c_(int* iatom);
#define tinker_f_elj2c elj2c_

// elj3.f
void elj3_();
#define tinker_f_elj3 elj3_
void elj3a_();
#define tinker_f_elj3a elj3a_
void elj3b_();
#define tinker_f_elj3b elj3b_
void elj3c_();
#define tinker_f_elj3c elj3c_
void elj3d_();
#define tinker_f_elj3d elj3d_
void elj3e_();
#define tinker_f_elj3e elj3e_

// embed.f
void embed_();
#define tinker_f_embed embed_
void kchiral_();
#define tinker_f_kchiral kchiral_
void triangle_();
#define tinker_f_triangle triangle_
void geodesic_();
#define tinker_f_geodesic geodesic_
void minpath_(int* root, double* upper, double* lower, int* start, int* stop, int* list);
#define tinker_f_minpath minpath_
void trifix_(int* p, int* q);
#define tinker_f_trifix trifix_
void grafic_(int* n, double* a, char* title, tinker_fchar_len_t title_cap);
inline void tinker_f_grafic(int* n, double* a, tinker_fchars title) {
    return grafic_(n, a, title.string, title.capacity);
}
void dstmat_(double* dmx);
#define tinker_f_dstmat dstmat_
void metric_(double* gmx, int* nneg);
#define tinker_f_metric metric_
void eigen_(double* evl, double* evc, double* gmx, int* valid);
#define tinker_f_eigen eigen_
void coords_(double* evl, double* evc);
#define tinker_f_coords coords_
void chksize_();
#define tinker_f_chksize chksize_
void majorize_(double* dmx);
#define tinker_f_majorize majorize_
void refine_(char* mode, double* fctval, double* grdmin, tinker_fchar_len_t mode_cap);
inline void tinker_f_refine(tinker_fchars mode, double* fctval, double* grdmin) {
    return refine_(mode.string, fctval, grdmin, mode.capacity);
}
void explore_(char* mode, int* nstep, double* dt, double* mass, double* temp_start, double* temp_stop, double* v, double* a, tinker_fchar_len_t mode_cap);
inline void tinker_f_explore(tinker_fchars mode, int* nstep, double* dt, double* mass, double* temp_start, double* temp_stop, double* v, double* a) {
    return explore_(mode.string, nstep, dt, mass, temp_start, temp_stop, v, a, mode.capacity);
}
void fracdist_(char* title, tinker_fchar_len_t title_cap);
inline void tinker_f_fracdist(tinker_fchars title) {
    return fracdist_(title.string, title.capacity);
}
void rmserror_(char* title, tinker_fchar_len_t title_cap);
inline void tinker_f_rmserror(tinker_fchars title) {
    return rmserror_(title.string, title.capacity);
}
void dmdump_(double* dmd);
#define tinker_f_dmdump dmdump_
double initerr_(double* xx, double* g);
#define tinker_f_initerr initerr_
double miderr_(double* xx, double* g);
#define tinker_f_miderr miderr_
double toterr_(double* xx, double* g);
#define tinker_f_toterr toterr_
double bnderr_(double* derivs);
#define tinker_f_bnderr bnderr_
double vdwerr_(double* derivs);
#define tinker_f_vdwerr vdwerr_
double locerr_(double* derivs);
#define tinker_f_locerr locerr_
double chirer_(double* derivs);
#define tinker_f_chirer chirer_
double torser_(double* derivs);
#define tinker_f_torser torser_

// emetal.f
void emetal_();
#define tinker_f_emetal emetal_

// emetal1.f
void emetal1_();
#define tinker_f_emetal1 emetal1_

// emetal2.f
void emetal2_(int* i);
#define tinker_f_emetal2 emetal2_

// emetal3.f
void emetal3_();
#define tinker_f_emetal3 emetal3_

// emm3hb.f
void emm3hb_();
#define tinker_f_emm3hb emm3hb_
void emm3hb0a_();
#define tinker_f_emm3hb0a emm3hb0a_
void emm3hb0b_();
#define tinker_f_emm3hb0b emm3hb0b_
void emm3hb0c_();
#define tinker_f_emm3hb0c emm3hb0c_

// emm3hb1.f
void emm3hb1_();
#define tinker_f_emm3hb1 emm3hb1_
void emm3hb1a_();
#define tinker_f_emm3hb1a emm3hb1a_
void emm3hb1b_();
#define tinker_f_emm3hb1b emm3hb1b_
void emm3hb1c_();
#define tinker_f_emm3hb1c emm3hb1c_

// emm3hb2.f
void emm3hb2_(int* iatom);
#define tinker_f_emm3hb2 emm3hb2_

// emm3hb3.f
void emm3hb3_();
#define tinker_f_emm3hb3 emm3hb3_
void emm3hb3a_();
#define tinker_f_emm3hb3a emm3hb3a_
void emm3hb3b_();
#define tinker_f_emm3hb3b emm3hb3b_
void emm3hb3c_();
#define tinker_f_emm3hb3c emm3hb3c_

// empole.f
void empole_();
#define tinker_f_empole empole_
void empole0a_();
#define tinker_f_empole0a empole0a_
void empole0b_();
#define tinker_f_empole0b empole0b_
void empole0c_();
#define tinker_f_empole0c empole0c_
void emreal0c_();
#define tinker_f_emreal0c emreal0c_
void empole0d_();
#define tinker_f_empole0d empole0d_
void emreal0d_();
#define tinker_f_emreal0d emreal0d_
void emrecip_();
#define tinker_f_emrecip emrecip_

// empole1.f
void empole1_();
#define tinker_f_empole1 empole1_
void empole1a_();
#define tinker_f_empole1a empole1a_
void empole1b_();
#define tinker_f_empole1b empole1b_
void empole1c_();
#define tinker_f_empole1c empole1c_
void emreal1c_();
#define tinker_f_emreal1c emreal1c_
void empole1d_();
#define tinker_f_empole1d empole1d_
void emreal1d_();
#define tinker_f_emreal1d emreal1d_
void emrecip1_();
#define tinker_f_emrecip1 emrecip1_

// empole2.f
void empole2_(int* i);
#define tinker_f_empole2 empole2_
void empole2a_(int* nlist, int* list);
#define tinker_f_empole2a empole2a_

// empole3.f
void empole3_();
#define tinker_f_empole3 empole3_
void empole3a_();
#define tinker_f_empole3a empole3a_
void empole3b_();
#define tinker_f_empole3b empole3b_
void empole3c_();
#define tinker_f_empole3c empole3c_
void emreal3c_();
#define tinker_f_emreal3c emreal3c_
void empole3d_();
#define tinker_f_empole3d empole3d_
void emreal3d_();
#define tinker_f_emreal3d emreal3d_
void emrecip3_();
#define tinker_f_emrecip3 emrecip3_

// energy.f
double energy_();
#define tinker_f_energy energy_

// eopbend.f
void eopbend_();
#define tinker_f_eopbend eopbend_

// eopbend1.f
void eopbend1_();
#define tinker_f_eopbend1 eopbend1_

// eopbend2.f
void eopbend2_(int* i);
#define tinker_f_eopbend2 eopbend2_
void eopbend2a_(int* i, double* de);
#define tinker_f_eopbend2a eopbend2a_

// eopbend3.f
void eopbend3_();
#define tinker_f_eopbend3 eopbend3_

// eopdist.f
void eopdist_();
#define tinker_f_eopdist eopdist_

// eopdist1.f
void eopdist1_();
#define tinker_f_eopdist1 eopdist1_

// eopdist2.f
void eopdist2_(int* i);
#define tinker_f_eopdist2 eopdist2_

// eopdist3.f
void eopdist3_();
#define tinker_f_eopdist3 eopdist3_

// epitors.f
void epitors_();
#define tinker_f_epitors epitors_

// epitors1.f
void epitors1_();
#define tinker_f_epitors1 epitors1_

// epitors2.f
void epitors2_(int* i);
#define tinker_f_epitors2 epitors2_
void epitors2a_(int* i, double* de);
#define tinker_f_epitors2a epitors2a_

// epitors3.f
void epitors3_();
#define tinker_f_epitors3 epitors3_

// epolar.f
void epolar_();
#define tinker_f_epolar epolar_
void epolar0a_();
#define tinker_f_epolar0a epolar0a_
void epolar0b_();
#define tinker_f_epolar0b epolar0b_
void epolar0c_();
#define tinker_f_epolar0c epolar0c_
void epreal0c_();
#define tinker_f_epreal0c epreal0c_
void epolar0d_();
#define tinker_f_epolar0d epolar0d_
void epreal0d_();
#define tinker_f_epreal0d epreal0d_
void epolar0e_();
#define tinker_f_epolar0e epolar0e_
void eprecip_();
#define tinker_f_eprecip eprecip_

// epolar1.f
void epolar1_();
#define tinker_f_epolar1 epolar1_
void epolar1a_();
#define tinker_f_epolar1a epolar1a_
void epolar1b_();
#define tinker_f_epolar1b epolar1b_
void epolar1c_();
#define tinker_f_epolar1c epolar1c_
void epreal1c_();
#define tinker_f_epreal1c epreal1c_
void epolar1d_();
#define tinker_f_epolar1d epolar1d_
void epreal1d_();
#define tinker_f_epreal1d epreal1d_
void epolar1e_();
#define tinker_f_epolar1e epolar1e_
void eprecip1_();
#define tinker_f_eprecip1 eprecip1_

// epolar2.f
void epolar2_(int* i);
#define tinker_f_epolar2 epolar2_
void epolar2a_(int* nlist, int* list, int* reinduce);
#define tinker_f_epolar2a epolar2a_

// epolar3.f
void epolar3_();
#define tinker_f_epolar3 epolar3_
void epolar3a_();
#define tinker_f_epolar3a epolar3a_
void epolar3b_();
#define tinker_f_epolar3b epolar3b_
void epolar3c_();
#define tinker_f_epolar3c epolar3c_
void epreal3c_();
#define tinker_f_epreal3c epreal3c_
void epolar3d_();
#define tinker_f_epolar3d epolar3d_
void epreal3d_();
#define tinker_f_epreal3d epreal3d_
void epolar3e_();
#define tinker_f_epolar3e epolar3e_
void eprecip3_();
#define tinker_f_eprecip3 eprecip3_

// erepel.f
void erepel_();
#define tinker_f_erepel erepel_
void erepel0a_();
#define tinker_f_erepel0a erepel0a_
void erepel0b_();
#define tinker_f_erepel0b erepel0b_

// erepel1.f
void erepel1_();
#define tinker_f_erepel1 erepel1_
void erepel1a_();
#define tinker_f_erepel1a erepel1a_
void erepel1b_();
#define tinker_f_erepel1b erepel1b_

// erepel2.f
void erepel2_(int* i);
#define tinker_f_erepel2 erepel2_
void erepel2a_(int* nlist, int* list);
#define tinker_f_erepel2a erepel2a_

// erepel3.f
void erepel3_();
#define tinker_f_erepel3 erepel3_
void erepel3a_();
#define tinker_f_erepel3a erepel3a_
void erepel3b_();
#define tinker_f_erepel3b erepel3b_

// erf.f
double erf_(double* x);
#define tinker_f_erf erf_
double erfc_(double* x);
#define tinker_f_erfc erfc_
void erfcore_(double* arg, double* result, int* mode);
#define tinker_f_erfcore erfcore_
double erfinv_(double* x);
#define tinker_f_erfinv erfinv_

// erxnfld.f
void erxnfld_();
#define tinker_f_erxnfld erxnfld_
void erfik_(int* ii, int* kk, int* i, int* k, double* rpi, double* rpk, double* eik);
#define tinker_f_erfik erfik_
void rfindex_(int* n, int* m, int* ind_x, int* ind_y, int* ind_z, int* p_s, int* p_e);
#define tinker_f_rfindex rfindex_
void ijkpts_();
#define tinker_f_ijkpts ijkpts_
double d1d2_(int* n, double* x1, double* y1, double* z1, double* x2, double* y2, double* z2, double* d, double* r1sq, double* r2sq, int* i, int* j, int* k, int* s, int* t, int* u);
#define tinker_f_d1d2 d1d2_

// erxnfld1.f
void erxnfld1_();
#define tinker_f_erxnfld1 erxnfld1_

// erxnfld2.f
void erxnfld2_(int* i);
#define tinker_f_erxnfld2 erxnfld2_

// erxnfld3.f
void erxnfld3_();
#define tinker_f_erxnfld3 erxnfld3_

// esolv.f
void esolv_();
#define tinker_f_esolv esolv_
void egb0a_();
#define tinker_f_egb0a egb0a_
void egb0b_();
#define tinker_f_egb0b egb0b_
void egb0c_();
#define tinker_f_egb0c egb0c_
void egk_();
#define tinker_f_egk egk_
void egk0a_();
#define tinker_f_egk0a egk0a_
void epb_();
#define tinker_f_epb epb_
void ediff_();
#define tinker_f_ediff ediff_
void pbempole_();
#define tinker_f_pbempole pbempole_
void enp_(double* ecav, double* edisp);
#define tinker_f_enp enp_
void ewca_(double* edisp);
#define tinker_f_ewca ewca_
void ewcax_(double* edisp);
#define tinker_f_ewcax ewcax_
void ehpmf_(double* ehp);
#define tinker_f_ehpmf ehpmf_

// esolv1.f
void esolv1_();
#define tinker_f_esolv1 esolv1_
void egb1a_();
#define tinker_f_egb1a egb1a_
void egb1b_();
#define tinker_f_egb1b egb1b_
void egb1c_();
#define tinker_f_egb1c egb1c_
void egk1_();
#define tinker_f_egk1 egk1_
void egk1a_();
#define tinker_f_egk1a egk1a_
void epb1_();
#define tinker_f_epb1 epb1_
void epb1a_();
#define tinker_f_epb1a epb1a_
void ediff1a_();
#define tinker_f_ediff1a ediff1a_
void ediff1b_();
#define tinker_f_ediff1b ediff1b_
void enp1_(double* ecav, double* edisp);
#define tinker_f_enp1 enp1_
void ewca1_(double* edisp);
#define tinker_f_ewca1 ewca1_
void ehpmf1_(double* ehp);
#define tinker_f_ehpmf1 ehpmf1_

// esolv2.f
void esolv2_(int* i);
#define tinker_f_esolv2 esolv2_
void esolv2a_(int* i);
#define tinker_f_esolv2a esolv2a_
void esolv2b_(int* nlist, int* list, int* reborn, int* reinduce);
#define tinker_f_esolv2b esolv2b_
void egb2a_(int* i);
#define tinker_f_egb2a egb2a_
void egb2b_(int* i);
#define tinker_f_egb2b egb2b_

// esolv3.f
void esolv3_();
#define tinker_f_esolv3 esolv3_
void egb3a_();
#define tinker_f_egb3a egb3a_
void egb3b_();
#define tinker_f_egb3b egb3b_
void egb3c_();
#define tinker_f_egb3c egb3c_
void egk3_();
#define tinker_f_egk3 egk3_
void egk3a_();
#define tinker_f_egk3a egk3a_
void epb3_();
#define tinker_f_epb3 epb3_
void ediff3_();
#define tinker_f_ediff3 ediff3_
void enp3_(double* ecav, double* aecav, double* edisp, double* aedisp);
#define tinker_f_enp3 enp3_
void ewca3_(double* edisp, double* aedisp);
#define tinker_f_ewca3 ewca3_
void ewca3x_(double* edisp, double* aedisp);
#define tinker_f_ewca3x ewca3x_
void ehpmf3_(double* ehp, int* nehp, double* aehp);
#define tinker_f_ehpmf3 ehpmf3_

// estrbnd.f
void estrbnd_();
#define tinker_f_estrbnd estrbnd_

// estrbnd1.f
void estrbnd1_();
#define tinker_f_estrbnd1 estrbnd1_

// estrbnd2.f
void estrbnd2_(int* iatom);
#define tinker_f_estrbnd2 estrbnd2_

// estrbnd3.f
void estrbnd3_();
#define tinker_f_estrbnd3 estrbnd3_

// estrtor.f
void estrtor_();
#define tinker_f_estrtor estrtor_

// estrtor1.f
void estrtor1_();
#define tinker_f_estrtor1 estrtor1_

// estrtor2.f
void estrtor2_(int* i);
#define tinker_f_estrtor2 estrtor2_

// estrtor3.f
void estrtor3_();
#define tinker_f_estrtor3 estrtor3_

// etors.f
void etors_();
#define tinker_f_etors etors_
void etors0a_();
#define tinker_f_etors0a etors0a_
void etors0b_();
#define tinker_f_etors0b etors0b_

// etors1.f
void etors1_();
#define tinker_f_etors1 etors1_
void etors1a_();
#define tinker_f_etors1a etors1a_
void etors1b_();
#define tinker_f_etors1b etors1b_

// etors2.f
void etors2_(int* i);
#define tinker_f_etors2 etors2_
void etors2a_(int* i);
#define tinker_f_etors2a etors2a_
void etors2b_(int* i);
#define tinker_f_etors2b etors2b_

// etors3.f
void etors3_();
#define tinker_f_etors3 etors3_
void etors3a_();
#define tinker_f_etors3a etors3a_
void etors3b_();
#define tinker_f_etors3b etors3b_

// etortor.f
void etortor_();
#define tinker_f_etortor etortor_
void chkttor_(int* ib, int* ic, int* id, double* sign, double* value1, double* value2);
#define tinker_f_chkttor chkttor_

// etortor1.f
void etortor1_();
#define tinker_f_etortor1 etortor1_

// etortor2.f
void etortor2_(int* i);
#define tinker_f_etortor2 etortor2_

// etortor3.f
void etortor3_();
#define tinker_f_etortor3 etortor3_

// eurey.f
void eurey_();
#define tinker_f_eurey eurey_

// eurey1.f
void eurey1_();
#define tinker_f_eurey1 eurey1_

// eurey2.f
void eurey2_(int* i);
#define tinker_f_eurey2 eurey2_

// eurey3.f
void eurey3_();
#define tinker_f_eurey3 eurey3_

// evcorr.f
void evcorr_(char* mode, double* elrc, tinker_fchar_len_t mode_cap);
inline void tinker_f_evcorr(tinker_fchars mode, double* elrc) {
    return evcorr_(mode.string, elrc, mode.capacity);
}
void evcorr1_(char* mode, double* elrc, double* vlrc, tinker_fchar_len_t mode_cap);
inline void tinker_f_evcorr1(tinker_fchars mode, double* elrc, double* vlrc) {
    return evcorr1_(mode.string, elrc, vlrc, mode.capacity);
}

// extra.f
void extra_();
#define tinker_f_extra extra_

// extra1.f
void extra1_();
#define tinker_f_extra1 extra1_

// extra2.f
void extra2_(int* i);
#define tinker_f_extra2 extra2_

// extra3.f
void extra3_();
#define tinker_f_extra3 extra3_

// fatal.f
void fatal_();
#define tinker_f_fatal fatal_

// fft3d.f
void fftsetup_();
#define tinker_f_fftsetup fftsetup_
void fftfront_();
#define tinker_f_fftfront fftfront_
void fftback_();
#define tinker_f_fftback fftback_
void fftclose_();
#define tinker_f_fftclose fftclose_

// fftpack.f
void cffti_(int* n, double* wsave, int* ifac);
#define tinker_f_cffti cffti_
void cffti1_(int* n, double* wa, int* ifac);
#define tinker_f_cffti1 cffti1_
void cfftf_(int* n, double* c, double* wsave, int* ifac);
#define tinker_f_cfftf cfftf_
void cfftf1_(int* n, double* c, double* ch, double* wa, int* ifac);
#define tinker_f_cfftf1 cfftf1_
void cfftb_(int* n, double* c, double* wsave, int* ifac);
#define tinker_f_cfftb cfftb_
void cfftb1_(int* n, double* c, double* ch, double* wa, int* ifac);
#define tinker_f_cfftb1 cfftb1_
void passf_(int* nac, int* ido, int* ip, int* l1, int* idl1, double* cc, double* c1, double* c2, double* ch, double* ch2, double* wa);
#define tinker_f_passf passf_
void passf2_(int* ido, int* l1, double* cc, double* ch, double* wa1);
#define tinker_f_passf2 passf2_
void passf3_(int* ido, int* l1, double* cc, double* ch, double* wa1, double* wa2);
#define tinker_f_passf3 passf3_
void passf4_(int* ido, int* l1, double* cc, double* ch, double* wa1, double* wa2, double* wa3);
#define tinker_f_passf4 passf4_
void passf5_(int* ido, int* l1, double* cc, double* ch, double* wa1, double* wa2, double* wa3, double* wa4);
#define tinker_f_passf5 passf5_
void passb_(int* nac, int* ido, int* ip, int* l1, int* idl1, double* cc, double* c1, double* c2, double* ch, double* ch2, double* wa);
#define tinker_f_passb passb_
void passb2_(int* ido, int* l1, double* cc, double* ch, double* wa1);
#define tinker_f_passb2 passb2_
void passb3_(int* ido, int* l1, double* cc, double* ch, double* wa1, double* wa2);
#define tinker_f_passb3 passb3_
void passb4_(int* ido, int* l1, double* cc, double* ch, double* wa1, double* wa2, double* wa3);
#define tinker_f_passb4 passb4_
void passb5_(int* ido, int* l1, double* cc, double* ch, double* wa1, double* wa2, double* wa3, double* wa4);
#define tinker_f_passb5 passb5_

// field.f
void field_();
#define tinker_f_field field_

// final.f
void final_();
#define tinker_f_final final_

// flatten.f
void flatten_();
#define tinker_f_flatten flatten_

// freeunit.f
int freeunit_();
#define tinker_f_freeunit freeunit_

// geometry.f
double geometry_(int* ia, int* ib, int* ic, int* id);
#define tinker_f_geometry geometry_

// getarc.f
void getarc_(int* iarc);
#define tinker_f_getarc getarc_

// getcart.f
void getcart_(int* ixyz);
#define tinker_f_getcart getcart_

// getdcd.f
void getdcd_(int* idcd);
#define tinker_f_getdcd getdcd_

// getint.f
void getint_();
#define tinker_f_getint getint_

// getkey.f
void getkey_();
#define tinker_f_getkey getkey_

// getmol.f
void getmol_();
#define tinker_f_getmol getmol_

// getmol2.f
void getmol2_();
#define tinker_f_getmol2 getmol2_

// getnumb.f
void getnumb_(char* string, int* number, int* next, tinker_fchar_len_t string_cap);
inline void tinker_f_getnumb(tinker_fchars string, int* number, int* next) {
    return getnumb_(string.string, number, next, string.capacity);
}

// getpdb.f
void getpdb_();
#define tinker_f_getpdb getpdb_

// getprm.f
void getprm_();
#define tinker_f_getprm getprm_

// getref.f
void getref_(int* iref);
#define tinker_f_getref getref_

// getstring.f
void getstring_(char* string, char* text, int* next, tinker_fchar_len_t string_cap, tinker_fchar_len_t text_cap);
inline void tinker_f_getstring(tinker_fchars string, tinker_fchars text, int* next) {
    return getstring_(string.string, text.string, next, string.capacity, text.capacity);
}

// gettext.f
void gettext_(char* string, char* text, int* next, tinker_fchar_len_t string_cap, tinker_fchar_len_t text_cap);
inline void tinker_f_gettext(tinker_fchars string, tinker_fchars text, int* next) {
    return gettext_(string.string, text.string, next, string.capacity, text.capacity);
}

// getword.f
void getword_(char* string, char* word, int* next, tinker_fchar_len_t string_cap, tinker_fchar_len_t word_cap);
inline void tinker_f_getword(tinker_fchars string, tinker_fchars word, int* next) {
    return getword_(string.string, word.string, next, string.capacity, word.capacity);
}

// getxyz.f
void getxyz_();
#define tinker_f_getxyz getxyz_

// ghmcstep.f
void ghmcstep_(int* istep, double* dt);
#define tinker_f_ghmcstep ghmcstep_
void ghmcterm_(int* istep, double* dt, double* alpha, double* beta);
#define tinker_f_ghmcterm ghmcterm_

// gradient.f
void gradient_(double* energy, double* derivs);
#define tinker_f_gradient gradient_

// gradrgd.f
void gradrgd_(double* energy, double* derivs);
#define tinker_f_gradrgd gradrgd_

// gradrot.f
void gradrot_(double* energy, double* derivs);
#define tinker_f_gradrot gradrot_

// groups.f
void groups_(int* proceed, double* weigh, int* ia, int* ib, int* ic, int* id, int* ie, int* ig);
#define tinker_f_groups groups_

// grpline.f
void grpline_();
#define tinker_f_grpline grpline_

// gyrate.f
void gyrate_(double* rg);
#define tinker_f_gyrate gyrate_

// hessian.f
void hessian_(double* h, int* hinit, int* hstop, int* hindex, double* hdiag);
#define tinker_f_hessian hessian_

// hessrgd.f
void hessrgd_(double* hrigid);
#define tinker_f_hessrgd hessrgd_

// hessrot.f
void hessrot_(char* mode, double* hrot, tinker_fchar_len_t mode_cap);
inline void tinker_f_hessrot(tinker_fchars mode, double* hrot) {
    return hessrot_(mode.string, hrot, mode.capacity);
}

// hybrid.f
void hybrid_();
#define tinker_f_hybrid hybrid_
void hatom_();
#define tinker_f_hatom hatom_
void hbond_();
#define tinker_f_hbond hbond_
void hangle_();
#define tinker_f_hangle hangle_
void hstrbnd_();
#define tinker_f_hstrbnd hstrbnd_
void himptor_();
#define tinker_f_himptor himptor_
void htors_();
#define tinker_f_htors htors_
void hstrtor_();
#define tinker_f_hstrtor hstrtor_
void hvdw_();
#define tinker_f_hvdw hvdw_
void hcharge_();
#define tinker_f_hcharge hcharge_
void hdipole_();
#define tinker_f_hdipole hdipole_

// image.f
void image_(double* xr, double* yr, double* zr);
#define tinker_f_image image_
void imager_(double* xr, double* yr, double* zr, int* i);
#define tinker_f_imager imager_
void imagen_(double* xr, double* yr, double* zr);
#define tinker_f_imagen imagen_

// impose.f
void impose_(int* n1, double* x1, double* y1, double* z1, int* n2, double* x2, double* y2, double* z2, double* rmsvalue);
#define tinker_f_impose impose_

// induce.f
void induce_();
#define tinker_f_induce induce_
void induce0a_();
#define tinker_f_induce0a induce0a_
void dfield0a_(double* field, double* fieldp);
#define tinker_f_dfield0a dfield0a_
void ufield0a_(double* field, double* fieldp);
#define tinker_f_ufield0a ufield0a_
void dfield0b_(double* field, double* fieldp);
#define tinker_f_dfield0b dfield0b_
void ufield0b_(double* field, double* fieldp);
#define tinker_f_ufield0b ufield0b_
void dfield0c_(double* field, double* fieldp);
#define tinker_f_dfield0c dfield0c_
void udirect1_(double* field);
#define tinker_f_udirect1 udirect1_
void udirect2a_(double* field, double* fieldp);
#define tinker_f_udirect2a udirect2a_
void udirect2b_(double* field, double* fieldp);
#define tinker_f_udirect2b udirect2b_
void ufield0c_(double* field, double* fieldp);
#define tinker_f_ufield0c ufield0c_
void umutual1_(double* field, double* fieldp);
#define tinker_f_umutual1 umutual1_
void umutual2a_(double* field, double* fieldp);
#define tinker_f_umutual2a umutual2a_
void umutual2b_(double* field, double* fieldp);
#define tinker_f_umutual2b umutual2b_
void induce0c_();
#define tinker_f_induce0c induce0c_
void dfield0d_(double* field, double* fieldp, double* fields, double* fieldps);
#define tinker_f_dfield0d dfield0d_
void ufield0d_(double* field, double* fieldp, double* fields, double* fieldps);
#define tinker_f_ufield0d ufield0d_
void induce0d_();
#define tinker_f_induce0d induce0d_
void dfield0e_(double* field, double* fieldp, double* fields, double* fieldps);
#define tinker_f_dfield0e dfield0e_
void ufield0e_(double* field, double* fieldp, double* fields, double* fieldps);
#define tinker_f_ufield0e ufield0e_
void ulspred_();
#define tinker_f_ulspred ulspred_
void uscale0a_(char* mode, double* rsd, double* rsdp, double* zrsd, double* zrsdp, tinker_fchar_len_t mode_cap);
inline void tinker_f_uscale0a(tinker_fchars mode, double* rsd, double* rsdp, double* zrsd, double* zrsdp) {
    return uscale0a_(mode.string, rsd, rsdp, zrsd, zrsdp, mode.capacity);
}
void uscale0b_(char* mode, double* rsd, double* rsdp, double* zrsd, double* zrsdp, tinker_fchar_len_t mode_cap);
inline void tinker_f_uscale0b(tinker_fchars mode, double* rsd, double* rsdp, double* zrsd, double* zrsdp) {
    return uscale0b_(mode.string, rsd, rsdp, zrsd, zrsdp, mode.capacity);
}

// inertia.f
void inertia_(int* mode);
#define tinker_f_inertia inertia_

// initatom.f
void initatom_();
#define tinker_f_initatom initatom_

// initial.f
void initial_();
#define tinker_f_initial initial_

// initprm.f
void initprm_();
#define tinker_f_initprm initprm_
void initmmff_();
#define tinker_f_initmmff initmmff_

// initres.f
void initres_();
#define tinker_f_initres initres_

// initrot.f
void initrot_();
#define tinker_f_initrot initrot_
int rotcheck_(int* base, int* partner);
#define tinker_f_rotcheck rotcheck_

// insert.f
void insert_(int* iatom);
#define tinker_f_insert insert_

// invbeta.f
double invbeta_(double* a, double* b, double* y);
#define tinker_f_invbeta invbeta_
double betai_(double* a, double* b, double* x);
#define tinker_f_betai betai_
double betacf_(double* a, double* b, double* x);
#define tinker_f_betacf betacf_
double gammln_(double* x);
#define tinker_f_gammln gammln_

// invert.f
void invert_(int* n, double* a);
#define tinker_f_invert invert_

// jacobi.f
void jacobi_(int* n, double* a, double* d, double* v);
#define tinker_f_jacobi jacobi_

// kangang.f
void kangang_();
#define tinker_f_kangang kangang_

// kangle.f
void kangle_();
#define tinker_f_kangle kangle_
void kanglem_();
#define tinker_f_kanglem kanglem_

// kangtor.f
void kangtor_();
#define tinker_f_kangtor kangtor_

// katom.f
void katom_();
#define tinker_f_katom katom_

// kbond.f
void kbond_();
#define tinker_f_kbond kbond_
void keneg_();
#define tinker_f_keneg keneg_
void kbondm_();
#define tinker_f_kbondm kbondm_

// kcharge.f
void kcharge_();
#define tinker_f_kcharge kcharge_
void kchargem_();
#define tinker_f_kchargem kchargem_

// kchgflx.f
void kchgflx_();
#define tinker_f_kchgflx kchgflx_

// kchgtrn.f
void kchgtrn_();
#define tinker_f_kchgtrn kchgtrn_

// kdipole.f
void kdipole_();
#define tinker_f_kdipole kdipole_

// kdisp.f
void kdisp_();
#define tinker_f_kdisp kdisp_

// kewald.f
void kewald_();
#define tinker_f_kewald kewald_
void ewaldcof_(double* alpha, double* cutoff);
#define tinker_f_ewaldcof ewaldcof_
void extent_(double* rmax);
#define tinker_f_extent extent_

// kexpol.f
void kexpol_();
#define tinker_f_kexpol kexpol_

// kextra.f
void kextra_();
#define tinker_f_kextra kextra_

// kgeom.f
void kgeom_();
#define tinker_f_kgeom kgeom_

// kimprop.f
void kimprop_();
#define tinker_f_kimprop kimprop_

// kimptor.f
void kimptor_();
#define tinker_f_kimptor kimptor_

// kinetic.f
void kinetic_(double* eksum, double* ekin, double* temp);
#define tinker_f_kinetic kinetic_
void kinaux_(double* temp_aux, double* temp_auxp);
#define tinker_f_kinaux kinaux_

// kmetal.f
void kmetal_();
#define tinker_f_kmetal kmetal_

// kmpole.f
void kmpole_();
#define tinker_f_kmpole kmpole_

// kopbend.f
void kopbend_();
#define tinker_f_kopbend kopbend_
void kopbendm_();
#define tinker_f_kopbendm kopbendm_

// kopdist.f
void kopdist_();
#define tinker_f_kopdist kopdist_

// korbit.f
void korbit_();
#define tinker_f_korbit korbit_

// kpitors.f
void kpitors_();
#define tinker_f_kpitors kpitors_

// kpolar.f
void kpolar_();
#define tinker_f_kpolar kpolar_
void polargrp_();
#define tinker_f_polargrp polargrp_

// krepel.f
void krepel_();
#define tinker_f_krepel krepel_

// ksolv.f
void ksolv_();
#define tinker_f_ksolv ksolv_
void ksa_();
#define tinker_f_ksa ksa_
void kgb_();
#define tinker_f_kgb kgb_
void kgk_();
#define tinker_f_kgk kgk_
void kpb_();
#define tinker_f_kpb kpb_
void knp_();
#define tinker_f_knp knp_
void khpmf_();
#define tinker_f_khpmf khpmf_
void setrad_(char* radtyp, tinker_fchar_len_t radtyp_cap);
inline void tinker_f_setrad(tinker_fchars radtyp) {
    return setrad_(radtyp.string, radtyp.capacity);
}

// kstrbnd.f
void kstrbnd_();
#define tinker_f_kstrbnd kstrbnd_
void kstrbndm_();
#define tinker_f_kstrbndm kstrbndm_

// kstrtor.f
void kstrtor_();
#define tinker_f_kstrtor kstrtor_

// ktors.f
void ktors_();
#define tinker_f_ktors ktors_
void ktorsm_();
#define tinker_f_ktorsm ktorsm_

// ktortor.f
void ktortor_();
#define tinker_f_ktortor ktortor_

// kurey.f
void kurey_();
#define tinker_f_kurey kurey_

// kvdw.f
void kvdw_();
#define tinker_f_kvdw kvdw_

// lattice.f
void lattice_();
#define tinker_f_lattice lattice_

// lbfgs.f
void lbfgs_(int* nvar, double* x0, double* minimum, double* grdmin, double (*fgvalue)(double*, double*), void (*optsave)(int*, double*, double*));
#define tinker_f_lbfgs lbfgs_

// lights.f
void lights_(double* cutoff, int* nsite, double* xsort, double* ysort, double* zsort, int* unique);
#define tinker_f_lights lights_

// lusolve.f
void lusolve_(int* nvar, double* a, double* b);
#define tinker_f_lusolve lusolve_

// makeint.f
void makeint_(int* mode);
#define tinker_f_makeint makeint_
int adjacent_(int* i1, int* i2, int* mode, int* more, int* iz0, int* iz1);
#define tinker_f_adjacent adjacent_

// makeref.f
void makeref_(int* iref);
#define tinker_f_makeref makeref_

// makexyz.f
void makexyz_();
#define tinker_f_makexyz makexyz_

// maxwell.f
double maxwell_(double* mass, double* temper);
#define tinker_f_maxwell maxwell_

// mdinit.f
void mdinit_(double* dt);
#define tinker_f_mdinit mdinit_

// mdrest.f
void mdrest_(int* istep);
#define tinker_f_mdrest mdrest_

// mdsave.f
void mdsave_(int* istep, double* dt, double* epot, double* eksum);
#define tinker_f_mdsave mdsave_

// mdstat.f
void mdstat_(int* istep, double* dt, double* etot, double* epot, double* ekin, double* temp, double* pres);
#define tinker_f_mdstat mdstat_

// mechanic.f
void mechanic_();
#define tinker_f_mechanic mechanic_

// merge.f
void merge_(int* iref);
#define tinker_f_merge merge_

// molecule.f
void molecule_();
#define tinker_f_molecule molecule_

// moments.f
void moments_(char* mode, tinker_fchar_len_t mode_cap);
inline void tinker_f_moments(tinker_fchars mode) {
    return moments_(mode.string, mode.capacity);
}

// mutate.f
void mutate_();
#define tinker_f_mutate mutate_
void altelec_();
#define tinker_f_altelec altelec_
void alttors_(int* ntbnd, int* itbnd);
#define tinker_f_alttors alttors_

// nblist.f
void nblist_();
#define tinker_f_nblist nblist_
void vlist_();
#define tinker_f_vlist vlist_
void vbuild_();
#define tinker_f_vbuild vbuild_
void vlight_();
#define tinker_f_vlight vlight_
void dlist_();
#define tinker_f_dlist dlist_
void dbuild_();
#define tinker_f_dbuild dbuild_
void dlight_();
#define tinker_f_dlight dlight_
void clist_();
#define tinker_f_clist clist_
void cbuild_();
#define tinker_f_cbuild cbuild_
void clight_();
#define tinker_f_clight clight_
void mlist_();
#define tinker_f_mlist mlist_
void mbuild_();
#define tinker_f_mbuild mbuild_
void mlight_();
#define tinker_f_mlight mlight_
void ulist_();
#define tinker_f_ulist ulist_
void ubuild_();
#define tinker_f_ubuild ubuild_
void ulight_();
#define tinker_f_ulight ulight_

// nextarg.f
void nextarg_(char* string, int* exist, tinker_fchar_len_t string_cap);
inline void tinker_f_nextarg(tinker_fchars string, int* exist) {
    return nextarg_(string.string, exist, string.capacity);
}

// nexttext.f
int nexttext_(char* string, tinker_fchar_len_t string_cap);
inline int tinker_f_nexttext(tinker_fchars string) {
    return nexttext_(string.string, string.capacity);
}

// nose.f
void nose_(int* istep, double* dt);
#define tinker_f_nose nose_
void hoover_(double* dt, double* press);
#define tinker_f_hoover hoover_

// nspline.f
void nspline_(int* n, double* x0, double* y0, double* s1, double* s2, double* h, double* g, double* dy, double* dla, double* dmu);
#define tinker_f_nspline nspline_

// number.f
int number_(char* string, tinker_fchar_len_t string_cap);
inline int tinker_f_number(tinker_fchars string) {
    return number_(string.string, string.capacity);
}

// numeral.f
void numeral_(int* number, char* string, int* size, tinker_fchar_len_t string_cap);
inline void tinker_f_numeral(int* number, tinker_fchars string, int* size) {
    return numeral_(number, string.string, size, string.capacity);
}

// numgrad.f
void numgrad_(double (*evalue)(), double* g, double* eps);
#define tinker_f_numgrad numgrad_

// ocvm.f
void ocvm_(int* nvar, double* x0, double* f0, double* grdmin, double (*fgvalue)(double*, double*), void (*optsave)(int*, double*, double*));
#define tinker_f_ocvm ocvm_

// openend.f
void openend_(int* iunit, char* name, tinker_fchar_len_t name_cap);
inline void tinker_f_openend(int* iunit, tinker_fchars name) {
    return openend_(iunit, name.string, name.capacity);
}

// optinit.f
void optinit_();
#define tinker_f_optinit optinit_

// optsave.f
void optsave_(int* ncycle, double* f, double* xx);
#define tinker_f_optsave optsave_

// orbital.f
void orbital_();
#define tinker_f_orbital orbital_
void piplane_();
#define tinker_f_piplane piplane_

// orient.f
void orient_();
#define tinker_f_orient orient_
void xyzrigid_();
#define tinker_f_xyzrigid xyzrigid_
void roteuler_(double* a, double* phi, double* theta, double* psi);
#define tinker_f_roteuler roteuler_
void rigidxyz_();
#define tinker_f_rigidxyz rigidxyz_

// orthog.f
void orthog_(int* m, int* n, double* a);
#define tinker_f_orthog orthog_

// overlap.f
void overlap_(int* atmnum1, int* atmnum2, double* rang, double* ovlap);
#define tinker_f_overlap overlap_
void slater_(int* na, int* la, double* za, int* nb, int* lb, double* zb, double* r, double* s);
#define tinker_f_slater slater_
void polyp_(double* c, int* ia, int* ib, int* max, double* d, int* iaa, int* ibb, int* n);
#define tinker_f_polyp polyp_
double cjkm_(int* j, int* k, int* m);
#define tinker_f_cjkm cjkm_
void aset_(double* alpha, int* n, double* a);
#define tinker_f_aset aset_
void bset_(double* beta, int* n, double* b);
#define tinker_f_bset bset_
double bmax_(double* beta, int* n);
#define tinker_f_bmax bmax_

// picalc.f
void picalc_();
#define tinker_f_picalc picalc_
void piscf_();
#define tinker_f_piscf piscf_
void pitilt_(double* povlap);
#define tinker_f_pitilt pitilt_
void pimove_(int* list, double* xr, double* yr, double* zr);
#define tinker_f_pimove pimove_
void pialter_();
#define tinker_f_pialter pialter_

// pmestuf.f
void getchunk_();
#define tinker_f_getchunk getchunk_
void moduli_();
#define tinker_f_moduli moduli_
void bspline_(double* x, int* n, double* c);
#define tinker_f_bspline bspline_
void dftmod_(double* bsmod, double* bsarray, int* nfft, int* order);
#define tinker_f_dftmod dftmod_
void bspline_fill_();
#define tinker_f_bspline_fill bspline_fill_
void bsplgen_(double* w, double* thetai);
#define tinker_f_bsplgen bsplgen_
void table_fill_();
#define tinker_f_table_fill table_fill_
void setchunk_(int* i, int* cid, int* off1, int* off2, int* off3);
#define tinker_f_setchunk setchunk_
void grid_pchg_();
#define tinker_f_grid_pchg grid_pchg_
void grid_mpole_(double* fmp);
#define tinker_f_grid_mpole grid_mpole_
void grid_uind_(double* fuind, double* fuinp);
#define tinker_f_grid_uind grid_uind_
void grid_disp_();
#define tinker_f_grid_disp grid_disp_
void adjust_(int* offset, int* nfft, int* nchk, int* amin, int* amax, int* cmin, int* cmax);
#define tinker_f_adjust adjust_
void fphi_pchg_(double* fphi);
#define tinker_f_fphi_pchg fphi_pchg_
void fphi_mpole_(double* fphi);
#define tinker_f_fphi_mpole fphi_mpole_
void fphi_uind_(double* fdip_phi1, double* fdip_phi2, double* fdip_sum_phi);
#define tinker_f_fphi_uind fphi_uind_
void cmp_to_fmp_(double* cmp, double* fmp);
#define tinker_f_cmp_to_fmp cmp_to_fmp_
void cart_to_frac_(double* ctf);
#define tinker_f_cart_to_frac cart_to_frac_
void fphi_to_cphi_(double* fphi, double* cphi);
#define tinker_f_fphi_to_cphi fphi_to_cphi_
void frac_to_cart_(double* ftc);
#define tinker_f_frac_to_cart frac_to_cart_

// pmpb.f
void apbsinitial_(int* dime, double* grid, double* gcent, double* cgrid, double* cgcent, double* fgrid, double* fgcent, double* pdie, double* sdie, double* srad, double* swin, double* sdens, double* kelvin, int* ionn, double* ionc, int* ionq, double* ionr, char* pbtyp, int* pbtyplen, char* pbsoln, int* pbsolnlen, char* bcfl, int* bcfllen, char* chgm, int* chgmlen, char* srfm, int* srfmlen, tinker_fchar_len_t pbtyp_cap, tinker_fchar_len_t pbsoln_cap, tinker_fchar_len_t bcfl_cap, tinker_fchar_len_t chgm_cap, tinker_fchar_len_t srfm_cap);
inline void tinker_f_apbsinitial(int* dime, double* grid, double* gcent, double* cgrid, double* cgcent, double* fgrid, double* fgcent, double* pdie, double* sdie, double* srad, double* swin, double* sdens, double* kelvin, int* ionn, double* ionc, int* ionq, double* ionr, tinker_fchars pbtyp, int* pbtyplen, tinker_fchars pbsoln, int* pbsolnlen, tinker_fchars bcfl, int* bcfllen, tinker_fchars chgm, int* chgmlen, tinker_fchars srfm, int* srfmlen) {
    return apbsinitial_(dime, grid, gcent, cgrid, cgcent, fgrid, fgcent, pdie, sdie, srad, swin, sdens, kelvin, ionn, ionc, ionq, ionr, pbtyp.string, pbtyplen, pbsoln.string, pbsolnlen, bcfl.string, bcfllen, chgm.string, chgmlen, srfm.string, srfmlen, pbtyp.capacity, pbsoln.capacity, bcfl.capacity, chgm.capacity, srfm.capacity);
}
void apbsempole_(int* n, double* pos, double* rsolv, double* pbpole, double* pbe, double* apbe, double* pbep, double* pbfp, double* pbtp);
#define tinker_f_apbsempole apbsempole_
void apbsinduce_(double* indpole, double* pbeuind);
#define tinker_f_apbsinduce apbsinduce_
void apbsnlinduce_(double* inppole, double* pbeuinp);
#define tinker_f_apbsnlinduce apbsnlinduce_
void pbdirectpolforce_(double* indpole, double* inppole, double* directf, double* directt);
#define tinker_f_pbdirectpolforce pbdirectpolforce_
void pbmutualpolforce_(double* indpole, double* inppole, double* mutualf);
#define tinker_f_pbmutualpolforce pbmutualpolforce_
void apbsfinal_();
#define tinker_f_apbsfinal apbsfinal_

// polymer.f
void polymer_();
#define tinker_f_polymer polymer_

// predict.f
void predict_();
#define tinker_f_predict predict_
void auxinit_();
#define tinker_f_auxinit auxinit_

// pressure.f
void pressure_(double* dt, double* epot, double* ekin, double* temp, double* pres, double* stress);
#define tinker_f_pressure pressure_
void pressure2_(double* epot, double* temp);
#define tinker_f_pressure2 pressure2_
void pmonte_(double* epot, double* temp);
#define tinker_f_pmonte pmonte_
void pscale_(double* dt, double* pres, double* stress);
#define tinker_f_pscale pscale_

// prmkey.f
void prmkey_(char* text, tinker_fchar_len_t text_cap);
inline void tinker_f_prmkey(tinker_fchars text) {
    return prmkey_(text.string, text.capacity);
}
void potoff_();
#define tinker_f_potoff potoff_
void valoff_();
#define tinker_f_valoff valoff_
void nbondoff_();
#define tinker_f_nbondoff nbondoff_

// promo.f
void promo_();
#define tinker_f_promo promo_

// prtarc.f
void prtarc_(int* iarc, int* first);
#define tinker_f_prtarc prtarc_
void prtarcf_(int* iarc);
#define tinker_f_prtarcf prtarcf_
void prtarcb_(int* idcd, int* first);
#define tinker_f_prtarcb prtarcb_

// prtdcd.f
void prtdcd_(int* idcd, int* first);
#define tinker_f_prtdcd prtdcd_

// prtdyn.f
void prtdyn_();
#define tinker_f_prtdyn prtdyn_

// prterr.f
void prterr_();
#define tinker_f_prterr prterr_

// prtint.f
void prtint_(int* izmt);
#define tinker_f_prtint prtint_

// prtmol2.f
void prtmol2_(int* imol2);
#define tinker_f_prtmol2 prtmol2_
void setmol2_(char* atmnam, char* atmtyp, double* atmchg, char* bndtyp, tinker_fchar_len_t atmnam_cap, tinker_fchar_len_t atmtyp_cap, tinker_fchar_len_t bndtyp_cap);
inline void tinker_f_setmol2(tinker_fchars atmnam, tinker_fchars atmtyp, double* atmchg, tinker_fchars bndtyp) {
    return setmol2_(atmnam.string, atmtyp.string, atmchg, bndtyp.string, atmnam.capacity, atmtyp.capacity, bndtyp.capacity);
}

// prtpdb.f
void prtpdb_(int* ipdb);
#define tinker_f_prtpdb prtpdb_

// prtprm.f
void prtprm_(int* itxt);
#define tinker_f_prtprm prtprm_

// prtseq.f
void prtseq_(int* iseq);
#define tinker_f_prtseq prtseq_

// prtxyz.f
void prtxyz_(int* ixyz);
#define tinker_f_prtxyz prtxyz_

// qrsolve.f
void qrfact_(int* n, int* m, double* a, int* pivot, int* ipvt, double* rdiag);
#define tinker_f_qrfact qrfact_
void qrsolve_(int* n, int* np, double* r, int* ipvt, double* diag, double* qtb, double* x, double* sdiag, double* xpvt);
#define tinker_f_qrsolve qrsolve_

// quatfit.f
void quatfit_(int* n1, double* x1, double* y1, double* z1, int* n2, double* x2, double* y2, double* z2);
#define tinker_f_quatfit quatfit_

// random.f
double random_();
#define tinker_f_random random_
double normal_();
#define tinker_f_normal normal_
void ranvec_(double* vector);
#define tinker_f_ranvec ranvec_
void sphere_(int* ndot, double* dot);
#define tinker_f_sphere sphere_

// rattle.f
void rattle_(double* dt, double* xold, double* yold, double* zold);
#define tinker_f_rattle rattle_
void rattle2_(double* dt);
#define tinker_f_rattle2 rattle2_
void shake_(double* xold, double* yold, double* zold);
#define tinker_f_shake shake_
void shake2_(double* derivs);
#define tinker_f_shake2 shake2_

// readcart.f
void readcart_(int* ixyz, int* first);
#define tinker_f_readcart readcart_

// readdcd.f
void readdcd_(int* idcd, int* first);
#define tinker_f_readdcd readdcd_

// readdyn.f
void readdyn_(int* idyn);
#define tinker_f_readdyn readdyn_

// readgau.f
void readgau_();
#define tinker_f_readgau readgau_
void readgarc_(int* igau, char* string, char* word, int* length, int* next, tinker_fchar_len_t string_cap, tinker_fchar_len_t word_cap);
inline void tinker_f_readgarc(int* igau, tinker_fchars string, tinker_fchars word, int* length, int* next) {
    return readgarc_(igau, string.string, word.string, length, next, string.capacity, word.capacity);
}

// readgdma.f
void readgdma_(int* idma);
#define tinker_f_readgdma readgdma_
void match1_(int* i, char* record, tinker_fchar_len_t record_cap);
inline void tinker_f_match1(int* i, tinker_fchars record) {
    return match1_(i, record.string, record.capacity);
}
void match2_(int* i, char* record, tinker_fchar_len_t record_cap);
inline void tinker_f_match2(int* i, tinker_fchars record) {
    return match2_(i, record.string, record.capacity);
}
void match3_(int* i, char* record, tinker_fchar_len_t record_cap);
inline void tinker_f_match3(int* i, tinker_fchars record) {
    return match3_(i, record.string, record.capacity);
}

// readint.f
void readint_(int* izmt);
#define tinker_f_readint readint_

// readmol.f
void readmol_(int* imdl);
#define tinker_f_readmol readmol_

// readmol2.f
void readmol2_(int* imol2);
#define tinker_f_readmol2 readmol2_

// readpdb.f
void readpdb_(int* ipdb);
#define tinker_f_readpdb readpdb_
void scanpdb_(int* ipdb);
#define tinker_f_scanpdb scanpdb_
void fixpdb_(char* resname, char* atmname, tinker_fchar_len_t resname_cap, tinker_fchar_len_t atmname_cap);
inline void tinker_f_fixpdb(tinker_fchars resname, tinker_fchars atmname) {
    return fixpdb_(resname.string, atmname.string, resname.capacity, atmname.capacity);
}

// readprm.f
void readprm_();
#define tinker_f_readprm readprm_

// readseq.f
void readseq_(int* iseq);
#define tinker_f_readseq readseq_

// readxyz.f
void readxyz_(int* ixyz);
#define tinker_f_readxyz readxyz_

// replica.f
void replica_(double* cutoff);
#define tinker_f_replica replica_

// respa.f
void respa_(int* istep, double* dt);
#define tinker_f_respa respa_
void gradfast_(double* energy, double* derivs);
#define tinker_f_gradfast gradfast_
void gradslow_(double* energy, double* derivs);
#define tinker_f_gradslow gradslow_

// rgdstep.f
void rgdstep_(int* istep, double* dt);
#define tinker_f_rgdstep rgdstep_
void rotrgd_(double* dfi, double* arot);
#define tinker_f_rotrgd rotrgd_
void linbody_(int* i, double* inert, double* wcp);
#define tinker_f_linbody linbody_

// rings.f
void rings_();
#define tinker_f_rings rings_

// rmsfit.f
double rmsfit_(double* x1, double* y1, double* z1, double* x2, double* y2, double* z2);
#define tinker_f_rmsfit rmsfit_

// rotlist.f
void rotlist_(int* base, int* partner);
#define tinker_f_rotlist rotlist_

// rotpole.f
void rotpole_(double* inpole, double* outpole);
#define tinker_f_rotpole rotpole_
void rotrpole_(double* inpole, double* outpole);
#define tinker_f_rotrpole rotrpole_
void rotmat_(int* i, double* a, int* planar);
#define tinker_f_rotmat rotmat_
void rotsite_(int* isite, double* a, int* planar, double* inpole, double* outpole);
#define tinker_f_rotsite rotsite_

// sdstep.f
void sdstep_(int* istep, double* dt);
#define tinker_f_sdstep sdstep_
void sdterm_(int* istep, double* dt, double* pfric, double* vfric, double* afric, double* prand, double* vrand);
#define tinker_f_sdterm sdterm_
void sdarea_(int* istep);
#define tinker_f_sdarea sdarea_

// search.f
void search_(int* nvar, double* f, double* g, double* x, double* p, double* f_move, double* angle, int* ncalls, double (*fgvalue)(double*, double*), char* status, tinker_fchar_len_t status_cap);
inline void tinker_f_search(int* nvar, double* f, double* g, double* x, double* p, double* f_move, double* angle, int* ncalls, double (*fgvalue)(double*, double*), tinker_fchars status) {
    return search_(nvar, f, g, x, p, f_move, angle, ncalls, fgvalue, status.string, status.capacity);
}

// server.f
void chksocket_(int* flag);
#define tinker_f_chksocket chksocket_
void createjvm_(int* flag);
#define tinker_f_createjvm createjvm_
void destroyjvm_();
#define tinker_f_destroyjvm destroyjvm_
void createserver_(int* flag);
#define tinker_f_createserver createserver_
void destroyserver_();
#define tinker_f_destroyserver destroyserver_
void createsystem_(int* n, int* nkey, int* flag);
#define tinker_f_createsystem createsystem_
void getmonitor_();
#define tinker_f_getmonitor getmonitor_
void releasemonitor_();
#define tinker_f_releasemonitor releasemonitor_
void createupdate_(int* n, int* mode, int* amoeba, int* flag);
#define tinker_f_createupdate createupdate_
void needupdate_(int* flag);
#define tinker_f_needupdate needupdate_
void setupdated_();
#define tinker_f_setupdated setupdated_
void setfile_(char* filename, tinker_fchar_len_t filename_cap);
inline void tinker_f_setfile(tinker_fchars filename) {
    return setfile_(filename.string, filename.capacity);
}
void setforcefield_(char* forcefield, tinker_fchar_len_t forcefield_cap);
inline void tinker_f_setforcefield(tinker_fchars forcefield) {
    return setforcefield_(forcefield.string, forcefield.capacity);
}
void setkeyword_(int* i, char* keyline, tinker_fchar_len_t keyline_cap);
inline void tinker_f_setkeyword(int* i, tinker_fchars keyline) {
    return setkeyword_(i, keyline.string, keyline.capacity);
}
void setatomtypes_(int* n, int* type);
#define tinker_f_setatomtypes setatomtypes_
void setatomic_(int* n, int* atomic);
#define tinker_f_setatomic setatomic_
void setmass_(int* n, double* mass);
#define tinker_f_setmass setmass_
void setcharge_(int* n, double* charge);
#define tinker_f_setcharge setcharge_
void setconnectivity_(int* n, int* b1, int* b2, int* b3, int* b4);
#define tinker_f_setconnectivity setconnectivity_
void setname_(int* i, char* name, tinker_fchar_len_t name_cap);
inline void tinker_f_setname(int* i, tinker_fchars name) {
    return setname_(i, name.string, name.capacity);
}
void setstory_(int* i, char* story, tinker_fchar_len_t story_cap);
inline void tinker_f_setstory(int* i, tinker_fchars story) {
    return setstory_(i, story.string, story.capacity);
}
void setcoordinates_(int* n, double* x, double* y, double* z);
#define tinker_f_setcoordinates setcoordinates_
void setstep_(int* ncycle);
#define tinker_f_setstep setstep_
void setmdtime_(double* time);
#define tinker_f_setmdtime setmdtime_
void setenergy_(double* energy);
#define tinker_f_setenergy setenergy_
void setgradients_(int* n, double* dx, double* dy, double* dz);
#define tinker_f_setgradients setgradients_
void setvelocity_(int* n, double* vx, double* vy, double* vz);
#define tinker_f_setvelocity setvelocity_
void setacceleration_(int* n, double* ax, double* ay, double* az);
#define tinker_f_setacceleration setacceleration_
void setinduced_(int* n, double* ux, double* uy, double* uz);
#define tinker_f_setinduced setinduced_

// setprm.f
void setprm_();
#define tinker_f_setprm setprm_

// shakeup.f
void shakeup_();
#define tinker_f_shakeup shakeup_
void chkangle_(int* ia, int* ib, int* ic);
#define tinker_f_chkangle chkangle_

// sigmoid.f
double sigmoid_(double* beta, double* x);
#define tinker_f_sigmoid sigmoid_

// simplex.f
void simplex_(int* nvar, int* iter, int* ntest, double* x0, double* y0, double* step, double* toler, double (*fvalue)(double*));
#define tinker_f_simplex simplex_

// sktstuf.f
void sktdyn_(int* istep, double* dt, double* epot);
#define tinker_f_sktdyn sktdyn_
void sktopt_(int* ncycle, double* eopt);
#define tinker_f_sktopt sktopt_
void sktinit_();
#define tinker_f_sktinit sktinit_
void sktkill_();
#define tinker_f_sktkill sktkill_

// sort.f
void sort_(int* n, int* list);
#define tinker_f_sort sort_
void sort2_(int* n, double* list, int* key);
#define tinker_f_sort2 sort2_
void sort3_(int* n, int* list, int* key);
#define tinker_f_sort3 sort3_
void sort4_(int* n, int* list);
#define tinker_f_sort4 sort4_
void sort5_(int* n, int* list, int* m);
#define tinker_f_sort5 sort5_
void sort6_(int* n, char* list, tinker_fchar_len_t list_cap);
inline void tinker_f_sort6(int* n, tinker_fchars list) {
    return sort6_(n, list.string, list.capacity);
}
void sort7_(int* n, char* list, int* key, tinker_fchar_len_t list_cap);
inline void tinker_f_sort7(int* n, tinker_fchars list, int* key) {
    return sort7_(n, list.string, key, list.capacity);
}
void sort8_(int* n, int* list);
#define tinker_f_sort8 sort8_
void sort9_(int* n, double* list);
#define tinker_f_sort9 sort9_
void sort10_(int* n, char* list, tinker_fchar_len_t list_cap);
inline void tinker_f_sort10(int* n, tinker_fchars list) {
    return sort10_(n, list.string, list.capacity);
}

// square.f
void square_(int* n, int* m, double* xlo, double* xhi, double* xc, double* fc, double* gc, double* fjac, double* grdmin, void (*rsdvalue)(int*, int*, double*, double*), void (*lsqwrite)(int*, int*, double*, double*, double*));
#define tinker_f_square square_
void lmstep_(int* n, int* m, double* ga, double* a, int* ipvt, double* xscale, double* qtf, double* stpmax, double* delta, double* amu, int* first, double* sa, int* gauss);
#define tinker_f_lmstep lmstep_
void trust_(int* n, int* m, double* xc, double* fcnorm, double* gc, double* a, int* ipvt, double* sc, double* sa, double* xscale, int* gauss, double* stpmax, double* delta, int* icode, double* xp, double* xpprev, double* fc, double* fp, double* fpnorm, double* fpprev, int* bigstp, int* ncalls, double* xlo, double* xhi, int* nactive, double* stpmin, double* rftol, double* faketol, void (*rsdvalue)(int*, int*, double*, double*));
#define tinker_f_trust trust_

// suffix.f
void suffix_(char* string, char* extension, char* status, tinker_fchar_len_t string_cap, tinker_fchar_len_t extension_cap, tinker_fchar_len_t status_cap);
inline void tinker_f_suffix(tinker_fchars string, tinker_fchars extension, tinker_fchars status) {
    return suffix_(string.string, extension.string, status.string, string.capacity, extension.capacity, status.capacity);
}

// surface.f
void surface_(double* total, double* area, double* radius, double* weight, double* probe);
#define tinker_f_surface surface_
void surface1_(double* total, double* area, double* darea, double* radius, double* weight, double* probe);
#define tinker_f_surface1 surface1_

// surfatom.f
void surfatom_(int* ir, double* area, double* radius);
#define tinker_f_surfatom surfatom_
void surfatom1_(int* ir, double* area, double* darea, double* radius);
#define tinker_f_surfatom1 surfatom1_

// switch.f
void switch_(char* mode, tinker_fchar_len_t mode_cap);
inline void tinker_f_switch(tinker_fchars mode) {
    return switch_(mode.string, mode.capacity);
}

// tcgstuf.f
void induce0b_();
#define tinker_f_induce0b induce0b_
void indtcga_();
#define tinker_f_indtcga indtcga_
void indtcgb_();
#define tinker_f_indtcgb indtcgb_
void tcg_alphaquad_(double* sum, double* a, double* b);
#define tinker_f_tcg_alphaquad tcg_alphaquad_
void tcg_resource_(int* order);
#define tinker_f_tcg_resource tcg_resource_
void tcg_alpha12_(double* source1, double* source2);
#define tinker_f_tcg_alpha12 tcg_alpha12_
void tcg_alpha22_(double* source1, double* source2, double* result1, double* result2);
#define tinker_f_tcg_alpha22 tcg_alpha22_
void tcg_dotprod_(double* sum, int* n, double* a, double* b);
#define tinker_f_tcg_dotprod tcg_dotprod_
void tcg_ufield_(double* ind, double* inp, double* v3d, double* v3p);
#define tinker_f_tcg_ufield tcg_ufield_
void tcg_t0_(double* ind, double* inp, double* v3d, double* v3p);
#define tinker_f_tcg_t0 tcg_t0_
void tcgswap_(double* uind1, double* uinp1, double* uind2, double* uinp2);
#define tinker_f_tcgswap tcgswap_
void tcg_update_(double* pvec, double* rvec, double* beta);
#define tinker_f_tcg_update tcg_update_

// temper.f
void temper_(double* dt, double* eksum, double* ekin, double* temp);
#define tinker_f_temper temper_
void temper2_(double* dt, double* temp);
#define tinker_f_temper2 temper2_

// tncg.f
void tncg_(char* mode, char* method, int* nvar, double* x0, double* minimum, double* grdmin, double (*fgvalue)(double*, double*), void (*hmatrix)(char*, double*, double*, int*, int*, int*, double*, tinker_fchar_len_t), void (*optsave)(int*, double*, double*), tinker_fchar_len_t mode_cap, tinker_fchar_len_t method_cap);
inline void tinker_f_tncg(tinker_fchars mode, tinker_fchars method, int* nvar, double* x0, double* minimum, double* grdmin, double (*fgvalue)(double*, double*), void (*hmatrix)(char*, double*, double*, int*, int*, int*, double*, tinker_fchar_len_t), void (*optsave)(int*, double*, double*)) {
    return tncg_(mode.string, method.string, nvar, x0, minimum, grdmin, fgvalue, hmatrix, optsave, mode.capacity, method.capacity);
}
void tnsolve_(char* mode, char* method, int* negtest, int* nvar, double* p, double* x0, double* g, double* h, int* h_init, int* h_stop, int* h_index, double* h_diag, int* cycle, int* iter_cg, int* fg_call, double (*fgvalue)(double*, double*), char* status, tinker_fchar_len_t mode_cap, tinker_fchar_len_t method_cap, tinker_fchar_len_t status_cap);
inline void tinker_f_tnsolve(tinker_fchars mode, tinker_fchars method, int* negtest, int* nvar, double* p, double* x0, double* g, double* h, int* h_init, int* h_stop, int* h_index, double* h_diag, int* cycle, int* iter_cg, int* fg_call, double (*fgvalue)(double*, double*), tinker_fchars status) {
    return tnsolve_(mode.string, method.string, negtest, nvar, p, x0, g, h, h_init, h_stop, h_index, h_diag, cycle, iter_cg, fg_call, fgvalue, status.string, mode.capacity, method.capacity, status.capacity);
}
void precond_(char* method, int* iter, int* nvar, double* s, double* r, double* h, int* h_init, int* h_stop, int* h_index, double* h_diag, tinker_fchar_len_t method_cap);
inline void tinker_f_precond(tinker_fchars method, int* iter, int* nvar, double* s, double* r, double* h, int* h_init, int* h_stop, int* h_index, double* h_diag) {
    return precond_(method.string, iter, nvar, s, r, h, h_init, h_stop, h_index, h_diag, method.capacity);
}

// torphase.f
void torphase_(int* ft, double* vt, double* st);
#define tinker_f_torphase torphase_

// torque.f
void torque_(int* i, double* trq, double* frcx, double* frcy, double* frcz, double* de);
#define tinker_f_torque torque_

// torsions.f
void torsions_();
#define tinker_f_torsions torsions_

// trimtext.f
int trimtext_(char* string, tinker_fchar_len_t string_cap);
inline int tinker_f_trimtext(tinker_fchars string) {
    return trimtext_(string.string, string.capacity);
}
void trimhead_(char* string, tinker_fchar_len_t string_cap);
inline void tinker_f_trimhead(tinker_fchars string) {
    return trimhead_(string.string, string.capacity);
}
void justify_(char* string, tinker_fchar_len_t string_cap);
inline void tinker_f_justify(tinker_fchars string) {
    return justify_(string.string, string.capacity);
}
void upcase_(char* string, tinker_fchar_len_t string_cap);
inline void tinker_f_upcase(tinker_fchars string) {
    return upcase_(string.string, string.capacity);
}
void lowcase_(char* string, tinker_fchar_len_t string_cap);
inline void tinker_f_lowcase(tinker_fchars string) {
    return lowcase_(string.string, string.capacity);
}

// unitcell.f
void unitcell_();
#define tinker_f_unitcell unitcell_

// verlet.f
void verlet_(int* istep, double* dt);
#define tinker_f_verlet verlet_

// version.f
void version_(char* string, char* status, tinker_fchar_len_t string_cap, tinker_fchar_len_t status_cap);
inline void tinker_f_version(tinker_fchars string, tinker_fchars status) {
    return version_(string.string, status.string, string.capacity, status.capacity);
}

// volume.f
void volume_(double* volume_tot, double* radius, double* exclude);
#define tinker_f_volume volume_
void volume1_(double* radius, double* probe, double* dex);
#define tinker_f_volume1 volume1_
void volume2_(int* iatom, double* radius, double* probe, double* xhess, double* yhess, double* zhess);
#define tinker_f_volume2 volume2_

// xyzatm.f
void xyzatm_(int* i, int* ia, double* bond, int* ib, double* angle1, int* ic, double* angle2, int* chiral);
#define tinker_f_xyzatm xyzatm_

// zatom.f
void zatom_(int* bionum, double* bond, double* angle, double* dihed, int* iz1, int* iz2, int* iz3, int* iz4);
#define tinker_f_zatom zatom_

#ifdef __cplusplus
}
#endif
