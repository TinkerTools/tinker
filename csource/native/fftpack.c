#include "f2c.h"


int cfftb_(n, c, wsave)
integer *n;
doublereal *c, *wsave;
{
    extern int cfftb1_();
    static integer iw1, iw2;

    --wsave;
    --c;

    if (*n == 1) {
	return 0;
    }

    iw1 = *n + *n + 1;
    iw2 = iw1 + *n + *n;
    cfftb1_(n, &c[1], &wsave[1], &wsave[iw1], &wsave[iw2]);

    return 0;
}


int cfftf_(n, c, wsave)
integer *n;
doublereal *c, *wsave;
{
    extern int cfftf1_();
    static integer iw1, iw2;

    --wsave;
    --c;

    if (*n == 1) {
	return 0;
    }

    iw1 = *n + *n + 1;
    iw2 = iw1 + *n + *n;
    cfftf1_(n, &c[1], &wsave[1], &wsave[iw1], &wsave[iw2]);

    return 0;
}


int cffti_(n, wsave)
integer *n;
doublereal *wsave;
{
    extern int cffti1_();
    static integer iw1, iw2;

    --wsave;

    if (*n == 1) {
	return 0;
    }

    iw1 = *n + *n + 1;
    iw2 = iw1 + *n + *n;
    cffti1_(n, &wsave[iw1], &wsave[iw2]);

    return 0;
}


int cfftb1_(n, c, ch, wa, ifac)
integer *n;
doublereal *c, *ch, *wa;
integer *ifac;
{
    integer i__1;

    static integer idot, i;
    extern /* Subroutine */ int passb_();
    static integer k1, l1, l2, n2;
    extern /* Subroutine */ int passb2_(), passb3_(), passb4_(), passb5_();
    static integer na, nf, ip, iw, ix2, ix3, ix4, nac, ido, idl1;

    --ifac;
    --wa;
    --ch;
    --c;

    nf = ifac[2];
    na = 0;
    l1 = 1;
    iw = 1;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	l2 = ip * l1;
	ido = *n / l2;
	idot = ido + ido;
	idl1 = idot * l1;
	if (ip != 4) {
	    goto L103;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	if (na != 0) {
	    goto L101;
	}
	passb4_(&idot, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L102;
L101:
	passb4_(&idot, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2], &wa[ix3]);
L102:
	na = 1 - na;
	goto L115;

L103:
	if (ip != 2) {
	    goto L106;
	}
	if (na != 0) {
	    goto L104;
	}
	passb2_(&idot, &l1, &c[1], &ch[1], &wa[iw]);
	goto L105;
L104:
	passb2_(&idot, &l1, &ch[1], &c[1], &wa[iw]);
L105:
	na = 1 - na;
	goto L115;

L106:
	if (ip != 3) {
	    goto L109;
	}
	ix2 = iw + idot;
	if (na != 0) {
	    goto L107;
	}
	passb3_(&idot, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L108;
L107:
	passb3_(&idot, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2]);
L108:
	na = 1 - na;
	goto L115;

L109:
	if (ip != 5) {
	    goto L112;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	ix4 = ix3 + idot;
	if (na != 0) {
	    goto L110;
	}
	passb5_(&idot, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
	goto L111;
L110:
	passb5_(&idot, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
L111:
	na = 1 - na;
	goto L115;

L112:
	if (na != 0) {
	    goto L113;
	}
	passb_(&nac, &idot, &ip, &l1, &idl1, &c[1], &c[1], &c[1], &ch[1], &ch[
		1], &wa[iw]);
	goto L114;
L113:
	passb_(&nac, &idot, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c[1], &
		c[1], &wa[iw]);
L114:
	if (nac != 0) {
	    na = 1 - na;
	}

L115:
	l1 = l2;
	iw += (ip - 1) * idot;
/* L116: */
    }
    if (na == 0) {
	return 0;
    }

    n2 = *n + *n;
    i__1 = n2;
    for (i = 1; i <= i__1; ++i) {
	c[i] = ch[i];
/* L117: */
    }

    return 0;
}


int cfftf1_(n, c, ch, wa, ifac)
integer *n;
doublereal *c, *ch, *wa;
integer *ifac;
{
    integer i__1;

    static integer idot, i;
    extern /* Subroutine */ int passf_();
    static integer k1, l1, l2, n2;
    extern /* Subroutine */ int passf2_(), passf3_(), passf4_(), passf5_();
    static integer na, nf, ip, iw, ix2, ix3, ix4, nac, ido, idl1;

    --ifac;
    --wa;
    --ch;
    --c;

    nf = ifac[2];
    na = 0;
    l1 = 1;
    iw = 1;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	l2 = ip * l1;
	ido = *n / l2;
	idot = ido + ido;
	idl1 = idot * l1;
	if (ip != 4) {
	    goto L103;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	if (na != 0) {
	    goto L101;
	}
	passf4_(&idot, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L102;
L101:
	passf4_(&idot, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2], &wa[ix3]);
L102:
	na = 1 - na;
	goto L115;

L103:
	if (ip != 2) {
	    goto L106;
	}
	if (na != 0) {
	    goto L104;
	}
	passf2_(&idot, &l1, &c[1], &ch[1], &wa[iw]);
	goto L105;
L104:
	passf2_(&idot, &l1, &ch[1], &c[1], &wa[iw]);
L105:
	na = 1 - na;
	goto L115;

L106:
	if (ip != 3) {
	    goto L109;
	}
	ix2 = iw + idot;
	if (na != 0) {
	    goto L107;
	}
	passf3_(&idot, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L108;
L107:
	passf3_(&idot, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2]);
L108:
	na = 1 - na;
	goto L115;

L109:
	if (ip != 5) {
	    goto L112;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	ix4 = ix3 + idot;
	if (na != 0) {
	    goto L110;
	}
	passf5_(&idot, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
	goto L111;
L110:
	passf5_(&idot, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
L111:
	na = 1 - na;
	goto L115;

L112:
	if (na != 0) {
	    goto L113;
	}
	passf_(&nac, &idot, &ip, &l1, &idl1, &c[1], &c[1], &c[1], &ch[1], &ch[
		1], &wa[iw]);
	goto L114;
L113:
	passf_(&nac, &idot, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c[1], &
		c[1], &wa[iw]);
L114:
	if (nac != 0) {
	    na = 1 - na;
	}

L115:
	l1 = l2;
	iw += (ip - 1) * idot;
/* L116: */
    }
    if (na == 0) {
	return 0;
    }

    n2 = *n + *n;
    i__1 = n2;
    for (i = 1; i <= i__1; ++i) {
	c[i] = ch[i];
/* L117: */
    }

    return 0;
}


int cffti1_(n, wa, ifac)
integer *n;
doublereal *wa;
integer *ifac;
{

    static integer ntryh[4] = { 3,4,2,5 };
    static doublereal tpi = 6.28318530717958647692528676655900577;

    integer i__1, i__2, i__3;

    double cos(), sin();

    static doublereal argh;
    static integer idot, ntry, i, j;
    static doublereal argld;
    static integer i1, k1, l1, l2, ib;
    static doublereal fi;
    static integer ld, ii, nf, ip, nl, nq, nr;
    static doublereal arg;
    static integer ido, ipm;

    --ifac;
    --wa;

    nl = *n;
    nf = 0;
    j = 0;

L101:
    ++j;
    if (j <= 4) {
	ntry = ntryh[j - 1];
    }
    if (j > 4) {
	ntry += 2;
    }
L104:
    nq = nl / ntry;
    nr = nl - ntry * nq;
    if (nr != 0) {
	goto L101;
    }

/* L105: */
    ++nf;
    ifac[nf + 2] = ntry;
    nl = nq;
    if (ntry != 2) {
	goto L107;
    }
    if (nf == 1) {
	goto L107;
    }
    i__1 = nf;
    for (i = 2; i <= i__1; ++i) {
	ib = nf - i + 2;
	ifac[ib + 2] = ifac[ib + 1];
/* L106: */
    }
    ifac[3] = 2;

L107:
    if (nl != 1) {
	goto L104;
    }

    ifac[1] = *n;
    ifac[2] = nf;

    argh = tpi / (doublereal) (*n);
    i = 2;
    l1 = 1;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	ld = 0;
	l2 = l1 * ip;
	ido = *n / l2;
	idot = ido + ido + 2;
	ipm = ip - 1;

	i__2 = ipm;
	for (j = 1; j <= i__2; ++j) {
	    i1 = i;
	    wa[i - 1] = 1.;
	    wa[i] = 0.;
	    ld += l1;
	    fi = 0.;
	    argld = (doublereal) ld * argh;
	    i__3 = idot;
	    for (ii = 4; ii <= i__3; ii += 2) {
		i += 2;
		fi += 1.;
		arg = fi * argld;
		wa[i - 1] = cos(arg);
		wa[i] = sin(arg);
/* L108: */
	    }
	    if (ip <= 5) {
		goto L109;
	    }
	    wa[i1 - 1] = wa[i - 1];
	    wa[i1] = wa[i];
L109:
	    ;
	}

	l1 = l2;
/* L110: */
    }

    return 0;
}


int passb_(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
integer *nac, *ido, *ip, *l1, *idl1;
doublereal *cc, *c1, *c2, *ch, *ch2, *wa;
{
    integer cc_dim1, cc_dim2, cc_offset, c1_dim1, c1_dim2, c1_offset, c2_dim1,
	     c2_offset, ch_dim1, ch_dim2, ch_offset, ch2_dim1, ch2_offset, 
	    i__1, i__2, i__3;

    static integer idij, idlj, idot, ipph, i, j, k, l, jc, lc, ik, nt, idj, 
	    idl, inc, idp;
    static doublereal wai, war;
    static integer ipp2;

    --wa;
    ch2_dim1 = *idl1;
    ch2_offset = ch2_dim1 + 1;
    ch2 -= ch2_offset;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    c2_dim1 = *idl1;
    c2_offset = c2_dim1 + 1;
    c2 -= c2_offset;
    c1_dim1 = *ido;
    c1_dim2 = *l1;
    c1_offset = c1_dim1 * (c1_dim2 + 1) + 1;
    c1 -= c1_offset;
    cc_dim1 = *ido;
    cc_dim2 = *ip;
    cc_offset = cc_dim1 * (cc_dim2 + 1) + 1;
    cc -= cc_offset;

    idot = *ido / 2;
    nt = *ip * *idl1;
    ipp2 = *ip + 2;
    ipph = (*ip + 1) / 2;
    idp = *ip * *ido;

    if (*ido < *l1) {
	goto L106;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i = 1; i <= i__3; ++i) {
		ch[i + (k + j * ch_dim2) * ch_dim1] = cc[i + (j + k * cc_dim2)
			 * cc_dim1] + cc[i + (jc + k * cc_dim2) * cc_dim1];
		ch[i + (k + jc * ch_dim2) * ch_dim1] = cc[i + (j + k * 
			cc_dim2) * cc_dim1] - cc[i + (jc + k * cc_dim2) * 
			cc_dim1];
/* L101: */
	    }
/* L102: */
	}
/* L103: */
    }

    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 1; i <= i__2; ++i) {
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + (k * cc_dim2 + 1) * 
		    cc_dim1];
/* L104: */
	}
/* L105: */
    }
    goto L112;

L106:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i = 1; i <= i__2; ++i) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i + (k + j * ch_dim2) * ch_dim1] = cc[i + (j + k * cc_dim2)
			 * cc_dim1] + cc[i + (jc + k * cc_dim2) * cc_dim1];
		ch[i + (k + jc * ch_dim2) * ch_dim1] = cc[i + (j + k * 
			cc_dim2) * cc_dim1] - cc[i + (jc + k * cc_dim2) * 
			cc_dim1];
/* L107: */
	    }
/* L108: */
	}
/* L109: */
    }

    i__1 = *ido;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + (k * cc_dim2 + 1) * 
		    cc_dim1];
/* L110: */
	}
/* L111: */
    }

L112:
    idl = 2 - *ido;
    inc = 0;
    i__1 = ipph;
    for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	idl += *ido;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    c2[ik + l * c2_dim1] = ch2[ik + ch2_dim1] + wa[idl - 1] * ch2[ik 
		    + (ch2_dim1 << 1)];
	    c2[ik + lc * c2_dim1] = wa[idl] * ch2[ik + *ip * ch2_dim1];
/* L113: */
	}
	idlj = idl;
	inc += *ido;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    idlj += inc;
	    if (idlj > idp) {
		idlj -= idp;
	    }
	    war = wa[idlj - 1];
	    wai = wa[idlj];
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		c2[ik + l * c2_dim1] += war * ch2[ik + j * ch2_dim1];
		c2[ik + lc * c2_dim1] += wai * ch2[ik + jc * ch2_dim1];
/* L114: */
	    }
/* L115: */
	}
/* L116: */
    }

    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + ch2_dim1] += ch2[ik + j * ch2_dim1];
/* L117: */
	}
/* L118: */
    }

    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *idl1;
	for (ik = 2; ik <= i__2; ik += 2) {
	    ch2[ik - 1 + j * ch2_dim1] = c2[ik - 1 + j * c2_dim1] - c2[ik + 
		    jc * c2_dim1];
	    ch2[ik - 1 + jc * ch2_dim1] = c2[ik - 1 + j * c2_dim1] + c2[ik + 
		    jc * c2_dim1];
	    ch2[ik + j * ch2_dim1] = c2[ik + j * c2_dim1] + c2[ik - 1 + jc * 
		    c2_dim1];
	    ch2[ik + jc * ch2_dim1] = c2[ik + j * c2_dim1] - c2[ik - 1 + jc * 
		    c2_dim1];
/* L119: */
	}
/* L120: */
    }

    *nac = 1;
    if (*ido == 2) {
	return 0;
    }
    *nac = 0;

    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	c2[ik + c2_dim1] = ch2[ik + ch2_dim1];
/* L121: */
    }

    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    c1[(k + j * c1_dim2) * c1_dim1 + 1] = ch[(k + j * ch_dim2) * 
		    ch_dim1 + 1];
	    c1[(k + j * c1_dim2) * c1_dim1 + 2] = ch[(k + j * ch_dim2) * 
		    ch_dim1 + 2];
/* L122: */
	}
/* L123: */
    }

    if (idot > *l1) {
	goto L127;
    }
    idij = 0;
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	idij += 2;
	i__2 = *ido;
	for (i = 4; i <= i__2; i += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		c1[i - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i 
			- 1 + (k + j * ch_dim2) * ch_dim1] - wa[idij] * ch[i 
			+ (k + j * ch_dim2) * ch_dim1];
		c1[i + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i + (
			k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i - 1 + (
			k + j * ch_dim2) * ch_dim1];
/* L124: */
	    }
/* L125: */
	}
/* L126: */
    }
    return 0;

L127:
    idj = 2 - *ido;
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	idj += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = idj;
	    i__3 = *ido;
	    for (i = 4; i <= i__3; i += 2) {
		idij += 2;
		c1[i - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i 
			- 1 + (k + j * ch_dim2) * ch_dim1] - wa[idij] * ch[i 
			+ (k + j * ch_dim2) * ch_dim1];
		c1[i + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i + (
			k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i - 1 + (
			k + j * ch_dim2) * ch_dim1];
/* L128: */
	    }
/* L129: */
	}
/* L130: */
    }

    return 0;
}


int passb2_(ido, l1, cc, ch, wa1)
integer *ido, *l1;
doublereal *cc, *ch, *wa1;
{
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    static integer i, k;
    static doublereal ti2, tr2;

    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = cc_dim1 * 3 + 1;
    cc -= cc_offset;

    if (*ido > 2) {
	goto L102;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 + 1] + 
		cc[((k << 1) + 2) * cc_dim1 + 1];
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 
		+ 1] - cc[((k << 1) + 2) * cc_dim1 + 1];
	ch[(k + ch_dim2) * ch_dim1 + 2] = cc[((k << 1) + 1) * cc_dim1 + 2] + 
		cc[((k << 1) + 2) * cc_dim1 + 2];
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = cc[((k << 1) + 1) * cc_dim1 
		+ 2] - cc[((k << 1) + 2) * cc_dim1 + 2];
/* L101: */
    }
    return 0;

L102:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 2; i <= i__2; i += 2) {
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = cc[i - 1 + ((k << 1) + 1) * 
		    cc_dim1] + cc[i - 1 + ((k << 1) + 2) * cc_dim1];
	    tr2 = cc[i - 1 + ((k << 1) + 1) * cc_dim1] - cc[i - 1 + ((k << 1) 
		    + 2) * cc_dim1];
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + ((k << 1) + 1) * cc_dim1]
		     + cc[i + ((k << 1) + 2) * cc_dim1];
	    ti2 = cc[i + ((k << 1) + 1) * cc_dim1] - cc[i + ((k << 1) + 2) * 
		    cc_dim1];
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * ti2 + wa1[i]
		     * tr2;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * tr2 - 
		    wa1[i] * ti2;
/* L103: */
	}
/* L104: */
    }

    return 0;
}


int passb3_(ido, l1, cc, ch, wa1, wa2)
integer *ido, *l1;
doublereal *cc, *ch, *wa1, *wa2;
{
    static doublereal taur = -.5;
    static doublereal taui = .86602540378443864676372317075293618;

    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    static integer i, k;
    static doublereal ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;

    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = (cc_dim1 << 2) + 1;
    cc -= cc_offset;

/*     one half sqrt(3) = .866025.....  . */

    if (*ido != 2) {
	goto L102;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	tr2 = cc[(k * 3 + 2) * cc_dim1 + 1] + cc[(k * 3 + 3) * cc_dim1 + 1];
	cr2 = cc[(k * 3 + 1) * cc_dim1 + 1] + taur * tr2;
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 3 + 1) * cc_dim1 + 1] + tr2;
	ti2 = cc[(k * 3 + 2) * cc_dim1 + 2] + cc[(k * 3 + 3) * cc_dim1 + 2];
	ci2 = cc[(k * 3 + 1) * cc_dim1 + 2] + taur * ti2;
	ch[(k + ch_dim2) * ch_dim1 + 2] = cc[(k * 3 + 1) * cc_dim1 + 2] + ti2;
	cr3 = taui * (cc[(k * 3 + 2) * cc_dim1 + 1] - cc[(k * 3 + 3) * 
		cc_dim1 + 1]);
	ci3 = taui * (cc[(k * 3 + 2) * cc_dim1 + 2] - cc[(k * 3 + 3) * 
		cc_dim1 + 2]);
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr2 + ci3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ci2 + cr3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ci2 - cr3;
/* L101: */
    }
    return 0;

L102:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 2; i <= i__2; i += 2) {
	    tr2 = cc[i - 1 + (k * 3 + 2) * cc_dim1] + cc[i - 1 + (k * 3 + 3) *
		     cc_dim1];
	    cr2 = cc[i - 1 + (k * 3 + 1) * cc_dim1] + taur * tr2;
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = cc[i - 1 + (k * 3 + 1) * 
		    cc_dim1] + tr2;
	    ti2 = cc[i + (k * 3 + 2) * cc_dim1] + cc[i + (k * 3 + 3) * 
		    cc_dim1];
	    ci2 = cc[i + (k * 3 + 1) * cc_dim1] + taur * ti2;
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + (k * 3 + 1) * cc_dim1] + 
		    ti2;
	    cr3 = taui * (cc[i - 1 + (k * 3 + 2) * cc_dim1] - cc[i - 1 + (k * 
		    3 + 3) * cc_dim1]);
	    ci3 = taui * (cc[i + (k * 3 + 2) * cc_dim1] - cc[i + (k * 3 + 3) *
		     cc_dim1]);
	    dr2 = cr2 - ci3;
	    dr3 = cr2 + ci3;
	    di2 = ci2 + cr3;
	    di3 = ci2 - cr3;
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * di2 + wa1[i]
		     * dr2;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * dr2 - 
		    wa1[i] * di2;
	    ch[i + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * di3 + wa2[i] * 
		    dr3;
	    ch[i - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * dr3 - wa2[
		    i] * di3;
/* L103: */
	}
/* L104: */
    }

    return 0;
}


int passb4_(ido, l1, cc, ch, wa1, wa2, wa3)
integer *ido, *l1;
doublereal *cc, *ch, *wa1, *wa2, *wa3;
{
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    static integer i, k;
    static doublereal ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, 
	    tr2, tr3, tr4;

    --wa3;
    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = cc_dim1 * 5 + 1;
    cc -= cc_offset;

    if (*ido != 2) {
	goto L102;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ti1 = cc[((k << 2) + 1) * cc_dim1 + 2] - cc[((k << 2) + 3) * cc_dim1 
		+ 2];
	ti2 = cc[((k << 2) + 1) * cc_dim1 + 2] + cc[((k << 2) + 3) * cc_dim1 
		+ 2];
	tr4 = cc[((k << 2) + 4) * cc_dim1 + 2] - cc[((k << 2) + 2) * cc_dim1 
		+ 2];
	ti3 = cc[((k << 2) + 2) * cc_dim1 + 2] + cc[((k << 2) + 4) * cc_dim1 
		+ 2];
	tr1 = cc[((k << 2) + 1) * cc_dim1 + 1] - cc[((k << 2) + 3) * cc_dim1 
		+ 1];
	tr2 = cc[((k << 2) + 1) * cc_dim1 + 1] + cc[((k << 2) + 3) * cc_dim1 
		+ 1];
	ti4 = cc[((k << 2) + 2) * cc_dim1 + 1] - cc[((k << 2) + 4) * cc_dim1 
		+ 1];
	tr3 = cc[((k << 2) + 2) * cc_dim1 + 1] + cc[((k << 2) + 4) * cc_dim1 
		+ 1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = tr2 + tr3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = tr2 - tr3;
	ch[(k + ch_dim2) * ch_dim1 + 2] = ti2 + ti3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ti2 - ti3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = tr1 + tr4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = tr1 - tr4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ti1 + ti4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 2] = ti1 - ti4;
/* L101: */
    }
    return 0;

L102:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 2; i <= i__2; i += 2) {
	    ti1 = cc[i + ((k << 2) + 1) * cc_dim1] - cc[i + ((k << 2) + 3) * 
		    cc_dim1];
	    ti2 = cc[i + ((k << 2) + 1) * cc_dim1] + cc[i + ((k << 2) + 3) * 
		    cc_dim1];
	    ti3 = cc[i + ((k << 2) + 2) * cc_dim1] + cc[i + ((k << 2) + 4) * 
		    cc_dim1];
	    tr4 = cc[i + ((k << 2) + 4) * cc_dim1] - cc[i + ((k << 2) + 2) * 
		    cc_dim1];
	    tr1 = cc[i - 1 + ((k << 2) + 1) * cc_dim1] - cc[i - 1 + ((k << 2) 
		    + 3) * cc_dim1];
	    tr2 = cc[i - 1 + ((k << 2) + 1) * cc_dim1] + cc[i - 1 + ((k << 2) 
		    + 3) * cc_dim1];
	    ti4 = cc[i - 1 + ((k << 2) + 2) * cc_dim1] - cc[i - 1 + ((k << 2) 
		    + 4) * cc_dim1];
	    tr3 = cc[i - 1 + ((k << 2) + 2) * cc_dim1] + cc[i - 1 + ((k << 2) 
		    + 4) * cc_dim1];
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = tr2 + tr3;
	    cr3 = tr2 - tr3;
	    ch[i + (k + ch_dim2) * ch_dim1] = ti2 + ti3;
	    ci3 = ti2 - ti3;
	    cr2 = tr1 + tr4;
	    cr4 = tr1 - tr4;
	    ci2 = ti1 + ti4;
	    ci4 = ti1 - ti4;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * cr2 - 
		    wa1[i] * ci2;
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * ci2 + wa1[i]
		     * cr2;
	    ch[i - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * cr3 - wa2[
		    i] * ci3;
	    ch[i + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * ci3 + wa2[i] * 
		    cr3;
	    ch[i - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 1] * cr4 - 
		    wa3[i] * ci4;
	    ch[i + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 1] * ci4 + wa3[i]
		     * cr4;
/* L103: */
	}
/* L104: */
    }

    return 0;
}


int passb5_(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
integer *ido, *l1;
doublereal *cc, *ch, *wa1, *wa2, *wa3, *wa4;
{
    static doublereal tr11 = .30901699437494742410229341718281906;
    static doublereal ti11 = .95105651629515357211643933337938214;
    static doublereal tr12 = -.80901699437494742410229341718281906;
    static doublereal ti12 = .58778525229247312916870595463907277;

    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;
    doublereal d__1;

    static integer i, k;
    static doublereal ci2, ci3, ci4, ci5, di2, di3, di4, di5, cr2, cr3, cr4, 
	    cr5, dr2, dr3, dr4, dr5, ti2, ti3, ti4, ti5, tr2, tr3, tr4, tr5;

    --wa4;
    --wa3;
    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = cc_dim1 * 6 + 1;
    cc -= cc_offset;

/*     sin(pi/10) = .30901699....    . */
/*     cos(pi/10) = .95105651....    . */
/*     sin(pi/5 ) = .58778525....    . */
/*     cos(pi/5 ) = .80901699....    . */

    if (*ido != 2) {
	goto L102;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ti5 = cc[(k * 5 + 2) * cc_dim1 + 2] - cc[(k * 5 + 5) * cc_dim1 + 2];
	ti2 = cc[(k * 5 + 2) * cc_dim1 + 2] + cc[(k * 5 + 5) * cc_dim1 + 2];
	ti4 = cc[(k * 5 + 3) * cc_dim1 + 2] - cc[(k * 5 + 4) * cc_dim1 + 2];
	ti3 = cc[(k * 5 + 3) * cc_dim1 + 2] + cc[(k * 5 + 4) * cc_dim1 + 2];
	tr5 = cc[(k * 5 + 2) * cc_dim1 + 1] - cc[(k * 5 + 5) * cc_dim1 + 1];
	tr2 = cc[(k * 5 + 2) * cc_dim1 + 1] + cc[(k * 5 + 5) * cc_dim1 + 1];
	tr4 = cc[(k * 5 + 3) * cc_dim1 + 1] - cc[(k * 5 + 4) * cc_dim1 + 1];
	tr3 = cc[(k * 5 + 3) * cc_dim1 + 1] + cc[(k * 5 + 4) * cc_dim1 + 1];
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr2;
	ch[(k + ch_dim2) * ch_dim1 + 1] = d__1 + tr3;
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 2] + ti2;
	ch[(k + ch_dim2) * ch_dim1 + 2] = d__1 + ti3;
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr11 * tr2;
	cr2 = d__1 + tr12 * tr3;
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 2] + tr11 * ti2;
	ci2 = d__1 + tr12 * ti3;
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr12 * tr2;
	cr3 = d__1 + tr11 * tr3;
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 2] + tr12 * ti2;
	ci3 = d__1 + tr11 * ti3;
	cr5 = ti11 * tr5 + ti12 * tr4;
	ci5 = ti11 * ti5 + ti12 * ti4;
	cr4 = ti12 * tr5 - ti11 * tr4;
	ci4 = ti12 * ti5 - ti11 * ti4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci5;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 1] = cr2 + ci5;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ci2 + cr5;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ci3 + cr4;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr3 - ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = cr3 + ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 2] = ci3 - cr4;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 2] = ci2 - cr5;
/* L101: */
    }
    return 0;

L102:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 2; i <= i__2; i += 2) {
	    ti5 = cc[i + (k * 5 + 2) * cc_dim1] - cc[i + (k * 5 + 5) * 
		    cc_dim1];
	    ti2 = cc[i + (k * 5 + 2) * cc_dim1] + cc[i + (k * 5 + 5) * 
		    cc_dim1];
	    ti4 = cc[i + (k * 5 + 3) * cc_dim1] - cc[i + (k * 5 + 4) * 
		    cc_dim1];
	    ti3 = cc[i + (k * 5 + 3) * cc_dim1] + cc[i + (k * 5 + 4) * 
		    cc_dim1];
	    tr5 = cc[i - 1 + (k * 5 + 2) * cc_dim1] - cc[i - 1 + (k * 5 + 5) *
		     cc_dim1];
	    tr2 = cc[i - 1 + (k * 5 + 2) * cc_dim1] + cc[i - 1 + (k * 5 + 5) *
		     cc_dim1];
	    tr4 = cc[i - 1 + (k * 5 + 3) * cc_dim1] - cc[i - 1 + (k * 5 + 4) *
		     cc_dim1];
	    tr3 = cc[i - 1 + (k * 5 + 3) * cc_dim1] + cc[i - 1 + (k * 5 + 4) *
		     cc_dim1];
	    d__1 = cc[i - 1 + (k * 5 + 1) * cc_dim1] + tr2;
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = d__1 + tr3;
	    d__1 = cc[i + (k * 5 + 1) * cc_dim1] + ti2;
	    ch[i + (k + ch_dim2) * ch_dim1] = d__1 + ti3;
	    d__1 = cc[i - 1 + (k * 5 + 1) * cc_dim1] + tr11 * tr2;
	    cr2 = d__1 + tr12 * tr3;
	    d__1 = cc[i + (k * 5 + 1) * cc_dim1] + tr11 * ti2;
	    ci2 = d__1 + tr12 * ti3;
	    d__1 = cc[i - 1 + (k * 5 + 1) * cc_dim1] + tr12 * tr2;
	    cr3 = d__1 + tr11 * tr3;
	    d__1 = cc[i + (k * 5 + 1) * cc_dim1] + tr12 * ti2;
	    ci3 = d__1 + tr11 * ti3;
	    cr5 = ti11 * tr5 + ti12 * tr4;
	    ci5 = ti11 * ti5 + ti12 * ti4;
	    cr4 = ti12 * tr5 - ti11 * tr4;
	    ci4 = ti12 * ti5 - ti11 * ti4;
	    dr3 = cr3 - ci4;
	    dr4 = cr3 + ci4;
	    di3 = ci3 + cr4;
	    di4 = ci3 - cr4;
	    dr5 = cr2 + ci5;
	    dr2 = cr2 - ci5;
	    di5 = ci2 - cr5;
	    di2 = ci2 + cr5;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * dr2 - 
		    wa1[i] * di2;
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * di2 + wa1[i]
		     * dr2;
	    ch[i - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * dr3 - wa2[
		    i] * di3;
	    ch[i + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * di3 + wa2[i] * 
		    dr3;
	    ch[i - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 1] * dr4 - 
		    wa3[i] * di4;
	    ch[i + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 1] * di4 + wa3[i]
		     * dr4;
	    ch[i - 1 + (k + ch_dim2 * 5) * ch_dim1] = wa4[i - 1] * dr5 - wa4[
		    i] * di5;
	    ch[i + (k + ch_dim2 * 5) * ch_dim1] = wa4[i - 1] * di5 + wa4[i] * 
		    dr5;
/* L103: */
	}
/* L104: */
    }

    return 0;
}


int passf_(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
integer *nac, *ido, *ip, *l1, *idl1;
doublereal *cc, *c1, *c2, *ch, *ch2, *wa;
{
    integer cc_dim1, cc_dim2, cc_offset, c1_dim1, c1_dim2, c1_offset, c2_dim1,
	     c2_offset, ch_dim1, ch_dim2, ch_offset, ch2_dim1, ch2_offset, 
	    i__1, i__2, i__3;

    static integer idij, idlj, idot, ipph, i, j, k, l, jc, lc, ik, nt, idj, 
	    idl, inc, idp;
    static doublereal wai, war;
    static integer ipp2;

    --wa;
    ch2_dim1 = *idl1;
    ch2_offset = ch2_dim1 + 1;
    ch2 -= ch2_offset;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    c2_dim1 = *idl1;
    c2_offset = c2_dim1 + 1;
    c2 -= c2_offset;
    c1_dim1 = *ido;
    c1_dim2 = *l1;
    c1_offset = c1_dim1 * (c1_dim2 + 1) + 1;
    c1 -= c1_offset;
    cc_dim1 = *ido;
    cc_dim2 = *ip;
    cc_offset = cc_dim1 * (cc_dim2 + 1) + 1;
    cc -= cc_offset;

    idot = *ido / 2;
    nt = *ip * *idl1;
    ipp2 = *ip + 2;
    ipph = (*ip + 1) / 2;
    idp = *ip * *ido;

    if (*ido < *l1) {
	goto L106;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i = 1; i <= i__3; ++i) {
		ch[i + (k + j * ch_dim2) * ch_dim1] = cc[i + (j + k * cc_dim2)
			 * cc_dim1] + cc[i + (jc + k * cc_dim2) * cc_dim1];
		ch[i + (k + jc * ch_dim2) * ch_dim1] = cc[i + (j + k * 
			cc_dim2) * cc_dim1] - cc[i + (jc + k * cc_dim2) * 
			cc_dim1];
/* L101: */
	    }
/* L102: */
	}
/* L103: */
    }

    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 1; i <= i__2; ++i) {
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + (k * cc_dim2 + 1) * 
		    cc_dim1];
/* L104: */
	}
/* L105: */
    }
    goto L112;

L106:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i = 1; i <= i__2; ++i) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i + (k + j * ch_dim2) * ch_dim1] = cc[i + (j + k * cc_dim2)
			 * cc_dim1] + cc[i + (jc + k * cc_dim2) * cc_dim1];
		ch[i + (k + jc * ch_dim2) * ch_dim1] = cc[i + (j + k * 
			cc_dim2) * cc_dim1] - cc[i + (jc + k * cc_dim2) * 
			cc_dim1];
/* L107: */
	    }
/* L108: */
	}
/* L109: */
    }

    i__1 = *ido;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + (k * cc_dim2 + 1) * 
		    cc_dim1];
/* L110: */
	}
/* L111: */
    }

L112:
    idl = 2 - *ido;
    inc = 0;
    i__1 = ipph;
    for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	idl += *ido;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    c2[ik + l * c2_dim1] = ch2[ik + ch2_dim1] + wa[idl - 1] * ch2[ik 
		    + (ch2_dim1 << 1)];
	    c2[ik + lc * c2_dim1] = -wa[idl] * ch2[ik + *ip * ch2_dim1];
/* L113: */
	}
	idlj = idl;
	inc += *ido;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    idlj += inc;
	    if (idlj > idp) {
		idlj -= idp;
	    }
	    war = wa[idlj - 1];
	    wai = wa[idlj];
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		c2[ik + l * c2_dim1] += war * ch2[ik + j * ch2_dim1];
		c2[ik + lc * c2_dim1] -= wai * ch2[ik + jc * ch2_dim1];
/* L114: */
	    }
/* L115: */
	}
/* L116: */
    }

    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + ch2_dim1] += ch2[ik + j * ch2_dim1];
/* L117: */
	}
/* L118: */
    }

    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *idl1;
	for (ik = 2; ik <= i__2; ik += 2) {
	    ch2[ik - 1 + j * ch2_dim1] = c2[ik - 1 + j * c2_dim1] - c2[ik + 
		    jc * c2_dim1];
	    ch2[ik - 1 + jc * ch2_dim1] = c2[ik - 1 + j * c2_dim1] + c2[ik + 
		    jc * c2_dim1];
	    ch2[ik + j * ch2_dim1] = c2[ik + j * c2_dim1] + c2[ik - 1 + jc * 
		    c2_dim1];
	    ch2[ik + jc * ch2_dim1] = c2[ik + j * c2_dim1] - c2[ik - 1 + jc * 
		    c2_dim1];
/* L119: */
	}
/* L120: */
    }

    *nac = 1;
    if (*ido == 2) {
	return 0;
    }
    *nac = 0;

    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	c2[ik + c2_dim1] = ch2[ik + ch2_dim1];
/* L121: */
    }

    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    c1[(k + j * c1_dim2) * c1_dim1 + 1] = ch[(k + j * ch_dim2) * 
		    ch_dim1 + 1];
	    c1[(k + j * c1_dim2) * c1_dim1 + 2] = ch[(k + j * ch_dim2) * 
		    ch_dim1 + 2];
/* L122: */
	}
/* L123: */
    }

    if (idot > *l1) {
	goto L127;
    }
    idij = 0;
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	idij += 2;
	i__2 = *ido;
	for (i = 4; i <= i__2; i += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		c1[i - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i 
			- 1 + (k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i 
			+ (k + j * ch_dim2) * ch_dim1];
		c1[i + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i + (
			k + j * ch_dim2) * ch_dim1] - wa[idij] * ch[i - 1 + (
			k + j * ch_dim2) * ch_dim1];
/* L124: */
	    }
/* L125: */
	}
/* L126: */
    }
    return 0;

L127:
    idj = 2 - *ido;
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	idj += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = idj;
	    i__3 = *ido;
	    for (i = 4; i <= i__3; i += 2) {
		idij += 2;
		c1[i - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i 
			- 1 + (k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i 
			+ (k + j * ch_dim2) * ch_dim1];
		c1[i + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i + (
			k + j * ch_dim2) * ch_dim1] - wa[idij] * ch[i - 1 + (
			k + j * ch_dim2) * ch_dim1];
/* L128: */
	    }
/* L129: */
	}
/* L130: */
    }

    return 0;
}


int passf2_(ido, l1, cc, ch, wa1)
integer *ido, *l1;
doublereal *cc, *ch, *wa1;
{
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    static integer i, k;
    static doublereal ti2, tr2;

    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = cc_dim1 * 3 + 1;
    cc -= cc_offset;

    if (*ido > 2) {
	goto L102;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 + 1] + 
		cc[((k << 1) + 2) * cc_dim1 + 1];
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 
		+ 1] - cc[((k << 1) + 2) * cc_dim1 + 1];
	ch[(k + ch_dim2) * ch_dim1 + 2] = cc[((k << 1) + 1) * cc_dim1 + 2] + 
		cc[((k << 1) + 2) * cc_dim1 + 2];
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = cc[((k << 1) + 1) * cc_dim1 
		+ 2] - cc[((k << 1) + 2) * cc_dim1 + 2];
/* L101: */
    }
    return 0;

L102:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 2; i <= i__2; i += 2) {
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = cc[i - 1 + ((k << 1) + 1) * 
		    cc_dim1] + cc[i - 1 + ((k << 1) + 2) * cc_dim1];
	    tr2 = cc[i - 1 + ((k << 1) + 1) * cc_dim1] - cc[i - 1 + ((k << 1) 
		    + 2) * cc_dim1];
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + ((k << 1) + 1) * cc_dim1]
		     + cc[i + ((k << 1) + 2) * cc_dim1];
	    ti2 = cc[i + ((k << 1) + 1) * cc_dim1] - cc[i + ((k << 1) + 2) * 
		    cc_dim1];
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * ti2 - wa1[i]
		     * tr2;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * tr2 + 
		    wa1[i] * ti2;
/* L103: */
	}
/* L104: */
    }

    return 0;
}


int passf3_(ido, l1, cc, ch, wa1, wa2)
integer *ido, *l1;
doublereal *cc, *ch, *wa1, *wa2;
{
    static doublereal taur = -.5;
    static doublereal taui = -.86602540378443864676372317075293618;

    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    static integer i, k;
    static doublereal ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;

    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = (cc_dim1 << 2) + 1;
    cc -= cc_offset;

    if (*ido != 2) {
	goto L102;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	tr2 = cc[(k * 3 + 2) * cc_dim1 + 1] + cc[(k * 3 + 3) * cc_dim1 + 1];
	cr2 = cc[(k * 3 + 1) * cc_dim1 + 1] + taur * tr2;
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 3 + 1) * cc_dim1 + 1] + tr2;
	ti2 = cc[(k * 3 + 2) * cc_dim1 + 2] + cc[(k * 3 + 3) * cc_dim1 + 2];
	ci2 = cc[(k * 3 + 1) * cc_dim1 + 2] + taur * ti2;
	ch[(k + ch_dim2) * ch_dim1 + 2] = cc[(k * 3 + 1) * cc_dim1 + 2] + ti2;
	cr3 = taui * (cc[(k * 3 + 2) * cc_dim1 + 1] - cc[(k * 3 + 3) * 
		cc_dim1 + 1]);
	ci3 = taui * (cc[(k * 3 + 2) * cc_dim1 + 2] - cc[(k * 3 + 3) * 
		cc_dim1 + 2]);
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr2 + ci3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ci2 + cr3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ci2 - cr3;
/* L101: */
    }
    return 0;

L102:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 2; i <= i__2; i += 2) {
	    tr2 = cc[i - 1 + (k * 3 + 2) * cc_dim1] + cc[i - 1 + (k * 3 + 3) *
		     cc_dim1];
	    cr2 = cc[i - 1 + (k * 3 + 1) * cc_dim1] + taur * tr2;
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = cc[i - 1 + (k * 3 + 1) * 
		    cc_dim1] + tr2;
	    ti2 = cc[i + (k * 3 + 2) * cc_dim1] + cc[i + (k * 3 + 3) * 
		    cc_dim1];
	    ci2 = cc[i + (k * 3 + 1) * cc_dim1] + taur * ti2;
	    ch[i + (k + ch_dim2) * ch_dim1] = cc[i + (k * 3 + 1) * cc_dim1] + 
		    ti2;
	    cr3 = taui * (cc[i - 1 + (k * 3 + 2) * cc_dim1] - cc[i - 1 + (k * 
		    3 + 3) * cc_dim1]);
	    ci3 = taui * (cc[i + (k * 3 + 2) * cc_dim1] - cc[i + (k * 3 + 3) *
		     cc_dim1]);
	    dr2 = cr2 - ci3;
	    dr3 = cr2 + ci3;
	    di2 = ci2 + cr3;
	    di3 = ci2 - cr3;
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * di2 - wa1[i]
		     * dr2;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * dr2 + 
		    wa1[i] * di2;
	    ch[i + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * di3 - wa2[i] * 
		    dr3;
	    ch[i - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * dr3 + wa2[
		    i] * di3;
/* L103: */
	}
/* L104: */
    }

    return 0;
}


int passf4_(ido, l1, cc, ch, wa1, wa2, wa3)
integer *ido, *l1;
doublereal *cc, *ch, *wa1, *wa2, *wa3;
{
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    static integer i, k;
    static doublereal ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, 
	    tr2, tr3, tr4;

    --wa3;
    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = cc_dim1 * 5 + 1;
    cc -= cc_offset;

    if (*ido != 2) {
	goto L102;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ti1 = cc[((k << 2) + 1) * cc_dim1 + 2] - cc[((k << 2) + 3) * cc_dim1 
		+ 2];
	ti2 = cc[((k << 2) + 1) * cc_dim1 + 2] + cc[((k << 2) + 3) * cc_dim1 
		+ 2];
	tr4 = cc[((k << 2) + 2) * cc_dim1 + 2] - cc[((k << 2) + 4) * cc_dim1 
		+ 2];
	ti3 = cc[((k << 2) + 2) * cc_dim1 + 2] + cc[((k << 2) + 4) * cc_dim1 
		+ 2];
	tr1 = cc[((k << 2) + 1) * cc_dim1 + 1] - cc[((k << 2) + 3) * cc_dim1 
		+ 1];
	tr2 = cc[((k << 2) + 1) * cc_dim1 + 1] + cc[((k << 2) + 3) * cc_dim1 
		+ 1];
	ti4 = cc[((k << 2) + 4) * cc_dim1 + 1] - cc[((k << 2) + 2) * cc_dim1 
		+ 1];
	tr3 = cc[((k << 2) + 2) * cc_dim1 + 1] + cc[((k << 2) + 4) * cc_dim1 
		+ 1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = tr2 + tr3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = tr2 - tr3;
	ch[(k + ch_dim2) * ch_dim1 + 2] = ti2 + ti3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ti2 - ti3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = tr1 + tr4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = tr1 - tr4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ti1 + ti4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 2] = ti1 - ti4;
/* L101: */
    }
    return 0;

L102:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 2; i <= i__2; i += 2) {
	    ti1 = cc[i + ((k << 2) + 1) * cc_dim1] - cc[i + ((k << 2) + 3) * 
		    cc_dim1];
	    ti2 = cc[i + ((k << 2) + 1) * cc_dim1] + cc[i + ((k << 2) + 3) * 
		    cc_dim1];
	    ti3 = cc[i + ((k << 2) + 2) * cc_dim1] + cc[i + ((k << 2) + 4) * 
		    cc_dim1];
	    tr4 = cc[i + ((k << 2) + 2) * cc_dim1] - cc[i + ((k << 2) + 4) * 
		    cc_dim1];
	    tr1 = cc[i - 1 + ((k << 2) + 1) * cc_dim1] - cc[i - 1 + ((k << 2) 
		    + 3) * cc_dim1];
	    tr2 = cc[i - 1 + ((k << 2) + 1) * cc_dim1] + cc[i - 1 + ((k << 2) 
		    + 3) * cc_dim1];
	    ti4 = cc[i - 1 + ((k << 2) + 4) * cc_dim1] - cc[i - 1 + ((k << 2) 
		    + 2) * cc_dim1];
	    tr3 = cc[i - 1 + ((k << 2) + 2) * cc_dim1] + cc[i - 1 + ((k << 2) 
		    + 4) * cc_dim1];
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = tr2 + tr3;
	    cr3 = tr2 - tr3;
	    ch[i + (k + ch_dim2) * ch_dim1] = ti2 + ti3;
	    ci3 = ti2 - ti3;
	    cr2 = tr1 + tr4;
	    cr4 = tr1 - tr4;
	    ci2 = ti1 + ti4;
	    ci4 = ti1 - ti4;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * cr2 + 
		    wa1[i] * ci2;
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * ci2 - wa1[i]
		     * cr2;
	    ch[i - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * cr3 + wa2[
		    i] * ci3;
	    ch[i + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * ci3 - wa2[i] * 
		    cr3;
	    ch[i - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 1] * cr4 + 
		    wa3[i] * ci4;
	    ch[i + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 1] * ci4 - wa3[i]
		     * cr4;
/* L103: */
	}
/* L104: */
    }

    return 0;
}


int passf5_(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
integer *ido, *l1;
doublereal *cc, *ch, *wa1, *wa2, *wa3, *wa4;
{
    static doublereal tr11 = .30901699437494742410229341718281906;
    static doublereal ti11 = -.95105651629515357211643933337938214;
    static doublereal tr12 = -.80901699437494742410229341718281906;
    static doublereal ti12 = -.58778525229247312916870595463907277;

    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;
    doublereal d__1;

    static integer i, k;
    static doublereal ci2, ci3, ci4, ci5, di2, di3, di4, di5, cr2, cr3, cr4, 
	    cr5, dr2, dr3, dr4, dr5, ti2, ti3, ti4, ti5, tr2, tr3, tr4, tr5;

    --wa4;
    --wa3;
    --wa2;
    --wa1;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = cc_dim1 * 6 + 1;
    cc -= cc_offset;

    if (*ido != 2) {
	goto L102;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ti5 = cc[(k * 5 + 2) * cc_dim1 + 2] - cc[(k * 5 + 5) * cc_dim1 + 2];
	ti2 = cc[(k * 5 + 2) * cc_dim1 + 2] + cc[(k * 5 + 5) * cc_dim1 + 2];
	ti4 = cc[(k * 5 + 3) * cc_dim1 + 2] - cc[(k * 5 + 4) * cc_dim1 + 2];
	ti3 = cc[(k * 5 + 3) * cc_dim1 + 2] + cc[(k * 5 + 4) * cc_dim1 + 2];
	tr5 = cc[(k * 5 + 2) * cc_dim1 + 1] - cc[(k * 5 + 5) * cc_dim1 + 1];
	tr2 = cc[(k * 5 + 2) * cc_dim1 + 1] + cc[(k * 5 + 5) * cc_dim1 + 1];
	tr4 = cc[(k * 5 + 3) * cc_dim1 + 1] - cc[(k * 5 + 4) * cc_dim1 + 1];
	tr3 = cc[(k * 5 + 3) * cc_dim1 + 1] + cc[(k * 5 + 4) * cc_dim1 + 1];
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr2;
	ch[(k + ch_dim2) * ch_dim1 + 1] = d__1 + tr3;
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 2] + ti2;
	ch[(k + ch_dim2) * ch_dim1 + 2] = d__1 + ti3;
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr11 * tr2;
	cr2 = d__1 + tr12 * tr3;
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 2] + tr11 * ti2;
	ci2 = d__1 + tr12 * ti3;
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr12 * tr2;
	cr3 = d__1 + tr11 * tr3;
	d__1 = cc[(k * 5 + 1) * cc_dim1 + 2] + tr12 * ti2;
	ci3 = d__1 + tr11 * ti3;
	cr5 = ti11 * tr5 + ti12 * tr4;
	ci5 = ti11 * ti5 + ti12 * ti4;
	cr4 = ti12 * tr5 - ti11 * tr4;
	ci4 = ti12 * ti5 - ti11 * ti4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci5;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 1] = cr2 + ci5;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ci2 + cr5;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ci3 + cr4;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr3 - ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = cr3 + ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 2] = ci3 - cr4;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 2] = ci2 - cr5;
/* L101: */
    }
    return 0;

L102:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i = 2; i <= i__2; i += 2) {
	    ti5 = cc[i + (k * 5 + 2) * cc_dim1] - cc[i + (k * 5 + 5) * 
		    cc_dim1];
	    ti2 = cc[i + (k * 5 + 2) * cc_dim1] + cc[i + (k * 5 + 5) * 
		    cc_dim1];
	    ti4 = cc[i + (k * 5 + 3) * cc_dim1] - cc[i + (k * 5 + 4) * 
		    cc_dim1];
	    ti3 = cc[i + (k * 5 + 3) * cc_dim1] + cc[i + (k * 5 + 4) * 
		    cc_dim1];
	    tr5 = cc[i - 1 + (k * 5 + 2) * cc_dim1] - cc[i - 1 + (k * 5 + 5) *
		     cc_dim1];
	    tr2 = cc[i - 1 + (k * 5 + 2) * cc_dim1] + cc[i - 1 + (k * 5 + 5) *
		     cc_dim1];
	    tr4 = cc[i - 1 + (k * 5 + 3) * cc_dim1] - cc[i - 1 + (k * 5 + 4) *
		     cc_dim1];
	    tr3 = cc[i - 1 + (k * 5 + 3) * cc_dim1] + cc[i - 1 + (k * 5 + 4) *
		     cc_dim1];
	    d__1 = cc[i - 1 + (k * 5 + 1) * cc_dim1] + tr2;
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = d__1 + tr3;
	    d__1 = cc[i + (k * 5 + 1) * cc_dim1] + ti2;
	    ch[i + (k + ch_dim2) * ch_dim1] = d__1 + ti3;
	    d__1 = cc[i - 1 + (k * 5 + 1) * cc_dim1] + tr11 * tr2;
	    cr2 = d__1 + tr12 * tr3;
	    d__1 = cc[i + (k * 5 + 1) * cc_dim1] + tr11 * ti2;
	    ci2 = d__1 + tr12 * ti3;
	    d__1 = cc[i - 1 + (k * 5 + 1) * cc_dim1] + tr12 * tr2;
	    cr3 = d__1 + tr11 * tr3;
	    d__1 = cc[i + (k * 5 + 1) * cc_dim1] + tr12 * ti2;
	    ci3 = d__1 + tr11 * ti3;
	    cr5 = ti11 * tr5 + ti12 * tr4;
	    ci5 = ti11 * ti5 + ti12 * ti4;
	    cr4 = ti12 * tr5 - ti11 * tr4;
	    ci4 = ti12 * ti5 - ti11 * ti4;
	    dr3 = cr3 - ci4;
	    dr4 = cr3 + ci4;
	    di3 = ci3 + cr4;
	    di4 = ci3 - cr4;
	    dr5 = cr2 + ci5;
	    dr2 = cr2 - ci5;
	    di5 = ci2 - cr5;
	    di2 = ci2 + cr5;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * dr2 + 
		    wa1[i] * di2;
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * di2 - wa1[i]
		     * dr2;
	    ch[i - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * dr3 + wa2[
		    i] * di3;
	    ch[i + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * di3 - wa2[i] * 
		    dr3;
	    ch[i - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 1] * dr4 + 
		    wa3[i] * di4;
	    ch[i + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 1] * di4 - wa3[i]
		     * dr4;
	    ch[i - 1 + (k + ch_dim2 * 5) * ch_dim1] = wa4[i - 1] * dr5 + wa4[
		    i] * di5;
	    ch[i + (k + ch_dim2 * 5) * ch_dim1] = wa4[i - 1] * di5 - wa4[i] * 
		    dr5;
/* L103: */
	}
/* L104: */
    }

    return 0;
}
