/* diagq.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"



/*     ################################################### */
/*     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ## */
/*     ##              All Rights Reserved              ## */
/*     ################################################### */

/*     ################################################################# */
/*     ##                                                             ## */
/*     ##  subroutine diagq  --  fast matrix diagonalization routine  ## */
/*     ##                                                             ## */
/*     ################################################################# */


/*     "diagq" is a matrix diagonalization routine which is derived */
/*     from the classical given, housec, and eigen algorithms with */
/*     several modifications to increase the efficiency and accuracy */

/*     variables and parameters: */

/*     n         logical dimension of the matrix to be diagonalized */
/*     np        physical dimension of the matrix storage area */
/*     nv        number of eigenvalues and eigenvectors desired */
/*     dd        upper triangle of the matrix to be diagonalized */
/*     ev        returned with the eigenvalues in ascending order */
/*     vec       returned with the eigenvectors of the matrix */
/*     a,b,p,w   some temporary work vectors of physical dimension */
/*     ta,tb,y   more temporary work vectors of physical dimension */

/*     literature reference: */

/*     adapted from an original program written by Bernie Brooks, */
/*     National Institutes of Health, Bethesda, MD */


/* Subroutine */ int diagq_(integer *n, integer *np, integer *nv, doublereal *
	dd, doublereal *ev, doublereal *vec, doublereal *a, doublereal *b, 
	doublereal *p, doublereal *w, doublereal *ta, doublereal *tb, 
	doublereal *y)
{
    /* System generated locals */
    integer vec_dim1, vec_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), d_mod(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal s, t, f0;
    static integer j1, ia, ii, ji, mi, mj, mk;
    static doublereal bx;
    static integer nn, mi1, mj1, mk1;
    static doublereal eta, epr, eps, rho;
    static integer nom;
    static doublereal sgn;
    static integer ipt;
    static doublereal sum1;
    static logical done;
    static integer iter;
    static doublereal zeta, temp, term, root;
    static integer ntot;
    static doublereal rand1, elim1, elim2, rpow1, gamma, delta, sigma, theta, 
	    trial, anorm, xnorm, rootx, factor, alimit;
    static integer nomtch;
    static doublereal rvalue, rpower;


#define vec_ref(a_1,a_2) vec[(a_2)*vec_dim1 + a_1]



/*     initialization of some necessary parameters */

    /* Parameter adjustments */
    --y;
    --tb;
    --ta;
    --w;
    --p;
    --b;
    --a;
    vec_dim1 = *np;
    vec_offset = 1 + vec_dim1;
    vec -= vec_offset;
    --ev;
    --dd;

    /* Function Body */
    eta = 1e-16;
    theta = 1e37;
    eps = eta * 100.;
    rho = eta / 100.;
/* Computing 2nd power */
    d__1 = eta;
    delta = d__1 * d__1 * 100.;
/* Computing 2nd power */
    d__1 = eta;
    gamma = d__1 * d__1 / 100.;
    zeta = 1e3 / theta;
    sigma = theta * delta / 1e3;
    rvalue = 4099.;
    rpower = 8388608.;
    rpow1 = rpower * .5;
    rand1 = rpower - 3.;

/*     get norm of the input matrix and scale */

    factor = 0.;
    ntot = *n * (*n + 1) / 2;
    i__1 = ntot;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = factor, d__3 = (d__1 = dd[i__], abs(d__1));
	factor = max(d__2,d__3);
    }
    if (factor == 0.) {
	return 0;
    }
    k = 0;
    anorm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    ++k;
/* Computing 2nd power */
	    d__1 = dd[k] / factor;
	    term = d__1 * d__1;
	    if (i__ == j) {
		term *= .5;
	    }
	    anorm += term;
	}
    }
    anorm = factor * sqrt(anorm * 2.);
    i__1 = ntot;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dd[i__] /= anorm;
    }

/*     compute the tridiagonalization of the matrix */

    nn = *n - 1;
    mi = 0;
    mi1 = *n - 1;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum1 = 0.;
	b[i__] = 0.;
	ji = i__ + 1;
	ipt = mi + i__;
	a[i__] = dd[ipt];
	++ipt;
	bx = dd[ipt];
	i__2 = *n;
	for (j = ji + 1; j <= i__2; ++j) {
	    ++ipt;
	    sum1 += dd[ipt] * dd[ipt];
	}
	if (sum1 < gamma) {
	    b[i__] = bx;
	    dd[mi + ji] = 0.;
	} else {
	    s = sqrt(sum1 + bx * bx);
	    sgn = 1.;
	    if (bx < 0.f) {
		sgn = -1.;
	    }
	    temp = abs(bx);
	    w[ji] = sqrt((temp / s + 1.) * .5);
	    ipt = mi + ji;
	    dd[ipt] = w[ji];
	    ii = i__ + 2;
	    if (ii <= *n) {
		temp = sgn / (w[ji] * 2. * s);
		i__2 = *n;
		for (j = ii; j <= i__2; ++j) {
		    ++ipt;
		    w[j] = temp * dd[ipt];
		    dd[ipt] = w[j];
		}
	    }
	    b[i__] = -sgn * s;
	    i__2 = *n;
	    for (j = ji; j <= i__2; ++j) {
		p[j] = 0.;
	    }
	    mk = mi + mi1;
	    mk1 = mi1 - 1;
	    i__2 = *n;
	    for (k = ji; k <= i__2; ++k) {
		ipt = mk + k;
		i__3 = *n;
		for (m = k; m <= i__3; ++m) {
		    bx = dd[ipt];
		    p[k] += bx * w[m];
		    if (k != m) {
			p[m] += bx * w[k];
		    }
		    ++ipt;
		}
		mk += mk1;
		--mk1;
	    }
	    term = 0.;
	    i__2 = *n;
	    for (k = ji; k <= i__2; ++k) {
		term += w[k] * p[k];
	    }
	    i__2 = *n;
	    for (k = ji; k <= i__2; ++k) {
		p[k] -= term * w[k];
	    }
	    mj = mi + mi1;
	    mj1 = mi1 - 1;
	    i__2 = *n;
	    for (j = ji; j <= i__2; ++j) {
		i__3 = *n;
		for (k = j; k <= i__3; ++k) {
		    dd[mj + k] -= (p[j] * w[k] + p[k] * w[j]) * 2.;
		}
		mj += mj1;
		--mj1;
	    }
	}
	mi += mi1;
	--mi1;
    }

/*     find the eigenvalues via the Sturm bisection method */

    a[*n] = dd[mi + *n];
    b[*n] = 0.;
    alimit = 1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = b[i__];
	b[i__] *= b[i__];
    }
    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ev[i__] = alimit;
    }
    root = -alimit;
    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rootx = alimit;
	i__2 = *nv;
	for (j = i__; j <= i__2; ++j) {
/* Computing MIN */
	    d__1 = rootx, d__2 = ev[j];
	    rootx = min(d__1,d__2);
	}
	ev[i__] = rootx;
	trial = (root + ev[i__]) * .5;
	while((d__1 = trial - root, abs(d__1)) >= eps && (d__2 = trial - ev[
		i__], abs(d__2)) >= eps) {
	    nomtch = *n;
	    j = 1;
	    while(j <= *n) {
		f0 = a[j] - trial;
		while(abs(f0) >= zeta) {
		    if (f0 >= 0.) {
			--nomtch;
		    }
		    ++j;
		    if (j > *n) {
			goto L10;
		    }
		    f0 = a[j] - trial - b[j - 1] / f0;
		}
		j += 2;
		--nomtch;
	    }
L10:
	    if (nomtch < i__) {
		root = trial;
	    } else {
		ev[i__] = trial;
		nom = min(*nv,nomtch);
		ev[nom] = trial;
	    }
	    trial = (root + ev[i__]) * .5;
	}
    }

/*     find the eigenvectors via a backtransformation step */

    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	root = ev[i__];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    y[j] = 1.;
	}
	++ia;
	if (i__ == 1 || (d__1 = ev[i__ - 1] - root, abs(d__1)) >= eps) {
	    ia = 0;
	}
	elim1 = a[1] - root;
	elim2 = w[1];
	i__2 = nn;
	for (j = 1; j <= i__2; ++j) {
	    if (abs(elim1) <= (d__1 = w[j], abs(d__1))) {
		ta[j] = w[j];
		tb[j] = a[j + 1] - root;
		p[j] = w[j + 1];
		temp = 1.;
		if ((d__1 = w[j], abs(d__1)) > zeta) {
		    temp = elim1 / w[j];
		}
		elim1 = elim2 - temp * tb[j];
		elim2 = -temp * w[j + 1];
	    } else {
		ta[j] = elim1;
		tb[j] = elim2;
		p[j] = 0.;
		temp = w[j] / elim1;
		elim1 = a[j + 1] - root - temp * elim2;
		elim2 = w[j + 1];
	    }
	    b[j] = temp;
	}
	ta[*n] = elim1;
	tb[*n] = 0.;
	p[*n] = 0.;
	p[nn] = 0.;
	iter = 1;
	if (ia != 0) {
	    goto L40;
	}
L20:
	m = *n + 1;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    --m;
	    done = FALSE_;
	    while(! done) {
		done = TRUE_;
		if (*n - m - 1 < 0) {
		    elim1 = y[m];
		} else {
		    if (*n - m - 1 == 0) {
			elim1 = y[m] - y[m + 1] * tb[m];
		    } else {
			elim1 = y[m] - y[m + 1] * tb[m] - y[m + 2] * p[m];
		    }
		}
		if (abs(elim1) <= sigma) {
		    temp = ta[m];
		    if (abs(temp) < delta) {
			temp = delta;
		    }
		    y[m] = elim1 / temp;
		} else {
		    i__3 = *n;
		    for (k = 1; k <= i__3; ++k) {
			y[k] /= sigma;
		    }
		    done = FALSE_;
		}
	    }
	}
	if (iter == 2) {
	    goto L50;
	}
	++iter;
L30:
	elim1 = y[1];
	i__2 = nn;
	for (j = 1; j <= i__2; ++j) {
	    if (ta[j] != w[j]) {
		y[j] = elim1;
		elim1 = y[j + 1] - elim1 * b[j];
	    } else {
		y[j] = y[j + 1];
		elim1 -= y[j + 1] * b[j];
	    }
	}
	y[*n] = elim1;
	goto L20;
L40:
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    d__1 = rvalue * rand1;
	    rand1 = d_mod(&d__1, &rpower);
	    y[j] = rand1 / rpow1 - 1.;
	}
	goto L20;
L50:
	if (ia != 0) {
	    i__2 = ia;
	    for (j1 = 1; j1 <= i__2; ++j1) {
		k = i__ - j1;
		temp = 0.;
		i__3 = *n;
		for (j = 1; j <= i__3; ++j) {
		    temp += y[j] * vec_ref(j, k);
		}
		i__3 = *n;
		for (j = 1; j <= i__3; ++j) {
		    y[j] -= temp * vec_ref(j, k);
		}
	    }
	}
	if (iter == 1) {
	    goto L30;
	}
	elim1 = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* Computing MAX */
	    d__2 = elim1, d__3 = (d__1 = y[j], abs(d__1));
	    elim1 = max(d__2,d__3);
	}
	temp = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    elim2 = y[j] / elim1;
	    temp += elim2 * elim2;
	}
	temp = 1. / (sqrt(temp) * elim1);
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    y[j] *= temp;
	    if ((d__1 = y[j], abs(d__1)) < rho) {
		y[j] = 0.;
	    }
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    vec_ref(j, i__) = y[j];
	}
    }

/*     normalization of the eigenvalues and eigenvectors */

    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    y[j] = vec_ref(j, i__);
	}
	mk = *n * (*n - 1) / 2 - 3;
	mk1 = 3;
	i__2 = *n - 2;
	for (j = 1; j <= i__2; ++j) {
	    t = 0.;
	    m = *n - j;
	    i__3 = *n;
	    for (k = m; k <= i__3; ++k) {
		t += dd[mk + k] * y[k];
	    }
	    i__3 = *n;
	    for (k = m; k <= i__3; ++k) {
		epr = t * dd[mk + k];
		y[k] -= epr * 2.;
	    }
	    mk -= mk1;
	    ++mk1;
	}
	t = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    t += y[j] * y[j];
	}
	xnorm = sqrt(t);
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    y[j] /= xnorm;
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    vec_ref(j, i__) = y[j];
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ev[i__] *= anorm;
    }
    return 0;
} /* diagq_ */

#undef vec_ref


