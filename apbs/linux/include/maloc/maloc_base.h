/*
 * ***************************************************************************
 * MALOC = < Minimal Abstraction Layer for Object-oriented C >
 * Copyright (C) 1994--2000  Michael Holst
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 675 Mass Ave, Cambridge, MA 02139, USA.
 * 
 * rcsid="$Id: maloc_base.h,v 1.23 2003/07/06 01:54:09 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     maloc_base.h
 *
 * Purpose:  The base (or foundation) header for MALOC.
 *
 * Notes:    This header sets things up correctly for using ISO/ANSI-C.
 *           The following macros affect the behavior of the header:
 *
 *    Inlining for speed:  (Normal C functions if VINLINE_XXX not defined.)
 *    -------------------
 *    -DVINLINE_VNM : Enables macro replacement of time-critical funcs in VNM.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _MALOC_BASE_H_
#define _MALOC_BASE_H_

/*
 * ***************************************************************************
 * Proper ISO-C header setup (with a slight "signals" tweek for setjmp)
 * ***************************************************************************
 */

/* Get the fifteen ISO-C headers (CAREFUL: POSIX/BSD flags delicate...) */

/* Some compilers don't set this for you; GCC does with -ansi */
/*
 * if !defined(__STDC__)
 *    define __STDC__ 1
 * endif
 */

/* Old Sparc compilers need this to give you prototypes */
/*
 * if !defined(__USE_FIXED_PROTOTYPES__)
 *     define __USE_FIXED_PROTOTYPES__
 * endif
 */

/* Include 14 of the 15 ISO-C headers (postponing setjmp.h) */
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <signal.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Fix to broken <time.h> on old SunOS */
#if !defined(CLOCKS_PER_SEC)
#   define CLOCKS_PER_SEC    60
#endif

/*
 * Problems using setjmp/longjmp for use in the MC-shell.
 *
 * Problem:  Some implementations of ISO-C "setjmp/longjmp"  do not return
 *           the interrupt mask to its pre-jump state after returning.
 *           The behavior this produces is for example you can capture a
 *           single CTRL-C, but that is it; the mask for CTRL-C is wiped
 *           after the first interrupt is handled.
 *
 * Solution: Use the "sigsetjmp/siglongjmp" extensions provided by most
 *           UNIX variants.  You just have to set an appropriate macro
 *           before including <setjmp.h> to get sigsetjmp rather than
 *           setjmp behavior.
 *
 * Notes:    You can run into trouble (e.g. some versions of Linux) if
 *           you set some of these special signal macros before some of
 *           the other ISO-C headers.  Therefore, we only set the macros
 *           just before including <setjmp.h> as the final ISO-C header.
 */
#define __FAVOR_BSD  /* Linux: uses sigsetjmp as the setjmp function */
#define _BSD_SIGNALS /* IRIX:  uses sigsetjmp as the setjmp function */

/* Now finally include the 15th header, setjmp.h */
#include <setjmp.h>

/*
 * ***************************************************************************
 * Setup so this include file (and subsequent) will work for both C and C++
 * ***************************************************************************
 */

#if defined(__cplusplus)
#   define VCXX
#   define VEXTERNC extern "C"
#else
#   define VCC
#   define VEXTERNC extern
#endif

/*
 * ***************************************************************************
 * Private and Public type modifier simulation
 * ***************************************************************************
 */

#define VPRIVATE static      /* Mimic C++ "Private" type modifier */
#define VPUBLIC  /*empty*/   /* Mimic C++ "Public" type modifier  */

/*
 * ***************************************************************************
 * Slick assertion macros
 * ***************************************************************************
 */

#define VWARN1(file, lineno) (fprintf(stderr,"VWARN: ASSERTION FAILURE! filename %s, line %u\n", (file), (lineno)), 0)
#define VASSERT1(file, lineno) (fprintf(stderr,"VASSERT: ASSERTION FAILURE! filename %s, line %u\n", (file), (lineno)), exit(1), 0)
#define VASSERT2(file, lineno) (fprintf(stderr,"VASSERT: ASSERTION FAILURE! filename %s, line %u\n", (file), (lineno)), abort(), 0)
#define VASSERT3(file, lineno, ex) (fprintf(stderr,"VASSERT: ASSERTION FAILURE!  filename %s, line %u, (%s)\n", (file), (lineno), (#ex)), abort(), 0)

#define VWARN(ex)   ((void) ((ex) ? 0 : VWARN1(__FILE__, __LINE__)))
#define VASSERT(ex) ((void) ((ex) ? 0 : VASSERT3(__FILE__, __LINE__, ex)))

/*
 * ***************************************************************************
 * A useful error handling macro
 * ***************************************************************************
 */

#define VJMPERR0(x) if (!(x)) goto VERROR0
#define VJMPERR1(x) if (!(x)) goto VERROR1
#define VJMPERR2(x) if (!(x)) goto VERROR2
#define VJMPERR3(x) if (!(x)) goto VERROR3
#define VJMPERR4(x) if (!(x)) goto VERROR4
#define VJMPERR5(x) if (!(x)) goto VERROR5
#define VJMPERR6(x) if (!(x)) goto VERROR6
#define VJMPERR7(x) if (!(x)) goto VERROR7
#define VJMPERR8(x) if (!(x)) goto VERROR8
#define VJMPERR9(x) if (!(x)) goto VERROR9

/*
 * ***************************************************************************
 * Global constants
 * ***************************************************************************
 */

#define VPI                3.14159265358979323846
#define VLARGE             1.0e+9   /* 1e9 just fits into 32-bit signed int */
#define VSMALL             1.0e-9
#define VPRTKEY            10000

#define VPTRSIZE           4
#define VMAX_ARGNUM        50
#define VMAX_ARGLEN        1024
#define VMAX_BUFSIZE       8192
#define VMAX_OBJECTS       16777216    /* (1<<24) = 2^24 */

#define VNULL              NULL
#define VINULL             -1
#define VTRUE              1
#define VFALSE             0
#define VSTDMODE           0600

#define VNULL_STRING       "\0"
#define VBLANK_STRING      " "
#define VNEWLINE_STRING    "\n"

#define VNULL_SYMBOL       '\0'
#define VBLANK_SYMBOL      ' '
#define VNEWLINE_SYMBOL    '\n'
#define VRDIN_SYMBOL       '<'
#define VRDOUT_SYMBOL      '>'
#define VPIPE_SYMBOL       '|'
#define VDELIM_SET         " ><|&"

/*
 * ***************************************************************************
 * Mathematical macros 
 * ***************************************************************************
 */

#define VABS(x)            ((x) >= 0 ? (x) : -(x))
#define VMIN2(x,y)         ((x) <= (y) ? (x) : (y))
#define VMAX2(x,y)         ((x) >= (y) ? (x) : (y))
#define VSIGN(x,y)         ((y) >= 0 ? (VABS(x)) : (-VABS(x)))

#define VODD(x)            ((x)&1)
#define VEVEN(x)           (!((x)&1))
#define VZERO(x)           ((x)==0)
#define VPOS(x)            ((x)>0)
#define VNEG(x)            ((x)<0)
#define VEVENP(x)          (VEVEN(x) && VPOS(x))
#define VEVENN(x)          (VEVEN(x) && VNEG(x))

#define VSQRT(x)           (sqrt(x))
#define VSQR(x)            ((x)*(x))
#define VSIN(x)            (sin(x))
#define VCOS(x)            (cos(x))
#define VTAN(x)            (tan(x))
#define VASIN(x)           (asin(x))
#define VACOS(x)           (acos(x))
#define VATAN(x)           (atan(x))
#define VSINH(x)           (sinh(x))
#define VCOSH(x)           (cosh(x))
#define VTANH(x)           (tanh(x))
#define VEXP(x)            (exp(x))
#define VPOW(x,y)          (pow(x,y))
#define VRINT(x)           ((int)(floor((x)+0.5)))

#define VRAND              (rand())
#define VRANDMAX           (RAND_MAX)

/*
 * ***************************************************************************
 * Inlining via macros for speed
 * ***************************************************************************
 */

#if 1
#   define VINLINE_MALOC
#endif

#endif /* _MALOC_BASE_H_ */

