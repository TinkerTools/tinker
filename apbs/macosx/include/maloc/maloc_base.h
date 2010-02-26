/**
 *  @file       maloc_base.h
 *  @brief      The base (or foundation) header for MALOC.
 *  @note       This header sets things up correctly for using ISO/ANSI-C.          
 *              The following macros affect the behavior of the header:     
    @verbatim
     Inlining for speed:  (Normal C functions if VINLINE_XXX not defined.)      
     ------------------
     -DVINLINE_VNM : Enables macro replacement of time-critical funcs in VNM.                                                
    @endverbatim                                                                                                      
 *                                                                                              
 *  @version    $Id: maloc_base.h,v 1.28 2008/02/05 00:08:31 fetk Exp $
 *  @author     Michael Holst
 *  
 *  @attention
 *  @verbatim
 *
 * MALOC = < Minimal Abstraction Layer for Object-oriented C >
 * Copyright (C) 1994--2008 Michael Holst
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 *  @endverbatim
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

/**
 *  @brief  Fix to broken <time.h> on old SunOS 
 */
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
#define __FAVOR_BSD  /**< @brief Linux: uses sigsetjmp as the setjmp function */
#define _BSD_SIGNALS /**< @brief IRIX:  uses sigsetjmp as the setjmp function */

/* Now finally include the 15th header, setjmp.h */
#include <setjmp.h>

#if defined(__cplusplus)
/** @brief Setup so this include file (and subsequent) will work for both C and C++ */
#   define VCXX
/** @brief Setup so this include file (and subsequent) will work for both C and C++ */
#   define VEXTERNC extern "C"
#else
/** @brief Setup so this include file (and subsequent) will work for both C and C++ */
#   define VCC
/** @brief Setup so this include file (and subsequent) will work for both C and C++ */
#   define VEXTERNC extern
#endif

/*
 * ***************************************************************************
 * Private and Public type modifier simulation
 * ***************************************************************************
 */

#define VPRIVATE static      /**< @brief Mimic C++ "Private" type modifier */
#define VPUBLIC  /*empty*/   /**< @brief Mimic C++ "Public" type modifier  */


/** @brief Slick assertion macro */
#define VWARN1(file, lineno) (fprintf(stderr,"VWARN: ASSERTION FAILURE! filename %s, line %u\n", (file), (lineno)), 0)
/** @brief Slick assertion macro */
#define VASSERT1(file, lineno) (fprintf(stderr,"VASSERT: ASSERTION FAILURE! filename %s, line %u\n", (file), (lineno)), exit(1), 0)
/** @brief Slick assertion macro */
#define VASSERT2(file, lineno) (fprintf(stderr,"VASSERT: ASSERTION FAILURE! filename %s, line %u\n", (file), (lineno)), abort(), 0)
/** @brief Slick assertion macro */
#define VASSERT3(file, lineno, ex) (fprintf(stderr,"VASSERT: ASSERTION FAILURE!  filename %s, line %u, (%s)\n", (file), (lineno), (#ex)), abort(), 0)

/** @brief Slick assertion macro */
#define VWARN(ex)   ((void) ((ex) ? 0 : VWARN1(__FILE__, __LINE__)))
/** @brief Slick assertion macro */
#define VASSERT(ex) ((void) ((ex) ? 0 : VASSERT3(__FILE__, __LINE__, ex)))


/** @brief A userful error handling macro */
#define VJMPERR0(x) if (!(x)) goto VERROR0
/** @brief A userful error handling macro */
#define VJMPERR1(x) if (!(x)) goto VERROR1
/** @brief A userful error handling macro */
#define VJMPERR2(x) if (!(x)) goto VERROR2
/** @brief A userful error handling macro */
#define VJMPERR3(x) if (!(x)) goto VERROR3
/** @brief A userful error handling macro */
#define VJMPERR4(x) if (!(x)) goto VERROR4
/** @brief A userful error handling macro */
#define VJMPERR5(x) if (!(x)) goto VERROR5
/** @brief A userful error handling macro */
#define VJMPERR6(x) if (!(x)) goto VERROR6
/** @brief A userful error handling macro */
#define VJMPERR7(x) if (!(x)) goto VERROR7
/** @brief A userful error handling macro */
#define VJMPERR8(x) if (!(x)) goto VERROR8
/** @brief A userful error handling macro */
#define VJMPERR9(x) if (!(x)) goto VERROR9

/*
 * ***************************************************************************
 * Global constants
 * ***************************************************************************
 */

/** @brief Global constant */
#define VPI                3.14159265358979323846
/** @brief Global constant. 1e9 just fits into 32-bit signed int */
#define VLARGE             1.0e+9   
/** @brief Global constant */
#define VSMALL             1.0e-9
/** @brief Global constant */
#define VVLARGE            1.0e+15
/** @brief Global constant */
#define VVSMALL            1.0e-15
/** @brief Global constant */
#define VPRTKEY            10000

/** @brief Global constant */
#define VPTRSIZE           4
/** @brief Global constant */
#define VMAX_ARGNUM        50
/** @brief Global constant */
#define VMAX_ARGLEN        1024
/** @brief Global constant */
#define VMAX_BUFSIZE       8192

/* 
 * #define  VMAX_OBJECTS      16777216     //(1<<24) = 2^24
 * #define  VBLOCK_POWER      12
 */

/** @brief Global constant */
#define VMAX_OBJECTS       1073741824    /* (1<<31) = 2^31 */
/** @brief Global constant */
#define VBLOCK_POWER       14

/** @brief Global constant */
#define VNULL              NULL
/** @brief Global constant */
#define VINULL             -1
/** @brief Global constant */
#define VTRUE              1
/** @brief Global constant */
#define VFALSE             0
/** @brief Global constant */
#define VSTDMODE           0600

/** @brief Global constant */
#define VNULL_STRING       "\0"
/** @brief Global constant */
#define VBLANK_STRING      " "
/** @brief Global constant */
#define VNEWLINE_STRING    "\n"

/** @brief Global constant */
#define VNULL_SYMBOL       '\0'
/** @brief Global constant */
#define VBLANK_SYMBOL      ' '
/** @brief Global constant */
#define VNEWLINE_SYMBOL    '\n'
/** @brief Global constant */
#define VRDIN_SYMBOL       '<'
/** @brief Global constant */
#define VRDOUT_SYMBOL      '>'
/** @brief Global constant */
#define VPIPE_SYMBOL       '|'
/** @brief Global constant */
#define VDELIM_SET         " ><|&"


/** @brief Mathematical macro */
#define VABS(x)            ((x) >= 0 ? (x) : -(x))
/** @brief Mathematical macro */
#define VMIN2(x,y)         ((x) <= (y) ? (x) : (y))
/** @brief Mathematical macro */
#define VMAX2(x,y)         ((x) >= (y) ? (x) : (y))
/** @brief Mathematical macro */
#define VSIGN(x,y)         ((y) >= 0 ? (VABS(x)) : (-VABS(x)))

/** @brief Mathematical macro */
#define VODD(x)            ((x)&1)
/** @brief Mathematical macro */
#define VEVEN(x)           (!((x)&1))
/** @brief Mathematical macro */
#define VZERO(x)           ((x)==0)
/** @brief Mathematical macro */
#define VPOS(x)            ((x)>0)
/** @brief Mathematical macro */
#define VNEG(x)            ((x)<0)
/** @brief Mathematical macro */
#define VEVENP(x)          (VEVEN(x) && VPOS(x))
/** @brief Mathematical macro */
#define VEVENN(x)          (VEVEN(x) && VNEG(x))

/** @brief Mathematical macro */
#define VSQRT(x)           (sqrt(x))
/** @brief Mathematical macro */
#define VSQR(x)            ((x)*(x))
/** @brief Mathematical macro */
#define VSIN(x)            (sin(x))
/** @brief Mathematical macro */
#define VCOS(x)            (cos(x))
/** @brief Mathematical macro */
#define VTAN(x)            (tan(x))
/** @brief Mathematical macro */
#define VASIN(x)           (asin(x))
/** @brief Mathematical macro */
#define VACOS(x)           (acos(x))
/** @brief Mathematical macro */
#define VATAN(x)           (atan(x))
/** @brief Mathematical macro */
#define VSINH(x)           (sinh(x))
/** @brief Mathematical macro */
#define VCOSH(x)           (cosh(x))
/** @brief Mathematical macro */
#define VTANH(x)           (tanh(x))
/** @brief Mathematical macro */
#define VEXP(x)            (exp(x))
/** @brief Mathematical macro */
#define VLOG(x)            (log(x))
/** @brief Mathematical macro */
#define VPOW(x,y)          (pow(x,y))
/** @brief Mathematical macro */
#define VRINT(x)           ((int)(floor((x)+0.5)))

/** @brief Mathematical macro */
#define VRAND              (rand())
/** @brief Mathematical macro */
#define VRANDMAX           (RAND_MAX)

/**
 * @brief Inlining via macros for speed
 */
#if 1
#   define VINLINE_MALOC
#endif

#endif /* _MALOC_BASE_H_ */

