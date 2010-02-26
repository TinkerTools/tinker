/** @defgroup Vcap Vcap class
 *  @brief  Collection of routines which cap certain exponential and hyperbolic
 *          functions
 *  @note   These routines are based on FORTRAN code by Mike Holst
 */

/**
 *  @file     vcap.h
 *  @ingroup  Vcap
 *  @brief    Contains declarations for class Vcap
 *  @version  $Id: vcap.h 1033 2006-12-29 17:08:22Z sobolevnrm $
 *  @author   Nathan A. Baker
 *
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Center for Computational Biology
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2002-2007.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2002.  The Regents of the University of
 * California.  
 * Portions Copyright (c) 1995.  Michael Holst.
 *
 * This file is part of APBS.
 *
 * APBS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * APBS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with APBS; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 * Linking APBS statically or dynamically with other modules is making a
 * combined work based on APBS. Thus, the terms and conditions of the GNU
 * General Public License cover the whole combination.
 * 
 * SPECIAL GPL EXCEPTION
 * In addition, as a special exception, the copyright holders of APBS
 * give you permission to combine the APBS program with free software
 * programs and libraries that are released under the GNU LGPL or with
 * code included in releases of ISIM, Ion Simulator Interface, PMV, PyMOL
 * SMOL, VMD, and Vision. Such combined software may be linked with APBS and 
 * redistributed together in original or modified form as mere aggregation
 * without requirement that the entire work be under the scope of the GNU 
 * General Public License. This special exception permission is also extended
 * to any software listed in the SPECIAL GPL EXCEPTION clauses by the PMG,
 * FEtk, MC, or MALOC libraries.
 * 
 * Note that people who make modified versions of APBS are not obligated
 * to grant this special exception for their modified versions; it is
 * their choice whether to do so. The GNU General Public License gives
 * permission to release a modified version without this exception; this
 * exception also makes it possible to release a modified version which
 * carries forward this exception.
 *
 * @endverbatim
 */

#ifndef _VCAP_H_
#define _VCAP_H_

/** @brief   Maximum argument for exp(), sinh(), or cosh()
 *  @ingroup Vcap
 */
#define EXPMAX  85.00

/** @brief   Minimum argument for exp(), sinh(), or cosh()
 *  @ingroup Vcap
 */
#define EXPMIN -85.00

#include "maloc/maloc.h"

/** @brief   Provide a capped exp() function
 *
 *           If the argument x of Vcap_exp() exceeds EXPMAX or EXPMIN, then we
 *           return exp(EXPMAX) or exp(EXPMIN) rather than exp(x).
 *
 *  @note    Original FORTRAN routine from PMG library by Mike Holst
 *           Original notes:
 *           to control overflow in the hyperbolic and exp functions, note
 *           that the following are the argument limits of the various 
 *           functions on various machines after which overflow occurs:
 *           Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
 *           sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
 *           dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
 *
 *  @ingroup Vcap
 *  @author  Nathan Baker (based on FORTRAN code by Mike Holst)
 *  @return  exp(x) or capped equivalent
 */
VEXTERNC double Vcap_exp(
        double x, /**< Argument to exp() */
        int *ichop /**< Set to 1 if function capped, 0 otherwise */
        );


/** @brief   Provide a capped sinh() function
 *
 *           If the argument x of Vcap_sinh() exceeds EXPMAX or EXPMIN, then we
 *           return sinh(EXPMAX) or sinh(EXPMIN) rather than sinh(x).
 *
 *  @note    Original FORTRAN routine from PMG library by Mike Holst
 *           Original notes:
 *           to control overflow in the hyperbolic and exp functions, note
 *           that the following are the argument limits of the various 
 *           functions on various machines after which overflow occurs:
 *           Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
 *           sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
 *           dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
 *
 *  @ingroup Vcap
 *  @author  Nathan Baker (based on FORTRAN code by Mike Holst)
 *  @return  sinh(x) or capped equivalent
 */
VEXTERNC double Vcap_sinh(
        double x, /**< Argument to sinh() */
        int *ichop /**< Set to 1 if function capped, 0 otherwise */
        );

/** @brief   Provide a capped cosh() function
 *
 *           If the argument x of Vcap_cosh() exceeds EXPMAX or EXPMIN, then we
 *           return cosh(EXPMAX) or cosh(EXPMIN) rather than cosh(x).
 *
 *  @note    Original FORTRAN routine from PMG library by Mike Holst
 *           Original notes:
 *           to control overflow in the hyperbolic and exp functions, note
 *           that the following are the argument limits of the various
 *           functions on various machines after which overflow occurs:
 *           Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
 *           sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
 *           dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
 *
 *  @ingroup Vcap
 *  @author  Nathan Baker (based on FORTRAN code by Mike Holst)
 *  @return  cosh(x) or capped equivalent
 */
VEXTERNC double Vcap_cosh(
        double x, /**< Argument to cosh() */
        int *ichop /**< Set to 1 if function capped, 0 otherwise */
        );

#endif    /* ifndef _VCAP_H_ */
