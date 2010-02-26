/**
 *  @file    vunit.h
 *  @ingroup Vunit
 *  @author  Nathan Baker
 *  @brief   Contains a collection of useful constants and conversion factors
 *  @author  Nathan A. Baker
 *  @version $Id: vunit.h 1033 2006-12-29 17:08:22Z sobolevnrm $
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

/** @defgroup Vunit Vunit class
 *  @brief    Collection of constants and conversion factors
 */

#ifndef _VUNIT_H_
#define _VUNIT_H_

/** @brief   Multiply by this to convert J to cal 
 *  @ingroup Vunit */
#define Vunit_J_to_cal		4.1840000e+00

/** @brief   Multiply by this to convert cal to J
 *  @ingroup Vunit */
#define Vunit_cal_to_J		2.3900574e-01

/** @brief   Multiply by this to convert amu to kg
 *  @ingroup Vunit */
#define Vunit_amu_to_kg 	1.6605402e-27

/** @brief   Multiply by this to convert kg to amu
 *  @ingroup Vunit */
#define Vunit_kg_to_amu 	6.0221367e+26

/** @brief   Multiply by this to convert ec to C
 *  @ingroup Vunit */
#define Vunit_ec_to_C		1.6021773e-19

/** @brief   Multiply by this to convert C to ec
 *  @ingroup Vunit */
#define Vunit_C_to_ec		6.2415065e+18

/** @brief   Charge of an electron in C
 *  @ingroup Vunit */
#define Vunit_ec		1.6021773e-19

/** @brief   Boltzmann constant
 *  @ingroup Vunit */
#define Vunit_kb		1.3806581e-23

/** @brief   Avogadro's number
 *  @ingroup Vunit */
#define Vunit_Na		6.0221367e+23

/** @brief   Pi
 *  @ingroup Vunit */
#define Vunit_pi		VPI

/** @brief   Vacuum permittivity
 *  @ingroup Vunit */
#define Vunit_eps0		8.8541878e-12

/** @brief \f${e_c}^2/\AA\f$ in ESU units => kcal/mol 
 *  @ingroup Vunit */
#define Vunit_esu_ec2A		3.3206364e+02

/** @brief \f$k_b\f$ in ESU units => kcal/mol 
 *  @ingroup Vunit */
#define Vunit_esu_kb            1.9871913e-03

#endif /* ifndef _VUNIT_H_ */
