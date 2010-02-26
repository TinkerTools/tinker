/** @defgroup FEMparm FEMparm class
 *  @brief    Parameter structure for FEM-specific variables from input files
 */

/**
 *  @file     femparm.h
 *  @ingroup  FEMparm
 *  @brief    Contains declarations for class FEMparm
 *  @version  $Id: femparm.h 1033 2006-12-29 17:08:22Z sobolevnrm $
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


#ifndef _FEMPARM_H_
#define _FEMPARM_H_

/* Generic header files */
#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/vstring.h"

/**
 * @brief  Adaptive refinment error estimate tolerance key
 * @ingroup FEMparm
 * @author  Nathan Baker
 */
enum eFEMparm_EtolType {
    FET_SIMP=0,  /**< per-simplex error tolerance */
    FET_GLOB=1,  /**< global error tolerance */
    FET_FRAC=2   /**< fraction of simplices we want to have refined */
};

/**
 * @brief  Declare FEparm_EtolType type
 * @ingroup  FEMparm
 * @author  Nathan Baker
 */
typedef enum eFEMparm_EtolType FEMparm_EtolType;

/**
 * @brief  Adaptive refinment error estimator method
 * @ingroup FEMparm
 * @note  Do not change these values; they correspond to settings in FEtk
 * @author  Nathan Baker
 */
enum eFEMparm_EstType {
    FRT_UNIF=0,  /**< Uniform refinement */
    FRT_GEOM=1,  /**< Geometry-based (i.e. surfaces and charges) refinement */
    FRT_RESI=2,  /**< Nonlinear residual estimate-based refinement */
    FRT_DUAL=3,  /**< Dual-solution weight nonlinear residual estimate-based
                  * refinement */ 
    FRT_LOCA=4  /**< Local problem error estimate-based refinement */
};

/**
 * @brief  Declare FEMparm_EstType type
 * @ingroup  FEMparm
 */
typedef enum eFEMparm_EstType FEMparm_EstType;

/**
 * @brief  Calculation type 
 * @ingroup FEMparm
 */
enum eFEMparm_CalcType {
    FCT_MANUAL,  /**< fe-manual */
	FCT_NONE  /**< unspecified */
};

/**
 * @brief  Declare FEMparm_CalcType type
 * @ingroup  FEMparm
 */
typedef enum eFEMparm_CalcType FEMparm_CalcType;

/**
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @brief   Parameter structure for FEM-specific variables from input files
 */
struct sFEMparm {
   
    int parsed;  /**< Flag:  Has this structure been filled with
                  * anything other than * the default values? (0 = no,
                  * 1 = yes) */
    FEMparm_CalcType type;  /**<  Calculation type */
    int settype;  /**< Boolean */
    double glen[3];  /**< Domain side lengths (in &Aring;) */
    int setglen;  /**< Boolean */
    double etol;  /**< Error tolerance for refinement; interpretation depends 
                   * on the adaptive refinement method chosen */
    int setetol;  /**< Boolean */
    FEMparm_EtolType ekey;  /**< Adaptive refinment interpretation of error 
                            * tolerance */
    int setekey;  /**< Boolean */
    FEMparm_EstType akeyPRE;  /**< Adaptive refinment error estimator method 
                               * for pre-solution refine.  Note, this should
                               * either be FRT_UNIF or FRT_GEOM.  */
    int setakeyPRE;  /**< Boolean */
    FEMparm_EstType akeySOLVE;  /**< Adaptive refinment error estimator method 
                               * for a posteriori solution-based refinement. */
    int setakeySOLVE;  /**< Boolean */
    int targetNum;    /**< Initial mesh will continue to be marked and refined
                        * with the method specified by akeyPRE until the mesh
                        * contains this many vertices or until targetRes is
                        * reached. */
    int settargetNum;  /**< Boolean */
    double targetRes; /**< Initial mesh will continue to be marked and refined
                        * with the method specified by akeyPRE until the mesh
                        * contains no markable simplices with longest edges
                        * above this size or until targetNum is reached. */
    int settargetRes;  /**< Boolean */
    int maxsolve;  /**< Maximum number of solve-estimate-refine cycles */
    int setmaxsolve;  /**< Boolean */
    int maxvert;  /**< Maximum number of vertices in mesh (ignored if less 
                   * than zero) */
    int setmaxvert;  /**< Boolean */
	int pkey;		/**< Boolean sets the pkey type for going into AM_Refine
					  * pkey = 0 for non-HB based methods
					  * pkey = 1 for HB based methods */

};

/** @typedef FEMparm
 *  @ingroup FEMparm
 *  @brief   Declaration of the FEMparm class as the FEMparm structure
 */
typedef struct sFEMparm FEMparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (nosh.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct FEMparm
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param  type  FEM calculation type
 *  @returns Newly allocated and initialized Vpmgp object
 */
VEXTERNC FEMparm* FEMparm_ctor(FEMparm_CalcType type);

/** @brief   FORTRAN stub to construct FEMparm
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param   thee Pointer to allocated FEMparm object
 *  @param  type  FEM calculation type
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int FEMparm_ctor2(FEMparm *thee, FEMparm_CalcType type);

/** @brief   Object destructor
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of FEMparm object
 */
VEXTERNC void FEMparm_dtor(FEMparm **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to FEMparm object
 */
VEXTERNC void FEMparm_dtor2(FEMparm *thee);

/** 
 * @brief   Consistency check for parameter values stored in object
 * @ingroup FEMparm
 * @author  Nathan Baker
 * @param   thee   FEMparm object
 * @returns 1 if OK, 0 otherwise
 */
VEXTERNC int FEMparm_check(FEMparm *thee);

/**	@brief	Copy target object into thee
	@ingroup	FEMparm
	@author	Nathan Baker
	@param	thee	Destination object
	@param	source	Source object
*/
VEXTERNC void FEMparm_copy(FEMparm *thee, FEMparm *source);

/** 
 * @brief   Parse an MG keyword from an input file
 * @ingroup MGparm
 * @author  Nathan Baker
 * @param   thee   MGparm object 
 * @param   tok    Token to parse
 * @param   sock   Stream for more tokens
 * @return   1 if matched and assigned; -1 if matched, but there's
 * some sort of error (i.e., too few args); 0 if not matched
 */
VEXTERNC int FEMparm_parseToken(FEMparm *thee, char tok[VMAX_BUFSIZE], 
  Vio *sock);

#endif 

