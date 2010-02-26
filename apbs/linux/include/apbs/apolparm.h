/** @defgroup APOLparm APOLparm class
 *  @brief    Parameter structure for APOL-specific variables from input files
 */

/**
 *  @file     femparm.h
 *  @ingroup  APOLparm
 *  @brief    Contains declarations for class APOLparm
 *  @version  $Id: femparm.h 907 2006-07-27 20:36:20Z sobolevnrm $
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


#ifndef _APOLPARM_H_
#define _APOLPARM_H_

/* Generic header files */
#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/vstring.h"
#include "apbs/vparam.h"

/**
* @ingroup APOLparm
 * @brief  Define energy calculation enumeration
 */
enum eAPOLparm_calcEnergy {
    ACE_NO=0, /**< Do not perform energy calculation */
    ACE_TOTAL=1, /**< Calculate total energy only */
    ACE_COMPS=2 /**< Calculate per-atom energy components */
};

/**
* @ingroup APOLparm
 * @brief  Define eAPOLparm_calcEnergy enumeration as APOLparm_calcEnergy
 */
typedef enum eAPOLparm_calcEnergy APOLparm_calcEnergy;

/**
* @ingroup APOLparm
 * @brief  Define force calculation enumeration
 */
enum eAPOLparm_calcForce {
    ACF_NO=0, /**< Do not perform force calculation */
    ACF_TOTAL=1, /**< Calculate total force only */
    ACF_COMPS=2 /**< Calculate per-atom force components */
};

/**
* @ingroup APOLparm
 * @brief  Define eAPOLparm_calcForce enumeration as APOLparm_calcForce
 */
typedef enum eAPOLparm_calcForce APOLparm_calcForce;

/**
* @ingroup APOLparm
 * @brief  Define force calculation enumeration
 */
enum eAPOLparm_doCalc {
    ACD_NO=0, /**< Do not perform calculation */
    ACD_YES=1, /**< Perform calculations */
	ACD_ERROR=2 /**< Error setting up calculation */
};

/**
* @ingroup APOLparm
 * @brief  Define eAPOLparm_calcForce enumeration as APOLparm_calcForce
 */
typedef enum eAPOLparm_doCalc APOLparm_doCalc;


/**
 *  @ingroup APOLparm
 *  @author  David Gohara
 *  @brief   Parameter structure for APOL-specific variables from input files
 */
struct sAPOLparm {
   
    int parsed;  /**< Flag:  Has this structure been filled with
                  * anything other than * the default values? (0 = no,
                  * 1 = yes) */
    int dime[3];  /**< Grid dimensions */
	int setdime;  /**< Flag, @see dime */
	
	double glen[3];  /**< Grid side lengths. */
	int setglen;  /**< Flag, @see fglen */
	
	int molid;  /**< Molecule ID to perform calculation on */
    int setmolid;  /**< Flag, @see molid */
	
	double bconc; /**< Vacc sphere density */
    int setbconc; /**< Flag, @see bconc */
	
	double sdens; /**< Vacc sphere density */
    int setsdens; /**< Flag, @see sdens */
	
	double dpos; /**< Atom position offset */
    int setdpos; /**< Flag, @see dpos */
	
	double press; /**< Solvent pressure */
	int setpress; /**< Flag, @see press */
	
	Vsurf_Meth srfm;  /**< Surface calculation method */
    int setsrfm;  /**< Flag, @see srfm */

    double srad;  /**< Solvent radius */
    int setsrad;  /**< Flag, @see srad */
    
	double swin;  /**< Cubic spline window */
    int setswin;  /**< Flag, @see swin */
	
	double temp;  /**< Temperature (in K) */
    int settemp;  /**< Flag, @see temp */
	
    double gamma;  /**< Surface tension for apolar energies/forces
					* (in kJ/mol/A^2) */
    int setgamma;  /**< Flag, @see gamma */
	
	APOLparm_calcEnergy calcenergy;  /**< Energy calculation flag */
    int setcalcenergy;  /**< Flag, @see calcenergy */
	
    APOLparm_calcForce calcforce;  /**< Atomic forces calculation */
    int setcalcforce;  /**< Flag, @see calcforce */
	
	double watsigma;  /**< Water oxygen Lennard-Jones radius (A) */
	double watepsilon;  /**< Water oxygen Lennard-Jones well depth (kJ/mol) */
	double sasa; /**< Solvent accessible surface area for this calculation */
	double sav;   /**< Solvent accessible volume for this calculation */
	double wcaEnergy; /** wcaEnergy */	
	double totForce[3]; /**< Total forces on x, y, z */
};

/** @typedef APOLparm
 *  @ingroup APOLparm
 *  @brief   Declaration of the APOLparm class as the APOLparm structure
 */
typedef struct sAPOLparm APOLparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (nosh.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct APOLparm
 *  @ingroup APOLparm
 *  @author  David Gohara
 *  @returns Newly allocated and initialized Vpmgp object
 */
VEXTERNC APOLparm* APOLparm_ctor();

/** @brief   FORTRAN stub to construct APOLparm
 *  @ingroup APOLparm
 *  @author  David Gohara
 *  @param   thee Pointer to allocated APOLparm object
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int APOLparm_ctor2(APOLparm *thee);

/** @brief   Object destructor
 *  @ingroup APOLparm
 *  @author  David Gohara
 *  @param   thee  Pointer to memory location of APOLparm object
 */
VEXTERNC void APOLparm_dtor(APOLparm **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup APOLparm
 *  @author  David Gohara
 *  @param   thee  Pointer to APOLparm object
 */
VEXTERNC void APOLparm_dtor2(APOLparm *thee);

/** 
 * @brief   Consistency check for parameter values stored in object
 * @ingroup APOLparm
 * @author  David Gohara
 * @param   thee   APOLparm object
 * @returns 1 if OK, 0 otherwise
 */
VEXTERNC int APOLparm_check(APOLparm *thee);

/**	@brief	Copy target object into thee
	@ingroup	APOLparm
	@author	Nathan Baker
	@param	thee	Destination object
	@param	source	Source object
*/
VEXTERNC void APOLparm_copy(APOLparm *thee, APOLparm *source);

/** 
 * @brief   Parse an MG keyword from an input file
 * @ingroup MGparm
 * @author  David Gohara
 * @param   thee   MGparm object 
 * @param   tok    Token to parse
 * @param   sock   Stream for more tokens
 * @return   1 if matched and assigned; -1 if matched, but there's
 * some sort of error (i.e., too few args); 0 if not matched
 */
VEXTERNC int APOLparm_parseToken(APOLparm *thee, char tok[VMAX_BUFSIZE], 
  Vio *sock);

#endif 

