/** @defgroup PBEparm PBEparm class
 *  @brief    Parameter structure for PBE variables independent of solver
 */

/**
 *  @file     pbeparm.h
 *  @ingroup  PBEparm
 *  @brief    Contains declarations for class PBEparm
 *  @version  $Id: pbeparm.h 1033 2006-12-29 17:08:22Z sobolevnrm $
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

#ifndef _PBEPARM_H_
#define _PBEPARM_H_

/* Generic headers */
#include "maloc/maloc.h"

/* Headers specific to this file */
#include "apbs/vhal.h"

/** @brief   Number of things that can be written out in a single calculation
 *  @ingroup PBEparm
 */
#define PBEPARM_MAXWRITE 10

/**
 * @ingroup PBEparm
 * @brief  Define energy calculation enumeration
 */
enum ePBEparm_calcEnergy {
    PCE_NO=0, /**< Do not perform energy calculation */
    PCE_TOTAL=1, /**< Calculate total energy only */
    PCE_COMPS=2 /**< Calculate per-atom energy components */
};

/**
 * @ingroup PBEparm
 * @brief  Define ePBEparm_calcEnergy enumeration as PBEparm_calcEnergy
 */
typedef enum ePBEparm_calcEnergy PBEparm_calcEnergy;

/**
 * @ingroup PBEparm
 * @brief  Define force calculation enumeration
 */
enum ePBEparm_calcForce {
    PCF_NO=0, /**< Do not perform force calculation */
    PCF_TOTAL=1, /**< Calculate total force only */
    PCF_COMPS=2 /**< Calculate per-atom force components */
};

/**
 * @ingroup PBEparm
 * @brief  Define ePBEparm_calcForce enumeration as PBEparm_calcForce
 */
typedef enum ePBEparm_calcForce PBEparm_calcForce;

/**
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @brief   Parameter structure for PBE variables from input files
 *  @note    If you add/delete/change something in this class, the member
 *           functions -- especially PBEparm_copy -- must be modified
 *           accordingly
 */
struct sPBEparm {

    int molid;  /**< Molecule ID to perform calculation on */
    int setmolid;  /**< Flag, @see molid */
    int useDielMap;  /**< Indicates whether we use an external
                      * dielectric maps (note plural) */
    int dielMapID;  /**< Dielectric map ID (if used) */
    int useKappaMap;  /**< Indicates whether we use an external
                       * kappa map */
    int kappaMapID;  /**< Kappa map ID (if used) */
    int useChargeMap;  /**< Indicates whether we use an external
                        * charge distribution map */
    int chargeMapID;  /**< Charge distribution map ID (if used) */
    Vhal_PBEType pbetype;  /**< Which version of the PBE are we solving? */
    int setpbetype;  /**< Flag, @see pbetype */
    Vbcfl bcfl;  /**< Boundary condition method */
    int setbcfl;  /**< Flag, @see bcfl */
    int nion;  /**< Number of counterion species */
    int setnion;  /**< Flag, @see nion */
    double ionq[MAXION];  /**< Counterion charges (in e) */
    double ionc[MAXION];  /**< Counterion concentrations (in M) */
    double ionr[MAXION];  /**< Counterion radii (in A) */
    int setion[MAXION];  /**< Flag, @see ionq */
    double pdie;  /**< Solute dielectric */
    int setpdie;  /**< Flag, @see pdie */
    double sdens; /**< Vacc sphere density */
    int setsdens; /**< Flag, @see sdens */
    double sdie;  /**< Solvent dielectric */
    int setsdie;  /**< Flag, @see sdie */
    Vsurf_Meth srfm;  /**< Surface calculation method */
    int setsrfm;  /**< Flag, @see srfm */
    double srad;  /**< Solvent radius */
    int setsrad;  /**< Flag, @see srad */
    double swin;  /**< Cubic spline window */
    int setswin;  /**< Flag, @see swin */
    double temp;  /**< Temperature (in K) */
    int settemp;  /**< Flag, @see temp */
    PBEparm_calcEnergy calcenergy;  /**< Energy calculation flag */
    int setcalcenergy;  /**< Flag, @see calcenergy */
    PBEparm_calcForce calcforce;  /**< Atomic forces calculation */
    int setcalcforce;  /**< Flag, @see calcforce */
    int numwrite;  /**< Number of write statements encountered */
    char writestem[PBEPARM_MAXWRITE][VMAX_ARGLEN]; /**< File stem to write 
                                                    * data to */
    Vdata_Type writetype[PBEPARM_MAXWRITE];  /**< What data to write */
    Vdata_Format writefmt[PBEPARM_MAXWRITE];  /**< File format to write data 
                                               * in */
    int writemat;  /**< Write out the operator matrix? 
                    * \li 0 => no 
                    * \li 1 => yes */
    int setwritemat;  /**< Flag, @see writemat */
    char writematstem[VMAX_ARGLEN];  /**< File stem to write mat */
    int writematflag;  /**< What matrix should we write:
                        * \li 0 => Poisson (differential operator)
                        * \li 1 => Poisson-Boltzmann operator linearized around
                        * solution (if applicable) */

    int parsed;  /**< Has this been filled with anything other
                  * than the default values? */
};

/** 
 *  @ingroup PBEparm
 *  @brief   Declaration of the PBEparm class as the PBEparm structure
 */
typedef struct sPBEparm PBEparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (mcsh.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Get charge (e) of specified ion species
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns Charge of ion species (e)
 */
VEXTERNC double PBEparm_getIonCharge(
        PBEparm *thee, /** PBEparm object */
        int iion  /** Ion species ID/index */
        );

/** @brief   Get concentration (M) of specified ion species
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns Concentration of ion species (M)
 */
VEXTERNC double PBEparm_getIonConc(
        PBEparm *thee, /** PBEparm object */
        int iion /** Ion species ID/index */
        );

/** @brief   Get radius (A) of specified ion species
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns Radius of ion species (A)
 */
VEXTERNC double PBEparm_getIonRadius(
        PBEparm *thee, /** PBEparm object */ 
        int iion /** Ion species ID/index */
        );


/** @brief   Construct PBEparm object
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns Newly allocated and initialized PBEparm object
 */
VEXTERNC PBEparm* PBEparm_ctor();

/** @brief   FORTRAN stub to construct PBEparm object
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns 1 if succesful, 0 otherwise
 */
VEXTERNC int PBEparm_ctor2(
        PBEparm *thee /** Memory location for object */
        );

/** @brief   Object destructor
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 */
VEXTERNC void PBEparm_dtor(
        PBEparm **thee /** Pointer to memory location of object */
        );

/** @brief   FORTRAN stub for object destructor
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 */
VEXTERNC void PBEparm_dtor2(
        PBEparm *thee /** Pointer to object to be destroyed */
        );

/** @brief   Consistency check for parameter values stored in object
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns 1 if OK, 0 otherwise
 */
VEXTERNC int PBEparm_check(
        PBEparm *thee /** Object to be checked */
        );

/** @brief   Copy PBEparm object into thee
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 */
VEXTERNC void PBEparm_copy(
        PBEparm *thee, /** Target for copy */
        PBEparm *parm /** Source for copy */
        );

/** @brief   Parse a keyword from an input file
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @return   1 if matched and assigned; -1 if matched, but there's some sort
 *            of error (i.e., too few args); 0 if not matched
 */
VEXTERNC int PBEparm_parseToken(
        PBEparm *thee, /** Parsing object */
        char tok[VMAX_BUFSIZE], /** Token to parse */ 
        Vio *sock /** Socket for additional tokens */
        );


#endif 

