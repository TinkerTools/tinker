/** @defgroup MGparm MGparm class
 *  @brief    Parameter which holds useful parameters for generic multigrid
 *            calculations
 */

/**
 *  @file     mgparm.h
 *  @ingroup  MGparm
 *  @brief    Contains declarations for class MGparm
 *  @version  $Id: mgparm.h 1033 2006-12-29 17:08:22Z sobolevnrm $
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


#ifndef _MGPARM_H_
#define _MGPARM_H_

#include "maloc/maloc.h"
#include "apbs/vhal.h"

/**
 * @brief  Calculation type 
 * @ingroup MGparm
 */
enum eMGparm_CalcType {
    MCT_MANUAL=0,  /**< mg-manual */
    MCT_AUTO=1,  /**< mg-auto */
    MCT_PARALLEL=2,  /**< mg-para */
    MCT_DUMMY=3,  /**< mg-dummy */
	MCT_NONE=4  /**< unspecified */
};

/**
 * @brief  Declare MGparm_CalcType type
 * @ingroup  MGparm
 */
typedef enum eMGparm_CalcType MGparm_CalcType;

/**
 * @brief  Centering method
 * @ingroup MGparm
 */
enum eMGparm_CentMeth {
    MCM_POINT=0, /**< Center on a point */ 
    MCM_MOLECULE=1,  /**< Center on a molecule */
	MCM_FOCUS=2  /**< Determined by focusing */
};

/**
 * @brief  Declare MGparm_CentMeth type
 * @ingroup  MGparm
 */
typedef enum eMGparm_CentMeth MGparm_CentMeth;
/**
 *  @ingroup MGparm
 *  @author  Nathan Baker and Todd Dolinsky
 *  @brief   Parameter structure for MG-specific variables from input files
 *  @note    If you add/delete/change something in this class, the member
 *           functions -- especially MGparm_copy -- must be modified
 *           accordingly
 */
struct sMGparm {

    MGparm_CalcType type;  /**< What type of MG calculation? */
    int parsed;  /**< Has this structure been filled? (0 = no, 1 = yes) */

    /* *** GENERIC PARAMETERS *** */
    int dime[3];  /**< Grid dimensions */
    int setdime;  /**< Flag, @see dime */
    Vchrg_Meth chgm;  /**< Charge discretization method */
    int setchgm;  /**< Flag, @see chgm */
    Vchrg_Src  chgs; /**< Charge source (Charge, Multipole, Induced Dipole, 
                      * NL Induced */

    /* *** TYPE 0 PARAMETERS (SEQUENTIAL MANUAL) *** */
    int nlev;  /**< Levels in multigrid hierarchy 
                *   @deprecated Just ignored now */
    int setnlev;  /**< Flag, @see nlev */
    double grid[3];  /**< Grid spacings */
    int setgrid;  /**< Flag, @see grid */
    double glen[3];  /**< Grid side lengths. */
    int setglen;  /**< Flag, @see glen */
    MGparm_CentMeth cmeth;  /**< Centering method */  
    double center[3];  /**< Grid center. If ispart = 0, then this is
                        * only meaningful if cmeth = 0.  However, if
                        * ispart = 1 and cmeth = MCM_PNT, then this is the
                        * center of the non-disjoint (overlapping)
                        * partition.  If ispart = 1 and cmeth = MCM_MOL, then
                        * this is the vector that must be added to the
                        * center of the molecule to give the center of
                        * the non-disjoint partition.  */
    int centmol;  /**< Particular molecule on which we want to center the grid.  
		This should be the appropriate index in an array of molecules, not the 
		positive definite integer specified by the user. */
    int setgcent;  /**< Flag, @see cmeth */

    /* ******** TYPE 1 & 2 PARAMETERS (SEQUENTIAL & PARALLEL AUTO-FOCUS) *** */
    double cglen[3];  /**< Coarse grid side lengths */
    int setcglen;  /**< Flag, @see cglen */
    double fglen[3];  /**< Fine grid side lengths */
    int setfglen;  /**< Flag, @see fglen */
    MGparm_CentMeth ccmeth;  /**< Coarse grid centering method */
    double ccenter[3];  /**< Coarse grid center.  */
    int ccentmol;  /**< Particular molecule on which we want to center the grid.  
		This should be the appropriate index in an array of molecules, not the 
		positive definite integer specified by the user. */
    int setcgcent;  /**< Flag, @see ccmeth */
    MGparm_CentMeth fcmeth;  /**< Fine grid centering method */
    double fcenter[3];  /**< Fine grid center.  */
    int fcentmol; /**< Particular molecule on which we want to center the grid.  
		This should be the appropriate index in an array of molecules, not the 
		positive definite integer specified by the user. */
    int setfgcent;  /**< Flag, @see fcmeth */


    /* ********* TYPE 2 PARAMETERS (PARALLEL AUTO-FOCUS) ******** */
    double partDisjCenter[3];  /**< This gives the center
                                     of the disjoint partitions */
    double partDisjLength[3];  /**< This gives the lengths of the disjoint 
                                * partitions */
    int partDisjOwnSide[6];  /**< Tells whether the boundary points are ours 
                              * (1) or not (0) */

    int pdime[3];  /**< Grid of processors to be used in calculation */
    int setpdime;  /**< Flag, @see pdime */
    int proc_rank;  /**< Rank of this processor */
    int setrank;  /**< Flag, @see proc_rank */
    int proc_size;  /**< Total number of processors */
    int setsize;  /**< Flag, @see proc_size */
    double ofrac;  /**< Overlap fraction between procs */
    int setofrac;  /**< Flag, @see ofrac */
    int async; /**< Processor ID for asynchronous calculation */
    int setasync; /**< Flag, @see asynch */
};

/** @typedef MGparm
 *  @ingroup MGparm
 *  @brief   Declaration of the MGparm class as the MGparm structure
 */
typedef struct sMGparm MGparm;

/** @brief   Get number of grid points in x direction
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Number of grid points in the x direction
 */
VEXTERNC int MGparm_getNx(MGparm *thee);

/** @brief   Get number of grid points in y direction
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Number of grid points in the y direction
 */
VEXTERNC int MGparm_getNy(MGparm *thee);

/** @brief   Get number of grid points in z direction
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Number of grid points in the z direction
 */
VEXTERNC int MGparm_getNz(MGparm *thee);

/** @brief   Get grid spacing in x direction (&Aring;)
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Grid spacing in the x direction
 */
VEXTERNC double MGparm_getHx(MGparm *thee);

/** @brief   Get grid spacing in y direction (&Aring;)
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Grid spacing in the y direction
 */
VEXTERNC double MGparm_getHy(MGparm *thee);

/** @brief   Get grid spacing in z direction (&Aring;)
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Grid spacing in the z direction
 */
VEXTERNC double MGparm_getHz(MGparm *thee);

/** @brief   Set center x-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @param   x     x-coordinate
 */
VEXTERNC void MGparm_setCenterX(MGparm *thee, double x);

/** @brief   Set center y-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @param   y     y-coordinate
 */  
VEXTERNC void MGparm_setCenterY(MGparm *thee, double y);

/** @brief   Set center z-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @param   z     z-coordinate
 */
VEXTERNC void MGparm_setCenterZ(MGparm *thee, double z);

/** @brief   Get center x-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  x-coordinate
 */
VEXTERNC double MGparm_getCenterX(MGparm *thee);

/** @brief   Get center y-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  y-coordinate
 */
VEXTERNC double MGparm_getCenterY(MGparm *thee);

/** @brief   Get center z-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  z-coordinate
 */
VEXTERNC double MGparm_getCenterZ(MGparm *thee);

/** @brief   Construct MGparm object
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   type Type of MG calculation
 *  @returns Newly allocated and initialized MGparm object
 */
VEXTERNC MGparm*  MGparm_ctor(MGparm_CalcType type);

/** @brief   FORTRAN stub to construct MGparm object
 *  @ingroup MGparm
 *  @author  Nathan Baker and Todd Dolinsky
 *  @param   thee Space for MGparm object
 *  @param   type Type of MG calculation
 *  @returns 1 if succesful, 0 otherwise
 */
VEXTERNC int      MGparm_ctor2(MGparm *thee, MGparm_CalcType type);

/** @brief   Object destructor
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of MGparm object
 */
VEXTERNC void     MGparm_dtor(MGparm **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to MGparm object
 */
VEXTERNC void     MGparm_dtor2(MGparm *thee);

/** @brief   Consistency check for parameter values stored in object
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee   MGparm object
 *  @returns 1 if OK, 0 otherwise
 */
VEXTERNC int      MGparm_check(MGparm *thee);

/** @brief   Copy MGparm object into thee
 *  @ingroup MGparm
 *  @author  Nathan Baker and Todd Dolinsky
 *  @param   thee   MGparm object (target for copy)
 *  @param   parm   MGparm object (source for copy)
 */
VEXTERNC void     MGparm_copy(MGparm *thee, MGparm *parm);

/** @brief   Parse an MG keyword from an input file
 *  @ingroup MGparm
 *  @author  Nathan Baker and Todd Dolinsky
 *  @param   thee   MGparm object 
 *  @param   tok    Token to parse
 *  @param   sock   Stream for more tokens
 *  @return   1 if matched and assigned; -1 if matched, but there's some sort
 *            of error (i.e., too few args); 0 if not matched
 */
VEXTERNC int      MGparm_parseToken(MGparm *thee, char tok[VMAX_BUFSIZE], 
                    Vio *sock);

#endif 

