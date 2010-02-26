/** @defgroup Vpee Vpee class
 *  @brief  This class provides some functionality for error esimation
 *          in parallel. 
 * 
 *    This class provides some functionality for error esimation in parallel.
 *    The purpose is to modulate the error returned by some external error
 *    estimator according to the partitioning of the mesh.  For example, the
 *    Bank/Holst parallel refinement routine essentially reduces the error
 *    outside the ``local" partition to zero.  However,  this leads to the need
 *    for a few final overlapping Schwarz solves to smooth out the errors near
 *    partition boundaries.  Supposedly, if the region in which we allow
 *    error-based refinement includes the ``local" partition and an external
 *    buffer zone approximately equal in size to the local region, then the
 *    solution will asymptotically approach the solution obtained via more
 *    typical methods.  This is essentially a more flexible parallel
 *    implementation of MC's AM_markRefine.
 */

/**
 *  @file     vpee.h
 *  @ingroup  Vpee
 *  @brief    Contains declarations for class Vpee
 *  @version  $Id: vpee.h 1033 2006-12-29 17:08:22Z sobolevnrm $
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

#ifndef _VPEE_H
#define _VPEE_H

/* Generic headers */
#include "maloc/maloc.h"
#include "mc/mc.h"

/**
 *  @ingroup Vpee
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vpee class/module
 */
struct sVpee {

  Gem *gm;                     /**< Grid manager */
  int localPartID;             /**< The local partition ID: i.e. the partition 
                                * whose boundary simplices we're keeping
                                * track of */
  double localPartCenter[3];   /**< The coordinates of the center of the local
                                * partition */
  double localPartRadius;      /**< The radius of the circle/sphere which
                                * circumscribes the local partition */
  int killFlag;                /**< A flag indicating the method we're using to
                                * artificially decrease the error esimate
                                * outside the local partition */
  double killParam;            /**< A parameter for the error estimate
                                * attenuation method */
  Vmem *mem;                   /**< Memory manager */

};

/** 
 *  @ingroup Vpee
 *  @brief   Declaration of the Vpee class as the Vpee structure
 */
typedef struct sVpee Vpee;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpee Inlineable methods 
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPEE)
#else /* if defined(VINLINE_VPEE) */
#endif /* if !defined(VINLINE_VPEE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpee: Non-Inlineable methods (vpee.c)
/////////////////////////////////////////////////////////////////////////// */

/** 
 * @brief   Construct the Vpee object
 * @ingroup Vpee
 * @author  Nathan Baker
 * @return   Newly constructed Vpee object
 */
VEXTERNC Vpee* Vpee_ctor(
        Gem *gm,  /** FEtk geometry manager object */
        int localPartID,  /** ID of the local partition (focus of refinement) */
        int killFlag,  /** A flag to indicate how error estimates are to be
                         attenuated outside the local partition:
                         \li 0:  no attenuation
                         \li 1:  all error outside the local partition set to
                               zero
                         \li 2:  all error is set to zero outside a sphere of
                               radius (killParam*partRadius), where
                               partRadius is the radius of the sphere
                               circumscribing the local partition
                         \li 3:  all error is set to zero except for the local
                               partition and its immediate neighbors */
        double killParam /** @see killFlag for usage */
        );

/** 
 * @brief  FORTRAN stub to construct the Vpee object
 * @ingroup  Vpee
 * @author  Nathan Baker
 * @return  1 if successful, 0 otherwise
 */
VEXTERNC int Vpee_ctor2(
        Vpee *thee,  /** The Vpee object */
        Gem *gm,  /** FEtk geometry manager object */
        int localPartID,  /** ID of the local partition (focus of refinement) */
        int killFlag,  /** A flag to indicate how error estimates are to be
                         attenuated outside the local partition:
                         \li 0:  no attenuation
                         \li 1:  all error outside the local partition set to
                               zero
                         \li 2:  all error is set to zero outside a sphere of
                               radius (killParam*partRadius), where
                               partRadius is the radius of the sphere
                               circumscribing the local partition
                         \li 3:  all error is set to zero except for the local
                               partition and its immediate neighbors */
        double killParam /** @see killFlag for usage */
        );

/** @brief   Object destructor
 *  @ingroup Vpee
 *  @author  Nathan Baker
 */
VEXTERNC void Vpee_dtor(
        Vpee **thee /** Pointer to memory location of the Vpee object */
        );

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vpee
 *  @author  Nathan Baker
 */
VEXTERNC void Vpee_dtor2(
        Vpee *thee /** Pointer to object to be destroyed */
        );

/** @brief   Mark simplices for refinement based on attenuated error estimates.
 *  
 *  A wrapper/reimplementation of AM_markRefine that allows for more flexible
 *  attenuation of error-based markings outside the local partition.  The error
 *  in each simplex is modified by the method (see killFlag) specified in the
 *  Vpee constructor.  This allows the user to confine refinement to an
 *  arbitrary area around the local partition.
 * 
 *  @ingroup Vpee
 *  @author  Nathan Baker and Mike Holst
 *  @note  This routine borrows very heavily from FEtk routines by Mike Holst.
 *  @return The number of simplices marked for refinement.
 *  @bug  This function is no longer up-to-date with FEtk and may not function
 *  properly
 */
VEXTERNC int Vpee_markRefine(
        Vpee *thee,  /** The Vpee object */
        AM *am,  /** The FEtk algebra manager currently used to solve the PB */
        int level,  /** The current level of the multigrid hierarchy */
        int akey,  /** The marking method: 
                      \li -1:  Reset markings  --> killFlag has no effect.
                      \li 0:  Uniform.
                      \li 1:  User defined (geometry-based).
                      \li >1:  A numerical estimate for the error has already been
                               set in am and should be attenuated according to
                               killFlag and used, in conjunction with etol, to mark
                               simplices for refinement. */
        int rcol, /** The ID of the main parition on which to mark (or -1 if
                    all partitions should be marked).  NOte that we shouldhave
                    (rcol == thee->localPartID) for (thee->killFlag == 2 or 3) */
        double etol,  /** The error tolerance criterion for marking */
        int bkey  /** How the error tolerance is interpreted:
                     \li 0:  Simplex marked if error > etol.
                     \li 1:  Simplex marked if error > 
                     sqrt(etol^2/L) where L$ is the number of simplices */
        );

/** @brief   Returns the number of simplices in the local partition
 *  @ingroup Vpee
 *  @author  Nathan Baker
 *  @return  Number of simplices in the local partition
 */
VEXTERNC int Vpee_numSS(
        Vpee *thee /** The Vpee object */
        );

#endif    /* ifndef _VPEE_H_ */
