/** @defgroup Vmgrid Vmgrid class
 *  @brief    Oracle for Cartesian mesh data
 */

/**
 *  @file    vmgrid.h
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @brief   Multiresolution oracle for Cartesian mesh data
 *  @version $Id: vmgrid.h 1033 2006-12-29 17:08:22Z sobolevnrm $
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

#ifndef _VMGRID_H_
#define _VMGRID_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

/* Headers specific to this file */
#include "apbs/vgrid.h"

/** @def VMGRIDMAX   The maximum number of levels in the grid hiearchy
 *  @ingroup Vmgrid
 */
#define VMGRIDMAX 20


/**
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @brief   Multiresoltion oracle for Cartesian mesh data
 */
struct sVmgrid {

    int ngrids;                /**< Number of grids in hiearchy */
    Vgrid *grids[VMGRIDMAX];   /**< Grids in hiearchy.  Our convention will be
                                *   to have the finest grid first, however,
                                *   this will not be enforced as it may be
                                *   useful to search multiple grids for
                                *   parallel datasets, etc. */
};

/** 
 *  @ingroup Vmgrid
 *  @brief   Declaration of the Vmgrid class as the Vgmrid structure
 */
typedef struct sVmgrid Vmgrid;

/** @brief   Construct Vmgrid object 
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @returns Newly allocated and initialized Vmgrid object
 */
VEXTERNC Vmgrid*  Vmgrid_ctor();

/** @brief   Initialize Vmgrid object 
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee Newly allocated Vmgrid object
 *  @returns Newly allocated and initialized Vmgrid object
 */
VEXTERNC int Vmgrid_ctor2(Vmgrid *thee);

/** @brief   Get potential value (from mesh or approximation) at a point
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee  Vmgrid obejct
 *  @param   x     Point at which to evaluate potential
 *  @param   value Value of data at point x
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vmgrid_value(Vmgrid *thee, double x[3], double *value);

/** @brief   Object destructor
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void Vmgrid_dtor(Vmgrid **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void Vmgrid_dtor2(Vmgrid *thee);

/** @brief   Add a grid to the hierarchy
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 *  @param   grid   Grid to be added.  As mentioned above, we would prefer to
 *           have the finest grid added first, next-finest second, ...,
 *           coarsest last -- this is how the grid will be searched when
 *           looking up values for points.  However, this is not enforced to
 *           provide flexibility for cases where the dataset is decomposed into
 *           disjoint partitions, etc.
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vmgrid_addGrid(Vmgrid *thee, Vgrid *grid);


/** @brief   Get second derivative values at a point
 *  @ingroup Vmgrid
 *  @author  Nathan Baker (wrapper for Vgrid routine by Steve Bond)
 *  @param   thee   Pointer to Vmgrid object
 *  @param   pt     Location to evaluate second derivative
 *  @param   cflag  
 *             \li  0:  Reduced Maximal Curvature
 *             \li  1:  Mean Curvature (Laplace)
 *             \li  2:  Gauss Curvature
 *             \li  3:  True Maximal Curvature
 *  @param   curv Specified curvature value
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vmgrid_curvature(Vmgrid *thee, double pt[3], int cflag, 
  double *curv);

/** @brief   Get first derivative values at a point
 *  @ingroup Vmgrid
 *  @author  Nathan Baker and Steve Bond
 *  @param   thee   Pointer to Vmgrid object
 *  @param   pt     Location to evaluate gradient
 *  @param   grad   Gradient
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vmgrid_gradient(Vmgrid *thee, double pt[3], double grad[3] );

/** @brief   Get specific grid in hiearchy
 *  @ingroup Vmgrid
 *  @author  Nathan Baker 
 *  @param   thee   Pointer to Vmgrid object
 *  @param   num    Number of grid in hiearchy 
 *  @return  Pointer to specified grid
 */
VEXTERNC Vgrid* Vmgrid_getGridByNum(Vmgrid *thee, int num);

/** @brief   Get grid in hiearchy which contains specified point or VNULL
 *  @ingroup Vmgrid
 *  @author  Nathan Baker 
 *  @param   thee   Pointer to Vmgrid object
 *  @param   pt     Point to check
 *  @return  Pointer to specified grid
 */
VEXTERNC Vgrid* Vmgrid_getGridByPoint(Vmgrid *thee, double pt[3]);

#endif

