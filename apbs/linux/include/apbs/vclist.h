/** @defgroup Vclist Vclist class
 *  @brief    Atom cell list
 */

/**
 *  @file     vclist.h
 *  @ingroup  Vclist
 *  @brief    Contains declarations for class Vclist
 *  @version  $Id: vclist.h 1033 2006-12-29 17:08:22Z sobolevnrm $
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

#ifndef _VCLIST_H_
#define _VCLIST_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

/* Headers specific to this file */
#include "apbs/valist.h"
#include "apbs/vatom.h"
#include "apbs/vunit.h"

/** 
 *  @ingroup Vclist
 *  @author  Nathan Baker
 *  @brief   Atom cell list domain setup mode
 */
enum eVclist_DomainMode {
    CLIST_AUTO_DOMAIN,  /**< Setup the cell list domain automatically to
                         * encompass the entire molecule */
    CLIST_MANUAL_DOMAIN   /**< Specify the cell list domain manually through
                           * the constructor */
};

/** 
 * @typedef Vclist_DomainMode
 * @ingroup Vclist
 * @brief Declaration of Vclist_DomainMode enumeration type
 */
typedef enum eVclist_DomainMode Vclist_DomainMode;

/**
 * @ingroup Vclist
 * @author Nathan Baker
 * @brief Atom cell list cell
 */
struct sVclistCell {
    Vatom **atoms;  /**< Array of atom objects associated with this cell */
    int natoms;  /**< Length of thee->atoms array */
};

/** 
 *  @ingroup Vclist
 *  @brief   Declaration of the VclistCell class as the VclistCell structure
 */
typedef struct sVclistCell VclistCell;

/**
 *  @ingroup Vclist
 *  @author  Nathan Baker
 *  @brief   Atom cell list
 */
struct sVclist {

  Vmem *vmem;  /**< Memory management object for this class */
  Valist *alist;  /**< Original Valist structure for list of atoms */
  Vclist_DomainMode mode;  /**< How the cell list was constructed */
  int npts[VAPBS_DIM];  /**< Hash table grid dimensions */
  int n;  /**< n = nx*nz*ny */
  double max_radius;  /**< Maximum probe radius */
  VclistCell *cells;  /**< Cell array of length thee->n */
  double lower_corner[VAPBS_DIM]; /**< Hash table grid corner */
  double upper_corner[VAPBS_DIM]; /**< Hash table grid corner */
  double spacs[VAPBS_DIM];  /**< Hash table grid spacings */

};

/** 
 *  @ingroup Vclist
 *  @brief   Declaration of the Vclist class as the Vclist structure
 */
typedef struct sVclist Vclist;

#if !defined(VINLINE_VCLIST)

    /** @brief   Get number of bytes in this object and its members
     *  @ingroup Vclist
     *  @author  Nathan Baker
     *  @returns Number of bytes allocated for object
     */
    VEXTERNC unsigned long int Vclist_memChk(
            Vclist *thee /** Object for memory check */
            );

    /**
     * @brief  Get the max probe radius value (in A) the cell list was
     *         constructed with
     * @ingroup Vclist
     * @author Nathan Baker
     * @returns Max probe radius (in A)
     */
    VEXTERNC double Vclist_maxRadius(
            Vclist *thee /** Cell list object */
            );

#else /* if defined(VINLINE_VCLIST) */

#   define Vclist_memChk(thee) (Vmem_bytes((thee)->vmem))
#   define Vclist_maxRadius(thee) ((thee)->max_radius)

#endif /* if !defined(VINLINE_VCLIST) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vclist: Non-Inlineable methods (vclist.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct the cell list object
 *  @ingroup Vclist
 *  @author  Nathan Baker
 *  @returns Newly allocated Vclist object */
VEXTERNC Vclist* Vclist_ctor(
        Valist *alist, /** Molecule for cell list queries */
        double max_radius, /** Max probe radius (&Aring;) to be queried */ 
        int npts[VAPBS_DIM], /** Number of in hash table points in each
                              * direction*/ 
        Vclist_DomainMode mode, /** Mode to construct table */
        double lower_corner[VAPBS_DIM],  /** Hash table lower corner for
                                           manual construction (see mode
                                           variable); ignored otherwise */
        double upper_corner[VAPBS_DIM]   /** Hash table upper corner for
                                           manual construction (see mode
                                           variable); ignored otherwise */
        );

/** @brief   FORTRAN stub to construct the cell list object
 *  @ingroup Vclist
 *  @author  Nathan Baker
 *  @returns 1 if successful, 0 otherwise */
VEXTERNC int Vclist_ctor2(
        Vclist *thee, /** Memory for Vclist objet */
        Valist *alist, /** Molecule for cell list queries */
        double max_radius, /** Max probe radius (&Aring;) to be queried */ 
        int npts[VAPBS_DIM], /** Number of in hash table points in each
                              * direction*/ 
        Vclist_DomainMode mode, /** Mode to construct table */
        double lower_corner[VAPBS_DIM],  /** Hash table lower corner for
                                           manual construction (see mode
                                           variable); ignored otherwise */
        double upper_corner[VAPBS_DIM]   /** Hash table upper corner for
                                           manual construction (see mode
                                           variable); ignored otherwise */
        );

/** @brief   Destroy object
 *  @ingroup Vclist
 *  @author  Nathan Baker
 */
VEXTERNC void Vclist_dtor(
        Vclist **thee /** Pointer to memory location of object */
        );

/** @brief   FORTRAN stub to destroy object
 *  @ingroup Vclist
 *  @author  Nathan Baker
 */
VEXTERNC void Vclist_dtor2(
        Vclist *thee /** Pointer to object */
        );

/** 
 * @brief  Return cell corresponding to specified position or return VNULL.
 * @ingroup Vclist
 * @author Nathan Baker
 * @returns Pointer to VclistCell object or VNULL if no cell available (away
 * from molecule).
 */
VEXTERNC VclistCell* Vclist_getCell(
        Vclist *thee, /** Pointer to Vclist cell list */
        double position[VAPBS_DIM] /** Position to evaluate */
        );

/**
 * @brief  Allocate and construct a cell list cell object
 * @ingroup Vclist
 * @author Nathan Baker
 * @returns Pointer to newly-allocated and constructed object.
 */
VEXTERNC VclistCell* VclistCell_ctor(
        int natoms  /** Number of atoms associated with this cell */
        );

/**
 * @brief  Construct a cell list object
 * @ingroup  Vclist
 * @author  Nathan Baker
 * @returns 1 if successful, 0 otherwise
 */
VEXTERNC int VclistCell_ctor2(
        VclistCell *thee,  /** Memory location for object */
        int natoms  /** Number of atoms associated with this cell */
        );

/** @brief   Destroy object
 *  @ingroup Vclist
 *  @author  Nathan Baker
 */
VEXTERNC void VclistCell_dtor(
        VclistCell **thee /** Pointer to memory location of object */
        );

/** @brief   FORTRAN stub to destroy object
 *  @ingroup Vclist
 *  @author  Nathan Baker
 */
VEXTERNC void VclistCell_dtor2(
        VclistCell *thee /** Pointer to object */
        );

#endif    /* ifndef _VCLIST_H_ */
