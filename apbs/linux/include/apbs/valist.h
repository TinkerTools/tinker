/** @defgroup Valist Valist class
 *  @brief    Container class for list of atom objects
 */

/** 
 *  @file     valist.h
 *  @ingroup  Valist
 *  @brief    Contains declarations for class Valist
 *  @version  $Id: valist.h 1046 2007-01-02 16:58:40Z sobolevnrm $
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

#ifndef _VALIST_H_
#define _VALIST_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

/* Headers specific to this file */
#include "apbs/vatom.h"
#include "apbs/vparam.h"

/** 
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @brief   Container class for list of atom objects
 */
struct sValist { 

  int number;         /**< Number of atoms in list */
  double center[3];   /**< Molecule center (xmin - xmax)/2, etc.*/
  double mincrd[3];   /**< Minimum coordinates */
  double maxcrd[3];   /**< Maximum coordinates */
  double maxrad;      /**< Maximum radius */
  double charge;      /**< Net charge */
  Vatom *atoms;       /**< Atom list */
  Vmem *vmem;         /**< Memory management object */

};

/**
 *  @ingroup Valist
 *  @brief Declaration of the Valist class as the Valist structure
 */
typedef struct sValist Valist;

#if !defined(VINLINE_VATOM)

/** 
 * @brief   Get actual array of atom objects from the list
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Array of atom objects 
 */
VEXTERNC Vatom* Valist_getAtomList(
        Valist *thee /**< Atom list object */
        );

/** @brief   Get x-coordinate of molecule center
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  X-coordinate of molecule center
 */
VEXTERNC double Valist_getCenterX(
        Valist *thee /**< Atom list object */
        );

/** @brief   Get y-coordinate of molecule center
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Y-coordinate of molecule center
 */ 
VEXTERNC double Valist_getCenterY(
        Valist *thee /**< Atom list object */
        );

/** @brief   Get z-coordinate of molecule center
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Z-coordinate of molecule center
 */ 
VEXTERNC double Valist_getCenterZ(
        Valist *thee /**< Atom list object */
        );

/** @brief   Get number of atoms in the list
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Number of atoms in list 
 */
VEXTERNC int Valist_getNumberAtoms(
        Valist *thee /**< Atom list object */
        );

/** @brief   Get pointer to particular atom in list
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Pointer to atom object i
 */
VEXTERNC Vatom* Valist_getAtom(
        Valist *thee, /**< Atom list object */
        int i /**< Index of atom in list */
        );

/** @brief   Get total memory allocated for this object and its members
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Total memory in bytes
 */
VEXTERNC unsigned long int Valist_memChk(
        Valist *thee /**< Atom list object */
        );

#else /* if defined(VINLINE_VATOM) */
#   define Valist_getAtomList(thee) ((thee)->atoms)
#   define Valist_getNumberAtoms(thee) ((thee)->number)
#   define Valist_getAtom(thee, i) (&((thee)->atoms[i]))
#   define Valist_memChk(thee) (Vmem_bytes((thee)->vmem))
#   define Valist_getCenterX(thee) ((thee)->center[0])
#   define Valist_getCenterY(thee) ((thee)->center[1])
#   define Valist_getCenterZ(thee) ((thee)->center[2])
#endif /* if !defined(VINLINE_VATOM) */

/** @brief   Construct the atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @returns Pointer to newly allocated (empty) atom list 
 */
VEXTERNC Valist* Valist_ctor();

/** @brief   FORTRAN stub to construct the atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Valist_ctor2(
        Valist *thee /**< Storage for new atom list */
        );

/** @brief   Destroys atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker
 */
VEXTERNC void Valist_dtor(
        Valist **thee /**< Pointer to storage for atom list */
        );

/** @brief   FORTRAN stub to destroy atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker
 */
VEXTERNC void Valist_dtor2(
        Valist *thee /**< Pointer to atom list object */
        );

/** 
 * @brief  Fill atom list with information from a PQR file
 * @ingroup Valist
 * @author  Nathan Baker
 * @return  1 if successful, 0 otherwise
 * @note  \li A PQR file has PDB structure with charge and radius in the last
 *            two columns instead of weight and occupancy
 *        \li We don't actually respect PDB format; instead recognize
 *            whitespace- or tab-delimited fields which allows us to deal with
 *            structures with coordinates > 999 or < -999.
 */
VEXTERNC int Valist_readPQR(
        Valist *thee, /**< Atom list object */
		Vparam *param, /**< A pre-initialized parameter object */
        Vio *sock /**< Socket reading for reading PQR file */
        );

/** 
 * @brief  Fill atom list with information from a PDB file
 * @ingroup Valist
 * @author  Nathan Baker, Todd Dolinsky
 * @return  1 if successful, 0 otherwise
 * @note  We don't actually respect PDB format; instead recognize whitespace-
 * or tab-delimited fields which allows us to deal with structures with
 * coordinates > 999 or < -999.
 */
VEXTERNC int Valist_readPDB(
        Valist *thee, /**< Atom list object */
        Vparam *param, /**< A pre-initialized parameter object */
        Vio *sock /**< Socket read for reading PDB file */
        );

/** 
 * @brief  Fill atom list with information from an XML file
 * @ingroup Valist
 * @author  Todd Dolinsky
 * @return  1 if successful, 0 otherwise
 * @note  \li The XML file must adhere to some guidelines, notably the 
 *            presence of an &lt;atom&gt; tag with all other useful information
 *            (x, y, z, charge, and radius) as nested elements.
 */
VEXTERNC int Valist_readXML(
        Valist *thee, /**< Atom list object */
		Vparam *param, /**< A pre-initialized parameter object */
        Vio *sock /**< Socket reading for reading PQR file */
        );

/** 
 * @brief   Load up Valist with various statistics
 * @ingroup Valist
 * @author  Nathan Baker
 * @return  1 if successful, 0 otherwise
 */
VEXTERNC int Valist_getStatistics(Valist *thee);


#endif /* ifndef _VALIST_H_ */
