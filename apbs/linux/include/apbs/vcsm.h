/** @defgroup Vcsm Vcsm class
 *  @brief  A charge-simplex map for evaluating integrals of delta functions
 *          in a finite element setting
 */

/**
 *  @file      vcsm.h
 *  @brief     Contains declarations for the Vcsm class
 *  @ingroup   Vcsm
 *  @version   $Id: vcsm.h 1033 2006-12-29 17:08:22Z sobolevnrm $
 *  @author    Nathan A. Baker
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

#ifndef _VCSM_H_
#define _VCSM_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/valist.h"

/* Specific headers */
#include "mc/mc.h"

/** @brief   External function for FEtk Gem class to use during mesh refinement
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 */
VEXTERNC void Gem_setExternalUpdateFunction(
        Gem *thee, /** The FEtk geometry managery */ 
        void (*externalUpdate)(SS **simps, int num) /** Function pointer for
                                                      call during mesh
                                                      refinement */
        );

/** @brief   Charge-simplex map class
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 */
struct sVcsm { 

  Valist *alist;      /**< Atom (charge) list */
  int natom;          /**< Size of thee->alist; redundant, but useful for
                       * convenience */
  Gem *gm;            /**< Grid manager (container class for master vertex
                       * and simplex lists as well as prolongation
                       * operator for updating after refinement ) */
  int **sqm;          /**< The map which gives the list charges associated with
                       * each simplex in gm->simplices.  The indices of
                       * the first dimension are associated with the
                       * simplex ID's in Vgm.  Each charge list (second 
                       * dimension) contains entries corresponding to
                       * indicies in thee->alist with lengths given in 
                       * thee->nsqm */
  int *nsqm;          /**< The length of the charge lists in thee->sqm */
  int nsimp;          /**< The _currently used) length of sqm, nsqm -- may not 
                       * always be up-to-date with Gem */
  int msimp;          /**< The maximum number of entries that can be 
                       * accomodated by sqm or nsqm  -- saves on realloc's */
  int **qsm;          /**< The inverse of sqm; the list of simplices
                       * associated with a given charge */
  int *nqsm;          /**< The length of the simplex lists in thee->qsm */
  int initFlag;       /**< Indicates whether the maps have been initialized
                       * yet */
  Vmem *vmem;         /**< Memory management object */

};

/** 
 *  @ingroup Vcsm
 *  @brief   Declaration of the Vcsm class as the Vcsm structure
 */
typedef struct sVcsm Vcsm;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Inlineable methods (vcsm.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VCSM)

    /** @brief   Get atom list
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Pointer to Valist atom list
     */
    VEXTERNC Valist* Vcsm_getValist(
            Vcsm *thee /** The Vcsm object */
            );

    /** @brief   Get number of atoms associated with a simplex
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Number of atoms associated with a simplex
     */
    VEXTERNC int Vcsm_getNumberAtoms(
            Vcsm *thee,  /** The Vcsm object */
            int isimp  /** Simplex ID */
            );

    /** @brief   Get particular atom associated with a simplex
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Array of atoms associated with a simplex
     */
    VEXTERNC Vatom* Vcsm_getAtom(
            Vcsm *thee,  /** The Vcsm object */
            int iatom,  /** Index of atom in Vcsm list ofr this simplex */
            int isimp  /** Simplex ID */
            );

    /** @brief   Get ID of particular atom in a simplex
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Index of atom in Valist object
     */
    VEXTERNC int Vcsm_getAtomIndex(
            Vcsm *thee,  /** The Vcsm object */
            int iatom,  /** Index of atom in Vcsm list for this simplex */
            int isimp  /** Simplex ID */
            );

    /** @brief   Get number of simplices associated with an atom
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Number of simplices associated with an atom
     */
    VEXTERNC int Vcsm_getNumberSimplices(
            Vcsm *thee,  /** The Vcsm object */
            int iatom  /** The Valist atom index */
            );

    /** @brief   Get particular simplex associated with an atom
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Pointer to simplex object
     */ 
    VEXTERNC SS* Vcsm_getSimplex(
            Vcsm *thee,  /** The Vcsm object */
            int isimp,  /** Index of simplex in Vcsm list */
            int iatom  /** Valist atom index */
            );

    /** @brief   Get index particular simplex associated with an atom
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Gem index of specified simplex
     */ 
    VEXTERNC int Vcsm_getSimplexIndex(
            Vcsm *thee,  /** THe Vcsm object */
            int isimp,  /** Index of simplex in Vcsm list */
            int iatom  /** Index of atom in Valist */
            );

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC unsigned long int Vcsm_memChk(
            Vcsm *thee /** The Vcsm object */
            );

#else /* if defined(VINLINE_VCSM) */
#   define Vcsm_getValist(thee) ((thee)->alist)
#   define Vcsm_getNumberAtoms(thee, isimp) ((thee)->nsqm[isimp])
#   define Vcsm_getAtom(thee, iatom, isimp) (Valist_getAtom((thee)->alist, ((thee)->sqm)[isimp][iatom]))
#   define Vcsm_getAtomIndex(thee, iatom, isimp) (((thee)->sqm)[isimp][iatom])
#   define Vcsm_getNumberSimplices(thee, iatom) (((thee)->nqsm)[iatom])
#   define Vcsm_getSimplex(thee, isimp, iatom) (Gem_SS((thee)->gm, ((thee)->qsm)[iatom][isimp]))
#   define Vcsm_getSimplexIndex(thee, isimp, iatom) (((thee)->qsm)[iatom][isimp])
#   define Vcsm_memChk(thee) (Vmem_bytes((thee)->vmem))
#endif /* if !defined(VINLINE_VCSM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Non-Inlineable methods (vcsm.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct Vcsm object
 *  @ingroup Vcsm 
 *  @author  Nathan Baker
 *  @note    \li The initial mesh must be sufficiently coarse for the assignment
 *             procedures to be efficient
 *           \li The map is not built until Vcsm_init is called
 *  @return  Pointer to newly allocated Vcsm object 
 */
VEXTERNC Vcsm* Vcsm_ctor(
        Valist *alist,  /** List of atoms */
        Gem *gm  /** FEtk geometry manager defining the mesh */
        );

/** @brief   FORTRAN stub to construct Vcsm object
 *  @ingroup Vcsm 
 *  @author  Nathan Baker
 *  @note    \li The initial mesh must be sufficiently coarse for the assignment
 *             procedures to be efficient
 *           \li The map is not built until Vcsm_init is called
 *  @return  1 if successful, 0 otherwise 
 */
VEXTERNC int Vcsm_ctor2(
        Vcsm *thee,  /** The Vcsm object */
        Valist *alist,  /** The list of atoms */
        Gem *gm  /** The FEtk geometry manager defining the mesh */
        );

/** @brief   Destroy Vcsm object
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 */
VEXTERNC void Vcsm_dtor(
        Vcsm **thee  /** Pointer to memory location for Vcsm object */
        );

/** @brief   FORTRAN stub to destroy Vcsm object
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 */
VEXTERNC void Vcsm_dtor2(
        Vcsm *thee /** Pointer to Vcsm object */
        );

/** @brief   Initialize charge-simplex map with mesh and atom data
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @note    The initial mesh must be sufficiently coarse for the assignment
 *            procedures to be efficient
 */
VEXTERNC void Vcsm_init(
        Vcsm *thee /** The Vcsm object */
        );

/** @brief   Update the charge-simplex and simplex-charge maps after
 *           refinement
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @return  1 if successful, 0 otherwise 
 */
VEXTERNC int Vcsm_update(
        Vcsm *thee, /** The Vcsm object */
        SS **simps, /** List of pointer to newly created (by refinement)
                      simplex objects.  The first simplex is expected to be
                      derived from the parent simplex and therefore have the
                      same ID.  The remaining simplices are the children and
                      should represent new entries in the charge-simplex map. */
        int num /** Number of simplices in simps list */
        );

#endif /* ifndef _VCSM_H_ */
