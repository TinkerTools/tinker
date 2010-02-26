/** @defgroup Vgrid Vgrid class
 *  @brief    Oracle for Cartesian mesh data
 */

/**
 *  @file    vgrid.h
 *  @ingroup Vgrid
 *  @author  Nathan Baker and Steve Bond
 *  @brief   Potential oracle for Cartesian mesh data
 *  @version $Id: vgrid.h 1033 2006-12-29 17:08:22Z sobolevnrm $
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

#ifndef _VGRID_H_
#define _VGRID_H_

#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/vstring.h"


/** @brief Number of decimal places for comparisons and formatting
 *  @ingroup Vgrid */
#define VGRID_DIGITS 6

/**
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @brief   Electrostatic potential oracle for Cartesian mesh data
 */
struct sVgrid {

    int nx;       /**< Number grid points in x direction */
    int ny;       /**< Number grid points in y direction */
    int nz;       /**< Number grid points in z direction */
    double hx;    /**< Grid spacing in x direction */
    double hy;    /**< Grid spacing in y direction */
    double hzed;  /**< Grid spacing in z direction */
    double xmin;  /**< x coordinate of lower grid corner */
    double ymin;  /**< y coordinate of lower grid corner */
    double zmin;  /**< z coordinate of lower grid corner */
    double xmax;  /**< x coordinate of upper grid corner */
    double ymax;  /**< y coordinate of upper grid corner */
    double zmax;  /**< z coordinate of upper grid corner */
    double *data; /**< nx*ny*nz array of data */
    int readdata; /**< flag indicating whether data was read from file */
    int ctordata; /**< flag indicating whether data was included at
                   *   construction */
    Vmem *mem;    /**< Memory manager object */
};

/** 
 *  @ingroup Vgrid
 *  @brief   Declaration of the Vgrid class as the sVgrid structure
 */
typedef struct sVgrid Vgrid;

#if !defined(VINLINE_VGRID)

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vgrid
     *  @author  Nathan Baker
     *  @param   thee  Vgrid object
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC unsigned long int Vgrid_memChk(Vgrid *thee);

#else /* if defined(VINLINE_VGRID) */

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vgrid
     *  @author  Nathan Baker
     *  @param   thee  Vgrid object
     *  @return  The memory used by this structure and its contents in bytes
     */
#   define Vgrid_memChk(thee) (Vmem_bytes((thee)->vmem))

#endif /* if !defined(VINLINE_VPMG) */

/** @brief   Construct Vgrid object with values obtained from Vpmg_readDX (for
 *           example)
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @param   nx    Number grid points in x direction
 *  @param   ny    Number grid points in y direction
 *  @param   nz    Number grid points in z direction
 *  @param   hx    Grid spacing in x direction
 *  @param   hy    Grid spacing in y direction
 *  @param   hzed  Grid spacing in z direction
 *  @param   xmin  x coordinate of lower grid corner
 *  @param   ymin  y coordinate of lower grid corner
 *  @param   zmin  z coordinate of lower grid corner
 *  @param   data  nx*ny*nz array of data.  This can be VNULL if you are
 *                 planning to read in data later with one of the read routines
 *  @returns Newly allocated and initialized Vgrid object
 */
VEXTERNC Vgrid*  Vgrid_ctor(int nx, int ny, int nz, 
                  double hx, double hy, double hzed,
                  double xmin, double ymin, double zmin,
                  double *data);

/** @brief   Initialize Vgrid object with values obtained from Vpmg_readDX (for
 *           example)
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @param   thee  Pointer to newly allocated Vgrid object
 *  @param   nx    Number grid points in x direction
 *  @param   ny    Number grid points in y direction
 *  @param   nz    Number grid points in z direction
 *  @param   hx    Grid spacing in x direction
 *  @param   hy    Grid spacing in y direction
 *  @param   hzed  Grid spacing in z direction
 *  @param   xmin  x coordinate of lower grid corner
 *  @param   ymin  y coordinate of lower grid corner
 *  @param   zmin  z coordinate of lower grid corner
 *  @param   data  nx*ny*nz array of data.  This can be VNULL if you are
 *                 planning to read in data later with one of the read routines
 *  @returns Newly allocated and initialized Vgrid object
 */
VEXTERNC int Vgrid_ctor2(Vgrid *thee, int nx, int ny, int nz, 
                  double hx, double hy, double hzed,
                  double xmin, double ymin, double zmin,
                  double *data);

/** @brief   Get potential value (from mesh or approximation) at a point
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @param   thee  Vgrid obejct
 *  @param   x     Point at which to evaluate potential
 *  @param   value Value of data at point x
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vgrid_value(Vgrid *thee, double x[3], double *value);

/** @brief   Object destructor
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void Vgrid_dtor(Vgrid **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void Vgrid_dtor2(Vgrid *thee);

/** @brief   Get second derivative values at a point
 *  @ingroup Vgrid
 *  @author  Steve Bond and Nathan Baker
 *  @param   thee   Pointer to Vgrid object
 *  @param   pt     Location to evaluate second derivative
 *  @param   cflag  
 *             \li  0:  Reduced Maximal Curvature
 *             \li  1:  Mean Curvature (Laplace)
 *             \li  2:  Gauss Curvature
 *             \li  3:  True Maximal Curvature
 *  @param   curv Specified curvature value
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vgrid_curvature(Vgrid *thee, double pt[3], int cflag, 
  double *curv);

/** @brief   Get first derivative values at a point
 *  @ingroup Vgrid
 *  @author  Nathan Baker and Steve Bond
 *  @param   thee   Pointer to Vgrid object
 *  @param   pt     Location to evaluate gradient
 *  @param   grad   Gradient
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vgrid_gradient(Vgrid *thee, double pt[3], double grad[3] );

/** @brief Write out the data in UHBD grid format 
 *  @note   \li The mesh spacing should be uniform
 *          \li Format changed from %12.6E to %12.5E
 * @ingroup Vgrid
 * @author  Nathan Baker
 * @param   thee   Grid object
 * @param   iodev  Output device type (FILE/BUFF/UNIX/INET)
 * @param   iofmt  Output device format (ASCII/XDR)
 * @param   thost  Output hostname (for sockets)
 * @param   fname  Output FILE/BUFF/UNIX/INET name
 * @param   title  Title to be inserted in grid file
 * @param   pvec   Partition weight (
 *                 if 1: point in current partition, 
 *                 if 0 point not in current partition
 *                 if > 0 && < 1 point on/near boundary )
 * @bug     This routine does not respect partition information
 */
VEXTERNC void Vgrid_writeUHBD(Vgrid *thee, const char *iodev, 
  const char *iofmt, const char *thost, const char *fname, char *title, 
  double *pvec);

/** @brief  Write out the data in OpenDX grid format 
 * @ingroup Vgrid
 * @author  Nathan Baker
 * @param   thee   Grid object
 * @param   iodev  Output device type (FILE/BUFF/UNIX/INET)
 * @param   iofmt  Output device format (ASCII/XDR)
 * @param   thost  Output hostname (for sockets)
 * @param   fname  Output FILE/BUFF/UNIX/INET name
 * @param   title  Title to be inserted in grid file
 * @param   pvec   Partition weight (
 *                 if 1: point in current partition, 
 *                 if 0 point not in current partition
 *                 if > 0 && < 1 point on/near boundary )
 */
VEXTERNC void Vgrid_writeDX(Vgrid *thee, const char *iodev, 
  const char *iofmt,  const char *thost, const char *fname, char *title, 
  double *pvec);

/** @brief   Read in data in OpenDX grid format
 *  @note    All dimension information is given in order: z, y, x
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @param   thee   Vgrid object
 *  @param   iodev  Input device type (FILE/BUFF/UNIX/INET)
 *  @param   iofmt  Input device format (ASCII/XDR)
 *  @param   thost  Input hostname (for sockets)
 *  @param   fname  Input FILE/BUFF/UNIX/INET name
 *  @returns 1 if sucessful, 0 otherwise
 */
VEXTERNC int Vgrid_readDX(Vgrid *thee, const char *iodev, const char *iofmt,
  const char *thost, const char *fname);

/**
 * @brief  Get the integral of the data
 * @ingroup  Vgrid
 * @author  Nathan Baker
 * @param  thee  Vgrid object
 * @returns  Integral of data */
VEXTERNC double Vgrid_integrate(Vgrid *thee);

/**
 * @ingroup  Vgrid
 * @author  Nathan Baker
 * @param  thee  Vgrid object
 * @returns  \f$L_1\f$ norm of data 
 * @brief  Get the \f$L_1\f$ norm of the data.  This returns the integral:
 *   \f[ \| u \|_{L_1} = \int_\Omega | u(x) | dx  \f] 
 */
VEXTERNC double Vgrid_normL1(Vgrid *thee);

/**
 * @ingroup  Vgrid
 * @author  Nathan Baker
 * @param  thee  Vgrid object
 * @returns  \f$L_2\f$ norm of data 
 * @brief  Get the \f$L_2\f$ norm of the data.  This returns the integral:
 *   \f[ \| u \|_{L_2} = \left( \int_\Omega | u(x) |^2 dx \right)^{1/2} \f] 
 */
VEXTERNC double Vgrid_normL2(Vgrid *thee);

/**
 * @ingroup  Vgrid
 * @author  Nathan Baker
 * @param  thee  Vgrid object
 * @returns  \f$L\infty\f$ norm of data 
 * @brief  Get the \f$L_\infty\f$ norm of the data.  This returns the integral:
 *   \f[ \| u \|_{L_\infty} = \sup_{x \in \Omega} | u(x) | \f]
 */
VEXTERNC double Vgrid_normLinf(Vgrid *thee);

/**
 * @ingroup  Vgrid
 * @author  Nathan Baker
 * @param  thee  Vgrid object
 * @returns  Integral of data 
 * @brief  Get the \f$H_1\f$ semi-norm of the data.  
 * This returns the integral:
 *   \f[ | u |_{H_1} = \left( \int_\Omega |\nabla u(x)|^2 dx \right)^{1/2} \f] 
 */
VEXTERNC double Vgrid_seminormH1(Vgrid *thee);

/**
 * @ingroup  Vgrid
 * @author  Nathan Baker
 * @param  thee  Vgrid object
 * @returns  Integral of data 
 * @brief  Get the \f$H_1\f$ norm (or energy norm) of the data.  
 * This returns the integral:
 *   \f[ \| u \|_{H_1} = \left( \int_\Omega |\nabla u(x)|^2 dx 
 *                     +        \int_\Omega |u(x)|^2 dx \right)^{1/2} \f] 
 */
VEXTERNC double Vgrid_normH1(Vgrid *thee);

#endif
