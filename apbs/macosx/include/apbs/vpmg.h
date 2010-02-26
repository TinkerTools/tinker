/** @defgroup Vpmg Vpmg class
 *  @brief  A wrapper for Mike Holst's PMG multigrid code.  
 *  @note   Many of the routines and macros are borrowed from the main.c driver
 *          (written by Mike Holst) provided with the PMG code.
 */

/**
 *  @file     vpmg.h
 *  @ingroup  Vpmg
 *  @brief    Contains declarations for class Vpmg
 *  @version  $Id: vpmg.h 1350 2009-02-12 00:38:48Z yhuang01 $
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
 * Copyright (c) 2002-2009, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2009.  Nathan A. Baker
 * Portions Copyright (c) 1999-2002.  The Regents of the University of California.
 * Portions Copyright (c) 1995.  Michael Holst
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  
 * 
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */


#ifndef _VPMG_H_
#define _VPMG_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

/* Headers specific to this file */
#include "apbs/vpmgp.h"
#include "apbs/vacc.h"
#include "apbs/vcap.h"
#include "apbs/vpbe.h"
#include "apbs/vgrid.h"
#include "apbs/mgparm.h"
#include "apbs/pbeparm.h"

/** @def VPMGMAXPART The maximum number of partitions the
 *                   mesh can be divided into 
 *  @ingroup Vpmg
 */
#define VPMGMAXPART 2000  

/** 
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vpmg class/module
 *      
 *  Many of the routines and macros are borrowed from the main.c driver 
 *  (written by Mike Holst) provided with the PMG code.
 *     
 */
struct sVpmg {

  Vmem *vmem;  /**< Memory management object for this class */
  Vpmgp *pmgp;  /**< Parameters */
  Vpbe *pbe;  /**< Information about the PBE system */

  double *epsx;  /**< X-shifted dielectric map */
  double *epsy;  /**< Y-shifted dielectric map */
  double *epsz;  /**< Y-shifted dielectric map */
  double *kappa;  /**< Ion accessibility map (0 <= kappa(x) <= 1) */
  double *charge;  /**< Charge map */

  int *iparm;  /**< Passing int parameters to FORTRAN */
  double *rparm;  /**< Passing real parameters to FORTRAN */
  int *iwork;  /**< Work array */
  double *rwork;  /**< Work array */
  double *a1cf;  /**< Operator coefficient values (a11) -- this array can be
                  * overwritten */
  double *a2cf;  /**< Operator coefficient values (a22) -- this array can be
                   overwritten */
  double *a3cf;  /**< Operator coefficient values (a33) -- this array can be
                   overwritten */
  double *ccf;  /**< Helmholtz term -- this array can be overwritten */
  double *fcf;  /**< Right-hand side -- this array can be overwritten */
  double *tcf;  /**< True solution */
  double *u;  /**< Solution */
  double *xf;  /**< Mesh point x coordinates */
  double *yf;  /**< Mesh point y coordinates */
  double *zf;  /**< Mesh point z coordinates */
  double *gxcf;  /**< Boundary conditions for x faces */
  double *gycf;  /**< Boundary conditions for y faces */
  double *gzcf;  /**< Boundary conditions for z faces */
  double *pvec;  /**< Partition mask array */
  double extDiEnergy;  /**< Stores contributions to the dielectric energy from
                        * regions outside the problem domain */
  double extQmEnergy;  /**< Stores contributions to the mobile ion energy from
                        * regions outside the problem domain */
  double extQfEnergy;  /**< Stores contributions to the fixed charge energy
                        * from regions outside the problem domain */
  double extNpEnergy;  /**< Stores contributions to the apolar energy from
                        * regions outside the problem domain */
  Vsurf_Meth surfMeth;  /**< Surface definition method */
  double splineWin;  /**< Spline window parm for surf defs */
  Vchrg_Meth chargeMeth;  /**< Charge discretization method */
  Vchrg_Src chargeSrc;  /**< Charge source */
  
  int filled;  /**< Indicates whether Vpmg_fillco has been called */

  int useDielXMap;  /**< Indicates whether Vpmg_fillco was called with an
                      external x-shifted dielectric map */
  Vgrid *dielXMap;  /**< External x-shifted dielectric map */
  int useDielYMap;  /**< Indicates whether Vpmg_fillco was called with an
                     * external y-shifted dielectric map */
  Vgrid *dielYMap;  /**< External y-shifted dielectric map */
  int useDielZMap;  /**< Indicates whether Vpmg_fillco was called with an
                     * external z-shifted dielectric map */
  Vgrid *dielZMap;  /**< External z-shifted dielectric map */
  int useKappaMap;  /**< Indicates whether Vpmg_fillco was called with an
                     * external kappa map */
  Vgrid *kappaMap;  /**< External kappa map */
  int useChargeMap;  /**< Indicates whether Vpmg_fillco was called with an
                      * external charge distribution map */
  Vgrid *chargeMap;  /**< External charge distribution map */
};

/** 
 *  @ingroup Vpmg
 *  @brief   Declaration of the Vpmg class as the Vpmg structure
 */
typedef struct sVpmg Vpmg;

/* /////////////////////////////////////////////////////////////////////////
/// Inlineable methods
//////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VPMG)

    /** @brief   Return the memory used by this structure (and its contents) 
     *           in bytes
     *  @ingroup Vpmg
     *  @author  Nathan Baker
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC unsigned long int Vpmg_memChk(
            Vpmg *thee  /** Object for memory check */
            );

#else /* if defined(VINLINE_VPMG) */

#   define Vpmg_memChk(thee) (Vmem_bytes((thee)->vmem))

#endif /* if !defined(VINLINE_VPMG) */

/* /////////////////////////////////////////////////////////////////////////
/// Non-inlineable methods
//////////////////////////////////////////////////////////////////////////// */
/** @brief   Constructor for the Vpmg class (allocates new memory)
 *  @author  Nathan Baker
 *  @ingroup Vpmg
 *  @returns Pointer to newly allocated Vpmg object 
 */
VEXTERNC Vpmg* Vpmg_ctor(
        Vpmgp *parms,  /** PMG parameter object */
        Vpbe *pbe,  /** PBE-specific variables */
        int focusFlag,  /** 1 for focusing, 0 otherwise */
        Vpmg *pmgOLD,  /** Old Vpmg object to use for boundary conditions */
        MGparm *mgparm,  /** MGparm parameter object for boundary conditions */
        PBEparm_calcEnergy energyFlag  /** What types of energies to calculate */
        );

/** 
 * @brief  FORTRAN stub constructor for the Vpmg class (uses
 *         previously-allocated memory)
 *  @author  Nathan Baker
 *  @ingroup Vpmg
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_ctor2(
        Vpmg *thee,  /** Memory location for object */
        Vpmgp *parms,  /** PMG parameter object */
        Vpbe *pbe,  /** PBE-specific variables */
        int focusFlag,  /** 1 for focusing, 0 otherwise */
        Vpmg *pmgOLD,  /** Old Vpmg object to use for boundary conditions (can
                         be VNULL if focusFlag = 0) */
        MGparm *mgparm,  /** MGparm parameter object for boundary 
                          * conditions (can be VNULL if focusFlag = 0) */
        PBEparm_calcEnergy energyFlag  /** What types of energies to 
                                        * calculate (ignored if focusFlag
                                        * = 0) */
        );

/** @brief   Object destructor
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 */
VEXTERNC void Vpmg_dtor(
        Vpmg **thee  /** Pointer to memory location of object to be 
                      * destroyed */
        );

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 */
VEXTERNC void Vpmg_dtor2(
        Vpmg *thee  /** Pointer to object to be destroyed */
        );

/** @brief  Fill the coefficient arrays prior to solving the equation
 *  @ingroup  Vpmg
 *  @author  Nathan Baker
 *  @returns  1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_fillco(
        Vpmg *thee,  /** Vpmg object */ 
        Vsurf_Meth surfMeth,  /** Surface discretization method */
        double splineWin,  /** Spline window (in A) for surfMeth = 
                            * VSM_SPLINE */
        Vchrg_Meth chargeMeth,  /** Charge discretization method */ 
        int useDielXMap,  /** Boolean to use dielectric map argument */
        Vgrid *dielXMap,  /** External dielectric map */
        int useDielYMap,  /** Boolean to use dielectric map argument */
        Vgrid *dielYMap,  /** External dielectric map */
        int useDielZMap,  /** Boolean to use dielectric map argument */
        Vgrid *dielZMap,  /** External dielectric map */
        int useKappaMap,  /** Boolean to use kappa map argument */
        Vgrid *kappaMap,  /** External kappa map */
        int useChargeMap,  /** Boolean to use charge map argument */
        Vgrid *chargeMap  /** External charge map */
        );

/** @brief   Solve the PBE using PMG
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @returns  1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_solve(
        Vpmg *thee  /** Vpmg object */
        );

/** @brief   Solve Poisson's equation with a homogeneous Laplacian operator
 *           using the solvent dielectric constant.  This solution is
 *           performed by a sine wave decomposition.
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @returns  1 if successful, 0 otherwise
 *  @note    This function is really only for testing purposes as the
 *           PMG multigrid solver can solve the homogeneous system much more
 *           quickly.  Perhaps we should implement an FFT version at some
 *           point...
 */
VEXTERNC int Vpmg_solveLaplace(
        Vpmg *thee  /** Vpmg object */
        );

/** @brief   Get the total electrostatic energy.
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The electrostatic energy in units of k_B T.
 */
VEXTERNC double Vpmg_energy(
        Vpmg *thee,  /** Vpmg object */
        int extFlag  /** If this was a focused calculation, include (1 -- for
                      * serial calculations) or ignore (0 -- for parallel
                      * calculations) energy contributions from outside the
                      * focusing domain */
        );

/** @brief   Get the "fixed charge" contribution to the electrostatic energy 
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the fixed charges
 *           with the potential: \f[ G = \sum_i q_i u(r_i) \f]
 *           and return the result in units of k_B T.  Clearly, no
 *           self-interaction terms are removed.  A factor a 1/2 has to be
 *           included to convert this to a real energy.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The fixed charge electrostatic energy in units of k_B T.
 */
VEXTERNC double Vpmg_qfEnergy(
        Vpmg *thee,  /** Vpmg object */
        int extFlag  /** If this was a focused calculation, include (1 -- for
                      * serial calculations) or ignore (0 -- for parallel
                      * calculations) energy contributions from outside the
                      * focusing domain */
        );

/** @brief   Get the per-atom "fixed charge" contribution to the electrostatic
 *           energy
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the fixed charges
 *           with the potential: \f[ G = q u(r), \f] where q$ is the
 *           charge and r is the location of the atom of interest.  The
 *           result is returned in units of k_B T.  Clearly, no
 *           self-interaction terms are removed.  A factor a 1/2 has to be
 *           included to convert this to a real energy.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The fixed charge electrostatic energy in units of k_B T.
 */
VEXTERNC double Vpmg_qfAtomEnergy(
        Vpmg *thee,  /** The Vpmg object */
        Vatom *atom  /** The atom for energy calculations */
        );

/** @brief Get the "mobile charge" contribution to the electrostatic energy.
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the mobile charges 
 *           with the potential: 
 *              \f[ G = \frac{1}{4 I_s} \sum_i c_i q_i^2 \int
 *              \kappa^2(x) e^{-q_i u(x)} dx \f]
 *           for the NPBE and
 *              \f[ G = \frac{1}{2} \int \overline{\kappa}^2(x) u^2(x) dx \f]
 *           for the LPBE.  Here i denotes the counterion species, 
 *           I_s is the bulk ionic strength, kappa^2(x)
 *           is the modified Debye-Huckel parameter, c_i is the 
 *           concentration of species i, q_i is the charge of
 *           species i, and u(x) is the dimensionless electrostatic
 *           potential.  The energy is scaled to units of k_b T.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The mobile charge electrostatic energy in units of k_B T.
 */
VEXTERNC double Vpmg_qmEnergy(
        Vpmg *thee,  /** Vpmg object */
        int extFlag  /** If this was a focused calculation, include (1 -- for
                      * serial calculations) or ignore (0 -- for parallel
                      * calculations) energy contributions from outside the
                      * focusing domain */
        );


/** @brief Get the "polarization" contribution to the electrostatic energy.
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the mobile charges 
 *           with the potential: 
 *              \f[ G = \frac{1}{2} \int \epsilon (\nabla u)^2 dx \f]
 *           where epsilon is the dielectric parameter and u(x) is
 *           the dimensionless electrostatic potential.  The energy is scaled
 *           to units of k_b T.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The polarization electrostatic energy in units of k_B T.
 */
VEXTERNC double Vpmg_dielEnergy(
        Vpmg *thee,  /** Vpmg object */
        int extFlag  /** If this was a focused calculation, include (1 -- for
                      * serial calculations) or ignore (0 -- for parallel
                      * calculations) energy contributions from outside the
                      * focusing domain */
        );


/** @brief Get the integral of the gradient of the dielectric function
 *
 *           Using the dielectric map at the finest mesh level, calculate the
 *           integral of the norm of the dielectric function gradient
 *           routines of Im et al (see Vpmg_dbForce for reference):
 *              \f[ \int \| \nabla \epsilon \| dx \f]
 *           where epsilon is the dielectric parameter.
 *           The integral is returned in units of A^2.
 * 
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The integral in units of A^2.
 */
VEXTERNC double Vpmg_dielGradNorm(
        Vpmg *thee  /** Vpmg object */
        );

/** @brief    Calculate the total force on the specified atom in units of
 *            k_B T/AA
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing.
 * @returns  1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_force(
        Vpmg *thee,  /** Vpmg object */
        double *force, /** 3*sizeof(double) space to hold the force in units
                         of k_B T/AA */
        int atomID,  /** Valist ID of desired atom */
        Vsurf_Meth srfm,  /** Surface discretization method */ 
        Vchrg_Meth chgm  /** Charge discretization method */
        );

/** @brief    Calculate the "charge-field" force on the specified atom in units
 *           of k_B T/AA
 * @ingroup  Vpmg
 * @author   Nathan Baker
 * @note     \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field 
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing. 
 * @returns  1 if sucessful, 0 otherwise
 */
VEXTERNC int Vpmg_qfForce(
        Vpmg *thee,  /** Vpmg object */
        double *force, /** 3*sizeof(double) space to hold the force in units
                         of k_B T/A */
        int atomID,  /** Valist ID of desired atom */
        Vchrg_Meth chgm  /** Charge discretization method */
        );

/** @brief   Calculate the dielectric boundary forces on the
 *           specified atom in units of k_B T/AA
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field 
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing. 
 * @returns  1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_dbForce(
        Vpmg *thee,  /** Vpmg object */
        double *dbForce, /** 3*sizeof(double) space to hold the dielectric
                           boundary force in units of k_B T/AA */
        int atomID,  /** Valist ID of desired atom */
        Vsurf_Meth srfm  /** Surface discretization method */ 
        );

/** @brief   Calculate the osmotic pressure on the specified atom in units of
 *           k_B T/AA
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field 
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing. 
 *  @returns  1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_ibForce(
        Vpmg *thee,  /** Vpmg object */
        double *force, /** 3*sizeof(double) space to hold the 
                           boundary force in units of k_B T/AA */
        int atomID,  /** Valist ID of desired atom */
        Vsurf_Meth srfm  /** Surface discretization method */ 
        );

/** @brief   Set partition information which restricts the calculation of
 *           observables to a (rectangular) subset of the problem domain
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 */
VEXTERNC void Vpmg_setPart(
        Vpmg *thee,  /** Vpmg object */
        double lowerCorner[3],  /** Partition lower corner */ 
        double upperCorner[3],  /** Partition upper corner */
        int bflags[6]  /** Booleans indicating whether a particular processor
                         is on the boundary with another partition.  0 if the
                         face is not bounded (next to) another partition, and
                         1 otherwise. */
        );

/** @brief  Remove partition restrictions
 *  @ingroup  Vpmg
 *  @author  Nathan Baker
 */
VEXTERNC void Vpmg_unsetPart(
        Vpmg *thee  /** Vpmg object */
        );

/** @brief  Fill the specified array with accessibility values 
 *  @ingroup  Vpmg
 *  @author  Nathan Baker
 *  @returns  1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_fillArray(
        Vpmg *thee,  /** Vpmg object */
        double *vec,  /** A nx*ny*nz*sizeof(double) array to contain the
                        values to be written */
        Vdata_Type type,  /** What to write */ 
        double parm,  /** Parameter for data type definition (if needed) */
        Vhal_PBEType pbetype  /** Parameter for PBE type (if needed) */
        );

/** @brief   Computes the field at an atomic center using a stencil based
 *           on the first derivative of a 5th order B-spline
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VPUBLIC void Vpmg_fieldSpline4(
             Vpmg *thee,     /** Vpmg object */
             int atomID,     /** Atom index */
             double field[3] /** The (returned) electric field */
             );

/** @brief   Computes the permanent multipole electrostatic hydration 
 *           energy (the polarization component of the hydration energy 
 *           currently computed in TINKER).
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 *  @returns The permanent multipole electrostatic hydration energy
 */
VEXTERNC double Vpmg_qfPermanentMultipoleEnergy(
             Vpmg *thee,     /** Vpmg object */
             int atomID      /** Atom index */
             );

/** @brief   Computes the q-Phi Force for permanent multipoles based on
 *           5th order B-splines
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_qfPermanentMultipoleForce(
             Vpmg *thee,      /** Vpmg object */
             int atomID,      /** Atom index */
             double force[3], /** (returned) force */ 
             double torque[3] /** (returned) torque */
             );

/** @brief   Compute the ionic boundary force for permanent multipoles.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_ibPermanentMultipoleForce( 
             Vpmg *thee,      /** Vpmg object */
             int atomID,      /** Atom index */
             double force[3]  /** (returned) force */ 
             );

/** @brief   Compute the dielectric boundary force for permanent multipoles.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */                                          
VEXTERNC void Vpmg_dbPermanentMultipoleForce(
             Vpmg *thee,      /** Vpmg object */
             int atomID,      /** Atom index */
             double force[3]  /** (returned) force */ 
             );
                                             
/** @brief   q-Phi direct polarization force between permanent multipoles and
 *           induced dipoles, which are induced by the sum of the permanent
 *           intramolecular field and the permanent reaction field. 
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */         
VEXTERNC void Vpmg_qfDirectPolForce( 
             Vpmg *thee,      /** Vpmg object */
             Vgrid *perm,     /** Permanent multipole potential */
             Vgrid *induced,  /** Induced dipole potential */
             int atomID,      /** Atom index */
             double force[3], /** (returned) force */ 
             double torque[3] /** (returned) torque */
             );

/** @brief   q-Phi direct polarization force between permanent multipoles and
 *           non-local induced dipoles based on 5th Order B-Splines.
 *           Keep in mind that the "non-local" induced dipooles are just
 *           a mathematical quantity that result from differentiation of
 *           the AMOEBA polarization energy. 
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */        
VEXTERNC void Vpmg_qfNLDirectPolForce( 
             Vpmg *thee,      /** Vpmg object */
             Vgrid *perm,     /** Permanent multipole potential */
             Vgrid *nlInduced,/** Non-local induced dipole potential */
             int atomID,      /** Atom index */
             double force[3], /** (returned) force */ 
             double torque[3] /** (returned) torque */
             );

/** @brief   Ionic boundary direct polarization force between permanent
 *           multipoles and induced dipoles, which are induced by the 
 *           sum of the permanent intramolecular field and the permanent 
 *           reaction field. 
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */  
VEXTERNC void Vpmg_ibDirectPolForce(
             Vpmg *thee,      /** Vpmg object */
             Vgrid *perm,     /** Permanent multipole potential */
             Vgrid *induced,  /** Induced dipole potential */
             int atomID,      /** Atom index */
             double force[3]  /** (returned) force */ 
             );

/** @brief   Ionic boundary direct polarization force between permanent
 *           multipoles and non-local induced dipoles based on 5th order 
 *           Keep in mind that the "non-local" induced dipooles are just
 *           a mathematical quantity that result from differentiation of
 *           the AMOEBA polarization energy. 
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */  
VEXTERNC void Vpmg_ibNLDirectPolForce( 
             Vpmg *thee,      /** Vpmg object */
             Vgrid *perm,     /** Permanent multipole potential */
             Vgrid *nlInduced,/** Induced dipole potential */
             int atomID,      /** Atom index */
             double force[3]  /** (returned) force */ 
             );

/** @brief   Dielectric boundary direct polarization force between permanent
 *           multipoles and induced dipoles, which are induced by the 
 *           sum of the permanent intramolecular field and the permanent 
 *           reaction field. 
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */  
VEXTERNC void Vpmg_dbDirectPolForce(
             Vpmg *thee,      /** Vpmg object */
             Vgrid *perm,     /** Permanent multipole potential */
             Vgrid *induced,  /** Induced dipole potential */
             int atomID,      /** Atom index */
             double force[3]  /** (returned) force */ 
             );

/** @brief   Dielectric bounday direct polarization force between
 *           permanent multipoles and non-local induced dipoles.
 *           Keep in mind that the "non-local" induced dipooles are just
 *           a mathematical quantity that result from differentiation of
 *           the AMOEBA polarization energy. 
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */  
VEXTERNC void Vpmg_dbNLDirectPolForce(
             Vpmg *thee,      /** Vpmg object */
             Vgrid *perm,     /** Permanent multipole potential */
             Vgrid *nlInduced,/** Non-local induced dipole potential */
             int atomID,      /** Atom index */
             double force[3]  /** (returned) force */ 
             );

/** @brief   Mutual polarization force for induced dipoles based on 5th
 *           order B-Splines. This force arises due to self-consistent
 *           convergence of the solute induced dipoles and reaction field.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */  
VEXTERNC void Vpmg_qfMutualPolForce(
             Vpmg *thee,      /** Vpmg object */
             Vgrid *induced,  /** Induced dipole potential */
             Vgrid *nlInduced,/** Non-local induced dipole potential */
             int atomID,      /** Atom index */
             double force[3]  /** (returned) force */ 
             );

/** @brief   Ionic boundary mutual polarization force for induced dipoles 
 *           based on 5th order B-Splines. This force arises due to 
 *           self-consistent convergence of the solute induced dipoles 
 *           and reaction field.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */ 
VEXTERNC void Vpmg_ibMutualPolForce(
             Vpmg *thee,      /** Vpmg object */
             Vgrid *induced,  /** Induced dipole potential */
             Vgrid *nlInduced,/** Non-local induced dipole potential */
             int atomID,      /** Atom index */
             double force[3]  /** (returned) force */ 
             );

/** @brief   Dielectric boundary mutual polarization force for induced dipoles 
 *           based on 5th order B-Splines. This force arises due to 
 *           self-consistent convergence of the solute induced dipoles 
 *           and reaction field.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */ 
VEXTERNC void Vpmg_dbMutualPolForce(
             Vpmg *thee,      /** Vpmg object */
             Vgrid *induced,  /** Induced dipole potential */
             Vgrid *nlInduced,/** Non-local induced dipole potential */
             int atomID,      /** Atom index */
             double force[3]  /** (returned) force */ 
             );

/** @brief   Print out a column-compressed sparse matrix in Harwell-Boeing
 *           format.  
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @bug  Can this path variable be replaced with a Vio socket?
 */
VEXTERNC void Vpmg_printColComp(
        Vpmg *thee,  /**  Vpmg object */
        char path[72],  /** The file to which the matrix is to be written */
        char title[72],  /** The title of the matrix */
        char mxtype[3],   /** The type of REAL-valued matrix, a 3-character
                            string of the form "R_A" where the '_' can be one
                            of:  
                            \li S:  symmetric matrix
                            \li U:  unsymmetric matrix
                            \li H:  Hermitian matrix
                            \li Z:  skew-symmetric matrix
                            \li R:  rectangular matrix */
        int flag  /** The operator to compress:
                    \li 0:  Poisson operator
                    \li 1:  Linearization of the full Poisson-Boltzmann
                            operator around the current solution */
        );

#endif    /* ifndef _VPMG_H_ */

