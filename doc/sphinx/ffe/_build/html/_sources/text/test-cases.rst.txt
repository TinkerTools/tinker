Test Cases & Examples
=====================

This section contains brief descriptions of the sample calculations found in the EXAMPLE subdirectory of the Tinker distribution. These examples exercise several of the current Tinker programs and are intended to provide a flavor of the capabilities of the package.

**ANION Test**

Computes an estimation of the free energy of hydration of Cl- anion vs. Br- anion via a 2 picosecond simulation on a "hybrid" anion in a box of water followed by a free energy perturbation calculation.

**ARGON Test**

Performs an initial energy minimization on a periodic box containing 150 argon atoms followed by 6 picoseconds of a molecular dynamics using a modified Beeman integration algorithm and a Bersedsen thermostat.

**CLUSTER Test**

Performs a set of 10 Gaussian density annealing (GDA) trials on a cluster of 13 argon atoms in an attempt to locate the global minimum energy structure.

**CRAMBIN Test**

Generates a Tinker file from a PDB file, followed by a single point energy computation and determination of the molecular volume and surface area.

**CYCLOHEX Test**

First approximately locates the transition state between chair and boat cyclohexane, followed by subsequent refinement of the transition state and a final vibrational analysis to show that a single negative frequency is associated with the saddle point.

**DHFR Test**

Performs 10 steps of molecular dynamics on a pre-equilibrated system of DHFR protein in a box or water using the AMOEBA force field. Note this test case is the so-called Joint Amber-CHARMM "JAC" benchmark containing 23558 total atoms.

**DIALANINE Test**

Finds all the local minima of alanine dipeptide via a potential energy surface scan using torsional modes to jump between the minima.

**ENKEPHALIN Test**

Produces coordinates from the met-enkephalin amino acid sequence and phi/psi angles, followed by truncated Newton energy minimization and determination of the lowest frequency normal mode.

**ETHANOL Test**

Performs fitting of torsional parameter values for the ethanol C-C-O-H bond based on relative quantum mechanical (G09) energies for rotating the C-O bond.

**FORMAMIDE Test**

Generates a unit cell from fractional coordinates, followed by full crystal energy minimization and determination of optimal carbonyl oxygen energy parameters from a fit to lattice energy and structure.

**GPCR Test**

Finds the lowest-frequency normal mode of bacteriorhodopsin using vibrational analysis via a sliding block iterative matrix diagonalization. Alter the gpcr.run script to save the file gpcr.001 for later viewing of the mode.

**HELIX Test**

Performs a rigid-body optimization of the packing of two idealized polyalanine helices using only van der Waals interactions.

**ICE Test**

Performs a short MD simulation of the monoclinic ice V crystal form using the iAMOEBA water model, pairwise neighbor lists and PME electrostatics.

**IFABP Test**

Generates three distance geometry structures for intestinal fatty acid binding protein from a set of NOE distance restraints and torsional restraints.

**METHANOL Test**

Processes distributed multipole analysis (DMA) output to extract coordinates and permanent multipoles, set local frames and polarization groups, remove intramolecular polarization, detect and average equivalent atomic sites.

**NITROGEN Test**

Calculates the self-diffusion constant and the N-N radial distribution function for liquid nitrogen via analysis of a 50ps MD trajectory.

**SALT Test**

Converts a sodium chloride assymetric unit to the corresponding unit cell, then runs a crystal minimization starting from the initial diffraction structure using Ewald summation to model the long-range electrostatic interactions.

**TETRAALA Test**

Generates capped alanine tetrapeptide in an extended conformation, then use Monte Carlo Minimization with random torsional moves to find the global minimum energy structure.

**WATER Test**

Fits the electrostatic potential around an AMOEBA water molecule to the QM-derived potential (MP2/aug-cc-pVTZ) on a grid of points outside the molecular surface.
