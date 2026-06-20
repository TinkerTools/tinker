Test Cases & Examples
=====================

This section contains brief descriptions of the sample calculations found in the TEST subdirectory of the Tinker distribution. These test cases exercise several of the current Tinker programs and are intended to provide an overview of the capabilities of the package. In addition, a number of input files for various molecular systems are provided in the EXAMPLE subdirectory.

**ANION Test**

Estimates the hydration free energy difference for Cl- vs. Br- anion via a 10 picosecond simulation of a "hybrid" anion in a box of water, followed by free energy perturbation

**ARGON Test**

Performs an initial energy minimization on a periodic box containing 150 argon atoms, then performs 25 picoseconds of molecular dynamics simulation, on a box with 150 argon atoms

**CATION Test**

Computes the hydration free energy difference for Rb+ vs. Cs+ cation via a 2 picosecond simulation of each cation in a box of water, followed by a BAR free energy calculation

**CLUSTER Test**

Performs a set of 10 Gaussian density annealing (GDA) trials on a cluster of 13 argon atoms to find the global minimum energy structure

**CRAMBIN Test**

Generates a Tinker XYZ file from a PDB file, followed by single point energy computation and determination of the molecular volume and surface area

**CYCLOHEXANE Test**

Locates the transition state between chair and boat cyclohexane via two methods: the Muller-Brown saddle point method, and path sampling using the Elber algorithm; vibrational analysis of the results shows the same TS with one negative frequency

**DHFR Test**

Runs 10 steps of molecular dynamics on an equilibrated system of DHFR protein in water using the AMOEBA force field; note this is the so-called Joint Amber-CHARMM "JAC" benchmark containing 23558 total atoms

**DIALANINE Test**

Finds all the local minima of alanine dipeptide via a potential energy surface scan using torsional modes to jump between minima

**ENKEPHALIN Test**

Builds coordinates for Met-enkephalin from amino acid sequence and phi/psi angles, followed by truncated Newton energy minimization and determination of the lowest frequency normal mode

**ETHANOL Test**

Fits torsional parameter values for the ethanol C-C-O-H bond based on relative quantum mechanical energies from Gaussian for rotating the C-O bond

**FORMAMIDE Test**

Generates a unit cell from fractional coordinates, followed by full crystal energy minimization and determination of optimal carbonyl oxygen parameters via a fit to lattice energy and structure

**GPCR Test**

Finds the lowest-frequency bacteriorhodopsin normal mode using a sliding block iterative diagonalization; alter the gpcr.run script to save the file gpcr.001 if you want view of the mode; this example can require up to an hour to complete

**HELIX Test**

Performs rigid-body optimization of the packing of two ideal polyalanine helices using only van der Waals interactions

**ICE Test**

Performs a short MD simulation of the monoclinic ice V crystal form using the iAMOEBA water model, pairwise neighbor lists and PME electrostatics

**IFABP Test**

Generates three distance geometry structures for intestinal fatty acid binding protein from a set of NOE distance restraints and torsional restraints.

**LIQUID Test**

Prints the system setup and computes the force field energy components for three small liquid water boxes using the AMOEBA, AMOEBA+ and HIPPO force fields

**METHANOL Test**

Processes distributed multipole analysis (DMA) output to extract coordinates and permanent multipoles, set local frames and polarization groups, modify the intramolecular polarization, detect and average equivalent atomic sites

**NITROGEN Test**

Calculates the self-diffusion constant and the N-N radial distribution function for liquid nitrogen via analysis of a 50 picosecond MD trajectory

**POLYALA Test**

Generates an extended conformation of capped alanine octapeptide, then uses Monte Carlo Minimization with torsion moves to find the 3/10 helix global minimum

**PYRIDINE Test**

Converts a simple XYZ file for pyridine to Tinker XYZ format using the BASIC force field, and computes the molecular mechanics energy

**SALT Test**

Converts a sodium chloride assymetric unit to the corresponding unit cell, then minimizes the crystal starting from the diffraction structure using Ewald summation to model long-range electrostatics

**SCORPION Test**

Converts a multi-model RCSB PDBx/mmCIF for scorpion toxin (1CHL) to Tinker XYZ and evaluates energies, then converts the XYZ to legacy PDB format and back to Tinker XYZ, and finally evaluates energies of the new XYZ file to compare with the prior values

**TIP4P Test**

Minimizes the energy to a small RMS gradient and enforces holonomic constraints for a cubic box of 1000 TIP4P water molecules, then performs a short MD simulation in the NVT ensemble at 298K

**VASOPRESSIN Test**

Compares analytical and finite difference numerical gradients over Cartesian and internal coordinates for vasopressin using the AMOEBA force field model

**WATER Test**

Fits the electrostatic potential for the TIP3P, AMOEBA and HIPPO water models to a QM-derived potential at the MP2/aug-cc-pVTZ level on a grid of points outside the molecular surface
