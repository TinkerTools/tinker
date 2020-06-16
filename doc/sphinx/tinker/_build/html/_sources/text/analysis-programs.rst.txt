Analysis & Utility Programs & Scripts
=====================================

This section of the manual contains a brief description of each of the Tinker structure manipulation, geometric calculation and auxiliary programs. A detailed example showing how to run each program is included in a later section. The programs listed below are all part of the main, supported distribution. Additional source code for various unsupported programs can be found in the /other directory of the Tinker distribution.

**ARCHIVE**

ARCHIVE is a program for concatenating Tinker cycle files into a single archive file; useful for storing the intermediate results of minimizations, dynamics trajectories, and so on. The program can also extract individual cycle files from a Tinker archive.

**BAR**

The BAR program computes a free energy from sampling of adjacent "lambda" windows using the Bennett acceptance ratio (BAR) algorithm. Input consists of trajectories or configurations sampled from the adjacent windows, as well as keyfiles and parameters used to define the states for the simulations. In a first phase, the BAR program computes the energies of all structures from both simulations under the control of both sets of potential energy parameters, i.e., four sets of numbers which are written to an intermediate .bar file. In its second phase, BAR reads a .bar file and uses the free energy perturbation (FEP) and Bennett acceptance ratio formula to compute the free energy, enthalpy and entropy between the two states.

**CORRELATE**

The CORRELATE program to compute time correlation functions from collections of Tinker cycle files. Its use requires a user supplied function property that computes the value of the property for which a time correlation is desired for two input structures. A sample routine is supplied that computes either a velocity autocorrelation function or an rms structural superposition as a function of time. The main body of the program organizes the overall computation in an efficient manner and outputs the final time correlation function.

**CRYSTAL**

CRYSTAL is a program for the manipulation of crystal structures including interconversion of fractional and Cartesian coordinates, generation of the unit cell from an asymmetric unit, and building of a crystalline block of specified size via replication of a single unit cell. The present version can handle about 25 of the most common space groups, others can easily be added as needed by modification of the routine symmetry.

**DIFFUSE**

DIFFUSE computes the self-diffusion constant for a homogeneous liquid via the Einstein equation. A previously saved dynamics trajectory is read in and "unfolded" to reverse translation of molecules due to use of periodic boundary conditions. The average motion over all molecules is then used to compute the self-diffusion constant. While the current program assumes a homogeneous system, it should be easy to modify the code to handle diffusion of individual molecules or other desired effects.

**DISTGEOM**

The DISTGEOM program performs distance geometry calculations using variations on the classic metric matrix method. A user specified number of structures consistent with keyfile input distance and dihedral restraints is generated. Bond length and angle restraints are derived from the input structure. Trial distances between the triangle smoothed lower and upper bounds can be chosen via any of several metrization methods, including a very effective partial random pairwise scheme. The correct radius of gyration of the structure is automatically maintained by choosing trial distances from Gaussian distributions of appropriate mean and width. The initial embedded structures can be further refined against a geometric restraint-only potential using either a sequential minimization protocol or simulated annealing.

**DOCUMENT**

The DOCUMENT program is provided as a minimal listing and documentation tool. It operates on the Tinker source code, either individual files or the complete source listing produced by the command script listing.make, to generate lists of routines, common blocks or valid keywords. In addition, the program has the ability to output a formatted parameter listing from the standard Tinker parameter files.

**FREEFIX**

FREEFIX is a small utility to compute the analytical enthalpy, entropy and free energy associated with the release of a flat-bottomed harmonic distance restraint between two sites within a simulation system.

**INTEDIT**

INTEDIT allows interactive inspection and alteration of the internal coordinate definitions and values of a Tinker structure. If the structure is altered, the user has the option to write out a new internal coordinates file upon exit.

**INTXYZ**

The INTXYZ program to convert a Tinker .int internal coordinates formatted file into a Tinker .xyz Cartesian coordinates formatted file.

**MOLXYZ**

MOLXYZ is a program for converting a MDL (Molecular Design Limited) MOL file into a Tinker XYZ Cartesian coordinate file. The current version of the program converts the MDL atoms types into Tinker "tiny force field" atom types based on atomic number and connectivity (i.e., a tetravalent carbon is type 64).

**MOL2XYZ**

The MOL2XYZ program converts a Tripos Sybyl MOL2 file into a Tinker XYZ Cartesian coordinate file. The current version of the program converts the Sybyl MOL2 atoms types into Tinker "tiny force field" atom types based on atomic number and connectivity (i.e., a tetravalent carbon is type 64).

**NUCLEIC**

The NUCLEIC program automates building of nucleic acid structures. Upon interactive input of a nucleotide sequence with optional phosphate backbone angles, the program builds internal and Cartesian coordinates. Standard bond lengths and angles are used. Both DNA and RNA sequences are supported as are A-, B- and Z-form structures. Double helixes of complementary sequence can be automatically constructed via a rigid docking of individual strands.

**PDBXYZ**

PDBXYZ is a program for converting a Brookhaven Protein Data Bank file (a PDB file) into a Tinker .xyz Cartesian coordinate file. If the PDB file contains only protein/peptide amino acid residues, then standard protein connectivity is assumed, and transferred to the .xyz file. For non-protein portions of the PDB file, atom connectivity is determined by the program based on interatomic distances. The program also has the ability to add or remove hydrogen atoms from a protein as required by the force field specified during the computation.

**POLARIZE**

POLARIZE is a simple program for computing molecular polarizability from an atom-based distributed model of polarizability. POLARIZE implements whichever damped interaction model is specified via keyfile and parameter settings. A Tinker .xyz file is required as input. The output consists of the overall polarizability tensor in the global coordinates and its eigenvalues.

**POLEDIT**

POLEDIT is a program for manipulating and processing polarizable atomic multipole models. Its primary use is to read a distributed multipole analysis (DMA) from output of the GDMA or Psi4 quantum chemistry programs. The program defines local coordinate frames, sets atomic polarizabilities, removes molecular mechanics polarization from the quantum DMA, averages over symmetrical atoms and outputs parameters in Tinker format. There are additional invocation options to only change local coordinate frame definitions or remove intramolecular polarization from an existing multipole model.

**POTENTIAL**

The POTENTIAL program performs electrostatic potential comparisons and fitting. POTENTIAL can compare two different force field electrostatic models via computing the RMS between the electrostatic potentials on a grid of points outside the molecular envelope. An electrostatic potential grid can also be generated from quantum chemistry output, and compare against a force field model. Finally, a flexible fitting of a force field model to an existing potential grid is available. The program can also take as model input a set of different molecules containing common types, and multiple conformations of a single molecule.

**PRMEDIT**

PRMEDIT is a program for formatting and renumbering Tinker force field parameter files. When atom types or classes are added to a parameter file, this utility program has the ability to renumber all the atom records sequentially, and alter type and class numbers in all other parameter entries to maintain consistency.

**PROTEIN**

The PROTEIN program automates building of peptide and protein structures. Upon interactive input of an amino acid sequence with optional phi/psi/omega/chi angles, D/L chirality, etc., the program builds internal and Cartesian coordinates. Standard bond lengths and angles are assumed for the peptide. The program will optionally convert the structure to a cyclic peptide, or add either or both N- and C-terminal capping groups. Atom type numbers are automatically assigned for the specified force field. The final coordinates and a sequence file are produced as the output.

**RADIAL**

The RADIAL program finds the pair radial distribution function between two atom types. The user supplies the two atom names for which the distribution function is to be computed, and the width of the distance bins for data analysis. A previously saved dynamics trajectory is read as input. The raw radial distribution and a spline smoothed version are then output from zero to a distance equal to half the minimum periodic box dimension. The atom names are matched to the atom name column of the Tinker .xyz file, independent of atom type.

**SPACEFILL**

The SPACEFILL program computes the volume and surface areas of molecules. Using a modified version of Connolly's original analytical description of the molecular surface, the program determines either the van der Waals, accessible or molecular (contact/reentrant) volume and surface area. Both surface area and volume are broken down into their geometric components, and surface area is decomposed into the convex contribution for each individual atom. The probe radius is input as a user option, and atomic radii can be set via the keyword file. If Tinker archive files are used as input, the program will compute the volume and surface area of each structure in the input file.

**SPECTRUM**

SPECTRUM is a program to compute a power spectrum from velocity autocorrelation data. As input, this program requires a velocity autocorrelation function as produced by the CORRELATE program. This data, along with a user input time step, are Fourier transformed to generate the spectral intensities over a wavelength range. The result is a power spectrum, and the positions of the bands are those predicted for an infrared or Raman spectrum. However, the data is not weighted by molecular dipole moment derivatives as would be required to produce correct IR intensities.

**SUPERPOSE**

The SUPERPOSE program is used to superimpose two molecular structures in 3-dimensions. A variety of options for input of the atom sets to be used during the superposition are presented interactively to the user. The superposition can be mass-weighted if desired, and the coordinates of the second structure superimposed on the first structure are optionally output. If Tinker archive files are used as input, the program will compute all pairwise superpositions between structures in the input files.

**TORSFIT**

TORSFIT is a program for setting force field parameters for torsional terms by fitting 1-fold to 6-fold torsional amplitudes to the difference between a quantum chemistry rotational profile and a force field rotational profile without any torsional terms.

**VALENCE**

VALENCE is a program for setting force field parameters for local valence terms, either from quantum chemistry data or from embedded empirical rules. [This program is still under development.]

**XYZEDIT**

XYZEDIT is a program to perform a variety of manipulations on an input Tinker .xyz Cartesian coordinates formatted file. The present version of the program has the following interactively selectable options: (1) Offset the Numbers of the Current Atoms, (2) Deletion of Individual Specified Atoms, (3) Deletion of Specified Types of Atoms, (4) Deletion of Atoms outside Cutoff Range, (5) Insertion of Individual Specified Atoms, (6) Replace Old Atom Type with a New Type, (7) Assign Connectivities based on Distance, (8) Convert Units from Bohrs to Angstroms, (9) Invert thru Origin to give Mirror Image, (10) Translate Center of Mass to the Origin, (11) Translate a Specified Atom to the Origin, (12) Translate and Rotate to Inertial Frame, (13) Move to Specified Rigid Body Coordinates, (14) Create and Fill a Periodic Boundary Box, (15) Soak Current Molecule in Box of Solvent, (16) Append another XYZ file to Current One. In most cases, multiply options can be applied sequentially to an input file. At the end of the editing process, a new version of the original .xyz file is written as output.

**XYZINT**

XYZINT converts a Tinker .xyz Cartesian coordinate formatted file into a Tinker .int internal coordinates formatted file. This program can optionally use an existing internal coordinates file as a template for the connectivity information.

**XYZMOL2**

XYZMOL2 is a program to convert a Tinker .xyz Cartesian coordinates file into a Tripos Sybyl MOL2 file. The conversion generates only the MOLECULE, ATOM, BOND and SUBSTRUCTURE record type in the MOL2 file. Generic Sybyl atom types are used in most cases; while these atom types may need to be altered in some cases, Sybyl is usually able to correctly display the resulting MOL2 file.

**XYZPDB**

The XYZPDB program converts a Tinker .xyz Cartesian coordinate file into a Brookhaven Protein Data Bank file (a PDB file). A Tinker .seq file with the biopolymer sequence must be present if the output PDB file is to be formatted as a protein or nucleic acid with a defined sequence.
