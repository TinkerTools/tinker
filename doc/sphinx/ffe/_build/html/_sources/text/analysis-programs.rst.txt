Analysis & Utility Programs & Scripts
=====================================

This section of the manual contains a brief description of each of the Tinker structure manipulation, geometric calculation and auxiliary programs. A detailed example showing how to run each program is included in a later section. The programs listed below are all part of the main, supported distribution. Additional source code for various unsupported programs can be found in the /other directory of the Tinker distribution.

**ARCHIVE**

A program for concatenating Tinker cycle files into a single archive file; useful for storing the intermediate results of minimizations, dynamics trajectories, and so on. The program can also extract individual cycle files from a Tinker archive.

**CORRELATE**

A program to compute time correlation functions from collections of Tinker cycle files. Its use requires a user supplied function property that computes the value of the property for which a time correlation is desired for two input structures. A sample routine is supplied that computes either a velocity autocorrelation function or an rms structural superposition as a function of time. The main body of the program organizes the overall computation in an efficient manner and outputs the final time correlation function.

**CRYSTAL**

A program for the manipulation of crystal structures including interconversion of fractional and Cartesian coordinates, generation of the unit cell from an asymmetric unit, and building of a crystalline block of specified size via replication of a single unit cell. The present version can handle about 25 of the most common space groups, others can easily be added as needed by modification of the routine symmetry.

**DIFFUSE**

A program to compute the self-diffusion constant for a homogeneous liquid via the Einstein equation. A previously saved dynamics trajectory is read in and ``unfolded'' to reverse translation of molecules due to use of periodic boundary conditions. The average motion over all molecules is then used to compute the self-diffusion constant. While the current program assumes a homogeneous system, it should be easy to modify the code to handle diffusion of individual molecules or other desired effects.

**DISTGEOM**

A program to perform distance geometry calculations using variations on the classic metric matrix method. A user specified number of structures consistent with keyfile input distance and dihedral restraints is generated. Bond length and angle restraints are derived from the input structure. Trial distances between the triangle smoothed lower and upper bounds can be chosen via any of several metrization methods, including a very effective partial random pairwise scheme. The correct radius of gyration of the structure is automatically maintained by choosing trial distances from Gaussian distributions of appropriate mean and width. The initial embedded structures can be further refined against a geometric restraint-only potential using either a sequential minimization protocol or simulated annealing.

**DOCUMENT**

The DOCUMENT program is provided as a minimal listing and documentation tool. It operates on the Tinker source code, either individual files or the complete source listing produced by the command script listing.make, to generate lists of routines, common blocks or valid keywords. In addition, the program has the ability to output a formatted parameter listing from the standard Tinker parameter files.

**INTEDIT**

A program to allow interactive inspection and alteration of the internal coordinate definitions and values of a Tinker structure. If the structure is altered, the user has the option to write out a new internal coordinates file upon exit.

**INTXYZ**

A program to convert a Tinker .int internal coordinates formatted file into a Tinker .xyz Cartesian coordinates formatted file.

**MOL2XYZ**

A program for converting a Tripos Sybyl MOL2 file into a Tinker XYZ Cartesian coordinate file. The current version of the program does not attempt to convert the Sybyl atoms types into the active Tinker force field types, i.e., all atoms types are simply set to zero.

**NUCLEIC**

A program for automated building of nucleic acid structures. Upon interactive input of a nucleotide sequence with optional phosphate backbone angles, the program builds internal and Cartesian coordinates. Standard bond lengths and angles are used. Both DNA and RNA sequences are supported as are A-, B- and Z-form structures. Double helixes of complementary sequence can be automatically constructed via a rigid docking of individual strands.

**PDBXYZ**

A program for converting a Brookhaven Protein Data Bank file (a PDB file) into a Tinker .xyz Cartesian coordinate file. If the PDB file contains only protein/peptide amino acid residues, then standard protein connectivity is assumed, and transferred to the .xyz file. For non-protein portions of the PDB file, atom connectivity is determined by the program based on interatomic distances. The program also has the ability to add or remove hydrogen atoms from a protein as required by the force field specified during the computation.

**POLARIZE**

A program for computing molecular polarizability from an atom-based distributed model of polarizability. A damped interaction model due to Thole is optionally via keyfile settings. A Tinker .xyz file is required as input. The output consists of the overall polarizability tensor in the global coordinates and its eigenvalues.

**PRMEDIT**

A program for formatting and renumbering Tinker force field parameter files. When atom types or classes are added to a parameter file, this utility program has the ability to renumber all the atom records sequentially, and alter type and class numbers in all other parameter entries to maintain consistency.

**PROTEIN**

A program for automated building of peptide and protein structures. Upon interactive input of an amino acid sequence with optional phi/psi/omega/chi angles, D/L chirality, etc., the program builds internal and Cartesian coordinates. Standard bond lengths and angles are assumed for the peptide. The program will optionally convert the structure to a cyclic peptide, or add either or both N- and C-terminal capping groups. Atom type numbers are automatically assigned for the specified force field. The final coordinates and a sequence file are produced as the output.

**RADIAL**

A program to compute the pair radial distribution function between two atom types. The user supplies the two atom names for which the distribution function is to be computed, and the width of the distance bins for data analysis. A previously saved dynamics trajectory is read as input. The raw radial distribution and a spline smoothed version are then output from zero to a distance equal to half the minimum periodic box dimension. The atom names are matched to the atom name column of the Tinker .xyz file, independent of atom type.

**SPACEFILL**

A program to compute the volume and surface areas of molecules. Using a modified version of Connolly's original analytical description of the molecular surface, the program determines either the van der Waals, accessible or molecular (contact/reentrant) volume and surface area. Both surface area and volume are broken down into their geometric components, and surface area is decomposed into the convex contribution for each individual atom. The probe radius is input as a user option, and atomic radii can be set via the keyword file. If Tinker archive files are used as input, the program will compute the volume and surface area of each structure in the input file.

**SPECTRUM**

A program to compute a power spectrum from velocity autocorrelation data. As input, this program requires a velocity autocorrelation function as produced by the CORRELATE program. This data, along with a user input time step, are Fourier transformed to generate the spectral intensities over a wavelength range. The result is a power spectrum, and the positions of the bands are those predicted for an infrared or Raman spectrum. However, the data is not weighted by molecular dipole moment derivatives as would be required to produce correct IR intensities.

**SUPERPOSE**

A program to superimpose two molecular structures in 3-dimensions. A variety of options for input of the atom sets to be used during the superposition are presented interactively to the user. The superposition can be mass-weighted if desired, and the coordinates of the second structure superimposed on the first structure are optionally output. If Tinker archive files are used as input, the program will compute all pairwise superpositions between structures in the input files.

**XYZEDIT**

A program that performs and of a variety of manipulations on an input Tinker .xyz Cartesian coordinates formatted file. The present version of the program has the following interactively selectable options: (1) Offset the Numbers of the Current Atoms, (2) Deletion of Individual Specified Atoms, (3) Deletion of Specified Types of Atoms, (4) Deletion of Atoms outside Cutoff Range, (5) Insertion of Individual Specified Atoms, (6) Replace Old Atom Type with a New Type, (7) Assign Connectivities based on Distance, (8) Convert Units from Bohrs to Angstroms, (9) Invert thru Origin to give Mirror Image, (10) Translate Center of Mass to the Origin, (11) Translate a Specified Atom to the Origin, (12) Translate and Rotate to Inertial Frame, (13) Move to Specified Rigid Body Coordinates, (14) Create and Fill a Periodic Boundary Box, (15) Soak Current Molecule in Box of Solvent, (16) Append another XYZ file to Current One. In most cases, multiply options can be applied sequentially to an input file. At the end of the editing process, a new version of the original .xyz file is written as output.

**XYZINT**

A program for converting a Tinker .xyz Cartesian coordinate formatted file into a Tinker .int internal coordinates formatted file. This program can optionally use an existing internal coordinates file as a template for the connectivity information.

**XYZMOL2**

A program to convert a Tinker .xyz Cartesian coordinates file into a Tripos Sybyl MOL2 file. The conversion generates only the MOLECULE, ATOM, BOND and SUBSTRUCTURE record type in the MOL2 file. Generic Sybyl atom types are used in most cases; while these atom types may need to be altered in some cases, Sybyl is usually able to correctly display the resulting MOL2 file.

**XYZPDB**

A program for converting a Tinker .xyz Cartesian coordinate file into a Brookhaven Protein Data Bank file (a PDB file).
