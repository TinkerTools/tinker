Types of Input & Output Files
=============================

This section describes the basic file types used by the Tinker package. Let's say you wish to perform a calculation on a particular small organic molecule. Assume that the file name chosen for our input and output files is sample. Then all of the Tinker files will reside on the computer under the name sample.xxx where .xxx is any of the several extension types to be described below.

**SAMPLE.XYZ**

The .xyz file is the basic Tinker Cartesian coordinates file type. It contains a title line followed by one line for each atom in the structure. Each line contains: the sequential number within the structure, an atomic symbol or name, X-, Y-, and Z-coordinates, the force field atom type number of the atom, and a list of the atoms connected to the current atom. Except for programs whose basic operation is in torsional space, all Tinker calculations are done from some version of the .xyz format.

**SAMPLE.INT**

The .int file contains an internal coordinates representation of the molecular structure. It consists of a title line followed by one line for each atom in the structure. Each line contains: the sequential number within the structure, an atomic symbol or name, the force field atom type number of the atom, and internal coordinates in the usual Z-matrix format. For each atom the internal coordinates consist of a distance to some previously defined atom, and either two bond angles or a bond angle and a dihedral angle to previous atoms. The length, angle and dihedral definitions do not have to represent real bonded interactions. Following the last atom definition are two optional blank line separated sets of atom number pairs. The first list contains pairs of atoms that are covalently bonded, but whose bond length was not used as part of the atom definitions. These pairs are typically used to close ring structures. The second list contains ``bonds'' that are to be broken, i.e., pairs of atoms that are not covalently bonded, but which were used to define a distance in the atom definitions.

**SAMPLE.KEY**

The keyword parameter file always has the extension .key and is optionally present during Tinker calculations. It contains values for any of a wide variety of switches and parameters that are used to change the course of the computation from the default. The detailed contents of this file is explained in a latter section of this User's Guide. If a molecular system specific keyfile, in this case sample.key, is not present, the the Tinker program will look in the same directory for a generic file named Tinker.key.

**SAMPLE.DYN**

The .dyn file contains values needed to restart a molecular or stochastic dynamics computation. It stores the current position, current velocity and current and previous accelerations for each atom, as well as the size and shape of any periodic box or crystal unit cell. This information can be used to start a new dynamics run from the final state of a previous run. Upon startup, the dynamics programs always check for the presence of a .dyn file and make use of it whenever possible. The .dyn file is updated concurrent with the saving of a new dynamics trajectory snapshot.

**SAMPLE.END**

The .end file type provides a mechanism to gracefully stop a running Tinker calculation. At appropriate checkpoints during a calculation, Tinker will test for the presence of a sample.end file, and if found will terminate the calculation after updating the output. The .end file can be created at any time during a computation, and will be detected when the next checkpoint is reached. The file may be of zero size, and its contents are unimportant. In the current version of Tinker, the .end mechanism is only available within dynamics-based programs.

**SAMPLE.001, SAMPLE.002, ....**

Several types of computations produce files containing a three or more digit extension (.001 as shown; or .002, .137, .5678, etc.). These are referred to as cycle files, and are used to store various types of output structures. The cycle files from a given computation are identical in internal structure to either the .xyz or .int files described above. For example, the vibrational analysis program can save the tenth normal mode in sample.010. A molecular dynamics-based program might save its tenth 0.1 picosecond frame (or an energy minimizer its tenth partially minimized intermediate) in a file of the same name.

**SAMPLE.LOG**

The Force Field Explorer interface to Tinker saves results of all calculations launched from the GUI to a log file with the .log suffix. Any output that would normally be directed to the screen after starting a program from the command line is appended to this log file by Force Field Explorer.

**SAMPLE.ARC**

A Tinker archive file is simply a series of .xyz Cartesian coordinate files appended together one after another. This file can be used to condense the results from intermediate stages of an optimization, frames from a molecular dynamics trajectory, or set of normal mode vibrations into a single file for storage. Tinker archive files can be displayed as sequential frame "movies" by the Force Field Explorer modeling program.

**SAMPLE.PDB**

This file type contains coordinate information in the PDB format developed by the Brookhaven Protein Data Bank for deposition of model structures based on macromolecular X-ray diffraction and NMR data. Although Tinker itself does not use .pdb files directly for input/output, auxiliary programs are provided with the system for interconverting .pdb files with the .xyz format described above.

**SAMPLE.SEQ**

This file type contains the primary sequence of a biopolymer in the standard one-letter code with 50 residues per line. The .seq file for a biopolymer is generated automatically when a PDB file is converted to Tinker .xyz format or when using the PROTEIN or NUCLEIC programs to build a structure from sequence It is required for the reverse conversion of a Tinker file back to PDB format..

**SAMPLE.FRAC**

The fractional coordinates corresponding to the asymmetric unit of a crystal unit cell are stored in the .frac file. The internal format of this file is identical to the .xyz file; except that the coordinates are fractional instead of in Angstrom units.

**SAMPLE.MOL2**

File conversion to and from the Tripos Sybyl MOL2 file format is supported by Tinker. The utility programs XYZMOL2 and MOL2XYZ transform a Tinker XYZ file to MOL2 format, and the reverse.

**PARAMETER FILES (*.PRM)**

The potential energy parameter files distributed with the Tinker package all end in the extension .prm, although this is not required by the programs themselves. Each of these files contains a definition of the potential energy functional forms for that force field as well as values for individual energy parameters. For example, the mm3pro.prm file contains the energy parameters and definitions needed for a protein-specific version of the MM3 force field.
