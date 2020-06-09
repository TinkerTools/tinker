Routines & Functions
====================

The distribution version of the Tinker package contains over 1000 separate programs, subroutines and functions. This section contains a brief description of the purpose of most of these code units. Further information can be found in the comments located at the top of each source code file.

**ACTIVE Subroutine**

"active" sets the list of atoms that are used during
each potential energy function calculation

**ADDBASE Subroutine**

"addbase" builds the Cartesian coordinates for a single nucleic
acid base; coordinates are read from the Protein Data Bank file
or found from internal coordinates, then atom types are assigned
and connectivity data generated

**ADDBOND Subroutine**

"addbond" adds entries to the attached atoms list in
order to generate a direct connection between two atoms

**ADDIONS Subroutine**

"addions" takes a currently defined solvated system and
places ions, with removal of solvent molecules

**ADDSIDE Subroutine**

"addside" builds the Cartesian coordinates for a single amino
acid side chain; coordinates are read from the Protein Data
Bank file or found from internal coordinates, then atom types
are assigned and connectivity data generated

**ADJACENT Function**

"adjacent" finds an atom connected to atom "i1" other than
atom "i2"; if no such atom exists, then the closest atom
in space is returned

**ADJUST Subroutine**

"adjust" modifies site bounds on the PME grid and returns
an offset into the B-spline coefficient arrays

**ALCHEMY Program**

"alchemy" computes the free energy difference corresponding
to a small perturbation by Boltzmann weighting the potential
energy difference over a number of sample states; current
version (incorrectly) considers the charge energy to be
intermolecular in finding the perturbation energies

**ALTELEC Subroutine**

"altelec" constructs mutated electrostatic parameters based
on the lambda mutation parameter "elambda"

**ALTERCHG Subroutine**

"alterchg" calculates the change in atomic partial charge or
monopole values due to bond and angle charge flux coupling

**ALTERPOL Subroutine**

"alterpol" finds an output set of atomic multipole parameters
which when used with an intergroup polarization model will
give the same electrostatic potential around the molecule as
the input set of multipole parameters with all atoms in one
polarization group

**ALTTORS Subroutine**

"alttors" constructs mutated torsional parameters based
on the lambda mutation parameter "tlambda"

**AMBERYZE Subroutine**

"amberyze" prints the force field parameters in a format needed
by the Amber setup protocol for using AMOEBA within Amber

**ANALYSIS Subroutine**

"analysis" calls the series of routines needed to calculate
the potential energy and perform energy partitioning analysis
in terms of type of interaction or atom number

**ANALYZE Program**

"analyze" computes and displays the total potential energy;
options are provided to display system and force field info,
partition the energy by atom or by potential function type,
show force field parameters by atom; output the large energy
interactions and find electrostatic and inertial properties

**ANGCHG Subroutine**

"angchg" computes modifications to atomic partial charges or
monopoles due to angle bending using a charge flux formulation

**ANGGUESS Function**

"angguess" sets approximate angle bend force constants based
on atom type and connected atoms

**ANGLES Subroutine**

"angles" finds the total number of bond angles and stores
the atom numbers of the atoms defining each angle; for
each angle to a trivalent central atom, the third bonded
atom is stored for use in out-of-plane bending

**ANNEAL Program**

"anneal" performs a simulated annealing protocol by means of
variable temperature molecular dynamics using either linear,
exponential or sigmoidal cooling schedules

**ANORM Function**

"anorm" finds the norm (length) of a vector; used as a
service routine by the Connolly surface area and volume
computation

**APBSEMPOLE Subroutine**

**APBSFINAL Subroutine**

**APBSINDUCE Subroutine**

**APBSINITIAL Subroutine**

**APBSNLINDUCE Subroutine**

**ARCHIVE Program**

"archive" is a utility program for coordinate files which
concatenates multiple coordinate sets into a new archive or
performs any of several manipulations on an existing archive

**ASET Subroutine**

"aset" computes by recursion the A functions used in the
evaluation of Slater-type (STO) overlap integrals

**ATOMYZE Subroutine**

"atomyze" prints the potential energy components broken
down by atom and to a choice of precision

**ATTACH Subroutine**

"attach" generates lists of 1-3, 1-4 and 1-5 connectivities
starting from the previously determined list of attached
atoms (ie, 1-2 connectivity)

**AUXINIT Subroutine**

"auxinit" initializes auxiliary variables and settings for
inertial extended Lagrangian induced dipole prediction

**AVGPOLE Subroutine**

"avgpole" condenses the number of multipole atom types based
upon atoms with equivalent attachments and additional user
specified sets of equivalent atoms

**BAOAB Subroutine**

"baoab" implements a constrained stochastic dynamics time
step using the geodesic BAOAB scheme

**BAR Program**

"bar" computes the free energy, enthalpy and entropy difference
between two states via Zwanzig free energy perturbation (FEP)
and Bennett acceptance ratio (BAR) methods

**BARCALC Subroutine**

**BASEFILE Subroutine**

"basefile" extracts from an input filename the portion
consisting of any directory name and the base filename;
also reads any keyfile and sets information level values

**BCUCOF Subroutine**

"bcucof" determines the coefficient matrix needed for bicubic
interpolation of a function, gradients and cross derivatives

**BCUINT Subroutine**

"bcuint" performs a bicubic interpolation of the function
value on a 2D spline grid

**BCUINT1 Subroutine**

"bcuint1" performs a bicubic interpolation of the function
value and gradient along the directions of a 2D spline grid

**BCUINT2 Subroutine**

"bcuint2" performs a bicubic interpolation of the function value,
gradient and Hessian along the directions of a 2D spline grid

**BEEMAN Subroutine**

"beeman" performs a single molecular dynamics time step
via the Beeman multistep recursion formula; uses original
coefficients or Bernie Brooks' "Better Beeman" values

**BETACF Function**

"betacf" computes a rapidly convergent continued fraction needed
by routine "betai" to evaluate the cumulative Beta distribution

**BETAI Function**

"betai" evaluates the cumulative Beta distribution function
as the probability that a random variable from a distribution
with Beta parameters "a" and "b" will be less than "x"

**BIGBLOCK Subroutine**

"bigblock" replicates the coordinates of a single unit cell
to give a larger unit cell as a block of repeated units

**BIOSORT Subroutine**

"biosort" renumbers and formats biotype parameters used to
convert biomolecular structure into force field atom types

**BITORS Subroutine**

"bitors" finds the total number of bitorsions as pairs
of adjacent torsional angles, and the numbers of the five
atoms defining each bitorsion

**BMAX Function**

"bmax" computes the maximum order of the B functions needed
for evaluation of Slater-type (STO) overlap integrals

**BNDCHG Subroutine**

"bndchg" computes modifications to atomic partial charges or
monopoles due to bond stretch using a charge flux formulation

**BNDERR Function**

"bnderr" is the distance bound error function and derivatives;
this version implements the original and Havel's normalized
lower bound penalty, the normalized version is preferred when
lower bounds are small (as with NMR NOE restraints), the
original penalty is needed if large lower bounds are present

**BNDGUESS Function**

"bndguess" sets approximate bond stretch force constants based
on atom type and connected atoms

**BONDS Subroutine**

"bonds" finds the total number of covalent bonds and
stores the atom numbers of the atoms defining each bond

**BORN Subroutine**

"born" computes the Born radius of each atom for use with
the various implicit solvation models

**BORN1 Subroutine**

"born1" computes derivatives of the Born radii with respect
to atomic coordinates and increments total energy derivatives
and virial components for potentials involving Born radii

**BOUNDS Subroutine**

"bounds" finds the center of mass of each molecule and
translates any stray molecules back into the periodic box

**BOXMIN Subroutine**

"boxmin" uses minimization of valence and vdw potential energy
to expand and refine a collection of solvent molecules in a
periodic box

**BOXMIN1 Function**

"boxmin1" is a service routine that computes the energy and
gradient during refinement of a periodic box

**BSET Subroutine**

"bset" computes by downward recursion the B functions used
in the evaluation of Slater-type (STO) overlap integrals

**BSPLGEN Subroutine**

"bsplgen" gets B-spline coefficients and derivatives for
a single PME atomic site along a particular direction

**BSPLINE Subroutine**

"bspline" calculates the coefficients for an n-th order
B-spline approximation

**BSPLINE_FILL Subroutine**

"bspline_fill" finds B-spline coefficients and derivatives
for PME atomic sites along the fractional coordinate axes

**BSSTEP Subroutine**

"bsstep" takes a single Bulirsch-Stoer step with monitoring
of local truncation error to ensure accuracy

**BUSSI Subroutine**

"bussi" performs a single molecular dynamics time step via
the Bussi-Parrinello isothermal-isobaric algorithm

**CALENDAR Subroutine**

"calendar" returns the current time as a set of integer values
representing the year, month, day, hour, minute and second

**CART_TO_FRAC Subroutine**

"cart_to_frac" computes a transformation matrix to convert
a multipole object in Cartesian coordinates to fractional

**CBUILD Subroutine**

"cbuild" performs a complete rebuild of the partial charge
electrostatic neighbor list for all sites

**CELLANG Subroutine**

"cellang" computes atomic coordinates and unit cell parameters
from fractional coordinates and lattice vectors

**CELLATOM Subroutine**

"cellatom" completes the addition of a symmetry related atom
to a unit cell by updating the atom type and attachment arrays

**CENTER Subroutine**

"center" moves the weighted centroid of each coordinate
set to the origin during least squares superposition

**CERROR Subroutine**

"cerror" is the error handling routine for the Connolly
surface area and volume computation

**CFFTB Subroutine**

"cfftb" computes the backward complex discrete Fourier
transform, the Fourier synthesis

**CFFTB1 Subroutine**

**CFFTF Subroutine**

"cfftf" computes the forward complex discrete Fourier
transform, the Fourier analysis

**CFFTF1 Subroutine**

**CFFTI Subroutine**

"cffti" initializes arrays used in both forward and backward
transforms; "ifac" is the prime factorization of "n", and
"wsave" contains a tabulation of trigonometric functions

**CFFTI1 Subroutine**

**CHIRER Function**

"chirer" computes the chirality error and its derivatives
with respect to atomic Cartesian coordinates as a sum the
squares of deviations of chiral volumes from target values

**CHKANGLE Subroutine**

"chkangle" tests angles to be constrained for their presence
in small rings and removes constraints that are redundant

**CHKAROM Function**

"chkatom" tests for the presence of a specified atom as a
member of an aromatic ring

**CHKPOLE Subroutine**

"chkpole" inverts atomic multipole moments as necessary
at sites with chiral local reference frame definitions

**CHKRING Subroutine**

"chkring" tests an atom or a set of connected atoms for
their presence within a single 3- to 6-membered ring

**CHKSIZE Subroutine**

"chksize" computes a measure of overall global structural
expansion or compaction from the number of excess upper
or lower bounds matrix violations

**CHKSOCKET Subroutine**

**CHKTREE Subroutine**

"chktree" tests a minimum energy structure to see if it
belongs to the correct progenitor in the existing map

**CHKTTOR Subroutine**

"chkttor" tests the attached atoms at a torsion-torsion central
site and inverts the angle values if the site is chiral

**CHKXYZ Subroutine**

"chkxyz" finds any pairs of atoms with identical Cartesian
coordinates, and prints a warning message

**CHOLESKY Subroutine**

"cholesky" uses a modified Cholesky method to solve the linear
system Ax = b, returning "x" in "b"; "A" is a real symmetric
positive definite matrix with its upper triangle (including the
diagonal) stored by rows

**CIRPLN Subroutine**

"cirpln" determines the points of intersection between a
specified circle and plane

**CJKM Function**

"cjkm" computes the coefficients of spherical harmonics
expressed in prolate spheroidal coordinates

**CLIGHT Subroutine**

"clight" performs a complete rebuild of the partial charge
pair neighbor list for all sites using the method of lights

**CLIMBER Subroutine**

**CLIMBRGD Subroutine**

**CLIMBROT Subroutine**

**CLIMBTOR Subroutine**

**CLIMBXYZ Subroutine**

**CLIST Subroutine**

"clist" performs an update or a complete rebuild of the
nonbonded neighbor lists for partial charges

**CLUSTER Subroutine**

"cluster" gets the partitioning of the system into groups
and stores a list of the group to which each atom belongs

**CMP_TO_FMP Subroutine**

"cmp_to_fmp" transforms the atomic multipoles from Cartesian
to fractional coordinates

**COLUMN Subroutine**

"column" takes the off-diagonal Hessian elements stored
as sparse rows and sets up indices to allow column access

**COMMAND Subroutine**

"command" uses the standard Unix-like iargc/getarg routines
to get the number and values of arguments specified on the
command line at program runtime

**COMPRESS Subroutine**

"compress" transfers only the non-buried tori from
the temporary tori arrays to the final tori arrays

**CONNECT Subroutine**

"connect" sets up the attached atom arrays
starting from a set of internal coordinates

**CONNOLLY Subroutine**

"connolly" uses the algorithms from the AMS/VAM programs of
Michael Connolly to compute the analytical molecular surface
area and volume of a collection of spherical atoms; thus
it implements Fred Richards' molecular surface definition as
a set of analytically defined spherical and toroidal polygons

**CONNYZE Subroutine**

"connyze" prints information onconnected atoms as lists
of all atom pairs that are 1-2 through 1-5 interactions

**CONTACT Subroutine**

"contact" constructs the contact surface, cycles and convex faces

**CONTROL Subroutine**

"control" gets initial values for parameters that determine
the output style and information level provided by Tinker

**COORDS Subroutine**

"coords" converts the three principal eigenvalues/vectors from
the metric matrix into atomic coordinates, and calls a routine
to compute the rms deviation from the bounds

**CORRELATE Program**

"correlate" computes the time correlation function of some
user-supplied property from individual snapshot frames taken
from a molecular dynamics or other trajectory

**CREATEJVM Subroutine**

**CREATESERVER Subroutine**

**CREATESYSTEM Subroutine**

**CREATEUPDATE Subroutine**

**CRYSTAL Program**

"crystal" is a utility which converts between fractional and
Cartesian coordinates, and can generate full unit cells from
asymmetric units

**CSPLINE Subroutine**

"cspline" computes the coefficients for a periodic interpolating
cubic spline

**CUTOFFS Subroutine**

"cutoffs" initializes and stores spherical energy cutoff
distance windows, Hessian element and Ewald sum cutoffs,
and allocates pairwise neighbor lists

**CYTSY Subroutine**

"cytsy" solves a system of linear equations for a cyclically
tridiagonal, symmetric, positive definite matrix

**CYTSYP Subroutine**

"cytsyp" finds the Cholesky factors of a cyclically tridiagonal
symmetric, positive definite matrix given by two vectors

**CYTSYS Subroutine**

"cytsys" solves a cyclically tridiagonal linear system
given the Cholesky factors

**D1D2 Function**

"d1d2" is a utility function used in computation of the
reaction field recursive summation elements

**DAMPDIR Subroutine**

"dampdir" generates coefficients for the direct field damping
function for powers of the interatomic distance

**DAMPEWALD Subroutine**

"dampewald" generates coefficients for Ewald error function
damping for powers of the interatomic distance

**DAMPMUT Subroutine**

"dampmut" generates coefficients for the mutual field damping
function for powers of the interatomic distance

**DAMPPOLAR Subroutine**

"damppolar" generates coefficients for the charge penetration
damping function used for polarization interactions

**DAMPPOLE Subroutine**

"damppole" generates coefficients for the charge penetration
damping function for powers of the interatomic distance

**DAMPPOT Subroutine**

"damppot" generates coefficients for the charge penetration
damping function used for the electrostatic potential

**DAMPREP Subroutine**

"damprep" generates coefficients for the Pauli repulsion
damping function for powers of the interatomic distance

**DAMPTHOLE Subroutine**

"dampthole" generates coefficients for the Thole damping
function for powers of the interatomic distance

**DBUILD Subroutine**

"dbuild" performs a complete rebuild of the damped dispersion
neighbor list for all sites

**DCFLUX Subroutine**

"dcflux" takes as input the electrostatic potential at each
atomic site and calculates gradient chain rule corrections due
to charge flux coupled with bond stretching and angle bending

**DEFLATE Subroutine**

"deflate" uses the power method with deflation to compute the
few largest eigenvalues and eigenvectors of a symmetric matrix

**DELETE Subroutine**

"delete" removes a specified atom from the Cartesian
coordinates list and shifts the remaining atoms

**DEPTH Function**

**DESTROYJVM Subroutine**

**DESTROYSERVER Subroutine**

**DFIELD0A Subroutine**

"dfield0a" computes the direct electrostatic field due to
permanent multipole moments via a double loop

**DFIELD0B Subroutine**

"dfield0b" computes the direct electrostatic field due to
permanent multipole moments via a pair list

**DFIELD0C Subroutine**

"dfield0c" computes the mutual electrostatic field due to
permanent multipole moments via Ewald summation

**DFIELD0D Subroutine**

"dfield0d" computes the direct electrostatic field due to
permanent multipole moments for use with with generalized
Kirkwood implicit solvation

**DFIELD0E Subroutine**

"dfield0e" computes the direct electrostatic field due to
permanent multipole moments for use with in Poisson-Boltzmann

**DFIELDI Subroutine**

"dfieldi" computes the electrostatic field due to permanent
multipole moments

**DFTMOD Subroutine**

"dftmod" computes the modulus of the discrete Fourier transform
of "bsarray" and stores it in "bsmod"

**DIAGBLK Subroutine**

"diagblk" performs diagonalization of the Hessian for a
block of atoms within a larger system

**DIAGQ Subroutine**

"diagq" is a matrix diagonalization routine which is derived
from the classical given, housec, and eigen algorithms with
several modifications to increase efficiency and accuracy

**DIFFEQ Subroutine**

"diffeq" performs the numerical integration of an ordinary
differential equation using an adaptive stepsize method to
solve the corresponding coupled first-order equations of the
general form dyi/dx = f(x,y1,...,yn) for yi = y1,...,yn

**DIFFUSE Program**

"diffuse" finds the self-diffusion constant for a homogeneous
liquid via the Einstein relation from a set of stored molecular
dynamics frames; molecular centers of mass are unfolded and mean
squared displacements are computed versus time separation

**DIST2 Function**

"dist2" finds the distance squared between two points; used
as a service routine by the Connolly surface area and volume
computation

**DISTGEOM Program**

"distgeom" uses a metric matrix distance geometry procedure to
generate structures with interpoint distances that lie within
specified bounds, with chiral centers that maintain chirality,
and with torsional angles restrained to desired values; the
user also has the ability to interactively inspect and alter
the triangle smoothed bounds matrix prior to embedding

**DLIGHT Subroutine**

"dlight" performs a complete rebuild of the damped dispersion
pair neighbor list for all sites using the method of lights

**DLIST Subroutine**

"dlist" performs an update or a complete rebuild of the
nonbonded neighbor lists for damped dispersion sites

**DMDUMP Subroutine**

"dmdump" puts the distance matrix of the final structure
into the upper half of a matrix, the distance of each atom
to the centroid on the diagonal, and the individual terms
of the bounds errors into the lower half of the matrix

**DOCUMENT Program**

"document" generates a formatted description of all the routines
and modules, an index of routines called by each source file, a
list of all valid keywords, a list of include file dependencies
as needed by a Unix-style Makefile, or a formatted force field
parameter summary

**DOT Function**

"dot" finds the dot product of two vectors

**DSTMAT Subroutine**

"dstmat" selects a distance matrix containing values between
the previously smoothed upper and lower bounds; the distance
values are chosen from uniform distributions, in a triangle
correlated fashion, or using random partial metrization

**DYNAMIC Program**

"dynamic" computes a molecular or stochastic dynamics trajectory
in one of the standard statistical mechanical ensembles and using
any of several possible integration methods

**EANGANG Subroutine**

"eangang" calculates the angle-angle potential energy

**EANGANG1 Subroutine**

"eangang1" calculates the angle-angle potential energy and
first derivatives with respect to Cartesian coordinates

**EANGANG2 Subroutine**

"eangang2" calculates the angle-angle potential energy
second derivatives with respect to Cartesian coordinates
using finite difference methods

**EANGANG2A Subroutine**

"eangang2a" calculates the angle-angle first derivatives for
a single interaction with respect to Cartesian coordinates;
used in computation of finite difference second derivatives

**EANGANG3 Subroutine**

"eangang3" calculates the angle-angle potential energy;
also partitions the energy among the atoms

**EANGLE Subroutine**

"eangle" calculates the angle bending potential energy;
projected in-plane angles at trigonal centers, special
linear or Fourier angle bending terms are optionally used

**EANGLE1 Subroutine**

"eangle1" calculates the angle bending potential energy and
the first derivatives with respect to Cartesian coordinates;
projected in-plane angles at trigonal centers, special linear
or Fourier angle bending terms are optionally used

**EANGLE2 Subroutine**

"eangle2" calculates second derivatives of the angle bending
energy for a single atom using a mixture of analytical and
finite difference methods; projected in-plane angles at trigonal
centers, special linear or Fourier angle bending terms are
optionally used

**EANGLE2A Subroutine**

"eangle2a" calculates bond angle bending potential energy
second derivatives with respect to Cartesian coordinates

**EANGLE2B Subroutine**

"eangle2b" computes projected in-plane bending first derivatives
for a single angle with respect to Cartesian coordinates;
used in computation of finite difference second derivatives

**EANGLE3 Subroutine**

"eangle3" calculates the angle bending potential energy, also
partitions the energy among the atoms; projected in-plane
angles at trigonal centers, spceial linear or Fourier angle
bending terms are optionally used

**EANGTOR Subroutine**

"eangtor" calculates the angle-torsion potential energy

**EANGTOR1 Subroutine**

"eangtor1" calculates the angle-torsion energy and first
derivatives with respect to Cartesian coordinates

**EANGTOR2 Subroutine**

"eangtor2" calculates the angle-torsion potential energy
second derivatives with respect to Cartesian coordinates

**EANGTOR3 Subroutine**

"eangtor3" calculates the angle-torsion potential energy;
also partitions the energy terms among the atoms

**EBOND Subroutine**

"ebond" calculates the bond stretching energy

**EBOND1 Subroutine**

"ebond1" calculates the bond stretching energy and
first derivatives with respect to Cartesian coordinates

**EBOND2 Subroutine**

"ebond2" calculates second derivatives of the bond
stretching energy for a single atom at a time

**EBOND3 Subroutine**

"ebond3" calculates the bond stretching energy; also
partitions the energy among the atoms

**EBUCK Subroutine**

"ebuck" calculates the Buckingham exp-6 van der Waals energy

**EBUCK0A Subroutine**

"ebuck0a" calculates the Buckingham exp-6 van der Waals energy
using a pairwise double loop

**EBUCK0B Subroutine**

"ebuck0b" calculates the Buckingham exp-6 van der Waals energy
using the method of lights

**EBUCK0C Subroutine**

"ebuck0c" calculates the Buckingham exp-6 van der Waals energy
using a pairwise neighbor list

**EBUCK0D Subroutine**

"ebuck0d" calculates the Buckingham exp-6 van der Waals energy
via a Gaussian approximation for potential energy smoothing

**EBUCK1 Subroutine**

"ebuck1" calculates the Buckingham exp-6 van der Waals energy
and its first derivatives with respect to Cartesian coordinates

**EBUCK1A Subroutine**

"ebuck1a" calculates the Buckingham exp-6 van der Waals energy
and its first derivatives using a pairwise double loop

**EBUCK1B Subroutine**

"ebuck1b" calculates the Buckingham exp-6 van der Waals energy
and its first derivatives using the method of lights

**EBUCK1C Subroutine**

"ebuck1c" calculates the Buckingham exp-6 van der Waals energy
and its first derivatives using a pairwise neighbor list

**EBUCK1D Subroutine**

"ebuck1d" calculates the Buckingham exp-6 van der Waals energy
and its first derivatives via a Gaussian approximation for
potential energy smoothing

**EBUCK2 Subroutine**

"ebuck2" calculates the Buckingham exp-6 van der Waals
second derivatives for a single atom at a time

**EBUCK2A Subroutine**

"ebuck2a" calculates the Buckingham exp-6 van der Waals second
derivatives using a double loop over relevant atom pairs

**EBUCK2B Subroutine**

"ebuck2b" calculates the Buckingham exp-6 van der Waals second
derivatives via a Gaussian approximation for use with potential
energy smoothing

**EBUCK3 Subroutine**

"ebuck3" calculates the Buckingham exp-6 van der Waals energy
and partitions the energy among the atoms

**EBUCK3A Subroutine**

"ebuck3a" calculates the Buckingham exp-6 van der Waals
energy and partitions the energy among the atoms using
a pairwise double loop

**EBUCK3B Subroutine**

"ebuck3b" calculates the Buckingham exp-6 van der Waals
energy and also partitions the energy among the atoms using
the method of lights

**EBUCK3C Subroutine**

"ebuck3c" calculates the Buckingham exp-6 van der Waals energy
and also partitions the energy among the atoms using a pairwise
neighbor list

**EBUCK3D Subroutine**

"ebuck3d" calculates the Buckingham exp-6 van der Waals energy
via a Gaussian approximation for potential energy smoothing

**ECHARGE Subroutine**

"echarge" calculates the charge-charge interaction energy

**ECHARGE0A Subroutine**

"echarge0a" calculates the charge-charge interaction energy
using a pairwise double loop

**ECHARGE0B Subroutine**

"echarge0b" calculates the charge-charge interaction energy
using the method of lights

**ECHARGE0C Subroutine**

"echarge0c" calculates the charge-charge interaction energy
using a pairwise neighbor list

**ECHARGE0D Subroutine**

"echarge0d" calculates the charge-charge interaction energy
using a particle mesh Ewald summation

**ECHARGE0E Subroutine**

"echarge0e" calculates the charge-charge interaction energy
using a particle mesh Ewald summation and the method of lights

**ECHARGE0F Subroutine**

"echarge0f" calculates the charge-charge interaction energy
using a particle mesh Ewald summation and a neighbor list

**ECHARGE0G Subroutine**

"echarge0g" calculates the charge-charge interaction energy
for use with potential smoothing methods

**ECHARGE1 Subroutine**

"echarge1" calculates the charge-charge interaction energy
and first derivatives with respect to Cartesian coordinates

**ECHARGE1A Subroutine**

"echarge1a" calculates the charge-charge interaction energy
and first derivatives with respect to Cartesian coordinates
using a pairwise double loop

**ECHARGE1B Subroutine**

"echarge1b" calculates the charge-charge interaction energy
and first derivatives with respect to Cartesian coordinates
using the method of lights

**ECHARGE1C Subroutine**

"echarge1c" calculates the charge-charge interaction energy
and first derivatives with respect to Cartesian coordinates
using a pairwise neighbor list

**ECHARGE1D Subroutine**

"echarge1d" calculates the charge-charge interaction energy
and first derivatives with respect to Cartesian coordinates
using a particle mesh Ewald summation

**ECHARGE1E Subroutine**

"echarge1e" calculates the charge-charge interaction energy
and first derivatives with respect to Cartesian coordinates
using a particle mesh Ewald summation and the method of lights

**ECHARGE1F Subroutine**

"echarge1f" calculates the charge-charge interaction energy
and first derivatives with respect to Cartesian coordinates
using a particle mesh Ewald summation and a neighbor list

**ECHARGE1G Subroutine**

"echarge1g" calculates the charge-charge interaction energy
and first derivatives with respect to Cartesian coordinates
for use with potential smoothing methods

**ECHARGE2 Subroutine**

"echarge2" calculates second derivatives of the
charge-charge interaction energy for a single atom

**ECHARGE2A Subroutine**

"echarge2a" calculates second derivatives of the charge-charge
interaction energy for a single atom using a pairwise loop

**ECHARGE2B Subroutine**

"echarge2b" calculates second derivatives of the charge-charge
interaction energy for a single atom using a neighbor list

**ECHARGE2C Subroutine**

"echarge2c" calculates second derivatives of the reciprocal
space charge-charge interaction energy for a single atom using
a particle mesh Ewald summation via numerical differentiation

**ECHARGE2D Subroutine**

"echarge2d" calculates second derivatives of the real space
charge-charge interaction energy for a single atom using a
pairwise loop

**ECHARGE2E Subroutine**

"echarge2e" calculates second derivatives of the real space
charge-charge interaction energy for a single atom using a
pairwise neighbor list

**ECHARGE2F Subroutine**

"echarge2f" calculates second derivatives of the charge-charge
interaction energy for a single atom for use with potential
smoothing methods

**ECHARGE2R Subroutine**

"echarge2r" computes reciprocal space charge-charge first
derivatives; used to get finite difference second derivatives

**ECHARGE3 Subroutine**

"echarge3" calculates the charge-charge interaction energy
and partitions the energy among the atoms

**ECHARGE3A Subroutine**

"echarge3a" calculates the charge-charge interaction energy
and partitions the energy among the atoms using a pairwise
double loop

**ECHARGE3B Subroutine**

"echarge3b" calculates the charge-charge interaction energy
and partitions the energy among the atoms using the method
of lights

**ECHARGE3C Subroutine**

"echarge3c" calculates the charge-charge interaction energy
and partitions the energy among the atoms using a pairwise
neighbor list

**ECHARGE3D Subroutine**

"echarge3d" calculates the charge-charge interaction energy
and partitions the energy among the atoms using a particle
mesh Ewald summation

**ECHARGE3E Subroutine**

"echarge3e" calculates the charge-charge interaction energy
and partitions the energy among the atoms using a particle
mesh Ewald summation and the method of lights

**ECHARGE3F Subroutine**

"echarge3f" calculates the charge-charge interaction energy
and partitions the energy among the atoms using a particle
mesh Ewald summation and a pairwise neighbor list

**ECHARGE3G Subroutine**

"echarge3g" calculates the charge-charge interaction energy
and partitions the energy among the atoms for use with
potential smoothing methods

**ECHGDPL Subroutine**

"echgdpl" calculates the charge-dipole interaction energy

**ECHGDPL1 Subroutine**

"echgdpl1" calculates the charge-dipole interaction energy
and first derivatives with respect to Cartesian coordinates

**ECHGDPL2 Subroutine**

"echgdpl2" calculates second derivatives of the
charge-dipole interaction energy for a single atom

**ECHGDPL3 Subroutine**

"echgdpl3" calculates the charge-dipole interaction energy;
also partitions the energy among the atoms

**ECHGTRN Subroutine**

"echgtrn" calculates the charge transfer potential energy

**ECHGTRN0A Subroutine**

"echgtrn0a" calculates the charge transfer interaction energy
using a double loop

**ECHGTRN0B Subroutine**

"echgtrn0b" calculates the charge transfer interaction energy
using the method of lights

**ECHGTRN0C Subroutine**

"echgtrn0c" calculates the charge transfer interaction energy
using a neighbor list

**ECHGTRN1 Subroutine**

"echgtrn1" calculates the charge transfer energy and first
derivatives with respect to Cartesian coordinates

**ECHGTRN1A Subroutine**

"echgtrn1a" calculates the charge transfer interaction energy
and first derivatives using a double loop

**ECHGTRN1B Subroutine**

"echgtrn1b" calculates the charge transfer energy and first
derivatives using a pairwise neighbor list

**ECHGTRN2 Subroutine**

"echgtrn2" calculates the second derivatives of the charge
transfer energy using a double loop over relevant atom pairs

**ECHGTRN3 Subroutine**

"echgtrn3" calculates the charge transfer energy; also partitions
the energy among the atoms

**ECHGTRN3A Subroutine**

"echgtrn3a" calculates the charge transfer interaction energy
and also partitions the energy among the atoms using a pairwise
double loop

**ECHGTRN3B Subroutine**

"echgtrn3b" calculates the charge transfer interaction energy
and also partitions the energy among the atoms using the method
of lights

**ECHGTRN3C Subroutine**

"echgtrn3c" calculates the charge transfer interaction energy
and also partitions the energy among the atoms using a pairwise
neighbor list

**ECRECIP Subroutine**

"ecrecip" evaluates the reciprocal space portion of the particle
mesh Ewald energy due to partial charges

**ECRECIP1 Subroutine**

"ecrecip1" evaluates the reciprocal space portion of the particle
mesh Ewald summation energy and gradient due to partial charges

**EDIFF Subroutine**

"ediff" calculates the energy of polarizing the vacuum induced
dipoles to their SCRF polarized values

**EDIFF1A Subroutine**

"ediff1a" calculates the energy and derivatives of polarizing
the vacuum induced dipoles to their SCRF polarized values using
a double loop

**EDIFF1B Subroutine**

"ediff1b" calculates the energy and derivatives of polarizing
the vacuum induced dipoles to their SCRF polarized values using
a neighbor list

**EDIFF3 Subroutine**

"ediff3" calculates the energy of polarizing the vacuum induced
dipoles to their generalized Kirkwood values with energy analysis

**EDIPOLE Subroutine**

"edipole" calculates the dipole-dipole interaction energy

**EDIPOLE1 Subroutine**

"edipole1" calculates the dipole-dipole interaction energy
and first derivatives with respect to Cartesian coordinates

**EDIPOLE2 Subroutine**

"edipole2" calculates second derivatives of the
dipole-dipole interaction energy for a single atom

**EDIPOLE3 Subroutine**

"edipole3" calculates the dipole-dipole interaction energy;
also partitions the energy among the atoms

**EDISP Subroutine**

"edisp" calculates the damped dispersion potential energy

**EDISP0A Subroutine**

"edisp0a" calculates the damped dispersion potential energy
using a pairwise double loop

**EDISP0B Subroutine**

"edisp0b" calculates the damped dispersion potential energy
using a pairwise neighbor list

**EDISP0C Subroutine**

"edisp0c" calculates the dispersion interaction energy using
particle mesh Ewald summation and a double loop

**EDISP0D Subroutine**

"edisp0d" calculates the dispersion interaction energy using
particle mesh Ewald summation and a neighbor list

**EDISP1 Subroutine**

"edisp1" calculates the damped dispersion energy and first
derivatives with respect to Cartesian coordinates

**EDISP1A Subroutine**

"edisp1a" calculates the damped dispersion energy and
derivatives with respect to Cartesian coordinates using
a pairwise double loop

**EDISP1B Subroutine**

"edisp1b" calculates the damped dispersion energy and
derivatives with respect to Cartesian coordinates using
a pairwise neighbor list

**EDISP1C Subroutine**

"edisp1c" calculates the damped dispersion energy and
derivatives with respect to Cartesian coordinates using
particle mesh Ewald summation and a double loop

**EDISP1D Subroutine**

"edisp1d" calculates the damped dispersion energy and
derivatives with respect to Cartesian coordinates using
particle mesh Ewald summation and a neighbor list

**EDISP2 Subroutine**

"edisp2" calculates the damped dispersion second derivatives
for a single atom at a time

**EDISP3 Subroutine**

"edisp3" calculates the dispersion energy; also partitions
the energy among the atoms

**EDISP3A Subroutine**

"edisp3a" calculates the dispersion potential energy and
also partitions the energy among the atoms using a pairwise
double loop

**EDISP3B Subroutine**

"edisp3b" calculates the damped dispersion potential energy
and also partitions the energy among the atomsusing a pairwise
neighbor list

**EDISP3C Subroutine**

"edisp3c" calculates the dispersion interaction energy using
particle mesh Ewald summation and a double loop

**EDISP3D Subroutine**

"edisp3d" calculates the damped dispersion energy and analysis
using particle mesh Ewald summation and a neighbor list

**EDREAL0C Subroutine**

"edreal0c" calculates the damped dispersion potential energy
using a particle mesh Ewald sum and pairwise double loop

**EDREAL0D Subroutine**

"edreal0d" evaluated the real space portion of the damped
dispersion energy using a neighbor list

**EDREAL1C Subroutine**

"edreal1c" evaluates the real space portion of the Ewald
summation energy and gradient due to damped dispersion
interactions via a double loop

**EDREAL1D Subroutine**

"edreal1d" evaluates the real space portion of the Ewald
summation energy and gradient due to damped dispersion
interactions via a neighbor list

**EDREAL3C Subroutine**

"edreal3c" calculates the real space portion of the damped
dispersion energy and analysis using Ewald and a double loop

**EDREAL3D Subroutine**

"edreal3d" evaluated the real space portion of the damped
dispersion energy and analysis using Ewald and a neighbor list

**EDRECIP Subroutine**

"edrecip" evaluates the reciprocal space portion of the particle
mesh Ewald energy due to damped dispersion

**EDRECIP1 Subroutine**

"edrecip1" evaluates the reciprocal space portion of particle
mesh Ewald energy and gradient due to damped dispersion

**EGAUSS Subroutine**

"egauss" calculates the Gaussian expansion van der Waals energy

**EGAUSS0A Subroutine**

"egauss0a" calculates the Gaussian expansion van der Waals
energy using a pairwise double loop

**EGAUSS0B Subroutine**

"egauss0b" calculates the Gaussian expansion van der Waals energy
using the method of lights

**EGAUSS0C Subroutine**

"egauss0c" calculates the Gaussian expansion van der Waals
energy using a pairwise neighbor list

**EGAUSS0D Subroutine**

"egauss0d" calculates the Gaussian expansion van der Waals
energy for use with potential energy smoothing

**EGAUSS1 Subroutine**

"egauss1" calculates the Gaussian expansion van der Waals
interaction energy and its first derivatives with respect
to Cartesian coordinates

**EGAUSS1A Subroutine**

"egauss1a" calculates the Gaussian expansion van der Waals
interaction energy and its first derivatives using a pairwise
double loop

**EGAUSS1B Subroutine**

"egauss1b" calculates the Gaussian expansion van der Waals
energy and its first derivatives with respect to Cartesian
coordinates using the method of lights

**EGAUSS1C Subroutine**

"egauss1c" calculates the Gaussian expansion van der Waals
energy and its first derivatives with respect to Cartesian
coordinates using a pairwise neighbor list

**EGAUSS1D Subroutine**

"egauss1d" calculates the Gaussian expansion van der Waals
interaction energy and its first derivatives for use with
potential energy smoothing

**EGAUSS2 Subroutine**

"egauss2" calculates the Gaussian expansion van der Waals
second derivatives for a single atom at a time

**EGAUSS2A Subroutine**

"egauss2a" calculates the Gaussian expansion van der Waals
second derivatives using a pairwise double loop

**EGAUSS2B Subroutine**

"egauss2b" calculates the Gaussian expansion van der Waals
second derivatives for use with potential energy smoothing

**EGAUSS3 Subroutine**

"egauss3" calculates the Gaussian expansion van der Waals
interaction energy and partitions the energy among the atoms

**EGAUSS3A Subroutine**

"egauss3a" calculates the Gaussian expansion van der Waals
energy and partitions the energy among the atoms using a
pairwise double loop

**EGAUSS3B Subroutine**

"egauss3b" calculates the Gaussian expansion van der Waals
energy and partitions the energy among the atoms using the
method of lights

**EGAUSS3C Subroutine**

"egauss3c" calculates the Gaussian expansion van der Waals
energy and partitions the energy among the atoms using a
pairwise neighbor list

**EGAUSS3D Subroutine**

"egauss3d" calculates the Gaussian expansion van der Waals
interaction energy and partitions the energy among the atoms
for use with potential energy smoothing

**EGB0A Subroutine**

"egb0a" calculates the generalized Born polarization energy
for the GB/SA solvation models using a pairwise double loop

**EGB0B Subroutine**

"egb0b" calculates the generalized Born polarization energy
for the GB/SA solvation models using a pairwise neighbor list

**EGB0C Subroutine**

"egb0c" calculates the generalized Born polarization energy
for the GB/SA solvation models for use with potential smoothing
methods via analogy to the smoothing of Coulomb's law

**EGB1A Subroutine**

"egb1a" calculates the generalized Born electrostatic energy
and first derivatives of the GB/SA solvation models using a
double loop

**EGB1B Subroutine**

"egb1b" calculates the generalized Born electrostatic energy
and first derivatives of the GB/SA solvation models using a
neighbor list

**EGB1C Subroutine**

"egb1c" calculates the generalized Born energy and first
derivatives of the GB/SA solvation models for use with
potential smoothing methods

**EGB2A Subroutine**

"egb2a" calculates second derivatives of the generalized
Born energy term for the GB/SA solvation models

**EGB2B Subroutine**

"egb2b" calculates second derivatives of the generalized
Born energy term for the GB/SA solvation models for use with
potential smoothing methods

**EGB3A Subroutine**

"egb3a" calculates the generalized Born electrostatic energy
for GB/SA solvation models using a pairwise double loop; also
partitions the energy among the atoms

**EGB3B Subroutine**

"egb3b" calculates the generalized Born electrostatic energy
for GB/SA solvation models using a pairwise neighbor list; also
partitions the energy among the atoms

**EGB3C Subroutine**

"egb3c" calculates the generalized Born electrostatic energy
for GB/SA solvation models for use with potential smoothing
methods via analogy to the smoothing of Coulomb's law; also
partitions the energy among the atoms

**EGEOM Subroutine**

"egeom" calculates the energy due to restraints on positions,
distances, angles and torsions as well as Gaussian basin and
spherical droplet restraints

**EGEOM1 Subroutine**

"egeom1" calculates the energy and first derivatives
with respect to Cartesian coordinates due to restraints
on positions, distances, angles and torsions as well as
Gaussian basin and spherical droplet restraints

**EGEOM2 Subroutine**

"egeom2" calculates second derivatives of restraints
on positions, distances, angles and torsions as well
as Gaussian basin and spherical droplet restraints

**EGEOM3 Subroutine**

"egeom3" calculates the energy due to restraints on positions,
distances, angles and torsions as well as Gaussian basin and
droplet restraints; also partitions energy among the atoms

**EGK Subroutine**

"egk" calculates the generalized Kirkwood electrostatic
solvation free energy for the GK/NP implicit solvation model

**EGK0A Subroutine**

"egk0a" calculates the electrostatic portion of the implicit
solvation energy via the generalized Kirkwood model

**EGK1 Subroutine**

"egk1" calculates the implicit solvation energy and derivatives
via the generalized Kirkwood plus nonpolar implicit solvation

**EGK1A Subroutine**

"egk1a" calculates the electrostatic portion of the implicit
solvation energy and derivatives via the generalized Kirkwood
model

**EGK3 Subroutine**

"egk3" calculates the generalized Kirkwood electrostatic
energy for GK/NP solvation models; also partitions the
energy among the atoms

**EGK3A Subroutine**

"egk3a" calculates the electrostatic portion of the implicit
solvation energy via the generalized Kirkwood model; also
partitions the energy among the atoms

**EHAL Subroutine**

"ehal" calculates the buffered 14-7 van der Waals energy

**EHAL0A Subroutine**

"ehal0a" calculates the buffered 14-7 van der Waals energy
using a pairwise double loop

**EHAL0B Subroutine**

"ehal0b" calculates the buffered 14-7 van der Waals energy
using the method of lights

**EHAL0C Subroutine**

"ehal0c" calculates the buffered 14-7 van der Waals energy
using a pairwise neighbor list

**EHAL1 Subroutine**

"ehal1" calculates the buffered 14-7 van der Waals energy and
its first derivatives with respect to Cartesian coordinates

**EHAL1A Subroutine**

"ehal1a" calculates the buffered 14-7 van der Waals energy and
its first derivatives with respect to Cartesian coordinates
using a pairwise double loop

**EHAL1B Subroutine**

"ehal1b" calculates the buffered 14-7 van der Waals energy and
its first derivatives with respect to Cartesian coordinates
using the method of lights

**EHAL1C Subroutine**

"ehal1c" calculates the buffered 14-7 van der Waals energy and
its first derivatives with respect to Cartesian coordinates
using a pairwise neighbor list

**EHAL2 Subroutine**

"ehal2" calculates the buffered 14-7 van der Waals second
derivatives for a single atom at a time

**EHAL3 Subroutine**

"ehal3" calculates the buffered 14-7 van der Waals energy
and partitions the energy among the atoms

**EHAL3A Subroutine**

"ehal3a" calculates the buffered 14-7 van der Waals energy
and partitions the energy among the atoms using a pairwise
double loop

**EHAL3B Subroutine**

"ehal3b" calculates the buffered 14-7 van der Waals energy
and also partitions the energy among the atoms using the
method of lights

**EHAL3C Subroutine**

"ehal3c" calculates the buffered 14-7 van der Waals energy
and also partitions the energy among the atoms using a
pairwise neighbor list

**EHPMF Subroutine**

"ehpmf" calculates the hydrophobic potential of mean force
energy using a pairwise double loop

**EHPMF1 Subroutine**

"ehpmf1" calculates the hydrophobic potential of mean force
energy and first derivatives using a pairwise double loop

**EHPMF3 Subroutine**

"ehpmf3" calculates the hydrophobic potential of mean force
nonpolar energy; also partitions the energy among the atoms

**EIGEN Subroutine**

"eigen" uses the power method to compute the largest eigenvalues
and eigenvectors of the metric matrix, "valid" is set true if the
first three eigenvalues are positive

**EIGENRGD Subroutine**

**EIGENROT Subroutine**

**EIGENROT Subroutine**

**EIGENTOR Subroutine**

**EIGENXYZ Subroutine**

**EIMPROP Subroutine**

"eimprop" calculates the improper dihedral potential energy

**EIMPROP1 Subroutine**

"eimprop1" calculates improper dihedral energy and its
first derivatives with respect to Cartesian coordinates

**EIMPROP2 Subroutine**

"eimprop2" calculates second derivatives of the improper
dihedral angle energy for a single atom

**EIMPROP3 Subroutine**

"eimprop3" calculates the improper dihedral potential
energy; also partitions the energy terms among the atoms

**EIMPTOR Subroutine**

"eimptor" calculates the improper torsion potential energy

**EIMPTOR1 Subroutine**

"eimptor1" calculates improper torsion energy and its
first derivatives with respect to Cartesian coordinates

**EIMPTOR2 Subroutine**

"eimptor2" calculates second derivatives of the improper
torsion energy for a single atom

**EIMPTOR3 Subroutine**

"eimptor3" calculates the improper torsion potential energy;
also partitions the energy terms among the atoms

**ELJ Subroutine**

"elj" calculates the Lennard-Jones 6-12 van der Waals energy

**ELJ0A Subroutine**

"elj0a" calculates the Lennard-Jones 6-12 van der Waals energy
using a pairwise double loop

**ELJ0B Subroutine**

"elj0b" calculates the Lennard-Jones 6-12 van der Waals energy
using the method of lights

**ELJ0C Subroutine**

"elj0c" calculates the Lennard-Jones 6-12 van der Waals energy
using a pairwise neighbor list

**ELJ0D Subroutine**

"elj0d" calculates the Lennard-Jones 6-12 van der Waals energy
via a Gaussian approximation for potential energy smoothing

**ELJ0E Subroutine**

"elj0e" calculates the Lennard-Jones 6-12 van der Waals energy
for use with stophat potential energy smoothing

**ELJ1 Subroutine**

"elj1" calculates the Lennard-Jones 6-12 van der Waals energy
and its first derivatives with respect to Cartesian coordinates

**ELJ1A Subroutine**

"elj1a" calculates the Lennard-Jones 6-12 van der Waals energy
and its first derivatives using a pairwise double loop

**ELJ1B Subroutine**

"elj1b" calculates the Lennard-Jones 6-12 van der Waals energy
and its first derivatives using the method of lights

**ELJ1C Subroutine**

"elj1c" calculates the Lennard-Jones 12-6 van der Waals energy
and its first derivatives using a pairwise neighbor list

**ELJ1D Subroutine**

"elj1d" calculates the Lennard-Jones 6-12 van der Waals energy
 and its first derivatives via a Gaussian approximation for
 potential energy smoothing

**ELJ1E Subroutine**

"elj1e" calculates the van der Waals interaction energy and its
first derivatives for use with stophat potential energy smoothing

**ELJ2 Subroutine**

"elj2" calculates the Lennard-Jones 6-12 van der Waals second
derivatives for a single atom at a time

**ELJ2A Subroutine**

"elj2a" calculates the Lennard-Jones 6-12 van der Waals second
derivatives using a double loop over relevant atom pairs

**ELJ2B Subroutine**

"elj2b" calculates the Lennard-Jones 6-12 van der Waals second
derivatives via a Gaussian approximation for use with potential
energy smoothing

**ELJ2C Subroutine**

"elj2c" calculates the Lennard-Jones 6-12 van der Waals second
derivatives for use with stophat potential energy smoothing

**ELJ3 Subroutine**

"elj3" calculates the Lennard-Jones 6-12 van der Waals energy
and also partitions the energy among the atoms

**ELJ3A Subroutine**

"elj3a" calculates the Lennard-Jones 6-12 van der Waals
energy and also partitions the energy among the atoms using
a pairwise double loop

**ELJ3B Subroutine**

"elj3b" calculates the Lennard-Jones 6-12 van der Waals
energy and also partitions the energy among the atoms using
the method of lights

**ELJ3C Subroutine**

"elj3c" calculates the Lennard-Jones van der Waals energy
and also partitions the energy among the atoms using a
pairwise neighbor list

**ELJ3D Subroutine**

"elj3d" calculates the Lennard-Jones 6-12 van der Waals energy
and also partitions the energy among the atoms via a Gaussian
approximation for potential energy smoothing

**ELJ3E Subroutine**

"elj3e" calculates the Lennard-Jones 6-12 van der Waals energy
and also partitions the energy among the atoms for use with
stophat potential energy smoothing

**EMBED Subroutine**

"embed" is a distance geometry routine patterned after the
ideas of Gordon Crippen, Irwin Kuntz and Tim Havel; it takes
as input a set of upper and lower bounds on the interpoint
distances, chirality restraints and torsional restraints,
and attempts to generate a set of coordinates that satisfy
the input bounds and restraints

**EMETAL Subroutine**

"emetal" calculates the transition metal ligand field energy

**EMETAL1 Subroutine**

"emetal1" calculates the transition metal ligand field energy
and its first derivatives with respect to Cartesian coordinates

**EMETAL2 Subroutine**

"emetal2" calculates the transition metal ligand field second
derivatives for a single atom at a time

**EMETAL3 Subroutine**

"emetal3" calculates the transition metal ligand field energy
and also partitions the energy among the atoms

**EMM3HB Subroutine**

"emm3hb" calculates the MM3 exp-6 van der Waals and directional
charge transfer hydrogen bonding energy

**EMM3HB0A Subroutine**

"emm3hb0a" calculates the MM3 exp-6 van der Waals and
directional charge transfer hydrogen bonding energy using
a pairwise double loop

**EMM3HB0B Subroutine**

"emm3hb0b" calculates the MM3 exp-6 van der Waals and
directional charge transfer hydrogen bonding energy using
the method of lights

**EMM3HB0C Subroutine**

"emm3hb0c" calculates the MM3 exp-6 van der Waals and
directional charge transfer hydrogen bonding energy using
a pairwise neighbor list

**EMM3HB1 Subroutine**

"emm3hb1" calculates the MM3 exp-6 van der Waals and directional
charge transfer hydrogen bonding energy with respect to Cartesian
coordinates

**EMM3HB1A Subroutine**

"emm3hb1a" calculates the MM3 exp-6 van der Waals and directional
charge transfer hydrogen bonding energy with respect to Cartesian
coordinates using a pairwise double loop

**EMM3HB1B Subroutine**

"emm3hb1b" calculates the MM3 exp-6 van der Waals and directional
charge transfer hydrogen bonding energy with respect to Cartesian
coordinates using the method of lights

**EMM3HB1C Subroutine**

"emm3hb1c" calculates the MM3 exp-6 van der Waals and directional
charge transfer hydrogen bonding energy with respect to Cartesian
coordinates using a pairwise neighbor list

**EMM3HB2 Subroutine**

"emm3hb2" calculates the MM3 exp-6 van der Waals and directional
charge transfer hydrogen bonding second derivatives for a single
atom at a time

**EMM3HB3 Subroutine**

"emm3hb3" calculates the MM3 exp-6 van der Waals and directional
charge transfer hydrogen bonding energy, and partitions the energy
among the atoms

**EMM3HB3A Subroutine**

"emm3hb3" calculates the MM3 exp-6 van der Waals and
directional charge transfer hydrogen bonding energy, and
partitions the energy among the atoms

**EMM3HB3B Subroutine**

"emm3hb3b" calculates the MM3 exp-6 van der Waals and
directional charge transfer hydrogen bonding energy using
the method of lights

**EMM3HB3C Subroutine**

"emm3hb3c" calculates the MM3 exp-6 van der Waals and
directional charge transfer hydrogen bonding energy using
a pairwise neighbor list

**EMPOLE Subroutine**

"empole" calculates the electrostatic energy due to atomic
multipole interactions

**EMPOLE0A Subroutine**

"empole0a" calculates the atomic multipole interaction energy
using a double loop

**EMPOLE0B Subroutine**

"empole0b" calculates the atomic multipole interaction energy
using a neighbor list

**EMPOLE0C Subroutine**

"empole0c" calculates the atomic multipole interaction energy
using particle mesh Ewald summation and a double loop

**EMPOLE0D Subroutine**

"empole0d" calculates the atomic multipole interaction energy
using particle mesh Ewald summation and a neighbor list

**EMPOLE1 Subroutine**

"empole1" calculates the atomic multipole energy and first
derivatives with respect to Cartesian coordinates

**EMPOLE1A Subroutine**

"empole1a" calculates the multipole energy and derivatives with
respect to Cartesian coordinates using a pairwise double loop

**EMPOLE1B Subroutine**

"empole1b" calculates the multipole energy and derivatives
with respect to Cartesian coordinates using a neighbor list

**EMPOLE1C Subroutine**

"empole1c" calculates the multipole energy and derivatives
with respect to Cartesian coordinates using particle mesh
Ewald summation and a double loop

**EMPOLE1D Subroutine**

"empole1d" calculates the multipole energy and derivatives
with respect to Cartesian coordinates using particle mesh Ewald
summation and a neighbor list

**EMPOLE2 Subroutine**

"empole2" calculates second derivatives of the multipole energy
for a single atom at a time

**EMPOLE2A Subroutine**

"empole2a" computes multipole first derivatives for a single
atom; used to get finite difference second derivatives

**EMPOLE3 Subroutine**

"empole3" calculates the electrostatic energy due to atomic
multipole interactions, and partitions the energy among atoms

**EMPOLE3A Subroutine**

"empole3a" calculates the atomic multipole interaction energy
using a double loop, and partitions the energy among atoms

**EMPOLE3B Subroutine**

"empole3b" calculates the atomic multipole interaction energy
using a neighbor list, and partitions the energy among the atoms

**EMPOLE3C Subroutine**

"empole3c" calculates the atomic multipole interaction energy
using a particle mesh Ewald summation and double loop, and
partitions the energy among the atoms

**EMPOLE3D Subroutine**

"empole3d" calculates the atomic multipole interaction energy
using particle mesh Ewald summation and a neighbor list, and
partitions the energy among the atoms

**EMREAL0C Subroutine**

"emreal0c" evaluates the real space portion of the Ewald sum
energy due to atomic multipoles using a double loop

**EMREAL0D Subroutine**

"emreal0d" evaluates the real space portion of the Ewald sum
energy due to atomic multipoles using a neighbor list

**EMREAL1C Subroutine**

"emreal1c" evaluates the real space portion of the Ewald
summation energy and gradient due to multipole interactions
via a double loop

**EMREAL1D Subroutine**

"emreal1d" evaluates the real space portion of the Ewald
summation energy and gradient due to multipole interactions
via a neighbor list

**EMREAL3C Subroutine**

"emreal3c" evaluates the real space portion of the Ewald sum
energy due to atomic multipole interactions and partitions
the energy among the atoms

**EMREAL3D Subroutine**

"emreal3d" evaluates the real space portion of the Ewald sum
energy due to atomic multipole interactions, and partitions
the energy among the atoms using a pairwise neighbor list

**EMRECIP Subroutine**

"emrecip" evaluates the reciprocal space portion of the particle
mesh Ewald energy due to atomic multipole interactions

**EMRECIP1 Subroutine**

"emrecip1" evaluates the reciprocal space portion of particle
mesh Ewald summation energy and gradient due to multipoles

**ENERGY Function**

"energy" calls the subroutines to calculate the potential
energy terms and sums up to form the total energy

**ENP Subroutine**

"enp" calculates the nonpolar implicit solvation energy
as a sum of cavity and dispersion terms

**ENP1 Subroutine**

"enp1" calculates the nonpolar implicit solvation energy
and derivatives as a sum of cavity and dispersion terms

**ENP3 Subroutine**

"enp3" calculates the nonpolar implicit solvation energy as
a sum of cavity and dispersion terms; also partitions the
energy among the atoms

**ENRGYZE Subroutine**

"enrgyze" is an auxiliary routine for the analyze program
that performs the energy analysis and prints the total and
intermolecular energies

**EOPBEND Subroutine**

"eopbend" computes the out-of-plane bend potential energy at
trigonal centers via a Wilson-Decius-Cross or Allinger angle

**EOPBEND1 Subroutine**

"eopbend1" computes the out-of-plane bend potential energy and
first derivatives at trigonal centers via a Wilson-Decius-Cross
or Allinger angle

**EOPBEND2 Subroutine**

"eopbend2" calculates second derivatives of the out-of-plane
bend energy via a Wilson-Decius-Cross or Allinger angle for
a single atom using finite difference methods

**EOPBEND2A Subroutine**

"eopbend2a" calculates out-of-plane bend first derivatives at
a trigonal center via a Wilson-Decius-Cross or Allinger angle;
used in computation of finite difference second derivatives

**EOPBEND3 Subroutine**

"eopbend3" computes the out-of-plane bend potential energy at
trigonal centers via a Wilson-Decius-Cross or Allinger angle;
also partitions the energy among the atoms

**EOPDIST Subroutine**

"eopdist" computes the out-of-plane distance potential
energy at trigonal centers via the central atom height

**EOPDIST1 Subroutine**

"eopdist1" computes the out-of-plane distance potential
energy and first derivatives at trigonal centers via
the central atom height

**EOPDIST2 Subroutine**

"eopdist2" calculates second derivatives of the out-of-plane
distance energy for a single atom via the central atom height

**EOPDIST3 Subroutine**

"eopdist3" computes the out-of-plane distance potential energy
at trigonal centers via the central atom height; also partitions
the energy among the atoms

**EPB Subroutine**

"epb" calculates the implicit solvation energy via the
Poisson-Boltzmann plus nonpolar implicit solvation

**EPB1 Subroutine**

"epb1" calculates the implicit solvation energy and derivatives
via the Poisson-Boltzmann plus nonpolar implicit solvation

**EPB1A Subroutine**

"epb1a" calculates the solvation energy and gradients for the
PB/NP solvation model

**EPB3 Subroutine**

"epb3" calculates the implicit solvation energy via the
Poisson-Boltzmann model; also partitions the energy among
the atoms

**EPITORS Subroutine**

"epitors" calculates the pi-system torsion potential energy

**EPITORS1 Subroutine**

"epitors1" calculates the pi-system torsion potential energy
and first derivatives with respect to Cartesian coordinates

**EPITORS2 Subroutine**

"epitors2" calculates the second derivatives of the pi-system
torsion energy for a single atom using finite difference methods

**EPITORS2A Subroutine**

"epitors2a" calculates the pi-system torsion first derivatives;
used in computation of finite difference second derivatives

**EPITORS3 Subroutine**

"epitors3" calculates the pi-system torsion potential energy;
also partitions the energy terms among the atoms

**EPOLAR Subroutine**

"epolar" calculates the polarization energy due to induced
dipole interactions

**EPOLAR0A Subroutine**

"epolar0a" calculates the induced dipole polarization energy
using a double loop, and partitions the energy among atoms

**EPOLAR0B Subroutine**

"epolar0b" calculates the induced dipole polarization energy
using a neighbor list

**EPOLAR0C Subroutine**

"epolar0c" calculates the dipole polarization energy with respect
to Cartesian coordinates using particle mesh Ewald summation and
a double loop

**EPOLAR0D Subroutine**

"epolar0d" calculates the dipole polarization energy with respect
to Cartesian coordinates using particle mesh Ewald summation and
a neighbor list

**EPOLAR0E Subroutine**

"epolar0e" calculates the dipole polarizability interaction
from the induced dipoles times the electric field

**EPOLAR1 Subroutine**

"epolar1" calculates the induced dipole polarization energy
and first derivatives with respect to Cartesian coordinates

**EPOLAR1A Subroutine**

"epolar1a" calculates the dipole polarization energy and
derivatives with respect to Cartesian coordinates using a
pairwise double loop

**EPOLAR1B Subroutine**

"epolar1b" calculates the dipole polarization energy and
derivatives with respect to Cartesian coordinates using a
neighbor list

**EPOLAR1C Subroutine**

"epolar1c" calculates the dipole polarization energy and
derivatives with respect to Cartesian coordinates using
particle mesh Ewald summation and a double loop

**EPOLAR1D Subroutine**

"epolar1d" calculates the dipole polarization energy and
derivatives with respect to Cartesian coordinates using
particle mesh Ewald summation and a neighbor list

**EPOLAR1E Subroutine**

"epolar1e" calculates the dipole polarizability interaction
from the induced dipoles times the electric field

**EPOLAR2 Subroutine**

"epolar2" calculates second derivatives of the dipole polarization
energy for a single atom at a time

**EPOLAR2A Subroutine**

"epolar2a" computes polarization first derivatives for a single
atom with respect to Cartesian coordinates; used to get finite
difference second derivatives

**EPOLAR3 Subroutine**

"epolar3" calculates the induced dipole polarization energy,
and partitions the energy among atoms

**EPOLAR3A Subroutine**

"epolar3a" calculates the induced dipole polarization energy
using a double loop, and partitions the energy among atoms

**EPOLAR3B Subroutine**

"epolar3b" calculates the induced dipole polarization energy
using a neighbor list, and partitions the energy among atoms

**EPOLAR3C Subroutine**

"epolar3c" calculates the polarization energy and analysis with
respect to Cartesian coordinates using particle mesh Ewald and
a double loop

**EPOLAR3D Subroutine**

"epolar3d" calculates the polarization energy and analysis with
respect to Cartesian coordinates using particle mesh Ewald and
a neighbor list

**EPOLAR3E Subroutine**

"epolar3e" calculates the dipole polarizability interaction
from the induced dipoles times the electric field

**EPREAL0C Subroutine**

"epreal0c" calculates the induced dipole polarization energy
using particle mesh Ewald summation and a double loop

**EPREAL0D Subroutine**

"epreal0d" calculates the induced dipole polarization energy
using particle mesh Ewald summation and a neighbor list

**EPREAL1C Subroutine**

"epreal1c" evaluates the real space portion of the Ewald
summation energy and gradient due to dipole polarization
via a double loop

**EPREAL1D Subroutine**

"epreal1d" evaluates the real space portion of the Ewald
summation energy and gradient due to dipole polarization
via a neighbor list

**EPREAL3C Subroutine**

"epreal3c" calculates the induced dipole polarization energy and
analysis using particle mesh Ewald summation and a double loop

**EPREAL3D Subroutine**

"epreal3d" calculates the induced dipole polarization energy
and analysis using particle mesh Ewald and a neighbor list

**EPRECIP Subroutine**

"eprecip" evaluates the reciprocal space portion of particle
mesh Ewald summation energy due to dipole polarization

**EPRECIP1 Subroutine**

"eprecip1" evaluates the reciprocal space portion of the particle
mesh Ewald summation energy and gradient due to dipole polarization

**EQUCLC Subroutine**

**EREPEL Subroutine**

"erepel" calculates the Pauli exchange repulsion energy

**EREPEL0A Subroutine**

"erepel0a" calculates the Pauli repulsion interaction energy
using a double loop

**EREPEL0B Subroutine**

"erepel0b" calculates the Pauli repulsion interaction energy
using a pairwise neighbor list

**EREPEL1 Subroutine**

"erepel1" calculates the Pauli repulsion energy and first
derivatives with respect to Cartesian coordinates

**EREPEL1A Subroutine**

"erepel1a" calculates the Pauli repulsion energy and first
derivatives with respect to Cartesian coordinates using a
pairwise double loop

**EREPEL1B Subroutine**

"erepel1b" calculates the Pauli repulsion energy and first
derivatives with respect to Cartesian coordinates using a
pariwise neighbor list

**EREPEL2 Subroutine**

"erepel2" calculates the second derivatives of the Pauli
repulsion energy

**EREPEL2A Subroutine**

"erepel2a" computes Pauli repulsion first derivatives for a
single atom via a double loop; used to get finite difference
second derivatives

**EREPEL3 Subroutine**

"erepel3" calculates the Pauli repulsion energy and partitions
the energy among the atoms

**EREPEL3A Subroutine**

"erepel3a" calculates the Pauli repulsion energy and also
partitions the energy among the atoms using a double loop

**EREPEL3B Subroutine**

"erepel3b" calculates the Pauli repulsion energy and also
partitions the energy among the atoms using a neighbor list

**ERF Function**

"erf" computes a numerical approximation to the value of
the error function via a Chebyshev approximation

**ERFC Function**

"erfc" computes a numerical approximation to the value of the
complementary error function via a Chebyshev approximation

**ERFCORE Subroutine**

"erfcore" evaluates erf(x) or erfc(x) for a real argument x;
when called with mode set to 0 it returns erf, a mode of 1
returns erfc; uses rational functions that approximate erf(x)
and erfc(x) to at least 18 significant decimal digits

**ERFIK Subroutine**

"erfik" compute the reaction field energy due to a single pair
of atomic multipoles

**ERFINV Function**

"erfinv" evaluates the inverse of the error function for
an argument in the range (-1,1) using a rational function
approximation followed by cycles of Newton-Raphson correction

**ERXNFLD Subroutine**

"erxnfld" calculates the macroscopic reaction field energy
arising from a set of atomic multipoles

**ERXNFLD1 Subroutine**

"erxnfld1" calculates the macroscopic reaction field energy
and derivatives with respect to Cartesian coordinates

**ERXNFLD2 Subroutine**

"erxnfld2" calculates second derivatives of the macroscopic
reaction field energy for a single atom at a time

**ERXNFLD3 Subroutine**

"erxnfld3" calculates the macroscopic reaction field energy,
and also partitions the energy among the atoms

**ESOLV Subroutine**

"esolv" calculates the implicit solvation energy for surface area,
generalized Born, generalized Kirkwood and Poisson-Boltzmann
solvation models

**ESOLV1 Subroutine**

"esolv1" calculates the implicit solvation energy and
first derivatives with respect to Cartesian coordinates
for surface area, generalized Born, generalized Kirkwood
and Poisson-Boltzmann solvation models

**ESOLV2 Subroutine**

"esolv2" calculates second derivatives of the implicit
solvation energy for surface area, generalized Born,
generalized Kirkwood and Poisson-Boltzmann solvation models

**ESOLV2A Subroutine**

"esolv2a" calculates second derivatives of the implicit solvation
potential energy by finite differences

**ESOLV2B Subroutine**

"esolv2b" finds implicit solvation gradients needed for
calculation of the Hessian matrix by finite differences

**ESOLV3 Subroutine**

"esolv3" calculates the implicit solvation energy for
surface area, generalized Born, generalized Kirkwood
and Poisson-Boltzmann solvation models; also partitions
the energy among the atoms

**ESTRBND Subroutine**

"estrbnd" calculates the stretch-bend potential energy

**ESTRBND1 Subroutine**

"estrbnd1" calculates the stretch-bend potential energy and
first derivatives with respect to Cartesian coordinates

**ESTRBND2 Subroutine**

"estrbnd2" calculates the stretch-bend potential energy
second derivatives with respect to Cartesian coordinates

**ESTRBND3 Subroutine**

"estrbnd3" calculates the stretch-bend potential energy;
also partitions the energy among the atoms

**ESTRTOR Subroutine**

"estrtor" calculates the stretch-torsion potential energy

**ESTRTOR1 Subroutine**

"estrtor1" calculates the stretch-torsion energy and first
derivatives with respect to Cartesian coordinates

**ESTRTOR2 Subroutine**

"estrtor2" calculates the stretch-torsion potential energy
second derivatives with respect to Cartesian coordinates

**ESTRTOR3 Subroutine**

"estrtor3" calculates the stretch-torsion potential energy;
also partitions the energy terms among the atoms

**ETORS Subroutine**

"etors" calculates the torsional potential energy

**ETORS0A Subroutine**

"etors0a" calculates the torsional potential energy
using a standard sum of Fourier terms

**ETORS0B Subroutine**

"etors0b" calculates the torsional potential energy
for use with potential energy smoothing methods

**ETORS1 Subroutine**

"etors1" calculates the torsional potential energy and first
derivatives with respect to Cartesian coordinates

**ETORS1A Subroutine**

"etors1a" calculates the torsional potential energy and first
derivatives with respect to Cartesian coordinates using a
standard sum of Fourier terms

**ETORS1B Subroutine**

"etors1b" calculates the torsional potential energy and first
derivatives with respect to Cartesian coordinates for use with
potential energy smoothing methods

**ETORS2 Subroutine**

"etors2" calculates the second derivatives of the torsional
energy for a single atom

**ETORS2A Subroutine**

"etors2a" calculates the second derivatives of the torsional
energy for a single atom using a standard sum of Fourier terms

**ETORS2B Subroutine**

"etors2b" calculates the second derivatives of the torsional
energy for a single atom for use with potential energy
smoothing methods

**ETORS3 Subroutine**

"etors3" calculates the torsional potential energy; also
partitions the energy among the atoms

**ETORS3A Subroutine**

"etors3a" calculates the torsional potential energy using
a standard sum of Fourier terms and partitions the energy
among the atoms

**ETORS3B Subroutine**

"etors3b" calculates the torsional potential energy for use
with potential energy smoothing methods and partitions the
energy among the atoms

**ETORTOR Subroutine**

"etortor" calculates the torsion-torsion potential energy

**ETORTOR1 Subroutine**

"etortor1" calculates the torsion-torsion energy and first
derivatives with respect to Cartesian coordinates

**ETORTOR2 Subroutine**

"etortor2" calculates the torsion-torsion potential energy
second derivatives with respect to Cartesian coordinates

**ETORTOR3 Subroutine**

"etortor3" calculates the torsion-torsion potential energy;
also partitions the energy terms among the atoms

**EUREY Subroutine**

"eurey" calculates the Urey-Bradley 1-3 interaction energy

**EUREY1 Subroutine**

"eurey1" calculates the Urey-Bradley interaction energy and
its first derivatives with respect to Cartesian coordinates

**EUREY2 Subroutine**

"eurey2" calculates second derivatives of the Urey-Bradley
interaction energy for a single atom at a time

**EUREY3 Subroutine**

"eurey3" calculates the Urey-Bradley energy; also
partitions the energy among the atoms

**EVCORR Subroutine**

"evcorr" computes the long range van der Waals correction
to the energy via numerical integration

**EVCORR1 Subroutine**

"evcorr1" computes the long range van der Waals correction
to the energy and virial via numerical integration

**EWALDCOF Subroutine**

"ewaldcof" finds an Ewald coefficient such that all terms
beyond the specified cutoff distance will have a value less
than a specified tolerance

**EWCA Subroutine**

"ewca" find the Weeks-Chandler-Andersen dispersion energy
of a solute using an HCT-like method

**EWCA1 Subroutine**

"ewca1" finds the Weeks-Chandler-Anderson dispersion energy
and derivatives of a solute

**EWCA3 Subroutine**

"ewca3" find the Weeks-Chandler-Andersen dispersion energy
of a solute; also partitions the energy among the atoms

**EWCA3X Subroutine**

"ewca3x" finds the Weeks-Chandler-Anderson dispersion energy
of a solute using a numerical "onion shell" method; also
partitions the energy among the atoms

**EWCAX Subroutine**

"ewcax" finds the Weeks-Chandler-Anderson dispersion energy
of a solute using a numerical "onion shell" method

**EXPLORE Subroutine**

"explore" uses simulated annealing on an initial crude
embedded distance geoemtry structure to refine versus the
bound, chirality, planarity and torsional error functions

**EXTENT Subroutine**

"extent" finds the largest interatomic distance in a system

**EXTRA Subroutine**

"extra" calculates any additional user defined potential
energy contribution

**EXTRA1 Subroutine**

"extra1" calculates any additional user defined potential
energy contribution and its first derivatives

**EXTRA2 Subroutine**

"extra2" calculates second derivatives of any additional
user defined potential energy contribution for a single
atom at a time

**EXTRA3 Subroutine**

"extra3" calculates any additional user defined potential
contribution and also partitions the energy among the atoms

**FATAL Subroutine**

"fatal" terminates execution due to a user request, a severe
error or some other nonstandard condition

**FFTBACK Subroutine**

"fftback" performs a 3-D FFT backward transform via a single
3-D transform or three separate 1-D transforms

**FFTCLOSE Subroutine**

"fftclose" does cleanup after performing a 3-D FFT by destroying
the FFTW plans for the forward and backward transforms

**FFTFRONT Subroutine**

"fftfront" performs a 3-D FFT forward transform via a single
3-D transform or three separate 1-D transforms

**FFTSETUP Subroutine**

"fftsetup" does initialization for a 3-D FFT to be computed
via either the FFTPACK or FFTW libraries

**FIELD Subroutine**

"field" sets the force field potential energy functions from
a parameter file and modifications specified in a keyfile

**FINAL Subroutine**

"final" performs any final program actions such as deallocation
of global memory, prints a status message, and then pauses if
necessary to avoid closing the execution window

**FINDATM Subroutine**

"findatm" locates a specific PDB atom name type within a
range of atoms from the PDB file, returns zero if the name
type was not found

**FITRSD Subroutine**

"fitrsd" computes residuals for electrostatic potential fitting
including total charge restraints, dipole and quadrupole moment
targets, and restraints to initial parameter values

**FITTORS Subroutine**

"fittors" refines torsion parameters based on a quantum
mechanical optimized energy surface

**FIXFRAME Subroutine**

"fixframe" is a service routine that alters the local frame
definition for specified atoms

**FIXPDB Subroutine**

"fixpdb" corrects problems with PDB files by converting residue
and atom names to the standard forms used by Tinker

**FIXPOLE Subroutine**

"fixpole" performs unit conversion of the multipole components,
rounds moments to desired precision, and enforces integer net
charge and traceless quadrupoles

**FLATTEN Subroutine**

"flatten" sets the type of smoothing method and the extent of
surface deformation for use with potential energy smoothing

**FPHI_MPOLE Subroutine**

"fphi_mpole" extracts the permanent multipole potential from
the particle mesh Ewald grid

**FPHI_TO_CPHI Subroutine**

"fphi_to_cphi" transforms the reciprocal space potential from
fractional to Cartesian coordinates

**FPHI_UIND Subroutine**

"fphi_uind" extracts the induced dipole potential from
the particle mesh Ewald grid

**FRACDIST Subroutine**

"fracdist" computes a normalized distribution of the pairwise
fractional distances between the smoothed upper and lower bounds

**FRAC_TO_CART Subroutine**

"frac_to_cart" computes a transformation matrix to convert
a multipole object in fraction coordinates to Cartesian

**FRAME13 Subroutine**

"frame13" finds local coordinate frame defining atoms in cases
where the use of 1-3 connected atoms is required

**FREEUNIT Function**

"freeunit" finds an unopened Fortran I/O unit and returns
its numerical value from 1 to 99; the units already assigned
to "input" and "iout" (usually 5 and 6) are skipped since
they have special meaning as the default I/O units

**GAMMLN Function**

"gammln" uses a series expansion due to Lanczos to compute
the natural logarithm of the Gamma function at "x" in [0,1]

**GAUSSJORDAN Subroutine**

"gaussjordan" solves a system of linear equations by using
the method of Gaussian elimination with partial pivoting

**GDA Program**

"gda" implements Gaussian Density Annealing (GDA) algorithm
for global optimization via simulated annealing

**GDA1 Subroutine**

**GDA2 Function**

**GDA3 Subroutine**

**GDASTAT Subroutine**

for a GDA integration step; also saves the coordinates

**GENDOT Subroutine**

"gendot" finds the coordinates of a specified number of surface
points for a sphere with the input radius and coordinate center

**GEODESIC Subroutine**

"geodesic" smooths the upper and lower distance bounds via
the triangle inequality using a sparse matrix version of a
shortest path algorithm

**GEOMETRY Function**

"geometry" finds the value of the interatomic distance, angle
or dihedral angle defined by two to four input atoms

**GETARC Subroutine**

"getarc" asks for a coordinate archive or trajectory file name,
then reads in the initial set of coordinates

**GETBASE Subroutine**

"getbase" finds the base heavy atoms for a single nucleotide
residue and copies the names and coordinates to the Protein
Data Bank file

**GETCHUNK Subroutine**

"getchunk" determines the number of grid point "chunks" used
along each axis of the PME grid for parallelization

**GETINT Subroutine**

"getint" asks for an internal coordinate file name, then reads
the internal coordinates and computes Cartesian coordinates

**GETKEY Subroutine**

"getkey" finds a valid keyfile and stores its contents as
line images for subsequent keyword parameter searching

**GETMOL Subroutine**

"getmol" asks for a MDL MOL molecule file name,
then reads the coordinates from the file

**GETMOL2 Subroutine**

"getmol2" asks for a Tripos MOL2 molecule file name,
then reads the coordinates from the file

**GETMONITOR Subroutine**

**GETNUCH Subroutine**

"getnuch" finds the nucleotide hydrogen atoms for a single
residue and copies the names and coordinates to the Protein
Data Bank file

**GETNUMB Subroutine**

"getnumb" searches an input string from left to right for an
integer and puts the numeric value in "number"; returns zero
with "next" unchanged if no integer value is found

**GETPDB Subroutine**

"getpdb" asks for a Protein Data Bank file name,
then reads in the coordinates file

**GETPRB Subroutine**

"getprb" tests for a possible probe position at the interface
between three neighboring atoms

**GETPRM Subroutine**

"getprm" finds the potential energy parameter file
and then opens and reads the parameters

**GETPROH Subroutine**

"getproh" finds the hydrogen atoms for a single amino acid
residue and copies the names and coordinates to the Protein
Data Bank file

**GETREF Subroutine**

"getref" copies structure information from the reference area
into the standard variables for the current system structure

**GETSEQ Subroutine**

"getseq" asks the user for the amino acid sequence
and torsional angle values needed to define a peptide

**GETSEQN Subroutine**

"getseqn" asks the user for the nucleotide sequence and
torsional angle values needed to define a nucleic acid

**GETSIDE Subroutine**

"getside" finds the side chain heavy atoms for a single amino
acid residue and copies the names and coordinates to the Protein
Data Bank file

**GETSTRING Subroutine**

"getstring" searches for a quoted text string within an input
character string; the region between the first and second
double quote is returned as the "text"; if the actual text is
too long, only the first part is returned

**GETTEXT Subroutine**

"gettext" searches an input string for the first string of
non-blank characters; the region from a non-blank character
to the first space or tab is returned as "text"; if the
actual text is too long, only the first part is returned

**GETTIME Subroutine**

"gettime" finds the elapsed wall clock and CPU times in seconds
since the last call to "settime"

**GETTOR Subroutine**

"gettor" tests for a possible torus position at the interface
between two atoms, and finds the torus radius, center and axis

**GETWORD Subroutine**

"getword" searches an input string for the first alphabetic
character (A-Z or a-z); the region from this first character
to the first blank space or separator is returned as a "word";
if the actual word is too long, only the first part is returned

**GETXYZ Subroutine**

"getxyz" asks for a Cartesian coordinate file name,
then reads in the coordinates file

**GHMCSTEP Subroutine**

"ghmcstep" performs a single stochastic dynamics time step via
the generalized hybrid Monte Carlo (GHMC) algorithm to ensure
exact sampling from the Boltzmann density

**GHMCTERM Subroutine**

"ghmcterm" finds the friction and fluctuation terms needed
to update velocities during GHMC stochastic dynamics

**GRADFAST Subroutine**

"gradfast" calculates the potential energy and first derivatives
for the fast-evolving local valence potential energy terms

**GRADIENT Subroutine**

"gradient" calls subroutines to calculate the potential energy
and first derivatives with respect to Cartesian coordinates

**GRADRGD Subroutine**

"gradrgd" calls subroutines to calculate the potential energy
and first derivatives with respect to rigid body coordinates

**GRADROT Subroutine**

"gradrot" calls subroutines to calculate the potential
energy and its torsional first derivatives

**GRADSLOW Subroutine**

"gradslow" calculates the potential energy and first derivatives
for the slow-evolving nonbonded potential energy terms

**GRAFIC Subroutine**

"grafic" outputs the upper & lower triangles and diagonal
of a square matrix in a schematic form for visual inspection

**GRID_DISP Subroutine**

"grid_disp" places the damped dispersion coefficients onto
the particle mesh Ewald grid

**GRID_MPOLE Subroutine**

"grid_mpole" places the fractional atomic multipoles onto
the particle mesh Ewald grid

**GRID_PCHG Subroutine**

"grid_pchg" places the fractional atomic partial charges onto
the particle mesh Ewald grid

**GRID_UIND Subroutine**

"grid_uind" places the fractional induced dipoles onto the
particle mesh Ewald grid

**GROUPS Subroutine**

"groups" tests a set of atoms to see if all are members of a
single atom group or a pair of atom groups; if so, then the
correct intra- or intergroup weight is assigned

**GRPLINE Subroutine**

"grpline" tests each atom group for linearity of the sites
contained in the group

**GSORT Subroutine**

"gsort" uses the Gram-Schmidt algorithm to build orthogonal
vectors for sliding block interative matrix diagonalization

**GYRATE Subroutine**

"gyrate" computes the radius of gyration of a molecular system
from its atomic coordinates; only active atoms are included

**HANGLE Subroutine**

"hangle" constructs hybrid angle bending parameters given
an initial state, final state and "lambda" value

**HATOM Subroutine**

"hatom" assigns a new atom type to each hybrid site

**HBOND Subroutine**

"hbond" constructs hybrid bond stretch parameters given
an initial state, final state and "lambda" value

**HCHARGE Subroutine**

"hcharge" constructs hybrid charge interaction parameters
given an initial state, final state and "lambda" value

**HDIPOLE Subroutine**

"hdipole" constructs hybrid dipole interaction parameters
given an initial state, final state and "lambda" value

**HESSBLK Subroutine**

"hessblk" calls subroutines to calculate the Hessian elements
for each atom in turn with respect to Cartesian coordinates

**HESSIAN Subroutine**

"hessian" calls subroutines to calculate the Hessian elements
for each atom in turn with respect to Cartesian coordinates

**HESSRGD Subroutine**

"hessrgd" computes the numerical Hessian elements with
respect to rigid body coordinates via 6*ngroup+1 gradient
evaluations

**HESSROT Subroutine**

"hessrot" computes numerical Hessian elements with respect
to torsional angles; either the diagonal or the full matrix
can be calculated; the full matrix needs nomega+1 gradient
evaluations while the diagonal needs just two evaluations

**HETATOM Subroutine**

"hetatom" translates water molecules and ions in Protein Data
Bank format to a Cartesian coordinate file and sequence file

**HIMPTOR Subroutine**

"himptor" constructs hybrid improper torsional parameters
given an initial state, final state and "lambda" value

**HOOVER Subroutine**

"hoover" applies a combined thermostat and barostat via a
Nose-Hoover chain algorithm

**HSTRBND Subroutine**

"hstrbnd" constructs hybrid stretch-bend parameters given
an initial state, final state and "lambda" value

**HSTRTOR Subroutine**

"hstrtor" constructs hybrid stretch-torsion parameters
given an initial state, final state and "lambda" value

**HTORS Subroutine**

"htors" constructs hybrid torsional parameters for a given
initial state, final state and "lambda" value

**HVDW Subroutine**

"hvdw" constructs hybrid van der Waals  parameters given
an initial state, final state and "lambda" value

**HYBRID Subroutine**

"hybrid" constructs the hybrid hamiltonian for a specified
initial state, final state and mutation parameter "lambda"

**IJKPTS Subroutine**

"ijkpts" stores a set of indices used during calculation
of macroscopic reaction field energetics

**IMAGE Subroutine**

"image" takes the components of pairwise distance between
two points in a periodic box and converts to the components
of the minimum image distance

**IMAGEN Subroutine**

"imagen" takes the components of pairwise distance between
two points and converts to the components of the minimum
image distance

**IMAGER Subroutine**

"imager" takes the components of pairwise distance between
two points in the same or neighboring periodic boxes and
converts to the components of the minimum image distance

**IMPOSE Subroutine**

"impose" performs the least squares best superposition
of two atomic coordinate sets via a quaternion method;
upon return, the first coordinate set is unchanged while
the second set is translated and rotated to give best fit;
the final root mean square fit is returned in "rmsvalue"

**INDTCGA Subroutine**

"indtcga" computes the induced dipoles and intermediates used
in polarization force calculation for the TCG method with dp
cross terms = true, initial guess mu0 = 0 and using a diagonal
preconditioner

**INDTCGB Subroutine**

"indtcgb" computes the induced dipoles and intermediates used
in polarization force calculation for the TCG method with dp
cross terms = true, initial guess mu0 = direct and using diagonal
preconditioner

**INDUCE Subroutine**

"induce" computes the induced dipole moments at polarizable
sites due to direct or mutual polarization

**INDUCE0A Subroutine**

"induce0a" computes the induced dipole moments at polarizable
sites using a preconditioned conjugate gradient solver

**INDUCE0B Subroutine**

"induce0b" computes and stores the induced dipoles via
the truncated conjugate gradient (TCG) method

**INDUCE0C Subroutine**

"induce0c" computes the induced dipole moments at polarizable
sites for generalized Kirkwood SCRF and vacuum environments

**INDUCE0D Subroutine**

"induce0d" computes the induced dipole moments at polarizable
sites for Poisson-Boltzmann SCRF and vacuum environments

**INEDGE Subroutine**

"inedge" inserts a concave edge into the
linked list for its temporary torus

**INERTIA Subroutine**

"inertia" computes the principal moments of inertia for the
system, and optionally translates the center of mass to the
origin and rotates the principal axes onto the global axes

**INITATOM Subroutine**

"initatom" sets the atomic symbol, standard atomic weight,
van der Waals radius and covalent radius for each element in
the periodic table

**INITERR Function**

"initerr" is the initial error function and derivatives for
a distance geometry embedding; it includes components from
the local geometry and torsional restraint errors

**INITIAL Subroutine**

"initial" sets up original values for some parameters and
variables that might not otherwise get initialized

**INITMMFF Subroutine**

"initmmff" initializes some parameter values for the Merck
Molecular force field

**INITPRM Subroutine**

"initprm" completely initializes a force field by setting all
parameters to zero and using defaults for control values

**INITRES Subroutine**

"initres" sets biopolymer residue names and biotype codes used
in PDB file conversion and automated generation of structures

**INITROT Subroutine**

"initrot" sets the torsional angles which are to be rotated
in subsequent computation, by default automatically selects
all rotatable single bonds; optionally makes atoms inactive
when they are not moved by any torsional rotation

**INSERT Subroutine**

"insert" adds the specified atom to the Cartesian
coordinates list and shifts the remaining atoms

**INTEDIT Program**

"intedit" allows the user to extract information from
or alter the values within an internal coordinates file

**INTERPOL Subroutine**

"interpol" computes intergroup induced dipole moments for use
during removal of intergroup polarization

**INTXYZ Program**

"intxyz" takes as input an internal coordinates file,
converts to and then writes out Cartesian coordinates

**INVBETA Function**

"invbeta" computes the inverse Beta distribution function
via a combination of Newton iteration and bisection search

**INVERT Subroutine**

"invert" inverts a matrix using the Gauss-Jordan method

**IPEDGE Subroutine**

"ipedge" inserts convex edge into linked list for atom

**JACOBI Subroutine**

"jacobi" performs a matrix diagonalization of a real
symmetric matrix by the method of Jacobi rotations

**JUSTIFY Subroutine**

"justify" converts a text string to right justified format
with leading blank spaces

**KANGANG Subroutine**

"kangang" assigns the parameters for angle-angle cross term
interactions and processes new or changed parameter values

**KANGLE Subroutine**

"kangle" assigns the force constants and ideal angles for
the bond angles; also processes new or changed parameters

**KANGLEM Subroutine**

"kanglem" assigns the force constants and ideal angles for
bond angles according to the Merck Molecular Force Field (MMFF)

**KANGTOR Subroutine**

"kangtor" assigns parameters for angle-torsion interactions
and processes new or changed parameter values

**KATOM Subroutine**

"katom" assigns an atom type definitions to each atom in
the structure and processes any new or changed values

**KBOND Subroutine**

"kbond" assigns a force constant and ideal bond length
to each bond in the structure and processes any new or
changed parameter values

**KBONDM Subroutine**

"kbondm" assigns a force constant and ideal bond length to
each bond according to the Merck Molecular Force Field (MMFF)

**KCHARGE Subroutine**

"kcharge" assigns partial charges to the atoms within
the structure and processes any new or changed values

**KCHARGEM Subroutine**

"kchargem" assigns partial charges to the atoms according to
the Merck Molecular Force Field (MMFF)

**KCHGFLX Subroutine**

"kchgflx" assigns a force constant and ideal bond length
to each bond in the structure and processes any new or
changed parameter values

**KCHGTRN Subroutine**

"kchgtrn" assigns charge magnitude and damping parameters for
charge transfer interactions and processes any new or changed
values for these parameters

**KCHIRAL Subroutine**

"kchiral" determines the target value for each chirality
and planarity restraint as the signed volume of the
parallelpiped spanned by vectors from a common atom to
each of three other atoms

**KDIPOLE Subroutine**

"kdipole" assigns bond dipoles to the bonds within
the structure and processes any new or changed values

**KDISP Subroutine**

"kdisp" assigns C6 coefficients and damping parameters for
dispersion interactions and processes any new or changed
values for these parameters

**KENEG Subroutine**

"keneg" applies primary and secondary electronegativity bond
length corrections to applicable bond parameters

**KEWALD Subroutine**

"kewald" assigns particle mesh Ewald parameters and options
for a periodic system

**KEXTRA Subroutine**

"kextra" assigns parameters to any additional user defined
potential energy contribution

**KGB Subroutine**

"kgb" initializes parameters needed for the generalized
Born implicit solvation models

**KGEOM Subroutine**

"kgeom" asisgns parameters for geometric restraint terms
to be included in the potential energy calculation

**KGK Subroutine**

"kgk" initializes parameters needed for the generalized
Kirkwood implicit solvation model

**KHPMF Subroutine**

"khpmf" initializes parameters needed for the hydrophobic
potential of mean force nonpolar implicit solvation model

**KIMPROP Subroutine**

"kimprop" assigns potential parameters to each improper
dihedral in the structure and processes any changed values

**KIMPTOR Subroutine**

"kimptor" assigns torsional parameters to each improper
torsion in the structure and processes any changed values

**KINAUX Subroutine**

"kinaux" computes the total kinetic energy and temperature
for auxiliary dipole variables used in iEL polarization

**KINETIC Subroutine**

"kinetic" computes the total kinetic energy and kinetic energy
contributions to the pressure tensor by summing over velocities

**KMETAL Subroutine**

"kmetal" assigns ligand field parameters to transition metal
atoms and processes any new or changed parameter values

**KMPOLE Subroutine**

"kmpole" assigns atomic multipole moments to the atoms of
the structure and processes any new or changed values

**KNP Subroutine**

"knp" initializes parameters needed for the cavity-plus-
dispersion nonpolar implicit solvation model

**KONVEC Subroutine**

"konvec" finds a Hessian-vector product via finite-difference
evaluation of the gradient based on atomic displacements

**KOPBEND Subroutine**

"kopbend" assigns the force constants for out-of-plane bends
at trigonal centers via Wilson-Decius-Cross or Allinger angles;
also processes any new or changed parameter values

**KOPBENDM Subroutine**

"kopbendm" assigns the force constants for out-of-plane bends
according to the Merck Molecular Force Field (MMFF)

**KOPDIST Subroutine**

"kopdist" assigns the force constants for out-of-plane
distance at trigonal centers via the central atom height;
also processes any new or changed parameter values

**KORBIT Subroutine**

"korbit" assigns pi-orbital parameters to conjugated systems
and processes any new or changed parameters

**KPB Subroutine**

"kpb" assigns parameters needed for the Poisson-Boltzmann
implicit solvation model implemented via APBS

**KPITORS Subroutine**

"kpitors" assigns pi-system torsion parameters to torsions
needing them, and processes any new or changed values

**KPOLAR Subroutine**

"kpolar" assigns atomic dipole polarizabilities to the atoms
within the structure and processes any new or changed values

**KREPEL Subroutine**

"krepel" assigns the size values, exponential parameter and
number of valence electrons for Pauli repulsion interactions
and processes any new or changed values for these parameters

**KSA Subroutine**

"ksa" initializes parameters needed for surface area-based
implicit solvation models including ASP and SASA

**KSOLV Subroutine**

"ksolv" assigns implicit solvation energy parameters for
the surface area, generalized Born, generalized Kirkwood,
Poisson-Boltzmann, cavity-dispersion and HPMF models

**KSTRBND Subroutine**

"kstrbnd" assigns parameters for stretch-bend interactions
and processes new or changed parameter values

**KSTRBNDM Subroutine**

"kstrbndm" assigns parameters for stretch-bend interactions
according to the Merck Molecular Force Field (MMFF)

**KSTRTOR Subroutine**

"kstrtor" assigns stretch-torsion parameters to torsions
needing them, and processes any new or changed values

**KTORS Subroutine**

"ktors" assigns torsional parameters to each torsion in
the structure and processes any new or changed values

**KTORSM Subroutine**

"ktorsm" assigns torsional parameters to each torsion according
to the Merck Molecular Force Field (MMFF)

**KTORTOR Subroutine**

"ktortor" assigns torsion-torsion parameters to adjacent
torsion pairs and processes any new or changed values

**KUREY Subroutine**

"kurey" assigns the force constants and ideal distances
for the Urey-Bradley 1-3 interactions; also processes any
new or changed parameter values

**KVDW Subroutine**

"kvdw" assigns the parameters to be used in computing the
van der Waals interactions and processes any new or changed
values for these parameters

**LATTICE Subroutine**

"lattice" stores the periodic box dimensions and sets angle
values to be used in computing fractional coordinates

**LBFGS Subroutine**

"lbfgs" is a limited memory BFGS quasi-newton nonlinear
optimization routine

**LIGASE Subroutine**

"ligase" translates a nucleic acid structure in Protein Data
Bank format to a Cartesian coordinate file and sequence file

**LIGHTS Subroutine**

"lights" computes the set of nearest neighbor interactions
using the method of lights algorithm

**LINBODY Subroutine**

"linbody" finds the angular velocity of a linear rigid body
given the inertia tensor and angular momentum

**LMSTEP Subroutine**

"lmstep" computes a Levenberg-Marquardt step during a nonlinear
least squares calculation using ideas from the MINPACK LMPAR
routine and the internal doubling strategy of Dennis and Schnabel

**LOCALMIN Subroutine**

"localmin" is used during normal mode local search to
perform a Cartesian coordinate energy minimization

**LOCALRGD Subroutine**

"localrgd" is used during the PSS local search procedure
to perform a rigid body energy minimization

**LOCALROT Subroutine**

"localrot" is used during the PSS local search procedure
to perform a torsional space energy minimization

**LOCALXYZ Subroutine**

"localxyz" is used during the potential smoothing and search
procedure to perform a local optimization at the current
smoothing level

**LOCERR Function**

"locerr" is the local geometry error function and derivatives
including the 1-2, 1-3 and 1-4 distance bound restraints

**LOWCASE Subroutine**

"lowcase" converts a text string to all lower case letters

**MAJORIZE Subroutine**

"majorize" refines the projected coordinates by attempting to
minimize the least square residual between the trial distance
matrix and the distances computed from the coordinates

**MAKEBAR Subroutine**

**MAKEBOX Subroutine**

"makebox" builds a periodic box of a desired size by randomly
copying a specified number of monomers into a target box size,
followed by optional excluded volume refinement

**MAKEINT Subroutine**

"makeint" converts Cartesian to internal coordinates where
selection of internal coordinates is controlled by "mode"

**MAKEPDB Subroutine**

"makepdb" cconstructs a Protein Data Bank file from a set
of Cartesian coordinates with special handling for systems
consisting of biopolymer chains, ligands and water molecules

**MAKEREF Subroutine**

"makeref" copies the information contained in the "xyz" file
of the current structure into corresponding reference areas

**MAKEXYZ Subroutine**

"makexyz" generates a complete set of Cartesian coordinates
for a full structure from the internal coordinate values

**MAPCHECK Subroutine**

"mapcheck" checks the current minimum energy structure
for possible addition to the master list of local minima

**MATCH1 Subroutine**

"match1" finds and stores the first multipole component found
on a line of output from Stone's GDMA program

**MATCH2 Subroutine**

"match2" finds and stores the second multipole component found
on a line of output from Stone's GDMA program

**MATCH3 Subroutine**

"match3" finds and stores the third multipole component found
on a line of output from Stone's GDMA program

**MAXWELL Function**

"maxwell" returns a speed in Angstroms/picosecond randomly
selected from a 3-D Maxwell-Boltzmann distribution for the
specified particle mass and system temperature

**MBUILD Subroutine**

"mbuild" performs a complete rebuild of the atomic multipole
electrostatic neighbor list for all sites

**MCM1 Function**

"mcm1" is a service routine that computes the energy and
gradient for truncated Newton optimization in Cartesian
coordinate space

**MCM2 Subroutine**

"mcm2" is a service routine that computes the sparse matrix
Hessian elements for truncated Newton optimization in Cartesian
coordinate space

**MCMSTEP Function**

"mcmstep" implements the minimization phase of an MCM step
via Cartesian minimization following a Monte Carlo step

**MDINIT Subroutine**

"mdinit" initializes the velocities and accelerations
for a molecular dynamics trajectory, including restarts

**MDREST Subroutine**

"mdrest" finds and removes any translational or rotational
kinetic energy of the overall system center of mass

**MDSAVE Subroutine**

"mdsave" writes molecular dynamics trajectory snapshots and
auxiliary files with velocity, force or induced dipole data;
also checks for user requested termination of a simulation

**MDSTAT Subroutine**

"mdstat" is called at each molecular dynamics time step to
form statistics on various average values and fluctuations,
and to periodically save the state of the trajectory

**MEASFN Subroutine**

**MEASFQ Subroutine**

**MEASFS Subroutine**

**MEASPM Subroutine**

"measpm" computes the volume of a single prism section of
the full interior polyhedron

**MECHANIC Subroutine**

"mechanic" sets up needed parameters for the potential energy
calculation and reads in many of the user selectable options

**MERGE Subroutine**

"merge" combines the reference and current structures into
a single new "current" structure containing the reference
atoms followed by the atoms of the current structure

**METRIC Subroutine**

"metric" takes as input the trial distance matrix and computes
the metric matrix of all possible dot products between the atomic
vectors and the center of mass using the law of cosines and the
following formula for the distances to the center of mass:

**MIDERR Function**

"miderr" is the secondary error function and derivatives
for a distance geometry embedding; it includes components
from the distance bounds, local geometry, chirality and
torsional restraint errors

**MINIMIZ1 Function**

"minimiz1" is a service routine that computes the energy and
gradient for a low storage BFGS optimization in Cartesian
coordinate space

**MINIMIZE Program**

"minimize" performs energy minimization in Cartesian coordinate
space using a low storage BFGS nonlinear optimization

**MINIROT Program**

"minirot" performs an energy minimization in torsional
angle space using a low storage BFGS nonlinear optimization

**MINIROT1 Function**

"minirot1" is a service routine that computes the energy
and gradient for a low storage BFGS nonlinear optimization
in torsional angle space

**MINPATH Subroutine**

"minpath" is a routine for finding the triangle smoothed upper
and lower bounds of each atom to a specified root atom using a
sparse variant of the Bellman-Ford shortest path algorithm

**MINRIGID Program**

"minrigid" performs an energy minimization of rigid body atom
groups using a low storage BFGS nonlinear optimization

**MINRIGID1 Function**

"minrigid1" is a service routine that computes the energy
and gradient for a low storage BFGS nonlinear optimization
of rigid bodies

**MLIGHT Subroutine**

"mlight" performs a complete rebuild of the atomic multipole
pair neighbor list for all sites using the method of lights

**MLIST Subroutine**

"mlist" performs an update or a complete rebuild of the
nonbonded neighbor lists for atomic multipoles

**MMID Subroutine**

"mmid" implements a modified midpoint method to advance the
integration of a set of first order differential equations

**MODECART Subroutine**

**MODERGD Subroutine**

**MODEROT Subroutine**

**MODESRCH Subroutine**

**MODETORS Subroutine**

**MODULI Subroutine**

"moduli" sets the moduli of the inverse discrete Fourier
transform of the B-splines

**MOL2XYZ Program**

"mol2xyz" takes as input a Tripos MOL2 coordinates file,
converts to and then writes out Cartesian coordinates

**MOLECULE Subroutine**

"molecule" counts the molecules, assigns each atom to
its molecule and computes the mass of each molecule

**MOLMERGE Subroutine**

"molmerge" connects fragments and removes duplicate atoms
during generation of a unit cell from an asymmetric unit

**MOLSETUP Subroutine**

"molsetup" generates trial parameters needed to perform
polarizable multipole calculations on a structure read
from distributed multipole analysis output

**MOLUIND Subroutine**

"moluind" computes the molecular induced dipole components
in the presence of an external electric field

**MOLXYZ Program**

"molxyz" takes as input a MDL MOL coordinates file,
converts to and then writes out Cartesian coordinates

**MOMENTS Subroutine**

"moments" computes the total electric charge, dipole and
quadrupole moments for the active atoms as a sum over the
partial charges, bond dipoles and atomic multipole moments

**MOMFULL Subroutine**

"momfull" computes the electric moments for the full system
as a sum over the partial charges, bond dipoles and atomic
multipole moments

**MOMYZE Subroutine**

"momyze" finds and prints the total charge, dipole moment
components, radius of gyration and moments of inertia

**MONTE Program**

"monte" performs a Monte Carlo-Minimization conformational
search using Cartesian single atom or torsional move sets

**MUTATE Subroutine**

"mutate" constructs the hybrid hamiltonian for a specified
initial state, final state and mutation parameter "lambda"

**NBLIST Subroutine**

"nblist" builds and maintains nonbonded pair neighbor lists
for vdw, dispersion, electrostatic and polarization terms

**NEARBY Subroutine**

"nearby" finds all of the through-space neighbors of each
atom for use in surface area and volume calculations

**NEEDUPDATE Subroutine**

**NEWATM Subroutine**

"newatm" creates and defines an atom needed for the
Cartesian coordinates file, but which may not present
in the original Protein Data Bank file

**NEWTON Program**

"newton" performs an energy minimization in Cartesian
coordinate space using a truncated Newton method

**NEWTON1 Function**

"newton1" is a service routine that computes the energy
and gradient for truncated Newton optimization in Cartesian
coordinate space

**NEWTON2 Subroutine**

"newton2" is a service routine that computes the sparse
matrix Hessian elements for truncated Newton optimization
in Cartesian coordinate space

**NEWTROT Program**

"newtrot" performs an energy minimization in torsional angle
space using a truncated Newton conjugate gradient method

**NEWTROT1 Function**

"newtrot1" is a service routine that computes the energy
and gradient for truncated Newton conjugate gradient
optimization in torsional angle space

**NEWTROT2 Subroutine**

"newtrot2" is a service routine that computes the sparse
matrix Hessian elements for truncated Newton optimization
in torsional angle space

**NEXTARG Subroutine**

"nextarg" finds the next unused command line argument
and returns it in the input character string

**NEXTTEXT Function**

"nexttext" finds and returns the location of the first
non-blank character within an input text string; zero
is returned if no such character is found

**NORMAL Function**

"normal" generates a random number from a normal Gaussian
distribution with a mean of zero and a variance of one

**NOSE Subroutine**

"nose" performs a single molecular dynamics time step via
a Nose-Hoover extended system isothermal-isobaric algorithm

**NSPLINE Subroutine**

"nspline" computes coefficients for an nonperiodic cubic spline
with natural boundary conditions where the first and last second
derivatives are already known

**NUCBASE Subroutine**

"nucbase" builds the side chain for a single nucleotide base
in terms of internal coordinates

**NUCCHAIN Subroutine**

"nucchain" builds up the internal coordinates for a nucleic
acid sequence from the sugar type, backbone and glycosidic
torsional values

**NUCLEIC Program**

"nucleic" builds the internal and Cartesian coordinates
of a polynucleotide from nucleic acid sequence and torsional
angle values for the nucleic acid backbone and side chains

**NUMBER Function**

"number" converts a text numeral into an integer value;
the input string must contain only numeric characters

**NUMERAL Subroutine**

"numeral" converts an input integer number into the
corresponding right- or left-justified text numeral

**NUMGRAD Subroutine**

"numgrad" computes the gradient of the objective function
"fvalue" with respect to Cartesian coordinates of the atoms
via a one-sided or two-sided numerical differentiation

**OCVM Subroutine**

"ocvm" is an optimally conditioned variable metric nonlinear
optimization routine without line searches

**OLDATM Subroutine**

"oldatm" get the Cartesian coordinates for an atom from
the Protein Data Bank file, then assigns the atom type
and atomic connectivities

**OPBGUESS Function**

"opbguess" sets approximate out-of-plane bend force constants
based on atom type and connected atoms

**OPENEND Subroutine**

"openend" opens a file on a Fortran unit such that the position
is set to the bottom for appending to the end of the file

**OPREP Subroutine**

"oprep" sets up the frictional and random terms needed to
update positions and velocities for the BAOAB integrator

**OPTFIT Function**

**OPTIMIZ1 Function**

"optimiz1" is a service routine that computes the energy and
gradient for optimally conditioned variable metric optimization
in Cartesian coordinate space

**OPTIMIZE Program**

"optimize" performs energy minimization in Cartesian coordinate
space using an optimally conditioned variable metric method

**OPTINIT Subroutine**

"optinit" initializes values and keywords used by multiple
structure optimization methods

**OPTIROT Program**

"optirot" performs an energy minimization in torsional angle
space using an optimally conditioned variable metric method

**OPTIROT1 Function**

"optirot1" is a service routine that computes the energy and
gradient for optimally conditioned variable metric optimization
in torsional angle space

**OPTRIGID Program**

"optrigid" performs an energy minimization of rigid body atom
groups using an optimally conditioned variable metric method

**OPTRIGID1 Function**

"optrigid1" is a service routine that computes the energy
and gradient for optimally conditioned variable metric
optimization of rigid bodies

**OPTSAVE Subroutine**

"optsave" is used by the optimizers to write imtermediate
coordinates and other relevant information; also checks for
user requested termination of an optimization

**ORBITAL Subroutine**

"orbital" finds and organizes lists of atoms in a pisystem,
bonds connecting pisystem atoms and torsions whose central
atoms are both pisystem atoms

**ORIENT Subroutine**

"orient" computes a set of reference Cartesian coordinates
in standard orientation for each rigid body atom group

**ORTHOG Subroutine**

"orthog" performs an orthogonalization of an input matrix
via the modified Gram-Schmidt algorithm

**OVERLAP Subroutine**

"overlap" computes the overlap for two parallel p-orbitals
given the atomic numbers and distance of separation

**PARAMYZE Subroutine**

"paramyze" prints the force field parameters used in the
computation of each of the potential energy terms

**PARTYZE Subroutine**

"partyze" prints the energy component and number of
interactions for each of the potential energy terms

**PASSB Subroutine**

**PASSB2 Subroutine**

**PASSB3 Subroutine**

**PASSB4 Subroutine**

**PASSB5 Subroutine**

**PASSF Subroutine**

**PASSF2 Subroutine**

**PASSF3 Subroutine**

**PASSF4 Subroutine**

**PASSF5 Subroutine**

**PATH Program**

"path" locates a series of structures equally spaced along
a conformational pathway connecting the input reactant and
product structures; a series of constrained optimizations
orthogonal to the path is done via Lagrangian multipliers

**PATH1 Function**

**PATHPNT Subroutine**

"pathpnt" finds a structure on the synchronous transit path
with the specified path value "tpath"

**PATHSCAN Subroutine**

"pathscan" makes a scan of a synchronous transit pathway by
computing structures and energies for specific path values

**PATHVAL Subroutine**

"pathval" computes the synchronous transit path value for
the specified structure

**PAULING Subroutine**

"pauling" uses a rigid body optimization to approximately
pack multiple polypeptide chains

**PAULING1 Function**

"pauling1" is a service routine that computes the energy
and gradient for optimally conditioned variable metric
optimization of rigid bodies

**PBDIRECTPOLFORCE Subroutine**

**PBEMPOLE Subroutine**

"pbempole" calculates the permanent multipole PB energy,
field, forces and torques

**PBMUTUALPOLFORCE Subroutine**

**PDBATOM Subroutine**

"pdbatom" adds an atom to the Protein Data Bank file

**PDBXYZ Program**

"pdbxyz" takes as input a Protein Data Bank file and then
converts to and writes out a Cartesian coordinates file and,
for biopolymers, a sequence file

**PIALTER Subroutine**

"pialter" modifies bond lengths and force constants according
to the "planar" P-P-P bond order values; also alters 2-fold
torsional parameters based on the "nonplanar" bond orders

**PICALC Subroutine**

"picalc" performs a modified Pariser-Parr-Pople molecular
orbital calculation for each conjugated pisystem

**PIMOVE Subroutine**

"pimove" rotates the vector between atoms "list(1)" and
"list(2)" so that atom 1 is at the origin and atom 2 along
the x-axis; the atoms defining the respective planes are
also moved and their bond lengths normalized

**PIPLANE Subroutine**

"piplane" selects the three atoms which specify the plane
perpendicular to each p-orbital; the current version will
fail in certain situations, including ketenes, allenes,
and isolated or adjacent triple bonds

**PISCF Subroutine**

"piscf" performs an SCF molecular orbital calculation for a
pisystem to determine bond orders used in parameter scaling

**PITILT Subroutine**

"pitilt" calculates for each pibond the ratio of the
actual p-orbital overlap integral to the ideal overlap
if the same orbitals were perfectly parallel

**PLACE Subroutine**

"place" finds the probe sites by putting the probe sphere
tangent to each triple of neighboring atoms

**PMONTE Subroutine**

"pmonte" implements a Monte Carlo barostat via random trial
changes in the periodic box volume and shape

**POLARGRP Subroutine**

"polargrp" generates members of the polarization group of
each atom and separate lists of the 1-2, 1-3 and 1-4 group
connectivities

**POLARIZE Program**

"polarize" computes the molecular polarizability by applying
an external field along each axis followed by diagonalization
of the resulting polarizability tensor

**POLEDIT Program**

"poledit" provides for the modification and manipulation
of polarizable atomic multipole electrostatic models

**POLESORT Subroutine**

"polesort" sorts a set of atomic multipole parameters based
on the atom types of centers involved

**POLYMER Subroutine**

"polymer" tests for the presence of an infinite polymer
extending across periodic boundaries

**POLYP Subroutine**

"polyp" is a polynomial product routine that multiplies two
algebraic forms

**POTENTIAL Program**

"potential" calculates the electrostatic potential for a
molecule at a set of grid points; optionally compares to a
target potential or optimizes electrostatic parameters

**POTGRID Subroutine**

"potgrid" generates electrostatic potential grid points in
radially distributed shells based on the molecular surface

**POTNRG Function**

**POTOFF Subroutine**

"potoff" clears the forcefield definition by turning off
the use of each of the potential energy functions

**POTPOINT Subroutine**

"potpoint" calculates the electrostatic potential at a grid
point "i" as the total electrostatic interaction energy of
the system with a positive charge located at the grid point

**POTSTAT Subroutine**

"potstat" computes and prints statistics for the electrostatic
potential over a set of grid points

**POTWRT Subroutine**

**PRECONBLK Subroutine**

"preconblk" applies a preconditioner to an atom block section
of the Hessian matrix

**PRECOND Subroutine**

"precond" solves a simplified version of the Newton equations
Ms = r, and uses the result to precondition linear conjugate
gradient iterations on the full Newton equations in "tnsolve"

**PRESSURE Subroutine**

"pressure" uses the internal virial to find the pressure
in a periodic box and maintains a constant desired pressure
via a barostat method

**PRESSURE2 Subroutine**

"pressure2" applies a box size and velocity correction at
the half time step as needed for the Monte Carlo barostat

**PRIORITY Function**

"priority" decides which of a set of connected atoms should
have highest priority in construction of a local coordinate
frame and returns its atom number; if all atoms are of equal
priority then zero is returned

**PRMEDIT Program**

"prmedit" reformats an existing parameter file, and revises
type and class numbers based on the "atom" parameter ordering

**PRMFORM Subroutine**

"prmform" formats each individual parameter record to conform
to a consistent text layout

**PRMKEY Subroutine**

"prmkey" parses a text string to extract keywords related to
force field potential energy functional forms and constants

**PRMORDER Subroutine**

"prmorder" places a list of atom type or class numbers into
canonical order for potential energy parameter definitions

**PRMSORT Subroutine**

"prmsort" places a list of atom type or class numbers into
canonical order for potential energy parameter definitions

**PRMVAR Subroutine**

"prmvar" determines the optimization values from the
corresponding electrostatic potential energy parameters

**PRMVAR Subroutine**

"prmvar" determines the optimization values from the
corresponding valence potential energy parameters

**PROCHAIN Subroutine**

"prochain" builds up the internal coordinates for an amino
acid sequence from the phi, psi, omega and chi values

**PROJCT Subroutine**

**PROJECT Subroutine**

"project" reads locked vectors from a binary file and projects
them out of the components of the set of trial eigenvectors
using the relation Y = X - U * U^T * X

**PROJECTK Subroutine**

"projectk" reads locked vectors from a binary file and projects
them out of the components of the set of trial eigenvectors
using the relation Y = X - U * U^T * X

**PROMO Subroutine**

"promo" writes a banner message containing information
about the Tinker version, release date and copyright notice

**PROPERTY Function**

"property" takes two input snapshot frames and computes the
value of the property for which the correlation function is
being accumulated

**PROSIDE Subroutine**

"proside" builds the side chain for a single amino acid
residue in terms of internal coordinates

**PROTEIN Program**

"protein" builds the internal and Cartesian coordinates
of a polypeptide from amino acid sequence and torsional
angle values for the peptide backbone and side chains

**PRTARC Subroutine**

"prtarc" writes out a set of Cartesian coordinates for
all active atoms in the Tinker XYZ archive format

**PRTDYN Subroutine**

"prtdyn" writes out the information needed to restart a
molecular dynamics trajectory to an external disk file

**PRTERR Subroutine**

"prterr" writes out a set of coordinates to a disk
file prior to aborting on a serious error

**PRTFIT Subroutine**

"prtfit" makes a key file containing results from fitting a
charge or multipole model to an electrostatic potential grid

**PRTINT Subroutine**

"prtint" writes out a set of Z-matrix internal
coordinates to an external disk file

**PRTMOD Subroutine**

"prtmod" writes out a set of modified Cartesian coordinates
with an optional atom number offset to an external disk file

**PRTMOL2 Program**

"prtmol2" writes out a set of coordinates in Tripos MOL2
format to an external disk file

**PRTPDB Subroutine**

"prtpdb" writes out a set of Protein Data Bank coordinates
to an external disk file

**PRTPOLE Subroutine**

"prtpole" creates a coordinates file, and a key file with
atomic multipoles corrected for intergroup polarization

**PRTPRM Subroutine**

"prtprm" writes out a formatted listing of the default
set of potential energy parameters for a force field

**PRTSEQ Subroutine**

"prtseq" writes out a biopolymer sequence to an external
disk file with 15 residues per line and distinct chains
separated by blank lines

**PRTVAL Subroutine**

"prtval" writes the final valence parameter results to the
standard output and appends the values to a key file

**PRTVIB Subroutine**

"prtvib" writes to an external disk file a series of
coordinate sets representing motion along a vibrational
normal mode

**PRTXYZ Subroutine**

"prtxyz" writes out a set of Cartesian coordinates
to an external disk file

**PSCALE Subroutine**

"pscale" implements a Berendsen barostat by scaling the
coordinates and box dimensions via coupling to an external
constant pressure bath

**PSS Program**

"pss" implements the potential smoothing plus search method
for global optimization in Cartesian coordinate space with
local searches performed in Cartesian or torsional space

**PSS1 Function**

"pss1" is a service routine that computes the energy
and gradient during PSS global optimization in Cartesian
coordinate space

**PSS2 Subroutine**

"pss2" is a service routine that computes the sparse
matrix Hessian elements during PSS global optimization
in Cartesian coordinate space

**PSSRGD1 Function**

"pssrgd1" is a service routine that computes the energy and
gradient during PSS global optimization over rigid bodies

**PSSRIGID Program**

"pssrigid" implements the potential smoothing plus search method
for global optimization for a set of rigid bodies

**PSSROT Program**

"pssrot" implements the potential smoothing plus search method
for global optimization in torsional space

**PSSROT1 Function**

"pssrot1" is a service routine that computes the energy and
gradient during PSS global optimization in torsional space

**PSSWRITE Subroutine**

**PTEST Subroutine**

"ptest" determines the numerical virial tensor, and compares
analytical to numerical values for dE/dV and isotropic pressure

**PTINCY Function**

**PZEXTR Subroutine**

"pzextr" is a polynomial extrapolation routine used during
Bulirsch-Stoer integration of ordinary differential equations

**QIROTMAT Subroutine**

"qirotmat" finds a rotation matrix that describes the
interatomic vector

**QONVEC Subroutine**

"qonvec" is a vector utility routine used during sliding
block iterative matrix diagonalization

**QRFACT Subroutine**

"qrfact" computes the QR factorization of an m by n matrix a
via Householder transformations with optional column pivoting;
the routine determines an orthogonal matrix q, a permutation
matrix p, and an upper trapezoidal matrix r with diagonal
elements of nonincreasing magnitude, such that a*p = q*r; the
Householder transformation for column k, k = 1,2,...,min(m,n),
is of the form:

**QRSOLVE Subroutine**

"qrsolve" solves a*x = b and d*x = 0 in the least squares sense;
used with routine "qrfact" to solve least squares problems

**QUATFIT Subroutine**

"quatfit" uses a quaternion-based method to achieve the best
fit superposition of two sets of coordinates

**RADIAL Program**

"radial" finds the radial distribution function for a specified
pair of atom types via analysis of a set of coordinate frames

**RANDOM Function**

"random" generates a random number on [0,1] via a long
period generator due to L'Ecuyer with Bays-Durham shuffle

**RANVEC Subroutine**

"ranvec" generates a unit vector in 3-dimensional
space with uniformly distributed random orientation

**RATTLE Subroutine**

"rattle" implements the first portion of the RATTLE algorithm
by correcting atomic positions and half-step velocities to
maintain interatomic distance and absolute spatial constraints

**RATTLE2 Subroutine**

"rattle2" implements the second portion of the RATTLE algorithm
by correcting the full-step velocities in order to maintain
interatomic distance constraints

**READBLK Subroutine**

"readblk" reads in a set of snapshot frames and transfers
the values to internal arrays for use in the computation
of time correlation functions

**READDYN Subroutine**

"readdyn" get the positions, velocities and accelerations
for a molecular dynamics restart from an external disk file

**READGARC Subroutine**

"readgarc" reads data from Gaussian archive section; each
entry is terminated with a backslash symbol

**READGAU Subroutine**

"readgau" reads an ab initio optimized structure, forces,
Hessian and frequencies from a Gaussian 09 output file

**READGDMA Subroutine**

"readgdma" takes the DMA output in spherical harmonics from
the GDMA program and converts to Cartesian multipoles in
the global coordinate frame

**READINT Subroutine**

"readint" gets a set of Z-matrix internal coordinates
from an external file

**READMOL Subroutine**

"readmol" gets a set of MDL MOL coordinates from
an external disk file

**READMOL2 Subroutine**

"readmol2" gets a set of Tripos MOL2 coordinates from an
external disk file

**READPDB Subroutine**

"readpdb" gets a set of Protein Data Bank coordinates
from an external disk file

**READPOT Subroutine**

"readpot" gets a set of grid points and target electrostatic
potential values from an external disk file

**READPRM Subroutine**

"readprm" processes the potential energy parameter file
in order to define the default force field parameters

**READSEQ Subroutine**

"readseq" gets a biopolymer sequence containing one or more
separate chains from an external file; all lines containing
sequence must begin with the starting sequence number, the
actual sequence is read from subsequent nonblank characters

**READXYZ Subroutine**

"readxyz" gets a set of Cartesian coordinates from
an external disk file

**REFINE Subroutine**

"refine" performs minimization of the atomic coordinates
of an initial crude embedded distance geometry structure versus
the bound, chirality, planarity and torsional error functions

**RELEASEMONITOR Subroutine**

**REPLICA Subroutine**

"replica" decides between images and replicates for generation
of periodic boundary conditions, and sets the cell replicate
list if the replicates method is to be used

**RESPA Subroutine**

"respa" performs a single multiple time step molecular dynamics
step using the reversible reference system propagation algorithm
(r-RESPA) via a Verlet core with the potential split into fast-
and slow-evolving portions

**RFINDEX Subroutine**

"rfindex" finds indices for each multipole site for use
in computing reaction field energetics

**RGDSTEP Subroutine**

"rgdstep" performs a single molecular dynamics time step
via a rigid body integration algorithm

**RIBOSOME Subroutine**

"ribosome" translates a polypeptide structure in Protein Data
Bank format to a Cartesian coordinate file and sequence file

**RIGIDXYZ Subroutine**

"rigidxyz" computes Cartesian coordinates for a rigid body
group via rotation and translation of reference coordinates

**RINGS Subroutine**

"rings" searches the structure for small rings and stores
their constituent atoms, and optionally reduces large rings
into their component smaller rings

**RMSERROR Subroutine**

"rmserror" computes the maximum absolute deviation and the
rms deviation from the distance bounds, and the number and
rms value of the distance restraint violations

**RMSFIT Function**

"rmsfit" computes the rms fit of two coordinate sets

**ROTANG Function**

**ROTCHECK Function**

"rotcheck" tests a specified candidate rotatable bond for
the disallowed case where inactive atoms are found on both
sides of the candidate bond

**ROTEULER Subroutine**

"roteuler" computes a set of Euler angle values consistent
with an input rotation matrix

**ROTFRAME Subroutine**

"rotframe" takes the global multipole moments and rotates them
into the local coordinate frame defined at each atomic site

**ROTLIST Subroutine**

"rotlist" generates the minimum list of all the atoms lying
to one side of a pair of directly bonded atoms; optionally
finds the minimal list by choosing the side with fewer atoms

**ROTMAT Subroutine**

"rotmat" finds the rotation matrix that rotates the local
coordinate system into the global frame at a multipole site

**ROTPOLE Subroutine**

"rotpole" constructs the set of atomic multipoles in the global
frame by applying the correct rotation matrix for each site

**ROTRGD Subroutine**

"rotrgd" finds the rotation matrix for a rigid body due
to a single step of dynamics

**ROTSITE Subroutine**

"rotsite" rotates the local frame atomic multipoles at a
specified site into the global coordinate frame by applying
a rotation matrix

**SADDLE Program**

"saddle" finds a transition state between two conformational
minima using a combination of ideas from the synchronous transit
(Halgren-Lipscomb) and quadratic path (Bell-Crighton) methods

**SADDLE1 Function**

"saddle1" is a service routine that computes the energy and
gradient for transition state optimization

**SADDLES Subroutine**

"saddles" constructs circles, convex edges and saddle faces

**SAVEYZE Subroutine**

"saveyze" prints the atomic forces and/or the induced dipoles
to separate external disk files

**SBGUESS Subroutine**

"sbguess" sets approximate stretch-bend force constants based
on atom type and connected atoms

**SCAN Program**

"scan" attempts to find all the local minima on a potential
energy surface via an iterative series of local searches along
normal mode directions

**SCAN1 Function**

"scan1" is a service routine that computes the energy and
gradient during exploration of a potential energy surface
via iterative local search

**SCAN2 Subroutine**

"scan2" is a service routine that computes the sparse matrix
Hessian elements during exploration of a potential energy
surface via iterative local search

**SCANPDB Subroutine**

"scanpdb" reads the first model in a Protein Data Bank file and
sets chains, alternate sites and insertion records to be used

**SDAREA Subroutine**

"sdarea" optionally scales the atomic friction coefficient
of each atom based on its accessible surface area

**SDSTEP Subroutine**

"sdstep" performs a single stochastic dynamics time step
via the velocity Verlet integration algorithm

**SDTERM Subroutine**

"sdterm" finds the frictional and random terms needed to
update positions and velocities during stochastic dynamics

**SEARCH Subroutine**

"search" is a unidimensional line search based upon parabolic
extrapolation and cubic interpolation using both function and
gradient values

**SETACCELERATION Subroutine**

**SETATOMIC Subroutine**

**SETATOMTYPES Subroutine**

**SETCHARGE Subroutine**

**SETCHUNK Subroutine**

"setchunk" marks a chunk in the PME spatial table which is
overlapped by the B-splines for a site

**SETCONNECTIVITY Subroutine**

**SETCOORDINATES Subroutine**

**SETELECT Subroutine**

"setelect" assigns partial charge, bond dipole and atomic
multipole parameters for the current structure, as needed
for computation of the electrostatic potential

**SETENERGY Subroutine**

**SETFILE Subroutine**

**SETFORCEFIELD Subroutine**

**SETFRAME Subroutine**

"setframe" assigns a local coordinate frame at each atomic
multipole site using high priority connected atoms along axes

**SETGRADIENTS Subroutine**

**SETINDUCED Subroutine**

**SETKEYWORD Subroutine**

**SETMASS Subroutine**

**SETMDTIME Subroutine**

**SETMOL2 Program**

"setmol2" assigns MOL2 atom names/types/charges and bond types
based upon atomic numbers and connectivity

**SETNAME Subroutine**

**SETPAIR Program**

"setpair" is a service routine that assigns flags, sets cutoffs
and allocates arrays used by different pairwise neighbor methods

**SETPOLAR Subroutine**

"setpolar" assigns atomic polarizabilities, Thole damping or
charge penetration parameters, and polarization groups with
user modification of these values

**SETSTEP Subroutine**

**SETSTORY Subroutine**

**SETTIME Subroutine**

"settime" initializes the wall clock and elapsed CPU times

**SETUPDATED Subroutine**

**SETVELOCITY Subroutine**

**SHAKE Subroutine**

"shake" implements the SHAKE algorithm by correcting atomic
positions to maintain interatomic distance and absolute spatial
constraints

**SHAKEF Subroutine**

"shakef" modifies the gradient to remove components along any
holonomic distance contraints using a variant of SHAKE

**SHAKEUP Subroutine**

"shakeup" initializes any holonomic constraints for use with
the SHAKE and RATTLE algorithms

**SHROTMAT Subroutine**

"shrotmat" finds the rotation matrix that converts spherical
harmonic quadrupoles from the local to the global frame given
the required dipole rotation matrix

**SHROTSITE Subroutine**

"shrotsite" converts spherical harmonic multipoles from the
local to the global frame given required rotation matrices

**SIGMOID Function**

"sigmoid" implements a normalized sigmoidal function on the
interval [0,1]; the curves connect (0,0) to (1,1) and have
a cooperativity controlled by beta, they approach a straight
line as beta -> 0 and get more nonlinear as beta increases

**SIMPLEX Subroutine**

"simplex" is a general multidimensional Nelder-Mead simplex
optimization routine requiring only repeated evaluations of
the objective function

**SIMPLEX1 Function**

"simplex1" is a service routine used only by the Nelder-Mead
simplex optimization method

**SKTDYN Subroutine**

"sktdyn" sends the current dynamics info via a socket

**SKTINIT Subroutine**

"sktinit" sets up socket communication with the graphical
user interface by starting a Java virtual machine, initiating
a server, and loading an object with system information

**SKTKILL Subroutine**

"sktkill" closes the server and Java virtual machine

**SKTOPT Subroutine**

"sktopt" sends the current optimization info via a socket

**SLATER Subroutine**

"slater" is a general routine for computing the overlap
integrals between two Slater-type orbitals

**SNIFFER Program**

"sniffer" performs a global energy minimization using a
discrete version of Griewank's global search trajectory

**SNIFFER1 Function**

"sniffer1" is a service routine that computes the energy
and gradient for the Sniffer global optimization method

**SOAK Subroutine**

"soak" takes a currently defined solute system and places
it into a solvent box, with removal of any solvent molecules
that overlap the solute

**SORT Subroutine**

"sort" takes an input list of integers and sorts it
into ascending order using the Heapsort algorithm

**SORT10 Subroutine**

"sort10" takes an input list of character strings and sorts
it into alphabetical order using the Heapsort algorithm,
duplicate values are removed from the final sorted list

**SORT2 Subroutine**

"sort2" takes an input list of reals and sorts it
into ascending order using the Heapsort algorithm;
it also returns a key into the original ordering

**SORT3 Subroutine**

"sort3" takes an input list of integers and sorts it
into ascending order using the Heapsort algorithm;
it also returns a key into the original ordering

**SORT4 Subroutine**

"sort4" takes an input list of integers and sorts it into
ascending absolute value using the Heapsort algorithm

**SORT5 Subroutine**

"sort5" takes an input list of integers and sorts it
into ascending order based on each value modulo "m"

**SORT6 Subroutine**

"sort6" takes an input list of character strings and sorts
it into alphabetical order using the Heapsort algorithm

**SORT7 Subroutine**

"sort7" takes an input list of character strings and sorts it
into alphabetical order using the Heapsort algorithm; it also
returns a key into the original ordering

**SORT8 Subroutine**

"sort8" takes an input list of integers and sorts it into
ascending order using the Heapsort algorithm, duplicate
values are removed from the final sorted list

**SORT9 Subroutine**

"sort9" takes an input list of reals and sorts it into
ascending order using the Heapsort algorithm, duplicate
values are removed from the final sorted list

**SPACEFILL Program**

"spacefill" computes the surface area and volume of
a structure; the van der Waals, accessible-excluded,
and contact-reentrant definitions are available

**SPECTRUM Program**

"spectrum" computes a power spectrum over a wavelength range
from the velocity autocorrelation as a function of time

**SPHERE Subroutine**

"sphere" finds a specified number of uniformly distributed
points on a sphere of unit radius centered at the origin

**SQUARE Subroutine**

"square" is a nonlinear least squares routine derived from the
IMSL BCLSF routine and the MINPACK LMDER routine; the Jacobian
is estimated by finite differences and bounds can be specified
for the variables to be refined

**SUFFIX Subroutine**

"suffix" checks a filename for the presence of an extension,
and appends an extension and version if none is found

**SUPERPOSE Program**

"superpose" takes pairs of structures and superimposes them
in the optimal least squares sense; it will attempt to match
all atom pairs or only those specified by the user

**SURFACE Subroutine**

"surface" performs an analytical computation of the weighted
solvent accessible surface area of each atom and the first
derivatives of the area with respect to Cartesian coordinates

**SURFACE1 Subroutine**

"surface1" performs an analytical computation of the weighted
solvent accessible surface area of each atom and the first
derivatives of the area with respect to Cartesian coordinates

**SURFATOM Subroutine**

"surfatom" performs an analytical computation of the surface
area of a specified atom; a simplified version of "surface"

**SURFATOM1 Subroutine**

"surfatom1" performs an analytical computation of the surface
area and first derivatives with respect to Cartesian coordinates
of a specified atom

**SWITCH Subroutine**

"switch" sets the coeffcients used by the fifth and seventh
order polynomial switching functions for spherical cutoffs

**SYMMETRY Subroutine**

"symmetry" applies symmetry operators to the fractional
coordinates of the asymmetric unit in order to generate
the symmetry related atoms of the full unit cell

**SYSTYZE Subroutine**

"systyze" is an auxiliary routine for the analyze program
that prints general information about the molecular system
and the force field model

**TABLE_FILL Subroutine**

"table_fill" constructs an array which stores the spatial
regions of the particle mesh Ewald grid with contributions
from each site

**TANGENT Subroutine**

"tangent" finds the projected gradient on the synchronous
transit path for a point along the transit pathway

**TCGSWAP Subroutine**

"tcgswap" switches two sets of induced dipole quantities for
use with the TCG induced dipole solver

**TCG_ALPHA12 Subroutine**

"tcg_alpha12" computes source1 = alpha*source1 and
source2 = alpha*source2

**TCG_ALPHA22 Subroutine**

"tcg_alpha22" computes result1 = alpha*source1 and
result2 = alpha*source2

**TCG_ALPHAQUAD Subroutine**

"tcg_alphaquad" computes the quadratic form, <a*alpha*b>,
where alpha is the diagonal atomic polarizability matrix

**TCG_DOTPROD Subroutine**

"tcg_dotprod" computes the dot product of two vectors
of length n elements

**TCG_RESOURCE Subroutine**

"tcg_resource" sets the number of mutual induced dipole
pairs based on the passed argument

**TCG_T0 Subroutine**

"tcg_t0" applies T matrix to ind/p, and returns v3d/p
T = 1/alpha + Tu

**TCG_UFIELD Subroutine**

"tcg_ufield" applies -Tu to ind/p and returns v3d/p

**TCG_UPDATE Subroutine**

"tcg_update" computes pvec = alpha*rvec + beta*pvec;
if the preconditioner is not used, then alpha = identity

**TEMPER Subroutine**

"temper" computes the instantaneous temperature and applies a
thermostat via Berendsen or Bussi-Parrinello velocity scaling,
Andersen stochastic collisions or Nose-Hoover chains; also uses
Berendsen scaling for any iEL induced dipole variables

**TEMPER2 Subroutine**

"temper2" applies a velocity correction at the half time step
as needed for the Nose-Hoover thermostat

**TESTGRAD Program**

"testgrad" computes and compares the analytical and numerical
gradient vectors of the potential energy function with respect
to Cartesian coordinates

**TESTHESS Program**

"testhess" computes and compares the analytical and numerical
Hessian matrices of the potential energy function with respect
to Cartesian coordinates

**TESTPAIR Program**

"testpair" performs a set of timing tests to compare the
evaluation of potential energy and energy/gradient using
different methods for finding pairwise neighbors

**TESTPOL Program**

"testpol" compares the induced dipoles from direct polarization,
mutual SCF iterations, perturbation theory extrapolation (OPT),
and truncated conjugate gradient (TCG) solvers

**TESTROT Program**

"testrot" computes and compares the analytical and numerical
gradient vectors of the potential energy function with respect
to rotatable torsional angles

**TESTVIR Program**

"testvir" computes the analytical internal virial and compares
it to a numerical virial derived from the finite difference
derivative of the energy with respect to lattice vectors

**TIMER Program**

"timer" measures the CPU time required for file reading and
parameter assignment, potential energy computation, energy
and gradient computation, and Hessian matrix evaluation

**TIMEROT Program**

"timerot" measures the CPU time required for file reading
and parameter assignment, potential energy computation,
energy and gradient over torsions, and torsional angle
Hessian matrix evaluation

**TNCG Subroutine**

"tncg" implements a truncated Newton optimization algorithm
in which a preconditioned linear conjugate gradient method is
used to approximately solve Newton's equations; special features
include use of an explicit sparse Hessian or finite-difference
gradient-Hessian products within the PCG iteration; the exact
Newton search directions can be used optionally; by default the
algorithm checks for negative curvature to prevent convergence
to a stationary point having negative eigenvalues; if a saddle
point is desired this test can be removed by disabling "negtest"

**TNSOLVE Subroutine**

"tnsolve" uses a linear conjugate gradient method to find
an approximate solution to the set of linear equations
represented in matrix form by Hp = -g (Newton's equations)

**TORFIT1 Function**

"torfit1" is a service routine that computes the energy and
gradient for a low storage BFGS optimization in Cartesian
coordinate space

**TORGUESS Subroutine**

"torguess" set approximate torsion amplitude parameters based
on atom type and connected atoms

**TORPHASE Subroutine**

"torphase" sets the n-fold amplitude and phase values
for each torsion via sorting of the input parameters

**TORQUE Subroutine**

"torque" takes the torque values on a single site defined by
a local coordinate frame and converts to Cartesian forces on
the original site and sites specifying the local frame, also
gives the x,y,z-force components needed for virial computation

**TORSER Function**

"torser" computes the torsional error function and its first
derivatives with respect to the atomic Cartesian coordinates
based on the deviation of specified torsional angles from
desired values, the contained bond angles are also restrained
to avoid a numerical instability

**TORSFIT Program**

"torsfit" refines torsional force field parameters based on
a quantum mechanical potential surface and analytical gradient

**TORSIONS Subroutine**

"torsions" finds the total number of torsional angles and
the numbers of the four atoms defining each torsional angle

**TORUS Subroutine**

"torus" sets a list of all of the temporary torus positions
by testing for a torus between each atom and its neighbors

**TOTERR Function**

"toterr" is the error function and derivatives for a distance
geometry embedding; it includes components from the distance
bounds, hard sphere contacts, local geometry, chirality and
torsional restraint errors

**TRANSFORM Subroutine**

"transform" diagonalizes the current basis vectors to produce
trial roots for sliding block iterative matrix diagonalization

**TRANSIT Function**

"transit" evaluates the synchronous transit function and
gradient; linear and quadratic transit paths are available

**TRBASIS Subroutine**

"trbasis" forms translation and rotation basis vectors used
during vibrational analysis via block iterative diagonalization

**TRIANGLE Subroutine**

"triangle" smooths the upper and lower distance bounds via
the triangle inequality using a full-matrix variant of the
Floyd-Warshall shortest path algorithm; this routine is
usually much slower than the sparse matrix shortest path
methods in "geodesic" and "trifix", and should be used only
for comparison with answers generated by those routines

**TRIFIX Subroutine**

"trifix" rebuilds both the upper and lower distance bound
matrices following tightening of one or both of the bounds
between a specified pair of atoms, "p" and "q", using a
modification of Murchland's shortest path update algorithm

**TRIGGER Subroutine**

"trigger" constructs a set of initial trial vectors for
use during sliding block iterative matrix diagonalization

**TRIMHEAD Subroutine**

"trimhead" removes blank spaces before the first non-blank
character in a text string by shifting the string to the left

**TRIMTEXT Function**

"trimtext" finds and returns the location of the last
non-blank character before the first null character in
an input text string; the function returns zero if no
such character is found

**TRIPLE Function**

"triple" finds the triple product of three vectors; used as
a service routine by the Connolly surface area and volume
computation

**TRUST Subroutine**

"trust" updates the model trust region for a nonlinear least
squares calculation based on ideas found in NL2SOL and Dennis
and Schnabel's book

**UBUILD Subroutine**

"ubuild" performs a complete rebuild of the polarization
preconditioner neighbor list for all sites

**UDIRECT1 Subroutine**

"udirect1" computes the reciprocal space contribution of the
permanent atomic multipole moments to the field

**UDIRECT2A Subroutine**

"udirect2a" computes the real space contribution of the permanent
atomic multipole moments to the field via a double loop

**UDIRECT2B Subroutine**

"udirect2b" computes the real space contribution of the permanent
atomic multipole moments to the field via a neighbor list

**UFIELD0A Subroutine**

"ufield0a" computes the mutual electrostatic field due to
induced dipole moments via a double loop

**UFIELD0B Subroutine**

"ufield0b" computes the mutual electrostatic field due to
induced dipole moments via a pair list

**UFIELD0C Subroutine**

"ufield0c" computes the mutual electrostatic field due to
induced dipole moments via Ewald summation

**UFIELD0D Subroutine**

"ufield0d" computes the mutual electrostatic field due to
induced dipole moments for use with with generalized Kirkwood
implicit solvation

**UFIELD0E Subroutine**

"ufield0e" computes the mutual electrostatic field due to
induced dipole moments via a Poisson-Boltzmann solver

**UFIELDI Subroutine**

"ufieldi" computes the electrostatic field due to intergroup
induced dipole moments

**ULIGHT Subroutine**

"ulight" performs a complete rebuild of the polarization
preconditioner pair neighbor list for all sites using the
method of lights

**ULIST Subroutine**

"ulist" performs an update or a complete rebuild of the
neighbor lists for the polarization preconditioner

**ULSPRED Subroutine**

"ulspred" uses standard extrapolation or a least squares fit
to set coefficients of an induced dipole predictor polynomial

**UMUTUAL1 Subroutine**

"umutual1" computes the reciprocal space contribution of the
induced atomic dipole moments to the field

**UMUTUAL2A Subroutine**

"umutual2a" computes the real space contribution of the induced
atomic dipole moments to the field via a double loop

**UMUTUAL2B Subroutine**

"umutual2b" computes the real space contribution of the induced
atomic dipole moments to the field via a neighbor list

**UNITCELL Subroutine**

"unitcell" gets the periodic boundary box size and related
values from an external keyword file

**UPCASE Subroutine**

"upcase" converts a text string to all upper case letters

**URYGUESS Function**

"uryguess" sets approximate Urey-Bradley force constants
based on atom type and connected atoms

**USCALE0A Subroutine**

"uscale0a" builds and applies a preconditioner for the conjugate
gradient induced dipole solver using a double loop

**USCALE0B Subroutine**

"uscale0b" builds and applies a preconditioner for the conjugate
gradient induced dipole solver using a neighbor pair list

**VALENCE Program**

"valence" refines force field parameters for valence terms based
on a quantum mechanical optimized structure and frequencies

**VALFIT1 Function**

"valfit1" is a service routine that computes the RMS error
and gradient for valence parameters fit to QM results

**VALGUESS Subroutine**

"valguess" sets approximate valence parameter values based on
quantum mechanical structure and frequency data

**VALMIN1 Function**

"valmin1" is a service routine that computes the molecular
energy and gradient during valence parameter optimization

**VALRMS Function**

"valrms" evaluates a valence parameter goodness-of-fit error
function based on comparison of forces, frequencies, bond
lengths and angles to QM results

**VAM Subroutine**

"vam" takes the analytical molecular surface defined
as a collection of spherical and toroidal polygons
and uses it to compute the volume and surface area

**VARPRM Subroutine**

"varprm" copies the current optimization values into the
corresponding electrostatic potential energy parameters

**VARPRM Subroutine**

"varprm" copies the current optimization values into the
corresponding valence potential energy parameters

**VBUILD Subroutine**

"vbuild" performs a complete rebuild of the van der Waals
pair neighbor list for all sites

**VCROSS Subroutine**

"vcross" finds the cross product of two vectors

**VDWERR Function**

"vdwerr" is the hard sphere van der Waals bound error function
and derivatives that penalizes close nonbonded contacts,
pairwise neighbors are generated via the method of lights

**VDWGUESS Subroutine**

"vdwguess" sets initial VDW parameters based on atom type
and connected atoms

**VECANG Function**

"vecang" finds the angle between two vectors handed with respect
to a coordinate axis; returns an angle in the range [0,2*pi]

**VERLET Subroutine**

"verlet" performs a single molecular dynamics time step
via the velocity Verlet multistep recursion formula

**VERSION Subroutine**

"version" checks the name of a file about to be opened; if
if "old" status is passed, the name of the highest current
version is returned; if "new" status is passed the filename
of the next available unused version is generated

**VIBBIG Program**

"vibbig" performs large-scale vibrational mode analysis using
only vector storage and gradient evaluations; preconditioning
is via an approximate inverse from a block diagonal Hessian,
and a sliding block method is used to converge any number of
eigenvectors starting from either lowest or highest frequency

**VIBRATE Program**

"vibrate" performs a vibrational normal mode analysis; the
Hessian matrix of second derivatives is determined and then
diagonalized both directly and after mass weighting; output
consists of the eigenvalues of the force constant matrix as
well as the vibrational frequencies and displacements

**VIBROT Program**

"vibrot" computes the eigenvalues and eigenvectors of the
torsional Hessian matrix

**VIRIYZE Subroutine**

"propyze" finds and prints the internal virial, the dE/dV value
and an estimate of the pressure

**VLIGHT Subroutine**

"vlight" performs a complete rebuild of the van der Waals
pair neighbor list for all sites using the method of lights

**VLIST Subroutine**

"vlist" performs an update or a complete rebuild of the
nonbonded neighbor lists for vdw sites

**VNORM Subroutine**

"vnorm" normalizes a vector to unit length; used as a
service routine by the Connolly surface area and volume
computation

**VOLUME Subroutine**

"volume" calculates the excluded volume via the Connolly
analytical volume and surface area algorithm

**VOLUME1 Subroutine**

"volume1" calculates first derivatives of the total excluded
volume with respect to the Cartesian coordinates of each atom

**VOLUME2 Subroutine**

"volume2" calculates second derivatives of the total excluded
volume with respect to the Cartesian coordinates of the atoms

**WATSON Subroutine**

"watson" uses a rigid body optimization to approximately
align the paired strands of a nucleic acid double helix

**WATSON1 Function**

"watson1" is a service routine that computes the energy
and gradient for optimally conditioned variable metric
optimization of rigid bodies

**WIGGLE Subroutine**

"wiggle" applies a random perturbation to the atomic coordinates
to avoid numerical instabilities for various linear, planar and
symmetric structures

**XTALERR Subroutine**

"xtalerr" computes an error function value derived from
lattice energies, dimer intermolecular energies and the
gradient with respect to structural parameters

**XTALFIT Program**

"xtalfit" determines optimized van der Waals and electrostatic
parameters by fitting to crystal structures, lattice energies,
and dimer structures and interaction energies

**XTALMIN Program**

"xtalmin" performs a full crystal energy minimization by
optimizing over fractional atomic coordinates and the six
lattice lengths and angles

**XTALMIN1 Function**

"xtalmin1" is a service routine that computes the energy and
gradient with respect to fractional coordinates and lattice
dimensions for a crystal energy minimization

**XTALMOVE Subroutine**

"xtalmove" converts fractional to Cartesian coordinates for
rigid molecules during optimization of force field parameters

**XTALPRM Subroutine**

"xtalprm" stores or retrieves a molecular structure; used to
make a previously stored structure the active structure, or to
store a structure for later use

**XTALWRT Subroutine**

"xtalwrt" prints intermediate results during fitting of
force field parameters to structures and energies

**XYZATM Subroutine**

"xyzatm" computes the Cartesian coordinates of a single
atom from its defining internal coordinate values

**XYZEDIT Program**

"xyzedit" provides for modification and manipulation
of the contents of Cartesian coordinates files

**XYZINT Program**

"xyzint" takes as input a Cartesian coordinates file, then
converts to and writes out an internal coordinates file

**XYZMOL2 Program**

"xyzmol2" takes as input a Cartesian coordinates file,
converts to and then writes out a Tripos MOL2 file

**XYZPDB Program**

"xyzpdb" takes as input a Cartesian coordinates file,
then converts to and writes out a Protein Data Bank file

**XYZRIGID Subroutine**

"xyzrigid" computes the center of mass and Euler angle rigid
body coordinates for each atom group in the system

**ZATOM Subroutine**

"zatom" adds an atom to the end of the current Z-matrix
and then increments the atom counter; atom type, defining
atoms and internal coordinates are passed as arguments

**ZHELP Subroutine**

"zhelp" prints the general information and instructions
for the Z-matrix editing program

**ZVALUE Subroutine**

"zvalue" gets user supplied values for selected coordinates
as needed by the internal coordinate editing program
