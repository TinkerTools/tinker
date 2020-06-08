Modules & Global Variables
==========================

The Fortran modules found in the Tinker package are listed below along with a brief description of the variables associated with each module. Each individual module contains a set of globally allocated variables available to any program unit upon inclusion of that module. A source listing containing each of the Tinker functions and subroutines and its included modules can be produced by running the "listing.make" script found in the distribution.

**MODULE ACTION        total number of each energy term type

.. code-block:: text

 neb             number of bond stretch energy terms computed
 nea             number of angle bend energy terms computed
 neba            number of stretch-bend energy terms computed
 neub            number of Urey-Bradley energy terms computed
 neaa            number of angle-angle energy terms computed
 neopb           number of out-of-plane bend energy terms computed
 neopd           number of out-of-plane distance energy terms computed
 neid            number of improper dihedral energy terms computed
 neit            number of improper torsion energy terms computed
 net             number of torsional energy terms computed
 nept            number of pi-system torsion energy terms computed
 nebt            number of stretch-torsion energy terms computed
 neat            number of angle-torsion energy terms computed
 nett            number of torsion-torsion energy terms computed
 nev             number of van der Waals energy terms computed
 ner             number of Pauli repulsion energy terms computed
 nedsp           number of dispersion energy terms computed
 nec             number of charge-charge energy terms computed
 necd            number of charge-dipole energy terms computed
 ned             number of dipole-dipole energy terms computed
 nem             number of multipole energy terms computed
 nep             number of polarization energy terms computed
 nect            number of charge transfer energy terms computed
 new             number of Ewald summation energy terms computed
 nerxf           number of reaction field energy terms computed
 nes             number of solvation energy terms computed
 nelf            number of metal ligand field energy terms computed
 neg             number of geometric restraint energy terms computed
 nex             number of extra energy terms computed

**MODULE ALIGN         information for structure superposition

.. code-block:: text

 nfit            number of atoms to use in superimposing two structures
 ifit            atom numbers of pairs of atoms to be superimposed
 wfit            weights assigned to atom pairs during superposition

**MODULE ANALYZ        energy components partitioned to atoms

.. code-block:: text

 aesum           total potential energy partitioned over atoms
 aeb             bond stretch energy partitioned over atoms
 aea             angle bend energy partitioned over atoms
 aeba            stretch-bend energy partitioned over atoms
 aeub            Urey-Bradley energy partitioned over atoms
 aeaa            angle-angle energy partitioned over atoms
 aeopb           out-of-plane bend energy partitioned over atoms
 aeopd           out-of-plane distance energy partitioned over atoms
 aeid            improper dihedral energy partitioned over atoms
 aeit            improper torsion energy partitioned over atoms
 aet             torsional energy partitioned over atoms
 aept            pi-system torsion energy partitioned over atoms
 aebt            stretch-torsion energy partitioned over atoms
 aeat            angle-torsion energy partitioned over atoms
 aett            torsion-torsion energy partitioned over atoms
 aev             van der Waals energy partitioned over atoms
 aer             Pauli repulsion energy partitioned over atoms
 aedsp           damped dispersion energy partitioned over atoms
 aec             charge-charge energy partitioned over atoms
 aecd            charge-dipole energy partitioned over atoms
 aed             dipole-dipole energy partitioned over atoms
 aem             multipole energy partitioned over atoms
 aep             polarization energy partitioned over atoms
 aect            charge transfer energy partitioned over atoms
 aerxf           reaction field energy partitioned over atoms
 aes             solvation energy partitioned over atoms
 aelf            metal ligand field energy partitioned over atoms
 aeg             geometric restraint energy partitioned over atoms
 aex             extra energy term partitioned over atoms

**MODULE ANGANG        angle-angles in current structure

.. code-block:: text

 nangang         total number of angle-angle interactions
 iaa             angle numbers used in each angle-angle term
 kaa             force constant for angle-angle cross terms

**MODULE ANGBND        bond angle bends in current structure

.. code-block:: text

 nangle          total number of angle bends in the system
 iang            numbers of the atoms in each angle bend
 ak              harmonic angle force constant (kcal/mole/rad**2)
 anat            ideal bond angle or phase shift angle (degrees)
 afld            periodicity for Fourier angle bending term

**MODULE ANGPOT        angle bend functional form details

.. code-block:: text

 angunit         convert angle bending energy to kcal/mole
 stbnunit        convert stretch-bend energy to kcal/mole
 aaunit          convert angle-angle energy to kcal/mole
 opbunit         convert out-of-plane bend energy to kcal/mole
 opdunit         convert out-of-plane distance energy to kcal/mole
 cang            cubic coefficient in angle bending potential
 qang            quartic coefficient in angle bending potential
 pang            quintic coefficient in angle bending potential
 sang            sextic coefficient in angle bending potential
 copb            cubic coefficient in out-of-plane bend potential
 qopb            quartic coefficient in out-of-plane bend potential
 popb            quintic coefficient in out-of-plane bend potential
 sopb            sextic coefficient in out-of-plane bend potential
 copd            cubic coefficient in out-of-plane distance potential
 qopd            quartic coefficient in out-of-plane distance potential
 popd            quintic coefficient in out-of-plane distance potential
 sopd            sextic coefficient in out-of-plane distance potential
 opbtyp          type of out-of-plane bend potential energy function
 angtyp          type of angle bending function for each bond angle

**MODULE ANGTOR        angle-torsions in current structure

.. code-block:: text

 nangtor         total number of angle-torsion interactions
 iat             torsion and angle numbers used in angle-torsion
 kant            1-, 2- and 3-fold angle-torsion force constants

**MODULE ARGUE         command line arguments at run time

.. code-block:: text

 maxarg          maximum number of command line arguments
 narg            number of command line arguments to the program
 listarg         flag to mark available command line arguments
 arg             strings containing the command line arguments

**MODULE ASCII         selected ASCII character code values

.. code-block:: text

 null            decimal value of ASCII code for null (0)
 tab             decimal value of ASCII code for tab (9)
 linefeed        decimal value of ASCII code for linefeed (10)
 formfeed        decimal value of ASCII code for formfeed (12)
 carriage        decimal value of ASCII code for carriage return (13)
 escape          decimal value of ASCII code for escape (27)
 space           decimal value of ASCII code for blank space (32)
 exclamation     decimal value of ASCII code for exclamation (33)
 quote           decimal value of ASCII code for double quote (34)
 pound           decimal value of ASCII code for pound sign (35)
 dollar          decimal value of ASCII code for dollar sign (36)
 percent         decimal value of ASCII code for percent sign (37)
 ampersand       decimal value of ASCII code for ampersand (38)
 apostrophe      decimal value of ASCII code for single quote (39)
 asterisk        decimal value of ASCII code for asterisk (42)
 plus            decimal value of ASCII code for plus sign (43)
 comma           decimal value of ASCII code for comma (44)
 minus           decimal value of ASCII code for minus sign (45)
 period          decimal value of ASCII code for period (46)
 frontslash      decimal value of ASCII codd for frontslash (47)
 colon           decimal value of ASCII code for colon (58)
 semicolon       decimal value of ASCII code for semicolon (59)
 equal           decimal value of ASCII code for equal sign (61)
 question        decimal value of ASCII code for question mark (63)
 atsign          decimal value of ASCII code for at sign (64)
 backslash       decimal value of ASCII code for backslash (92)
 caret           decimal value of ASCII code for caret (94)
 underbar        decimal value of ASCII code for underbar (95)
 vertical        decimal value of ASCII code for vertical bar (124)
 tilde           decimal value of ASCII code for tilde (126)

**MODULE ATMLST        bond and angle local geometry indices

.. code-block:: text

 bndlist         numbers of the bonds involving each atom
 anglist         numbers of the angles centered on each atom
 balist          numbers of the bonds comprising each angle

**MODULE ATOMID        atomic properties for current atoms

.. code-block:: text

 tag             integer atom labels from input coordinates file
 class           atom class number for each atom in the system
 atomic          atomic number for each atom in the system
 valence         valence number for each atom in the system
 mass            atomic weight for each atom in the system
 name            atom name for each atom in the system
 story           descriptive type for each atom in system

**MODULE ATOMS         number, position and type of atoms

.. code-block:: text

 n               total number of atoms in the current system
 type            atom type number for each atom in the system
 x               current x-coordinate for each atom in the system
 y               current y-coordinate for each atom in the system
 z               current z-coordinate for each atom in the system

**MODULE BATH          thermostat and barostat control values

.. code-block:: text

 maxnose         maximum length of Nose-Hoover thermostat chain
 voltrial        mean number of steps between Monte Carlo moves
 kelvin          target value for the system temperature (K)
 atmsph          target value for the system pressure (atm)
 tautemp         time constant for Berendsen thermostat (psec)
 taupres         time constant for Berendsen barostat (psec)
 compress        isothermal compressibility of medium (atm-1)
 collide         collision frequency for Andersen thermostat
 eta             velocity value for Bussi-Parrinello barostat
 volmove         maximum volume move for Monte Carlo barostat (Ang**3)
 vbar            velocity of log volume for Nose-Hoover barostat
 qbar            mass of the volume for Nose-Hoover barostat
 gbar            force for the volume for Nose-Hoover barostat
 vnh             velocity of each chained Nose-Hoover thermostat
 qnh             mass for each chained Nose-Hoover thermostat
 gnh             force for each chained Nose-Hoover thermostat
 isothermal      logical flag governing use of temperature control
 isobaric        logical flag governing use of pressure control
 anisotrop       logical flag governing use of anisotropic pressure
 thermostat      choice of temperature control method to be used
 barostat        choice of pressure control method to be used
 volscale        choice of scaling method for Monte Carlo barostat

**MODULE BITOR         bitorsions in the current structure

.. code-block:: text

 nbitor          total number of bitorsions in the system
 ibitor          numbers of the atoms in each bitorsion

**MODULE BNDPOT        bond stretch functional form details

.. code-block:: text

 cbnd            cubic coefficient in bond stretch potential
 qbnd            quartic coefficient in bond stretch potential
 bndunit         convert bond stretch energy to kcal/mole
 bndtyp          type of bond stretch potential energy function

**MODULE BNDSTR        bond stretches in the current structure

.. code-block:: text

 nbond           total number of bond stretches in the system
 ibnd            numbers of the atoms in each bond stretch
 bk              bond stretch force constants (kcal/mole/Ang**2)
 bl              ideal bond length values in Angstroms

**MODULE BOUND         periodic boundary condition controls

.. code-block:: text

 polycut         cutoff distance for infinite polymer nonbonds
 polycut2        square of infinite polymer nonbond cutoff
 use_bounds      flag to use periodic boundary conditions
 use_replica     flag to use replicates for periodic system
 use_polymer     flag to mark presence of infinite polymer

**MODULE BOXES         periodic boundary condition parameters

.. code-block:: text

 xbox            length of a-axis of periodic box in Angstroms
 ybox            length of b-axis of periodic box in Angstroms
 zbox            length of c-axis of periodic box in Angstroms
 alpha           angle between b- and c-axes of box in degrees
 beta            angle between a- and c-axes of box in degrees
 gamma           angle between a- and b-axes of box in degrees
 xbox2           half of the a-axis length of periodic box
 ybox2           half of the b-axis length of periodic box
 zbox2           half of the c-axis length of periodic box
 box34           three-fourths axis length of truncated octahedron
 volbox          volume in Ang**3 of the periodic box
 beta_sin        sine of the beta periodic box angle
 beta_cos        cosine of the beta periodic box angle
 gamma_sin       sine of the gamma periodic box angle
 gamma_cos       cosine of the gamma periodic box angle
 beta_term       term used in generating triclinic box
 gamma_term      term used in generating triclinic box
 lvec            real space lattice vectors as matrix rows
 recip           reciprocal lattice vectors as matrix columns
 orthogonal      flag to mark periodic box as orthogonal
 monoclinic      flag to mark periodic box as monoclinic
 triclinic       flag to mark periodic box as triclinic
 octahedron      flag to mark box as truncated octahedron
 spacegrp        space group symbol for the unit cell type

**MODULE CELL          replicated cell periodic boundaries

.. code-block:: text

 ncell           total number of cell replicates for periodic boundaries
 icell           offset along axes for each replicate periodic cell
 xcell           length of the a-axis of the complete replicated cell
 ycell           length of the b-axis of the complete replicated cell
 zcell           length of the c-axis of the complete replicated cell
 xcell2          half the length of the a-axis of the replicated cell
 ycell2          half the length of the b-axis of the replicated cell
 zcell2          half the length of the c-axis of the replicated cell

**MODULE CFLUX         charge flux terms in current system

.. code-block:: text

 bflx            bond stretching charge flux constant (electrons/Ang)
 aflx            angle bending charge flux constant (electrons/radian)
 abflx           asymmetric stretch charge flux constant (electrons/Ang)

**MODULE CHARGE        partial charges in current structure

.. code-block:: text

 nion            total number of partial charges in system
 iion            number of the atom site for each partial charge
 jion            neighbor generation site for each partial charge
 kion            cutoff switching site for each partial charge
 pchg            current atomic partial charge values (e-)
 pchg0           original partial charge values for charge flux

**MODULE CHGPEN        charge penetration in current structure

.. code-block:: text

 ncp             total number of charge penetration sites in system
 pcore           number of core electrons at each multipole site
 pval            number of valence electrons at each multipole site
 pval0           original number of valence electrons for charge flux
 palpha          charge penetration damping at each multipole site

**MODULE CHGPOT        charge-charge functional form details

.. code-block:: text

 electric        energy factor in kcal/mole for current force field
 dielec          dielectric constant for electrostatic interactions
 ebuffer         electrostatic buffering constant added to distance
 c1scale         factor by which 1-1 charge interactions are scaled
 c2scale         factor by which 1-2 charge interactions are scaled
 c3scale         factor by which 1-3 charge interactions are scaled
 c4scale         factor by which 1-4 charge interactions are scaled
 c5scale         factor by which 1-5 charge interactions are scaled
 neutnbr         logical flag governing use of neutral group neighbors
 neutcut         logical flag governing use of neutral group cutoffs

**MODULE CHGTRN        charge transfer for current structure

.. code-block:: text

 nct             total number of dispersion sites in the system
 chgct           charge for charge transfer at each multipole site
 dmpct           charge transfer damping factor at each multipole site

**MODULE CHRONO        clock time values for current program

.. code-block:: text

 twall           current processor wall clock time in seconds
 tcpu            elapsed cpu time from start of program in seconds

**MODULE CHUNKS        PME grid spatial decomposition values

.. code-block:: text

 nchunk          total number of spatial regions for PME grid
 nchk1           number of spatial regions along the a-axis
 nchk2           number of spatial regions along the b-axis
 nchk3           number of spatial regions along the c-axis
 ngrd1           number of grid points per region along a-axis
 ngrd2           number of grid points per region along b-axis
 ngrd3           number of grid points per region along c-axis
 nlpts           PME grid points to the left of center point
 nrpts           PME grid points to the right of center point
 grdoff          offset for index into B-spline coefficients
 pmetable        PME grid spatial regions involved for each site

**MODULE COUPLE        atom neighbor connectivity lists

.. code-block:: text

 n12             number of atoms directly bonded to each atom
 n13             number of atoms in a 1-3 relation to each atom
 n14             number of atoms in a 1-4 relation to each atom
 n15             number of atoms in a 1-5 relation to each atom
 i12             atom numbers of atoms 1-2 connected to each atom
 i13             atom numbers of atoms 1-3 connected to each atom
 i14             atom numbers of atoms 1-4 connected to each atom
 i15             atom numbers of atoms 1-5 connected to each atom

**MODULE CTRPOT        charge transfer functional form details

.. code-block:: text

 ctrntyp         type of charge transfer term (SEPARATE or COMBINED)

**MODULE DERIV         Cartesian coord derivative components

.. code-block:: text

 desum           total energy Cartesian coordinate derivatives
 deb             bond stretch Cartesian coordinate derivatives
 dea             angle bend Cartesian coordinate derivatives
 deba            stretch-bend Cartesian coordinate derivatives
 deub            Urey-Bradley Cartesian coordinate derivatives
 deaa            angle-angle Cartesian coordinate derivatives
 deopb           out-of-plane bend Cartesian coordinate derivatives
 deopd           out-of-plane distance Cartesian coordinate derivatives
 deid            improper dihedral Cartesian coordinate derivatives
 deit            improper torsion Cartesian coordinate derivatives
 det             torsional Cartesian coordinate derivatives
 dept            pi-system torsion Cartesian coordinate derivatives
 debt            stretch-torsion Cartesian coordinate derivatives
 deat            angle-torsion Cartesian coordinate derivatives
 dett            torsion-torsion Cartesian coordinate derivatives
 dev             van der Waals Cartesian coordinate derivatives
 der             Pauli repulsion Cartesian coordinate derivatives
 dedsp           damped dispersion Cartesian coordinate derivatives
 dec             charge-charge Cartesian coordinate derivatives
 decd            charge-dipole Cartesian coordinate derivatives
 ded             dipole-dipole Cartesian coordinate derivatives
 dem             multipole Cartesian coordinate derivatives
 dep             polarization Cartesian coordinate derivatives
 dect            charge transfer Cartesian coordinate derivatives
 derxf           reaction field Cartesian coordinate derivatives
 des             solvation Cartesian coordinate derivatives
 delf            metal ligand field Cartesian coordinate derivatives
 deg             geometric restraint Cartesian coordinate derivatives
 dex             extra energy term Cartesian coordinate derivatives

**MODULE DIPOLE        bond dipoles in current structure

.. code-block:: text

 ndipole         total number of dipoles in the system
 idpl            numbers of atoms that define each dipole
 bdpl            magnitude of each of the dipoles (Debye)
 sdpl            position of each dipole between defining atoms

**MODULE DISGEO        distance geometry bounds & parameters

.. code-block:: text

 vdwmax          maximum value of hard sphere sum for an atom pair
 compact         index of local distance compaction on embedding
 pathmax         maximum value of upper bound after smoothing
 dbnd            distance geometry upper and lower bounds matrix
 georad          hard sphere radii for distance geometry atoms
 use_invert      flag to use enantiomer closest to input structure
 use_anneal      flag to use simulated annealing refinement

**MODULE DISP          damped dispersion for current structure

.. code-block:: text

 ndisp           total number of dispersion sites in the system
 idisp           number of the atom for each dispersion site
 csixpr          pairwise sum of C6 dispersion coefficients
 csix            C6 dispersion coefficient value at each site
 adisp           alpha dispersion damping value at each site

**MODULE DMA           QM spherical harmonic multipole moments

.. code-block:: text

 mp              atomic monopole charge values from DMA
 dpx             atomic dipole moment x-component from DMA
 dpy             atomic dipole moment y-component from DMA
 dpz             atomic dipole moment z-component from DMA
 q20             atomic Q20 quadrupole component from DMA (zz)
 q21c            atomic Q21c quadrupole component from DMA (xz)
 q21s            atomic Q21s quadrupole component from DMA (yz)
 q22c            atomic Q22c quadrupole component from DMA (xx-yy)
 q22s            atomic Q22s quadrupole component from DMA (xy)

**MODULE DOMEGA        derivative components over torsions

.. code-block:: text

 tesum           total energy derivatives over torsions
 teb             bond stretch derivatives over torsions
 tea             angle bend derivatives over torsions
 teba            stretch-bend derivatives over torsions
 teub            Urey-Bradley derivatives over torsions
 teaa            angle-angle derivatives over torsions
 teopb           out-of-plane bend derivatives over torsions
 teopd           out-of-plane distance derivatives over torsions
 teid            improper dihedral derivatives over torsions
 teit            improper torsion derivatives over torsions
 tet             torsional derivatives over torsions
 tept            pi-system torsion derivatives over torsions
 tebt            stretch-torsion derivatives over torsions
 teat            angle-torsion derivatives over torsions
 tett            torsion-torsion derivatives over torsions
 tev             van der Waals derivatives over torsions
 ter             Pauli repulsion derivatives over torsions
 tedsp           dampled dispersion derivatives over torsions
 tec             charge-charge derivatives over torsions
 tecd            charge-dipole derivatives over torsions
 ted             dipole-dipole derivatives over torsions
 tem             atomic multipole derivatives over torsions
 tep             polarization derivatives over torsions
 tect            charge transfer derivatives over torsions
 terxf           reaction field derivatives over torsions
 tes             solvation derivatives over torsions
 telf            metal ligand field derivatives over torsions
 teg             geometric restraint derivatives over torsions
 tex             extra energy term derivatives over torsions

**MODULE DSPPOT        dispersion interaction scale factors

.. code-block:: text

 dsp2scale       scale factor for 1-2 dispersion energy interactions
 dsp3scale       scale factor for 1-3 dispersion energy interactions
 dsp4scale       scale factor for 1-4 dispersion energy interactions
 dsp5scale       scale factor for 1-5 dispersion energy interactions
 use_dcorr       flag to use long range dispersion correction

**MODULE ENERGI        individual potential energy components

.. code-block:: text

 esum            total potential energy of the system
 eb              bond stretch potential energy of the system
 ea              angle bend potential energy of the system
 eba             stretch-bend potential energy of the system
 eub             Urey-Bradley potential energy of the system
 eaa             angle-angle potential energy of the system
 eopb            out-of-plane bend potential energy of the system
 eopd            out-of-plane distance potential energy of the system
 eid             improper dihedral potential energy of the system
 eit             improper torsion potential energy of the system
 et              torsional potential energy of the system
 ept             pi-system torsion potential energy of the system
 ebt             stretch-torsion potential energy of the system
 eat             angle-torsion potential energy of the system
 ett             torsion-torsion potential energy of the system
 ev              van der Waals potential energy of the system
 er              Pauli repulsion potential energy of the system
 edsp            dampled dispersion potential energy of the system
 ec              charge-charge potential energy of the system
 ecd             charge-dipole potential energy of the system
 ed              dipole-dipole potential energy of the system
 em              atomic multipole potential energy of the system
 ep              polarization potential energy of the system
 ect             charge transfer potential energy of the system
 erxf            reaction field potential energy of the system
 es              solvation potential energy of the system
 elf             metal ligand field potential energy of the system
 eg              geometric restraint potential energy of the system
 ex              extra term potential energy of the system

**MODULE EWALD         Ewald summation parameters and options

.. code-block:: text

 aewald          current value of Ewald convergence coefficient
 aeewald         Ewald convergence coefficient for electrostatics
 apewald         Ewald convergence coefficient for polarization
 adewald         Ewald convergence coefficient for dispersion
 boundary        Ewald boundary condition; none, tinfoil or vacuum

**MODULE FACES         Connolly area and volume variables

.. code-block:: text

 maxcls          maximum number of neighboring atom pairs
 maxtt           maximum number of temporary tori
 maxt            maximum number of total tori
 maxp            maximum number of probe positions
 maxv            maximum number of vertices
 maxen           maximum number of concave edges
 maxfn           maximum number of concave faces
 maxc            maximum number of circles
 maxeq           maximum number of convex edges
 maxfs           maximum number of saddle faces
 maxfq           maximum number of convex faces
 maxcy           maximum number of cycles
 mxcyeq          maximum number of convex edge cycles
 mxfqcy          maximum number of convex face cycles

**MODULE FFT           Fast Fourier transform control values

.. code-block:: text

 maxprime        maximum number of prime factors of FFT dimension
 iprime          prime factorization of each FFT dimension (FFTPACK)
 planf           pointer to forward transform data structure (FFTW)
 planb           pointer to backward transform data structure (FFTW)
 ffttable        intermediate array used by the FFT routine (FFTPACK)
 ffttyp          type of FFT package; currently FFTPACK or FFTW

**MODULE FIELDS        molecular mechanics force field type

.. code-block:: text

 maxbio          maximum number of biopolymer atom definitions
 biotyp          force field atom type of each biopolymer type
 forcefield      string used to describe the current forcefield

**MODULE FILES         name & number of current structure file

.. code-block:: text

 nprior          number of previously existing cycle files
 ldir            length in characters of the directory name
 leng            length in characters of the base filename
 filename        base filename used by default for all files
 outfile         output filename used for intermediate results

**MODULE FRACS         distances to molecular center of mass

.. code-block:: text

 xfrac           fractional coordinate along a-axis of center of mass
 yfrac           fractional coordinate along b-axis of center of mass
 zfrac           fractional coordinate along c-axis of center of mass

**MODULE FREEZE        definition of holonomic constraints

.. code-block:: text

 nrat            number of holonomic distance constraints to apply
 nratx           number of atom group holonomic constraints to apply
 iratx           group number of group in a holonomic constraint
 kratx           spatial constraint type (1=plane, 2=line, 3=point)
 irat            atom numbers of atoms in a holonomic constraint
 rateps          convergence tolerance for holonomic constraints
 krat            ideal distance value for holonomic constraint
 use_rattle      logical flag to set use of holonomic contraints
 ratimage        flag to use minimum image for holonomic constraint

**MODULE GKSTUF        generalized Kirkwood solvation values

.. code-block:: text

 gkc             tuning parameter exponent in the f(GB) function
 gkr             generalized Kirkwood cavity radii for atom types

**MODULE GROUP         partitioning of system into atom groups

.. code-block:: text

 ngrp            total number of atom groups in the system
 kgrp            contiguous list of the atoms in each group
 grplist         number of the group to which each atom belongs
 igrp            first and last atom of each group in the list
 grpmass         total mass of all the atoms in each group
 wgrp            weight for each set of group-group interactions
 use_group       flag to use partitioning of system into groups
 use_intra       flag to include only intragroup interactions
 use_inter       flag to include only intergroup interactions

**MODULE HESCUT        cutoff for Hessian matrix elements

.. code-block:: text

 hesscut         magnitude of smallest allowed Hessian element

**MODULE HESSN         Cartesian Hessian elements for one atom

.. code-block:: text

 hessx           Hessian elements for x-component of current atom
 hessy           Hessian elements for y-component of current atom
 hessz           Hessian elements for z-component of current atom

**MODULE HPMF          hydrophobic potential of mean force term

.. code-block:: text

 rcarbon         radius of a carbon atom for use with HPMF
 rwater          radius of a water molecule for use with HPMF
 acsurf          surface area of a hydrophobic carbon atom
 safact          constant for calculation of atomic surface area
 tgrad           tanh slope (set very steep, default=100)
 toffset         shift the tanh plot along the x-axis (default=6)
 hpmfcut         cutoff distance for pairwise HPMF interactions
 hd1             hd2,hd3  hydrophobic PMF well depth parameter
 hc1             hc2,hc3  hydrophobic PMF well center point
 hw1             hw2,hw3  reciprocal of the hydrophobic PMF well width
 npmf            number of hydrophobic carbon atoms in the system
 ipmf            number of the atom for each HPMF carbon atom site
 rpmf            radius of each atom for use with hydrophobic PMF
 acsa            SASA value for each hydrophobic PMF carbon atom

**MODULE IELSCF        extended Lagrangian induced dipoles

.. code-block:: text

 nfree_aux       total degrees of freedom for auxiliary dipoles
 tautemp_aux     time constant for auliliary Berendsen thermostat
 kelvin_aux      target system temperature for auxiliary dipoles
 uaux            auxiliary induced dipole value at each site
 upaux           auxiliary shadow induced dipoles at each site
 vaux            auxiliary induced dipole velocity at each site
 vpaux           auxiliary shadow dipole velocity at each site
 aaux            auxiliary induced dipole acceleration at each site
 apaux           auxiliary shadow dipole acceleration at each site
 use_ielscf      flag to use inertial extended Lagrangian method

**MODULE IMPROP        improper dihedrals in current structure

.. code-block:: text

 niprop          total number of improper dihedral angles in the system
 iiprop          numbers of the atoms in each improper dihedral angle
 kprop           force constant values for improper dihedral angles
 vprop           ideal improper dihedral angle value in degrees

**MODULE IMPTOR        improper torsions in current structure

.. code-block:: text

 nitors          total number of improper torsional angles in the system
 iitors          numbers of the atoms in each improper torsional angle
 itors1          1-fold amplitude and phase for each improper torsion
 itors2          2-fold amplitude and phase for each improper torsion
 itors3          3-fold amplitude and phase for each improper torsion

**MODULE INFORM        program I/O and flow control values

.. code-block:: text

 maxask          maximum number of queries for interactive input
 digits          decimal places output for energy and coordinates
 iprint          steps between status printing (0=no printing)
 iwrite          steps between coordinate saves (0=no saves)
 isend           steps between socket communication (0=no sockets)
 silent          logical flag to turn off all information printing
 verbose         logical flag to turn on extra information printing
 debug           logical flag to turn on full debug printing
 holdup          logical flag to wait for carriage return on exit
 abort           logical flag to stop execution at next chance

**MODULE INTER         sum of intermolecular energy components

.. code-block:: text

 einter          total intermolecular potential energy

**MODULE IOUNIT        Fortran input/output unit numbers

.. code-block:: text

 input           Fortran I/O unit for main input (default=5)
 iout            Fortran I/O unit for main output (default=6)

**MODULE KANANG        angle-angle term forcefield parameters

.. code-block:: text

 anan            angle-angle cross term parameters for each atom class

**MODULE KANGS         bond angle bend forcefield parameters

.. code-block:: text

 maxna           maximum number of harmonic angle bend parameter entries
 maxna5          maximum number of 5-membered ring angle bend entries
 maxna4          maximum number of 4-membered ring angle bend entries
 maxna3          maximum number of 3-membered ring angle bend entries
 maxnap          maximum number of in-plane angle bend parameter entries
 maxnaf          maximum number of Fourier angle bend parameter entries
 acon            force constant parameters for harmonic angle bends
 acon5           force constant parameters for 5-ring angle bends
 acon4           force constant parameters for 4-ring angle bends
 acon3           force constant parameters for 3-ring angle bends
 aconp           force constant parameters for in-plane angle bends
 aconf           force constant parameters for Fourier angle bends
 ang             bond angle parameters for harmonic angle bends
 ang5            bond angle parameters for 5-ring angle bends
 ang4            bond angle parameters for 4-ring angle bends
 ang3            bond angle parameters for 3-ring angle bends
 angp            bond angle parameters for in-plane angle bends
 angf            phase shift angle and periodicity for Fourier bends
 ka              string of atom classes for harmonic angle bends
 ka5             string of atom classes for 5-ring angle bends
 ka4             string of atom classes for 4-ring angle bends
 ka3             string of atom classes for 3-ring angle bends
 kap             string of atom classes for in-plane angle bends
 kaf             string of atom classes for Fourier angle bends

**MODULE KANTOR        angle-torsion forcefield parameters

.. code-block:: text

 maxnat          maximum number of angle-torsion parameter entries
 atcon           torsional amplitude parameters for angle-torsion
 kat             string of atom classes for angle-torsion terms

**MODULE KATOMS        atom definition forcefield parameters

.. code-block:: text

 atmcls          atom class number for each of the atom types
 atmnum          atomic number for each of the atom types
 ligand          number of atoms to be attached to each atom type
 weight          average atomic mass of each atom type
 symbol          modified atomic symbol for each atom type
 describe        string identifying each of the atom types

**MODULE KBONDS        bond stretching forcefield parameters

.. code-block:: text

 maxnb           maximum number of bond stretch parameter entries
 maxnb5          maximum number of 5-membered ring bond stretch entries
 maxnb4          maximum number of 4-membered ring bond stretch entries
 maxnb3          maximum number of 3-membered ring bond stretch entries
 maxnel          maximum number of electronegativity bond corrections
 bcon            force constant parameters for harmonic bond stretch
 bcon5           force constant parameters for 5-ring bond stretch
 bcon4           force constant parameters for 4-ring bond stretch
 bcon3           force constant parameters for 3-ring bond stretch
 blen            bond length parameters for harmonic bond stretch
 blen5           bond length parameters for 5-ring bond stretch
 blen4           bond length parameters for 4-ring bond stretch
 blen3           bond length parameters for 3-ring bond stretch
 dlen            electronegativity bond length correction parameters
 kb              string of atom classes for harmonic bond stretch
 kb5             string of atom classes for 5-ring bond stretch
 kb4             string of atom classes for 4-ring bond stretch
 kb3             string of atom classes for 3-ring bond stretch
 kel             string of atom classes for electronegativity corrections

**MODULE KCHRGE        partial charge forcefield parameters

.. code-block:: text

 chg             partial charge parameters for each atom type

**MODULE KCPEN         charge penetration forcefield parameters

.. code-block:: text

 cpele           valence electron magnitude for each atom class
 cpalp           alpha charge penetration parameter for each atom class

**MODULE KCTRN         charge transfer forcefield parameters

.. code-block:: text

 ctchg           charge transfer magnitude for each atom class
 ctdmp           alpha charge transfer parameter for each atom class

**MODULE KDIPOL        bond dipole forcefield parameters

.. code-block:: text

 maxnd           maximum number of bond dipole parameter entries
 maxnd5          maximum number of 5-membered ring dipole entries
 maxnd4          maximum number of 4-membered ring dipole entries
 maxnd3          maximum number of 3-membered ring dipole entries
 dpl             dipole moment parameters for bond dipoles
 dpl5            dipole moment parameters for 5-ring dipoles
 dpl4            dipole moment parameters for 4-ring dipoles
 dpl3            dipole moment parameters for 3-ring dipoles
 pos             dipole position parameters for bond dipoles
 pos5            dipole position parameters for 5-ring dipoles
 pos4            dipole position parameters for 4-ring dipoles
 pos3            dipole position parameters for 3-ring dipoles
 kd              string of atom classes for bond dipoles
 kd5             string of atom classes for 5-ring dipoles
 kd4             string of atom classes for 4-ring dipoles
 kd3             string of atom classes for 3-ring dipoles

**MODULE KDSP          damped dispersion forcefield parameters

.. code-block:: text

 dspsix          C6 dispersion coefficient for each atom class
 dspdmp          alpha dispersion parameter for each atom class

**MODULE KEYS          contents of the keyword control file

.. code-block:: text

 maxkey          maximum number of lines in the keyword file
 nkey            number of nonblank lines in the keyword file
 keyline         contents of each individual keyword file line

**MODULE KHBOND        H-bonding term forcefield parameters

.. code-block:: text

 maxnhb          maximum number of hydrogen bonding pair entries
 radhb           radius parameter for hydrogen bonding pairs
 epshb           well depth parameter for hydrogen bonding pairs
 khb             string of atom types for hydrogen bonding pairs

**MODULE KIPROP        improper dihedral forcefield parameters

.. code-block:: text

 maxndi          maximum number of improper dihedral parameter entries
 dcon            force constant parameters for improper dihedrals
 tdi             ideal dihedral angle values for improper dihedrals
 kdi             string of atom classes for improper dihedral angles

**MODULE KITORS        improper torsion forcefield parameters

.. code-block:: text

 maxnti          maximum number of improper torsion parameter entries
 ti1             torsional parameters for improper 1-fold rotation
 ti2             torsional parameters for improper 2-fold rotation
 ti3             torsional parameters for improper 3-fold rotation
 kti             string of atom classes for improper torsional parameters

**MODULE KMULTI        atomic multipole forcefield parameters

.. code-block:: text

 maxnmp          maximum number of atomic multipole parameter entries
 multip          atomic monopole, dipole and quadrupole values
 mpaxis          type of local axis definition for atomic multipoles
 kmp             string of atom types for atomic multipoles

**MODULE KOPBND        out-of-plane bend forcefield parameters

.. code-block:: text

 maxnopb         maximum number of out-of-plane bending entries
 opbn            force constant parameters for out-of-plane bending
 kopb            string of atom classes for out-of-plane bending

**MODULE KOPDST        out-of-plane distance forcefield params

.. code-block:: text

 maxnopd         maximum number of out-of-plane distance entries
 opds            force constant parameters for out-of-plane distance
 kopd            string of atom classes for out-of-plane distance

**MODULE KORBS         pisystem orbital forcefield parameters

.. code-block:: text

 maxnpi          maximum number of pisystem bond parameter entries
 maxnpi5         maximum number of 5-membered ring pibond entries
 maxnpi4         maximum number of 4-membered ring pibond entries
 sslope          slope for bond stretch vs. pi-bond order
 sslope5         slope for 5-ring bond stretch vs. pi-bond order
 sslope4         slope for 4-ring bond stretch vs. pi-bond order
 tslope          slope for 2-fold torsion vs. pi-bond order
 tslope5         slope for 5-ring 2-fold torsion vs. pi-bond order
 tslope4         slope for 4-ring 2-fold torsion vs. pi-bond order
 electron        number of pi-electrons for each atom class
 ionize          ionization potential for each atom class
 repulse         repulsion integral value for each atom class
 kpi             string of atom classes for pisystem bonds
 kpi5            string of atom classes for 5-ring pisystem bonds
 kpi4            string of atom classes for 4-ring pisystem bonds

**MODULE KPITOR        pi-system torsion forcefield parameters

.. code-block:: text

 maxnpt          maximum number of pi-system torsion parameter entries
 ptcon           force constant parameters for pi-system torsions
 kpt             string of atom classes for pi-system torsion terms

**MODULE KPOLR         polarizability forcefield parameters

.. code-block:: text

 pgrp            connected types in polarization group of each atom type
 polr            dipole polarizability parameters for each atom type
 athl            Thole polarizability damping value for each atom type
 ddir            direct polarization damping value for each atom type

**MODULE KREPL         Pauli repulsion forcefield parameters

.. code-block:: text

 prsiz           Pauli repulsion size value for each atom class
 prdmp           alpha Pauli repulsion parameter for each atom class
 prele           number of valence electrons for each atom class

**MODULE KSTBND        stretch-bend forcefield parameters

.. code-block:: text

 maxnsb          maximum number of stretch-bend parameter entries
 stbn            force constant parameters for stretch-bend terms
 ksb             string of atom classes for stretch-bend terms

**MODULE KSTTOR        stretch-torsion forcefield parameters

.. code-block:: text

 maxnbt          maximum number of stretch-torsion parameter entries
 btcon           torsional amplitude parameters for stretch-torsion
 kbt             string of atom classes for stretch-torsion terms

**MODULE KTORSN        torsional angle forcefield parameters

.. code-block:: text

 maxnt           maximum number of torsional angle parameter entries
 maxnt5          maximum number of 5-membered ring torsion entries
 maxnt4          maximum number of 4-membered ring torsion entries
 t1              torsional parameters for standard 1-fold rotation
 t2              torsional parameters for standard 2-fold rotation
 t3              torsional parameters for standard 3-fold rotation
 t4              torsional parameters for standard 4-fold rotation
 t5              torsional parameters for standard 5-fold rotation
 t6              torsional parameters for standard 6-fold rotation
 t15             torsional parameters for 1-fold rotation in 5-ring
 t25             torsional parameters for 2-fold rotation in 5-ring
 t35             torsional parameters for 3-fold rotation in 5-ring
 t45             torsional parameters for 4-fold rotation in 5-ring
 t55             torsional parameters for 5-fold rotation in 5-ring
 t65             torsional parameters for 6-fold rotation in 5-ring
 t14             torsional parameters for 1-fold rotation in 4-ring
 t24             torsional parameters for 2-fold rotation in 4-ring
 t34             torsional parameters for 3-fold rotation in 4-ring
 t44             torsional parameters for 4-fold rotation in 4-ring
 t54             torsional parameters for 5-fold rotation in 4-ring
 t64             torsional parameters for 6-fold rotation in 4-ring
 kt              string of atom classes for torsional angles
 kt5             string of atom classes for 5-ring torsions
 kt4             string of atom classes for 4-ring torsions

**MODULE KTRTOR        torsion-torsion forcefield parameters

.. code-block:: text

 maxntt          maximum number of torsion-torsion parameter entries
 maxtgrd         maximum dimension of torsion-torsion spline grid
 maxtgrd2        maximum number of torsion-torsion spline grid points
 tnx             number of columns in torsion-torsion spline grid
 tny             number of rows in torsion-torsion spline grid
 ttx             angle values for first torsion of spline grid
 tty             angle values for second torsion of spline grid
 tbf             function values at points on spline grid
 tbx             gradient over first torsion of spline grid
 tby             gradient over second torsion of spline grid
 tbxy            Hessian cross components over spline grid
 ktt             string of torsion-torsion atom classes

**MODULE KURYBR        Urey-Bradley term forcefield parameters

.. code-block:: text

 maxnu           maximum number of Urey-Bradley parameter entries
 ucon            force constant parameters for Urey-Bradley terms
 dst13           ideal 1-3 distance parameters for Urey-Bradley terms
 ku              string of atom classes for Urey-Bradley terms

**MODULE KVDWPR        special vdw term forcefield parameters

.. code-block:: text

 maxnvp          maximum number of special van der Waals pair entries
 radpr           radius parameter for special van der Waals pairs
 epspr           well depth parameter for special van der Waals pairs
 kvpr            string of atom classes for special van der Waals pairs

**MODULE KVDWS         van der Waals term forcefield parameters

.. code-block:: text

 rad             van der Waals radius parameter for each atom type
 eps             van der Waals well depth parameter for each atom type
 rad4            van der Waals radius parameter in 1-4 interactions
 eps4            van der Waals well depth parameter in 1-4 interactions
 reduct          van der Waals reduction factor for each atom type

**MODULE LIGHT         method of lights pair neighbors indices

.. code-block:: text

 nlight          total number of sites for method of lights calculation
 kbx             low index of neighbors of each site in the x-sorted list
 kby             low index of neighbors of each site in the y-sorted list
 kbz             low index of neighbors of each site in the z-sorted list
 kex             high index of neighbors of each site in the x-sorted list
 key             high index of neighbors of each site in the y-sorted list
 kez             high index of neighbors of each site in the z-sorted list
 locx            maps the x-sorted list into original interaction list
 locy            maps the y-sorted list into original interaction list
 locz            maps the z-sorted list into original interaction list
 rgx             maps the original interaction list into x-sorted list
 rgy             maps the original interaction list into y-sorted list
 rgz             maps the original interaction list into z-sorted list

**MODULE LIMITS        interaction taper & cutoff distances

.. code-block:: text

 vdwcut          cutoff distance for van der Waals interactions
 repcut          cutoff distance for Pauli repulsion interactions
 dispcut         cutoff distance for dispersion interactions
 chgcut          cutoff distance for charge-charge interactions
 dplcut          cutoff distance for dipole-dipole interactions
 mpolecut        cutoff distance for atomic multipole interactions
 ctrncut         cutoff distance for charge transfer interactions
 vdwtaper        distance at which van der Waals switching begins
 reptaper        distance at which Pauli repulsion switching begins
 disptaper       distance at which dispersion switching begins
 chgtaper        distance at which charge-charge switching begins
 dpltaper        distance at which dipole-dipole switching begins
 mpoletaper      distance at which atomic multipole switching begins
 ctrntaper       distance at which charge transfer switching begins
 ewaldcut        cutoff distance for real space Ewald electrostatics
 dewaldcut       cutoff distance for real space Ewald dispersion
 usolvcut        cutoff distance for dipole solver preconditioner
 use_ewald       logical flag governing use of electrostatic Ewald
 use_dewald      logical flag governing use of dispersion Ewald
 use_lights      logical flag governing use of method of lights
 use_list        logical flag governing use of any neighbor lists
 use_vlist       logical flag governing use of van der Waals list
 use_dlist       logical flag governing use of dispersion list
 use_clist       logical flag governing use of charge list
 use_mlist       logical flag governing use of multipole list
 use_ulist       logical flag governing use of preconditioner list

**MODULE LINMIN        line search minimization parameters

.. code-block:: text

 stpmin          minimum step length in current line search direction
 stpmax          maximum step length in current line search direction
 cappa           stringency of line search (0=tight < cappa < 1=loose)
 slpmax          projected gradient above which stepsize is reduced
 angmax          maximum angle between search direction and -gradient
 intmax          maximum number of interpolations during line search

**MODULE MATH          mathematical and geometrical constants

.. code-block:: text

 pi              numerical value of the geometric constant
 elog            numerical value of the natural logarithm base
 radian          conversion factor from radians to degrees
 logten          numerical value of the natural log of ten
 twosix          numerical value of the sixth root of two
 sqrtpi          numerical value of the square root of Pi
 sqrttwo         numerical value of the square root of two
 sqrtthree       numerical value of the square root of three

**MODULE MDSTUF        molecular dynamics trajectory controls

.. code-block:: text

 nfree           total number of degrees of freedom for a system
 irest           steps between removal of COM motion (0=no removal)
 bmnmix          mixing coefficient for use with Beeman integrator
 arespa          inner time step for use with RESPA integrator
 dorest          logical flag to remove center of mass motion
 integrate       type of molecular dynamics integration algorithm

**MODULE MERCK         MMFF-specific force field parameters

.. code-block:: text

 nlignes         number of atom pairs having MMFF Bond Type 1
 bt_1            atom pairs having MMFF Bond Type 1
 eqclass         table of atom class equivalencies used to find
 default         parameters if explicit values are missing
 see             J. Comput. Chem., 17, 490-519, '95, Table IV)
 crd             number of attached neighbors    |
 val             valency value                   |  see T. A. Halgren,
 pilp            if 0, no lone pair              |  J. Comput. Chem.,
 if              1, one or more lone pair(s)  |  17, 616-645 (1995)
 mltb            multibond indicator             |
 arom            aromaticity indicator           |
 lin             linearity indicator             |
 sbmb            single- vs multiple-bond flag   |
 mmffarom        aromatic rings parameters
 mmffaromc       cationic aromatic rings parameters
 mmffaroma       anionic aromatic rings parameters

**MODULE MINIMA        general parameters for minimizations

.. code-block:: text

 fctmin          value below which function is deemed optimized
 hguess          initial value for the H-matrix diagonal elements
 maxiter         maximum number of iterations during optimization
 nextiter        iteration number to use for the first iteration

**MODULE MOLCUL        individual molecules in current system

.. code-block:: text

 nmol            total number of separate molecules in the system
 imol            first and last atom of each molecule in the list
 kmol            contiguous list of the atoms in each molecule
 molcule         number of the molecule to which each atom belongs
 totmass         total weight of all the molecules in the system
 molmass         molecular weight for each molecule in the system

**MODULE MOLDYN        MD trajectory velocity & acceleration

.. code-block:: text

 v               current velocity of each atom along the x,y,z-axes
 a               current acceleration of each atom along x,y,z-axes
 aalt            alternate acceleration of each atom along x,y,z-axes

**MODULE MOMENT        electric multipole moment components

.. code-block:: text

 netchg          net electric charge for the total system
 netdpl          dipole moment magnitude for the total system
 netqdp          diagonal quadrupole (Qxx, Qyy, Qzz) for system
 xdpl            dipole vector x-component in the global frame
 ydpl            dipole vector y-component in the global frame
 zdpl            dipole vector z-component in the global frame
 xxqdp           quadrupole tensor xx-component in global frame
 xyqdp           quadrupole tensor xy-component in global frame
 xzqdp           quadrupole tensor xz-component in global frame
 yxqdp           quadrupole tensor yx-component in global frame
 yyqdp           quadrupole tensor yy-component in global frame
 yzqdp           quadrupole tensor yz-component in global frame
 zxqdp           quadrupole tensor zx-component in global frame
 zyqdp           quadrupole tensor zy-component in global frame
 zzqdp           quadrupole tensor zz-component in global frame

**MODULE MPLPOT        multipole functional form details

.. code-block:: text

 m2scale         scale factor for 1-2 multipole energy interactions
 m3scale         scale factor for 1-3 multipole energy interactions
 m4scale         scale factor for 1-4 multipole energy interactions
 m5scale         scale factor for 1-5 multipole energy interactions
 use_chgpen      flag to use charge penetration damped potential
 pentyp          type of penetration damping (NONE, GORDON1, GORDON2)

**MODULE MPOLE         atomic multipoles in current structure

.. code-block:: text

 maxpole         max components (monopole=1,dipole=4,quadrupole=13)
 npole           total number of multipole sites in the system
 ipole           number of the atom for each multipole site
 polsiz          number of multipole components at each atom
 pollist         multipole site for each atom (0=no multipole)
 zaxis           number of the z-axis defining atom for each site
 xaxis           number of the x-axis defining atom for each site
 yaxis           number of the y-axis defining atom for each site
 pole            traceless Cartesian multipoles in the local frame
 rpole           traceless Cartesian multipoles in the global frame
 spole           spherical harmonic multipoles in the local frame
 srpole          spherical harmonic multipoles in the global frame
 mono0           original atomic monopole values for charge flux
 polaxe          local axis type for each multipole site

**MODULE MRECIP        reciprocal PME for permanent multipoles

.. code-block:: text

 vmxx            scalar sum xx-component of virial due to multipoles
 vmyy            scalar sum yy-component of virial due to multipoles
 vmzz            scalar sum zz-component of virial due to multipoles
 vmxy            scalar sum xy-component of virial due to multipoles
 vmxz            scalar sum xz-component of virial due to multipoles
 vmyz            scalar sum yz-component of virial due to multipoles
 cmp             Cartesian permenent multipoles as polytensor vector
 fmp             fractional permanent multipoles as polytensor vector
 cphi            Cartesian permanent multipole potential and field
 fphi            fractional permanent multipole potential and field

**MODULE MUTANT        free energy calculation hybrid atoms

.. code-block:: text

 nmut            number of atoms mutated from initial to final state
 vcouple         van der Waals lambda type (0=decouple, 1=annihilate)
 imut            atomic sites differing in initial and final state
 type0           atom type of each atom in the initial state system
 class0          atom class of each atom in the initial state system
 type1           atom type of each atom in the final state system
 class1          atom class of each atom in the final state system
 lambda          generic weighting between initial and final states
 tlambda         state weighting value for torsional potential
 vlambda         state weighting value for van der Waals potentials
 elambda         state weighting value for electrostatic potentials
 scexp           scale factor for soft core buffered 14-7 potential
 scalpha         scale factor for soft core buffered 14-7 potential
 mut             true if an atom is to be mutated, false otherwise

**MODULE NEIGH         pairwise neighbor list indices & storage

.. code-block:: text

 maxvlst         maximum size of van der Waals pair neighbor lists
 maxelst         maximum size of electrostatic pair neighbor lists
 maxulst         maximum size of dipole preconditioner pair lists
 nvlst           number of sites in list for each vdw site
 vlst            site numbers in neighbor list of each vdw site
 nelst           number of sites in list for each electrostatic site
 elst            site numbers in list of each electrostatic site
 nulst           number of sites in list for each preconditioner site
 ulst            site numbers in list of each preconditioner site
 lbuffer         width of the neighbor list buffer region
 pbuffer         width of the preconditioner list buffer region
 lbuf2           square of half the neighbor list buffer width
 pbuf2           square of half the preconditioner list buffer width
 vbuf2           square of van der Waals cutoff plus the list buffer
 vbufx           square of van der Waals cutoff plus 2X list buffer
 dbuf2           square of dispersion cutoff plus the list buffer
 dbufx           square of dispersion cutoff plus 2X list buffer
 cbuf2           square of charge cutoff plus the list buffer
 cbufx           square of charge cutoff plus 2X list buffer
 mbuf2           square of multipole cutoff plus the list buffer
 mbufx           square of multipole cutoff plus 2X list buffer
 ubuf2           square of preconditioner cutoff plus the list buffer
 ubufx           square of preconditioner cutoff plus 2X list buffer
 xvold           x-coordinate at last vdw/dispersion list update
 yvold           y-coordinate at last vdw/dispersion list update
 zvold           z-coordinate at last vdw/dispersion list update
 xeold           x-coordinate at last electrostatic list update
 yeold           y-coordinate at last electrostatic list update
 zeold           z-coordinate at last electrostatic list update
 xuold           x-coordinate at last preconditioner list update
 yuold           y-coordinate at last preconditioner list update
 zuold           z-coordinate at last preconditioner list update
 dovlst          logical flag to rebuild vdw neighbor list
 dodlst          logical flag to rebuild dispersion neighbor list
 doclst          logical flag to rebuild charge neighbor list
 domlst          logical flag to rebuild multipole neighbor list
 doulst          logical flag to rebuild preconditioner neighbor list

**MODULE NONPOL        nonpolar cavity & dispersion parameters

.. code-block:: text

 epso            water oxygen eps for implicit dispersion term
 epsh            water hydrogen eps for implicit dispersion term
 rmino           water oxygen Rmin for implicit dispersion term
 rminh           water hydrogen Rmin for implicit dispersion term
 awater          water number density at standard temp & pressure
 slevy           enthalpy-to-free energy scale factor for dispersion
 solvprs         limiting microscopic solvent pressure value
 surften         limiting macroscopic surface tension value
 spcut           starting radius for solvent pressure tapering
 spoff           cutoff radius for solvent pressure tapering
 stcut           starting radius for surface tension tapering
 stoff           cutoff radius for surface tension tapering
 rcav            atomic radius of each atom for cavitation energy
 rdisp           atomic radius of each atom for dispersion energy
 cdisp           maximum dispersion energy for each atom

**MODULE NUCLEO        parameters for nucleic acid structure

.. code-block:: text

 pucker          sugar pucker, either 2=2'-endo or 3=3'-endo
 glyco           glycosidic torsional angle for each nucleotide
 bkbone          phosphate backbone angles for each nucleotide
 dblhlx          flag to mark system as nucleic acid double helix
 deoxy           flag to mark deoxyribose or ribose sugar units
 hlxform         helix form (A, B or Z) of polynucleotide strands

**MODULE OMEGA         torsional space dihedral angle values

.. code-block:: text

 nomega          number of dihedral angles allowed to rotate
 iomega          numbers of two atoms defining rotation axis
 zline           line number in Z-matrix of each dihedral angle
 dihed           current value in radians of each dihedral angle

**MODULE OPBEND        out-of-plane bends in current structure

.. code-block:: text

 nopbend         total number of out-of-plane bends in the system
 iopb            bond angle numbers used in out-of-plane bending
 opbk            force constant values for out-of-plane bending

**MODULE OPDIST        out-of-plane distances in structure

.. code-block:: text

 nopdist         total number of out-of-plane distances in the system
 iopd            numbers of the atoms in each out-of-plane distance
 opdk            force constant values for out-of-plane distance

**MODULE OPENMP        OpenMP processor and thread values

.. code-block:: text

 nproc           number of processors available to OpenMP
 nthread         number of threads to be used with OpenMP

**MODULE ORBITS        conjugated pisystem orbital energies

.. code-block:: text

 qorb            number of pi-electrons contributed by each atom
 worb            ionization potential of each pisystem atom
 emorb           repulsion integral for each pisystem atom

**MODULE OUTPUT        output file format control parameters

.. code-block:: text

 archive         logical flag to save structures in an archive
 noversion       logical flag governing use of filename versions
 overwrite       logical flag to overwrite intermediate files inplace
 cyclesave       logical flag to mark use of numbered cycle files
 velsave         logical flag to save velocity vector components
 frcsave         logical flag to save force vector components
 uindsave        logical flag to save induced atomic dipoles
 coordtype       selects Cartesian, internal, rigid body or none

**MODULE PARAMS        force field parameter file contents

.. code-block:: text

 maxprm          maximum number of lines in the parameter file
 nprm            number of nonblank lines in the parameter file
 prmline         contents of each individual parameter file line

**MODULE PATHS         Elber reaction path method parameters

.. code-block:: text

 pnorm           length of the reactant-product vector
 acoeff          transformation matrix 'A' from Elber algorithm
 pc0             reactant Cartesian coordinates as variables
 pc1             product Cartesian coordinates as variables
 pvect           vector connecting the reactant and product
 pstep           step per cycle along reactant-product vector
 pzet            current projection on reactant-product vector
 gc              gradient of the path constraints

**MODULE PBSTUF        Poisson-Boltzmann solvation parameters

.. code-block:: text

 APBS            configuration parameters (see APBS documentation for details)
 In              the column on the right are possible values for each variable,
 with            default values given in brackets. Only a subset of the APBS
 options         are supported and/or are appropriate for use with AMOEBA
 pbtyp           lpbe
 At              some point AMOEBA with the non-linear PBE could be supported,
 but             this is only worked out for energies (no gradients)
 pbsoln          mg-auto, [mg-manual]
 Currently       there is only limited support for focusing calculations,
 which           is a powerful feature of APBS. At present, all energies and
 forces          must all be calculated using the finest solution
 bcfl            boundary conditions              zero, sdh, [mdh]
 chgm            multipole discretization         spl4
 other           charge discretization methods are not appropriate for AMOEBA
 srfm            surface method                   mol, smol, [spl4]
 spl4            is required for forces calculations, although mol is useful
 for             comparison with generalized Kirkwood
 dime            number of grid points            [65, 65, 65]
 grid            grid spacing (mg-manual)         fxn of "dime"
 cgrid           coarse grid spacing              fxn of "dime"
 fgrid           fine grid spacing                cgrid / 2
 stable          results require grid spacing to be fine enough to keep
 multipoles      inside the dielectric boundary (2.5 * grid < PBR)
 gcent           grid center (mg-manual)          center of mass
 cgcent          coarse grid center               center of mass
 fgcent          fine grid center                 center of mass
 pdie            solute/homogeneous dieletric     [1.0]
 sdie            solvent dieletric                [78.3]
 ionn            number of ion species            [0]
 ionc            ion concentration (M)            [0.0]
 ionq            ion charge (electrons)           [1.0]
 ionr            ion radius (A)                   [2.0]
 srad            solvent probe radius (A)         [1.4]
 swin            surface spline window width      [0.3]
 sdens           density of surface points        [10.0]
 additional      parameter to facilitate default grid setup
 smin            minimum distance between an      [10.0]
 atom            and the grid boundary (A)
 pbe             Poisson-Boltzmann permanent multipole solvation energy
 apbe            Poisson-Boltzmann permanent multipole energy over atoms
 pbr             Poisson-Boltzmann cavity radii for atom types
 pbep            Poisson-Boltzmann energies on permanent multipoles
 pbfp            Poisson-Boltzmann forces on permanent multipoles
 pbtp            Poisson-Boltzmann torques on permanent multipoles
 pbeuind         Poisson-Boltzmann field due to induced dipoles
 pbeuinp         Poisson-Boltzmann field due to non-local induced dipoles

**MODULE PDB           Protein Data Bank structure definition

.. code-block:: text

 npdb            number of atoms stored in Protein Data Bank format
 nres            number of residues stored in Protein Data Bank format
 resnum          number of the residue to which each atom belongs
 resatm          number of first and last atom in each residue
 npdb12          number of atoms directly bonded to each CONECT atom
 ipdb12          atom numbers of atoms connected to each CONECT atom
 pdblist         list of the Protein Data Bank atom number of each atom
 xpdb            x-coordinate of each atom stored in PDB format
 ypdb            y-coordinate of each atom stored in PDB format
 zpdb            z-coordinate of each atom stored in PDB format
 altsym          string with PDB alternate locations to be included
 pdbres          Protein Data Bank residue name assigned to each atom
 pdbatm          Protein Data Bank atom name assigned to each atom
 pdbtyp          Protein Data Bank record type assigned to each atom
 chnsym          string with PDB chain identifiers to be included
 instyp          string with PDB insertion records to be included

**MODULE PHIPSI        phi-psi-omega-chi angles for protein

.. code-block:: text

 chiral          chirality of each amino acid residue (1=L, -1=D)
 disulf          residue joined to each residue via a disulfide link
 phi             value of the phi angle for each amino acid residue
 psi             value of the psi angle for each amino acid residue
 omg             value of the omega angle for each amino acid residue
 chi             values of the chi angles for each amino acid residue

**MODULE PIORBS        conjugated system in current structure

.. code-block:: text

 norbit          total number of pisystem orbitals in the system
 nconj           total number of separate conjugated piystems
 reorbit         number of evaluations between orbital updates
 nbpi            total number of bonds affected by the pisystem
 ntpi            total number of torsions affected by the pisystem
 iorbit          numbers of the atoms containing pisystem orbitals
 iconj           first and last atom of each pisystem in the list
 kconj           contiguous list of atoms in each pisystem
 piperp          atoms defining a normal plane to each orbital
 ibpi            bond and piatom numbers for each pisystem bond
 itpi            torsion and pibond numbers for each pisystem torsion
 pbpl            pi-bond orders for bonds in "planar" pisystem
 pnpl            pi-bond orders for bonds in "nonplanar" pisystem
 listpi          atom list indicating whether each atom has an orbital

**MODULE PISTUF        bond order-related pisystem parameters

.. code-block:: text

 bkpi            bond stretch force constants for pi-bond order of 1.0
 blpi            ideal bond length values for a pi-bond order of 1.0
 kslope          rate of force constant decrease with bond order decrease
 lslope          rate of bond length increase with a bond order decrease
 torsp2          2-fold torsional energy barrier for pi-bond order of 1.0

**MODULE PITORS        pi-system torsions in current structure

.. code-block:: text

 npitors         total number of pi-system torsional interactions
 ipit            numbers of the atoms in each pi-system torsion
 kpit            2-fold pi-system torsional force constants

**MODULE PME           values for particle mesh Ewald summation

.. code-block:: text

 nfft1           current number of PME grid points along a-axis
 nfft2           current number of PME grid points along b-axis
 nfft3           current number of PME grid points along c-axis
 nefft1          number of grid points along electrostatic a-axis
 nefft2          number of grid points along electrostatic b-axis
 nefft3          number of grid points along electrostatic c-axis
 ndfft1          number of grid points along dispersion a-axis
 ndfft2          number of grid points along dispersion b-axis
 ndfft3          number of grid points along dispersion c-axis
 bsorder         current order of the PME B-spline values
 bseorder        order of the electrostatic PME B-spline values
 bsporder        order of the polarization PME B-spline values
 bsdorder        order of the dispersion PME B-spline values
 igrid           initial Ewald grid values for B-spline
 bsmod1          B-spline moduli along the a-axis direction
 bsmod2          B-spline moduli along the b-axis direction
 bsmod3          B-spline moduli along the c-axis direction
 bsbuild         B-spline derivative coefficient temporary storage
 thetai1         B-spline coefficients along the a-axis
 thetai2         B-spline coefficients along the b-axis
 thetai3         B-spline coefficients along the c-axis
 qgrid           values on the particle mesh Ewald grid
 qfac            prefactors for the particle mesh Ewald grid

**MODULE POLAR         induced dipole moments & polarizability

.. code-block:: text

 npolar          total number of polarizable sites in the system
 ipolar          number of the multipole for each polarizable site
 polarity        dipole polarizability for each multipole site (Ang**3)
 thole           Thole polarizability damping value for each site
 dirdamp         direct polarization damping value for each site
 pdamp           value of polarizability scale factor for each site
 udir            direct induced dipole components at each multipole site
 udirp           direct induced dipoles in field used for energy terms
 udirs           direct GK or PB induced dipoles at each multipole site
 udirps          direct induced dipoles in field used for GK or PB energy
 uind            mutual induced dipole components at each multipole site
 uinp            mutual induced dipoles in field used for energy terms
 uinds           mutual GK or PB induced dipoles at each multipole site
 uinps           mutual induced dipoles in field used for GK or PB energy
 uexact          exact SCF induced dipoles to full numerical precision
 douind          flag to allow induced dipoles at each atomic site

**MODULE POLGRP        polarization group connectivity lists

.. code-block:: text

 maxp11          maximum number of atoms in a polarization group
 maxp12          maximum number of atoms in groups 1-2 to an atom
 maxp13          maximum number of atoms in groups 1-3 to an atom
 maxp14          maximum number of atoms in groups 1-4 to an atom
 np11            number of atoms in polarization group of each atom
 np12            number of atoms in groups 1-2 to each atom
 np13            number of atoms in groups 1-3 to each atom
 np14            number of atoms in groups 1-4 to each atom
 ip11            atom numbers of atoms in same group as each atom
 ip12            atom numbers of atoms in groups 1-2 to each atom
 ip13            atom numbers of atoms in groups 1-3 to each atom
 ip14            atom numbers of atoms in groups 1-4 to each atom

**MODULE POLOPT        induced dipoles for OPT extrapolation

.. code-block:: text

 maxopt          maximum order for OPT induced dipole extrapolation
 optorder        highest coefficient order for OPT dipole extrapolation
 optlevel        current OPT order for reciprocal potential and field
 copt            coefficients for OPT total induced dipole moments
 copm            coefficients for OPT incremental induced dipole moments
 uopt            OPT induced dipole components at each multipole site
 uoptp           OPT induced dipoles in field used for energy terms
 uopts           OPT GK or PB induced dipoles at each multipole site
 uoptps          OPT induced dipoles in field used for GK or PB energy
 fopt            OPT fractional reciprocal potentials at multipole sites
 foptp           OPT fractional reciprocal potentials for energy terms

**MODULE POLPCG        induced dipoles via the PCG solver

.. code-block:: text

 mindex          index into preconditioner inverse for PCG solver
 pcgpeek         value of acceleration factor for PCG peek step
 minv            preconditioner inverse for induced dipole PCG solver
 pcgprec         flag to allow use of preconditioner with PCG solver
 pcgguess        flag to use initial PCG based on direct field

**MODULE POLPOT        polarization functional form details

.. code-block:: text

 politer         maximum number of induced dipole SCF iterations
 poleps          induced dipole convergence criterion (rms Debye/atom)
 p2scale         scale factor for 1-2 polarization energy interactions
 p3scale         scale factor for 1-3 polarization energy interactions
 p4scale         scale factor for 1-4 polarization energy interactions
 p5scale         scale factor for 1-5 polarization energy interactions
 p2iscale        scale factor for 1-2 intragroup polarization energy
 p3iscale        scale factor for 1-3 intragroup polarization energy
 p4iscale        scale factor for 1-4 intragroup polarization energy
 p5iscale        scale factor for 1-5 intragroup polarization energy
 d1scale         scale factor for intra-group direct induction
 d2scale         scale factor for 1-2 group direct induction
 d3scale         scale factor for 1-3 group direct induction
 d4scale         scale factor for 1-4 group direct induction
 u1scale         scale factor for intra-group mutual induction
 u2scale         scale factor for 1-2 group mutual induction
 u3scale         scale factor for 1-3 group mutual induction
 u4scale         scale factor for 1-4 group mutual induction
 w2scale         scale factor for 1-2 induced dipole interactions
 w3scale         scale factor for 1-3 induced dipole interactions
 w4scale         scale factor for 1-4 induced dipole interactions
 w5scale         scale factor for 1-5 induced dipole interactions
 udiag           acceleration factor for induced dipole SCF iterations
 dpequal         flag to set dscale values equal to pscale values
 use_thole       flag to use Thole damped polarization interactions
 use_dirdamp     flag to use damped direct polarization interactions
 poltyp          type of polarization (MUTUAL, DIRECT, OPT or TCG)

**MODULE POLTCG        induced dipoles via the TCG solver

.. code-block:: text

 tcgorder        total number of TCG iterations to be used
 tcgnab          number of mutual induced dipole components
 tcgpeek         value of acceleration factor for TCG peek step
 uad             left-hand side mutual induced d-dipoles
 uap             left-hand side mutual induced p-dipoles
 ubd             right-hand side mutual induced d-dipoles
 ubp             right-hand side mutual induced p-dipoles
 tcgguess        flag to use initial TCG based on direct field

**MODULE POTENT        usage of potential energy components

.. code-block:: text

 use_bond        logical flag governing use of bond stretch potential
 use_angle       logical flag governing use of angle bend potential
 use_strbnd      logical flag governing use of stretch-bend potential
 use_urey        logical flag governing use of Urey-Bradley potential
 use_angang      logical flag governing use of angle-angle cross term
 use_opbend      logical flag governing use of out-of-plane bend term
 use_opdist      logical flag governing use of out-of-plane distance
 use_improp      logical flag governing use of improper dihedral term
 use_imptor      logical flag governing use of improper torsion term
 use_tors        logical flag governing use of torsional potential
 use_pitors      logical flag governing use of pi-system torsion term
 use_strtor      logical flag governing use of stretch-torsion term
 use_angtor      logical flag governing use of angle-torsion term
 use_tortor      logical flag governing use of torsion-torsion term
 use_vdw         logical flag governing use of vdw der Waals potential
 use_repuls      logical flag governing use of Pauli repulsion term
 use_disp        logical flag governing use of dispersion potential
 use_charge      logical flag governing use of charge-charge potential
 use_chgdpl      logical flag governing use of charge-dipole potential
 use_dipole      logical flag governing use of dipole-dipole potential
 use_mpole       logical flag governing use of multipole potential
 use_polar       logical flag governing use of polarization term
 use_chgtrn      logical flag governing use of charge transfer term
 use_chgflx      logical flag governing use of charge flux term
 use_rxnfld      logical flag governing use of reaction field term
 use_solv        logical flag governing use of continuum solvation term
 use_metal       logical flag governing use of ligand field term
 use_geom        logical flag governing use of geometric restraints
 use_extra       logical flag governing use of extra potential term
 use_born        logical flag governing use of Born radii values
 use_orbit       logical flag governing use of pisystem computation

**MODULE POTFIT        values for electrostatic potential fit

.. code-block:: text

 nconf           total number of configurations to be analyzed
 namax           maximum number of atoms in the largest configuration
 ngatm           total number of atoms with active potential grid points
 nfatm           total number of atoms in electrostatic potential fit
 npgrid          total number of electrostatic potential grid points
 ipgrid          atom associated with each potential grid point
 resp            weight used to restrain parameters to original values
 xdpl0           target x-component of the molecular dipole moment
 ydpl0           target y-component of the molecular dipole moment
 zdpl0           target z-component of the molecular dipole moment
 xxqdp0          target xx-component of the molecular quadrupole moment
 xyqdp0          target xy-component of the molecular quadrupole moment
 xzqdp0          target xz-component of the molecular quadrupole moment
 yyqdp0          target yy-component of the molecular quadrupole moment
 yzqdp0          target yz-component of the molecular quadrupole moment
 zzqdp0          target zz-component of the molecular quadrupole moment
 fit0            initial value of each parameter used in potential fit
 fchg            partial charges by atom type during potential fit
 fpol            atomic multipoles by atom type during potential fit
 pgrid           Cartesian coordinates of potential grid points
 epot            values of electrostatic potential at grid points
 use_dpl         flag to include molecular dipole in potential fit
 use_qdp         flag to include molecular quadrupole in potential fit
 fit_mpl         flag for atomic monopoles to vary in potential fit
 fit_dpl         flag for atomic dipoles to vary in potential fit
 fit_qdp         flag for atomic quadrupoles to vary in potential fit
 fitchg          flag marking atom types for use in partial charge fit
 fitpol          flag marking atom types for use in atomic multipole fit
 gatm            flag to use potential grid points around each atom
 fatm            flag to use each atom in electrostatic potential fit

**MODULE PTABLE        symbols and info for chemical elements

.. code-block:: text

 maxele          maximum number of elements from periodic table
 atmass          standard atomic weight for each chemical element
 vdwrad          van der Waals radius for each chemical element
 covrad          covalent radius for each chemical element
 elemnt          atomic symbol for each chemical element

**MODULE REFER         reference atomic coordinate storage

.. code-block:: text

 nref            total number of atoms in each reference system
 refltitle       length in characters of each reference title line
 refleng         length in characters of each reference filename
 reftyp          atom types of the atoms in each reference system
 n12ref          number of atoms bonded to each reference atom
 i12ref          atom numbers of atoms 1-2 connected to each atom
 xboxref         reference a-axis length of periodic box
 yboxref         reference b-axis length of periodic box
 zboxref         reference c-axis length of periodic box
 alpharef        reference angle between b- and c-axes of box
 betaref         reference angle between a- and c-axes of box
 gammaref        reference angle between a- and b-axes of box
 xref            reference x-coordinates for atoms in each system
 yref            reference y-coordinates for atoms in each system
 zref            reference z-coordinates for atoms in each system
 refnam          atom names of the atoms in each reference system
 reffile         base filename for each reference system
 reftitle        title used to describe each reference system

**MODULE REPEL         Pauli repulsion for current structure

.. code-block:: text

 nrep            total number of repulsion sites in the system
 sizpr           Pauli repulsion size parameter value at each site
 dmppr           Pauli repulsion alpha damping value at each site
 elepr           Pauli repulsion valence electrons at each site

**MODULE REPPOT        repulsion interaction scale factors

.. code-block:: text

 r2scale         scale factor for 1-2 repulsion energy interactions
 r3scale         scale factor for 1-3 repulsion energy interactions
 r4scale         scale factor for 1-4 repulsion energy interactions
 r5scale         scale factor for 1-5 repulsion energy interactions

**MODULE RESDUE        amino acid & nucleotide residue names

.. code-block:: text

 maxamino        maximum number of amino acid residue types
 maxnuc          maximum number of nucleic acid residue types
 ntyp            biotypes for mid-chain peptide backbone N atoms
 catyp           biotypes for mid-chain peptide backbone CA atoms
 ctyp            biotypes for mid-chain peptide backbone C atoms
 hntyp           biotypes for mid-chain peptide backbone HN atoms
 otyp            biotypes for mid-chain peptide backbone O atoms
 hatyp           biotypes for mid-chain peptide backbone HA atoms
 cbtyp           biotypes for mid-chain peptide backbone CB atoms
 nntyp           biotypes for N-terminal peptide backbone N atoms
 cantyp          biotypes for N-terminal peptide backbone CA atoms
 cntyp           biotypes for N-terminal peptide backbone C atoms
 hnntyp          biotypes for N-terminal peptide backbone HN atoms
 ontyp           biotypes for N-terminal peptide backbone O atoms
 hantyp          biotypes for N-terminal peptide backbone HA atoms
 nctyp           biotypes for C-terminal peptide backbone N atoms
 cactyp          biotypes for C-terminal peptide backbone CA atoms
 cctyp           biotypes for C-terminal peptide backbone C atoms
 hnctyp          biotypes for C-terminal peptide backbone HN atoms
 octyp           biotypes for C-terminal peptide backbone O atoms
 hactyp          biotypes for C-terminal peptide backbone HA atoms
 o5typ           biotypes for nucleotide backbone and sugar O5' atoms
 c5typ           biotypes for nucleotide backbone and sugar C5' atoms
 h51typ          biotypes for nucleotide backbone and sugar H5' atoms
 h52typ          biotypes for nucleotide backbone and sugar H5'' atoms
 c4typ           biotypes for nucleotide backbone and sugar C4' atoms
 h4typ           biotypes for nucleotide backbone and sugar H4' atoms
 o4typ           biotypes for nucleotide backbone and sugar O4' atoms
 c1typ           biotypes for nucleotide backbone and sugar C1' atoms
 h1typ           biotypes for nucleotide backbone and sugar H1' atoms
 c3typ           biotypes for nucleotide backbone and sugar C3' atoms
 h3typ           biotypes for nucleotide backbone and sugar H3' atoms
 c2typ           biotypes for nucleotide backbone and sugar C2' atoms
 h21typ          biotypes for nucleotide backbone and sugar H2' atoms
 o2typ           biotypes for nucleotide backbone and sugar O2' atoms
 h22typ          biotypes for nucleotide backbone and sugar H2'' atoms
 o3typ           biotypes for nucleotide backbone and sugar O3' atoms
 ptyp            biotypes for nucleotide backbone and sugar P atoms
 optyp           biotypes for nucleotide backbone and sugar OP atoms
 h5ttyp          biotypes for nucleotide backbone and sugar H5T atoms
 h3ttyp          biotypes for nucleotide backbone and sugar H3T atoms
 amino           three-letter abbreviations for amino acids types
 nuclz           three-letter abbreviations for nucleic acids types
 amino1          one-letter abbreviations for amino acids types
 nuclz1          one-letter abbreviations for nucleic acids types

**MODULE RESTRN        parameters for geometrical restraints

.. code-block:: text

 npfix           number of position restraints to be applied
 ndfix           number of distance restraints to be applied
 nafix           number of angle restraints to be applied
 ntfix           number of torsional restraints to be applied
 ngfix           number of group distance restraints to be applied
 nchir           number of chirality restraints to be applied
 ipfix           atom number involved in each position restraint
 kpfix           flags to use x-, y-, z-coordinate position restraints
 idfix           atom numbers defining each distance restraint
 iafix           atom numbers defining each angle restraint
 itfix           atom numbers defining each torsional restraint
 igfix           group numbers defining each group distance restraint
 ichir           atom numbers defining each chirality restraint
 depth           depth of shallow Gaussian basin restraint
 width           exponential width coefficient of Gaussian basin
 rwall           radius of spherical droplet boundary restraint
 xpfix           x-coordinate target for each restrained position
 ypfix           y-coordinate target for each restrained position
 zpfix           z-coordinate target for each restrained position
 pfix            force constant and flat-well range for each position
 dfix            force constant and target range for each distance
 afix            force constant and target range for each angle
 tfix            force constant and target range for each torsion
 gfix            force constant and target range for each group distance
 chir            force constant and target range for chiral centers
 use_basin       logical flag governing use of Gaussian basin
 use_wall        logical flag governing use of droplet boundary

**MODULE RGDDYN        rigid body MD velocities and momenta

.. code-block:: text

 xcmo            x-component from each atom to center of rigid body
 ycmo            y-component from each atom to center of rigid body
 zcmo            z-component from each atom to center of rigid body
 vcm             current translational velocity of each rigid body
 wcm             current angular velocity of each rigid body
 lm              current angular momentum of each rigid body
 vc              half-step translational velocity for kinetic energy
 wc              half-step angular velocity for kinetic energy
 linear          logical flag to mark group as linear or nonlinear

**MODULE RIGID         rigid body coordinates for atom groups

.. code-block:: text

 xrb             rigid body reference x-coordinate for each atom
 yrb             rigid body reference y-coordinate for each atom
 zrb             rigid body reference z-coordinate for each atom
 rbc             current rigid body coordinates for each group
 use_rigid       flag to mark use of rigid body coordinate system

**MODULE RING          number and location of ring structures

.. code-block:: text

 nring3          total number of 3-membered rings in the system
 nring4          total number of 4-membered rings in the system
 nring5          total number of 5-membered rings in the system
 nring6          total number of 6-membered rings in the system
 nring7          total number of 7-membered rings in the system
 iring3          numbers of the atoms involved in each 3-ring
 iring4          numbers of the atoms involved in each 4-ring
 iring5          numbers of the atoms involved in each 5-ring
 iring6          numbers of the atoms involved in each 6-ring
 iring7          numbers of the atoms involved in each 7-ring

**MODULE ROTBND        molecule partitions for bond rotation

.. code-block:: text

 nrot            total number of atoms moving when bond rotates
 rot             atom numbers of atoms moving when bond rotates
 use_short       logical flag governing use of shortest atom list

**MODULE RXNFLD        reaction field matrix and indices

.. code-block:: text

 ijk             indices into the reaction field element arrays
 b1              first reaction field matrix element array
 b2              second reaction field matrix element array

**MODULE RXNPOT        reaction field functional form details

.. code-block:: text

 rfsize          radius of reaction field sphere centered at origin
 rfbulkd         bulk dielectric constant of reaction field continuum
 rfterms         number of terms to use in reaction field summation

**MODULE SCALES        optimization parameter scale factors

.. code-block:: text

 scale           multiplicative factor for each optimization parameter
 set_scale       logical flag to show if scale factors have been set

**MODULE SEQUEN        sequence information for biopolymer

.. code-block:: text

 nseq            total number of residues in biopolymer sequences
 nchain          number of separate biopolymer sequence chains
 ichain          first and last residue in each biopolymer chain
 seqtyp          residue type for each residue in the sequence
 seq             three-letter code for each residue in the sequence
 chnnam          one-letter identifier for each sequence chain
 chntyp          contents of each chain (GENERIC, PEPTIDE or NUCLEIC)

**MODULE SHUNT         polynomial switching function values

.. code-block:: text

 off             distance at which the potential energy goes to zero
 off2            square of distance at which the potential goes to zero
 cut             distance at which switching of the potential begins
 cut2            square of distance at which the switching begins
 c0              zeroth order coefficient of multiplicative switch
 c1              first order coefficient of multiplicative switch
 c2              second order coefficient of multiplicative switch
 c3              third order coefficient of multiplicative switch
 c4              fourth order coefficient of multiplicative switch
 c5              fifth order coefficient of multiplicative switch
 f0              zeroth order coefficient of additive switch function
 f1              first order coefficient of additive switch function
 f2              second order coefficient of additive switch function
 f3              third order coefficient of additive switch function
 f4              fourth order coefficient of additive switch function
 f5              fifth order coefficient of additive switch function
 f6              sixth order coefficient of additive switch function
 f7              seventh order coefficient of additive switch function

**MODULE SIZES         parameters to set array dimensions

"sizes" sets values for critical array dimensions used
throughout the software; these parameters fix the size of
the largest systems that can be handled

.. code-block:: text

 parameter       maximum allowed number of:
 maxatm          atoms in the molecular system
 maxtyp          force field atom type definitions
 maxclass        force field atom class definitions
 maxval          atoms directly bonded to an atom
 maxref          stored reference molecular systems
 maxgrp          user-defined groups of atoms
 maxres          residues in the macromolecule
 maxfix          geometric constraints and restraints

**MODULE SOCKET        socket communication control parameters

.. code-block:: text

 skttyp          socket information type (1=DYN, 2=OPT)
 cstep           current dynamics or optimization step number
 cdt             current dynamics cumulative simulation time
 cenergy         current potential energy from simulation
 sktstart        logical flag to indicate socket initialization
 sktstop         logical flag to indicate socket shutdown
 use_socket      logical flag governing use of external sockets

**MODULE SOLUTE        continuum solvation model parameters

.. code-block:: text

 doffset         dielectric offset to continuum solvation atomic radii
 p1              single-atom scale factor for analytical Still radii
 p2              1-2 interaction scale factor for analytical Still radii
 p3              1-3 interaction scale factor for analytical Still radii
 p4              nonbonded scale factor for analytical Still radii
 p5              soft cutoff parameter for analytical Still radii
 rsolv           atomic radius of each atom for continuum solvation
 asolv           atomic surface area solvation parameters
 rborn           Born radius of each atom for GB/SA solvation
 drb             solvation derivatives with respect to Born radii
 drbp            GK polarization derivatives with respect to Born radii
 drobc           chain rule term for Onufriev-Bashford-Case radii
 gpol            polarization self-energy values for each atom
 shct            overlap scale factors for Hawkins-Cramer-Truhlar radii
 aobc            alpha values for Onufriev-Bashford-Case radii
 bobc            beta values for Onufriev-Bashford-Case radii
 gobc            gamma values for Onufriev-Bashford-Case radii
 vsolv           atomic volume of each atom for use with ACE
 wace            "omega" values for atom class pairs for use with ACE
 s2ace           "sigma^2" values for atom class pairs for use with ACE
 uace            "mu" values for atom class pairs for use with ACE
 solvtyp         type of continuum solvation energy model in use
 borntyp         method to be used for the Born radius computation

**MODULE STODYN        SD trajectory frictional coefficients

.. code-block:: text

 friction        global frictional coefficient for exposed particle
 fgamma          atomic frictional coefficients for each atom
 use_sdarea      logical flag to use surface area friction scaling

**MODULE STRBND        stretch-bends in current structure

.. code-block:: text

 nstrbnd         total number of stretch-bend interactions
 isb             angle and bond numbers used in stretch-bend
 sbk             force constants for stretch-bend terms

**MODULE STRTOR        stretch-torsions in current structure

.. code-block:: text

 nstrtor         total number of stretch-torsion interactions
 ist             torsion and bond numbers used in stretch-torsion
 kst             1-, 2- and 3-fold stretch-torsion force constants

**MODULE SYNTRN        synchronous transit path definition

.. code-block:: text

 tpath           value of the path coordinate (0=reactant, 1=product)
 ppath           path coordinate for extra point in quadratic transit
 xmin1           reactant coordinates as array of optimization variables
 xmin2           product coordinates as array of optimization variables
 xm              extra coordinate set for quadratic synchronous transit

**MODULE TARRAY        store dipole-dipole matrix elements

.. code-block:: text

 ntpair          number of stored dipole-dipole matrix elements
 tindex          index into stored dipole-dipole matrix values
 tdipdip         stored dipole-dipole matrix element values

**MODULE TITLES        title for current molecular system

.. code-block:: text

 ltitle          length in characters of the nonblank title string
 title           title used to describe the current structure

**MODULE TORPOT        torsional functional form details

.. code-block:: text

 idihunit        convert improper dihedral energy to kcal/mole
 itorunit        convert improper torsion amplitudes to kcal/mole
 torsunit        convert torsional parameter amplitudes to kcal/mole
 ptorunit        convert pi-system torsion energy to kcal/mole
 storunit        convert stretch-torsion energy to kcal/mole
 atorunit        convert angle-torsion energy to kcal/mole
 ttorunit        convert torsion-torsion energy to kcal/mole

**MODULE TORS          torsional angles in current structure

.. code-block:: text

 ntors           total number of torsional angles in the system
 itors           numbers of the atoms in each torsional angle
 tors1           1-fold amplitude and phase for each torsional angle
 tors2           2-fold amplitude and phase for each torsional angle
 tors3           3-fold amplitude and phase for each torsional angle
 tors4           4-fold amplitude and phase for each torsional angle
 tors5           5-fold amplitude and phase for each torsional angle
 tors6           6-fold amplitude and phase for each torsional angle

**MODULE TORTOR        torsion-torsions in current structure

.. code-block:: text

 ntortor         total number of torsion-torsion interactions
 itt             atoms and parameter indices for torsion-torsion

**MODULE TREE          potential smoothing search tree levels

.. code-block:: text

 maxpss          maximum number of potential smoothing levels
 nlevel          number of levels of potential smoothing used
 etree           energy reference value at the top of the tree
 ilevel          smoothing deformation value at each tree level

**MODULE UNITS         physical constants and unit conversions

D. B. Newell, F. Cabiati, J. Fischer, K. Fujii, S. G. Karshenboim,
S. Margolis, E. de Mirandes, P. J. Mohr, F. Nez, K. Pachucki,
T. J. Quinn, N. Taylor, M. Wang, B. M. Wood and Z. Zhang, "The
CODATA 2017 Values of h, e, k, and Na for the Revision of the SI",
Metrologia, 55, L13-L16 (2018)

P. J. Mohr, D. B. Newell and B. N. Taylor, "CODATA Recommended
Values of the Fundamental Physical Constants: 2014", Journal of
Physical and Chemical Reference Data, 45, 043102 (2016)

Where available, values are from the 2017 CODATA adjustment
based on exact physical constants for the revised SI

Other values are from the 2014 CODATA reference constants; also
available online from the National Institute of Standards and
Technology at http://physics.nist.gov/cuu/Constants/index.html/

The conversion from calorie to Joule is the definition of the
thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992)

The "coulomb" energy conversion factor is found by dimensional
analysis of Coulomb's Law, ie, by dividing the square of the
elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is
the permittivity of vacuum (the "electric constant"); note that
eps0 is typically given in F/m, equivalent to C**2/(J-m)

The approximate value used for the Debye, 3.33564 x 10-30 C-m,
is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997)

The value of "prescon" is based on definition of 1 atmosphere
as 101325 Pa set by the 10th Conference Generale des Poids et
Mesures (1954), where a Pascal (Pa) is equal to a J/m**3

.. code-block:: text

 avogadro        Avogadro's number (N) in particles/mole
 lightspd        speed of light in vacuum (c) in cm/ps
 boltzmann       Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K
 gasconst        ideal gas constant (R) in kcal/mole/K
 elemchg         elementary charge of a proton in Coulombs
 vacperm         vacuum permittivity (electric constant, eps0) in F/m
 emass           mass of an electron in atomic mass units
 planck          Planck's constant (h) in J-s
 joule           conversion from calorie to joule
 ekcal           conversion from kcal to g*Ang**2/ps**2
 bohr            conversion from Bohr to Angstrom
 hartree         conversion from Hartree to kcal/mole
 evolt           conversion from Hartree to electron-volt
 efreq           conversion from Hartree to cm-1
 coulomb         conversion from electron**2/Ang to kcal/mole
 debye           conversion from electron-Ang to Debye
 prescon         conversion from kcal/mole/Ang**3 to Atm

**MODULE UPRIOR        previous values of induced dipoles

.. code-block:: text

 maxpred         maximum number of predictor induced dipoles to save
 nualt           number of sets of prior induced dipoles in storage
 maxualt         number of sets of induced dipoles needed for predictor
 gear            coefficients for Gear predictor binomial method
 aspc            coefficients for always stable predictor-corrector
 bpred           coefficients for induced dipole predictor polynomial
 bpredp          coefficients for predictor polynomial in energy field
 bpreds          coefficients for predictor for PB/GK solvation
 bpredps         coefficients for predictor in PB/GK energy field
 udalt           prior values for induced dipoles at each site
 upalt           prior values for induced dipoles in energy field
 usalt           prior values for induced dipoles for PB/GK solvation
 upsalt          prior values for induced dipoles in PB/GK energy field
 use_pred        flag to control use of induced dipole prediction
 polpred         type of predictor polynomial (GEAR, ASPC or LSQR)

**MODULE UREY          Urey-Bradley interactions in structure

.. code-block:: text

 nurey           total number of Urey-Bradley terms in the system
 iury            numbers of the atoms in each Urey-Bradley interaction
 uk              Urey-Bradley force constants (kcal/mole/Ang**2)
 ul              ideal 1-3 distance values in Angstroms

**MODULE URYPOT        Urey-Bradley functional form details

.. code-block:: text

 cury            cubic coefficient in Urey-Bradley potential
 qury            quartic coefficient in Urey-Bradley potential
 ureyunit        convert Urey-Bradley energy to kcal/mole

**MODULE USAGE         atoms active during energy computation

.. code-block:: text

 nuse            total number of active atoms in energy calculation
 iuse            numbers of the atoms active in energy calculation
 use             true if an atom is active, false if inactive

**MODULE VALFIT        valence term parameter fitting values

.. code-block:: text

 fit_bond        logical flag to fit bond stretch parameters
 fit_angle       logical flag to fit angle bend parameters
 fit_strbnd      logical flag to fit stretch-bend parameters
 fit_urey        logical flag to fit Urey-Bradley parameters
 fit_opbend      logical flag to fit out-of-plane bend parameters
 fit_tors        logical flag to fit torsional parameters
 fit_struct      logical flag to structure-fit valence parameters
 fit_force       logical flag to force-fit valence parameters

**MODULE VDW           van der Waals terms in current structure

.. code-block:: text

 nvdw            total number van der Waals active sites in the system
 ivdw            number of the atom for each van der Waals active site
 jvdw            type or class index into vdw parameters for each atom
 ired            attached atom from which reduction factor is applied
 kred            value of reduction factor parameter for each atom
 xred            reduced x-coordinate for each atom in the system
 yred            reduced y-coordinate for each atom in the system
 zred            reduced z-coordinate for each atom in the system
 radmin          minimum energy distance for each atom class pair
 epsilon         well depth parameter for each atom class pair
 radmin4         minimum energy distance for 1-4 interaction pairs
 epsilon4        well depth parameter for 1-4 interaction pairs
 radhbnd         minimum energy distance for hydrogen bonding pairs
 epshbnd         well depth parameter for hydrogen bonding pairs

**MODULE VDWPOT        van der Waals functional form details

.. code-block:: text

 igauss          coefficients of Gaussian fit to vdw potential
 ngauss          number of Gaussians used in fit to vdw potential
 abuck           value of "A" constant in Buckingham vdw potential
 bbuck           value of "B" constant in Buckingham vdw potential
 cbuck           value of "C" constant in Buckingham vdw potential
 ghal            value of "gamma" in buffered 14-7 vdw potential
 dhal            value of "delta" in buffered 14-7 vdw potential
 v2scale         factor by which 1-2 vdw interactions are scaled
 v3scale         factor by which 1-3 vdw interactions are scaled
 v4scale         factor by which 1-4 vdw interactions are scaled
 v5scale         factor by which 1-5 vdw interactions are scaled
 use_vcorr       flag to use long range van der Waals correction
 vdwindex        indexing mode (atom type or class) for vdw parameters
 vdwtyp          type of van der Waals potential energy function
 radtyp          type of parameter (sigma or R-min) for atomic size
 radsiz          atomic size provided as radius or diameter
 radrule         combining rule for atomic size parameters
 epsrule         combining rule for vdw well depth parameters
 gausstyp        type of Gaussian fit to van der Waals potential

**MODULE VIBS          iterative vibrational analysis components

.. code-block:: text

 rho             trial vectors for iterative vibrational analysis
 rhok            alternate vectors for iterative vibrational analysis
 rwork           temporary work array for eigenvector transformation

**MODULE VIRIAL        components of internal virial tensor

.. code-block:: text

 vir             total internal virial Cartesian tensor components
 use_virial      logical flag governing use of virial computation

**MODULE WARP          potential surface smoothing parameters

.. code-block:: text

 deform          value of the smoothing deformation parameter
 difft           diffusion coefficient for torsional potential
 diffv           diffusion coefficient for van der Waals potential
 diffc           diffusion coefficient for charge-charge potential
 m2              second moment of the GDA gaussian for each atom
 use_smooth      flag to use a potential energy smoothing method
 use_dem         flag to use diffusion equation method potential
 use_gda         flag to use gaussian density annealing potential
 use_tophat      flag to use analytical tophat smoothed potential
 use_stophat     flag to use shifted tophat smoothed potential

**MODULE XTALS         structures used for parameter fitting

.. code-block:: text

 maxlsq          maximum number of least squares variables
 maxrsd          maximum number of residual functions
 nxtal           number of molecular structures to be stored
 nvary           number of potential parameters to optimize
 ivary           index for the types of potential parameters
 iresid          structure to which each residual function refers
 vary            atom numbers involved in potential parameters
 e0_lattice      ideal lattice energy for the current crystal
 vartyp          type of each potential parameter to be optimized
 rsdtyp          experimental variable for each of the residuals

**MODULE ZCLOSE        Z-matrix ring openings and closures

.. code-block:: text

 nadd            number of added bonds between Z-matrix atoms
 ndel            number of bonds between Z-matrix bonds to delete
 iadd            numbers of the atom pairs defining added bonds
 idel            numbers of the atom pairs defining deleted bonds

**MODULE ZCOORD        Z-matrix internal coordinate values

.. code-block:: text

 iz              defining atom numbers for each Z-matrix atom
 zbond           bond length used to define each Z-matrix atom
 zang            bond angle used to define each Z-matrix atom
 ztors           angle or torsion used to define Z-matrix atom
