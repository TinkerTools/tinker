Acknowledgments
===============

The Tinker package has developed over a period of many years, very slowly during the late-1980s, and more rapidly since the mid-1990s in Jay Ponder's research group at the Washington University School of Medicine in Saint Louis. Many people have played significant roles in the development of the package into its current form. The major contributors are listed below:

**Stew Rubenstein**

coordinate interconversions; original optimization methods and torsional angle manipulation

**Craig Kundrot**

molecular surface area & volume and their derivatives

**Shawn Huston**

original Amber/OPLS implementation; free energy calculations; time correlation functions

**Mike Dudek**

initial multipole models for peptides and proteins

**Yong "Mike" Kong**

multipole electrostatics; dipole polarization; reaction field treatment; TINKER water model

**Reece Hart**

potential smoothing methodology; Scheraga's DEM, Straub's GDA and extensions

**Mike Hodsdon**

extension of the Tinker distance geometry program and its application to NMR NOE structure determination

**Rohit Pappu**

potential smoothing methodology and PSS algorithms; rigid body optimization; GB/SA solvation derivatives

**Wijnand Mooij**

MM3 directional hydrogen bonding term; crystal lattice minimization code

**Gerald Loeffler**

stochastic/Langevin dynamics implementation

**Marina Vorobieva & Nina Sokolova**

nucleic acid building module and parameter translation

**Peter Bagossi**

AMOEBA force field parameters for alkanes and diatomics

**Pengyu Ren**

Ewald summation for polarizable atomic multipoles; AMOEBA force field for water, organics and peptides

**Anders Carlsson**

original ligand field potential energy term for transition metals

**Andrey Kutepov**

integrator for rigid-body dynamics trajectories

**Tom Darden**

Particle Mesh Ewald (PME) code, and development of PME for the AMOEBA force field

**Alan Grossfield**

Monte Carlo minimization; tophat potential smoothing

**Michael Schnieders**

Force Field Explorer GUI for Tinker; neighbor lists for nonbonded interactions

**Chuanjie Wu**

solvation free energy calculations; AMOEBA nucleic acid force field; parameterization tools for Tinker

**Justin Xiang**

angular overlap and valence bond potential models for transition metals

**David Gohara**

OpenMP parallelization of energy terms including PME, and parallel neighbor lists

**Chao Lu**

derivatives of potential energy with respect to lambda for metadynamics and similar methods

**Aaron Gordon**

enthalpy and entropy estimates as an adjunct to BAR free energy calculation

**Zhi Wang**

Bennett acceptance ratio (BAR) for free energy calculations

**Josh Rackers & Rose Silva**

implementation of the HIPPO force field for water and general organic molecules

It is critically important that Tinker's distributed force field parameter sets exactly reproduce the intent of the original force field authors. We would like to thank Julian Tirado-Rives (OPLS-AA), Alex MacKerell (CHARMM27), Wilfred van Gunsteren (GROMOS), and Adrian Roitberg and Carlos Simmerling (AMBER) for their help in testing Tinker's results against those given by the authentic programs and parameter sets. Lou Allinger provided updated parameters for MM2 and MM3 on several occasions. His very successful methods provided the original inspiration for the development of Tinker.

Still other workers have devoted considerable time in developing code that will hopefully be incorporated into future Tinker versions; for example, Jim Kress (UFF implementation) and Michael Sheets (numerous code optimizations, thermodynamic integration). Finally, we wish to thank the many users of the Tinker package for their suggestions and comments, praise and criticism, which have resulted in a variety of improvements.
