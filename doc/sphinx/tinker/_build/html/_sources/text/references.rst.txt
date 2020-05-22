References
==========

This section contains a list of the references to general theory, algorithms and implementation details which have been of use during the development of the TINKER package. Methods described in some of the references have been implemented in detail within the TINKER source code. Other references contain useful background information although the algorithms themselves are now obsolete. Still other papers contain ideas or extensions planned for future inclusion in TINKER. References for specific force field parameter sets are provided in an earlier section of this User's Guide. This list is heavily skewed toward biomolecules in general and proteins in particular. This bias reflects our group's major interests; however an attempt has been made to include methods which should be generally applicable.

Partial List of Molecular Mechanics Software Packages
-----------------------------------------------------

.. code-block:: text

 AMBER          Peter Kollman, University of California, San Francisco
 AMMP           Rob Harrison, Thomas Jefferson University, Philadelphia
 ARGOS          Andy McCammon, University of California, San Diego
 BOSS           William Jorgensen, Yale University
 BRUGEL         Shoshona Wodak, Free University of Brussels
 CFF            Shneior Lifson, Weizmann Institute
 CHARMM         Martin Karplus, Harvard University
 CHARMM/GEMM    Bernard Brooks, National Institutes of Health, Bethesda
 DELPHI         Bastian van de Graaf, Delft University of Technology
 DISCOVER       Molecular Simulations Inc., San Diego
 DL_POLY        W. Smith & T. Forester, CCP5, Daresbury Laboratory
 ECEPP          Harold Scheraga, Cornell University
 ENCAD          Michael Levitt, Stanford University
 FANTOM         Werner Braun, University of Texas, Galveston
 FEDER/2        Nobuhiro Go, Kyoto University
 GROMACS        Herman Berendsen, University of Groningen
 GROMOS         Wilfred van Gunsteren, BIOMOS and ETH, Zurich
 IMPACT         Ronald Levy, Rutgers University
 MACROMODEL     Schodinger, Inc., Jersey City, New Jersey
 MM2/MM3/MM4    N. Lou Allinger, University of Georgia
 MMC            Cliff Dykstra, Indiana Univ.-Purdue Univ. at Indianapolis
 MMFF           Tom Halgren, Merck Research Laboratories, Rahway
 MMTK           Konrad Hinsen, Inst. of Structural Biology, Grenoble
 MOIL           Ron Elber, Cornell University
 MOLARIS        Arieh Warshal, University of Southern California
 MOLDY          Keith Refson, Oxford University
 MOSCITO        Dietmar Paschek & Alfons Geiger, Universit‰t Dortmund
 NAMD           Klaus Schulten, University of Illinois, Urbana
 OOMPAA         Andy McCammon, University of California, San Diego
 ORAL           Karel Zimmerman, INRA, Jouy-en-Josas, France
 ORIENT         Anthony Stone, Cambridge University
 PCMODEL        Kevin Gilbert, Serena Software, Bloomington, Indiana
 PEFF           Jan Dillen, University of Pretoria, South Africa
 Q              Johan Aqvist, Uppsala University
 SIBFA          Nohad Gresh, INSERM, CNRS, Paris
 SIGMA          Jan Hermans, University of North Carolina
 SPASIBA        Gerard Vergoten, UniversitÈ de Lille
 SPASMS         David Spellmeyer and the Kollman Group, UCSF
 TINKER         Jay Ponder, Washington University, St. Louis
 XPLOR/CNS      Axel Brunger, Stanford University
 YAMMP          Stephen Harvey, University of Alabama, Birmingham
 YASP           Florian Mueller-Plathe, ETH Zentrum, Zurich
 YETI           Angelo Vedani, Biografik-Labor 3R, Basel

**AMBER**     D. A Pearlman, D. A. Case, J. W. Caldwell, W. S. Ross, T. E. Cheatham III, S. DeBolt, D. Ferguson, G. Seibel and P. Kollman, AMBER, a Package of Computer Programs for Applying Molecular Mechanics, Normal Mode Analysis, Molecular Dynamics and Free Energy Calculations to Simulate the Structural and Energetic Properties of Molecules, Comp. Phys. Commun., 91, 1-41 (1995)

**ARGOS**     T. P. Straatsma and J. A. McCammon, ARGOS, a Vectorized General Molecular Dynamics Program, J. Comput. Chem., 11, 943-951 (1990)

**CHARMM**     B. R. Brooks, R. E. Bruccoleri, B. D. Olafson, D. J. States, S. Swaminathan and M. Karplus, CHARMM: A Program for Macromolecular Energy, Minimization, and Dynamics Calculations, J. Comput. Chem., 4, 187-217 (1983)

**ENCAD**     M. Levitt, M. Hirshberg, R. Sharon and V. Daggett, Potential Energy Function and Parameters for Simulations for the Molecular Dynamics of Proteins and Nucleic Acids in Solution, Comp. Phys. Commun., 91, 215-231 (1995)

**FANTOM**     T. Schaumann, W. Braun and K. Wurtrich, The Program FANTOM for Energy Refinement of Polypeptides and Proteins Using a Newton-Raphson Minimizer in Torsion Angle Space, Biopolymers, 29, 679-694 (1990)

**FEDER/2**     H. Wako, S. Endo, K. Nagayama and N. Go, FEDER/2: Program for Static and Dynamic Conformational Energy Analysis of Macro-molecules in Dihedral Angle Space, Comp. Phys. Commun., 91, 233-251 (1995)

**GROMACS**     E. Lindahl, B. Hess and D. van der Spoel, GROMACS 3.0: A Package for Molecular Simulation and Trajectory Analysis, J. Mol. Mol., 7, 306-317 (2001)

**GROMOS**     W. R. P. Scott, P. H. Hunenberger , I. G. Tironi, A. E. Mark, S. R. Billeter, J. Fennen, A. E. Torda, T. Huber, P. Kruger, W. F. van Gunsteren, The GROMOS Biomolecular Simulation Program Package, J. Phys. Chem. A, 103, 3596-3607 (1999)

**IMPACT**     D. B. Kitchen, F. Hirata, J. D. Westbrook, R. Levy, D. Kofke and M. Yarmush, Conserving Energy during Molecular Dynamics Simulations of Water, Proteins, and Proteins in Water, J. Comput. Chem., 10, 1169-1180 (1990)

**MACROMODEL**     F. Mahamadi, N. G. J. Richards, W. C. Guida, R. Liskamp, M. Lipton, C. Caufield, G. Chang, T. Hendrickson and W. C. Still, MacroModel: An Integrated Software System for Modeling Organic and Bioorganic Molecules Using Molecular Mechanics, J. Comput. Chem., 11, 440-467 (1990)

**MM2**     N. L. Allinger, Conformational Analysis. 130. MM2. A Hydrocarbon Force Field Utilizing V1 and V2 Torsional Terms, J. Am. Chem. Soc., 99, 8127-8134 (1977)

**MM3**     N. L. Allinger, Y. H. Yuh and J.-H. Lii, Molecular Mechanics. The MM3 Force Field for Hydrocarbons, J. Am. Chem. Soc., 111, 8551-8566 (1989)

**MM4**     N. L. Allinger, K. Chen and J.-H. Lii, An Improved Force Field (MM4) for Saturated Hydrocarbons, J. Comput. Chem., 17, 642-668 (1996)

**MMC**     C. E. Dykstra, Molecular Mechanics for Weakly Interacting Assemblies of Rare Gas Atoms and Small Molecules, J. Am. Chem. Soc., 111, 6168-6174 (1989)

**MMFF**     T. A. Halgren, Merck Molecular Force Field. I. Basis, Form, Scope, Parameterization, and Performance of MMFF94, J. Comput. Chem., 17, 490-516 (1996)

**MOIL**     R. Elber, A. Roitberg, C. Simmerling, R. Goldstein, H. Li, G. Verkhiver, C. Keasar, J. Zhang and A. Ulitsky, MOIL: A Program for Simulations of Macromolecules, Comp. Phys. Commun., 91, 159-189 (1995)

**MOSCITO**     See the web site at http:/ganter.chemie.uni-dortmund.de/~pas/moscito.html

**NAMD**     L. KalÈ, R. Skeel, M. Bhandarkar, R. Brunner, A. Gursoy, N. Krawetz, J. Phillips, A. Shinozaki, K. Varadarajan and K. Schulten, NAMD2: Greater Scalability for Parallel Molecular Dynamics, J. Comput. Phys., 151, 283-312 (1999)

**OOMPAA**     G. A. Huber and J. A. McCammon, OOMPAA: Object-oriented Model for Probing Assemblages of Atoms, J. Comput. Phys., 151, 264-282 (1999)

**ORAL**     K. Zimmermann, ORAL: All Purpose Molecular Mechanics Simulator and Energy Minimizer, J. Comput. Chem., 12, 310-319 (1991)

**PCMODEL**     See the web site at http:/www.serenasoft.com

**PEFF**     J. L. M. Dillen, PEFF: A Program for the Development of Empirical Force Fields, J. Comput. Chem., 13, 257-267 (1992)

**Q**     See the web site at http://aqvist.bmc.uu.se/Q

**SIBFA**     N. Gresh, Inter- and Intramolecular Interactions. Inception and Refinements of the SIBFA, Molecular Mechanics (SMM) Procedure, a Separable, Polarizable Methodology Grounded on ab Initio SCF/MP2 Computations. Examples of Applications to Molecular Recognition Problems, J. Chim. Phys. PCB, 94, 1365-1416 (1997)

**SIGMA**     See the web site at http://femto.med.unc.edu/SIGMA

**SPASIBA**     P. Derreumaux and G. Vergoten, A New Spectroscopic Molecular Mechanics Force-Field - Parameters For Proteins, J. Chem. Phys., 102, 8586-8605 (1995)

**TINKER**     See the web site at http://dasher.wustl.edu/tinker

**YAMMP**     R. K.-Z. Tan and S. C. Harvey, Yammp: Development of a Molecular Mechanics Program Using the Modular Programming Method, J. Comput. Chem., 14, 455-470 (1993)

**YETI**     A. Vedani, YETI: An Interactive Molecular Mechanics Program for Small-Molecule Protein Complexes, J. Comput. Chem., 9, 269-280 (1988)

Molecular Mechanics
-------------------

U. Burkert and N. L. Allinger, Molecular Mechanics, American Chemical Society, Washington, D.C., 1982

P. Comba and T. W. Hambley, Molecular Modeling of Inorganic Compounds, 2nd Ed., Wiley-VCH, New York, 2001

K. Machida, Principles of Molecular Mechanics, Kodansha/John Wiley & Sons, Tokyo/New York, 1999

A. K. Rappe and C. J. Casewit, Molecular Mechanics across Chemistry, University Science Books, Sausalito, CA, 1997

K. Rasmussen, Potential Energy Functions in Conformational Analysis (Lecture Notes in Chemistry, Vol. 27), Springer-Verlag, Berlin, 1985

Computer Simulation Methods
---------------------------

M. P. Allen and D. J. Tildesley, Computer Simulation of Liquids, Oxford University Press, Oxford, 1987

C. J. Cramer, Essentials of Computational Chemistry: Theories and Models, John Wiley and Sons, New York, 2002

M. J. Field, A Practical Introduction to the Simulation of Molecular Systems, Cambridge Univ. Press, Cambridge, 1999

D. Frankel and B. Smit, Understanding Molecular Simulation: From Algorithms to Applications, 2nd Ed., Academic Press, San Diego, CA, 2001

J. M. Haile, Molecular Dynamics Simulation: Elementary Methods, John Wiley and Sons, New York, 1992

F. Jensen, Introduction to Computational Chemistry, John Wiley and Sons, New York, 1998

A. R. Leach, Molecular Modelling: Principles and Applications, 2nd Ed., Addison Wesley Longman, Essex, England, 2001

D. C. Rapaport, The Art of Molecular Dynamics Simulation, 2nd Ed., Cambridge University Press, Cambridge, 2004

T. Schlick, Molecular Modeling and Simulation, Springer-Verlag, New York, 2002

Modeling of Biological Macromolecules
-------------------------------------

O. M. Becker, A. D. MacKerell, Jr., B. Roux and M. Watanabe, Eds., Computational Biochemistry and Biophysics, Marcel Dekker, New York, 2001

C. L. Brooks III, M. Karplus and B. M. Pettitt, Proteins: A Theoretical Perspective of Dynamics, Structure, and Thermodynamics, John Wiley and Sons, New York, 1988

V. Daggett, Ed., Protein Simulations (Advances in Protein Chemistry, Vol. 66), Academic Press/Elsevier, New York, 2003

J. A. McCammon and S. Harvey, Dynamics of Proteins and Nucleic Acids, Cambridge University Press, Cambridge, 1987

W. F. van Gunsteren, P. K. Weiner and A. J. Wilkinson, Computer Simulation of Biomolecular Systems, Vol. 1-3, Kluwer Academic Publishers, Dordrecht, 1989-1997

Conjugate Gradient and Quasi-Newton Optimization
------------------------------------------------

J. Nocedal and S. J. Wright, Numerical Optimization, Springer-Verlag, New York, 1999

S. G. Nash and A. Sofer, Linear and Nonlinear Programming, McGraw-Hill, New York, 1996

R. Fletcher, Practical Methods of Optimization, John Wiley & Sons Ltd., Chichester, 1987

D. G. Luenberger, Linear and Nonlinear Programming, 2nd Ed., Addison-Wesley, Reading, MA, 1984

P. E. Gill, W. Murray and M. H. Wright, Practical Optimization, Academic Press, New York, 1981

J. Nocedal, Updating Quasi-Newton Matrices with Limited Storage, Math. Comp., 773-782 (1980)

S. J. Watowich, E. S. Meyer, R. Hagstrom and R. Josephs, A Stable, Rapidly Converging Conjugate Gradient Method for Energy Minimization, J. Comput. Chem., 9, 650-661 (1988)

W. C. Davidon, Optimally Conditioned Optimization Algorithms without Line Searches, Math. Prog., 9, 1-30 (1975)

Truncated Newton Optimization
-----------------------------

J. W. Ponder and F. M. Richards, An Efficient Newton-like Method for Molecular Mechanics Energy Minimization of Large Molecules, J. Comput. Chem., 8, 1016-1024 (1987)

R. S. Dembo and T. Steihaug, Truncated-Newton Algorithms for Large-Scale Unconstrained Optimization, Math. Prog., 26, 190-212 (1983)

S. C. Eisenstat and H. F. Walker, Choosing the Forcing Terms in an Inexact Newton Method, SIAM J. Sci. Comput., 17, 16-32 (1996)

T. Schlick and M. Overton, A Powerful Truncated Newton Method for Potential Energy Minimization, J. Comput. Chem., 8, 1025-1039 (1987)

D. S. Kershaw, The Incomplete Cholesky-Conjugate Gradient Method for the Iterative Solution of Systems of Linear Equations, J. Comput. Phys., 26, 43-65 (1978)

T. A. Manteuffel, An Incomplete Factorization Technique for Positive Definite Linear Systems, Math. Comp., 34, 473-497 (1980)

P. Derreumaux, G. Zhang and T. Schlick and B. R. Brooks, A Truncated Newton Minimizer Adapted for CHARMM and Biomolecular Applications, J. Comput. Chem., 15, 532-552 (1994)

I. S. Duff, A. M. Erisman and J. K. Reid, Direct Methods for Sparse Matrices, Oxford University Press, Oxford, 1986

Potential Energy Smoothing
--------------------------

R. V. Pappu, R. K. Hart and J. W. Ponder, Analysis and Application of Potential Energy Smoothing Methods for Global Optimization, J. Phys. Chem. B, 102, 9725-9742 (1998)

L. Piela, J. Kostrowicki and H. A. Scheraga, The Multiple-Minima Problem in the Conformational Analysis of Molecules. Deformation of the Potential Energy Hypersurface by the Diffusion Equation Method, J. Phys. Chem., 93, 3339-3346 (1989)

J. Ma and J. E. Straub, Simulated Annealing Using the Classical Density Distribution, J. Chem. Phys., 101, 533-541 (1994)

C. Tsoo and C. L. Brooks, Cluster Structure Determination Using Gaussian Density Distribution Global Minimization Methods, J. Chem. Phys., 101, 6405-6411 (1994)

S. Nakamura, H. Hirose, M. Ikeguchi and J. Doi, Conformational Energy Minimization Using a Two-Stage Method, J. Phys. Chem., 99, 8374-8378 (1995)

T. Huber, A. E. Torda and W. F. van Gunsteren, Structure Optimization Combining Soft-Core Interaction Functions, the Diffusion Equation Method, and Molecular Dynamics, J. Phys. Chem. A, 101, 5926-5930 (1997)

S. Schelstraete and H. Verschelde, Finding Minimum-Energy Configurations of Lennard-Jones Clusters Using an Effective Potential, J. Phys. Chem. A, 101, 310-315 (1998)

I. Andricioaei and J. E. Straub, Global Optimization Using Bad Derivatives: Derivative-Free Method for Molecular Energy Minimization, J. Comput. Chem., 19, 1445-1455 (1998)

L. Piela, Search for the Most Stable Structures on Potential Energy Surfaces, Coll. Czech. Chem. Commun., 63, 1368-1380 (1998)

"Sniffer" Global Optimization
-----------------------------

A. O. Griewank, Generalized Descent for Global Optimization, J. Opt. Theor. Appl., 34, 11-39 (1981)

R. A. R. Butler and E. E. Slaminka, An Evaluation of the Sniffer Global Optimization Algorithm Using Standard Test Functions, J. Comput. Phys., 99, 28-32 (1993)

J. W. Rogers and R. A. Donnelly, Potential Transformation Methods for Large-Scale Global Optimization, SIAM J. Optim., 5, 871-891 (1995)

Integration Methods for Molecular Dynamics
------------------------------------------

D. Beeman, Some Multistep Methods for Use in Molecular Dynamics Calculations, J. Comput. Phys., 20, 130-139 (1976)

M. Levitt and H. Meirovitch, Integrating the Equations of Motion, J. Mol. Biol., 168, 617-620 (1983)

J. Aqvist, W. F. van Gunsteren, M. Leijonmarck and O. Tapia, A Molecular Dynamics Study of the C-Terminal Fragment of the L7/L12 Ribosomal Protein, J. Mol. Biol., 183, 461-477 (1985)

W. C. Swope, H. C. Andersen, P. H. Berens and K. R. Wilson, A Computer Simulation Method for the Calculation of Equilibrium Constants for the Formation of Physical Clusters of Molecules: Application to Small Water Clusters, J. Chem. Phys., 76, 637-649 (1982)

Constraint Dynamics
-------------------

W. F. van Gunsteren and H. J. C. Berendsen, Algorithms for Macromolecular Dynamics and Constraint Dynamics, Mol. Phys., 34, 1311-1327 (1977)

G. Ciccotti, M. Ferrario and J.-P. Ryckaert, Molecular Dynamics of Rigid Systems in Cartesian Coordinates: A General Formulation, Mol. Phys., 47, 1253-1264 (1982)

H. C. Andersen, Rattle: A "Velocity" Version of the Shake Algorithm for Molecular Dynamics Calculations, J. Comput. Phys., 52, 24-34 (1983)

R. Kutteh, RATTLE Recipe for General Holonomic Constraints: Angle and Torsion Constraints, CCP5 Newsletter, 46, 9-17 (1998) [available from the web site at http://www.dl.ac.uk/CCP/CCP5/newsletter_index.html]

B. J. Palmer, Direct Application of SHAKE to the Velocity Verlet Algorithm, J. Comput. Phys., 104, 470-472 (1993)

S. Miyamoto and P. A. Kollman, SETTLE: An Analytical Version of the SHAKE and RATTLE Algorithm for Rigid Water Models, J. Comput. Chem., 13, 952-962 (1992)

B. Hess, H. Bekker, H. J. C. Berendsen and J. G. E. M. Fraaije, LINCS: A Linear Constraint Solver for Molecular Simulations, J. Comput. Chem., 18, 1463-1472 (1997)

J. T. Slusher and P. T. Cummings, Non-Iterative Constraint Dynamics using Velocity-Explicit Verlet Methods, Mol. Simul., 18, 213-224 (1996)

Langevin, Brownian and Stochastic Dynamics
------------------------------------------

M. P. Allen, Brownian Dynamics Simulation of a Chemical Reaction in Solution, Mol. Phys., 40, 1073-1087 (1980)

W. F. van Gunsteren and H. J. C. Berendsen, Algorithms for Brownian Dynamics, Mol. Phys., 45, 637-647 (1982)

F. Guarnieri and W. C. Still, A Rapidly Convergent Simulation Method: Mixed Monte Carlo/Stochastic Dynamics, J. Comput. Chem., 15, 1302-1310 (1994)

M. G. Paterlini and D. M. Ferguson, Constant Temperature Simulations using the Langevin Equation with Velocity Verlet Integration, Chem. Phys., 236, 243-252 (1998)

Constant Temperature and Pressure Dynamics
------------------------------------------

H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren, A. DiNola and J. R. Haak, Molecular Dynamics with Coupling to an External Bath, J. Chem. Phys., 81, 3684-3690 (1984)

W. G. Hoover, Canonical Dynamics: Equilibrium Phase-space Distributions, Phys. Rev. A, 31, 1695-1697 (1985)

J. J. Morales, S. Toxvaerd and L. F. Rull, Computer Simulation of a Phase Transition at Constant Temperature and Pressure, Phys. Rev. A, 34, 1495-1498 (1986)

B. R. Brooks, Algorithms for Molecular Dynamics at Constant Temperature and Pressure, Internal Report of Division of Computer Research and Technology, National Institutes of Health, 1988.

M. Levitt, Molecular Dynamics of Native Protein: Computer Simulation of Trajectories, J. Mol. Biol., 168, 595-620 (1983)

Out-of-Plane Deformation Terms
------------------------------

J. R. Maple, U. Dinar and A. T. Hagler, Derivation of Force Fields for Molecular Mechanics and Dynamics from ab initio Energy Surfaces, Proc. Natl. Acad. Sci. USA, 85, 5350-5354 (1988)

S.-H. Lee, K. Palmo and S. Krimm, New Out-of-Plane Angle and Bond Angle Internal Coordinates and Related Potential Energy Functions for Molecular Mechanics and Dynamics Simulations, J. Comput. Chem., 20, 1067-1084 (1999)

Analytical Derivatives of Potential Functions
---------------------------------------------

K. J. Miller, R. J. Hinde and J. Anderson, First and Second Derivative Matrix Elements for the Stretching, Bending, and Torsional Energy, J. Comput. Chem., 10, 63-76 (1989)

D. H. Faber and C. Altona, UTAH5: A Versatile Programme Package for the Calculation of Molecular Properties by Force Field Methods, Computers & Chemistry, 1, 203-213 (1977)

W. C. Swope and D. M. Ferguson, Alternative Expressions for Energies and Forces Due to Angle Bending and Torsional Energy, Report G320-3561, J. Comput. Chem., 13, 585-594 (1992)

A. Blondel and M. Karplus, New Formulation for Derivatives of Torsion Angles and Improper Torsion Angles in Molecular Mechanics: Elimination of Singularities, J. Comput. Chem., 17, 1132-1141 (1996)

R. E. Tuzun, D. W. Noid and B. G. Sumpter, Efficient Treatment of Out-of-Plane Bend and Improper Torsion Interactions in MM2, MM3, and MM4 Molecular Mechanics Calculations, J. Comput. Chem., 18, 1804-1811 (1997)

Torsional Space Derivatives and Normal Modes
--------------------------------------------

M. Levitt, C. Sander and P. S. Stern, Protein Normal-mode Dynamics:  Trypsin Inhibitor, Crambin, Ribonuclease and Lysozyme, J. Mol. Biol., 181, 423-447 (1985)

M. Levitt, Protein Folding by Restrained Energy Minimization and Molecular Dynamics, J. Mol. Biol., 170, 723-764 (1983)

H. Wako and N. Go, Algorithm for Rapid Calculation of Hessian of Conformational Energy Function of Proteins by Supercomputer, J. Comput. Chem., 8, 625-635 (1987)

H. Abe, W. Braun, T. Noguti and N. Go, Rapid Calculation of First and Second Derivatives of Conformational Energy with Respect to Dihedral Angles for Proteins: General Recurrent Equations, Computers & Chemistry, 8, 239-247 (1984)

T. Noguti and N. Go, A Method of Rapid Calculation of a Second Derivative Matrix of Conformational Energy for Large Molecules, J. Phys. Soc. Japan, 52, 3685-3690 (1983)

Analytical Surface Area and Volume
----------------------------------

M. L. Connolly, Analytical Molecular Surface Calculation, J. Appl. Cryst., 16, 548-558 (1983)

M. L. Connolly, Computation of Molecular Volume, J. Am. Chem. Soc., 107, 1118-1124 (1985)

M. L. Connolly, Molecular Surfaces: A Review, available from the web site at http://www.netsci.org/Science/Compchem/feature14.html

C. E. Kundrot, J. W. Ponder and F. M. Richards, Algorithms for Calculating Excluded Volume and Its Derivatives as a Function of Molecular Conformation and Their Use in Energy Minimization, J. Comput. Chem., 12, 402-409 (1991)

T. J. Richmond, Solvent Accessible Surface Area and Excluded Volume in Proteins, J. Mol. Biol., 178, 63-89 (1984)

L. Wesson and D. Eisenberg, Atomic Solvation Parameters Applied to Molecular Dynamics of Proteins in Solution, Protein Science, 1, 227-235 (1992)

V. Gononea and E. Osawa, Implementation of Solvent Effect in Molecular Mechanics, Part 3. The First- and Second-order Analytical Derivatives of Excluded Volume, J. Mol. Struct. (Theochem), 311 305-324 (1994)

K. D. Gibson and H. A. Scheraga, Exact Calculation of the Volume and Surface Area of Fused Hard-sphere Molecules with Unequal Atomic Radii, Mol. Phys., 62, 1247-1265 (1987)

K. D. Gibson and H. A. Scheraga, Surface Area of the Intersection of Three Spheres with Unequal Radii: A Simplified Analytical Formula, Mol. Phys., 64, 641-644 (1988)

S. Sridharan, A. Nichols and K. A. Sharp, A Rapid Method for Calculating Derivatives of Solvent Accessible Surface Areas of Molecules, J. Comput, Chem., 16, 1038-1044 (1995)

Approximate Surface Area and Volume
-----------------------------------

S. J. Wodak and J. Janin, Analytical Approximation to the Accessible Surface Area of Proteins, Proc. Natl. Acad. Sci. USA, 77, 1736-1740 (1980)

W. Hasel, T. F. Hendrickson and W. C. Still, A Rapid Approximation to the Solvent Accessible Surface Areas of Atoms, Tetrahedron Comput. Method., 1, 103-116 (1988)

J. Weiser, P. S. Shenkin and W. C. Still, Approximate Solvent-Accessible Surface Areas from Tetrahedrally Directed Neighber Densities, Biopolymers, 50, 373-380 (1999)

Boundary Conditions and Neighbor Methods
----------------------------------------

W. F. van Gunsteren, H. J. C. Berendsen, F. Colonna, D. Perahia, J. P. Hollenberg and D. Lellouch, On Searching Neighbors in Computer Simulations of Macromolecular Systems, J. Comput. Chem., 5, 272-279  (1984)

F. Sullivan, R. D. Mountain and J. O'Connell, Molecular Dynamics on Vector Computers, J. Comput. Phys., 61, 138-153 (1985)

J. Boris, A Vectorized "Near Neighbors" Algorithm of Order N Using a Monotonic Logical Grid, J. Comput. Phys., 66, 1-20 (1986)

S. G. Lambrakos and J. P. Boris, Geometric Properties of the Monotonic Lagrangian Grid Algorithm for Near Neighbors Calculations, J. Comput. Phys., 73, 183-202 (1987)

T. A. Andrea, W. C. Swope and H. C. Andersen, The Role of Long Ranged Forces in Determining the Structure and Properties of Liquid Water, J. Chem. Phys., 79, 4576-4584 (1983)

D. N. Theodorou and U. W. Suter, Geometrical Considerations in Model Systems with Periodic Boundary Conditions, J. Chem. Phys., 82, 955-966 (1985)

J. Barnes and P. Hut, A Hierarchical O(NlogN) Force-calculation Algorithm, Nature, 234, 446-449 (1986)

Cutoff and Truncation Methods
-----------------------------

P. J. Steinbach and B. R. Brooks, New Spherical-Cutoff Methods for Long-Range Forces in Macromolecular Simulation, J. Comput. Chem., 15, 667-683 (1993)

R. J. Loncharich and B. R. Brooks, The Effects of Truncating Long-Range Forces on Protein Dynamics, Proteins, 6, 32-45 (1989)

C. L. Brooks III, B. M. Pettitt and M. Karplus, Structural and Energetic Effects of Truncating Long Ranged Interactions in Ionic and Polar Fluids, J. Chem. Phys., 83, 5897-5908 (1985)

Ewald Summation Techniques
--------------------------

A. Y. Toukmaji and J. A. Board, Jr., Ewald Summation Techniques in Perspective: A Survey, Comp. Phys. Commun., 95, 73-92 (1996)

T. Darden, L. Perera, L. Li and L. Pedersen, New Tricks for Modelers from the Crystallography Toolkit: The Particle Mesh Ewald Algorithm and its Use in Nucleic Acid Simulations, Structure, 7, R550-R60 (1999)

T. Darden, D. York and L. G. Pedersen, Particle Mesh Ewald: An Nlog(N) Method for Ewald Sums in Large Systems, J. Chem. Phys., 98, 10089-10092 (1993)

U. Essmann, L. Perera, M. L. Berkowitz, T. Darden, H. Lee and L. G. Pedersen, A Smooth Particle Mesh Ewald Method, J. Chem. Phys., 103, 8577-8593 (1995)

W. Smith, Point Multipoles in the Ewald Summation (Revisited), CCP5 Newsletter, 46, 18-30 (1998)  [available from http://www.dl.ac.uk/CCP/CCP5/newsletter_index.html]

S. E. Feller, R. W. Pastor, A. Rojnuckarin, S. Bogusz and B. R. Brooks, Effect of Electrostatic Force Truncation on Interfacial and Transport Properties of Water, J. Phys. Chem., 100, 17011-17020 (1996)

W. Weber, P. H. H¸nenberger and J. A. McCammon, Molecular Dynamics Simulations of a Polyalanine Octapeptide under Ewald Boundary Conditions: Influence of Artificial Periodicity on Peptide Conformation, J. Phys. Chem. B, 104, 3668-3675 (2000)

Conjugated and Aromatic Systems
-------------------------------

N. L. Allinger, F. Li, L. Yan and J. C. Tai, Molecular Mechanics (MM3) Calculations on Conjugated Hydrocarbons, J. Comput. Chem., 11, 868-895 (1990)

J. T. Sprague, J. C. Tai, Y. Yuh and N. L. Allinger, The MMP2 Calculational Method, J. Comput. Chem., 8, 581-603 (1987)

J. Kao, A Molecular Orbital Based Molecular Mechanics Approach to Study Conjugated Hydrocarbons, J. Am. Chem. Soc., 109, 3818-3829 (1987)

J. Kao and N. L. Allinger, Conformational Analysis: Heats of Formation of Conjugated Hydrocarbons by the Force Field Method, J. Am. Chem. Soc., 99, 975-986 (1977)

D. H. Lo and M. A. Whitehead, Accurate Heats of Atomization and Accurate Bond Lengths: Benzenoid Hydrocarbons, Can. J. Chem., 46, 2027-2040 (1968)

G. D. Zeiss and M. A. Whitehead, Hetero-atomic Molecules: Semi-empirical Molecular Orbital Calculations and Prediction of Physical Properties, J. Chem. Soc. A, 1727-1738 (1971)

Free Energy Simulation Methods
------------------------------

P. Kollman, Free Energy Calculations: Applications to Chemical and Biochemical Phenomena, Chem. Rev., 93, 2395-2417 (1993)

B. L. Tembe and J. A. McCammon, Ligand-Receptor Interactions, Computers & Chemistry, 8, 281-283 (1984)

W. L. Jorgensen and C. Ravimohan, Monte Carlo Simulation of Differences in Free Energy of Hydration, J. Chem. Phys., 83, 3050-3054 (1985)

W. L. Jorgensen, J. K. Buckner, S. Boudon and J. Tirado-Rives, Efficient Computation of Absolute Free Energies of Binding by Computer Simulations:  Application to the Methane Dimer in Water, J. Chem. Phys., 89, 3742-3746 (1988)

S. H. Fleischman and C. L. Brooks III, Thermodynamics of Aqueous Solvation:  Solution Properties of Alcohols and Alkanes, J. Chem. Phys., 87, 3029-3037 (1987)

U. C. Singh, F. K. Brown, P. A. Bash and P. A. Kollman, An Approach to the Application of Free Energy Perturbation Methods Using Molecular Dynamics, J. Am. Chem. Soc., 109, 1607-1614 (1987)

D. A. Pearlman and P. A. Kollman, A New Method for Carrying out Free Energy Perturbation Calculations: Dynamically Modified Windows, J. Chem. Phys., 90, 2460-2470 (1989)

T. P. Straatsma, H. J. C. Berendsen and J. P. M. Postma, Free Energy of Hydrophobic Hydration:  A Molecular Dynamics Study of Noble Gases in Water, J. Chem. Phys., 85, 6720-6727 (1986)

T. P. Straatsma and H. J. C. Berendsen, Free Energy of Ionic Hydration:  Analysis of a Thermodynamic Integration Technique to Evaluate Free Energy Differences by Molecular Dynamics Simulations, J. Chem. Phys., 89, 5876-5886 (1988)

M. Mezei, The Finite Difference Thermodynamic Integration, Tested on Calculating the Hydration Free Energy Difference between Acetone and Dimethylamine in Water, J. Chem. Phys., 86, 7084-7088 (1987)

A. E. Mark and W. F. van Gunsteren, Decomposition of the Free Energy of a System in Terms of Specific Interactions, J. Mol. Biol., 240, 167-176 (1994)

S. Boresch and M. Karplus, The Meaning of Copmponent Analysis: Decomposition of the Free Energy in Terms of Specific Interactions, J. Mol. Biol., 254, 801-807 (1995)

Methods for Parameter Determination
-----------------------------------

N. L. Allinger, X. Zhou and J. Bergsma, Molecular Mechanics Parameters, J. Mol. Struct. (THEOCHEM), 312, 69-83 (1994)

A. J. Pertsin and A. I. Kitaigorodsky, The Atom-Atom Potential Method: Application to Organic Molecular Solids, Springer-Verlag, Berlin, 1987

D. E. Williams, Transferable Empirical Nonbonded Potential Functions, in Crystal Cohesion and Conformational Energies, Ed. by R. M. Metzger, Springer-Verlag, Berlin, 1981

A. T. Hagler and S. Lifson, A Procedure for Obtaining Energy Parameters from Crystal Packing, Acta Cryst., B30, 1336-1341 (1974)

A. T. Hagler, S. Lifson and P. Dauber, Consistent Force Field Studies of Intermolecular Forces in Hydrogen-Bonded Crystals:  A Benchmark for the Objective Comparison of Alternative Force Fields, J. Am. Chem. Soc., 101, 5122-5130 (1979)

W. L. Jorgensen, J. D. Madura and C. J. Swenson, Optimized Intermolecular Potential Functions for Liquid Hydrocarbons, J. Am. Chem. Soc., 106, 6638-6646 (1984)

W. L. Jorgensen and C. J. Swenson, Optimized Intermolecular Potential Functions for Amides and Peptides: Structure and Properties of Liquid Amides, J. Am. Chem. Soc., 107, 569-578 (1985)

J. R. Maple, U. Dinur and A. T. Hagler, Derivation of Force Fields for Molecular Mechanics and Dynamics from ab Initio Surfaces, Proc. Nat. Acad. Sci. USA, 85, 5350-5354 (1988)

U. Dinur and A. T. Hagler, Direct Evaluation of Nonbonding Interactions from ab Initio Calculations, J. Am. Chem. Soc., 111, 5149-5151 (1989)

Electrostatic Interactions
--------------------------

S. L. Price, Towards More Accurate Model Intermolecular Potentials for Organic Molecules, Rev. Comput. Chem., 14, 225-289 (2000)

C. H. Faerman and S. L. Price, A Transferable Distributed Multipole Model for the Electrostatic Interactions of Peptides and Amides, J. Am. Chem. Soc., 112, 4915-4926 (1990)

C. E. Dykstra, Electrostatic Interaction Potentials in Molecular Force Fields, Chem. Rev., 93, 2339-2353 (1993)

M. J. Dudek and J. W. Ponder, Accurate Modeling of the Intramolecular Electrostatic Energy of Proteins, J. Comput. Chem., 16, 791-816 (1995)

U. Koch and E. Egert, An Improved Description of the Molecular Charge Density in Force Fields with Atomic Multipole Moments, J. Comput. Chem., 16, 937-944 (1995)

D. E. Williams, Representation of the Molecular Electrostatic Potential by Atomic Multipole and Bond Dipole Models, J. Comput. Chem., 9, 745-763 (1988)

F. Colonna, E. Evleth and J. G. Angyan, Critical Analysis of Electric Field Modeling: Formamide, J. Comput. Chem., 13, 1234-1245 (1992)

Polarization Effects
--------------------

S. Kuwajima and A. Warshel, Incorporating Electric Polarizabilities in Water-Water Interaction Potentials, J. Phys. Chem., 94, 460-466 (1990)

J. W. Caldwell and P. A. Kollman, Structure and Properties of Neat Liquids Using Nonadditive Molecular Dynamics: Water, Methanol, and N-Methylacetamide, J. Phys. Chem., 99, 6208-6219 (1995)

D. N. Bernardo, Y. Ding, K. Kroegh-Jespersen and R. M. Levy, An Anisotropic Polarizable Water Model: Incorporation of All-Atom Polarizabilities into Molecular Mechanics Force Fields, J. Phys. Chem., 98, 4180-4187 (1994)

P. T. van Duijnen and M. Swart, Molecular and Atomic Polarizabilities: Thole's Model Revisited, J. Phys. Chem. A, 102, 2399-2407 (1998)

K. J. Miller, Calculation of the Molecular Polarizability Tensor, J. Am. Chem. Soc., 112, 8543-8551 (1990)

J. Applequist, J. R. Carl and K.-K. Fung, An Atom Dipole Interaction Model for Molecular Polarizability. Application to Polyatomic Molecules and Determination of Atom Polarizabilities, J. Am. Chem. Soc., 94, 2952-2960 (1972)

J. Applequist, Atom Charge Transfer in Molecular Polarizabilities. Application of the Olson-Sundberg Model to Aliphatic and Aromatic Hydrocarbons, J. Phys. Chem., 97, 6016-6023 (1993)

A. J. Stone, Distributed Polarizabilities, Mol. Phys., 56, 1065-1082 (1985)

J. M. Stout and C. E. Dykstra, A Distributed Model of the Electrical Response of Organic Molecules, J. Phys. Chem. A, 102, 1576-1582 (1998)

Macroscopic Treatment of Solvent
--------------------------------

C. J. Cramer and D. G. Truhlar, Continuum Solvation Models: Classical and Quantum Mechanical Implementations, Rev. Comput. Chem., 6, 1-72 (1995)

B. Roux and T. Simonson, Implicit Solvation Models, Biophys. Chem., 78, 1-20 (1999)

M. K. Gilson, Introduction to Continuum Electrostatics with Molecular Applications, available from http://gilsonlab.umbi.umd.edu

Surface Area-Based Solvation Models
-----------------------------------

D. Eisenberg and A. D. McLachlan, Solvation Energy in Protein Folding and Binding, Nature, 319, 199-203 (1986)

L. Wesson and D. Eisenberg, Atomic Solvation Parameters Applied to Molecular Dynamics of Proteins in Solution, Prot. Sci., 1, 227-235 (1992)

T. Ooi, M. Oobatake, G. Nemethy and H. A. Scheraga, Accessible Surface Areas as a Measure of the Thermodynamic Parameters of Hydration of Peptides, Proc. Natl. Acad. Sci. USA, 84, 3086-3090 (1987)

J. D. Augspurger and H. A. Scheraga, An Efficient, Differentiable Hydration Potential for Peptides and Proteins, J. Comput. Chem., 17, 1549-1558 (1996)

Generalized Born Solvation Models
---------------------------------

W. C. Still, A. Tempczyk, R. C. Hawley and T. Hendrickson, A Semiempirical Treatment of Solvation for Molecular Mechanics and Dynamics, J. Am. Chem. Soc., 112, 6127-6129 (1990)

D. Qiu, P. S. Shenkin, F. P. Hollinger and W. C. Still, The GB/SA Continuum Model for Solvation. A Fast Analytical Method for the Calculation of Approximate Born Radii, J. Phys. Chem. A, 101, 3005-3014 (1997)

G. D. Hawkins, C. J. Cramer and D. G. Truhlar, Pairwise Solute Descreening of Solute Charges from a Dielectric Medium, Chem. Phys. Lett., 246, 122-129 (1995)

G. D. Hawkins, C. J. Cramer and D. G. Truhlar, Parametrized Models of Aqueous Free Energies of Solvation Based on Pairwise Descreening of Solute Atomic Charges from a Dielectric Medium, J. Phys. Chem., 100, 19824-19839 (1996)

A. Onufriev, D. Bashford and D. A. Case, Modification of the Generalized Born Model Suitable for Macromolecules, J. Phys. Chem. B, 104, 3712-3720 (2000)

M. Schaefer and M. Karplus, A Comprehensive Analytical Treatment of Continuum Electrostatics, J. Phys. Chem., 100, 1578-1599 (1996)

M. Schaefer, C. Bartels and M. Karplus, Solution Conformations and Thermodynamics of Structured Peptides: Molecular Dynamics Simulation with an Implicit Solvation Model, J. Mol. Biol., 284, 835-848 (1998)

Superposition of Coordinate Sets
--------------------------------

S. J. Kearsley, An Algorithm for the Simultaneous Superposition of a Structural Series, J. Comput. Chem., 11, 1187-1192 (1990)

R. Diamond, A Note on the Rotational Superposition Problem, Acta Cryst., A44, 211-216 (1988)

A. D. McLachlan, Rapid Comparison of Protein Structures, Acta Cryst., A38, 871-873 (1982)

S. C. Nyburg, Some Uses of a Best Molecular Fit Routine, Acta Cryst., B30, 251-253 (1974)

Location of Transition States
-----------------------------

R. Czerminski and R. Elber, Reaction Path Study of Conformational Transitions and Helix Formation in a Tetrapeptide, Proc. Nat. Acad. Sci. USA, 86, 6963 (1989)

R. S. Berry, H. L. Davis and T. L. Beck, Finding Saddles on Multidimensional Potential Surfaces, Chem. Phys. Lett., 147, 13 (1988)

K. Muller, Reaction Paths on Multidimensional Energy Hypersurfaces, Ang. Chem. Int. Ed. Engl., 19, 1-13 (1980)

S. Bell and J. S. Crighton, Locating Transition States, J. Chem. Phys., 80, 2464-2475 (1984)

S. Fischer and M. Karplus, Conjugate Peak Refinement: An Algorithm for Finding Reaction Paths and Accurate Transition States in Systems with Many Degrees of Freedom, Chem. Phys. Lett., 194, 252-261 (1992)

J. E. Sinclair and R. Fletcher, A New Method of Saddle-Point Location for the Calculation of Defect Migration Energies, J. Phys. C, 7, 864-870 (1974)

R. Elber and M. Karplus, A Method for Determining Reaction Paths in Large Molecules:  Application to Myoglobin, Chem. Phys. Lett., 139, 375-380 (1987)

D. T. Nguyen and D. A. Case, On Finding Stationary States on Large-Molecule Potential Energy Surfaces, J. Phys. Chem., 89, 4020-4026 (1985)

T. A. Halgren and W. N. Lipscomb, The Synchronous-Transit Method for Determining Reaction Pathways and Locating Molecular Transition States, Chem. Phys. Lett., 49, 225-232 (1977)

G. T. Barkema and N. Mousseau, Event-Based Relaxation of Continuous Disordered Systems, Phys. Rev. Lett., 77, 4358-4361 (1996)
