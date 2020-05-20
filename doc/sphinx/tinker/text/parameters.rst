Force Field Parameter Sets
==========================

The Tinker package is distributed with several force field parameter sets, implementing a selection of widely used literature force fields as well as the Tinker force field currently under construction in the Ponder lab. We try to exactly reproduce the intent of the original authors of our distributed, third-party force fields. In all cases the parameter sets have been validated against literature reports, results provided by the original developers, or calculations made with the authentic programs. With the few exceptions noted below, Tinker calculations can be treated as authentic results from the genuine force fields. A brief description of each parameter set, including some still in preparation and not distributed with the current version, is provided below with lead literature references for the force field:

**AMOEBA.PRM**

Parameters for the AMOEBA polarizable atomic multipole force field. As of the current Tinker release, we have completed parametrization for a number of ions and small organic molecules. For further information, or if you are interested in developing or testing parameters for other small molecules, please contact the Ponder lab.

P. Ren and J. W. Ponder, A Consistent Treatment of Inter- and Intramolecular Polarization in Molecular Mechanics Calculations, J. Comput. Chem., 23, 1497-1506 (2002)

P. Ren and J. W. Ponder, Polarizable Atomic Multipole Water Model for Molecular Mechanics Simulation, J. Phys. Chem. B, 107, 5933-5947 (2003)

P. Ren and J. W. Ponder, Ion Solvation Thermodynamics from Simulation with a Polarizable Force Field, A. Grossfield, J. Am. Chem. Soc., 125, 15671-15682 (2003)

**AMOEBAPRO.PRM**

Preliminary protein parameters for the AMOEBA polarizable atomic multipole force field. While the distributed parameters are still subject to minor alteration as we continue validation, they are now stable enough for other groups to begin using them. For further information, or if you are interested in testing the protein parameter set, please contact the Ponder lab.

J. W. Ponder and D. A. Case, Force Fields for Protein Simulation, Adv. Prot. Chem., 66, 27-85 (2003)

P. Ren and J. W. Ponder, Polarizable Atomic Multipole-based Potential for Proteins: Model and Parameterization, in preparation

**AMBER94.PRM**

AMBER ff94 parameters for proteins and nucleic acids. Note that with their "Cornell" force field, the Kollman group has devised separate, fully independent partial charge values for each of the N- and C-terminal amino acid residues. At present, the terminal residue charges for Tinker's version maintain the correct formal charge, but redistributed somewhat at the alpha carbon atoms from the original Kollman group values. The total magnitude of the redistribution is less than 0.01 electrons in most cases.

W. D. Cornell, P. Cieplak, C. I. Bayly, I. R. Gould, K. M. Merz, Jr., D. M. Ferguson, D. C. Spellmeyer, T. Fox, J. W. Caldwell and P. A. Kollman, A Second Generation Force Field for the Simulation of Proteins, Nucleic Acids, and Organic Molecules, J. Am. Chem. Soc., 117, 5179-5197 (1995)  [ff94]

G. Moyna, H. J. Williams, R. J. Nachman and A. I. Scott, Conformation in Solution and Dynamics of a Structurally Constrained Linear Insect Kinin Pentapeptide Analogue, Biopolymers, 49, 403-413 (1999)  [AIB charges]

W. S. Ross and C. C. Hardin, Ion-Induced Stabilization of the G-DNA Quadruplex: Free Energy Perturbation Studies, J. Am. Chem. Soc., 116, 4363-4366 (1994)   [alkali metal ions]

J. Aqvist, Ion-Water Interaction Potentials Derived from Free Energy Perturbation Simulations, J. Phys. Chem., 94, 8021-8024, 1990  [alkaline earth Ions, radii adapted for Amber combining rule]

Current force field parameter values and suggested procedures for development of parameters for additional molecules are available from the Amber web site in the Case lab at Scripps, http://amber.scripps.edu/

**AMBER96.PRM**

AMBER ff96 parameters for proteins and nucleic acids. The only change from the ff94 parameter set is in the torsional parameters for the protein phi/psi angles. These values were altered to give better agreement with  changes of ff96 with LMP2 QM results from the Friesner lab on alanine dipeptide and tetrapeptide.

P. Kollman, R. Dixon, W. Cornell, T. Fox, C. Chipot and A. Pohorille, The Development/ Application of a 'Minimalist' Organic/Biochemical Molecular Mechanic Force Field using a Combination of ab Initio Calculations and Experimental Data, in Computer Simulation of Biomolecular Systems, W. F. van Gunsteren, P. K. Weiner, A. J. Wilkinson, eds., Volume 3, 83-96 (1997)  [ff96]

Current force field parameter values and suggested procedures for development of parameters for additional molecules are available from the Amber web site in the Case lab at Scripps, http://amber.scripps.edu/

**AMBER98.PRM**

AMBER ff98 parameters for proteins and nucleic acids. The only change from the ff94 parameter set is in the glycosidic torsional parameters that control sugar pucker.

T. E. Cheatham III, P. Cieplak and P. A. Kollman, A Modified Version of the Cornell et al. Force Field with Improved Sugar Pucker Phases and Helical Repeat, J. Biomol. Struct. Dyn., 16, 845-862 (1999)

Current force field parameter values and suggested procedures for development of parameters for additional molecules are available from the Amber web site in the Case lab at Scripps, http://amber.scripps.edu/

**AMBER99.PRM**

AMBER ff99 parameters for proteins and nucleic acids. The original partial charges from the ff94 parameter set are retained, but many of the bond, angle and torsional parameters have been revised to provide better general agreement with experiment.

J. Wang, P. Cieplak and P. A. Kollman, How Well Does a Restrained Electrostatic Potential (RESP) Model Perform in Calcluating Conformational Energies of Organic and Biological Molecules?, J. Comput. Chem., 21, 1049-1074 (2000)

Current force field parameter values and suggested procedures for development of parameters for additional molecules are available from the Amber web site in the Case lab at Scripps, http://amber.scripps.edu/

**CHARMM19.PRM**

CHARMM19 united-atom parameters for proteins. The nucleic acid parameter are not yet implemented. There are some differences between authentic CHARMM19 and the Tinker version due to replacement of CHARMM impropers by torsions for cases that involve atoms not bonded to the trigonal atom and Tinker's use of all possible torsions across a bond instead of a single torsion per bond.

E. Neria, S. Fischer and M. Karplus, Simulation of Activation Free Energies in Molecular Systems, J. Chem. Phys., 105, 1902-1921 (1996)

L. Nilsson and M. Karplus, Empirical Energy Functions for Energy Minimizations and Dynamics of Nucleic Acids, J. Comput. Chem., 7, 591-616 (1986)

W. E. Reiher III, Theoretical Studies of Hydrogen Bonding, Ph.D. Thesis, Department of Chemistry, Harvard University, Cambridge, MA, 1985

**CHARMM22.PRM**

CHARMM22 all-atom parameters for proteins and lipids. Most of the nucleic acid and small model compound parameters are not yet implemented. We plan to provide these additional parameters in due course.

N. Foloppe and A. D. MacKerell, Jr., All-Atom Empirical Force Field for Nucleic Acids: 1) Parameter Optimization Based on Small Molecule and Condensed Phase Macromolecular Target Data, J. Comput. Chem., 21, 86-104 (2000)  [CHARMM27]

N. Banavali and A. D. MacKerell, Jr., All-Atom Empirical Force Field for Nucleic Acids: 2) Application to Molecular Dynamics Simulations of DNA and RNA in Solution, J. Comput. Chem., 21, 105-120 (2000)

A. D. MacKerrell, Jr., et al., All-Atom Empirical Potential for Molecular Modeling and Dynamics Studies of Proteins, J. Phys. Chem. B, 102, 3586-3616 (1998)  [CHARMM22]

A. D. MacKerell, Jr., J. Wiorkeiwicz-Kuczera and M. Karplus, An All-Atom Empirical Energy Function for the Simulation of Nucleic Acids, J. Am. Chem. Soc., 117, 11946-11975 (1995)

S. E. Feller, D. Yin, R. W. Pastor and A. D. MacKerell, Jr., Molecular Dynamics Simulation of Unsaturated Lipids at Low Hydration: Parametrization and Comparison with Diffraction Studies, Biophysical Journal, 73, 2269-2279 (1997)  [alkenes]

R. H. Stote and M. Karplus, Zinc Binding in Proteins and Solution - A Simple but Accurate Nonbonded Representation, Proteins, 23, 12-31 (1995)  [zinc ion]

Current and legacy parameter values are available from the CHARMM force field web site on Alex MacKerell's  Research Interests page at the University of Maryland School of Pharmacy, https://rxsecure.umaryland.edu/research/amackere/research.html/

**DUDEK.PRM**

Protein-only parameters for the early 1990's Tinker force field with multipole values of Dudek and Ponder. The current file contains only the multipole values from the 1995 paper by Dudek and Ponder. This set is now superceeded by the more recent Tinker force field developed by Pengyu Ren (see WATER.PRM, below).

M. J. Dudek and J. W. Ponder, Accurate Electrostatic Modelling of the Intramolecular Energy of Proteins, J. Comput. Chem., 16, 791-816 (1995)

**ENCAD.PRM**

ENCAD parameters for proteins and nucleic acids.  (in preparation)

M. Levitt, M. Hirshberg, R. Sharon and V. Daggett, Potential Energy Function and Parameters for Simulations of the Molecular Dynamics of Protein and Nucleic Acids in Solution, Comp. Phys. Commun., 91, 215-231 (1995)

M. Levitt, M. Hirshberg, R. Sharon, K. E. Laidig and V. Daggett, Calibration and Testing of a Water Model for Simulation of the Molecular Dynamics of Protein and Nucleic Acids in Solution, J. Phys. Chem. B, 101, 5051-5061 (1997)  [F3C water]

**HOCH.PRM**

Simple NMR-NOE force field of Hoch and Stern.

J. C. Hoch and A. S. Stern, A Method for Determining Overall Protein Fold from NMR Distance Restraints, J. Biomol. NMR, 2, 535-543 (1992)

**MM2.PRM**

Full MM2(1991) parameters including ?-systems. The anomeric and electronegativity correction terms included in some later versions of MM2 are not implemented.

N. L. Allinger, Conformational Analysis. 130. MM2. A Hydrocarbon Force Field Utilizing V1 and V2 Torsional Terms, J. Am. Chem. Soc., 99, 8127-8134 (1977)

J. T. Sprague, J. C. Tai, Y. Yuh and N. L. Allinger, The MMP2 Calculational Method, J. Comput. Chem., 8, 581-603 (1987)

J. C. Tai and N. L. Allinger, Molecular Mechanics Calculations on Conjugated Nitrogen-Containing Heterocycles, J. Am. Chem. Soc., 110, 2050-2055 (1988)

J. C. Tai, J.-H. Lii and N. L. Allinger, A Molecular Mechanics (MM2) Study of Furan, Thiophene, and Related Compounds, J. Comput. Chem., 10, 635-647 (1989)

N. L. Allinger, R. A. Kok and M. R. Imam, Hydrogen Bonding in MM2, J. Comput. Chem., 9, 591-595 (1988)

L. Norskov-Lauritsen and N. L. Allinger, A Molecular Mechanics Treatment of the Anomeric Effect, J. Comput. Chem., 5, 326-335 (1984)

All parameters distributed with Tinker are from the "MM2 (1991) Parameter Set", as provided by N. L. Allinger, University of Georgia

**MM3.PRM**

Full MM3(2000) parameters including pi-systems. The directional hydrogen bonding term and electronegativity bond length corrections are implemented, but the anomeric and Bohlmann correction terms are not implemented.

N. L. Allinger, Y. H. Yuh and J.-H. Lii, Molecular Mechanics. The MM3 Force Field for Hydrocarbons. 1, J. Am. Chem. Soc., 111, 8551-8566 (1989)

J.-H. Lii and N. L. Allinger, Molecular Mechanics. The MM3 Force Field for Hydrocarbons. 2. Vibrational Frequencies and Thermodynamics, J. Am. Chem. Soc., 111, 8566-8575 (1989)

J.-H. Lii and N. L. Allinger, Molecular Mechanics. The MM3 Force Field for Hydrocarbons. 3. The van der Waals' Potentials and Crystal Data for Aliphatic and Aromatic Hydrocarbons, J. Am. Chem. Soc., 111, 8576-8582 (1989)

N. L. Allinger, H. J. Geise, W. Pyckhout, L. A. Paquette and J. C. Gallucci, Structures of Norbornane and Dodecahedrane by Molecular Mechanics Calculations (MM3), X-ray Crystallography, and Electron Diffraction, J. Am. Chem. Soc., 111, 1106-1114 (1989)  [stretch-torsion cross term]

N. L. Allinger, F. Li and L. Yan, Molecular Mechanics. The MM3 Force Field for Alkenes, J. Comput. Chem., 11, 848-867 (1990)

N. L. Allinger, F. Li, L. Yan and J. C. Tai, Molecular Mechanics (MM3) Calculations on Conjugated Hydrocarbons, J. Comput. Chem., 11, 868-895 (1990)

J.-H. Lii and N. L. Allinger, Directional Hydrogen Bonding in the MM3 Force Field. I, J. Phys. Org. Chem., 7, 591-609 (1994)

J.-H. Lii and N. L. Allinger, Directional Hydrogen Bonding in the MM3 Force Field. II, J. Comput. Chem., 19, 1001-1016 (1998)

All parameters distributed with Tinker are from the "MM3 (2000) Parameter Set", as provided by N. L. Allinger, University of Georgia, August 2000

**MM3PRO.PRM**

Protein-only version of the MM3 parameters.

J.-H. Lii and N. L. Allinger, The MM3 Force Field for Amides, Polypeptides and Proteins, J. Comput. Chem., 12, 186-199 (1991)

**OPLSUA.PRM**

Complete OPLS-UA with united-atom parameters for proteins and many classes of organic molecules. Explicit hydrogens on polar atoms and aromatic carbons.

W. L. Jorgensen and J. Tirado-Rives, The OPLS Potential Functions for Proteins. Energy Minimizations for Crystals of Cyclic Peptides and Crambin, J. Am. Chem. Soc., 110, 1657-1666 (1988)  [peptide and proteins]

W. L. Jorgensen and D. L. Severance, Aromatic-Aromatic Interactions: Free Energy Profiles for the Benzene Dimer in Water, Chloroform, and Liquid Benzene, J. Am. Chem. Soc., 112, 4768-4774 (1990)  [aromatic hydrogens]

S. J. Weiner, P. A. Kollman, D. A. Case, U. C. Singh, C. Ghio, G. Alagona, S. Profeta, Jr. and P. Weiner, A New Force Field for Molecular Mechanical Simulation of Nucleic Acids and Proteins, J. Am. Chem. Soc., 106, 765-784 (1984)  [united-atom "AMBER/OPLS" local geometry]

S. J. Weiner, P. A. Kollman, D. T. Nguyen and D. A. Case, An All Atom Force Field for Simulations of Proteins and Nucleic Acids, J. Comput. Chem., 7, 230-252 (1986)  [all-atom "AMBER/OPLS" local geometry]

L. X. Dang and B. M. Pettitt, Simple Intramolecular Model Potentials for Water, J. Phys. Chem., 91, 3349-3354 (1987)  [flexible TIP3P and SPC water]

W. L. Jorgensen, J. D. Madura and C. J. Swenson, Optimized Intermolecular Potential Functions for Liquid Hydrocarbons, J. Am. Chem. Soc., 106, 6638-6646 (1984)  [hydrocarbons]

W. L. Jorgensen, E. R. Laird, T. B. Nguyen and J. Tirado-Rives, Monte Carlo Simulations of Pure Liquid Substituted Benzenes with OPLS Potential Functions, J. Comput. Chem., 14, 206-215 (1993)  [substituted benzenes]

E. M. Duffy, P. J. Kowalczyk and W. L. Jorgensen, Do Denaturants Interact with Aromatic Hydrocarbons in Water?, J. Am. Chem. Soc., 115, 9271-9275 (1993)  [benzene, naphthalene, urea, guanidinium, tetramethyl ammonium]

W. L. Jorgensen and C. J. Swenson, Optimized Intermolecular Potential Functions for Amides and Peptides. Structure and Properties of Liquid Amides, J. Am. Chem. Soc., 106, 765-784 (1984)  [amides]

W. L. Jorgensen, J. M. Briggs and M. L. Contreras, Relative Partition Coefficients for Organic Solutes form Fluid Simulations, J. Phys. Chem., 94, 1683-1686 (1990)  [chloroform, pyridine, pyrazine, pyrimidine]

J. M. Briggs, T. B. Nguyen and W. L. Jorgensen, Monte Carlo Simulations of Liquid Acetic Acid and Methyl Acetate with the OPLS Potential Functions, J. Phys. Chem., 95, 3315-3322 (1991)  [acetic acid, methyl acetate]

H. Liu, F. Muller-Plathe and W. F. van Gunsteren, A Force Field for Liquid Dimethyl Sulfoxide and Physical Properties of Liquid Dimethyl Sulfoxide Calculated Using Molecular Dynamics Simulation, J. Am. Chem. Soc., 117, 4363-4366 (1995)  [dimethyl sulfoxide]

J. Gao, X. Xia and T. F. George, Importance of Bimolecular Interactions in Developing Empirical Potential Functions for Liquid Ammonia, J. Phys. Chem., 97, 9241-9246 (1993)  [ammonia]

J. Aqvist, Ion-Water Interaction Potentials Derived from Free Energy Perturbation Simulations, J. Phys. Chem., 94, 8021-8024 (1990)  [metal ions]

W. S. Ross and C. C. Hardin, Ion-Induced Stabilization of the G-DNA Quadruplex: Free Energy Perturbation Studies, J. Am. Chem. Soc., 116, 4363-4366 (1994)  [alkali metal ions]

J. Chandrasekhar, D. C. Spellmeyer and W. L. Jorgensen, Energy Component Analysis for Dilute Aqueous Solutions of Li+, Na+, F-, and Cl- Ions, J. Am. Chem. Soc., 106, 903-910 (1984)  [halide ions]

Most parameters distributed with Tinker are from "OPLS and OPLS-AA Parameters for Organic Molecules, Ions, and Nucleic Acids" as provided by W. L. Jorgensen, Yale University, October 1997

**OPLSAA.PRM**

OPLS-AA force field with all-atom parameters for proteins and many general classes of organic molecules.

W. L. Jorgensen, D. S. Maxwell and J. Tirado-Rives, Development and Testing of the OPLS All-Atom Force Field on Conformational Energetics and Properties of Organic Liquids, J. Am. Chem. Soc., 117, 11225-11236 (1996)

D. S. Maxwell, J. Tirado-Rives and W. L. Jorgensen, A Comprehensive Study of the Rotational Energy Profiles of Organic Systems by Ab Initio MO Theory, Forming a Basis for Peptide Torsional Parameters, J. Comput. Chem., 16, 984-1010 (1995)

W. L. Jorgensen and N. A. McDonald, Development of an All-Atom Force Field for Heterocycles. Properties of Liquid Pyridine and Diazenes, THEOCHEM-J. Mol. Struct., 424, 145-155 (1998)

N. A. McDonald and W. L. Jorgensen, Development of an All-Atom Force Field for Heterocycles. Properties of Liquid Pyrrole, Furan, Diazoles, and Oxazoles, J. Phys. Chem. B, 102, 8049-8059 (1998)

R. C. Rizzo and W. L. Jorgensen, OPLS All-Atom Model for Amines: Resolution of the Amine Hydration Problem, J. Am. Chem. Soc., 121, 4827-4836 (1999)

M. L. P. Price, D. Ostrovsky and W. L. Jorgensen, Gas-Phase and Liquid-State Properties of Esters, Nitriles, and Nitro Compounds with the OPLS-AA Force Field, J. Comput. Chem., 22, 1340-1352 (2001)

All parameters distributed with Tinker are from "OPLS and OPLS-AA Parameters for Organic Molecules, Ions, and Nucleic Acids" as provided by W. L. Jorgensen, Yale University, October 1997

**OPLSAAL.PRM**

An improved OPLS-AA parameter set for proteins in which the only change is a reworking of many of the backbone and sidechain torsional parameters to give better agreement with LMP2 QM calculations. This parameter set is also known as OPLS(2000).

G. A. Kaminsky, R. A. Friesner, J. Tirado-Rives and W. L. Jorgensen, Evaluation and Reparametrization of the OPLS-AA Force Field for Proteins via Comparison with Accurate Quantum Chemical Calculations on Peptides, J. Phys. Chem. B, 105, 6474-6487 (2001)

**SMOOTH.PRM**

Version of OPLS-UA for use with potential smoothing. Largely adapted largely from standard OPLS-UA parameters with modifications to the vdw and improper torsion terms.

R. V. Pappu, R. K. Hart and J. W. Ponder, Analysis and Application of Potential Energy Smoothing and Search Methods for Global Optimization, J. Phys, Chem. B, 102, 9725-9742 (1998)  [smoothing modifications]

**SMOOTHAA.PRM**

Version of OPLS-AA for use with potential smoothing. Largely adapted largely from standard OPLS-AA parameters with modifications to the vdw and improper torsion terms.

R. V. Pappu, R. K. Hart and J. W. Ponder, Analysis and Application of Potential Energy Smoothing and Search Methods for Global Optimization, J. Phys, Chem. B, 102, 9725-9742 (1998)  [smoothing modifications]

**WATER.PRM**

The AMOEBA water parameters for a polarizable atomic multipole electrostatics model. This model is equal or better to the best available water models for many bulk and cluster properties.

P. Ren and J. W. Ponder, A Polarizable Atomic Multipole Water Model for Molecular Mechanics Simulation, J. Phys. Chem. B, 107, 5933-5947 (2003) 

P. Ren and J. W. Ponder, Ion Solvation Thermodynamics from Simulation with a Polarizable Force Field, A. Grossfield, J. Am. Chem. Soc., 125, 15671-15682 (2003)

P. Ren and J. W. Ponder, Temperature and Pressure Dependence of the AMOEBA Water Model, J. Phys. Chem. B, 108, 13427-13437 (2004)

An earlier version the AMOEBA water model is described in: Yong Kong, Multipole Electrostatic Methods for Protein Modeling with Reaction Field Treatment, Biochemistry & Molecular Biophysics, Washington University, St. Louis, August, 1997 [available from http://dasher.wustl.edu/ponder/]
