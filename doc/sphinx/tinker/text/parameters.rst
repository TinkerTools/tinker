Force Field Parameter Sets
==========================

The Tinker package is distributed with several force field parameter sets, implementing a selection of widely used literature force fields as well as the Tinker force field currently under construction in the Ponder lab. We try to exactly reproduce the intent of the original authors of our distributed, third-party force fields. In all cases the parameter sets have been validated against literature reports, results provided by the original developers, or calculations made with the authentic programs. With the few exceptions noted below, Tinker calculations can be treated as authentic results from the genuine force fields. A brief description of each parameter set, including some still in preparation and not distributed with the current version, is provided below with lead literature references for the force field:

**AMBER94.PRM**

Amber ff94 parameters for proteins and nucleic acids. Note that with their "Cornell" force field, the Kollman group has devised separate, fully independent partial charge values for each of the N- and C-terminal amino acid residues. At present, the terminal residue charges for Tinker's version maintain the correct formal charge, but redistributed somewhat at the alpha carbon atoms from the original Kollman group values. The total magnitude of the redistribution is less than 0.01 electrons in most cases.

W. D. Cornell, P. Cieplak, C. I. Bayly, I. R. Gould, K. M. Merz, Jr., D. M. Ferguson, D. C. Spellmeyer, T. Fox, J. W. Caldwell and P. A. Kollman, A Second Generation Force Field for the Simulation of Proteins, Nucleic Acids, and Organic Molecules, J. Am. Chem. Soc., 117, 5179-5197 (1995)  [ff94]

G. Moyna, H. J. Williams, R. J. Nachman and A. I. Scott, Conformation in Solution and Dynamics of a Structurally Constrained Linear Insect Kinin Pentapeptide Analogue, Biopolymers, 49, 403-413 (1999)  [AIB charges]

W. S. Ross and C. C. Hardin, Ion-Induced Stabilization of the G-DNA Quadruplex: Free Energy Perturbation Studies, J. Am. Chem. Soc., 116, 4363-4366 (1994)   [alkali metal ions]

J. Aqvist, Ion-Water Interaction Potentials Derived from Free Energy Perturbation Simulations, J. Phys. Chem., 94, 8021-8024, 1990  [alkaline earth Ions, radii adapted for Amber combining rule]

Current force field parameter values and suggested procedures for development of parameters for additional molecules are available from the Amber web site in the Case lab at Scripps, http://amber.scripps.edu/

**AMBER96.PRM**

Amber ff96 parameters for proteins and nucleic acids. The only change from the ff94 parameter set is in the torsional parameters for the protein phi/psi angles. These values were altered to give better agreement with  changes of ff96 with LMP2 QM results from the Friesner lab on alanine dipeptide and tetrapeptide.

P. Kollman, R. Dixon, W. Cornell, T. Fox, C. Chipot and A. Pohorille, The Development/ Application of a 'Minimalist' Organic/Biochemical Molecular Mechanic Force Field using a Combination of ab Initio Calculations and Experimental Data, in Computer Simulation of Biomolecular Systems, W. F. van Gunsteren, P. K. Weiner, A. J. Wilkinson, eds., Volume 3, 83-96 (1997)  [ff96]

Current force field parameter values and suggested procedures for development of parameters for additional molecules are available from the Amber web site in the Case lab at Scripps, http://amber.scripps.edu/

**AMBER98.PRM**

Amber ff98 parameters for proteins and nucleic acids. The only change from the ff94 parameter set is in the glycosidic torsional parameters that control sugar pucker.

T. E. Cheatham III, P. Cieplak and P. A. Kollman, A Modified Version of the Cornell et al. Force Field with Improved Sugar Pucker Phases and Helical Repeat, J. Biomol. Struct. Dyn., 16, 845-862 (1999)

Current force field parameter values and suggested procedures for development of parameters for additional molecules are available from the Amber web site in the Case lab at Scripps, http://amber.scripps.edu/

**AMBER99.PRM**

Amber ff99 parameters for proteins and nucleic acids. The original partial charges from the ff94 parameter set are retained, but many of the bond, angle and torsional parameters have been revised to provide better general agreement with experiment.

J. Wang, P. Cieplak and P. A. Kollman, How Well Does a Restrained Electrostatic Potential (RESP) Model Perform in Calcluating Conformational Energies of Organic and Biological Molecules?, J. Comput. Chem., 21, 1049-1074 (2000)

Current force field parameter values and suggested procedures for development of parameters for additional molecules are available from the Amber web site in the Case lab at Scripps, http://amber.scripps.edu/

**AMBER99SB.PRM**

Amber ff99SB parameters for proteins and nucleic acids from the Simmerling lab at Stonybrook University. Modifies the phi/psi and phi'/psi' torsional parameters of the ff99 protein force field.

V. Hornak, R. Abel, A. Okur, B. Strockbine, A. Roitberg and C. Simmerling, Comparison of Multiple Amber Force Fields and Development of Improved Protein Backbone Parameters, PROTEINS, 65, 712-725 (2006)  [PARM99SB]

**AMBER14SB.PRM**

Update Amber ff99SB protein parameters with improved backbone and sidechain torsional parameters. Also contains the updated OL15 DNA and OL3 RNA parameters values for nucleic acids.

J. A. Maier, C. Martinez, K. Kasavajhala, L. Wickstrom, K. E. House and C. Simmerling, ff14SB: Improving the Accuracy of Protein Side Chain
and Backbone Parameters from ff99SB, J. Chem. Theory Comput., 11, 2696-2713 (2015)  [Protein ff14SB]

R. Galindo-Murillo, J. C. Robertson, M. Zgarbovic, J. Sponer, M. Otyepka, P. Jureska, T. E. Cheatham, Assessing the Current State of Amber Force Field Modifications for DNA. J. Chem. Theory Comput., 12, 4114–4127 (2016)  [DNA OL15]

M. Zgarbova, M. Otyepka, J. Sponer, A. Mladek, P. Banas, T. E. Cheatham, P. Jurecka, Refinement of the Cornell et al. Nucleic Acids Force Field Based on Reference Quantum Chemical Calculations of Glycosidic Torsion Profiles, J. Chem. Theory Comput., 7, 2886–2902 (2011)  [RNA OL3]

**AMBER19SB.PRM**

Update Amber ff14SB protein parameters with amino acid-specific backbone CMAP torsion-torsion coupling values. Also contains the updated OL21 DNA and OL3 RNA parameters values for nucleic acids. Intended for use with the rigid four-point OPC water model.

C. Tian, K. Kasavajhala, K. Belfon, L. Raguette, H. Huang, A. Migues, J. Bickel, Y. Wang, J. Pincay, Q. Wu and C. Simmerling. ff19SB: Amino-Acid-Specific Protein Backbone Parameters Trained against Quantum Mechanics Energy Surfaces in Solution. J. Chem. Theory Comput., 16, 528–552 (2020)  [Protein ff19SB]

M. Zgarbova, J. Sponer and P. Jurecka, Z-DNA as a Touchstone for Additive Empirical Force Fields and a Refinement of the Alpha/Gamma DNA Torsions for AMBER, J. Chem. Theory Comput., 17, 6292–6301 (2017)  [DNA OL21]

M. Zgarbova, M. Otyepka, J. Sponer, A. Mladek, P. Banas, T. E. Cheatham and P. Jurecka, Refinement of the Cornell et al. Nucleic Acids Force Field Based on Reference Quantum Chemical Calculations of Glycosidic Torsion Profiles, J. Chem. Theory Comput., 7, 2886–2902 (2011)  [RNA OL3]

S. Izadi, R. Anandakrishnan and A. V. Onufriev, Building Water Models: A Different Approach, J. Phys. Chem. Lett., 5, 3863-3871 (2014)  [OPC Water]

A. Sengupta, Z. Li, L. F. Song, P. Li and K. M. Merz Jr., Parameterization of Monovalent Ions for the OPC3, OPC, TIP3P-FB, and TIP4P-FB Water Models, J. Chem. Inf. Model., 61, 869-880 (2021)  [Monovalent Ions]

Z. Li, L. F. Song, P. Li and K. M. Merz Jr., Systematic Parameterization of Divalent Metal Ions for the OPC3, OPC, TIP3P-FB, and TIP4P-FB Water Models, J. Chem. Theory Comput., 16, 4429-4442 (2020)  [Divalent Cations]

**AMOEBA09.PRM**

Parameters for the AMOEBA polarizable atomic multipole force field for a number of small organic molecules. These values were collected from 2009 and before, with values for each molecule developed individually via comparison against quantum mechanical and bulk phase data.

P. Ren and J. W. Ponder, A Consistent Treatment of Inter- and Intramolecular Polarization in Molecular Mechanics Calculations, J. Comput. Chem., 23, 1497-1506 (2002)

P. Ren and J. W. Ponder, Polarizable Atomic Multipole Water Model for Molecular Mechanics Simulation, J. Phys. Chem. B, 107, 5933-5947 (2003)

P. Ren and J. W. Ponder, Ion Solvation Thermodynamics from Simulation with a Polarizable Force Field, A. Grossfield, J. Am. Chem. Soc., 125, 15671-15682 (2003)

**AMOEBABIO09.PRM**

A combined biomolecular parameter file with 2009 AMOEBA biopolymer force field, with modifications to the original AMOEBA protein parameters made during 2009, and preliminary nucleic acid parameters due to Chuanjie Wu

A. Grossfield, P. Ren, J. W. Ponder, Ion Solvation Thermodynamics from Simulation with a Polarizable Force Field, J. Am. Chem. Soc., 125, 15671-15682 (2003)

P. Ren and J. W. Ponder, Polarizable Atomic Multipole Water Model for Molecular Mechanics Simulation, J. Phys. Chem. B, 107, 5933-5947 (2003)

**AMOEBABIO18.PRM**

An updated AMOEBA biomolecular parameter set based on the published 2013 protein parameters and 2018 nucleic acid parameters. Ion values taken from work by Zhi Wang, and the revised Generalized Kirkwood implicit solvent parameterization from Mike Schnieders group.

Y. Shi, Z. Xia, J. Zhang, R. Best, J. W. Ponder and P. Ren, Polarizable Atomic Multipole-Based AMOEBA Force Field for Proteins, J. Chem. Theory Comput., 9, 4046-4063 (2013)

C. Zhang, C. Lu, Z. Jing, C. Wu, J.-P. Piquemal, J. W. Ponder and P. Ren, AMOEBA Polarizable Atomic Multipole Force Field for Nucleic Acids, J. Chem. Theory Comput., 14, 2084-2108 (2018)

R. A. Corrigan, G. Qi, A. C. Thiel, J. R. Lynn, B. D. Walker, T. L. Casavant, L. Lagardere, J.-P. Piquemal, J. W. Ponder, P. Ren and M. J. Schnieders, Implicit Solvents for the Polarizable Atomic Multipole AMOEBA Force Field, J. Chem. Theory Comput., 17, 2323-2341 (2021)

**AMOEBABIO-HFC23.PRM**

Modified AMOEBA protein parameters with high field corrections (HFC) from Sameer Varma's group at the University of South Florida.

V. Wineman-Fisher, Y. Al-Hamdani, I. Addou, A. Tkatchenko and S. Varma, Ion-Hydroxyl Interactions: From High-Level Quantum Benchmarks to Transferable Polarizable Force Fields, J. Chem. Theory Comput., 15, 2444-2453 (2019)

V. Wineman-Fisher, J. M. Delgado, P. Nagy, E. Jakobsson, S. A. Pandit and S. Varma, Transferable Interactions of Li+ and Mg2+ Ions in Polarizable Models, J. Chem. Phys., 153, 104113 (2020)

V. Wineman-Fisher, Y. Al-Hamdani, P. Nagy, A. Tkatchenko and S. Varma, Improved Description of Ligand Polarization Enhances Transferability of Ion-Ligand Interactions, J. Chem. Phys., 153, 094115 (2020)

J. M. Delgado, V. Wineman-Fisher, S. A. Pandit and S. Varma, Inclusion of High-Field Target Data in AMOEBA’s Calibration Improves Predictions of Protein-Ion Interactions, J. Chem. Inf. Model., 62, 4713-4726 (2022)

J. M. Delgado, P. Nagy and S. Varma, Polarizable AMOEBA Model for Simulating Mg2+-Protein-Nucleotide Complexes, J. Chem. Inf. Model., 64, 378-392 (2024)

**AMOEBAIL22.PRM**

AMOEBA ionic liquid parameters from Andres Cisneros' group at the University of Texas at Dallas as provided in late 2022, and updated in May 2023.

**AMOEBAPRO04.PRM**

Initial protein parameters for the AMOEBA polarizable atomic multipole force field as determined in 2004.

J. W. Ponder and D. A. Case, Force Fields for Protein Simulation, Adv. Prot. Chem., 66, 27-85 (2003)

P. Ren and J. W. Ponder, Polarizable Atomic Multipole Water Model for Molecular Mechanics Simulation, J. Phys. Chem. B, 107, 5933-5947 (2003)

P. Ren, C. Wu and J. W. Ponder, Polarizable Atomic Multipole-Based Molecular Mechanics for Organic Molecules, J. Chem. Theory Comput., 7, 3143-3161 (2011)

**AMOEBAPRO13.prm**

Y. Shi, Z. Xia, J. Zhang, R. Best, J. W. Ponder and P. Ren, Polarizable Atomic Multipole-Based AMOEBA Force Field for Proteins, J. Chem. Theory Comput., 9, 4046-4063 (2013)

P. Ren, C. Wu and J. W. Ponder, Polarizable Atomic Multipole-Based Molecular Mechanics for Organic Molecules, J. Chem. Theory Comput., 7, 3143-3161 (2011)

J. C. Wu, J.-P. Piquemal, R. Chaudret, P. Reinhardt and P. Ren, Polarizable Molecular Dynamics Simulation of Zn(II) in Water Using the AMOEBA Force Field, J. Chem. Theory Comput., 6, 2059-2070 (2010)

A. Grossfield, P. Ren, J. W. Ponder, Ion Solvation Thermodynamics from Simulation with a Polarizable Force Field, J. Am. Chem. Soc., 125, 15671-15682 (2003)

P. Ren and J. W. Ponder, Polarizable Atomic Multipole Water Model for Molecular Mechanics Simulation, J. Phys. Chem. B, 107, 5933-5947 (2003)

**CHARMM19.PRM**

CHARMM19 united-atom parameters for proteins. The nucleic acid parameter are not yet implemented. There are some small differences between authentic CHARMM19 and the Tinker version due to: (1) replacement of CHARMM impropers by torsions for cases that involve atoms not bonded to the trigonal atom, and (2) Tinker's use of all possible torsions across a bond instead of a single torsion per bond.

E. Neria, S. Fischer and M. Karplus, Simulation of Activation Free Energies in Molecular Systems, J. Chem. Phys., 105, 1902-1921 (1996)

L. Nilsson and M. Karplus, Empirical Energy Functions for Energy Minimizations and Dynamics of Nucleic Acids, J. Comput. Chem., 7, 591-616 (1986)

W. E. Reiher III, Theoretical Studies of Hydrogen Bonding, Ph.D. Thesis, Department of Chemistry, Harvard University, Cambridge, MA, 1985

**CHARMM22.PRM**

CHARMM22 all-atom parameters for proteins and lipids. Many of the nucleic acid and small model compound parameters are not currently implemented.

N. Banavali and A. D. MacKerell, Jr., All-Atom Empirical Force Field for Nucleic Acids: 2) Application to Molecular Dynamics Simulations of DNA and RNA in Solution, J. Comput. Chem., 21, 105-120 (2000)

A. D. MacKerrell, Jr., et al., All-Atom Empirical Potential for Molecular Modeling and Dynamics Studies of Proteins, J. Phys. Chem. B, 102, 3586-3616 (1998)  [CHARMM22 Protein]

N. Foloppe and A. D. MacKerell, Jr., All-Atom Empirical Force Field for Nucleic Acids: 1) Parameter Optimization Based on Small Molecule and Condensed Phase Macromolecular Target Data, J. Comput. Chem., 21, 86-104 (2000)  [CHARMM Nucleic Acid]

A. D. MacKerell, Jr., J. Wiorkeiwicz-Kuczera and M. Karplus, An All-Atom Empirical Energy Function for the Simulation of Nucleic Acids, J. Am. Chem. Soc., 117, 11946-11975 (1995)

S. E. Feller, D. Yin, R. W. Pastor and A. D. MacKerell, Jr., Molecular Dynamics Simulation of Unsaturated Lipids at Low Hydration: Parametrization and Comparison with Diffraction Studies, Biophys. J., 73, 2269-2279 (1997)  [alkenes]

R. H. Stote and M. Karplus, Zinc Binding in Proteins and Solution - A Simple but Accurate Nonbonded Representation, Proteins, 23, 12-31 (1995)  [zinc ion]

Current parameter values are available from the CHARMM parameter site maintained by Alex MacKerell at the University of Maryland, Baltimore, http://mackerell.umaryland.edu/CHARMM_ff_params.html

**CHARMM27.PRM**

CHARMM27 is a modification of CHARMM22 to include CMAP torsion-torsion coupling parameters for the protein phi-psi backbone torsions.

A. D. MacKerell, Jr, M. Feig, C. L. Brooks III, "Extending the Treatment of Backbone Energetics in Protein Force Fields: Limitations of Gas-Phase Quantum Mechanics in Reproducing Protein Conformational Distributions in Molecular Dynamics Simulations", J. Comput. Chem., 25, 1400-1415 (2004)  [CMAP Torsion-Torsion Correction]

Current parameter values are available from the CHARMM parameter site maintained by Alex MacKerell at the University of Maryland, Baltimore, http://mackerell.umaryland.edu/CHARMM_ff_params.html

**CHARMM36M.PRM**

CHARMM36m is a major revision of the CHARMM parameters for proteins and nucleic acids. It includes the "36m" modifications intended to provide a balanced model for folded and disordered protein structures.

R. B. Best, X. Zhu, J. Shim, P. E. M. Lopes, J. Mittal, M. Feig and A. D. MacKerell Jr., "Optimization of the Additive CHARMM All-Atom Protein Force Field Targeting Improved Sampling of the Backbone Phi, Psi and Side-Chain Chi1 and Chi2 Dihedral Angles", J. Chem. Theory Comput., 8, 3257-3273 (2012)  [CHARMM36 Protein]

J. Huang and A. D. MacKerell Jr., "CHARMM36 All-Atom Additive Protein Force Field: Validation Based on Comparison to NMR Data", J. Comput. Chem., 34, 2135-2145 (2013)  [CHARMM36 Protein Validation]

J. Huang, S. Rauscher, G. Nawrocki, T. Ran, M. Feig, B. L. de Groot, H. Grubmuller and A. D MacKerell Jr., "CHARMM36m: An Improved Force Field for Folded and Intrinsically Disordered Proteins", Nature Methods, 14, 71-73 (2017)  [CHARMM36m Protein]

K. Hart, N. Foloppe, C. M. Baker, E. J. Denning, L. Nilsson and A. D. MacKerell Jr., "Optimization of the CHARMM Additive Force Field for DNA: Improved Treatment of the BI/BII Conformational Equilibrium", J. Chem. Theory Comput., 8, 348-362 (2012)  [CHARMM36 DNA]

E. J. Denning, U. D. Priyakumar, L. Nilsson and A. D. MacKerell Jr., "Impact of 2'-Hydroxyl Sampling on the Conformational Properties of RNA: Update of the CHARMM All-Atom Additive Force Field for RNA", J. Comput. Chem., 32, 1929-1943 (2011)  [CHARMM36 RNA]

Current parameter values are available from the CHARMM parameter site maintained by Alex MacKerell at the University of Maryland, Baltimore, http://mackerell.umaryland.edu/charmm-ff.shtml

**DANG.PRM**

L. X. Dang,  Development of Nonadditive Intermolecular Potentials Using Molecular Dynamics: Solvation of Li+ and F- Ions in Polarizable Water, J. Chem. Phys., 96, 6970-6977 (1992)

D. E. Smith and L. X. Dang, Interionic Potentials of Mean Force for SrCl2 in Polarizable Water: A Computer Simulation Study, Chem. Phys. Lett., 230, 209-214 (1994)

L. X. Dang and T.-M. Chang, Molecular Dynamics Study of Water Clusters, Liquid, and Liquid-Vapor Interface of Water with Many-Body Potentials, J. Chem. Phys., 106, 8149-8159 (1997)

T.-M. Chang and L. X. Dang, Ion Solvation in Polarizable Chloroform: A Molecular Dynamics Study, J. Phys. Chem. B, 101, 10518-10526 (1997)

T.-M. Chang and L. X. Dang, Detailed Study of Potassium Solvation Using Molecular Dynamics Techniques, J. Phys. Chem. B, 103, 4714-4720 (1999)

L. X. Dang, Computer Simulation Studies of Ion Transport Across a Liquid/Liquid Interface, J. Phys. Chem. B, 103, 8195-8200 (1999)

L. X. Dang, Intermolecular Interactions of Liquid Dichloromethane and Equilibrium Properties of Liquid-Vapor and Liquid-Liquid Interfaces: A Molecular Dynamics Study, J. Chem. Phys., 110, 10113-10122 (1999)

L. X. Dang and D. Feller, Molecular Dynamics Study of Water-Benzene Interactions at the Liquid/Vapor Interface of Water, J. Phys. Chem. B, 104, 4403-4407 (2000)

L. X. Dang, Molecular Dynamics Study of Benzene-Benzene and Benzene-Potassium Ion Interactions Using Polarizable Potential Models, J. Chem. Phys., 113, 266-273 (2000)

L. X. Dang, Computational Study of Ion Binding to the Liquid Interface of Water, J. Phys. Chem. B, 106, 10388-10394 (2002)

T.-M. Chang and L. X. Dang, On Rotational Dynamics of an NH4+ in Water, J. Chem. Phys., 118, 8813-8820 (2003)

L. X. Dang, Solvation of the Hydronium Ion at the Water Liquid/Vapor Interface, J. Chem. Phys., 119, 6351-6353 (2003)

L. X. Dang and T.-M. Chang, Many-Body Interactions in Liquid Methanol and its Liquid/Vapor Interface: A Molecular Dynamics Study, J. Chem. Phys., 119, 9851-9857 (2003)

L. X. Dang and B. C. Garrett, Molecular Mechanism of Water and Ammonia Uptake by the Liquid/Vapor Interface of Water, Chem. Phys. Lett., 385, 309-313 (2004)

L. X. Dang, Ions at the Liquid/Vapor Interface of Methanol, J. Phys. Chem. A, 108, 9014-9017 (2004)

M. Roeselova, J. Vieceli, L. X. Dang, B. C. Garrett and D. J. Tobias, Hydroxyl Radical at the Air-Water Interface, J. Am. Chem. Soc., 126, 16308-16309 (2004)

L. X. Dang, G. K. Schenter, V.-A. Glezakou and J. L. Fulton, Molecular Simulation Analysis and X-ray Absorption Measurement of Ca2+, K+ and Cl- Ions in Solution, J. Phys. Chem. B, 110, 23644-23654 (2006)

J. L. Thomas, M. Roeselova, L. X. Dang and D. J. Tobias, Molecular Dynamics Simulations of the Solution-Air Interface of Aqueous Sodium Nitrate, J. Phys. Chem. A, 111, 3091-3098 (2007)

C. D. Wick and L. X. Dang, Molecular Dynamics Study of Ion Transfer and Distribution at the Interface of Water and 1,2-Dichloroethane, J. Chem. Phys. C, 112, 647-649 (2008)

C. D. Wick and L. X. Dang, "Investigating Hydroxide Anion Interfacial Activity by Classical and Multistate Empirical Valence Bond Molecular Dynamics Simulations, J. Phys. Chem. A, 113, 6356-6364 (2009)

X. Sun, T.-M. Chang, Y. Cao, S. Niwayama, W. L. Hase and L. X. Dang, Solvation of Dimethyl Succinate in a Sodium Hydroxide Aqueous Solution. A Computational Study, J. Phys. Chem. B, 113, 6473-6477 (2009)

M. Baer, C. J. Mundy, T.-M. Chang, F.-M. Tao and L. X. Dang, Interpreting Vibrational Sum-Frequency Spectra of Sulfur Dioxide at the Air/Water Interface: A Comprehensive Molecular Dynamics Study, J. Phys. Chem. B, 114, 7245-7249 (2010)

L. X. Dang, T. B. Troung and B. Ginovska-Pangovska, Interionic Potentials of Mean Force for Ca+2-Cl- in Polarizable Water, J. Chem. Phys., 136, 126101 (2012)

**HIPPO19.PRM**

Initial parameters for the HIPPO force field for organic molecules, as derived starting in 2019 by Roseane dos Reis Silva in the Ponder lab.

J. A. Rackers and J. W. Ponder, Classical Pauli Repulsion: An Anisotropic, Atomic Multipole Model, J. Chem. Phys., 150, 084104 (2019)

J. A. Rackers, C. Liu, P. Ren and J. W. Ponder, A Physically Grounded Damped Dispersion Model with Particle Mesh Ewald Summation, J. Chem. Phys., 149, 084115 (2018)

J. A. Rackers, Q. Wang, C. Liu, J.-P. Piquemal, P. Ren and J. W. Ponder, An Optimized Charge Penetration Model for Use with the AMOEBA Force Field,  Phys. Chem. Chem. Phys., 19, 276-291 (2017)

**HOCH.PRM**

Simple NMR-NOE force field of Hoch and Stern.

J. C. Hoch and A. S. Stern, A Method for Determining Overall Protein Fold from NMR Distance Restraints, J. Biomol. NMR, 2, 535-543 (1992)

**IWATER.PRM**

Parameters for an AMOEBA-like polarizable multipole model for water including only "direct" (non-iterative) polarization.

L.-P. Wang, T. L. Head-Gordon, J. W. Ponder, P. Ren, J. D. Chodera, P. K. Eastman, T. J. Martinez and V. S. Pande, Systematic Improvement of a Classical Model of Water, J. Phys. Chem. B, 117, 9956-9972 (2013)

**MM2.PRM**

Full MM2(1991) parameters including pi-systems. The anomeric and electronegativity correction terms included in some later versions of MM2 are not implemented.

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

**MMFF94.PRM**

Nearly complete implementation of the MMFF94 force field, missing only the parameterization method for charged heteroaromatic ring systems.

T. A. Halgren, Merck Molecular Force Field. I. Basis, Form, Scope, Parametrization, and Performance of MMFF94, J. Comput. Chem., 17, 490-519 (1995)

T. A. Halgren, Merck Molecular Force Field. II. MMFF94 van der Waals and Electrostatic Parameters for Intermolecular Interactions, J. Comput. Chem., 17, 520-552 (1995)

T. A. Halgren, Merck Molecular Force Field. III. Molecular Geometries and Vibrational Frequencies for MMFF94, J. Comput. Chem., 17, 553-586 (1995)

T. A. Halgren and R. B. Nachbar, Merck Molecular Force Field. IV. Conformational Energies and Geometries for MMFF94, J. Comput. Chem., 17, 587-615 (1995)

T. A. Halgren, Merck Molecular Force Field. V. Extension of MMFF94 Using Experimental Data, Additional Computational Data, and Empirical Rules, J. Comput. Chem., 17, 616-641 (1995)

T. A. Halgren, MMFF VI. MMFF94s Option for Energy Minimization Studies, J. Comput. Chem., 20, 720-729 (1999)

T. A. Halgren, MMFF VII. Characterization of MMFF94, MMFF94s, and Other Widely Available Force Fields for Conformational Energies and for Intermolecular-Interaction Energies and Geometries, J. Comput. Chem., 20, 730-748 (1999)

**MMFF94S.PRM**

The MMFF94s parameter set enforces greater planarity at certain trigonal nitrogen atoms, including amides, anilines, ureas and others. It is meant to mimic the average bulk phase structure at these nitrogen centers. The changes from MMFF94 are only in out-of-plane and torsional parameters.

T. A. Halgren, MMFF VI. MMFF94s Option for Energy Minimization Studies, J. Comput. Chem., 20, 720-729 (1999)

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

**OPLSAA08.PRM**

OPLS-UA and OPLS-AA force fields with parameters for proteins and many general classes of organic molecules.

W. L. Jorgensen, D. S. Maxwell and J. Tirado-Rives, Development and Testing of the OPLS All-Atom Force Field on Conformational Energetics and Properties of Organic Liquids, J. Am. Chem. Soc., 117, 11225-11236 (1996)

D. S. Maxwell, J. Tirado-Rives and W. L. Jorgensen, A Comprehensive Study of the Rotational Energy Profiles of Organic Systems by Ab Initio MO Theory, Forming a Basis for Peptide Torsional Parameters, J. Comput. Chem., 16, 984-1010 (1995)

W. L. Jorgensen and N. A. McDonald, Development of an All-Atom Force Field for Heterocycles. Properties of Liquid Pyridine and Diazenes, THEOCHEM-J. Mol. Struct., 424, 145-155 (1998)

N. A. McDonald and W. L. Jorgensen, Development of an All-Atom Force Field for Heterocycles. Properties of Liquid Pyrrole, Furan, Diazoles, and Oxazoles, J. Phys. Chem. B, 102, 8049-8059 (1998)

R. C. Rizzo and W. L. Jorgensen, OPLS All-Atom Model for Amines: Resolution of the Amine Hydration Problem, J. Am. Chem. Soc., 121, 4827-4836 (1999)

M. L. P. Price, D. Ostrovsky and W. L. Jorgensen, Gas-Phase and Liquid-State Properties of Esters, Nitriles, and Nitro Compounds with the OPLS-AA Force Field, J. Comput. Chem., 22, 1340-1352 (2001)

The parameters distributed with Tinker are from "OPLS All-Atom Parameters for Organic Molecules, Ions, Peptides & Nucleic Acids, July 2008" as provided by W. L. Jorgensen, Yale University during June 2009

**OPLSAAL.PRM**

An improved OPLS-AA parameter set for proteins in which the only change is a reworking of many of the backbone and sidechain torsional parameters to give better agreement with LMP2 QM calculations. This parameter set is also known as OPLS(2000).

G. A. Kaminsky, R. A. Friesner, J. Tirado-Rives and W. L. Jorgensen, Evaluation and Reparametrization of the OPLS-AA Force Field for Proteins via Comparison with Accurate Quantum Chemical Calculations on Peptides, J. Phys. Chem. B, 105, 6474-6487 (2001)

**SMOOTH.PRM**

Version of OPLS-UA for use with potential smoothing. Largely adapted largely from standard OPLS-UA parameters with modifications to the vdw and improper torsion terms.

R. V. Pappu, R. K. Hart and J. W. Ponder, Analysis and Application of Potential Energy Smoothing and Search Methods for Global Optimization, J. Phys, Chem. B, 102, 9725-9742 (1998)  [smoothing modifications]

**SMOOTHAA.PRM**

Version of OPLS-AA for use with potential smoothing. Largely adapted largely from standard OPLS-AA parameters with modifications to the vdw and improper torsion terms.

R. V. Pappu, R. K. Hart and J. W. Ponder, Analysis and Application of Potential Energy Smoothing and Search Methods for Global Optimization, J. Phys, Chem. B, 102, 9725-9742 (1998)  [smoothing modifications]

**UWATER.PRM**

An AMOEBA-like water model with only a single polarizable multipole site centered on the water oxygen atom.

R. Qi, L.-P. Wang, Q. Wang, V. S. Pande, P. Ren, "United Polarizable Multipole Water Model for Molecular Mechanics Simulations", J. Chem. Phys., 143, 014504 (2015)

**WATER03.PRM**

The AMOEBA water parameters for a polarizable atomic multipole electrostatics model. This model from 2003 gives good agreement with ab initio and experimental data for many bulk and cluster properties.

P. Ren and J. W. Ponder, A Polarizable Atomic Multipole Water Model for Molecular Mechanics Simulation, J. Phys. Chem. B, 107, 5933-5947 (2003) 

P. Ren and J. W. Ponder, Ion Solvation Thermodynamics from Simulation with a Polarizable Force Field, A. Grossfield, J. Am. Chem. Soc., 125, 15671-15682 (2003)

P. Ren and J. W. Ponder, Temperature and Pressure Dependence of the AMOEBA Water Model, J. Phys. Chem. B, 108, 13427-13437 (2004)

An earlier version the AMOEBA water model is described in: Yong Kong, Multipole Electrostatic Methods for Protein Modeling with Reaction Field Treatment, Biochemistry & Molecular Biophysics, Washington University, St. Louis, August, 1997 [available from http://dasher.wustl.edu/ponder/]

**WATER14.PRM**

Revision of the AMOEBA 2003 water model based on fitting to high-level QM data on clusters via use of the Force Balance method of Leeping Wang at the University of California, Davis.

M. L. Laury, L.-P. Wang, V. S. Pande, T. Head-Gordon and J. W. Ponder, Revised Parameters for the AMOEBA Polarizable Atomic Multipole Water Model,  J. Phys. Chem. B., 119, 9423-9437 (2015)

**WATER21.PRM**

Water parameters for the HIPPO (Hydrogen-like Intermolecular Polarizable Potential) force field. It accounts for charge penetration, damped dispersion, and anisotropic repulsion as per the published HIPPO theory papers.

J. A. Rackers, R. R. Silva and J. W. Ponder, A Polarizable Water Model Derived from a Model Electron Density, J. Chem. Theory Comput., 17, 7056-7084 (2021)

J. A. Rackers and J. W. Ponder, Classical Pauli Repulsion: An Anisotropic, Atomic Multipole Model, J. Chem. Phys., 150, 084104 (2019)

J. A. Rackers, C. Liu, P. Ren and J. W. Ponder, A Physically Grounded Damped Dispersion Model with Particle Mesh Ewald Summation, J. Chem. Phys., 149, 084115 (2018)

J. A. Rackers, Q. Wang, C. Liu, J.-P. Piquemal, P. Ren and J. W. Ponder, An Optimized Charge Penetration Model for Use with the AMOEBA Force Field, Phys. Chem. Chem. Phys., 19, 276-291 (2017)

**WATER22.PRM**

An improved AMOEBA water model derived by Chengwen Liu at Univ. of Texas at Austin. Tunable parameters are vdw and bond/angle constants, the model provides
excellent agreement with experiment and high-quality QM data
