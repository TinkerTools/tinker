Benchmark Results
=================

The tables in this section provide CPU benchmarks for basic Tinker energy and derivative evaluations, vibrational analysis and molecular dynamics. All times are in seconds and were measured with Tinker executables dimensioned to maxatm of 10000 and maxhess of 1000000 in the source file sizes.i. All calculations were run twice in rapid succession on a quiet machine. The times reported for each benchmark are the results from the second run. If you have built Tinker on an alternative machine type and are able to run the benchmarks on the additional machine type, please send the results for inclusion in a future listing.

BENCHMARK #1:  Calmodulin Energy Evaluation

The system is an isolated molecule of the 148-residue protein calmodulin with 2264 atoms using the Amber ff94 force field. All interactions are computed with no use of cutoffs. Times listed are for calculation setup followed by a single energy, energy/gradient and Hessian evaluation.

MACHINE-OS-COMPILER TYPE	MHz	SETUP	ENERGY	GRAD	HESS

 Athlon XP 2400+ (RH 8.0, Intel)	2000	0.13	0.28	0.60	2.96
 Athlon XP 2400+ (RH 8.0, PGI)	2000	0.16	0.31	0.70	3.60
 Athlon XP 2400+ (RH 8.0, g77 3.2)	2000	0.17	0.28	0.66	3.67
 Athlon Thunderbird (RH 8.0, Intel)	1400	0.22	0.41	0.86	5.15
 Athlon Thunderbird (RH 8.0, PGI)	1400	0.21	0.44	1.00	5.92
 Athlon Thunderbird (RH 8.0, g77 3.2)	1400	0.19	0.40	0.94	5.81
 Athlon Classic (RH 8.0, Intel)	950	0.30	0.64	1.42	7.07
 Athlon Classic (RH 8.0, PGI)	950	0.30	0.69	1.65	7.96
 Athlon Classic (RH 8.0, g77 3.2)	950	0.31	0.63	1.57	7.94
 Compaq Evo N610c P4 (RH 8.0, Intel)	2000	0.18	0.45	0.87	3.08
 Compaq Evo N610c P4 (RH 8.0, PGI)	2000	0.22	0.44	1.06	4.27
 Compaq Evo N610c P4 (RH 8.0, Absoft)	2000	0.17	0.52	1.06	3.95
 Compaq Evo N610c P4 (RH 8.0, g77 3.2)	2000	0.19	0.41	1.07	4.41
 Compaq Evo N610c P4 (WinXP, CVF 6.6)	2000	0.16	0.38	0.98	3.54
 Compaq Evo N610c P4 (WinXP, g77 3.2)	2000	0.16	0.40	1.08	4.45
 Apple Power Mac G4 (OSX 10.2, Absoft)	733	0.41	2.96	5.12	17.83
 Apple Power Mac G4 (OSX 10.2, g77 3.3)	733	0.37	1.98	3.79	14.48
 Compaq AlphaServer DS10 (Tru64 5.0)	466	0.35	1.33	1.93	8.40
 SGI IndigoII R10K (Irix 6.5, MIPS)	195	1.17	3.49	6.35	23.03

BENCHMARK #2:  Crambin Crystal Energy Evaluation

The system is a unit cell of the 46-residue protein crambin containing 2 polypeptide chains, 2 ethanol and 178 water molecules for a total of 1360 atoms using the OPLS-UA force field. Periodic boundaries are used with particle mesh Ewald for electrostatics and a 9.0 Angstrom cutoff for vdW interactions. Times listed are for calculation setup followed by a single energy, energy/ gradient and Hessian evaluation.

MACHINE-OS-COMPILER TYPE	MHz	SETUP	ENERGY	GRAD	HESS

 Athlon XP 2400+ (RH 8.0, Intel)	2000	0.12	0.12	0.21	0.66
 Athlon XP 2400+ (RH 8.0, PGI)	2000	0.13	0.14	0.24	0.63
 Athlon XP 2400+ (RH 8.0, g77 3.2)	2000	0.14	0.13	0.28	0.81
 Athlon Thunderbird (RH 8.0, Intel)	1400	0.19	0.17	0.30	0.91
 Athlon Thunderbird (RH 8.0, PGI)	1400	0.18	0.18	0.32	0.91
 Athlon Thunderbird (RH 8.0, g77 3.2)	1400	0.17	0.17	0.38	1.11
 Athlon Classic (RH 8.0, Intel)	950	0.26	0.25	0.47	1.46
 Athlon Classic (RH 8.0, PGI)	950	0.29	0.27	0.50	1.42
 Athlon Classic (RH 8.0, g77 3.2)	950	0.27	0.27	0.56	1.70
 Compaq Evo N610c P4 (RH 8.0, Intel)	2000	0.15	0.14	0.27	0.64
 Compaq Evo N610c P4 (RH 8.0, PGI)	2000	0.22	0.19	0.33	0.88
 Compaq Evo N610c P4 (RH 8.0, Absoft)	2000	0.14	0.22	0.39	0.84
 Compaq Evo N610c P4 (RH 8.0, g77 3.2)	2000	0.15	0.20	0.45	1.13
 Compaq Evo N610c P4 (WinXP, CVF 6.6)	2000	0.14	0.17	0.33	0.83
 Compaq Evo N610c P4 (WinXP, g77 3.2)	2000	0.12	0.22	0.52	1.16
 Apple Power Mac G4 (OSX 10.2, Absoft)	733	0.32	0.58	1.09	3.11
 Apple Power Mac G4 (OSX 10.2, g77 3.3)	733	0.31	0.42	0.79	2.37
 Compaq AlphaServer DS10 (Tru64 5.0)	466	0.29	0.38	0.64	1.95
 SGI IndigoII R10K (Irix 6.5, MIPS)	195	0.92	0.74	1.41	3.89

BENCHMARK #3:  Peptide Normal Mode Calculation

The system is a minimum energy conformation of a 20-residue peptide containing one of each of the standard amino acids for a total of 328 atoms using the OPLS-AA force field without cutoffs. The time reported is for computation of the Hessian and calculation of the normal modes of the Hessian matrix and the vibration frequencies requiring two separate matrix diagonalization steps.

MACHINE-OS-COMPILER TYPE	MHz	NORMAL MODES

 Athlon XP 2400+ (RH 8.0, Intel)	2000	22
 Athlon XP 2400+ (RH 8.0, PGI)	2000	26
 Athlon XP 2400+ (RH 8.0, g77 3.2)	2000	24
 Athlon Thunderbird (RH 8.0, Intel)	1400	31
 Athlon Thunderbird (RH 8.0, PGI)	1400	34
 Athlon Thunderbird (RH 8.0, g77 3.2)	1400	33
 Athlon Classic (RH 8.0, Intel)	950	46
 Athlon Classic (RH 8.0, PGI)	950	51
 Athlon Classic (RH 8.0, g77 3.2)	950	48
 Compaq Evo N610c P4 (RH 8.0, Intel)	2000	19
 Compaq Evo N610c P4 (RH 8.0, PGI)	2000	19
 Compaq Evo N610c P4 (RH 8.0, Absoft)	2000	20
 Compaq Evo N610c P4 (RH 8.0, g77 3.2)	2000	19
 Compaq Evo N610c P4 (WinXP, CVF 6.6)	2000	19
 Compaq Evo N610c P4 (WinXP, g77 3.2)	2000	20
 Apple Power Mac G4 (OSX 10.2, Absoft)	733	67
 Apple Power Mac G4 (OSX 10.2, g77 3.3)	733	62
 Compaq AlphaServer DS10 (Tru64 5.0)	466	39
 SGI IndigoII R10K (Irix 6.5, MIPS)	195	144

BENCHMARK #4:  TIP3P Water Box Molecular Dynamics

The system consists of 216 rigid TIP3P water molecules in a 18.643 Angstrom periodic box, 9.0 Angstrom shifted energy switch cutoffs for nonbonded interactions. The time reported is for 1000 dynamics steps of 1.0 fs each using the modified Beeman integrator and Rattle constraints on all bond lengths.

MACHINE-OS-COMPILER TYPE	MHz	DYNAMICS

 Athlon XP 2400+ (RH 8.0, Intel)	2000	37
 Athlon XP 2400+ (RH 8.0, PGI)	2000	34
 Athlon XP 2400+ (RH 8.0, g77 3.2)	2000	45
 Athlon Thunderbird (RH 8.0, Intel)	1400	52
 Athlon Thunderbird (RH 8.0, PGI)	1400	47
 Athlon Thunderbird (RH 8.0, g77 3.2)	1400	63
 Athlon Classic (RH 8.0, Intel)	950	77
 Athlon Classic (RH 8.0, PGI)	950	71
 Athlon Classic (RH 8.0, g77 3.2)	950	96
 Compaq Evo N610c P4 (RH 8.0, Intel)	2000	53
 Compaq Evo N610c P4 (RH 8.0, PGI)	2000	54
 Compaq Evo N610c P4 (RH 8.0, Absoft)	2000	55
 Compaq Evo N610c P4 (RH 8.0, g77 3.2)	2000	91
 Compaq Evo N610c P4 (WinXP, CVF 6.6)	2000	63
 Compaq Evo N610c P4 (WinXP, g77 3.2)	2000	94
 Apple Power Mac G4 (OSX 10.2, Absoft)	733	209
 Apple Power Mac G4 (OSX 10.2, g77 3.3)	733	170
 Compaq AlphaServer DS10 (Tru64 5.0)	466	106
 SGI IndigoII R10K (Irix 6.5, MIPS)	195	280

BENCHMARK #5:  Tinker Water Box Molecular Dynamics

The system consists of 216 AMOEBA flexible polarizable atomic multipole water molecules in a 18.643 Angstrom periodic box using regular Ewald summation for the electrostatics and a 12.0 Angstrom switched cutoff for vdW interactions. The time reported is for 100 dynamics steps of 1.0 fs each using the modified Beeman integrator and 0.01 Debye rms convergence for induced dipole moments.

MACHINE-OS-COMPILER TYPE	MHz	DYNAMICS

 Athlon XP 2400+ (RH 8.0, Intel)	2000	108
 Athlon XP 2400+ (RH 8.0, PGI)	2000	104
 Athlon XP 2400+ (RH 8.0, g77 3.2)	2000	128
 Athlon Thunderbird (RH 8.0, Intel)	1400	165
 Athlon Thunderbird (RH 8.0, PGI)	1400	158
 Athlon Thunderbird (RH 8.0, g77 3.2)	1400	183
 Athlon Classic (RH 8.0, Intel)	950	282
 Athlon Classic (RH 8.0, PGI)	950	261
 Athlon Classic (RH 8.0, g77 3.2)	950	307
 Compaq Evo N610c P4 (RH 8.0, Intel)	2000	156
 Compaq Evo N610c P4 (RH 8.0, PGI)	2000	191
 Compaq Evo N610c P4 (RH 8.0, Absoft)	2000	226
 Compaq Evo N610c P4 (RH 8.0, g77 3.2)	2000	243
 Compaq Evo N610c P4 (WinXP, CVF 6.6)	2000	176
 Compaq Evo N610c P4 (WinXP, g77 3.2)	2000	263
 Apple Power Mac G4 (OSX 10.2, Absoft)	733	680
 Apple Power Mac G4 (OSX 10.2, g77 3.3)	733	479
 Compaq AlphaServer DS10 (Tru64 5.0)	466	358
 SGI IndigoII R10K (Irix 6.5, MIPS)	195	868
