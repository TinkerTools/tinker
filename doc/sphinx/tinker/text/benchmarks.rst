Benchmark Results
=================

The tables in this section provide CPU benchmarks for basic Tinker energy and derivative evaluations, vibrational analysis and molecular dynamics simulation. All times are in seconds and were measured with Tinker executables dimensioned to a maximum of 1000000 atoms. Each benchmark was run on an unloaded machine and is the fastest time reported for that particular machine. The first five benchmarks are run serial on a single thread, while the last four benchmarks reflect OpenMP parallel performance. If you have built Tinker on an alternative machine type and are able to run the benchmarks on the additional machine type, please send the results for inclusion in a future listing.

Calmodulin Energy Evaluation (Serial)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gas-Phase Calmodulin Molecule, 2264 Atoms, Amber ff94 Force Field, No Nonbonded Cutoffs, 100 Evaluations

**Machine Type (OS/Compiler)                           CPU        Energy       Gradient       Hessian**

.. code-block:: text

 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150        26.2         50.9       149.9
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462        18.2         34.2       106.2
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650        13.8         34.6        97.4
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670         7.2         16.6        50.4
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ         7.3         16.7        51.5
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H         5.4         10.4        33.2

Crambin Crystal Energy Evaluation (Serial)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Crambin Unit Cell, 1360 Atoms in Periodic Unit Cell, OPLS-UA Force Field with PME Electrostatics, 9.0 Ang vdw Cutoff, 1000 Evaluations

**Machine Type (OS/Compiler)                           CPU        Energy       Gradient       Hessian**

.. code-block:: text

 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150        36.6         54.1       220.5
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462        32.0         45.4       212.7
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650        30.8         45.1       171.2
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670        19.1         26.0        85.0
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ        18.8         26.6        87.8
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H        15.6         19.9        67.9

Crambin Normal Mode Calculation (Serial)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Hessian Eigenvalues, Normal Modes and Vibrational Freqencies for the 42-Amino Acid, 642-Atom Protein Crambin, CHARMM-22 Force Field with Cutoffs

**Machine Type (OS/Compiler)                           CPU                   Seconds**

.. code-block:: text

 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150                41.1
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462                39.6
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650                29.4
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670                15.9
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ                15.3
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H                11.0

Water Box Molecular Dynamics using TIP3P (Serial)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MD run of 10000 Steps for 216 TIP3P Waters in 18.643 Ang Periodic Box, 9.0 Ang Shifted & Switched Cutoffs, Rattle for Rigid TIP3P, 1.0 fs Time Step with Modified Beeman Integrator

**Machine Type (OS/Compiler)                           CPU                   Seconds**

.. code-block:: text

 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150               214.0
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462               156.7
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650               163.2
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670                80.4
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ                82.1
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H                55.1

Water Box Molecular Dynamics using AMOEBA (Serial)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MD run of 1000 Steps for 216 AMOEBA Waters in a 18.643 Ang Box, Neighbor Lists, PME with a 20x20x20 FFT and 7.0 Ang Real-Space Cutoff, 9.0 Ang vdW Cutoff with Correction, 1.0 fs Time Step with Modified Beeman Integrator, and 0.00001 RMS Induced Dipole Convergence

**Machine Type (OS/Compiler)                           CPU                   Seconds**

.. code-block:: text

 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150               112.8
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462                97.4
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650                86.6
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670                46.9
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ                48.7
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H                37.0

MD on DHFR in Water using CHARMM (OpenMP Parallel)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MD run of 100 Steps for CHARMM DHFR in Water (23558 Atoms, 62.23 Ang Box), Neighbor Lists, PME with a 64x64x64 FFT and 7.0 Ang Real-Space Cutoff, 9.0 Ang vdW Cutoff, 1.0 fs Time Step with Modified Beeman Integrator; OpenMP timings as "wall clock" time, with parallel speedup in parentheses

**Machine Type (OS/Compiler)                           CPU       Core/Thread          Seconds**

.. code-block:: text

 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150           4/1       105.8 (1.00)
 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150           4/4        68.1 (1.55)
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462           8/1        86.3 (1.00)
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462           8/8        49.0 (1.76)
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650          12/1        79.1 (1.00)
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650         12/24        36.8 (2.15)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/1        44.7 (1.00)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/2        32.5 (1.38)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/4        25.0 (1.79)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/1        45.4 (1.00)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/2        34.3 (1.32)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/4        26.9 (1.69)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/8        23.8 (1.91)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/1        34.5 (1.00)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/2        23.6 (1.46)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/4        16.8 (2.05)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/6        14.5 (2.38)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H          6/12        14.1 (2.45)

MD on DHFR in Water using AMOEBA (OpenMP Parallel)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MD run of 100 Steps for AMOEBA DHFR in Water (23558 Atoms, 62.23 Ang Box), Neighbor Lists, PME with a 64x64x64 FFT and 7.0 Ang Real-Space Cutoff, 9.0 Ang vdW Cutoff with Correction, 1.0 fs Time Step with Modified Beeman Integrator, and 0.00001 RMS Induced Dipole Convergence; OpenMP timings reported as "wall clock" time, with parallel speedup in parentheses

**Machine Type (OS/Compiler)                           CPU       Core/Thread          Seconds**

.. code-block:: text

 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150           4/1       507.5 (1.00)
 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150           4/4       246.3 (2.06)
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462           8/1       423.4 (1.00)
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462           8/8       161.1 (2.63)
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650          12/1       384.8 (1.00)
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650         12/24       122.0 (3.15)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/1       202.6 (1.00)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/2       143.5 (1.41)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/4        92.2 (2.20)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/1       219.1 (1.00)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/2       159.8 (1.37)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/4        98.6 (2.22)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/8        85.0 (2.58)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/1       153.3 (1.00)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/2       108.0 (1.42)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/4        69.7 (2.20)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/6        57.1 (2.68)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H          6/12        59.5 (2.58)

MD on COX-2 in Water using OPLS-AA (OpenMP Parallel)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MD run of 100 Steps for OPLS-AA COX-2 in Water (174219 Atoms, 120.0 Ang Box), Neighbor Lists, PME with a 128x128x128 FFT and 7.0 Ang Real-Space Cutoff, 9.0 Ang vdW Cutoff, 1.0 fs Time Step with Modified Beeman Integrator; RATTLE for all X-H bonds and rigid TIP3P Water; OpenMP timings reported as "wall clock" time, with parallel speedup in parentheses

**Machine Type (OS/Compiler)                           CPU       Core/Thread          Seconds**

.. code-block:: text

 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150           4/1        798.6 (1.00)
 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150           4/4        487.2 (1.65)
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462           8/1        666.6 (1.00)
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462           8/8        367.1 (1.82)
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650          12/1        531.9 (1.00)
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650         12/24        267.0 (1.99)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/1        344.2 (1.00)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/2        270.9 (1.27)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/4        196.8 (1.75)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/1        347.6 (1.00)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/2        276.1 (1.26)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/4        193.7 (1.79)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/8        173.2 (2.01)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/1        262.3 (1.00)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/2        202.4 (1.30)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/4        143.2 (1.83)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/6        127.2 (2.06)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H          6/12        120.9 (2.17)

MD on COX-2 in Water using AMOEBA (OpenMP Parallel)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MD run of 100 Steps for AMOEBA COX-2 in Water (174219 Atoms, 120.0 Ang Box), Neighbor Lists, PME with a 128x128x128 FFT and 7.0 Ang Real-Space Cutoff, 9.0 Ang vdW Cutoff with Correction, 1.0 fs Time Step with Modified Beeman Integrator, and 0.00001 RMS Induced Dipole Convergence; OpenMP timings reported as "wall clock" time, with parallel speedup in parentheses

**Machine Type (OS/Compiler)                           CPU       Core/Thread          Seconds**

.. code-block:: text

 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150           4/1      5427.4 (1.00)
 Mac Pro 1.1 (MacOS 10.11, GNU 7.1)             5150           4/4      2369.3 (2.29)
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462           8/1      4838.9 (1.00)
 Mac Pro 3.1 (MacOS 10.13, GNU 8.1)            E5462           8/8      1661.1 (2.91)
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650          12/1      3686.8 (1.00)
 Mac Pro 5.1 (MacOS 10.13, GNU 8.1)            X5650         12/24       933.3 (3.95)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/1      2240.2 (1.00)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/2      1509.2 (1.00)
 iMac 14.2 (MacOS 10.13, GNU 8.1)               4670           4/4       916.8 (1.00)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/1      2279.8 (1.00)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/2      1494.0 (1.00)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/4       897.6 (1.00)
 MacBook Pro 11.3 (MacOS 10.13, GNU 8.1)      4960HQ           4/8       763.5 (1.00)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/1      1621.2 (1.00)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/2      1114.9 (1.00)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/4       701.3 (1.00)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H           6/6       577.4 (1.00)
 Razer Blade 15 (Ubuntu 18.04, GNU 7.5)        9750H          6/12       545.2 (1.00)
