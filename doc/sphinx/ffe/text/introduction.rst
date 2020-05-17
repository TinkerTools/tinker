Introduction
============

What is Force Field Explorer?
-----------------------------

Force Field Explorer (FFE) is a graphical user interface to the Tinker suite of software tools for molecular modeling. The goal is to create a robust environment for computational chemistry, biomolecular modeling and molecular engineering applications. Current features of FFE include:

* Setup and real-time visualization of modeling routines, including local/global optimizations, molecular dynamics, sampling methods, superposition, docking, etc. 
* Basic 3D molecular visualization (wireframe, tube, spacefilled, ball & stick, etc., and several coloring schemes) with special features for force field-specific parameters such as Lennard-Jones radii and partial charges, induced dipole vectors, and atomic force and velocity vectors. 
* A comprehensive editor for Tinker keyword files containing parameters and options that control modeling algorithms. 
* Playback of Tinker molecular dynamics trajectories and structure archives.

Current Release
---------------

The current release of Force Field Explorer is denoted as Version 8.8 and dated July 2020 to mirror the corresponding Tinker 8.8 version. FFE is self-contained in that all dependencies, listed below, are included in install packages available for Linux, MacOS and Microsoft Windows operating systems.

* Tinker version 8.8
* Java Runtime Environment (JRE) 1.8 (Java 8) or later 
* Java3D version 1.6

Major improvements over the previous releases of FFE include:

* Self-Extracting Installers for Linux, MacOS and Microsoft Windows. The installer distributions include FFE-enabled Tinker binaries, full source code, FFE and Tinker documentation, and example files. Private JRE and Java3D libraries are bundled with all versions of FFE.
* An intuitive, consistent selection mechanism. 
* Improvements to the keyword and modeling panels. 
* Support for Tinker internal coordinate, restart and induced dipole file formats, plus of molecules and biopolymers directly from the PubChem, NCI or PDB databases,  or from disk files.
* Automatic conversion of downloaded small molecules to Tinker's "tiny" force field. Also supported is convertion of PDB files to Tinker format with automatic atom typing for many standard biomolecular force fields.

Future Plans
------------

The most important Force Field Explorer development is to add further molecular editing and visualization features, as well as tighter integration with Tinker. Currently, Java3D is used as the graphics toolkit underlying FFE's visualization engine. Java3D was originally developed by Sun, but is now maintained by an open source consortium and is freely available. Other planned enhancements include support for additional molecule file types and special capabilities for molecular complexes, docking and binding applications.
