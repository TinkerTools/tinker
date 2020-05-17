Introduction
============

What is Force Field Explorer?
-----------------------------

Force Field Explorer (FFE) is a graphical user interface to the Tinker suite of molecular modeling tools. The goal is to create a robust environment for computational chemistry and molecular engineering applications. Current features include:

* Setup and real-time visualization of modeling routines, including local/global optimizations, sampling methods, superposition, etc. 
* Basic 3D molecular visualization (wireframe, tube, spacefilling, ball & stick, etc.) with special features for force field specific parameters such as Lennard-Jones radii and partial charge magnitudes. 
* A comprehensive editor for Tinker keyword files, which contain parameters that control modeling algorithms. 
* Tinker trajectory playback. 

Current Release
---------------

The current release of Force Field Explorer is denoted as Version 8.8 and dated July 2020 to mirror the corresponding Tinker version. Dependencies, listed below, are now included in self-extracting installers for Linux, MacOS and Microsoft Windows.

* Tinker version 8.8
* Java Runtime Environment (JRE) 1.8 (Java 8) or later 
* Java3D version 1.6

Major improvements over the previous releases of FFE include:

* Self-Extracting Installers for Linux, MacOS and Microsoft Windows. These distributions include Tinker binaries, source code, documentation, and example files. Private JRE and Java3D libraries are bundled with all versions of FFE.
* An intuitive, consistent selection mechanism. 
* Improved keyword and modeling panels. 
* Support for Tinker internal coordinate, restart and induced dipole file formats, plus loading of Protein Databank files directly from the PDB or from disk. Conversion of PDB files to Tinker format (with automatic atom typing for many biomolecules) is also supported.

Future Plans
------------

The most important Force Field Explorer development goal is to add further molecular editing and visualization features. Currently, Java3D is used as the graphics toolkit underlying FFE. Java3D was originally developed by Sun, but is now maintained by an open source consortium and is freely available. Other goals include support for other file types and special capabilities for molecular complexes, docking and binding.
