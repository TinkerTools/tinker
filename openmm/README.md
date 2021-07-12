# Tinker-OpenMM: Tinker Interface for OpenMM

<H2><B>Introduction</B></H2>

Tinker-OpenMM is an interface between Tinker and OpenMM. It provides an enhanced set of executables using Tinker as the "front end" while utilizing OpenMM as the "back end" to perform fast energy-force eveluations and molecular dynamics simulations, particularly on NVIDIA GPUs. Tinker-OpenMM is primarily intended for use on Linux systems, and the instructions below support that operating system. The same instructions can be used for Macintosh computers running macOS. However, due to continuing disfunction in the relationship between Apple and NVIDIA, the last macOS version supported is High Sierra (10.13.6), the last CUDA version supported is 10.2, and NVIDIA drivers are only available for GPUs through the Pascal (GTX 10xx) series. We have no experience with native builds of either OpenMM or Tinker-OpenMM on Windows-based systems.

<H2><B>The Tinker Package</B></H2>

The Tinker molecular modeling software is a complete and general package for molecular mechanics and dynamics, with some special features for biopolymers. Tinker has the ability to use any of several common parameter sets, such as Amber (ff94, ff96, ff98, ff99, ff99SB), CHARMM (19, 22, 22/CMAP), Allinger MM (MM2-1991 and MM3-2000), OPLS (OPLS-UA, OPLS-AA), Merck Molecular Force Field (MMFF), Liam Dang's polarizable model, and the AMOEBA, AMOEBA+ and HIPPO polarizable atomic multipole force fields. Updated parameter sets and those for other widely-used force fields will be released in due course.

<H2><B>Hardware & Software Prerequisites</B></H2>

Prior to building Tinker-OpenMM, you must have available a Linux system with standard development tools, including the GNU compiler suite V7 or later with gcc, g++ and gfortran. We routinely test on Ubuntu systems, as well as various versions of CentOS, but most current Linux distributions should be fine. The system needs to contain a functioning, fairly recent NVIDIA GPU with a working NVIDIA graphics driver. Commodity gaming GPUs from the 700 series onward should work, as will most of the Tesla series commercial/research GPU models. Next, a CUDA toolkit providing the full developemnt environment must be present. As of late 2020, we recommend installation of CUDA 11, though older versions back through CUDA 9.2 may work as well. We recommend installing the CUDA toolkit from the corresponding "run" file obtained directly from NVIDIA, instead of getting CUDA through the package manager for your Linux distribution. Installation from the run file will place CUDA into a /usr/local/cuda directory. Various of the setup files for building Tinker-OpenMM will look for CUDA under /usr/local.

<H2><B>Building & Installing OpenMM</B></H2>

Tinker-OpenMM requires several dynamic object libraries from the OpenMM package. Prior to building Tinker-OpenMM, you must first build the OpenMM version obtained from the current master repository at https://github.com/openmm/openmm. Download or clone this version of OpenMM, and then follow the standard instructions for building and installing OpenMM that are provided with that package. The basic procedure involves using cmake or ccmake to configure a Makefile, followed by "make" and "make install". The final installation should be in the default location, /usr/local/openmm. You can also try running "make tests", which will run a lengthy series of small examples and test cases to verify the OpenMM installation and its ability to access your GPU. Once OpenMM is behaving correctly, you can proceed to the final stage of building the Tinker-OpenMM executables to provide a Tinker interface to OpenMM.

<H2><B>Building Tinker-OpenMM from Source Code</B></H2>

Tinker is provided as a complete distribution available from https://github.com/TinkerTools/tinker or from the Ponder lab website at https://dasher.wustl.edu/tinker. After unpacking the release, you can build the Tinker-OpenMM executables following the steps below:

<B>(1)</B> We will assume you have unpacked Tinker in your home directory, ~/tinker. Other locations may require corresponding modification to the Makefile or other setup scripts.

<B>(2)</B> First, build the FFTW Fourier transform package that is included with the Tinker distribution. Move to the tinker/fftw directory, then follow the instructions in the 0README file. After "make" and "make install", there will be static libraries for FFTW (libfftw3.a, etc.) in tinker/fftw/lib. 

<B>(3)</B> Move into the top-level ~/tinker directory, and create a /build subdirectory.

<B>(4)</B> Move into the new /build subdirectory, and copy source and related files via the commands "cp ../source/<B>**</B>.f ."</B> and "cp ../openmm/* ."

<B>(5)</B> Check the directory environment variables near the top of the Makefile. As distributed, the section to build on Linux using the GNU compilers is activated.

<B>(6)</B> Run the "make" command to build executables. When the build is finished, there will be three executables produced in the build directory, "analyze_omm.x", "bar_omm.x" and "dynamic_omm.x".

These Tinker-OpenMM executables are GPU-enabled versions of the corresponding Tinker programs, and use OpenMM libraries to perform energy and force evaluations, as well as to take MD steps. Tinker-OpenMM supports the standard Tinker coordinate and parameter files, as well as many of the usual keyfile options. The executables can be renamed or moved elsewhere, but are dynamically linked against CUDA and OpenMM, which must remain available in the same location as when the build was performed.
