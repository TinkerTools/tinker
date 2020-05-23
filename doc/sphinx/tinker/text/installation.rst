Installation on Your Computer
=============================

How to Obtain a Copy of Tinker
------------------------------

The Tinker package is distributed on the Internet at the Ponder lab's Tinker web site located at https://dasher.wustl.edu/tinker/, or via download from the Github site for the TinkerTools organization at https://github.com/TinkerTools/Tinker/. After unpacking the distribution, you can build a set of Tinker executables on almost any machine with a Fortran compiler. Makefiles, a CMakeLists.txt file for cmake, as well as standalone scripts to compile, build object libraries, and link executables on a wide variety of machine-CPU-operating system combinations are provided.

Prebuilt Tinker Executables
---------------------------

Pre-built Tinker executables for Linux, MacOS, and Windows are also available for download from the sites mentioned above. They should run on most recent vintage machines using the above operating systems, and can handle a maximum of 1 million atoms provided sufficient memory is available. The Linux executables require at least glibc-2.6 or later. Note starting with Tinker 8, we no longer provide pre-built executables for any 32-bit operating systems.

The provided executables are OpenMP capable, but do not support APBS or the Tinker-FFE interface. You will still need to have a copy of the complete Tinker distribution as it contains the parameter sets, examples, benchmarks, test files and documentation required to use the package.

Building your Own Executables
-----------------------------

The compilation and building of the Tinker executables should be easy for most of the common Linux, MacOS and Windows computers. We provide in the /make area of the distribution a Makefile that with minor modification can be used to build Tinker on any of these machines. As an alternative to Makefiles, we also provide machine-specific directories with three separate shell scripts to compile the source, build an object library, and link binary executables.

The first step in building Tinker using the script files is to run the appropriate compile.make script for your operating system and compiler version. Next you must use a library.make script to create an archive of object code modules. Finally, run a link.make script to produce the complete set of Tinker executables. The executables can be renamed and moved to wherever you like by editing and running the "rename" script. These steps will produce executables that can run from the command line, but without the capability to interact with the FFE GUI. Building FFE-enabled Tinker executables involves replacing the sockets.f source file with sockets.c, and included the object from the C code in the Tinker object library. Then executables must be linked against Java libraries in addition to the usual resources. Sample compgui.make and linkgui.make scripts are provided for systems capable of building GUI-enabled executables.

Regardless of your target machine, only a few small pieces of code can possibly require attention prior to building. The most common source alterations are to the master array dimensions found in the source file sizes.f. The basic limit is on the number of atoms allowed, "maxatm". This parameter can be set to 1000000 or more on most workstations. Personal computers with minimal memory may need a lower limit, depending on available memory, swap space and other resources. A description of the other parameter values is contained in the header of the file.

Tinker-FFE (Force Field Explorer)
---------------------------------

Tinker-FFE, formerly Force Field Explorer, is a Java-based GUI for the Tinker package. It provides visualization for Tinker molecule files, as well as launching of Tinker calculations from a graphical interface. Tinker-FFE for Linux, MacOS and Windows can be downloaded from the Ponder lab Tinker web site as "installation kits" containing the FFE GUI and an FFE-enabled version of Tinker. Tinker-FFE requires a 64-bit CPU and operating system, as 32-bit systems are no longer supported.

Integration with Tinker, including the ability to interactively run Tinker calculations, and to access molecule downloads from the PubChem, NCI and PDB databases make Tinker-FFE a useful tool in classroom teaching environments. For research work, we recommend using the latest command line version of Tinker for numerical calculations, and using FFE or another visualization program to view results. Several other visualization programs (including VMD, Avogadro, Jmol, MOLDEN, WebMO, some PyMOL versions, etc.) can display Tinker structure and MD trajectory files.

The Tinker-FFE Installer for Linux is provided as a gzipped shell script. Uncompress the the .gz archive to produce an .sh script, and then run the script. The script must have the "executable" attribute, set via "chmod +x installer-file-name.sh", prior to being run.

The Tinker-FFE Installer for MacOS is provided as a .dmg disk image file. Double-click on the file to run the installer. MacOS 10.8 and later contains a security feature called Gatekeeper that keeps applications not obtained via the App Store or Apple-approved developers from being opened. Gatekeeper is enabled by default, and may result in the (incorrect!) error message: "Tinker-FFE Installer.app is damaged and can't be opened." To turn off Gatekeeper, go to the panel System Preferences > Security & Privacy > General, and set "Allow apps downloaded from:" to "Anywhere". This will require an Administrator account, and must be done before invoking the FFE installer. Once FFE is installed and launched for the first time, you can return the System Preference to its prior value. On Sierra (10.12) and later, the "Anywhere" option has been removed. In most cases the Security & Privacy panel will open and permit the user to run the installer. Alternatively, the "Anywhere" option can be restored by running the command "sudo spctl --master-disable" in a Terminal window.

The Tinker-FFE Installer for Windows is provided as a zipped executable. First, unzip the .zip file, then run the resulting executable .exe file. In order to perform minimizations or molecular dynamics from within FFE, some environment variables and symbolic links must be set prior to using the program. A batch file named "FFESetupWin.bat" is installed in the main Tinker-FFE directory, which by default resides in the user's home directory. To complete the setup of FFE, this batch file should be run from a Command Prompt window following installation. It is only necessary to invoke this batch file once, as the settings should persist between logins.

For those wishing to modify the FFE GUI or build a version from source, we provide a complete development package for Tinker-FFE. This is a large download which contains the code for all components, including the Java source for FFE itself and the many required Java libraries. This package allows building Tinker-FFE on all three supported operating systems from a common code base. External requirements are the GNU compiler suite with gcc, g++ and gfortran (on Windows use MinGW-w64 compilers under Cygwin), and the Install4j Java installer builder. Note Install4j is a commercial product; only the compiler is needed, not the full Install4j GUI interface.

Documentation and Other Information
-----------------------------------

The documentation for the Tinker programs, including the guide you are currently reading, is located in the /doc subdirectory of the distribution. The documentation was prepared using the Sphinx documentation generator. Portable versions of the documentation are provided as PDF files and in HTML format for web display. Please read and return by mail the Tinker license. In particular, we note that Tinker is not an Open Source package as users are prohibited from redistribution of original or modified Tinker source code or binaries to other parties. While our intent is to distribute the Tinker code to anyone who wants it, the Ponder Lab would keep track of researchers using the package. The returned license forms also help us justify further development of Tinker. When new modules and capabilities become available, and when the inevitable bugs are uncovered, we will attempt to notify those who have returned a license form. Finally, we remind you that this software is copyrighted, and ask that it not be redistributed in any form.

Where to Direct Questions
-------------------------

Specific questions about the building or use of the Tinker package should be directed to ponder@dasher.wustl.edu. Tinker related questions or comments of more general interest can be posted on Twitter @TINKERtoolsMD. The Tinker developers monitor this account and will respond to the site or the individual poster as appropriate. 
