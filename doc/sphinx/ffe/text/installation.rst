Installation
============

Use of Self-Extracting Kits
---------------------------

Force Field Explorer and Tinker are distributed via self-extracting installers for Linux, Apple MacOS and Microsoft Windows. These install kits are available from the Tinker web site maintained by the Ponder group at https://dasher.wustl.edu/tinker/. Included within the installers are binary Tinker executables, source code, example/test files, documentation and the Java environment needed by Force Field Explorer. FFE contains a private JRE version that will not affect any other JREs installed system-wide.
    
For each platform, a launch script/exe located in the top level installation directory named “ffe” will launch Force Field Explorer.
    
On some platforms, double clicking the ffe.jar file in the “tinker/jar” directory will also launch Force Field Explorer. In this case, your system JRE will be invoked instead of the private JRE bundled with the installation kit. This should behave correctly, as long as your system JRE is a 1.8 version (Java 8) and includes Java3D extensions.

Specific installation instructions for operating systems are provided below:

* Linux

	The Tinker-FFE Installer for Linux is provided as a gzipped shell script. Uncompress the the .gz archive to produce an .sh script, and then run the script. The script must have the "executable" attribute, set via "chmod +x ffe-linux.sh", prior to being run.

	After installation, the top-level directory defaults to “$HOME/Tinker-FFE”. The FFE-enabled Tinker binaries are located in the “$HOME/Tinker-FFE/tinker/bin” directory.

* Apple MacOS

	The Tinker-FFE Installer for MacOS is provided as a .dmg disk image file. Double-click on the file to run the installer. Once installed, Force Field Explorer is launched by double clicking on the Desktop icon.

	Recent versions of MacOS contains a security feature called Gatekeeper that keeps applications not obtained via the App Store or Apple-approved developers from being opened. Gatekeeper is enabled by default, and may result in the (incorrect!) error message: "Tinker-FFE Installer.app is damaged and can't be opened." To turn off Gatekeeper, go to the panel System Preferences > Security & Privacy > General, and set "Allow apps downloaded from:" to "Anywhere". This will require an Administrator account, and must be done before invoking the FFE installer. After FFE is installed and launched for the first time, you can return the System Preference to its prior value. On Sierra (MacOS 10.12) and later, the "Anywhere" option has been removed. In most cases the Security & Privacy panel will open and permit the user to run the installer. Alternatively, the "Anywhere" option can be restored by running the command "sudo spctl --master-disable" in a Terminal window with Administrator privilege.

* Microsoft Windows
   
	The Tinker-FFE Installer for Windows is provided as a zipped executable. First, unzip the .zip file. Then run the resulting executable .exe file by double clicking. Following installation, Force Field Explorer is launched via the FFE desktop icon or from the Start Menu.

	In order to perform minimizations or molecular dynamics from within FFE, some environment variables and symbolic links must be set prior to using the program. A batch file named "FFESetupWin.bat" is installed in the main Tinker-FFE directory, which by default resides in the user's home directory. To complete the setup of FFE, this batch file should be run from a Command Prompt window following installation. It is only necessary to invoke this batch file once, as the settings should persist between logins.

Building FFE from Source
------------------------

For those wishing to modify the FFE GUI or build installation kits from source, we provide a complete development package for Tinker-FFE. The full source distribution for Force Field Explorer is available from the TinkerTools organization public repository on Github located at https://github/TinkerTools/FFE/. The 0README file in the top-level of the distribution has instructions for building installers for Linux, MacOS and Windows.

Additional requirements for building FFE are the GNU compiler suite with gcc, g++ and gfortran and the Install4j Java installer builder. On Windows systems we recommend use of the MinGW-w64 compilers under Cygwin. Note that Install4j is a commercial product; only the command line "install4jc" compiler is needed, and not the full Install4j GUI interface.
