Installation
============

Force Field Explorer and Tinker are distributed via self-extracting installers for Linux, Apple MacOS and Microsoft Windows. Included within the installers are binary Tinker executables, source code, example/test files, documentation and the Java environment needed by Force Field Explorer. FFE contains a private JRE version that will not affect any other JREs installed system wide.
    
For each platform, a launch script/exe located in the top level installation directory named “ffe” will launch Force Field Explorer.
    
On some platforms, double clicking the ffe.jar file in the “tinker/jar” directory will also launch Force Field Explorer. In this case, your system JRE will be invoked instead of the private JRE bundled with the installation kit. This should behave correctly, as long as your system JRE is a 1.8 version (Java 8) and includes Java3D extensions.

Specific installation instructions for operating systems are provided below:

* Linux

	After downloading the installer, execute it from a shell window. A script to run Force Field Explorer is located in the top-level installation directory, which defaults to “$HOME/Tinker-FFE”. In this case, Tinker binaries would be located in the “$HOME/Tinker-FFE/tinker/bin” directory.

* Apple MacOS

	After downloading the installer, double click to start installation. Force Field Explorer can then be launched by double clicking on the Desktop icon.

* Microsoft Windows
   
	After downloading the installer, double click to start installation. Force Field Explorer can then be launched from an icon or the Start Menu.
