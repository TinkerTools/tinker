Here are a few modifications to TINKER that provide hooks into PLUMED2. In the interest of collaboration, I am sharing this in a very ugly (but working) state. 

I have not yet run an enhanced free energy simulation to convergence, but I am running a 2D metdynamics simulation of the alanine dipeptide in water (AMOEBA-BIO-2013 force field), using the Nose-Hoover NPT integrator from TINKER. The files required to do so are in the folder aladi/

INSTALLATION:


1. Install Plumed2 (https://github.com/plumed/plumed2)


2. The file Plumed.inc will be generated.  Copy it into the source/ directory of (this version of) TINKER.


3. (Go back to TINKER directory) Make your machine-specific edits to the Makefile in make/, and copy Makefile into /source



4. From source, 'make.'

Have fun and report errors/bugs!
