This directory adds iEL-SCF functionality to Tinker.

iEL-SCF is a hybrid extended-Lagrangian/self-consistent field scheme originally developed by Anders Niklasson and colleauges to solve for the ground state density matrix in AIMD simulations.  iEL-SCF is the extension of this method to the solution of mutually induced dipoles.  'i' denotes 'interially-restrained', a way to limit the interia of the kinetic modes of the 'EL' part of iEL-SCF, which is necessary to maintain stability.  Several methods are possible for this, each taking the form of a familiar temperature control scheme from typical MD simulations.  More details will be available in an upcoming publication which will be listed here.


---
The use of iEL-SCF is invoked by using the 'iel-scf' keyword in a simulation keyfile.  The iEL-SCF specific keyworkds are listed below along with their definitions and uses:

'iel-scf' - Directs the TINKER MD simulation to use EL-SCF for solution of the mutually induced dipoles.

'auxstat [STRING]' - Informs the simulation to use inertial restraint using a [STRING] thermostat.  The currently available types for [STRING] are 'NOSE-HOOVER', 'NOSE-HOOVER1', 'BERENDSEN', and 'RESCALE'.  The default is 'NONE'.  'NOSE-HOOVER' uses a four-deep Nose-Hoover chain, 'NOSE-HOOVER1' uses a single Nose-Hoover chain, 'BERENDSEN' uses Berendsen rescaling, and 'RESCALE' uses straight periodic rescaling.  'NOSE-HOOVER' is recommended.

'aux-temp [REAL]' - Sets the auxiliary induced dipole 'temperature' at [REAL] in units of (e-Ang/ps)^2.  The default is 10^5 (e-Ang/ps)^2.

'aux-tau-temp [REAL]' - Sets the auxiliary 'thermostat' coupling timescale to [REAL] in units of ps.  Similar to 'tau-temperature' in default TINKER.  The default value is 0.1 ps.

---


A makefile is included in 'tinker/source' which includes added iEL-SCF files.  Change the default directories in that to your own and it should be good to go.

8/3/15:
This version should be good for water and ions.  There is some discrepancy between uinp and uind that may cause problems for more connected systems.  I'm currently working on this.

Email questions/comments/concerns to aalbaugh@berkeley.edu.


A. Albaugh, O. Demerdash, T. Head-Gordon
UC Berkeley
