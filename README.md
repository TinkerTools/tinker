This code is capable of resonance-free multiple step (MTS) dynamics via the stochastic isokinetic integrator detailed in Leimkuhler, Margul, and Tuckerman Mol. Phys. 111, 22-23, 3579-3594 (2013).  This paper details the application to a fixed-charge model, but this code supports polarizable force fields.  It has been tested with the AMOEBA-WATER-2003 model, simulating a periodic cube of 512 H2O molecules.

With properly chosen keywords, this code can generate canonical position distributions functions with well over an order of magnitude of computation time savings. Consider the Oxygen-Oxygen, Oxygen-Hydrogen, and Hydrogen-Hydrogen distribution functions of liquid water. Reference NVT data are shown in blue, and MTS data with a large time step of 105 fs are in red. Run the simulation yourself in xi-respa-test/:

![image](https://cloud.githubusercontent.com/assets/4325206/8332087/f47d509c-1a59-11e5-89be-9405b48b3039.png)


The following keywords are required to run a stochastic isokinetic MTS simulation:

 INTEGRATOR   STOCH-RESPA
 
 NRES-BOND [integer]
 
 NRES-TORS   [integer]
 
 NRES-SHORT [integer]
 
- For a single large outer time step such as 105 fs: 
  - The long-range non-bonded interactions are computed once. 
  - The short-range non-bonded interactions are computed NRES-SHORT times.  
  - The torsional terms are computed (NRES-SHORT x NRES-TORS) times.  
  - The bonded interactions are computed (NRES-SHORT x NRES-TORS x NRES-BOND) times. 
- The default value for all three NRES- integers is 1.

NEIGHBOR-LIST

RESPA-CUTOFF  [value, Å]

RESPA-TAPER    [value, Å]

 - RESPA-CUTOFF is the maximum distance between a pair of interacting particles which the user would classify as within "short-range."  This applies to the van der Waals interaction and either fixed charged PME or AMOEBA-like polarizable electrostatics.  For more control, the user may specify different short-range cutoffs for different interactions with the keywords:
  - RESPA-VDW-CUTOFF [value, Å] for van der Waals interactions 
  - RESPA-CHG-CUTOFF [value, Å] for fixed-charge electrostatics 
  - RESPA-MPOLE-CUTOFF [value, Å] for polarizable electrostatics
 - RESPA-TAPER represents the minimum distance between a pair of interacting particles for which the user would like to apply the quintic force switching function in Morrone, Zhou, and Berne, J. Chem. Theory Comput. 6, 1798-1804 (2010).  For more control, the user may specify this value for different interactions with the keywords:
  - RESPA-VDW-TAPER [value, Å] for van der Waals interactions 
  - RESPA-CHG-TAPER [value, Å] for fixed-charge electrostatics 
  - RESPA-MPOLE-TAPER [value, Å] for polarizable electrostatics

NHC-LENGTH <integer>
 - In Leimkuhler et al., this is L, the number of thermostat particles per degree of freedom. It has a default value of 4. For XI-RESPA (see below), the user may find that L = 1 is sufficient. 

XO-RESPA
 - By default, the integrator will follow an XI-RESPA scheme, where the Nose-like update to thermostat velocities occurs during the same time step level as the bonded interactions.  The user may enter the keyword XO-RESPA instead, so that this occurs during the time step level of the long-range nonbonded interactions.  This will speed up the code with some accuracy loss, and accuracy can be regained with larger values of 

TAU-TEMPERATURE [value, ps]
 - This is a characteristic timescale of the physical system of interest, and it is used to assign "masses" to the extended phase space thermostat particles. This variable is native to TINKER for use in Nose-like thermostats and barostats, and has a default value of 0.2 ps.  For stochastic isokinetic dynamics of liquid water, TAU-TEMPERATURE of 0.01 ps is recommended. 

STOCH-GAMMA [value, ps^-1]
 - This is the Langevin friction coefficient $\gamma$ in Leimkuhler et al. It has a default value of 100 ps^-1.

RESPA-THERM-NC   [integer]

RESPA-THERM-NSY [integer]

 - The user may fine tune the Suzuki-Yoshida decomposition of the Nose-like thermostat velocity update. RESPA-THERM-NC is the number of factorizations per update, and RESPA-THERM-NSY is the number of Suzuki-Yoshida time steps per factorization. 
 - By default, RESPA-THERM-NC = 5. For XI-RESPA, RESPA-THERM-NC = 2 may be sufficient. 
 - RESPA-THERM-NSY has a default value of 3, and will accept 3, 5, 7, or 15.

