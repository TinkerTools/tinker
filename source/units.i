c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  units.i  --  physical constants and unit conversions  ##
c     ##                                                        ##
c     ############################################################
c
c
c     literature references:
c
c     P. J. Mohr, B. N. Taylor and D.B. Newell, "The 2010 CODATA
c     Recommended Values of the Fundamental Physical Constants",
c     National Institute of Standards & Technology, Gaithersburg, MD
c     (Web Version 6.0, developed jointly by J. Baker, M. Douma and
c     S. Kotochigova, available at http://physics.nist.gov/constants)
c
c     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended
c     Values of the Fundamental Physical Constants: 2006", Reviews of
c     Modern Physics, 80, 633-730 (2008)
c
c     Most values below are derived from 2010 CODATA reference values
c
c     The conversion from calorie to Joule is the definition of the
c     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992)
c
c     The "coulomb" energy conversion factor is found by dimensional
c     analysis of Coulomb's Law, ie, by dividing the square of the
c     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is
c     the permittivity of vacuum (the "electric constant"); note that
c     eps0 is typically given in F/m, equivalent to C**2/(J-m)
c
c     The approximate value used for the Debye, 3.33564 x 10-30 C-m,
c     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997)
c
c     The value of "prescon" is based on definition of 1 atmosphere
c     as 101325 Pa set by the 10th Conference Generale des Poids et
c     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3
c
c     avogadro    Avogadro's number (N) in particles/mole
c     lightspd    speed of light in vacuum (c) in cm/ps
c     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K
c     gasconst    ideal gas constant (R) in kcal/mole/K
c     emass       mass of an electron in atomic mass units
c     joule       conversion from calories to joules
c     convert     conversion from kcal to g*Ang**2/ps**2
c     bohr        conversion from Bohrs to Angstroms
c     hartree     conversion from Hartree to kcal/mole
c     evolt       conversion from Hartree to electron-volts
c     efreq       conversion from Hartree to cm-1
c     coulomb     conversion from electron**2/Ang to kcal/mole
c     debye       conversion from electron-Ang to Debyes
c     prescon     conversion from kcal/mole/Ang**3 to Atm
c
c
      real*8 avogadro,lightspd
      real*8 boltzmann,gasconst
      real*8 emass,joule
      real*8 convert,bohr
      real*8 hartree,evolt
      real*8 efreq,coulomb
      real*8 debye,prescon
      parameter (avogadro=6.02214129d+23)
      parameter (lightspd=2.99792458d-2)
      parameter (boltzmann=0.831446215d0)
      parameter (gasconst=1.98720415d-3)
      parameter (emass=5.485799095d-4)
      parameter (joule=4.1840d0)
      parameter (convert=4.1840d+2)
      parameter (bohr=0.52917721092d0)
      parameter (hartree=627.5094743d0)
      parameter (evolt=27.21138503d0)
      parameter (efreq=2.194746313708d+5)
      parameter (coulomb=332.063714d0)
      parameter (debye=4.80321d0)
      parameter (prescon=6.85684112d+4)
