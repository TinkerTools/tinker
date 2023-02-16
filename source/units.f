c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module units  --  physical constants and unit conversions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     literature references:
c
c     M. Stock, R. Davis, E. de Mirandes and M. J. T. Milton, "The
c     Revision of the SI - The Result of Three Decades of Progress
c     in Metrology", Metrologia, 56, 022001 (2019)
c
c     P. J. Mohr, D. B. Newell and B. N. Taylor, "CODATA Recommended
c     Values of the Fundamental Physical Constants: 2014", Journal of
c     Physical and Chemical Reference Data, 45, 043102 (2016)
c
c     Where appropriate, values are from the November 2018 revision
c     of SI units to fixed values by the 26th General Conference on
c     Weights and Measures; other values are taken from 2014 CODATA
c     reference constants or are described below
c
c     The conversion from calorie to Joule is the definition of the
c     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992)
c
c     The "coulomb" energy conversion factor is found by dimensional
c     analysis of Coulomb's Law, that is by dividing the square of the
c     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is
c     the permittivity of vacuum (the "electric constant"); note that
c     eps0 is typically given in F/m, equivalent to C**2/(J-m)
c
c     The approximate value used for the Debye, 3.33564 x 10-30 C-m,
c     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997)
c
c     The value of "prescon" is based on definition of 1 atmosphere
c     as 101325 Pa set by the 10th Conference Generale des Poids et
c     Mesures (Paris, 1954), where a Pascal (Pa) is equal to a J/m**3
c
c     avogadro    Avogadro's number (N) in particles/mole
c     lightspd    speed of light in vacuum (c) in cm/ps
c     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K
c     gasconst    ideal gas constant (R) in kcal/mole/K
c     elemchg     elementary charge of a proton in Coulombs
c     vacperm     vacuum permittivity (electric constant, eps0) in F/m
c     emass       mass of an electron in atomic mass units
c     planck      Planck's constant (h) in J-s
c     joule       conversion from calorie to joule
c     ekcal       conversion from kcal to g*Ang**2/ps**2
c     bohr        conversion from Bohr to Angstrom
c     hartree     conversion from Hartree to kcal/mole
c     evolt       conversion from Hartree to electron-volt
c     efreq       conversion from Hartree to cm-1
c     coulomb     conversion from electron**2/Ang to kcal/mole
c     elefield    conversion from electron**2/Ang to megavolt/cm
c     debye       conversion from electron-Ang to Debye
c     prescon     conversion from kcal/mole/Ang**3 to Atm
c
c
      module units
      implicit none
      real*8 avogadro
      real*8 lightspd
      real*8 boltzmann
      real*8 gasconst
      real*8 elemchg
      real*8 vacperm
      real*8 emass,planck
      real*8 joule,ekcal
      real*8 bohr,hartree
      real*8 evolt,efreq
      real*8 coulomb
      real*8 elefield
      real*8 debye,prescon
      parameter (avogadro=6.02214076d+23)
      parameter (lightspd=2.99792458d-2)
      parameter (boltzmann=0.831446262d0)
      parameter (gasconst=1.987204259d-3)
      parameter (elemchg=1.602176634d-19)
      parameter (vacperm=8.854187817d-12)
      parameter (emass=5.4857990907d-4)
      parameter (planck=6.62607015d-34)
      parameter (joule=4.1840d0)
      parameter (ekcal=4.1840d+2)
      parameter (bohr=0.52917721067d0)
      parameter (hartree=627.5094736d0)
      parameter (evolt=27.21138602d0)
      parameter (efreq=2.194746314d+5)
      parameter (coulomb=332.063713d0)
      parameter (elefield=1439.96455d0)
      parameter (debye=4.80321d0)
      parameter (prescon=6.85684112d+4)
      save
      end
