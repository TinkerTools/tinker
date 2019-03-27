c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine initatom  --  setup atoms in periodic table  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "initatom" sets the atomic symbol, standard atomic weight,
c     van der Waals radius and covalent radius for each element in
c     the periodic table
c
c     literature references:
c
c     J. Emsley, The Elements, 3rd Edition, Oxford University Press,
c     (1999)  [relative atomic masses]
c
c     J. Meija, T. B. Coplen, M. Berglund, W. A. Brand, P. De Bievre,
c     M. Groning, N. E. Holden, J. Irrgeher, R. D. Loss, T. Walczyk and
c     R. Prohaska, Atomic Weights of the Elements 2013, Pure and Applied
c     Chemistry, 88, 265-291 (2016)  [standard atomic weights]
c
c     A. Bondi, van der Waals Volumes and Radii, Journal of Physical
c     Chemistry, 68, 441-451 (1964)  [original vdw radii; not used]
c
c     S. Alvarez, "A Cartography of the van der Waals Territories",
c     Dalton Transactions, 42, 8617-8636 (2013)  [vdw radii for most
c     elements 1-99]
c
c     B. Cordero, V. Gomez. A. E. Platero-Prats, M. Reves,
c     J. Echeverria, E. Cremades, F. Barragan and S. Alverez,
c     "Covalent Radii Revisited", Dalton Transactions, 2832-2838 (2008)
c     [covalent radii for elements 1-96]
c
c     P. Pyykko and M. Atsumi, "Molecular Single-Bond Covalent Radii
c     for Elements 1-118", Chemistry- A European Journal, 15, 187-197
c     (2009)  [covalent radii for elements 97-112]
c
c
      subroutine initatom
      use ptable
      implicit none
      integer i
      real*8 amas(maxele)
      real*8 vrad(maxele)
      real*8 crad(maxele)
      character*3 asym(maxele)
c
c     atomic symbol for each element
c
      data asym  / 'H  ', 'He ', 'Li ', 'Be ', 'B  ', 'C  ', 'N  ',
     &             'O  ', 'F  ', 'Ne ', 'Na ', 'Mg ', 'Al ', 'Si ',
     &             'P  ', 'S  ', 'Cl ', 'Ar ', 'K  ', 'Ca ', 'Sc ',
     &             'Ti ', 'V  ', 'Cr ', 'Mn ', 'Fe ', 'Co ', 'Ni ',
     &             'Cu ', 'Zn ', 'Ga ', 'Ge ', 'As ', 'Se ', 'Br ',
     &             'Kr ', 'Rb ', 'Sr ', 'Y  ', 'Zr ', 'Nb ', 'Mo ',
     &             'Tc ', 'Ru ', 'Rh ', 'Pd ', 'Ag ', 'Cd ', 'In ',
     &             'Sn ', 'Sb ', 'Te ', 'I  ', 'Xe ', 'Cs ', 'Ba ',
     &             'La ', 'Ce ', 'Pr ', 'Nd ', 'Pm ', 'Sm ', 'Eu ',
     &             'Gd ', 'Tb ', 'Dy ', 'Ho ', 'Er ', 'Tm ', 'Yb ',
     &             'Lu ', 'Hf ', 'Ta ', 'W  ', 'Re ', 'Os ', 'Ir ',
     &             'Pt ', 'Au ', 'Hg ', 'Tl ', 'Pb ', 'Bi ', 'Po ',
     &             'At ', 'Rn ', 'Fr ', 'Ra ', 'Ac ', 'Th ', 'Pa ',
     &             'U  ', 'Np ', 'Pu ', 'Am ', 'Cm ', 'Bk ', 'Cf ',
     &             'Es ', 'Fm ', 'Md ', 'No ', 'Lr ', 'Rf ', 'Db ',
     &             'Sg ', 'Bh ', 'Hs ', 'Mt ', 'Ds ', 'Rg ', 'Cn ' /
c
c     standard atomic weight for each element
c
      data amas  /  1.008d0,   4.003d0,   6.941d0,   9.012d0,  10.811d0,
     &             12.011d0,  14.007d0,  15.999d0,  18.998d0,  20.180d0,
     &             22.990d0,  24.305d0,  26.982d0,  28.086d0,  30.974d0,
     &             32.066d0,  35.453d0,  39.948d0,  39.098d0,  40.078d0,
     &             44.956d0,  47.867d0,  50.942d0,  51.996d0,  54.938d0,
     &             55.845d0,  58.933d0,  58.693d0,  63.546d0,  65.380d0,
     &             69.723d0,  72.630d0,  74.922d0,  78.971d0,  79.904d0,
     &             83.798d0,  85.468d0,  87.620d0,  88.906d0,  91.224d0,
     &             92.906d0,  95.950d0,  98.906d0, 101.070d0, 102.910d0,
     &            106.420d0, 107.870d0, 112.410d0, 114.820d0, 118.710d0,
     &            121.760d0, 127.600d0, 126.900d0, 131.290d0, 132.910d0,
     &            137.330d0, 138.910d0, 140.120d0, 140.910d0, 144.240d0,
     &            144.913d0, 150.360d0, 151.960d0, 157.250d0, 158.930d0,
     &            162.500d0, 164.930d0, 167.260d0, 168.930d0, 173.050d0,
     &            174.970d0, 178.490d0, 180.950d0, 183.840d0, 186.210d0,
     &            190.230d0, 192.220d0, 195.080d0, 196.970d0, 200.590d0,
     &            204.383d0, 207.200d0, 208.980d0, 208.982d0, 209.987d0,
     &            222.017d0, 223.020d0, 226.025d0, 227.027d0, 232.038d0,
     &            231.036d0, 238.029d0, 237.048d0, 244.064d0, 243.061d0,
     &            247.070d0, 247.070d0, 251.080d0, 252.083d0, 257.095d0,
     &            258.098d0, 259.101d0, 262.110d0, 267.122d0, 270.131d0,
     &            269.129d0, 270.133d0, 270.134d0, 278.156d0, 281.165d0,
     &            281.166d0, 285.177d0 /
c
c     van der Waals radius for each element (Angstroms)
c
      data vrad  / 1.20d0, 1.43d0, 2.12d0, 1.98d0, 1.91d0, 1.77d0,
     &             1.66d0, 1.50d0, 1.46d0, 1.58d0, 2.50d0, 2.51d0,
     &             2.25d0, 2.19d0, 1.90d0, 1.89d0, 1.82d0, 1.83d0,
     &             2.73d0, 2.62d0, 2.58d0, 2.46d0, 2.42d0, 2.45d0,
     &             2.45d0, 2.44d0, 2.40d0, 2.40d0, 2.38d0, 2.39d0,
     &             2.32d0, 2.29d0, 1.88d0, 1.82d0, 1.86d0, 2.25d0,
     &             3.21d0, 2.84d0, 2.75d0, 2.52d0, 2.56d0, 2.45d0,
     &             2.44d0, 2.46d0, 2.44d0, 2.15d0, 2.53d0, 2.49d0,
     &             2.43d0, 2.42d0, 2.47d0, 1.99d0, 2.04d0, 2.06d0,
     &             3.48d0, 3.03d0, 2.98d0, 2.88d0, 2.92d0, 2.95d0,
     &             0.00d0, 2.90d0, 2.87d0, 2.83d0, 2.79d0, 2.87d0,
     &             2.81d0, 2.83d0, 2.79d0, 2.80d0, 2.74d0, 2.63d0,
     &             2.53d0, 2.57d0, 2.49d0, 2.48d0, 2.41d0, 2.29d0,
     &             2.32d0, 2.45d0, 2.47d0, 2.60d0, 2.54d0, 0.00d0,
     &             0.00d0, 0.00d0, 0.00d0, 0.00d0, 2.80d0, 2.93d0,
     &             2.88d0, 2.71d0, 2.82d0, 2.81d0, 2.83d0, 3.05d0,
     &             3.40d0, 3.05d0, 2.70d0, 0.00d0, 0.00d0, 0.00d0,
     &             0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,
     &             0.00d0, 0.00d0, 0.00d0, 0.00d0 /
c
c     covalent radius for each element (Angstroms)
c
      data crad  / 0.31d0, 0.28d0, 1.28d0, 0.96d0, 0.84d0, 0.76d0,
     &             0.71d0, 0.66d0, 0.57d0, 0.58d0, 1.66d0, 1.41d0,
     &             1.21d0, 1.11d0, 1.07d0, 1.05d0, 1.02d0, 1.06d0,
     &             2.03d0, 1.76d0, 1.70d0, 1.60d0, 1.53d0, 1.39d0,
     &             1.39d0, 1.32d0, 1.26d0, 1.24d0, 1.32d0, 1.22d0,
     &             1.22d0, 1.20d0, 1.19d0, 1.20d0, 1.20d0, 1.16d0,
     &             2.20d0, 1.95d0, 1.90d0, 1.75d0, 1.64d0, 1.54d0,
     &             1.47d0, 1.46d0, 1.42d0, 1.39d0, 1.45d0, 1.44d0,
     &             1.42d0, 1.39d0, 1.39d0, 1.38d0, 1.39d0, 1.40d0,
     &             2.44d0, 2.15d0, 2.07d0, 2.04d0, 2.03d0, 2.01d0,
     &             1.99d0, 1.98d0, 1.98d0, 1.96d0, 1.94d0, 1.92d0,
     &             1.92d0, 1.89d0, 1.90d0, 1.87d0, 1.87d0, 1.75d0,
     &             1.70d0, 1.62d0, 1.51d0, 1.44d0, 1.41d0, 1.36d0,
     &             1.36d0, 1.32d0, 1.45d0, 1.46d0, 1.48d0, 1.40d0,
     &             1.50d0, 1.50d0, 2.60d0, 2.21d0, 2.15d0, 2.06d0,
     &             2.00d0, 1.96d0, 1.90d0, 1.87d0, 1.80d0, 1.69d0,
     &             1.68d0, 1.68d0, 1.65d0, 1.67d0, 1.73d0, 1.76d0,
     &             1.61d0, 1.57d0, 1.49d0, 1.43d0, 1.41d0, 1.34d0,
     &             1.29d0, 1.28d0, 1.21d0, 1.22d0 /
c
c
c     set the symbol, weight and radii for each element
c
      do i = 1, maxele
         atmass(i) = amas(i)
         elemnt(i) = asym(i)
         if (vrad(i) .eq. 0.0d0)  vrad(i) = 2.0d0
         vdwrad(i) = vrad(i)
         covrad(i) = crad(i)
      end do
      return
      end
