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
c     "initatom" sets the atomic symbol, van der Waals and covalent
c     radii for each element in the periodic table
c
c     literature references:
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
      real*8 vrad(maxele)
      real*8 crad(maxele)
      character*3 asym(maxele)
c
c     atomic symbols for the elements
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
c     van der Waals radii for the elements (Angstroms)
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
c     covalent radii for the elements (Angstroms)
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
c     set atomic symbol and covalent radius for each element
c
      do i = 1, maxele
         if (vrad(i) .eq. 0.0d0)  vrad(i) = 2.0d0
         vdwrad(i) = vrad(i)
         covrad(i) = crad(i)
         elemnt(i) = asym(i)
      end do
      return
      end
