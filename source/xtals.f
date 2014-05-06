c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module xtals  --  crystals used for parameter fitting  ##
c     ##                                                         ##
c     #############################################################
c
c
c     maxlsq       maximum number of least squares variables
c     maxrsd       maximum number of least squares residuals
c
c     nxtal        number of crystal structures to be stored
c     nvary        number of potential parameters to optimize
c     ivary        index for the types of potential parameters
c     vary         atom numbers involved in potential parameters
c     iresid       crystal structure to which each residual refers
c     e0_lattice   ideal lattice energy for the current crystal
c     moment_0     ideal dipole moment for monomer from crystal
c     rsdtyp       experimental variable for each of the residuals
c     vartyp       type of potential parameter to be optimized
c
c
      module xtals
      implicit none
      integer maxlsq,maxrsd
      parameter (maxlsq=50)
      parameter (maxrsd=100)
      integer nxtal,nvary
      integer ivary(maxlsq)
      integer vary(2,maxlsq)
      integer iresid(maxrsd)
      real*8 e0_lattice
      real*8 moment_0
      character*16 rsdtyp(maxrsd)
      character*16 vartyp(maxlsq)
      save
      end
