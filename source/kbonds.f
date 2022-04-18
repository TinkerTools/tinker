c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module kbonds  --  bond stretching forcefield parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxnb   maximum number of bond stretch parameter entries
c     maxnb5  maximum number of 5-membered ring bond stretch entries
c     maxnb4  maximum number of 4-membered ring bond stretch entries
c     maxnb3  maximum number of 3-membered ring bond stretch entries
c     maxnel  maximum number of electronegativity bond corrections
c
c     bcon    force constant parameters for harmonic bond stretch
c     bcon5   force constant parameters for 5-ring bond stretch
c     bcon4   force constant parameters for 4-ring bond stretch
c     bcon3   force constant parameters for 3-ring bond stretch
c     blen    bond length parameters for harmonic bond stretch
c     blen5   bond length parameters for 5-ring bond stretch
c     blen4   bond length parameters for 4-ring bond stretch
c     blen3   bond length parameters for 3-ring bond stretch
c     dlen    electronegativity bond length correction parameters
c     kb      string of atom classes for harmonic bond stretch
c     kb5     string of atom classes for 5-ring bond stretch
c     kb4     string of atom classes for 4-ring bond stretch
c     kb3     string of atom classes for 3-ring bond stretch
c     kel     string of atom classes for electronegativity corrections
c
c
      module kbonds
      implicit none
      integer maxnb
      integer maxnb5
      integer maxnb4
      integer maxnb3
      integer maxnel
      real*8, allocatable :: bcon(:)
      real*8, allocatable :: bcon5(:)
      real*8, allocatable :: bcon4(:)
      real*8, allocatable :: bcon3(:)
      real*8, allocatable :: blen(:)
      real*8, allocatable :: blen5(:)
      real*8, allocatable :: blen4(:)
      real*8, allocatable :: blen3(:)
      real*8, allocatable :: dlen(:)
      character*8, allocatable :: kb(:)
      character*8, allocatable :: kb5(:)
      character*8, allocatable :: kb4(:)
      character*8, allocatable :: kb3(:)
      character*12, allocatable :: kel(:)
      save
      end
