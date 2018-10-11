c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2018  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module polopt  --  induced dipoles via OPT extrapolation  ##        
c     ##                                                            ##
c     ################################################################
c
c
c     maxopt    maximum order for OPT induced dipole extrapolation
c
c     coptmax   maximum coefficient order for OPT dipole extrapolation
c     optlevel  current OPT order for reciprocal potential and field
c     copt      coefficients for OPT total induced dipole moments
c     copm      coefficients for OPT incremental induced dipole moments
c     uopt      OPT induced dipole components at each multipole site
c     uoptp     OPT induced dipoles in field used for energy terms
c     uopts     OPT GK or PB induced dipoles at each multipole site
c     uoptps    OPT induced dipoles in field used for GK or PB energy
c     fopt      OPT fractional reciprocal potentials at multipole sites
c     foptp     OPT fractional reciprocal potentials for energy terms
c
c
      module polopt
      implicit none
      integer maxopt
      parameter (maxopt=4)
      integer coptmax
      integer optlevel
      real*8, allocatable :: copt(:)
      real*8, allocatable :: copm(:)
      real*8, allocatable :: uopt(:,:,:)
      real*8, allocatable :: uoptp(:,:,:)
      real*8, allocatable :: uopts(:,:,:)
      real*8, allocatable :: uoptps(:,:,:)
      real*8, allocatable :: fopt(:,:,:)
      real*8, allocatable :: foptp(:,:,:)
      save
      end
