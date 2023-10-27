c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module potfit  --  values for electrostatic potential fit  ##
c     ##                                                             ##
c     #################################################################
c
c
c     nconf       total number of configurations to be analyzed
c     namax       maximum number of atoms in the largest configuration
c     ngatm       total atom number with active potential grid points
c     nfatm       total atom number in electrostatic potential fit
c     npgrid      total number of electrostatic potential grid points
c     ipgrid      atom associated with each potential grid point
c     wresp       weight used to restrain electrostatic parameters
c     xdpl0       target x-component of total dipole moment
c     ydpl0       target y-component of total dipole moment
c     zdpl0       target z-component of total dipole moment
c     xxqpl0      target xx-component of total quadrupole moment
c     xyqpl0      target xy-component of total quadrupole moment
c     xzqpl0      target xz-component of total quadrupole moment
c     yyqpl0      target yy-component of total quadrupole moment
c     yzqpl0      target yz-component of total quadrupole moment
c     zzqpl0      target zz-component of total quadrupole moment
c     fit0        initial value of each parameter used in potential fit      
c     fchg        partial charges by atom type during potential fit
c     fpol        atomic multipoles by atom type during potential fit
c     fcpen       charge penetration by atom type during potential fit
c     pgrid       Cartesian coordinates of potential grid points
c     epot        values of electrostatic potential at grid points
c     use_dpl     flag to include total dipole in potential fit
c     use_qpl     flag to include total quadrupole in potential fit
c     fit_mpl     flag for atomic monopoles to vary in potential fit
c     fit_dpl     flag for atomic dipoles to vary in potential fit
c     fit_qpl     flag for atomic quadrupoles to vary in potential fit
c     fit_chgpen  flag for atomic quadrupoles to vary in potential fit
c     fitchg      flag marking atom types used in partial charge fit
c     fitpol      flag marking atom types used in atomic multipole fit
c     fitcpen     flag marking atom types used in charge penetration
c     gatm        flag to use potential grid points around each atom
c     fatm        flag to use each atom in electrostatic potential fit
c     fxdpl       flag to use each atom x-dipole in electrostatic fit
c     fydpl       flag to use each atom y-dipole in electrostatic fit
c     fzdpl       flag to use each atom z-dipole in electrostatic fit
c     vchg        flag for partial charge at each atom in fitting
c     vpol        flag for atomic multipoles at each atom in fitting
c     vcpen       flag for charge penetration at each atom in fitting
c     resptyp     electrostatic restraint target (ORIG, ZERO or NONE)
c     varpot      descriptive name for each variable in potential fit
c
c
      module potfit
      use sizes
      implicit none
      integer nconf,namax
      integer ngatm,nfatm
      integer npgrid(maxref)
      integer, allocatable :: ipgrid(:,:)
      real*8 wresp
      real*8 xdpl0(maxref)
      real*8 ydpl0(maxref)
      real*8 zdpl0(maxref)
      real*8 xxqpl0(maxref)
      real*8 xyqpl0(maxref)
      real*8 xzqpl0(maxref)
      real*8 yyqpl0(maxref)
      real*8 yzqpl0(maxref)
      real*8 zzqpl0(maxref)
      real*8, allocatable :: fit0(:)
      real*8, allocatable :: fchg(:)
      real*8, allocatable :: fpol(:,:)
      real*8, allocatable :: fcpen(:)
      real*8, allocatable :: pgrid(:,:,:)
      real*8, allocatable :: epot(:,:,:)
      logical use_dpl,use_qpl
      logical fit_mpl,fit_dpl
      logical fit_qpl,fit_chgpen
      logical, allocatable :: fitchg(:)
      logical, allocatable :: fitpol(:)
      logical, allocatable :: fitcpen(:)
      logical, allocatable :: gatm(:)
      logical, allocatable :: fatm(:)
      logical, allocatable :: fxdpl(:)
      logical, allocatable :: fydpl(:)
      logical, allocatable :: fzdpl(:)
      logical, allocatable :: vchg(:,:)
      logical, allocatable :: vpol(:,:,:)
      logical, allocatable :: vcpen(:,:)
      character*4 resptyp
      character*6, allocatable :: varpot(:)
      save
      end
