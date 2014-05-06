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
c     nconf     total number of conformations to be analyzed
c     ngatm     total number of atoms with active potential grid points
c     nfatm     total number of atoms in electrostatic potential fit
c     npgrid    total number of electrostatic potential grid points
c     ipgrid    atom associated with each potential grid point
c     xdpl0     target x-component of the molecular dipole moment
c     ydpl0     target y-component of the molecular dipole moment
c     zdpl0     target z-component of the molecular dipole moment
c     xxqdp0    target xx-component of the molecular quadrupole moment
c     xyqdp0    target xy-component of the molecular quadrupole moment
c     xzqdp0    target xz-component of the molecular quadrupole moment
c     yyqdp0    target yy-component of the molecular quadrupole moment
c     yzqdp0    target yz-component of the molecular quadrupole moment
c     zzqdp0    target zz-component of the molecular quadrupole moment
c     pgrid     Cartesian coordinates of potential grid points
c     epot      values of electrostatic potential at grid points
c     use_dpl   flag to include molecular dipole in potential fit
c     use_qdp   flag to include molecular quadrupole in potential fit
c     fit_mpl   flag for atomic monopoles to vary in potential fit
c     fit_dpl   flag for atomic dipoles to vary in potential fit
c     fit_qdp   flag for atomic quadrupoles to vary in potential fit
c     gatm      flag to use potential grid points around each atom
c     fatm      flag to use each atom in electrostatic potential fit
c
c
      module potfit
      use sizes
      implicit none
      integer nconf
      integer ngatm,nfatm
      integer npgrid(maxref)
      integer, allocatable :: ipgrid(:,:)
      real*8 xdpl0(maxref)
      real*8 ydpl0(maxref)
      real*8 zdpl0(maxref)
      real*8 xxqdp0(maxref)
      real*8 xyqdp0(maxref)
      real*8 xzqdp0(maxref)
      real*8 yyqdp0(maxref)
      real*8 yzqdp0(maxref)
      real*8 zzqdp0(maxref)
      real*8, allocatable :: pgrid(:,:,:)
      real*8, allocatable :: epot(:,:,:)
      logical use_dpl,use_qdp
      logical fit_mpl,fit_dpl
      logical fit_qdp
      logical, allocatable :: gatm(:)
      logical, allocatable :: fatm(:)
      save
      end
