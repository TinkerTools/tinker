c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  potfit.i  --  values for electrostatic potential fitting  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxpgrd   maximum dimension of electrostatic potential grid
c
c     pgrid     Cartesian coordinates of potential grid points
c     epot      values of electrostatic potential at grid points
c     xdpl0     target x-component of the molecular dipole moment
c     ydpl0     target y-component of the molecular dipole moment
c     zdpl0     target z-component of the molecular dipole moment
c     xxqdp0    target xx-component of the molecular quadrupole moment
c     xyqdp0    target xy-component of the molecular quadrupole moment
c     xzqdp0    target xz-component of the molecular quadrupole moment
c     yyqdp0    target yy-component of the molecular quadrupole moment
c     yzqdp0    target yz-component of the molecular quadrupole moment
c     zzqdp0    target zz-component of the molecular quadrupole moment
c     nconf     total number of conformations to be analyzed
c     npgrid    total number of electrostatic potential grid points
c     ipgrid    atom associated with each potential grid point
c     ngatm     total number of atoms with active potential grid points
c     nfatm     total number of atoms in electrostatic potential fit
c     gatm      flag to use potential grid points around each atom
c     fatm      flag to use each atom in electrostatic potential fit
c     use_dpl   flag to include molecular dipole in potential fit
c     use_qdp   flag to include molecular quadrupole in potential fit
c     fit_mpl   flag for atomic monopoles to vary in potential fit
c     fit_dpl   flag for atomic dipoles to vary in potential fit
c     fit_qdp   flag for atomic quadrupoles to vary in potential fit
c
c
      integer maxpgrd
      parameter (maxpgrd=100000)
      integer nconf
      integer npgrid,ipgrid
      integer ngatm,nfatm
      real*8 pgrid,epot
      real*8 xdpl0,ydpl0,zdpl0
      real*8 xxqdp0,xyqdp0,xzqdp0
      real*8 yyqdp0,yzqdp0,zzqdp0
      logical gatm,fatm
      logical use_dpl,use_qdp
      logical fit_mpl,fit_dpl
      logical fit_qdp
      common /potfit/ pgrid(3,maxpgrd,maxref),epot(2,maxpgrd,maxref),
     &                xdpl0(maxref),ydpl0(maxref),zdpl0(maxref),
     &                xxqdp0(maxref),xyqdp0(maxref),xzqdp0(maxref),
     &                yyqdp0(maxref),yzqdp0(maxref),zzqdp0(maxref),
     &                nconf,npgrid(maxref),ipgrid(maxpgrd,maxref),ngatm,
     &                nfatm,gatm(maxatm),fatm(maxatm),use_dpl,use_qdp,
     &                fit_mpl,fit_dpl,fit_qdp
