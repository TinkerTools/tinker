c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  omega.i  --  dihedrals for torsional space computations  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     dihed    current value in radians of each dihedral angle
c     nomega   number of dihedral angles allowed to rotate
c     iomega   numbers of two atoms defining rotation axis
c     zline    line number in Z-matrix of each dihedral angle
c
c
      integer nomega,iomega,zline
      real*8 dihed
      common /omega/ dihed(maxrot),nomega,iomega(2,maxrot),
     &               zline(maxrot)
