c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module ktrtor  --  torsion-torsion forcefield parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxntt    maximum number of torsion-torsion parameter entries
c     maxtgrd   maximum dimension of torsion-torsion spline grid
c     maxtgrd2  maximum number of torsion-torsion spline grid points
c
c     tnx       number of columns in torsion-torsion spline grid
c     tny       number of rows in torsion-torsion spline grid
c     ttx       angle values for first torsion of spline grid
c     tty       angle values for second torsion of spline grid
c     tbf       function values at points on spline grid
c     tbx       gradient over first torsion of spline grid
c     tby       gradient over second torsion of spline grid
c     tbxy      Hessian cross components over spline grid
c     ktt       string of torsion-torsion atom classes
c
c
      module ktrtor
      implicit none
      integer maxntt
      integer maxtgrd
      integer maxtgrd2
      integer, allocatable :: tnx(:)
      integer, allocatable :: tny(:)
      real*8, allocatable :: ttx(:,:)
      real*8, allocatable :: tty(:,:)
      real*8, allocatable :: tbf(:,:)
      real*8, allocatable :: tbx(:,:)
      real*8, allocatable :: tby(:,:)
      real*8, allocatable :: tbxy(:,:)
      character*20, allocatable :: ktt(:)
      save
      end
