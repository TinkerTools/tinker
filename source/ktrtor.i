c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  ktrtor.i  --  forcefield parameters for torsion-torsions  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxntt    maximum number of torsion-torsion parameter entries
c     maxtgrd   maximum dimension of torsion-torsion spline grid
c     maxtgrd2  maximum number of torsion-torsion spline grid points
c
c     ttx       angle values for first torsion of spline grid
c     tty       angle values for second torsion of spline grid
c     tbf       function values at points on spline grid
c     tbx       gradient over first torsion of spline grid
c     tby       gradient over second torsion of spline grid
c     tbxy      Hessian cross components over spline grid
c     tnx       number of columns in torsion-torsion spline grid
c     tny       number of rows in torsion-torsion spline grid
c     ktt       string of torsion-torsion atom classes
c
c
      integer maxntt,maxtgrd,maxtgrd2
      parameter (maxntt=100)
      parameter (maxtgrd=30)
      parameter (maxtgrd2=maxtgrd*maxtgrd)
      integer tnx,tny
      real*8 ttx,tty,tbf
      real*8 tbx,tby,tbxy
      character*20 ktt
      common /ktrtor/ ttx(maxtgrd,maxntt),tty(maxtgrd,maxntt),
     &                tbf(maxtgrd2,maxntt),tbx(maxtgrd2,maxntt),
     &                tby(maxtgrd2,maxntt),tbxy(maxtgrd2,maxntt),
     &                tnx(maxntt),tny(maxntt),ktt(maxntt)
