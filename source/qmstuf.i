c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2009 by Chuanjie Wu and Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  qmstuf.i  --  quantum data from Gaussian 03 calculation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     gx       x-coordinate of each atom in the QM data file
c     gy       y-coordinate of each atom in the QM data file
c     gz       z-coordinate of each atom in the QM data file
c     gforce   force components on each atom from QM data
c     gh       Hessian maxtrix elements from QM data
c     gfreq    calculated vibrational frequencies from QM data
c     ngatom   number of atoms in the QM data file
c
c
      integer ngatom
      real*8  gx,gy,gz
      real*8  gforce,gh,gfreq
      common /qmstuf/ gx(maxatm),gy(maxatm),gz(maxatm),gforce(3,maxatm),
     &                gh(maxhess),gfreq(maxvib),ngatom
