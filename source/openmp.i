c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  openmp.i  --  system parameters for OpenMP computation  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nproc     number of processors available to OpenMP
c     nthread   number of threads to be used with OpenMP 
c
c
      integer nproc,nthread
      common /openmp/ nproc,nthread
