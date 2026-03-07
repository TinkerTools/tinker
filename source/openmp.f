c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module openmp  --  OpenMP processor and thread values  ##
c     ##                                                         ##
c     #############################################################
c
c
c     nproc     number of processors available to OpenMP
c     nthread   number of threads to be used with OpenMP
c     nnest     number of nested active parallel regions
c
c
      module openmp
      implicit none
      integer nproc
      integer nthread
      integer nnest
      save
      end
