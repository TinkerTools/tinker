c
c
c     #########################################################
c     ##  COPYRIGHT (C) 2024 by Moses Chung & Jay W. Ponder  ##
c     ##                 All Rights Reserved                 ##
c     #########################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module alfmol  --  AlphaMol area and volume information  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     alfthread  number of OpenMP threads to use for AlphaMol2
c     delcxeps   eps value for AlphaMol Delaunay triangulation
c     alfhydro   logical flag to include hydrogen atoms in AlphaMol
c     alfsosgmp  logical flag governing use of SoS GMP in AlphaMol
c     alfmethod  set AlphaMol1 or AlphaMol2 (SINGLE or MULTI threaded)
c     alfsort    set sort method (NONE, SORT3D, BRIO, SPLIT, KDTREE)
c
c
      module alfmol
      implicit none
      integer alfthread
      real*8 delcxeps
      logical alfhydro
      logical alfsosgmp
      character*6 alfmethod
      character*6 alfsort
      save
      end
